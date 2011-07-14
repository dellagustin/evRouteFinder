#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "../reuse/ge_points.h"
#include "../reuse/rt_array_index_shuffler.h"

using namespace std;

// switches
#define _LOGGING_

class CWayPoint : public CGePoint2d
{
public:
	CWayPoint();
	// ~CWayPoint();

	// Enums
public:
    enum IdTypes
    {
        idUnknown = 0x0000,
    };

    // Types
public:
    typedef unsigned int Id;

    // Member Access
public:
    Id id() const;

    // Operator
public:
    bool operator==(const CWayPoint&) const;

    // Members
public:
	Id m_Id;
};

CWayPoint::CWayPoint()
{
    m_Id = idUnknown;
}

CWayPoint::Id CWayPoint::id() const
{
    return m_Id;
}

bool CWayPoint::operator==(const CWayPoint& otherWayPoint) const
{
    // for a while lets use just the id for comparison
    if(id() != otherWayPoint.id())
        return false;

    return true;
}

typedef vector<CWayPoint> CWayPointArray;

// colection of waypoints
class CRoute : public vector<CWayPoint>
{
public:
	CRoute();

	// Query Methods
public:
    bool findWayPoint(CWayPoint::Id wayPointId, unsigned int& uiWayPointIndex) const;
	double linearLength() const;
	unsigned int maxWayPoints() const;

    // Editing Methods
public:
	void setToMaxWayPoints();
	void setMaxWayPoints(unsigned int newMaxWayPoints);
	void copyPropertiesFrom(const CRoute& otherRoute);

	// Operators
public:
    bool operator==(const CRoute&) const;

private:
	unsigned int m_uiMaxWayPoints;
};

CRoute::CRoute()
{
    m_uiMaxWayPoints = 0;
}

bool CRoute::findWayPoint(CWayPoint::Id wayPointId, unsigned int& uiWayPointIndex) const
{
    unsigned int i;

    for(i = 0; i < size(); i++)
    {
        if((*this)[i].id() == wayPointId)
        {
            uiWayPointIndex = i;
            return true;
        }
    }

    return false;
}

double CRoute::linearLength() const
{
	unsigned int i;
	double fLinearLength = 0;
	CWayPoint lastValidWayPoint;

	for(i = 0; i < size(); i++)
	{
	    // skip invalid waypoints
	    if((*this)[i].id())
	    {
	        if(lastValidWayPoint.id())
	        {
                fLinearLength += (*this)[i].distanceTo(lastValidWayPoint);
	        }

            // keep track of the last valid waypoint
	        lastValidWayPoint = (*this)[i];
	    }
	}

	return fLinearLength;
}

unsigned int CRoute::maxWayPoints() const
{
    return m_uiMaxWayPoints;
}

void CRoute::setToMaxWayPoints()
{
    resize(m_uiMaxWayPoints );
}

void CRoute::setMaxWayPoints(unsigned int newMaxWayPoints)
{
    m_uiMaxWayPoints = newMaxWayPoints;
}

void CRoute::copyPropertiesFrom(const CRoute& otherRoute)
{
    setMaxWayPoints(otherRoute.maxWayPoints());
    resize(otherRoute.size());
}

bool CRoute::operator==(const CRoute& otherRoute) const
{
    if(m_uiMaxWayPoints != otherRoute.m_uiMaxWayPoints)
        return false;

    unsigned int i;
    for(i = 0; i < size(); i++)
    {
        if(!((*this)[i] == otherRoute[i]))
            return false;
    }

    return true;
}

// CRouteArray {{

class CRouteArray : public vector<CRoute>
{
public:
	CRouteArray();
	// ~CRouteArray();

    // Query Methods
public:
    //! Verify if this route array is compatible with route array
    //! otherRouteArray
    bool isCompatibleWith(const CRouteArray& otherRouteArray) const;
	double linearLength() const;
	double totalTravelDistanceCached() const;
	bool findWayPoint(CWayPoint::Id id, unsigned int &uiRouteIndex, unsigned int &uiWayPointIndex) const;
	bool indexValidity(unsigned int uiRouteIndex, unsigned int uiWayPointIndex) const;

	// Tool Methods
public:
    bool incrementIndexes(unsigned int &uiRouteIndex, unsigned int &uiWayPointIndexAtRoute) const;

    // Edit Methods
public:
    bool setWayPoint(const CWayPoint& wayPoint, unsigned int uiRouteIndex, unsigned int uiWayPointIndex, bool bOverride = false);
    void copyPropertiesFrom(const CRouteArray& routeArray1);
    double updateTotalTravelDistanceCache();

public:
    bool operator< (const CRouteArray&) const;
    bool operator== (const CRouteArray&) const;

private:
    double m_fTotalTravelDistanceCache;
    // bool m_bUseCachedTotalTravelDistance;
};

typedef vector<CRouteArray> CRouteArrayArray;

CRouteArray::CRouteArray()
{
    m_fTotalTravelDistanceCache = 0;
}

bool CRouteArray::isCompatibleWith(const CRouteArray& otherRouteArray) const
{
    if(size() != otherRouteArray.size())
        return false;

    unsigned int i;
    for(i = 0; i < size(); i++)
    {
        if((*this)[i].size() != otherRouteArray[i].size())
            return false;
    }

    return true;
}

double CRouteArray::linearLength() const
{
	double fRoutesLinearLength = 0.0;

	unsigned int i;
	for(i = 0; i < size(); i++)
	{
		fRoutesLinearLength += (*this)[i].linearLength();
	}

	return fRoutesLinearLength;
}

double CRouteArray::updateTotalTravelDistanceCache()
{
    m_fTotalTravelDistanceCache = linearLength();
    return m_fTotalTravelDistanceCache;
}

double CRouteArray::totalTravelDistanceCached() const
{
    return m_fTotalTravelDistanceCache;
}

bool CRouteArray::incrementIndexes(unsigned int &uiRouteIndex, unsigned int &uiWayPointIndexAtRoute) const
{
    if (uiRouteIndex >= size())
        return false;

    if (uiWayPointIndexAtRoute >= (*this)[uiRouteIndex].size())
        return false;

    uiWayPointIndexAtRoute++;

    if (uiWayPointIndexAtRoute >= (*this)[uiRouteIndex].size())
    {
        // reset
        uiWayPointIndexAtRoute = 0;

        uiRouteIndex++;

        if (uiWayPointIndexAtRoute >= size())
        {
            uiRouteIndex = 0;
        }
    }

    return true;
}

bool CRouteArray::findWayPoint(CWayPoint::Id id, unsigned int &uiRouteIndex, unsigned int &uiWayPointIndex) const
{
    unsigned int i;
    for(i = 0; i < size(); i++)
    {
        if((*this)[i].findWayPoint(id, uiWayPointIndex))
        {
            uiRouteIndex = i;
            return true;
        }
    }

    return false;
}

bool CRouteArray::indexValidity(unsigned int uiRouteIndex, unsigned int uiWayPointIndex) const
{
    if(uiRouteIndex >= size())
        return false;

    if(uiWayPointIndex >= (*this)[uiRouteIndex].size())
        return false;

    return true;
}

bool CRouteArray::setWayPoint(const CWayPoint& wayPoint, unsigned int uiRouteIndex, unsigned int uiWayPointIndex, bool bOverride)
{
    if(!indexValidity(uiRouteIndex, uiWayPointIndex))
        return false;

    if((*this)[uiRouteIndex][uiWayPointIndex].id() && (!bOverride))
        return false;

    (*this)[uiRouteIndex][uiWayPointIndex] = wayPoint;

    return true;
}

void CRouteArray::copyPropertiesFrom(const CRouteArray& otherRouteArray)
{
    resize(otherRouteArray.size());

    unsigned int i;
    for(i = 0; i < size(); i++)
    {
        (*this)[i].copyPropertiesFrom(otherRouteArray[i]);
    }
}

bool CRouteArray::operator< (const CRouteArray& otherRouteArray) const
{
    return totalTravelDistanceCached() < otherRouteArray.totalTravelDistanceCached();
}

bool CRouteArray::operator== (const CRouteArray& otherRouteArray) const
{
    unsigned int i;
    for(i = 0; i < size(); i++)
    {
        if(!((*this)[i] == otherRouteArray[i]))
        {
            return false;
        }
    }

    return true;
}

// }} CRouteArray

void filter_best_samples(CRouteArrayArray& routeArrayArray, unsigned int uiSamples)
{
    /*
    bool bContinue;

    // simple bubble sort
    do
    {
        bContinue = false;

        unsigned int i;
        for(i = 0; i < routeArrayArray.size()-1; i++)
        {
            if(routeArrayArray[i].totalTravelDistanceCached() > routeArrayArray[i+1].totalTravelDistanceCached())
            {
                bContinue = true;
                swap(routeArrayArray[i], routeArrayArray[i+1]);
                // routeArrayArray[i].swap(routeArrayArray[i+1]);
            }
        }
    }
    while(bContinue);

    routeArrayArray.resize(uiSamples);
    */
    sort(routeArrayArray.begin(), routeArrayArray.end());

    // remove repeated samples
    unsigned int i;
    for(i = 0; i < routeArrayArray.size()-1; i++)
    {
        if(routeArrayArray[i].totalTravelDistanceCached() == routeArrayArray[i+1].totalTravelDistanceCached())
        {
            routeArrayArray.erase(routeArrayArray.begin() + i + 1);
            i--;
#ifdef _LOGGING_
            fprintf(stderr, "Twin removed.\n");
#endif
        }
    }

    routeArrayArray.resize(uiSamples);
}

int create_child(const CWayPointArray& wayPointArray, const CRouteArray& routeArray1, const CRouteArray& routeArray2, CRouteArray& childRouteArray)
{
    unsigned int uiWayPointIndex;

    if(!routeArray1.isCompatibleWith(routeArray2))
        return 0;

    CArrayIndexShuffler arrayIndexShuffler(wayPointArray.size());

    // clear and set properties of the child route array
    childRouteArray = CRouteArray();
    childRouteArray.copyPropertiesFrom(routeArray1);

    while((uiWayPointIndex = arrayIndexShuffler.nextIndex()) != (unsigned int)(-1))
    {
        bool bRetValue;

        unsigned int uiRouteIndex, uiWayPointIndexAtRoute;
        unsigned int uiOldRouteIndex, uiOldWayPointIndexAtRoute;
        const CRouteArray *pCurParentRouteAr = rand() & 0x01 ? &routeArray2 : &routeArray1;

        bRetValue = pCurParentRouteAr->findWayPoint(wayPointArray[uiWayPointIndex].id(), uiRouteIndex, uiWayPointIndexAtRoute);

        if(!bRetValue)
        {
            // this is a critical error
            childRouteArray = CRouteArray();
            return 0;
        }

        do
        {
            // try setting the current way point to its desired position
            bRetValue = childRouteArray.setWayPoint(wayPointArray[uiWayPointIndex], uiRouteIndex, uiWayPointIndexAtRoute);

            if(!bRetValue)
            {
                // no good, the place was already filled

                // keep the original values to avoid an infinite loop
                uiOldRouteIndex = uiRouteIndex;
                uiOldWayPointIndexAtRoute = uiWayPointIndexAtRoute;

                // go to the next position available an try again
                if(!childRouteArray.incrementIndexes(uiRouteIndex, uiWayPointIndexAtRoute))
                {
                    // out of bounds, should not happen, but is a critical error
                    childRouteArray = CRouteArray();
                    return 0;
                }
            }
        }
        while(!bRetValue || (uiOldRouteIndex == uiRouteIndex && uiOldWayPointIndexAtRoute == uiWayPointIndexAtRoute));
    }

    return 1;
}

int create_childs(const CWayPointArray& wayPointArray, CRouteArrayArray& routeArrayArray, unsigned int nChildsPerCombination = 1)
{
    unsigned int i, j, k;
    unsigned int uiOldRouteArrayArrayLength;
    // CRouteArrayArray childRoutesArrayArray;

    uiOldRouteArrayArrayLength = routeArrayArray.size();

    for(i = 0; i < uiOldRouteArrayArrayLength; i++)
    {
        for(j = i+1; j < uiOldRouteArrayArrayLength; j++)
        {
            for(k = 0; k < nChildsPerCombination; k++)
            {
                CRouteArray childRouteArray;

                if(create_child(wayPointArray, routeArrayArray[i], routeArrayArray[j], childRouteArray))
                {
                    // childRoutesArrayArray.push_back(childRouteArray);
                    routeArrayArray.push_back(childRouteArray);
                }
            }
        }
    }

    // routeArrayArray.push_back(childRoutesArrayArray);
    return 1;
}

int generate_random_route_array(const CWayPointArray& wayPointArray, CRouteArray& outRouteArray)
{
    CArrayIndexShuffler arrayIndexShuffler(wayPointArray.size());
    unsigned int uiIndex;

    for(uiIndex = 0; uiIndex < outRouteArray.size(); uiIndex++)
    {
        outRouteArray[uiIndex].setToMaxWayPoints();
    }

    while((uiIndex = arrayIndexShuffler.nextIndex()) != (unsigned int)(-1))
    {
        unsigned int uiRouteIndex;
        unsigned int uiWaypointIndex;
        // bool bPlacedOk = false;

        do
        {
            uiRouteIndex = rand() % outRouteArray.size();
            uiWaypointIndex = rand() % outRouteArray[uiRouteIndex].size();
        }
        while(!outRouteArray.setWayPoint(wayPointArray[uiIndex], uiRouteIndex, uiWaypointIndex));
    }

    return true;
}

int generate_random_route_arrays(const CWayPointArray& wayPointArray, const CRouteArray& propertySource, CRouteArrayArray& routeArrayArray,
                                    unsigned int uiNumberOfRouteArrays)
{
    unsigned int i;
    for(i = 0; i < uiNumberOfRouteArrays; i++)
    {
        CRouteArray bufferRouteArray;

        bufferRouteArray.copyPropertiesFrom(propertySource);

        if(generate_random_route_array(wayPointArray, bufferRouteArray))
        {
            routeArrayArray.push_back(bufferRouteArray);
        }
    }

    return true;
}

int average_distance(const CRouteArrayArray& routeArrayArray, double &fAverageDistance)
{
    unsigned int i;
    double fAverageBuffer = 0.0;

    for(i = 0; i < routeArrayArray.size(); i++)
    {
        fAverageBuffer += routeArrayArray[i].totalTravelDistanceCached();
    }

    fAverageDistance = fAverageBuffer / (double)routeArrayArray.size();

    return true;
}

int main(int argc, const char *argv[])
{
	unsigned int uiWinnerSamples;
	unsigned int uiMaxIterations;
	unsigned int uiMaxSamples;
	unsigned int uiIteration;
	// unsigned int uiRoutes.

	CWayPointArray wayPointArray;
	CRouteArrayArray routeArrayArray;
	CRouteArray propertySourceArray;

    // cfg values
	uiWinnerSamples = 50;
    uiMaxIterations = 100;
    uiMaxSamples = 100;

    // set route array properties
    propertySourceArray.resize(2);

    propertySourceArray[0].setMaxWayPoints(20);
    propertySourceArray[1].setMaxWayPoints(20);

    // set way point array
    unsigned int i;
    for(i = 0; i < 30; i++)
    {
        CWayPoint bufferWayPoint;

        bufferWayPoint.m_fx = rand() % 1000;
        bufferWayPoint.m_fy = rand() % 1000;
        bufferWayPoint.m_Id = i+1;

        wayPointArray.push_back(bufferWayPoint);

        fprintf(stdout, "WayPoint:%u\t%f\t%f\n", wayPointArray[i].id(), bufferWayPoint.m_fx, bufferWayPoint.m_fy);
    }

    // generate_random_route_arrays(wayPointArray, propertySourceArray, routeArrayArray, 60);

    /*
    for(i = 0; i < routeArrayArray.size(); i++)
    {
        // fprintf(stdout, "RouteArray:%u\n", i);

        unsigned int j;
        for(j = 0; j < routeArrayArray[i].size(); j++)
        {
            unsigned int k;
            for(k = 0 ; k < routeArrayArray[i][j].size(); k++)
            {
                fprintf(stdout, "RA: %4u, R: %2u, WP: %4u, %f, %f\n", i, j, routeArrayArray[i][j][k].id(), routeArrayArray[i][j][k].m_fx, routeArrayArray[i][j][k].m_fy);
            }
        }
    }
    */


    for(uiIteration = 0; uiIteration < uiMaxIterations; uiIteration++)
    {
        // generate the random routes needed
        if(routeArrayArray.size() < uiMaxSamples)
            generate_random_route_arrays(wayPointArray, propertySourceArray, routeArrayArray, uiMaxSamples - routeArrayArray.size());

        // generate childs
        create_childs(wayPointArray, routeArrayArray, 1);

        // obtain distance
    	unsigned int j;
    	for(j = 0; j < routeArrayArray.size(); j++)
    	{
            routeArrayArray[j].updateTotalTravelDistanceCache();
    	}

    	// get average distance
		double fAverageDistance;
		average_distance(routeArrayArray, fAverageDistance);
		fprintf(stdout, "BF: %u\t%f\n", uiIteration, fAverageDistance);

		// filter best samples
		filter_best_samples(routeArrayArray, uiWinnerSamples);

		// get average distance
		average_distance(routeArrayArray, fAverageDistance);
        fprintf(stdout, "AF: %u\t%f\t%f\n", uiIteration, fAverageDistance, routeArrayArray.size() > 0 ? routeArrayArray[0].totalTravelDistanceCached(): 0.0);
    }

    return 0;
}
