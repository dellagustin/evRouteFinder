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
// #define _LOGGING_

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

    // sub structures
public:
    class CCompleteIndex
    {
    public:
        CCompleteIndex();

        // Query methods
    public:
        unsigned int routeIndex() const;
        unsigned int wayPointIndex() const;

        // operators
    public:
        CCompleteIndex& operator=(const CCompleteIndex&);
        bool operator==(const CCompleteIndex&) const;

        // members
    public:
        unsigned int m_uiRouteIndex;
        unsigned int m_uiWayPointIndex;
    };

    // Query Methods
public:
    //! Verify if this route array is compatible with route array
    //! otherRouteArray
    bool isCompatibleWith(const CRouteArray& otherRouteArray) const;
	double linearLength() const;
	double totalTravelDistanceCached() const;
	bool findWayPoint(CWayPoint::Id id, CCompleteIndex& wayPointCompleteIndex) const;
	bool findWayPoint(CWayPoint::Id id, unsigned int &uiRouteIndex, unsigned int &uiWayPointIndex) const;
	bool indexValidity(const CCompleteIndex& wayPointCompleteIndex) const;
	bool indexValidity(unsigned int uiRouteIndex, unsigned int uiWayPointIndex) const;
	const CWayPoint& wayPoint(const CCompleteIndex& wayPointIndex) const;

	// Tool Methods
public:
    bool randomIndex(CCompleteIndex&) const;
    bool incrementIndexes(CCompleteIndex&) const;
    bool incrementIndexes(unsigned int &uiRouteIndex, unsigned int &uiWayPointIndexAtRoute) const;


    // Edit Methods
public:
    bool setWayPoint(const CWayPoint& wayPoint, const CCompleteIndex &wayPointCompleteIndex, bool bOverride = false);
    bool setWayPoint(const CWayPoint& wayPoint, unsigned int uiRouteIndex, unsigned int uiWayPointIndex, bool bOverride = false);
    bool swapWayPoints(const CCompleteIndex&, const CCompleteIndex&);
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

// CCompleteIndex {{
CRouteArray::CCompleteIndex::CCompleteIndex()
{
    m_uiRouteIndex = 0;
    m_uiWayPointIndex = 0;
}

unsigned int CRouteArray::CCompleteIndex::routeIndex() const
{
    return m_uiRouteIndex;
}

unsigned int CRouteArray::CCompleteIndex::wayPointIndex() const
{
    return m_uiWayPointIndex;
}

CRouteArray::CCompleteIndex& CRouteArray::CCompleteIndex::operator=(const CCompleteIndex& otherIndex)
{
    m_uiRouteIndex = otherIndex.m_uiRouteIndex;
    m_uiWayPointIndex = otherIndex.m_uiWayPointIndex;
    return *this;
}

bool CRouteArray::CCompleteIndex::operator== (const CCompleteIndex& otherIndex) const
{
    return m_uiRouteIndex == otherIndex.m_uiRouteIndex && m_uiWayPointIndex == otherIndex.m_uiWayPointIndex;
}

// }} CCompleteIndex

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

bool CRouteArray::randomIndex(CCompleteIndex& completeIndex) const
{
    if(!size())
        return false;

    completeIndex.m_uiRouteIndex = rand() % size();
    completeIndex.m_uiWayPointIndex = rand() % (*this)[completeIndex.m_uiRouteIndex].size();

    return true;
}

bool CRouteArray::incrementIndexes(CCompleteIndex& wayPointCompleteIndex) const
{
    return incrementIndexes(wayPointCompleteIndex.m_uiRouteIndex, wayPointCompleteIndex.m_uiWayPointIndex);
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

        if (uiRouteIndex >= size())
        {
            uiRouteIndex = 0;
        }
    }

    return true;
}

bool CRouteArray::findWayPoint(CWayPoint::Id id, CCompleteIndex& wayPointCompleteIndex) const
{
    return findWayPoint(id, wayPointCompleteIndex.m_uiRouteIndex, wayPointCompleteIndex.m_uiWayPointIndex);
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

bool CRouteArray::indexValidity(const CCompleteIndex& wayPointCompleteIndex) const
{
    return indexValidity(wayPointCompleteIndex.routeIndex(), wayPointCompleteIndex.wayPointIndex());
}

bool CRouteArray::indexValidity(unsigned int uiRouteIndex, unsigned int uiWayPointIndex) const
{
    if(uiRouteIndex >= size())
        return false;

    if(uiWayPointIndex >= (*this)[uiRouteIndex].size())
        return false;

    return true;
}

const CWayPoint& CRouteArray::wayPoint(const CCompleteIndex& wayPointCompleteIndex) const
{
    if(!indexValidity(wayPointCompleteIndex))
        return *((CWayPoint*)NULL);

    return (*this)[wayPointCompleteIndex.routeIndex()][wayPointCompleteIndex.wayPointIndex()];
}

bool CRouteArray::setWayPoint(const CWayPoint& wayPoint, const CCompleteIndex &wayPointCompleteIndex, bool bOverride)
{
    return setWayPoint(wayPoint, wayPointCompleteIndex.routeIndex(), wayPointCompleteIndex.wayPointIndex(), bOverride);
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

bool CRouteArray::swapWayPoints(const CCompleteIndex& wayPointCompleteIndex1, const CCompleteIndex& wayPointCompleteIndex2)
{
    if(!indexValidity(wayPointCompleteIndex1))
        return false;

    if(!indexValidity(wayPointCompleteIndex2))
        return false;

    ::swap((*this)[wayPointCompleteIndex1.routeIndex()][wayPointCompleteIndex1.wayPointIndex()],
        (*this)[wayPointCompleteIndex2.routeIndex()][wayPointCompleteIndex2.wayPointIndex()]);

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

// CRouteArrayPtrContainer {{

class CRouteArrayPtrContainer
{
public:
    CRouteArrayPtrContainer();
    CRouteArrayPtrContainer(CRouteArray* pRouteArray);
    ~CRouteArrayPtrContainer();

public:
    void setRouteArrayPtr(CRouteArray* pRouteArray, bool bDeletePtr = true);
    CRouteArray* routeArrayPtr();
    const CRouteArray* routeArrayPtr() const;
    void deletePtr();

public:
    bool operator< (const CRouteArrayPtrContainer&) const;

private:
    CRouteArray* m_pRouteArray;
};

typedef vector<CRouteArrayPtrContainer> CRouteArrayPtrContainerArray;

CRouteArrayPtrContainer::CRouteArrayPtrContainer()
{
    m_pRouteArray = NULL;
}

CRouteArrayPtrContainer::CRouteArrayPtrContainer(CRouteArray* pRouteArray)
{
    m_pRouteArray = pRouteArray;
}

CRouteArrayPtrContainer::~CRouteArrayPtrContainer()
{
}

void CRouteArrayPtrContainer::setRouteArrayPtr(CRouteArray* pRouteArray, bool bDeletePtr)
{
    if(bDeletePtr)
    {
        deletePtr();
    }

    m_pRouteArray = pRouteArray;
}

CRouteArray* CRouteArrayPtrContainer::routeArrayPtr()
{
    return m_pRouteArray;
}

const CRouteArray* CRouteArrayPtrContainer::routeArrayPtr() const
{
    return m_pRouteArray;
}

void CRouteArrayPtrContainer::deletePtr()
{
    if(m_pRouteArray)
        delete m_pRouteArray;

    m_pRouteArray = NULL;
}

bool CRouteArrayPtrContainer::operator< (const CRouteArrayPtrContainer& otherRouteArrayPtrContainer) const
{
    if(m_pRouteArray && otherRouteArrayPtrContainer.m_pRouteArray)
    {
        return *(m_pRouteArray) < *(otherRouteArrayPtrContainer.m_pRouteArray);
    }
    else if(otherRouteArrayPtrContainer.m_pRouteArray)
    {
        return true;
    }

    return false;
}

// }} CRouteArrayPtrContainer

struct struct_sim_settings
{
    unsigned int uiWinnerSamples;
	unsigned int uiMaxIterations;
	unsigned int uiMaxSamples;
	unsigned int uiChildsPerPair;
	unsigned int uiMaxMutationLevel;
	double fMutationChance;

	// output handlers
	FILE *pStdErr;
	FILE *pStdOut;
};

void filter_best_samples(CRouteArrayArray& routeArrayArray, unsigned int uiSamples)
{
    CRouteArrayArray otherRouteArrayArray;
    CRouteArrayPtrContainerArray routeArrayPtrContainerArray;
    unsigned int i;

    routeArrayPtrContainerArray.reserve(routeArrayArray.size());

    for(i = 0; i < routeArrayArray.size(); i++)
    {
        CRouteArrayPtrContainer routeArrayPtrContainer(&routeArrayArray[i]);
        routeArrayPtrContainerArray.push_back(routeArrayPtrContainer);
    }

    sort(routeArrayPtrContainerArray.begin(), routeArrayPtrContainerArray.end());

    // remove repeated samples

    otherRouteArrayArray.reserve(uiSamples);

    if(routeArrayPtrContainerArray.size() > 0)
    {
        otherRouteArrayArray.push_back(*(routeArrayPtrContainerArray[0].routeArrayPtr()));
    }

    for(i = 1; i < routeArrayPtrContainerArray.size() && otherRouteArrayArray.size() < uiSamples ; i++)
    {
        if(routeArrayPtrContainerArray[i].routeArrayPtr()->totalTravelDistanceCached() !=
            routeArrayPtrContainerArray[i-1].routeArrayPtr()->totalTravelDistanceCached())
        {
            otherRouteArrayArray.push_back(*(routeArrayPtrContainerArray[i].routeArrayPtr()));
        }
        else
        {
#ifdef _LOGGING_
            fprintf(stderr, "Twin detected.\n");
#endif
        }
    }

    routeArrayArray = otherRouteArrayArray;
}

int create_child(const struct_sim_settings& sim_settings, const CWayPointArray& wayPointArray, const CRouteArray& routeArray1, const CRouteArray& routeArray2, CRouteArray& childRouteArray)
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
        CRouteArray::CCompleteIndex wayPointCompleteIndex;
        CRouteArray::CCompleteIndex oldWayPointCompleteIndex;
        const CRouteArray *pCurParentRouteAr = rand() & 0x01 ? &routeArray2 : &routeArray1;

        bRetValue = pCurParentRouteAr->findWayPoint(wayPointArray[uiWayPointIndex].id(), wayPointCompleteIndex);

        if(!bRetValue)
        {
            // this is a critical error
            childRouteArray = CRouteArray();
            return false;
        }

        // keep the original values to avoid an infinite loop
        oldWayPointCompleteIndex = wayPointCompleteIndex;

        do
        {
            // try setting the current way point to its desired position
            bRetValue = childRouteArray.setWayPoint(wayPointArray[uiWayPointIndex], wayPointCompleteIndex);

            if(!bRetValue)
            {
                // no good, the place was already filled

                // go to the next position available an try again
                if(!childRouteArray.incrementIndexes(wayPointCompleteIndex))
                {
                    // out of bounds, should not happen, but is a critical error
                    childRouteArray = CRouteArray();
                    return 0;
                }
            }
        }
        while(!bRetValue && !(oldWayPointCompleteIndex == wayPointCompleteIndex));
    }

    // introduces mutations
    if(sim_settings.uiMaxMutationLevel)
    {
        if(((rand() % 10000) / 10000) < sim_settings.fMutationChance)
        {
            unsigned int uiMutationLevel = sim_settings.uiMaxMutationLevel == 1 ? 1 : 1 + (rand() % (sim_settings.uiMaxMutationLevel-1));

            unsigned int i;
            for(i = 0; i < uiMutationLevel; i++)
            {
                CRouteArray::CCompleteIndex wayPointCompleteIndex1;
                CRouteArray::CCompleteIndex wayPointCompleteIndex2;
                CRouteArray::CCompleteIndex oldWayPointCompleteIndex;

                childRouteArray.randomIndex(wayPointCompleteIndex1);
                oldWayPointCompleteIndex = wayPointCompleteIndex1;

                // ensure that the first wayPoint is not empty
                while(!childRouteArray.wayPoint(wayPointCompleteIndex1).id())
                {
                    if(!childRouteArray.incrementIndexes(wayPointCompleteIndex1))
                    {
                        // critical error!
                        childRouteArray = CRouteArray();
                        return false;
                    }

                    // avoid infinite loop
                    if(wayPointCompleteIndex1 == oldWayPointCompleteIndex)
                    {
                        // critical error!
                        childRouteArray = CRouteArray();
                        return false;
                    }
                }

                do
                {
                    childRouteArray.randomIndex(wayPointCompleteIndex2);
                }
                while(wayPointCompleteIndex1 == wayPointCompleteIndex2);

                childRouteArray.swapWayPoints(wayPointCompleteIndex1, wayPointCompleteIndex2);
            }
        }
    }

    return 1;
}

int create_childs(const struct_sim_settings& sim_settings, const CWayPointArray& wayPointArray, CRouteArrayArray& routeArrayArray)
{
    unsigned int i, j, k;
    unsigned int uiOldRouteArrayArrayLength;
    // CRouteArrayArray childRoutesArrayArray;

    uiOldRouteArrayArrayLength = routeArrayArray.size();

    routeArrayArray.reserve(sim_settings.uiChildsPerPair*routeArrayArray.size()*routeArrayArray.size()/2);

    for(i = 0; i < uiOldRouteArrayArrayLength; i++)
    {
        for(j = i+1; j < uiOldRouteArrayArrayLength; j++)
        {
            for(k = 0; k < sim_settings.uiChildsPerPair; k++)
            {
                CRouteArray childRouteArray;

                if(create_child(sim_settings, wayPointArray, routeArrayArray[i], routeArrayArray[j], childRouteArray))
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
        CRouteArray::CCompleteIndex wayPointCompleteIndex;
        // bool bPlacedOk = false;

        do
        {
            outRouteArray.randomIndex(wayPointCompleteIndex);
        }
        while(!outRouteArray.setWayPoint(wayPointArray[uiIndex], wayPointCompleteIndex));
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
	unsigned int uiIteration;
	struct_sim_settings sim_settings;

	CWayPointArray wayPointArray;
	CRouteArrayArray routeArrayArray;
	// CRouteArrayPtrContainerArray routeArrayPtrContainerArray;
	CRouteArray propertySourceArray;

    // cfg values
	sim_settings.uiWinnerSamples = 80;
    sim_settings.uiMaxIterations = 1200;
    sim_settings.uiMaxSamples = 100;
    sim_settings.uiChildsPerPair = 2;
    sim_settings.uiMaxMutationLevel = 5;
    sim_settings.fMutationChance = 0.3;

    // output handlers
    sim_settings.pStdErr = stderr;
    sim_settings.pStdOut = stdout;

    // set route array properties
    // propertySourceArray.resize(2);
    propertySourceArray.resize(1);

    propertySourceArray[0].setMaxWayPoints(100);
    // propertySourceArray[1].setMaxWayPoints(20);
    // propertySourceArray[2].setMaxWayPoints(40);
    // propertySourceArray[3].setMaxWayPoints(20);
    // propertySourceArray[4].setMaxWayPoints(20);
    // propertySourceArray[5].setMaxWayPoints(20);

    // set way point array
    unsigned int i;
    for(i = 0; i < 100; i++)
    {
        CWayPoint bufferWayPoint;

        // bufferWayPoint.m_fx = rand() % 1000;
        bufferWayPoint.m_fx = 0.0;
        // bufferWayPoint.m_fy = rand() % 1000;
        bufferWayPoint.m_fy = (double)i;
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

    for(uiIteration = 0; uiIteration < sim_settings.uiMaxIterations; uiIteration++)
    {
        fprintf(sim_settings.pStdErr, "%4u: Starting Iteration.\n", uiIteration);

        // generate the random routes needed
        if(routeArrayArray.size() < sim_settings.uiMaxSamples)
        {
            fprintf(sim_settings.pStdErr, "%4u: Generating %u random routes... ", uiIteration, (unsigned int)(sim_settings.uiMaxSamples - routeArrayArray.size()));
            generate_random_route_arrays(wayPointArray, propertySourceArray, routeArrayArray, sim_settings.uiMaxSamples - routeArrayArray.size());
            fprintf(sim_settings.pStdErr, "Done!\n");
        }

        // generate childs
        fprintf(sim_settings.pStdErr, "%4u: Generating childs... ", uiIteration);
        create_childs(sim_settings, wayPointArray, routeArrayArray);
        fprintf(sim_settings.pStdErr, "Done!\n");

        // obtain distance
        fprintf(sim_settings.pStdErr, "%4u: Updating distance cache... ", uiIteration);

    	unsigned int j;
    	for(j = 0; j < routeArrayArray.size(); j++)
    	{
            routeArrayArray[j].updateTotalTravelDistanceCache();
    	}

    	fprintf(sim_settings.pStdErr, "Done!\n");

    	// get average distance
		double fAverageDistance;
		// average_distance(routeArrayArray, fAverageDistance);
		// fprintf(stdout, "%4u: BF: %f\n", uiIteration, fAverageDistance);

		// filter best samples
		fprintf(sim_settings.pStdErr, "%4u: Filtering best samples... ", uiIteration);
		filter_best_samples(routeArrayArray, sim_settings.uiWinnerSamples);
		fprintf(sim_settings.pStdErr, "Done!\n");

		// get average distance
		average_distance(routeArrayArray, fAverageDistance);
        // fprintf(stdout, "%4u AF: %f\t%f\n", uiIteration, fAverageDistance, routeArrayArray.size() > 0 ? routeArrayArray[0].totalTravelDistanceCached() : 0.0);
        // fprintf(stdout, "%4u AF: %f\t%f\n", uiIteration, fAverageDistance, routeArrayArray.size() > 0 ? routeArrayArray[0].totalTravelDistanceCached() : 0.0);
        fprintf(stdout, "%4u\t%f\t%f\n", uiIteration, fAverageDistance, routeArrayArray.size() > 0 ? routeArrayArray[0].totalTravelDistanceCached() : 0.0);

        if(routeArrayArray[0].totalTravelDistanceCached() < 30.0)
            break;
    }

#if 0
    // print the best route
    for(i = 0; i < routeArrayArray[0].size(); i++)
    {
        fprintf(stdout, "Route: %u\n", i);

        unsigned int j;
        for(j = 0; j < routeArrayArray[0][i].size(); j++)
        {
            if(routeArrayArray[0][i][j].id())
            {
                fprintf(stdout, "WayPoint - index: %u, id: %u, pos: %f, %f\n", j, routeArrayArray[0][i][j].id(), routeArrayArray[0][i][j].m_fx, routeArrayArray[0][i][j].m_fy);
            }
        }
    }
#endif

    return 0;
}
