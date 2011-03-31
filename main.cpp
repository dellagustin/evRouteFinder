#include <iostream>
#include "../reuse/ge_points.h"

using namespace std;

class CWayPoint : public CGePoint2d
{
private:
	bool m_bIsFixed;
};

class CRoute // : vector<CWayPoint>
{
	// members
private:
	int m_nWayPoints;
};

int main(int argc, const char *argv[])
{
    cout << "Hello world!" << endl;
    return 0;
}
