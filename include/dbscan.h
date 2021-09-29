#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

using namespace std;


class Point{
public:
    double x,y,z;
    int clusterID;

    Point(double _x, double _y, double _z){
        this->clusterID = -1;
        this->x = _x;
        this->y = _y;
        this->z = _z;
    }

    Point(std::vector<double>& pts){
        this->clusterID = -1;
        this->x = pts[0];
        this->y = pts[1];
        this->z = pts[2];
    }

};

class DBSCAN {
public:    
    DBSCAN(unsigned int minPts, double eps, std::vector<Point> points){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = points.size();
    }
    ~DBSCAN(){}

    int run();
    std::vector<int> calculateCluster(Point point);
    int expandCluster(Point point, int clusterID);
    inline double calculateDistance(Point& pointCore, Point& pointTarget);

    int getTotalPointSize() {return m_pointSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}
    
public:
    std::vector<Point> m_points;
    
private:    
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    double m_epsilon;
};

#endif // DBSCAN_H
