#ifndef __GEOMETRYCALC_H__
#define __GEOMETRYCALC_H__

#include <stdlib.h>
#include <vector>
#include <list>
#include <stack>
#include <functional>
#include <assert.h>
#include <algorithm>
#include <xutility>
#include <array>

using std::pair;
using std::vector;
using std::list;
using std::array;

typedef array<int,2> tPointi;
typedef double tPointd[2];
typedef struct tPointStruct {
    int vnum;
    tPointd v;
    bool del;
    tPointStruct(int i, double x, double y, bool b) 
        :vnum(i), v{ x,y }, del(b)
    {}
    bool operator ==(const tPointStruct& a)
    {
        return sqrt(abs(v[0] - a.v[0])*abs(v[0] - a.v[0])
            + abs(v[1] - a.v[1])*abs(v[1] - a.v[1])) < 1e-12;
    }
    //operator bool() { return del; }
}tsPoint;

namespace gc
{
    class Point3d
    {
    public:
        Point3d(double a = 0.0, double b = 0.0, double c = 0.0)
            :m_data{ a,b,c }
        {}
        Point3d(const double p[3])
            :m_data{ p[0], p[1], p[2] }
        {}
        Point3d(const Point3d&) = default;
        Point3d& operator =(const Point3d&) = default;
        double x() const{
            return m_data[0];
        }
        double y() const{
            return m_data[1];
        }
        double z() const{
            return m_data[2];
        }
        const double & operator[](int id) const{
            assert(!(id < 0 || id > 2));
            return m_data[id];
        }
        double & operator[](int id){
            assert(!(id < 0 || id > 2));
            return m_data[id];
        }
        Point3d operator-(const Point3d& r) const
        {
            return Point3d(m_data[0] - r.m_data[0],
                m_data[1] - r.m_data[1],
                m_data[2] - r.m_data[2]);
        }
        double dot(const Point3d& r) const;
        Point3d crs(const Point3d& r) const;
        double len() const;
        operator const double*() const{
            return m_data;
        }
    private:
        double m_data[3];
    };

}


double Area2(const tPointd a, const tPointd b, const tPointd c);

void CrossProduct(const double a[3], const double b[3], double c[3]);

double GetPolygonArea(const tPointd P[], int num);

bool Left(const tPointd a, const tPointd b, const tPointd c);

bool LeftOn(const tPointd a, const tPointd b, const tPointd c);

bool Collinear(const tPointd a, const tPointd b, const tPointd c);

bool IntersectProp(const tPointd a, const tPointd b, const tPointd c, const tPointd d);

//class Cmp
//{
//public:
//    static std::vector<tsPoint>* ptr;
//    static int Compare(const void * tpi, const void * tpj)
//    {
//        assert(ptr->size() > 0);
//        tsPoint* pi = const_cast<tsPoint*>(static_cast<const tsPoint*>(tpi));
//        tsPoint* pj = const_cast<tsPoint*>(static_cast<const tsPoint*>(tpj));
//        double a = Area2((*ptr)[0].v, pi->v, pj->v);
//        if(abs(a) < 1e-6)
//        {
//            double x = abs(pi->v[0] - (*ptr)[0].v[0]) - abs(pj->v[0] - (*ptr)[0].v[0]);
//            double y = abs(pi->v[1] - (*ptr)[0].v[1]) - abs(pj->v[1] - (*ptr)[0].v[1]);
//            if (abs(x) < 1e-6 || abs(y) < 1e-6)
//            {
//                if (pi->vnum > pj->vnum)
//                    pj->del = true;
//                else
//                    pi->del = true;
//                return 0;
//            }
//            else if ((x < 0) || (y < 0))
//            {
//                pi->del = true;
//                return -1;
//            }
//            else //((x > 0) || (y > 0))
//            {
//                pj->del = true;
//                return 1;
//            }
//
//        }
//        else if (a > 0)
//            return -1;
//        else
//            return 1;
//
//    }
//
//};
//std::vector<tsPoint>* Cmp::ptr;

int Compare(const void * tpi, const void * tpj, const tsPoint* ptr);

class GrahamConvexHull
{
public:
    GrahamConvexHull(const double arr[], int num);

    int GetConvexHull(std::vector<double>* hull);

private:
    void FindLowest();
    /*int Squash()
    {
    int i = 1, j = 1;
    while (i < m_list.size())
    {
    if (!m_list[i].del)
    {
    if (i != j)
    {
    m_list[j] = m_list[i];
    }
    ++j;
    }
    ++i;
    }
    return j;
    }*/

    int Graham();

    //std::list<tsPoint> m_list;
    std::list<tsPoint> m_list1;
};


vector<tPointi> lineBresenham(int xBegin, int yBegin, int xEnd, int yEnd);
vector<pair<int, int>> circlePixels(int xCenter, int yCenter, int radius);

double twoPointDist(const double p1[], const double p2[], int dim);

list<double> pointPairsDistance(const double * startPtCoords, const double * endPtCoords, int dim, int numSections);

double polylineLength(const list<gc::Point3d>&);
double polylineLength(const double* ptCoords, int dim, int numPts);
list<double> accumulatePolylineLength(const list<gc::Point3d>&);
list<double> accumulatePolylineLength(const double* coords, int dim, int num);
list<double> linspace(double a, double b, int num = 100); // evenly distribute num values in [a, b] interval
list<double> midPartition(const list<double>& iUList);
void subdivision(list<double>* ioUList);


#endif // !__GEOMETRYCALC_H__
