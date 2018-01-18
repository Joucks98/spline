#pragma once
#include <stdlib.h>
#include <vector>
#include <list>
#include <stack>
#include <functional>
#include <assert.h>
#include <algorithm>

typedef double tPointd[2];
typedef struct tPointStruct {
    int vnum;
    tPointd v;
    bool del;
    tPointStruct(int i, double x, double y, bool b) :vnum(i), v{x,y}, del(b)
    {}
    bool operator ==(const tPointStruct& a)
    {
        return sqrt(abs(v[0] - a.v[0])*abs(v[0] - a.v[0]) 
            + abs(v[1] - a.v[1])*abs(v[1] - a.v[1])) < 1e-6;
    }
    //operator bool() { return del; }
}tsPoint;

static double Area2(const tPointd a, const tPointd b, const tPointd c)
{
    return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]);
}

static double GetPolygonArea(const tPointd P[], int num)
{
    if (num < 3)
        return 0;
    double sumArea = 0.0;
    const tPointd &a = P[0];
    for (int i = 1; i < num - 1; ++i)
    {
        const tPointd &b = P[i];
        const tPointd &c = P[i+1];
        sumArea += Area2(a, b, c);
    }
    return sumArea;
}

bool Left(const tPointd a, const tPointd b, const tPointd c)
{
    return Area2(a, b, c) > 0;
}

bool LeftOn(const tPointd a, const tPointd b, const tPointd c)
{
    double area = Area2(a, b, c);
    return  area > 0 || abs(area) < 1e-6;
}

bool Collinear(const tPointd a, const tPointd b, const tPointd c)
{
    return  abs(Area2(a, b, c)) < 1e-6;
}

bool IntersectProp(const tPointd a, const tPointd b, const tPointd c, const tPointd d)
{
    if (Collinear(a, b, c) ||
        Collinear(a, b, d) ||
        Collinear(c, d, a) ||
        Collinear(c, d, b))
        return false;
    return (Area2(a, b, c)*Area2(a, b, d) < 0) &&
        (Area2(c, d, a)*Area2(c, d, b) < 0);
}

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

static int Compare(const void * tpi, const void * tpj, const tsPoint* ptr)
{
    tsPoint* pi = const_cast<tsPoint*>(static_cast<const tsPoint*>(tpi));
    tsPoint* pj = const_cast<tsPoint*>(static_cast<const tsPoint*>(tpj));
    double a = Area2((*ptr).v, pi->v, pj->v);
    if (abs(a) < 1e-6)
    {
        double x = abs(pi->v[0] - (*ptr).v[0]) - abs(pj->v[0] - (*ptr).v[0]);
        double y = abs(pi->v[1] - (*ptr).v[1]) - abs(pj->v[1] - (*ptr).v[1]);
        if (abs(x) < 1e-6 || abs(y) < 1e-6)
        {
            if (pi->vnum > pj->vnum)
                pj->del = true;
            else
                pi->del = true;
            return 0;
        }
        else if ((x < 0) || (y < 0))
        {
            pi->del = true;
            return -1;
        }
        else //((x > 0) || (y > 0))
        {
            pj->del = true;
            return 1;
        }

    }
    else if (a > 0)
        return -1;
    else
        return 1;

}
class GrahamConvexHull
{
public:
    
    GrahamConvexHull(const double arr[], int num)
    {
        for (int i = 0; i < num; ++i)
        {
            m_list1.emplace_back(i, arr[i*2 + 0], arr[i*2 + 1], false
            /*tsPoint{ i, arr[i*DIM + X], arr[i*DIM + Y], false }*/);
        }
    }

    int GetConvexHull(std::vector<double>* hull)
    {
        if (hull == nullptr || Graham() != 0)
            return 1;

        hull->resize(m_list1.size() * 2);

        /*for (int i = 0; i < m_list1.size(); ++i)
        {
            tmp[i*2] = m_list[i].v[0];
            tmp[i*2 + 1] = m_list[i].v[1];
        }*/
        int i = 0;
        for_each(m_list1.begin(), m_list1.end(), [hull, &i](const tsPoint& val) 
        {
            /*tmp.push_back(val.v[0]);
            tmp.push_back(val.v[1]);*/
            (*hull)[2 * i] = val.v[0];
            (*hull)[2 * i + 1] = val.v[1];
            ++i;
        });
        return 0;
    }

private:
    void FindLowest()
    {
        auto itr = m_list1.begin();
        for (auto &p :m_list1)
        {
            if ((p.v[1] < itr->v[1]) ||
                ((abs(p.v[1] - itr->v[1]) < 1e-6) && (p.v[0] > itr->v[0])))
                std::swap(*itr, p);
        }
    }


    

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

    int Graham()
    {
        if (m_list1.size() < 3) // check colinearity when size 3
            return -1;

        FindLowest();
        
        //m_list.assign(m_list1.begin(), m_list1.end());
        //Cmp::ptr = &m_list;
        //qsort(&m_list[1], m_list.size() - 1, sizeof(tsPoint), Cmp::Compare); // ?
        //m_list1.assign(m_list.begin(), m_list.end());

        //int n = Squash();
        //m_list.erase(m_list.begin() + n, m_list.end());

        tsPoint reserve = m_list1.front();
        m_list1.pop_front();
        m_list1.sort([&reserve](const tsPoint& a, const tsPoint& b) { 
            return (Compare(&a, &b, &reserve) < 0 ? true : false); 
        });
        //std::sort(m_list1.begin(), m_list1.end(), Cmp());
        m_list1.remove_if([](const tsPoint& a) {return a.del; });
        m_list1.push_front(reserve);
        int n = (int)m_list1.size();

        auto iter = m_list1.begin();
        auto preIter2 = iter;
        auto preIter1 = ++iter;
        auto curIter = ++iter;
        while (curIter != m_list1.end())
        {
            if (Area2(preIter2->v, preIter1->v, curIter->v) > 0)
            {
                ++curIter;
                ++preIter1;
                ++preIter2;
            }
            else
            {
                auto iter = m_list1.erase(preIter1);
                curIter = iter;
                preIter1 = --iter;
                preIter2 = --iter;
            }
        }
        /*std::stack<tsPoint, std::vector<tsPoint>> P;
        P.push(m_list[0]);
        P.push(m_list[1]);
        int i = 2;
        while (i < n)
        {
            const tsPoint& p1 = *(P._Get_container().end()-2);
            const tsPoint& p2 = *(P._Get_container().end()-1);
            if (Area2(p1.v, p2.v, m_list[i].v) > 0)
            {
                P.push(m_list[i]);
                ++i;
            }
            else
            {
                P.pop();
            }
        }
        std::vector<tsPoint> tmp = P._Get_container();
        m_list.swap(tmp);*/
        return 0;
    }
    //std::list<tsPoint> m_list;
    std::list<tsPoint> m_list1;
};
