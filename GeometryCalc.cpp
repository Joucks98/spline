#define _SCL_SECURE_NO_WARNINGS
#include <numeric> // std::inner_product
#include "GeometryCalc.h"
using namespace std;


namespace
{
        
}

double Area2(const tPointd a, const tPointd b, const tPointd c)
{
    return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]);
}

void CrossProduct(const double a[3], const double b[3], double c[3])
{
    for (int i = 0; i < 3; ++i)
    {
        int i1 = (i + 1) % 3;
        int i2 = (i + 2) % 3;
        c[i] = a[i1] * b[i2] - a[i2] * b[i1];
    }
}

double GetPolygonArea(const tPointd P[], int num)
{
    if (num < 3)
        return 0;
    double sumArea = 0.0;
    const tPointd &a = P[0];
    for (int i = 1; i < num - 1; ++i)
    {
        const tPointd &b = P[i];
        const tPointd &c = P[i + 1];
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

int Compare(const void * tpi, const void * tpj, const tsPoint * ptr)
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


vector<tPointi> lineBresenham(int xBegin, int yBegin, int xEnd, int yEnd)
{
    int dx = abs(xEnd - xBegin);
    int dy = abs(yEnd - yBegin);

    vector<tPointi> output(max(dx, dy) + 1);

    int p = 2 * dy - dx;
    int two_dy = p + dx;
    int two_dy_minus_dx = p - dx;
    int x = 0, y = 0, ydir = 0;
    bool reversed = xBegin > xEnd;
    if (reversed)
    {
        ydir = yBegin - yEnd;
        x = xEnd;
        y = yEnd;
        xEnd = xBegin;
        yEnd = yBegin;
    }
    else
    {
        ydir = yEnd - yBegin;
        x = xBegin;
        y = yBegin;
    }
    int i = 0;
    //(*output).resize(dx + 1, {0,0});  // resize error, but vector<tPointi> o(m) success.
    while (x < xEnd)
    {
        output[i][0] = x;
        output[i][1] = y;

        if (p < 0)
        {
            p += two_dy;
        }
        else
        {
            p += two_dy_minus_dx;
            if (ydir > 0) ++y;
            else if (ydir < 0) --y;
        }
        ++x;
        ++i;
    }
    output[i][0] = xEnd;
    output[i][1] = yEnd;
    if (reversed)
    {
        reverse(output.begin(), output.end());
    }
    return output;
}


vector<pair<int, int>> circlePixels(int xCenter, int yCenter, int radius)
{
    int x = 0;
    int r = radius, y = r;
    int d = 3 - (r << 1);

    vector<pair<int, int>> tmp;
    while (x < y)
    {
        tmp.emplace_back(x, y);
        if (d >= 0)
        {
            y -= 1;
            d += ((x - y) << 2) + 10;
        }
        else
        {
            d += (x << 2) + 6;
        }
        x += 1;
        //DrawPoints(200, 200, x, y);
    }
    tmp.emplace_back(x, y);

    auto copyNum = static_cast<int>(tmp.size() - 1);
    vector<pair<int, int>> out(8 * copyNum);
    // {0b0001, 0b0100, 0b0110, 0b0011, 0b1011, 0b1110, 0b1100, 0b1001};
    int octaveQuadrant[] = { 1,4,6,3,11,14,12,9 };
    for (int ii = 0, index = 0; ii < 8; ++ii, index += copyNum)
    {
        bool swapXY = !(bool)(octaveQuadrant[ii] & 1);
        auto start = out.begin() + index;
        auto end = start + copyNum;
        if (ii & 1)
        {
            std::copy_n(tmp.rbegin(), copyNum, start);
        }
        else
        {
            std::copy_n(tmp.begin(), copyNum, start);
            /*std::transform(sIter, eIter, out->begin() + index, [&](auto& p) {
            auto q = p;
            if (swapXY) swap(q.first, q.second);
            if (minus1) q.first = ~q.first + 1;
            if (minus2) q.second = ~q.second + 1;
            return q;
            });*/
        }

        if (swapXY)
        {
            for (auto itr = start; itr != end; ++itr)
            {
                swap((*itr).first, (*itr).second);
            }
        }
        if (octaveQuadrant[ii] & 0b0010) // minus1
        {
            for (auto itr = start; itr != end; ++itr)
            {
                (*itr).first = ~(*itr).first + 1;
            }
        }
        if (octaveQuadrant[ii] & 0b1000) // minus2
        {
            for (auto itr = start; itr != end; ++itr)
            {
                (*itr).second = ~(*itr).second + 1;
            }
        }
    }

    if (xCenter != 0 || yCenter != 0)
    {
        for (auto& v : out)
        {
            v.first += xCenter;
            v.second += yCenter;
        }
    }

    return out;
}

double twoPointDist(const double p1[], const double p2[], int dim)
{
    if (p1 == nullptr || p2 == nullptr)
        return -1;

    auto op2 = [](double a, double b) {return (a - b)*(a - b); };
    return sqrt(std::inner_product(p1, p1 + dim, p2, 0.0, std::plus<double>(), op2));
}

vector<double> pointPairsDistance(const double * startPtCoords, const double * endPtCoords, int dim, int numSections)
{
    assert(startPtCoords != nullptr);
    assert(endPtCoords != nullptr);
    assert(numSections > 0);

    vector<double> re(numSections,0.0);
    for (int i = 0, base = 0; i < numSections; ++i, base += dim)
    {
        re[i] = twoPointDist(startPtCoords + base, endPtCoords + base, dim);
    }
    return re;

}

double polylineLength(const double * startPtCoords, const double * endPtCoords, int dim, int numSections)
{
    assert(startPtCoords != nullptr);
    assert(endPtCoords != nullptr);
    assert(numSections > 0);

    double sum = 0.0;
    for (int i = 0, base = 0; i < numSections; ++i, base += dim)
    {
        sum += twoPointDist(startPtCoords + base, endPtCoords + base, dim);
    }
    return sum;
}

double polylineLength(const double * ptCoords, int dim, int numPts)
{
    assert(ptCoords != nullptr);
    assert(numPts > 0);
    if (numPts < 2)
        return 0.0;

    double sum = 0.0;
    const double* cur = &ptCoords[0];
    for (int i = 0; i < numPts - 1; ++i, cur += dim)
    {
        sum += twoPointDist(cur, cur + dim, dim);        
    }
    return sum;
    //vector<double>tmp(pointPairsDistance(ptCoords, ptCoords + dim, dim, numPts - 1));
    //return std::accumulate(tmp.begin(), tmp.end(), 0.0);
}



vector<double> accumulatePolylineLength(const double * coords, int dim, int num)
{
    assert(coords != nullptr);
    assert(num > 0);

    if (num == 1)
        return{ 0.0 };

    //vector<double> re(num, 0.0);
    //vector<double> tmp(pointPairsDistance(coords, coords + dim, dim, num - 1));
    //std::partial_sum(tmp.begin(), tmp.end(), re.begin()+1);
    vector<double> re(num, 0.0);
    const double* cur = &coords[0];
    for (int i = 0; i < num - 1; ++i, cur += dim)
    {
        re[i + 1] = re[i] + twoPointDist(cur, cur + dim, dim);     
    }
    return re;
}

vector<double> linspace(double a, double b, int num/* = 100*/)
{
    assert(num > 0);
    assert(b > a);
    if (num == 1)
    {
        return{ b };
    }
    vector<double> series(num, 0.0);
    std::iota(series.begin(), series.end(), 0.0);
    transform(series.begin(), series.end(), series.begin(), [a, b, step = 1.0/(num - 1)](auto& it) {
        return a + (b - a) * it * step; });
    return series;
}

vector<double> midPartition(const double * uArr, int num)
{
    vector<double> tmp(num);
    std::adjacent_difference(uArr, uArr + num, tmp.begin(), [](auto&a, auto& b) {
        return (a + b)*.5;
    });
    tmp.erase(tmp.begin());
    //vector<double> re(2*num - 1);
    //std::merge(uArr, uArr + num, tmp.begin() + 1, tmp.end(), re.begin());
    return tmp;
}

vector<double> subdivision(const double* uArr, int num)
{
    auto interUSeries = midPartition(uArr, num);
    vector<double> tmp((num << 1) - 1);
    std::merge(uArr, uArr+num, interUSeries.begin(), interUSeries.end(), tmp.begin());
    return tmp;
}


GrahamConvexHull::GrahamConvexHull(const double arr[], int num)
{
    for (int i = 0; i < num; ++i)
    {
        m_list1.emplace_back(i, arr[i * 2 + 0], arr[i * 2 + 1], false
        /*tsPoint{ i, arr[i*DIM + X], arr[i*DIM + Y], false }*/);
    }
}


int GrahamConvexHull::GetConvexHull(std::vector<double>* hull)
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

void GrahamConvexHull::FindLowest()
{
    auto itr = m_list1.begin();
    for (auto &p : m_list1)
    {
        if ((p.v[1] < itr->v[1]) ||
            ((abs(p.v[1] - itr->v[1]) < 1e-6) && (p.v[0] > itr->v[0])))
            std::swap(*itr, p);
    }
}

int GrahamConvexHull::Graham()
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
