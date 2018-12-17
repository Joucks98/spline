#define _SCL_SECURE_NO_WARNINGS
#include <Windows.h>
#include <gl/GLU.h>
#include <Eigen/Dense>
#include <iterator>
#include <numeric> // std::inner_product
#include <functional> // std::plus<>
#include "NurbsBase.h"
#include "InterpolationCurve.h"
#include "GeometryCalc.h"


#define EPS 1e-12

using namespace Eigen;
using namespace Interpolation;

double Interpolation::norm2(const double v[], int dim)
{
    if (v == nullptr)
        return -1;
    return sqrt(std::inner_product(v, v+dim, v, 0.0));
}


InterpolationCurve::InterpolationCurve()
    : BCurve(),
    readyFlag(0), isAppend(1),
    derivateIsSet(0), inFocus(0), closeState(0), m_offsetLen(0)
{}

InterpolationCurve::InterpolationCurve(int deg, int dim)
    : BCurve(dim, deg),
    readyFlag(0),
    isAppend(1),
    derivateIsSet(0),
    inFocus(0),
    closeState(0),
    m_offsetLen(0)

{

    //ptrNurbs = gluNewNurbsRenderer();//创建NURBS对象ptrNurbs
    //gluNurbsProperty(ptrNurbs, GLU_SAMPLING_TOLERANCE, 25);
    //gluNurbsProperty(ptrNurbs, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);//把表面渲染为多边形 
}

void InterpolationCurve::setDerivEnd(bool setFlag, const stlDVec* derVec)
{
    derivateIsSet = setFlag;
    if (!derivateIsSet)
        return;

    if (!m_endDerVec.empty())
        m_endDerVec.clear();
    if (derVec == nullptr)
    {
        stlDVec D0(m_dimension), DN(m_dimension);
        for (int k = 0; k < m_dimension; ++k)
        {
            D0[k] = m_interPointCoordVec[m_dimension + k] - m_interPointCoordVec[k];
            DN[k] = m_interPointCoordVec[m_interPointCoordVec.size() - m_dimension + k] -
                m_interPointCoordVec[m_interPointCoordVec.size() - 2 * m_dimension + k];
        }
        double length = chordPolyLineLength();
        if (closeState)
        {
            for (int k = 0; k < m_dimension; ++k)
            {
                D0[k] += DN[k];
                DN[k] = D0[k];
            }
            double norm_0n = norm2(&D0[0], m_dimension); // ||D0+DN||
            double lamda0n = .5*length / norm_0n;
            for (int k = 0; k < m_dimension; ++k)
            {
                D0[k] *= lamda0n;
                DN[k] *= lamda0n;
            }
        }
        else
        {
            double norm_0 = norm2(&D0[0], m_dimension);
            double norm_n = norm2(&DN[0], m_dimension);         
            double lamda0 = length / norm_0;
            double lamdan = length / norm_n;
            for (int k = 0; k < m_dimension; ++k)
            {
                D0[k] *= lamda0;
                DN[k] *= lamdan;
            }
        }

        m_endDerVec.resize(2 * m_dimension);
        std::copy(D0.begin(), D0.end(), m_endDerVec.begin());
        std::copy(DN.begin(), DN.end(), m_endDerVec.begin() + m_dimension);
    }
    else
    {
        // set m_endDerVec by derVec
        m_endDerVec.resize(derVec->size());
        std::copy(derVec->begin(), derVec->end(), m_endDerVec.begin());
    }    
}

CURVESTATE InterpolationCurve::setInterPointCoords(stlDVec && a)
{
    CURVESTATE flag = UNCHANGE;
    if (a.empty())
        return flag;
    int num = (int)a.size() / m_dimension;
    
    for (int i = 0; i < num - 1; ++i)
    {
        if (twoPointDist(&a[i*m_dimension], &a[(i + 1)*m_dimension], m_dimension) < EPS)
        {
            flag = CRVERROR;
            break;
        }
    }
  
    //m_interPointCoordVec = std::move(a);
    if (flag != CRVERROR)
    {
        m_interPointCoordVec.swap(a);
        flag = APPEND;
    }
        
    return flag;
}

CURVESTATE InterpolationCurve::addInterPointCoords(double coord[], int num)
{
    if (coord == nullptr || num != m_dimension)
        return CRVERROR;
    if (!m_interPointCoordVec.empty())
    {
        std::vector<double> tmp(m_interPointCoordVec.end() - m_dimension, m_interPointCoordVec.end());
        if (twoPointDist(&tmp[0], coord, m_dimension) < EPS)
            return CRVERROR;
    }
    m_interPointCoordVec.insert(m_interPointCoordVec.end(), &coord[0], &coord[0] + m_dimension);
    return APPEND;
}

CURVESTATE InterpolationCurve::insertInterPointCoords(double coord[], int idx)
{
    assert(idx > 0 && idx < getInterPointNum());
    if (idx == 0) //  can not add an end point
        return UNCHANGE;
    m_interPointCoordVec.insert(m_interPointCoordVec.begin() + idx*m_dimension, coord, coord + m_dimension);
    return APPEND;
}

CURVESTATE InterpolationCurve::removeInterPointCoords(int idx)
{
    assert(idx >= 0 && idx < getInterPointNum());
    m_interPointCoordVec.erase(m_interPointCoordVec.begin() + idx*m_dimension, m_interPointCoordVec.begin() + (idx + 1)*m_dimension);
    return APPEND;
}

CURVESTATE InterpolationCurve::modifyInerPointCoords(int index, const double a[])
{
    if (a == nullptr || index < 0 || index >= getInterPointNum())
        return CRVERROR;
    else
    {
        for (int i = 0; i < m_dimension; ++i)
        {            
            m_interPointCoordVec[index*m_dimension + i] = a[i];
            if (closeState && (index == 0 || index == getInterPointNum() - 1))
            {
                m_interPointCoordVec[(getInterPointNum() - 1 )*m_dimension + i] = a[i];
            }
        }
        if (index == 0 || index == 1 ||
            index == getInterPointNum() - 2 ||
            index == getInterPointNum() - 1)
            return ENDMODIFY;
        else
            return MIDMODIFY;
    }
}

void InterpolationCurve::showInterPoints(int modeType)
{
    glColor3f(0.5f, 0.4f, 0.9f);
    glPointSize(8.0);
    glEnable(GL_POINT_SMOOTH);
    if (modeType)
    {
        glPointSize(12.0);
        glDisable(GL_POINT_SMOOTH);
    }

    std::vector<double> offsetPtVec(m_interPointCoordVec);
    if (m_offsetLen != 0)
        getOffsetPt(m_offsetLen, &m_uParam[0], (int)m_uParam.size(), &offsetPtVec);
    
    glBegin(GL_POINTS);
    drawPoint(offsetPtVec, m_dimension);
    glEnd();
}

void InterpolationCurve::update(CURVESTATE flag)
{
    if (flag == CRVERROR || flag == UNCHANGE)
        return;
        
    if (flag & APPEND)
    {
        int qNum = getInterPointNum();
        if (qNum == 1) m_degree = 0;
        else if (qNum == 2) m_degree = 1;
        else if (qNum == 3) m_degree = 2;
        else m_degree = 3;
    }
    
    if (((flag & APPEND) || (flag & ENDMODIFY)))
    {
        setDerivEnd(getInterPointNum() >= 3);
    }

    NurbsBase nurbsTool;
    if (nurbsTool.generateUparameter(this) < 0 ||
        nurbsTool.generateCrvKnots(this) < 0 ||
        nurbsTool.generateCrvControlPoints(this, 0) < 0)
    {
        std::cerr << "update curve error" << std::endl;
        return;
    }
}

void InterpolationCurve::showControlPoints(int modeType)
{
    glColor3f(1.0f, 0.0f, 0.f);
    glPointSize(5.0);
    glEnable(GL_POINT_SMOOTH);

    /*std::vector<double> offsetPtVec(m_controlPointCoordVec);
    if (m_offsetLen != 0)
        getOffsetPt(m_offsetLen, &m_uParam[0], m_uParam.size(), &offsetPtVec);*/
    glBegin(GL_POINTS);
    drawPoint(m_controlPointCoordVec, m_dimension);
    glEnd();

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(3, 0x1111);
    glLineWidth(1.5);
    glBegin(GL_LINE_STRIP);
    drawPoint(m_controlPointCoordVec, m_dimension);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}

int InterpolationCurve::display(int modeType)
{
    glColor3f(1.0f, 1.0f, .5f);
    glLineWidth(2.0f);
    if (modeType)
    {
        glColor3f(0.8f, 0.8f, 0.5f);
        //glLineWidth(4.0f);
        if (inFocus)
        {
            glColor3f(1.0f, 1.0f, .5f/*1.0f, 0.5f, 0.f*/);
            glEnable(GL_LINE_STIPPLE);
            //glLineStipple(3, 0x0101);
            glLineStipple(1, 0x1111);
        }
            
    }
    
    if (m_offsetLen != 0)
    {
        std::vector<double> offsetPtVec;
        // counterclockwise, reverse sign of m_offsetLen 
        // to maintain the behavior:positive offsetlen indicates outside contraction
        ///m_offsetLen *= chordPolygonArea() > 0 ? -1 : 1;
        auto uSeries = linspace(0, 1, 1000);
        getOffsetPt((chordPolygonArea() > 0?-1:1)*m_offsetLen, &uSeries[0], (int)uSeries.size(), &offsetPtVec);
        if (inFocus)
            glBegin(GL_LINE_STRIP);
        else
            glBegin(GL_LINES);
        drawPoint(offsetPtVec, m_dimension);
        glEnd();
    }
    else
    {
        NurbsBase nurbsTool;
        nurbsTool.plotNurbs(*this);
    }
    

    if (modeType && inFocus)
    {
        glDisable(GL_LINE_STIPPLE);
    }
    return 0;
}

void InterpolationCurve::clear()
{
    BCurve::clear();
    m_controlPointCoordVec.clear();
    m_endDerVec.clear();
    readyFlag = false;
    isAppend = true;
    derivateIsSet = false;
    inFocus = false;
    m_offsetLen = 0;
}

double InterpolationCurve::chordPolyLineLength() const
{    
    double sum = 0.0;
    if(m_interPointCoordVec.size() < 2 * m_dimension)
        return 0.0;
    
    for (int i = 0; i < getInterPointNum() - 1; ++i)
    {
        sum += twoPointDist(&m_interPointCoordVec[i*m_dimension], &m_interPointCoordVec[(i + 1)*m_dimension], m_dimension);
    }
    return sum;
}

double InterpolationCurve::chordPolygonArea() const
{
    if (m_dimension != 2)
        return 0;
    double sumArea = 0.0;
    tPointd a = { m_interPointCoordVec[0], m_interPointCoordVec[1] };
    for (int i = 1; i < getInterPointNum() - 1; ++i)
    {
        tPointd b = { m_interPointCoordVec[2 * i], m_interPointCoordVec[2 * i + 1] };
        tPointd c = { m_interPointCoordVec[2 * (i + 1)], m_interPointCoordVec[2 * (i + 1) + 1] };
        sumArea += Area2(a, b, c);
    }
    return sumArea;
}

double InterpolationCurve::curveLength(double a, double b, stlDVec* polylineCoords) const
{
    assert(b > a);
    assert(a >= m_uParam[0]);
    assert(b <= m_uParam.back());

    int num = 1001;
    auto uSeries = linspace(a, b, num);
    stlDVec coords = evaluate(uSeries);
    double pre = -1.;
    double cur = polylineLength(&coords[0], m_dimension, num);
    while (abs(cur - pre) > 1e-6)
    {
        pre = cur;

        auto interUSeries = subdivide(&uSeries[0], num);
        auto interCoords = evaluate(interUSeries);
        double s1 = polylineLength(&coords[0], &interCoords[0], m_dimension, num - 1);
        double s2 = polylineLength(&interCoords[0], &coords[m_dimension], m_dimension, num - 1);
        cur = s1 + s2;


        stlDVec tmp((num << 1) - 1);
        std::merge(uSeries.begin(), uSeries.end(), interUSeries.begin(), interUSeries.end(), tmp.begin());
        uSeries.swap(tmp);
        stlDVec tmpCoords(((num << 1) - 1 )*m_dimension);
        int base = 0, k = 0;
        for (int i = 0; i < num - 1; ++i)
        {
            std::copy(&coords[base], &coords[base] + m_dimension, &tmpCoords[k]);
            std::copy(&interCoords[base], &interCoords[base] + m_dimension, &tmpCoords[k + m_dimension]);
            base += m_dimension;
            k += (m_dimension << 1);
        }
        std::copy(&coords[base], &coords[base] + m_dimension, &tmpCoords[k]);
        coords.swap(tmpCoords);
        num = (num << 1) - 1;

        //num = (num<<1) - 1; // duplicate interval number, by subdividing old intervals
        //coords = evaluate(linspace(a, b, num));
        //cur = polylineLength(&coords[0], m_dimension, num);
    }
    if (polylineCoords)
    {
        *polylineCoords = std::move(coords);
    }
    return cur;
}

//InterpolationCurve::stlDVec InterpolationCurve::evaluate(double u)
//{
//    assert(u >= 0 && u <= 1);
//    if(!getReadyFlag())
//        return stlDVec();
//
//    int idx = NurbsBase::findSpan(m_degree, u, m_knotVec);
//    stlDVec N(m_degree + 1);
//    NurbsBase::basisFuns(idx, m_degree, u, m_knotVec, &N);
//    stlDVec re(m_dimension);
//    for (int i = 0; i < m_dimension; ++i)
//    {
//        for (int j = 0; j < m_degree + 1; ++j)
//        {
//            re[i] += N[j] * m_controlPointCoordVec[2 * (idx - m_degree + j) + i];
//        }
//    }
//    return re;
//}

InterpolationCurve::stlDVec InterpolationCurve::evaluate(const stlDVec & uSeries) const
{
    stlDVec re(uSeries.size()*m_dimension);
    for (int i = 0, base = 0; i < uSeries.size(); ++i, base+=m_dimension)
    {
        auto tmp = NurbsBase::deBoor(*this, uSeries[i]);
        std::copy(tmp.begin(), tmp.end(), re.begin() + base);
    }
    return re;
}

InterpolationCurve::stlDVec InterpolationCurve::linspacePoints(int num) const
{
    assert(num > 0);
    if (num == 1)
    {
        return{ m_interPointCoordVec.end() - m_dimension, m_interPointCoordVec.end() };
    }
    stlDVec polylineCoords;
    double len = curveLength(m_uParam[0], m_uParam.back(), &polylineCoords);
    stlDVec tmp = accumulatePolylineLength(&polylineCoords[0], m_dimension, (int)polylineCoords.size() / m_dimension);
    //double step = len / (num - 1);
    //double stake = step;
    //stlDVec re(num*m_dimension);
    //std::copy(&polylineCoords[0], &polylineCoords[0] + m_dimension, &re[0]);
    //for (int i = 1, k = 1; i < tmp.size() - 1; ++i) // the head and tail points will be copied outside the loop.
    //{
    //    if (tmp[i] >= stake)
    //    {
    //        std::copy(&polylineCoords[i*m_dimension], &polylineCoords[i*m_dimension] + m_dimension, &re[(k++)*m_dimension]);
    //        stake += step;
    //    }
    //}    
    //std::copy(polylineCoords.end() - m_dimension, polylineCoords.end(), re.end() - m_dimension);
    stlDVec re(num*m_dimension);
    auto stakes = linspace(0, len, num);
    for (int k = 0; k < num-1; ++k)
    {
        auto itr = std::lower_bound(tmp.begin(), tmp.end(), stakes[k]);
        auto j = distance(tmp.begin(), itr);
        std::copy(&polylineCoords[j*m_dimension], &polylineCoords[j*m_dimension] + m_dimension, &re[k*m_dimension]);
    }
    std::copy(polylineCoords.end() - m_dimension, polylineCoords.end(), re.end() - m_dimension);
    return re;
}

void InterpolationCurve::getDerNorEndPts(double u, stlDVec* derPts, stlDVec* norPts) const
{
    assert(derPts != nullptr || norPts != nullptr);
    assert(u >= 0 && u <= 1);
    if (derPts != nullptr && derPts->size() != 2 * m_dimension)
        derPts->resize(2 * m_dimension);
    if (norPts != nullptr && norPts->size() != 2 * m_dimension)
        norPts->resize(2 * m_dimension);

    stlDVec der;
    NurbsBase nurbsTool;
    nurbsTool.curveDer_1(*this, u, 1, &der);
    double len = norm2(&der[m_dimension], m_dimension);
    if (norPts != nullptr)
    {
        stlDVec nor(der); // normal  
                          /*counterclock rotate 90*/
        nor[nor.size() - 1] *= -1;
        std::swap(nor[nor.size() - 1], nor[nor.size() - 2]);

        for (int i = 0; i < m_dimension; ++i)
        {
            nor[m_dimension + i] /= len;
            nor[m_dimension + i] += nor[i];
        }
        std::copy(nor.begin(), nor.end(), norPts->begin());
    }

    if (derPts != nullptr)
    {
        for (int i = 0; i < m_dimension; ++i)
        {
            der[m_dimension + i] /= len;
            der[m_dimension + i] += der[i];
        }
        std::copy(der.begin(), der.end(), derPts->begin());
    }
}

void InterpolationCurve::getDerNorEndPts(const double u[], int num, stlDVec * allDerPts, stlDVec * allNorPts) const
{
    assert(allDerPts != nullptr || allNorPts != nullptr);
    int n = 2 * m_dimension;
    if ((allDerPts != nullptr) && (allDerPts->size() != n*num))
        allDerPts->resize(n*num);
    if ((allNorPts != nullptr) && (allNorPts->size() != n*num))
        allNorPts->resize(n*num);
    stlDVec derPts, norPts;
    for (int i = 0; i < num; ++i)
    {
        getDerNorEndPts(u[i], &derPts, &norPts);
        if(allDerPts != nullptr) 
            std::copy(derPts.begin(), derPts.end(), allDerPts->begin()+i*n);
        if(allNorPts != nullptr) 
            std::copy(norPts.begin(), norPts.end(), allNorPts->begin()+i*n);
    }
}

void InterpolationCurve::showDerivates()
{
    if (m_controlPointCoordVec.empty())
        return;
    stlDVec allDerPts, allNorPts;
    
    stlDVec uTmp;
    int num = 100;
    double step = 1.0 / num;
    for (int i = 0; i <= num; ++i)
    {
        uTmp.push_back((i*step > 1 ? 1 : i*step));
    }
    getDerNorEndPts(&uTmp[0], num+1/*&m_uParam[0], m_uParam.size()*/, &allDerPts, &allNorPts);
    

    glBegin(GL_LINES);
    //drawPoint(allDerPts, m_dimension);
    drawPoint(allNorPts, m_dimension);
    glEnd();
}

void InterpolationCurve::showHull()
{
    if (getControlPointNum() < 3)
        return;
    GrahamConvexHull tmp(&getControlPointCoords()[0], getControlPointNum());
    stlDVec convex;
    tmp.GetConvexHull(&convex);
    glBegin(GL_LINE_LOOP);
    drawPoint(convex, m_dimension);
    glEnd();
}


static double dot(const double v1[], const double v2[], int dim)
{
    double result = 0;
    for (int i = 0; i < dim; ++i)
        result += v1[i] * v2[i];
    return result;
}

static double NewtonMethod(const InterpolationCurve& crv, double u, const double pt[])
{
    assert(u >= 0 && u <= 1);
    NurbsBase nurbsTool;
    int num = 500;

    int n = crv.dimension();

    VectorXd Q(n); 
    for (int i = 0; i < n; ++i)
        Q(i) = pt[i];

    VectorXd P(n), T(n), uT(n), dT(n), P_Q(n);
    std::vector<double> der;
    double fk_1 = 0, dfk_1 = 0;
    double uPre = 0, uCur = u;
    double rou = 0.01;
    int type = 0;
    int i = 0;
    for (;i < num; ++i)
    {
        uPre = uCur;
        nurbsTool.curveDer_1(crv, uPre, 2, &der);        
        for (int j = 0; j < n; ++j)
        {
            P(j) = der[j];
            T(j) = der[n + j];
            dT(j) = der[2 * n + j];
        }
        P_Q = P - Q;
        fk_1 = T.dot(P_Q); 

        double lamda = 1 / sqrt(T.dot(T));
        uT = T * lamda;
        double tmp = dT.dot(uT); // this not zero! :)
        dfk_1 = (dT-(dT.dot(uT))*uT).dot(P_Q) + T.dot(T);

        uCur = uPre - fk_1 / dfk_1;       
        
        double tm = sqrt(T.dot(T))*sqrt(P_Q.dot(P_Q));
        if (abs(fk_1) / tm < EPS)
        {
            type = 1;
            break;
        }
        if (abs(uCur - uPre) < EPS)
        {
            type = 2;
            break;
        }
        if (uCur < 0)
        {
            uCur = uPre - rou*fk_1 / dfk_1;
            //uCur = uCur < 0 ? 0 : uCur;
        }
        if (uCur > 1)
        {
            uCur = uPre - rou*fk_1 / dfk_1;
            //uCur = uCur > 1 ? 1 : uCur;
        }
        if (abs(uCur - u) > 0.01)
        {
            type = 3;
            uCur = u + rou*(uCur - u); 
        }
        uCur = uCur > 1 ? 1 : uCur;
        uCur = uCur < 0 ? 0 : uCur;
    }
    if (type == 0)
    {
        std::cerr << "newton error:0" << std::endl;
    }
    if (type == 3)
    {
        std::cerr << "newton error:3" << std::endl;
    }
    return uCur;
}

int InterpolationCurve::FindNearestCurvePoint(const double Q[], stlDVec * crvPt) const
{
    NurbsBase nurbsTool;
    stlDVec P;
    double minDist = LONG_MAX;
    int num = 50;
    double uStep = 1.0/num, u = 0, minu = -1;
    for (int i = 0; i <= num; ++i)
    {
        if (nurbsTool.evaluate(*this, min(u,1), &P) != 0)
        {
            continue;
            //break; // reserve for go into
        }

        double distPQ = twoPointDist(&P[0], Q, m_dimension);
        if (minDist > distPQ)
        {
            minDist = distPQ;
            minu = min(u,1);
        }            
        u += uStep;
    }

    if (minu >= 0)
    {
        minu = NewtonMethod(*this, minu, Q);
        return nurbsTool.evaluate(*this, minu, crvPt);
    }
    else
    {
        return 1;
    }
}

CURVESTATE InterpolationCurve::SetClose(bool b)
{
    CURVESTATE flag = UNCHANGE;
    if (closeState != b)
    {
        closeState = b;
        if (getInterPointNum() < 3)
            return flag;
        if (closeState)
        {
            setAppendFlag(false);
            //error: m_interPointCoordVec address changed when push_back
            //flag = addInterPointCoords(&m_interPointCoordVec[0], m_dimension); 
            stlDVec frontPt(m_interPointCoordVec.begin(), m_interPointCoordVec.begin() + m_dimension);
            flag = addInterPointCoords(&frontPt[0], m_dimension);
            return flag;
        }
        else
        {
            setAppendFlag(true);
            flag =removeInterPointCoords(getInterPointNum() - 1);
            return flag;
        }
    }  
    return flag;    
}

void InterpolationCurve::getOffsetPt(double offsetRatio, double u, stlDVec * offsetPt) const
{
    assert(offsetPt != nullptr);
    assert(u >= 0 && u <= 1);
    stlDVec der;
    NurbsBase nurbsTool;
    nurbsTool.curveDer_1(*this, u, 1, &der);    
    
    stlDVec normal(der.begin()+m_dimension, der.end());
    if (m_dimension == 2)
    {
        /*counterclock rotate 90*/
        normal[1] *= -1;
        std::swap(normal[0], normal[1]);
        double len = norm2(&normal[0], 2);
        double lambda = offsetRatio / len;
        offsetPt->assign({ der[0] + normal[0] * lambda, der[1] + normal[1] * lambda });
    }
    else
    {
        //...
    }
    
    return;
}

void InterpolationCurve::getOffsetPt(double offsetRatio, const double u[], int num, stlDVec * offsetPts) const
{
    assert(offsetPts != nullptr);
    //assert(sizeof(u) / sizeof(double) == num);
    offsetPts->clear();
    offsetPts->reserve(m_dimension*num);
    for (int i = 0; i < num; ++i)
    {
        std::vector<double> tmp;
        getOffsetPt(offsetRatio, u[i], &tmp);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(*offsetPts));
    }
}

void InterpolationCurve::setOffsetLength(double l)
{
    m_offsetLen = l;
    return;
}

BCurve InterpolationCurve::bCurve() const
{
    return BCurve(m_dimension, m_degree, &m_knotVec[0], &m_controlPointCoordVec[0], getControlPointNum());
}

void InterpolationCurve::drawPoint(stlDVec & vec, int dim)
{
    /*if (!vec.empty())
    {
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(dim, GL_FLOAT, 0, &vec[0]);
        glDrawArrays(GL_POINTS, 0, vec.size() / dim);
        glDisableClientState(GL_VERTEX_ARRAY);
    }*/

    
    for (int i = 0; i < vec.size() / dim; ++i)
    {
    
        if (2 == dim)
            glVertex2d(vec[2 * i], vec[2 * i + 1]);
        else
            glVertex3d(vec[3 * i], vec[3 * i + 1], vec[3 * i + 2]);
    }
}
