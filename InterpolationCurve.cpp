#include <Windows.h>
#include <gl/GLU.h>
#include "NurbsBase.h"
#include "InterpolationCurve.h"
#include "GeometryCalc.h"
#include <Eigen/Dense>
#include <iterator>
#define EPS 1e-6
using namespace Eigen;
using namespace Interpolation;

InterpolationCurve::stlDVec InterpolationCurve::uniformVec;
bool InterpolationCurve::_init = InterpolationCurve::init();

double Interpolation::norm2(const double v[], int dim)
{
    if (v == nullptr)
        return -1;
    double tmp = 0.0;
    for (int i = 0; i < dim; ++i)
    {
        tmp += v[i] * v[i];
    }
    return sqrt(tmp);
}

double Interpolation::twoPointDist(const double p1[], const double p2[], int dim)
{
    if (p1 == nullptr || p2 == nullptr)
        return -1;
    std::vector<double> d(dim);
    for (int i = 0; i < dim; ++i)
    {
        d[i] = p1[i] - p2[i];
    }
    return norm2(&d[0], dim);
}
void InterpolationCurve::setDerivEnd(bool setFlag, const std::vector<double>* derVec)
{
    derivateIsSet = setFlag;
    if (!derivateIsSet)
        return;

    if (!endDerVec.empty())
        endDerVec.clear();
    if (derVec == nullptr)
    {
        stlDVec D0(dimension), DN(dimension);
        for (int k = 0; k < dimension; ++k)
        {
            D0[k] = interPointCoordVec[dimension + k] - interPointCoordVec[k];
            DN[k] = interPointCoordVec[interPointCoordVec.size() - dimension + k] -
                interPointCoordVec[interPointCoordVec.size() - 2 * dimension + k];
        }
        double length = getPolyLineLen();
        if (closeState)
        {
            for (int k = 0; k < dimension; ++k)
            {
                D0[k] += DN[k];
                DN[k] = D0[k];
            }
            double norm_0n = norm2(&D0[0], dimension); // ||D0+DN||
            double lamda0n = .5*length / norm_0n;
            for (int k = 0; k < dimension; ++k)
            {
                D0[k] *= lamda0n;
                DN[k] *= lamda0n;
            }
        }
        else
        {
            double norm_0 = norm2(&D0[0], dimension);
            double norm_n = norm2(&DN[0], dimension);         
            double lamda0 = length / norm_0;
            double lamdan = length / norm_n;
            for (int k = 0; k < dimension; ++k)
            {
                D0[k] *= lamda0;
                DN[k] *= lamdan;
            }
        }

        endDerVec.resize(2 * dimension);
        std::copy(D0.begin(), D0.end(), endDerVec.begin());
        std::copy(DN.begin(), DN.end(), endDerVec.begin() + dimension);
    }
    else
    {
        // set endDerVec by derVec
        endDerVec.resize(derVec->size());
        std::copy(derVec->begin(), derVec->end(), endDerVec.begin());
    }    
}

void InterpolationCurve::setDegree(int d)
{
    assert(d > 0);
    degree = d;
    // clear nurbs data to restart
    controlPointCoordVec.clear();
    knotVec.clear();
    readyFlag = false;
    isAppend = true;
}

CURVESTATE InterpolationCurve::setInterPointCoords(stlDVec && a)
{
    CURVESTATE flag = UNCHANGE;
    if (a.empty())
        return flag;
    int num = (int)a.size() / dimension;
    
    for (int i = 0; i < num - 1; ++i)
    {
        if (twoPointDist(&a[i*dimension], &a[(i + 1)*dimension], dimension) < EPS)
        {
            flag = CRVERROR;
            break;
        }
    }
  
    //interPointCoordVec = std::move(a);
    if (flag != CRVERROR)
    {
        interPointCoordVec.swap(a);
        flag = APPEND;
    }
        
    return flag;
}

CURVESTATE InterpolationCurve::addInterPointCoords(double coord[], int num)
{
    if (coord == nullptr || num != dimension)
        return CRVERROR;
    if (!interPointCoordVec.empty())
    {
        std::vector<double> tmp(interPointCoordVec.end() - dimension, interPointCoordVec.end());
        if (twoPointDist(&tmp[0], coord, dimension) < EPS)
            return CRVERROR;
    }    
    for (int i = 0; i < dimension; ++i)
    {
        interPointCoordVec.push_back(coord[i]);
    }
    return APPEND;
}

CURVESTATE InterpolationCurve::insertInterPointCoords(double coord[], int idx)
{
    assert(idx > 0 && idx < getInterPointNum());
    if (idx == 0) //  can not add an end point
        return UNCHANGE;
    interPointCoordVec.insert(interPointCoordVec.begin() + idx*dimension, coord, coord + dimension);
    return APPEND;
}

CURVESTATE InterpolationCurve::removeInterPointCoords(int idx)
{
    assert(idx >= 0 && idx < getInterPointNum());
    interPointCoordVec.erase(interPointCoordVec.begin() + idx*dimension, interPointCoordVec.begin() + (idx + 1)*dimension);
    return APPEND;
}

CURVESTATE InterpolationCurve::modifyInerPointCoords(int index, const double a[])
{
    if (a == nullptr || index < 0 || index >= getInterPointNum())
        return CRVERROR;
    else
    {
        for (int i = 0; i < dimension; ++i)
        {            
            interPointCoordVec[index*dimension + i] = a[i];
            if (closeState && (index == 0 || index == getInterPointNum() - 1))
            {
                interPointCoordVec[(getInterPointNum() - 1 )*dimension + i] = a[i];
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

    std::vector<double> offsetPtVec(interPointCoordVec);
    if (offsetLen != 0)
        getOffsetPt(offsetLen, &uParam[0], (int)uParam.size(), &offsetPtVec);
    
    glBegin(GL_POINTS);
    drawPoint(offsetPtVec, dimension);
    glEnd();
}

void InterpolationCurve::update(CURVESTATE flag)
{
    if (flag == CRVERROR || flag == UNCHANGE)
        return;
        
    if (flag & APPEND)
    {
        int qNum = getInterPointNum();
        if (qNum == 1) degree = 0;
        else if (qNum == 2) degree = 1;
        else if (qNum == 3) degree = 2;
        else degree = 3;
    }
    
    if (((flag & APPEND) || (flag & ENDMODIFY)))
    {
        if (getInterPointNum() >= 3)
            setDerivEnd(true);
        else
            setDerivEnd(false);
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

    /*std::vector<double> offsetPtVec(controlPointCoordVec);
    if (offsetLen != 0)
        getOffsetPt(offsetLen, &uParam[0], uParam.size(), &offsetPtVec);*/
    glBegin(GL_POINTS);
    drawPoint(controlPointCoordVec, dimension);
    glEnd();

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(3, 0x1111);
    glLineWidth(1.5);
    glBegin(GL_LINE_STRIP);
    drawPoint(controlPointCoordVec, dimension);
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
    
    if (offsetLen != 0)
    {
        std::vector<double> offsetPtVec;
        // counterclockwise, reverse sign of offsetLen 
        // to maintain the behavior:positive offsetlen indicates outside contraction
        ///offsetLen *= getPolygonArea() > 0 ? -1 : 1;
        getOffsetPt((getPolygonArea() > 0?-1:1)*offsetLen, &uniformVec[0], (int)uniformVec.size(), &offsetPtVec);
        if (inFocus)
            glBegin(GL_LINE_STRIP);
        else
            glBegin(GL_LINES);
        drawPoint(offsetPtVec, dimension);
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
    knotVec.clear();
    interPointCoordVec.clear();
    controlPointCoordVec.clear();
    endDerVec.clear();
    readyFlag = false;
    isAppend = true;
    derivateIsSet = false;
    inFocus = false;
    polylineLen = 0;
    offsetLen = 0;
}

double InterpolationCurve::getPolyLineLen() const
{
    double sum = 0.0;
    if(interPointCoordVec.size() < 2 * dimension)
        return 0.0;
    
    for (int i = 0; i < getInterPointNum() - 1; ++i)
    {
        sum += twoPointDist(&interPointCoordVec[i*dimension], &interPointCoordVec[(i + 1)*dimension], dimension);
    }
    return sum;
}

double InterpolationCurve::getPolygonArea() const
{
    if (dimension != 2)
        return 0;
    double sumArea = 0.0;
    tPointd a = { interPointCoordVec[0], interPointCoordVec[1] };
    for (int i = 1; i < getInterPointNum() - 1; ++i)
    {
        tPointd b = { interPointCoordVec[2 * i], interPointCoordVec[2 * i + 1] };
        tPointd c = { interPointCoordVec[2 * (i + 1)], interPointCoordVec[2 * (i + 1) + 1] };
        sumArea += Area2(a, b, c);
    }
    return sumArea;
}

InterpolationCurve::stlDVec InterpolationCurve::evaluate(double u)
{
    assert(u >= 0 && u <= 1);
    if(!getReadyFlag())
        return stlDVec();

    int idx = NurbsBase::findSpan(degree, u, knotVec);
    stlDVec N(degree + 1);
    NurbsBase::basisFuns(idx, degree, u, knotVec, &N);
    stlDVec result(dimension);
    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < degree + 1; ++j)
        {
            result[i] += N[j] * controlPointCoordVec[2 * (idx - degree + j) + i];
        }
    }
    return result;
}

int InterpolationCurve::evaluate(int num, stlDVec * coords)
{
    assert(num > 0);
    assert(coords != nullptr);
    coords->resize(dimension*(num + 1));
    double step = 1.0 / num;

    double v = 0;
    for (int i = 0; i <= num; ++i, v += step)
    {
        if (v > 1) v = 1;
        stlDVec tmp(evaluate(v));
        if (tmp.empty()) return -1;
        std::copy(tmp.begin(), tmp.end(), coords->begin() + i*dimension);
    }

    return 0;
}

void InterpolationCurve::getDerNorEndPts(double u, stlDVec* derPts, stlDVec* norPts) const
{
    assert(derPts != nullptr || norPts != nullptr);
    assert(u >= 0 && u <= 1);
    if (derPts != nullptr && derPts->size() != 2 * dimension)
        derPts->resize(2 * dimension);
    if (norPts != nullptr && norPts->size() != 2 * dimension)
        norPts->resize(2 * dimension);

    stlDVec der;
    NurbsBase nurbsTool;
    nurbsTool.curveDer_1(*this, u, 1, &der);
    double len = norm2(&der[dimension], dimension);
    if (norPts != nullptr)
    {
        stlDVec nor(der); // normal  
                          /*counterclock rotate 90*/
        nor[nor.size() - 1] *= -1;
        std::swap(nor[nor.size() - 1], nor[nor.size() - 2]);

        for (int i = 0; i < dimension; ++i)
        {
            nor[dimension + i] /= len;
            nor[dimension + i] += nor[i];
        }
        std::copy(nor.begin(), nor.end(), norPts->begin());
    }

    if (derPts != nullptr)
    {
        for (int i = 0; i < dimension; ++i)
        {
            der[dimension + i] /= len;
            der[dimension + i] += der[i];
        }
        std::copy(der.begin(), der.end(), derPts->begin());
    }
}

void InterpolationCurve::getDerNorEndPts(const double u[], int num, stlDVec * allDerPts, stlDVec * allNorPts) const
{
    assert(allDerPts != nullptr || allNorPts != nullptr);
    int n = 2 * dimension;
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
    if (controlPointCoordVec.empty())
        return;
    stlDVec allDerPts, allNorPts;
    
    stlDVec uTmp;
    int num = 100;
    double step = 1.0 / num;
    for (int i = 0; i <= num; ++i)
    {
        uTmp.push_back((i*step > 1 ? 1 : i*step));
    }
    getDerNorEndPts(&uTmp[0], num+1/*&uParam[0], uParam.size()*/, &allDerPts, &allNorPts);
    

    glBegin(GL_LINES);
    //drawPoint(allDerPts, dimension);
    drawPoint(allNorPts, dimension);
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
    drawPoint(convex, dimension);
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
    //double eps = 1e-6;
    int n = crv.getDimension();

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
        if (nurbsTool.evaluate(*this, (u > 1 ? 1 : u), &P) != 0)
        {
            continue;
            //break; // reserve for go into
        }

        double distPQ = twoPointDist(&P[0], Q, dimension);
        if (minDist > distPQ)
        {
            minDist = distPQ;
            minu = (u > 1 ? 1 : u);
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
            //error: interPointCoordVec address changed when push_back
            //flag = addInterPointCoords(&interPointCoordVec[0], dimension); 
            stlDVec frontPt(interPointCoordVec.begin(), interPointCoordVec.begin() + dimension);
            flag = addInterPointCoords(&frontPt[0], dimension);
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
    
    stlDVec normal(der.begin()+dimension, der.end());
    if (dimension == 2)
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
    offsetPts->reserve(dimension*num);
    for (int i = 0; i < num; ++i)
    {
        std::vector<double> tmp;
        getOffsetPt(offsetRatio, u[i], &tmp);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(*offsetPts));
    }
}

void InterpolationCurve::setOffsetLength(double l)
{
    offsetLen = l;
    return;
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

bool InterpolationCurve::init()
{
    int num = 1000;
    double step = 1.0 / num;
    double uVal = 0;
    uniformVec.resize(num + 1);
    for (int i = 1; i < num + 1; ++i)
    {
        uniformVec[i] = uniformVec[i - 1] + step;
    }
    uniformVec[num] = 1;

    return true;
}

