#pragma once
#ifndef __INTERPOLATIONCURVE__
#define __INTERPOLATIONCURVE__

#include <assert.h>
#include <vector>
#include <iostream>
using std::vector;

typedef enum
{
    CRVERROR = -1,
    UNCHANGE,
    APPEND,
    ENDMODIFY,
    MIDMODIFY = 4,
    CLOSED = 8
}CURVESTATE;

namespace Interpolation
{
    double norm2(const double v[], int dim);
    double twoPointDist(const double p1[], const double p2[], int dim);
}

class InterpolationCurve
{
    typedef std::vector<double> stlDVec;
public:
    InterpolationCurve(): degree(0), dimension(2), 
                          readyFlag(0), isAppend(1), 
                          derivateIsSet(0), inFocus(0), closeState(0),polylineLen(0), offsetLen(0)
    {}
    InterpolationCurve(int deg, int dim = 2) : readyFlag(0), 
                                               isAppend(1),
                                               derivateIsSet(0), 
                                               inFocus(0), 
                                               closeState(0),
                                               polylineLen(0),
                                               offsetLen(0)
                                               
    {
        degree = deg;
        dimension = dim;
        //ptrNurbs = gluNewNurbsRenderer();//创建NURBS对象ptrNurbs
        //gluNurbsProperty(ptrNurbs, GLU_SAMPLING_TOLERANCE, 25);
        //gluNurbsProperty(ptrNurbs, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);//把表面渲染为多边形 
    }
    ~InterpolationCurve()
    {
       /*if(ptrNurbs != nullptr)
       {
           gluDeleteNurbsRenderer(ptrNurbs);
            ptrNurbs = nullptr;
       }*/           
    }

    bool getReadyFlag()
    {
        //...
        /*if (!controlPointCoordVec.empty() && checkKnotNum())
            readyFlag = true;*/
        return readyFlag = !interPointCoordVec.empty() &&!controlPointCoordVec.empty() && checkKnotNum();

    }
    bool appendable()const  { return isAppend; }
    void setAppendFlag(bool a) { isAppend = a; }
    bool getDerivFlag()const { return derivateIsSet; }
    void setFocus(bool a) { inFocus = a; }
    bool isInFocus()const { return inFocus; }
    
    bool isClosed()const { return closeState; }

    const stlDVec& getEndDers() { return endDerVec; }
    void setDerivEnd(bool setFlag, const std::vector<double>* derVec = nullptr);

    // dimension operator
    /*void setDimension(int d)
    {
        assert(d > 0);
        dimension = d;
    }*/
    int getDimension() const  { return dimension;}

    // degree operator
    void setDegree(int d);
    int getDegree() const  { return degree;}    

    // uParam operator
    const stlDVec& getUparam() const { return uParam; }
    void setUparam(stlDVec&& u) { uParam = std::move(u); }

    // knots operator
    const stlDVec& getKnots() const { return knotVec; }
    void setKnots(stlDVec && a) { knotVec = std::move(a);}
    /*void setKnots(const stlDVec & a)
    {
        knotVec = a;
    }*/
    bool checkKnotNum() const
    {
        return knotVec.size() == getControlPointNum() + degree + 1;
    }
    
    bool isKnotReady() const
    {
        return knotVec.size() == getControlPointNum() + degree + 1;
    }



    // inter point operator
    const stlDVec& getInterPointCoords() const
    {
        return interPointCoordVec;
    }
    CURVESTATE setInterPointCoords(stlDVec && a);
    
    CURVESTATE addInterPointCoords(double coord[], int num); 
    CURVESTATE insertInterPointCoords(double coord[], int idx);
    CURVESTATE removeInterPointCoords(int idx);
    CURVESTATE modifyInerPointCoords(int index, const double a[]);
    
    int getInterPointNum() const
    {
        return (int)interPointCoordVec.size() / dimension;
    }
    void showInterPoints(int modeType = 0);    


    void update(CURVESTATE flag);

    // control point operator

    const stlDVec& getControlPointCoords() const
    {
        return controlPointCoordVec;
    }
    void setControlPointCoords(stlDVec && a)
    {
        controlPointCoordVec = std::move(a);
    }
    int getControlPointNum() const
    {
        return (int)controlPointCoordVec.size() / dimension;
    }
    void showControlPoints(int modeType = 0);
    int display(int modeType = 0);
    void clear();    
    double getPolyLineLen() const;
    double getPolygonArea() const;

    stlDVec evaluate(double u);
    int evaluate(int num, stlDVec* coords);

    void getDerNorEndPts(double u, stlDVec* derPts, stlDVec* norPts) const;
    void getDerNorEndPts(const double u[], int num, stlDVec* allDerPts, stlDVec* allNorPts) const;
    void showDerivates();
    void showHull();

    int FindNearestCurvePoint(const double Q[], stlDVec* crvPt) const;
    CURVESTATE SetClose(bool b);
    void getOffsetPt(double offsetRatio, double u, stlDVec* offsetPt) const;
    void getOffsetPt(double offsetRatio, const double u[], int num, stlDVec* offsetPts) const;
    void setOffsetLength(double l);
private:
    void drawPoint(stlDVec & vec, int dim);
    static bool init();

    int degree;
    int dimension;
    bool readyFlag, isAppend, derivateIsSet, inFocus, closeState;
    double polylineLen, offsetLen;
    stlDVec uParam;
    stlDVec knotVec;
    stlDVec interPointCoordVec;
    stlDVec controlPointCoordVec;
    stlDVec endDerVec;
    static stlDVec uniformVec;
    static bool _init;    
};

#endif // !__INTERPOLATIONCURVE__
