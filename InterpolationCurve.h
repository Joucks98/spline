#ifndef __INTERPOLATIONCURVE_H__
#define __INTERPOLATIONCURVE_H__

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
}

class InterpolationCurve
{
    typedef std::vector<double> stlDVec;
public:
    InterpolationCurve(): m_degree(0), m_dimension(2), 
                          readyFlag(0), isAppend(1), 
                          derivateIsSet(0), inFocus(0), closeState(0), m_offsetLen(0)
    {}
    InterpolationCurve(int deg, int dim = 2) : readyFlag(0), 
                                               isAppend(1),
                                               derivateIsSet(0), 
                                               inFocus(0), 
                                               closeState(0),
                                               m_offsetLen(0)
                                               
    {
        m_degree = deg;
        m_dimension = dim;
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
        /*if (!m_controlPointCoordVec.empty() && checkKnotNum())
            readyFlag = true;*/
        return readyFlag = !m_interPointCoordVec.empty() &&!m_controlPointCoordVec.empty() && checkKnotNum();

    }
    bool appendable()const  { return isAppend; }
    void setAppendFlag(bool a) { isAppend = a; }
    bool getDerivFlag()const { return derivateIsSet; }
    void setFocus(bool a) { inFocus = a; }
    bool isInFocus()const { return inFocus; }
    
    bool isClosed()const { return closeState; }

    const stlDVec& getEndDers() { return m_endDerVec; }
    void setDerivEnd(bool setFlag, const stlDVec* derVec = nullptr);

    // m_dimension operator
    /*void setDimension(int d)
    {
        assert(d > 0);
        m_dimension = d;
    }*/
    int getDimension() const  { return m_dimension;}

    // m_degree operator
    void setDegree(int d);
    int getDegree() const  { return m_degree;}    

    // m_uParam operator
    const stlDVec& getUparam() const { return m_uParam; }
    void setUparam(stlDVec&& u) { m_uParam = std::move(u); }

    // knots operator
    const stlDVec& getKnots() const { return m_knotVec; }
    void setKnots(stlDVec && a) { m_knotVec = std::move(a);}
    /*void setKnots(const stlDVec & a)
    {
        m_knotVec = a;
    }*/
    bool checkKnotNum() const
    {
        return m_knotVec.size() == getControlPointNum() + m_degree + 1;
    }
    
    bool isKnotReady() const
    {
        return m_knotVec.size() == getControlPointNum() + m_degree + 1;
    }



    // inter point operator
    const stlDVec& getInterPointCoords() const
    {
        return m_interPointCoordVec;
    }
    CURVESTATE setInterPointCoords(stlDVec && a);
    
    CURVESTATE addInterPointCoords(double coord[], int num); 
    CURVESTATE insertInterPointCoords(double coord[], int idx);
    CURVESTATE removeInterPointCoords(int idx);
    CURVESTATE modifyInerPointCoords(int index, const double a[]);
    
    int getInterPointNum() const
    {
        return (int)m_interPointCoordVec.size() / m_dimension;
    }
    void showInterPoints(int modeType = 0);    


    void update(CURVESTATE flag);

    // control point operator

    const stlDVec& getControlPointCoords() const
    {
        return m_controlPointCoordVec;
    }
    void setControlPointCoords(stlDVec && a)
    {
        m_controlPointCoordVec = std::move(a);
    }
    int getControlPointNum() const
    {
        return (int)m_controlPointCoordVec.size() / m_dimension;
    }
    void showControlPoints(int modeType = 0);
    int display(int modeType = 0);
    void clear();    
    double chordPolyLineLength() const;
    double chordPolygonArea() const;
    double curveLength(double a = 0.0, double b = 1.0, stlDVec* polylineCoords = nullptr) const;


    //stlDVec evaluate(double u);
    stlDVec evaluate(const stlDVec& uSeries) const;
    stlDVec linspacePoints(int num) const; // num : equally spaced points number. 

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
    //static bool init(int num);

    int m_degree;
    int m_dimension;
    bool readyFlag, isAppend, derivateIsSet, inFocus, closeState;
    double m_offsetLen;
    stlDVec m_uParam;
    stlDVec m_knotVec;
    stlDVec m_interPointCoordVec;
    stlDVec m_controlPointCoordVec;
    stlDVec m_endDerVec;
    //static stlDVec uniformVec;
    //static bool _init;    
};

#endif // !__INTERPOLATIONCURVE_H__
