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

class BCurve
{
public:
    typedef std::vector<double> stlDVec;
    BCurve(int dim = 2, int p = 0)
        :m_dimension(dim), m_degree(p)
    {}
    BCurve(int dim, int p, const double* uArr, const double* cptArr, int cptNum)
        :m_dimension(dim), m_degree(p)        
    {
        m_knotVec.assign(uArr, uArr + p + cptNum + 1);
        m_controlPointCoordVec.assign(cptArr, cptArr + dim * cptNum);
    }
    BCurve(int dim, int p, stlDVec&& U, stlDVec&& CP)
        :m_dimension(dim), m_degree(p), m_knotVec(U), m_controlPointCoordVec(CP)
    {}
    BCurve(const BCurve&) = default;
    BCurve& operator=(const BCurve&) = default;
    virtual ~BCurve() {}


    // m_dimension operator
    int dimension() const { return m_dimension; }
    // m_degree operator
    int p() const { return m_degree; }
    void changeP(int d)
    {
        assert(d >= 0);
        if (m_degree != d)
        {
            m_degree = d;
            clear();
        }
    }
    void clear()
    {
        m_knotVec.clear();
        m_controlPointCoordVec.clear();
    }
protected:
    int m_dimension;
    int m_degree;
    stlDVec m_knotVec;
    stlDVec m_controlPointCoordVec;
};

class InterpolationCurve : public BCurve
{
    //typedef std::vector<double> stlDVec;
public:
    InterpolationCurve(): BCurve(),
                          readyFlag(0), isAppend(1), 
                          derivateIsSet(0), inFocus(0), closeState(0), m_offsetLen(0)
    {}
    InterpolationCurve(int deg, int dim = 2) 
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

    void setDegree(int d) {
        changeP(d);
        if (d != m_degree)
        {
            readyFlag = false;
            isAppend = true;
        }
    }

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

    bool readyFlag, isAppend, derivateIsSet, inFocus, closeState;
    double m_offsetLen;
    stlDVec m_uParam;
    stlDVec m_interPointCoordVec;
    stlDVec m_endDerVec;
};

#endif // !__INTERPOLATIONCURVE_H__
