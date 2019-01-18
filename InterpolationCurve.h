#ifndef __INTERPOLATIONCURVE_H__
#define __INTERPOLATIONCURVE_H__

#include <assert.h>
#include <vector>
#include <iostream>
#include "BSpline.h"

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


class InterpolationCurve : public BSpline
{
public:
    InterpolationCurve();
    InterpolationCurve(int deg, int dim = 2);
    
    ~InterpolationCurve() {}

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

    // m_uParam operator
    const stlDVec& getUparam() const { return m_uParam; }
    void setUparam(stlDVec&& u) { m_uParam = std::move(u); }

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

    void update(CURVESTATE flag);

    //int display(int modeType = 0);
    virtual void clear() override;    
    double chordPolyLineLength() const;
    double chordPolygonArea() const;
    double curveLength(double a = 0.0, double b = 1.0, stlDVec* polylineCoords = nullptr) const;


    //stlDVec evaluate(double u);
    stlDVec evaluate(const stlDVec& uSeries) const;
    stlDVec linspacePoints(int num) const; // num : equally spaced points number. 

    void getDerNorEndPts(double u, stlDVec* derPts, stlDVec* norPts) const;
    void getDerNorEndPts(const double u[], int num, stlDVec* allDerPts, stlDVec* allNorPts) const;


    int FindNearestCurvePoint(const double Q[], stlDVec* crvPt = nullptr, double* u = nullptr) const;
    CURVESTATE SetClose(bool b);
    void getOffsetPt(double offsetRatio, double u, stlDVec* offsetPt) const;
    void getOffsetPt(double offsetRatio, const double u[], int num, stlDVec* offsetPts) const;
    
    void setOffsetLength(double l) { m_offsetLen = l; }
    double getOffsetLength() const { return m_offsetLen; };

    BSpline bSpline() const;
private:

    bool readyFlag, isAppend, derivateIsSet, inFocus, closeState;
    double m_offsetLen;
    stlDVec m_uParam;
    stlDVec m_interPointCoordVec;
    stlDVec m_endDerVec;
};

#endif // !__INTERPOLATIONCURVE_H__
