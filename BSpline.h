#ifndef __BCURVE_H__
#define __BCURVE_H__

#include <vector>
#include <cassert>

using std::vector;

class BSpline
{
public:
    typedef std::vector<double> stlDVec;
    BSpline(int dim = 2, int p = 0)
        :m_dimension(dim), m_degree(p)
    {}
    BSpline(int dim, int p, const double* uArr, const double* cptArr, int cptNum);

    BSpline(int dim, int p, stlDVec&& U, stlDVec&& CP)
        :m_dimension(dim), m_degree(p), m_knotVec(U), m_controlPointCoordVec(CP)
    {}
    BSpline(const BSpline&) = default;
    BSpline(BSpline&&) = default;
    BSpline& operator=(const BSpline&) = default;
    virtual ~BSpline() {}

    bool checkKnotNum() const
    {
        return m_knotVec.size() == getControlPointNum() + m_degree + 1;
    }

    // m_dimension operator
    int dimension() const { return m_dimension; }

    // m_degree operator
    int p() const { return m_degree; }
    void changeP(int d);

    virtual void clear();

    // knots operator
    const stlDVec& getKnots() const { return m_knotVec; }
    void setKnots(stlDVec && a) { m_knotVec = std::move(a); }

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

protected:

    int m_dimension;
    int m_degree;
    stlDVec m_knotVec;
    stlDVec m_controlPointCoordVec;
};


#endif // !__BCURVE_H__

