#ifndef __BCURVE_H__
#define __BCURVE_H__

#include <vector>
#include <cassert>

using std::vector;

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
    BCurve(BCurve&&) = default;
    BCurve& operator=(const BCurve&) = default;
    virtual ~BCurve() {}

    bool checkKnotNum() const
    {
        return m_knotVec.size() == getControlPointNum() + m_degree + 1;
    }

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
    virtual void clear()
    {
        m_knotVec.clear();
        m_controlPointCoordVec.clear();
    }

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

