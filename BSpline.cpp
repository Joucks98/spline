#include "BSpline.h"

BSpline::BSpline(int dim, int p, const double * uArr, const double * cptArr, int cptNum)
    :m_dimension(dim), m_degree(p)
{
    m_knotVec.assign(uArr, uArr + p + cptNum + 1);
    m_controlPointCoordVec.assign(cptArr, cptArr + dim * cptNum);
}

void BSpline::changeP(int d)
{
    assert(d >= 0);
    if (m_degree != d)
    {
        m_degree = d;
        clear();
    }
}

void BSpline::clear()
{
    m_knotVec.clear();
    m_controlPointCoordVec.clear();
}

