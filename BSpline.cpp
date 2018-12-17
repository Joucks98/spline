#include <Windows.h>
#include <gl/GLU.h>
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

void BSpline::showControlPoints(int modeType)
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

void BSpline::drawPoint(stlDVec & vec, int dim)
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
