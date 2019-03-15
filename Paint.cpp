#include <Windows.h>
#include <gl/GLU.h>
#include <memory>
#include "Paint.h"
#include "BSpline.h"
#include "InterpolationCurve.h"
#include "GeometryCalc.h"

using namespace std;


void paint::drawPoints(const double * pointCoords, int pointsNum, int dim)
{
    /*if (!vec.empty())
    {
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(dim, GL_FLOAT, 0, &vec[0]);
    glDrawArrays(GL_POINTS, 0, vec.size() / dim);
    glDisableClientState(GL_VERTEX_ARRAY);
    }*/
    if (pointCoords)
    {
        for (int i = 0; i < pointsNum; ++i)
        {

            if (2 == dim)
                glVertex2d(pointCoords[2 * i], pointCoords[2 * i + 1]);
            else
                glVertex3d(pointCoords[3 * i], pointCoords[3 * i + 1], pointCoords[3 * i + 2]);
        }
    }    
}

void paint::showControlPoints(const double * pointCoords, int pointsNum, int dim)
{
    glColor3f(1.0f, 0.0f, 0.f);
    glPointSize(5.0);
    glEnable(GL_POINT_SMOOTH);

    /*std::vector<double> offsetPtVec(m_controlPointCoordVec);
    if (m_offsetLen != 0)
    getOffsetPt(m_offsetLen, &m_uParam[0], m_uParam.size(), &offsetPtVec);*/
    glBegin(GL_POINTS);
    drawPoints(pointCoords, pointsNum, dim);
    glEnd();

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(3, 0x1111);
    glLineWidth(1.5);
    glBegin(GL_LINE_STRIP);
    drawPoints(pointCoords, pointsNum, dim);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}

void paint::showInterPoints(const double * pointCoords, int pointsNum, int dim, int modeType)
{
    glColor3f(0.5f, 0.4f, 0.9f);
    glPointSize(8.0);
    glEnable(GL_POINT_SMOOTH);
    if (modeType)
    {
        glPointSize(12.0);
        glDisable(GL_POINT_SMOOTH);
    }

    //std::vector<double> offsetPtVec(m_interPointCoordVec);
    //if (m_offsetLen != 0)
    //    getOffsetPt(m_offsetLen, &m_uParam[0], (int)m_uParam.size(), &offsetPtVec);

    glBegin(GL_POINTS);
    drawPoints(pointCoords, pointsNum, dim);
    glEnd();
}

void paint::showDerivates(const double * pointPairsArr, int linesNum, int dim)
{
    glBegin(GL_LINES);
    //drawPoints(allDerPts, m_dimension);
    drawPoints(pointPairsArr, 2* linesNum, dim);
    glEnd();
}

void paint::showHull(const double * pointCoords, int pointsNum, int dim)
{
    if (pointsNum < 3)
        return;    
    glBegin(GL_LINE_LOOP);
    drawPoints(pointCoords, pointsNum, dim);
    glEnd();
}

int paint::showInterpolationCurve(const InterpolationCurve &crv, int modeType)
{
    glColor3f(1.0f, 1.0f, .5f);
    glLineWidth(2.0f);
    if (modeType)
    {
        glColor3f(0.8f, 0.8f, 0.5f);
        //glLineWidth(4.0f);
        if (crv.isInFocus())
        {
            glColor3f(1.0f, 1.0f, .5f/*1.0f, 0.5f, 0.f*/);
            glEnable(GL_LINE_STIPPLE);
            //glLineStipple(3, 0x0101);
            glLineStipple(1, 0x1111);
        }

    }

    if (crv.getOffsetLength() != 0)
    {
        std::vector<double> offsetPtVec;
        // counterclockwise, reverse sign of m_offsetLen 
        // to maintain the behavior:positive offsetlen indicates outside contraction
        ///m_offsetLen *= chordPolygonArea() > 0 ? -1 : 1;
        auto uSeries = linspace(0, 1, 1000);
        crv.getOffsetPt((crv.chordPolygonArea() > 0 ? -1 : 1)*crv.getOffsetLength(), uSeries, &offsetPtVec);
        if (crv.isInFocus())
            glBegin(GL_LINE_STRIP);
        else
            glBegin(GL_LINES);
        drawPoints(&offsetPtVec[0], (int)offsetPtVec.size() / crv.dimension(), crv.dimension());
        glEnd();
    }
    else
    {
        plotNurbs(crv);
    }


    if (modeType && crv.isInFocus())
    {
        glDisable(GL_LINE_STIPPLE);
    }
    return 0;
}

int paint::plotNurbs(const BSpline & crv)
{
    GLUnurbsObj* ptrNurbs = gluNewNurbsRenderer();//创建NURBS对象ptrNurbs
    gluNurbsProperty(ptrNurbs, GLU_SAMPLING_TOLERANCE, 25);
    gluNurbsProperty(ptrNurbs, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);//把表面渲染为多边形

    if (!crv.checkKnotNum())
        return -1;
    unique_ptr<GLfloat[]> controlPoints(new GLfloat[crv.getControlPointNum() * 3]);
    for (int k = 0; k < crv.getControlPointNum(); ++k)
    {
        controlPoints[3 * k] = (GLfloat)crv.getControlPointCoords()[crv.dimension()*k + 0];
        controlPoints[3 * k + 1] = (GLfloat)crv.getControlPointCoords()[crv.dimension()*k + 1];
        controlPoints[3 * k + 2] = 1.0f;
    }
    unique_ptr<GLfloat[]> knots(new GLfloat[crv.getKnots().size()]);
    for (size_t k = 0; k < crv.getKnots().size(); ++k)
        knots[k] = (GLfloat)crv.getKnots()[k];

    gluBeginCurve(ptrNurbs);
    gluNurbsCurve(ptrNurbs, (GLint)crv.getKnots().size(), knots.get(), 3, controlPoints.get(), crv.p() + 1, GL_MAP1_VERTEX_3);
    gluEndCurve(ptrNurbs);

    gluDeleteNurbsRenderer(ptrNurbs);
    return 0;
}
