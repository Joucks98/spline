#ifndef __PAINT_H__
#define __PAINT_H__

class BSpline;
class InterpolationCurve;

namespace paint
{
    void drawPoints(const double* pointCoords, int pointsNum, int dim);
    void showControlPoints(const double* pointCoords, int pointsNum, int dim);
    void showInterPoints(const double* pointCoords, int pointsNum, int dim, int modeType = 0);
    void showDerivates(const double* lineCoords, int lineNum, int dim);
    void showHull(const double* pointCoords, int pointsNum, int dim);
    int  showInterpolationCurve(const InterpolationCurve&, int modeType = 0);

    int  plotNurbs(const BSpline& crv);

}
#endif // !__PAINT_H__

