#ifndef __EQUATIONSOLVER__
#define __EQUATIONSOLVER__

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace Eigen;

class EquationSolver
{
public:
    EquationSolver(const double* aElems, int rows, int cols,
        const double* bElems, int xCols) : A(rows, cols), B(rows, xCols)
    {
        /*for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
                A(i, j) = aElems[i*cols + j];
            for (int k = 0; k < xCols; ++k)
                B(i, k) = bElems[i*xCols + k];
        }*/
        int num = rows*cols;
        MatrixXd ATmp(cols, rows), BTmp(xCols, rows);
        for (int i = 0; i < num; ++i)
            ATmp(i) = aElems[i];
        A = ATmp.transpose();
        num = rows*xCols;
        for (int i = 0; i < num; ++i)
            BTmp(i) = bElems[i];
        B = BTmp.transpose();
    }
    EquationSolver(MatrixXd& a, MatrixXd& b) : A(a), B(b) {}
    ~EquationSolver() {}
    bool getSolution(double* xArr)
    {
        if (xArr == nullptr)
            return false;
        MatrixXd tmp(A.lu().solve(B));
        for (int i = 0; i < B.rows(); ++i)
            for (int j = 0; j < B.cols(); ++j)
                xArr[i*B.cols() + j] = tmp(i, j);
        std::cout << "A=\n " << A << std::endl;
        std::cout << "X=\n" << tmp << std::endl;
        std::cout << "B=\n" << B << std::endl;
        std::cout << "A*X=\n" << A*tmp << std::endl;
        return true;
    }
private:
    MatrixXd A, B;
};
#endif // !__EQUATIONSOLVER__
