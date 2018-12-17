#ifndef __NURBSBASE__
#define __NURBSBASE__

#include <vector>
#include <xutility>
#include<Windows.h>
#include<gl/GLU.h>

class BCurve;
class InterpolationCurve;
//class GLUnurbsObj;
class NurbsBase
{
public:
    NurbsBase();
    ~NurbsBase();
    // 根据Q0....Qm，通过向心法求解uk的值,
    // stride dim=3表示xyz三个坐标分量；dim=2表示xy;
    // num表示Q的个数, 也等于uArr数组的长度 num = n + 1；
    // return dist sum or error(-1)
    static double generateUparameter(const double* qArr, int num, int dim, double* uArr);
    int generateUparameter(InterpolationCurve* crv);
    


    // uArr is u parameter vector
    // len = n + 1 is the length of uArr
    // order = p + 1
    // make sure knotArr size = len + order 
    // or size = len + order + 2 if derivate ends are set.
    static int generateKnots(const double* uArr, int len, int order, double* knotArr, bool derivIsSet = false);
    int generateCrvKnots(InterpolationCurve* crv);
    

    // 此函数算法见<<the nurbs books>>  P74
    //输入：p表示p次样条,int类；m 表示节点区间长度length(U)-1，也就是u0-u10的角标10，int类;U[]表示整个节点区间，double类型;i表示第i个基函数，是与第i个控制点对应，int类型；u表示C(u)的自变量，从0-1的任意数值，double类型。
    /*
    ///OneBasicFuns函数实现小片段,若要在main函数中实现OneBasicFuns成员函数，则要先定义NURBSCLASS的一个具体的对象basicfun，并通过basicfun.OneBasicFuns进行函数调用，
    ///此外在main函数中只能调用NURBSClass中的public成员函数；
    int main()
    {
    double Nip;
    NURBS basicfun;
    int p = 1;
    int m = 10;
    double U[] = { 0,0,0,1,2,3,4,4,5,5,5 };
    int j = 3;
    double u = 2.5;
    Nip=basicfun.OneBasicFuns(p, m, U, j, u);
    cout << Nip;
    return 0;
    }
    */
    double oneBasicFuns(int p, int m, const double U[], int i, double u);
    static int findSpan(int p, double u, const std::vector<double>& U);
    static std::pair<int, int> findSpanMult(int p, double u, const std::vector<double>& U);
    static int basisFuns(int idx, int p, double u, const std::vector<double>& knotVec, std::vector<double>* N);
    // Make sure Q size is n+1
    // intput: Q, knotVec, D0, DN
    // output: P of n+3 elems
    int solveTridiagonal(const std::vector<double>& Q, const std::vector<double>& knotVec, double P1, double Pn1, std::vector<double>* P);

    int constructMatrix(const std::vector<double>& uVec, 
                        const std::vector<double>& knotArr, 
                        int p, std::vector<double>* A, bool extend = false);

    int generateCrvControlPoints(InterpolationCurve* crv, bool tri = true);
    int plotNurbs(const BCurve& crv);

    static int doolittleLU(double A[],int rowColnum);
    /*input :A is LU result matrix, B each column is b 
      output:x results are stored in B
      */
    static int luSolver(const double A[], int rowColnum, double B[], int colnum);

    /*DersBasisFun:求基函数的各级倒数，
    output:ders是(n+1)x(p+1)的矩阵，0-n阶的基函数的导数，*/
    static void dersBasisFuns(int i, double u, int p, int n, const std::vector<double>& U, double**ders);
    

    //==========================================================//
    /*求B样条曲线的各阶导数，
    output:CK是一个（d+1）x3的矩阵，第0行表示曲线上的点，第1行表示改点的一阶导数，以此类推*/
    static void curveDer_1(/*int n, */int p, const std::vector<double>& U, 
        const std::vector<double>& PC, double u, int d, int dim, std::vector<double>* CK);
    /*output: value and d derivates in parameter u*/
    void curveDer_1(const InterpolationCurve& crv, double u, int d, std::vector<double>* der);
    
    int evaluate(const InterpolationCurve& crv, double u, std::vector<double>* val);
    static std::vector<double> deBoor(const BCurve& crv, double u);
    static int curveKnotIns(
        const std::vector<double>& U, const std::vector<double>& CP, int dim, int p, double u, int h,
        std::vector<double>* UQ, std::vector<double>* Qw);
private:

    GLUnurbsObj* ptrNurbs;
};

#endif // !__NURBSBASE__
