#ifndef __NURBSBASE__
#define __NURBSBASE__

#include <vector>
#include <xutility>

using std::vector;

class BSpline;
class InterpolationCurve;
class NurbsBase
{
public:
    NurbsBase();
    ~NurbsBase();
    // ����Q0....Qm��ͨ�����ķ����uk��ֵ,
    // stride dim=3��ʾxyz�������������dim=2��ʾxy;
    // num��ʾQ�ĸ���, Ҳ����uArr����ĳ��� num = n + 1��
    // return dist sum or error(-1)
    static double generateUparameter(const double* qArr, int num, int dim, double* uArr, int type = 0);
    int generateUparameter(InterpolationCurve* crv);
    


    // uArr is u parameter vector
    // len = n + 1 is the length of uArr
    // order = p + 1
    // make sure knotArr size = len + order 
    // or size = len + order + 2 if derivate ends are set.
    static int generateKnots(const double* uArr, int len, int order, double* knotArr, bool derivIsSet = false);
    int generateCrvKnots(InterpolationCurve* crv);

    static vector<double> generateKonts(const double* uArr, int uLen, int h, int order);
    

    // �˺����㷨��<<the nurbs books>>  P74
    //���룺p��ʾp������,int�ࣻm ��ʾ�ڵ����䳤��length(U)-1��Ҳ����u0-u10�ĽǱ�10��int��;U[]��ʾ�����ڵ����䣬double����;i��ʾ��i���������������i�����Ƶ��Ӧ��int���ͣ�u��ʾC(u)���Ա�������0-1��������ֵ��double���͡�
    /*
    ///OneBasicFuns����ʵ��СƬ��,��Ҫ��main������ʵ��OneBasicFuns��Ա��������Ҫ�ȶ���NURBSCLASS��һ������Ķ���basicfun����ͨ��basicfun.OneBasicFuns���к������ã�
    ///������main������ֻ�ܵ���NURBSClass�е�public��Ա������
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
    static int findSpan(int p, double u, const vector<double>& U);
    static std::pair<int, int> findSpanMult(int p, double u, const vector<double>& U);
    static int basisFuns(int idx, int p, double u, const vector<double>& knotVec, vector<double>* N);
    static int basisFuns(int idx, int p, double u, const vector<double>& knotVec, double N[]);
    // Make sure Q size is n+1
    // intput: Q, knotVec, D0, DN
    // output: P of n+3 elems
    int solveTridiagonal(const vector<double>& Q, const vector<double>& knotVec, double P1, double Pn1, vector<double>* P);

    // store matrix elements row by row
    int constructMatrix(const vector<double>& uVec, 
                        const vector<double>& knotArr, 
                        int p, vector<double>* A, bool extend = false);
    int constructMatrixN(const vector<double>& uVec,
                         const vector<double>& knotArr,
                         int p, vector<double>* N);

    int generateCrvControlPoints(InterpolationCurve* crv, bool tri = true);
    //int plotNurbs(const BSpline& crv);

    static int doolittleLU(double A[],int rowColnum);
    /*input :A is LU result matrix, B each column is b 
      output:x results are stored in B
      */
    static int luSolver(const double A[], int rowColnum, double B[], int colnum);

    /*DersBasisFun:��������ĸ���������
    output:ders��(n+1)x(p+1)�ľ���0-n�׵Ļ������ĵ�����*/
    static void dersBasisFuns(int i, double u, int p, int n, const vector<double>& U, double**ders);
    

    //==========================================================//
    /*��B�������ߵĸ��׵�����
    output:CK��һ����d+1��x3�ľ��󣬵�0�б�ʾ�����ϵĵ㣬��1�б�ʾ�ĵ��һ�׵������Դ�����*/
    static void curveDer_1(/*int n, */int p, const vector<double>& U, 
        const vector<double>& PC, double u, int d, int dim, vector<double>* CK);
    /*output: value and d derivates in parameter u*/
    bool curveDer_1(const BSpline& crv, double u, int d, vector<double>* der);
    
    double curvature(const BSpline& crv, double u);
    int evaluate(const BSpline& crv, double u, vector<double>* val);
    static vector<double> deBoor(const BSpline& crv, double u);
    static int curveKnotIns(
        const vector<double>& U, const vector<double>& CP, int dim, int p, double u, int h,
        vector<double>* UQ, vector<double>* Qw);
    // ptsCoordVec: input fitting points coordinate array with len points, 
    //              ordered by encounter order along the spline that you sought after.
    // h+1: spline control points number
    // p: spline degree
    // make sure: npt > h >= p >= 1, if npt==h+1, function return an interpolation spline.
    BSpline splineFitting(const vector<double>& ptsCoordVec, int npt, int h, int p);
    

};

#endif // !__NURBSBASE__
