#include "NurbsBase.h"
#include "EquationSolver.h"
#include "InterpolationCurve.h"
//extern bool toShowInter;
using std::cout;
using std::ostream;
using std::endl;
using std::cerr;

NurbsBase::NurbsBase()
{
    ptrNurbs = gluNewNurbsRenderer();//����NURBS����ptrNurbs
    gluNurbsProperty(ptrNurbs, GLU_SAMPLING_TOLERANCE, 25);
    gluNurbsProperty(ptrNurbs, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);//�ѱ�����ȾΪ����� 
}

NurbsBase::~NurbsBase()
{
    if (ptrNurbs != nullptr)
    {
        gluDeleteNurbsRenderer(ptrNurbs);
        ptrNurbs = nullptr;
    }

}

double NurbsBase::generateUparameter(const double * qArr, int num, int dim, double * uArr)
{
    //cout << "toShowInter address:" << &toShowInter << endl;
    if (uArr == nullptr || qArr == nullptr || num < 2)
        return -1;

    uArr[0] = 0.0;
    uArr[num - 1] = 1.0;

    double distSum = 0.0;
    std::vector<double> distArr(num - 1, 0);

    for (int i = 0; i < num - 1; ++i)
    {
        std::vector<double> curQ(qArr + dim*i, qArr + dim*(i + 1));
        std::vector<double> nexQ(qArr + dim*(i + 1), qArr + dim*(i + 2));

        for (int k = 0; k < dim; ++k)
        {
            double tmp = curQ[k] - nexQ[k];
            distArr[i] += tmp * tmp;
        }
        distArr[i] = std::sqrt(distArr[i]);
        distArr[i] = std::sqrt(distArr[i]); // d_k = sqrt(|Q_k - Q_k+1|)
        distSum += distArr[i];
    }

    for (int j = 1; j < num - 1; ++j)
    {
        uArr[j] = uArr[j - 1] + distArr[j] / distSum;
    }

    return distSum;
}

int NurbsBase::generateUparameter(InterpolationCurve* crv)
{
    int num = (*crv).getInterPointNum();
    if (num < 2)
        return -1;
    std::vector<double> uVec(num);
    generateUparameter(&(crv->getInterPointCoords()[0]), num, crv->getDimension(), &uVec[0]);
    crv->setUparam(std::move(uVec));
    return 0;
}

int NurbsBase::generateKnots(const double * uArr, int len, int order, double * knotArr, bool derivIsSet)
{
    if (uArr == nullptr || knotArr == nullptr || len - order < 0) // n >= p
        return -1;
    int knotArrSize = len + order;
    if (derivIsSet)
        knotArrSize += 2;
    for (int i = 0; i < order; ++i)
    {
        knotArr[i] = 0;
        knotArr[(knotArrSize - 1) - i] = 1;
    }

    for (int j = 0; j <= len - order + 1; ++j)
    {
        if (j == 0 || j == len - order + 1)
        {
            if (!derivIsSet)
                continue;
        }
        double tmp = 0.0;
        for (int i = j; i < j + order - 1; ++i)
        {
            tmp += uArr[i];
        }

        if (derivIsSet)
            knotArr[order - 1 + j + 1] = tmp / (order - 1);
        else
            knotArr[order - 1 + j] = tmp / (order - 1);
    }


    return 0;
}

int NurbsBase::generateCrvKnots(InterpolationCurve* crv)
{
    const std::vector<double>& uVec = (*crv).getUparam();
    if (uVec.size() < 2)
        return -1;
    std::vector<double> knotVecTmp(uVec.size() + (*crv).getDegree() + 1);
    if ((*crv).getDerivFlag())
        knotVecTmp.resize(uVec.size() + (*crv).getDegree() + 3);
    if ((*crv).getDegree() == 3 && (*crv).getDerivFlag())
    {
        std::copy(uVec.begin(), uVec.end(), knotVecTmp.begin() + (*crv).getDegree());
        knotVecTmp[knotVecTmp.size() - 1] = knotVecTmp[knotVecTmp.size() - 2] = knotVecTmp[knotVecTmp.size() - 3] = 1;
    }
    else
    {
        if (generateKnots(&uVec[0], (int)uVec.size(), (*crv).getDegree() + 1, &knotVecTmp[0], (*crv).getDerivFlag()))
            return -1;
    }    

    (*crv).setKnots(std::move(knotVecTmp));
    return 0;
}

double NurbsBase::oneBasicFuns(int p, int m, const double U[], int i, double u)
{    
    double *N = new double[p + 1];
    double saved, temp;
    int j, k;
    double Uleft, Uright;
    if (((i == 0) && (u == U[0])) || ((i == m - p - 1) && (u == U[m]))) //special case,
    {
        return 1.0;
    }
    if (u < U[i] || u >= U[i + p + 1])  // local  property
    {
        return 0;
    }

    for (j = 0; j <= p; j++) // initialize zero'th degree functions
    {
        if (u >= U[i + j] && u<U[i + j + 1])
        {
            N[j] = 1.0;  //Ni,0=1.0
        }
        else
        {
            N[j] = 0.0;
        }
    }

    //computing the coefficients by de boor's algorithm,an efficient program

    for (k = 1; k <= p; k++)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = (u - U[i])*N[0] / (U[i + k] - U[i]);
            //printf("saved:\t%.9f\n",saved);
        }
        for (j = 0; j<p - k + 1; j++)
        {
            Uleft = U[i + j + 1];
            Uright = U[i + j + k + 1];
            if (N[j + 1] == 0)
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j + 1] / (Uright - Uleft);
                N[j] = saved + (Uright - u)*temp;
                saved = (u - Uleft)*temp;
            }
        }
    }
    temp = N[0];
    delete[] N;
    return temp;
}

int NurbsBase::findSpan(int p, double u, const std::vector<double>& knotVec)
{
    int n = (int)(knotVec.size() - 1 - p - 1);
    assert(n > 0);
    assert(u >= 0 && u <= 1);
    if (abs(u - knotVec[n + 1]) < 1e-6)
        return n;
    int low = p;
    int high = n + 1;
    int mid = (low + high) / 2;
    while (u < knotVec[mid] || u >= knotVec[mid + 1])
    {
        if (u < knotVec[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

int NurbsBase::basisFuns(int idx, int p, double u, const std::vector<double>& knotVec, std::vector<double>* N)
{
    assert(N != nullptr);
    assert(idx >= 0 && idx <= knotVec.size() - 1 - p - 1); // m-p-1 = n -2
    N->resize(p + 1);
    (*N)[0] = 1;
    std::vector<double> left(p + 1, 0);
    std::vector<double> right(p + 1, 0);

    for (int j = 1; j <= p; ++j)
    {
        left[j] = u - knotVec[idx + 1 - j];
        right[j] = knotVec[idx + j] - u;
        double saved = 0;
        for (int r = 0; r < j; ++r)
        {
            double temp = (*N)[r] / (right[r + 1] + left[j - r]);
            (*N)[r] = saved + right[r + 1]*temp;
            saved = left[j - r] * temp;
        }
        (*N)[j] = saved;
    }
  
    return 0;
}

int NurbsBase::solveTridiagonal(const std::vector<double>& Q, const std::vector<double>& knotVec, double P1, double Pn1, std::vector<double>* P)
{
    assert(P != nullptr);
    int n = int(Q.size() - 1);
    (*P).resize(n + 3);

    // (*P)[0] (*P)[1] (*P)[n+1] (*P)[n+2]
    (*P)[0] = Q[0];
    (*P)[1] = P1;
    (*P)[n + 1] = Pn1;
    (*P)[n + 2] = Q[n];
    std::vector<double> R(n + 1);
    for (int i = 3; i < n; ++i)
        R[i] = Q[i - 1];
    std::vector<double> abc;
    basisFuns(4, 3, knotVec[4], knotVec, &abc);
    double den = abc[1];    
    (*P)[2] = (Q[1] - abc[0] * (*P)[1]) / den;
    std::vector<double> dd(n + 1);
    for (int i = 3; i < n; ++i)
    {
        dd[i] = abc[2] / den;
        basisFuns(i + 2, 3, knotVec[i + 2], knotVec, &abc);
        den = abc[1] - abc[0] * dd[i];
        (*P)[i] = (R[i] - abc[0] * (*P)[i - 1]) / den;
    }
    dd[n] = abc[2] / den;
    basisFuns(n + 2, 3, knotVec[n + 2], knotVec, &abc);
    den = abc[1] - abc[0] * dd[n];
    (*P)[n] = (Q[n - 1] - abc[2] * (*P)[n + 1] - abc[0] * (*P)[n - 1]) / den;
    
    for (int i = n - 1; i >= 2; --i)
        (*P)[i] = (*P)[i] - dd[i + 1] * (*P)[i + 1];

    return 0;
}

int NurbsBase::constructMatrix(const std::vector<double>& uVec, 
                               const std::vector<double>& knotArr, 
                               int p, std::vector<double>* A, bool extend)
{
    if (A == nullptr)
        return -1;
    // check knotArr
    //...

    
    int sizeA = (int)uVec.size();
    if (extend)
    {
        sizeA += 2;
    }
    (*A).resize(sizeA*sizeA, 0);

    for (int i = 0; i < sizeA; ++i)
    {
        if (extend)
        {
            if (i > 1 && i < sizeA - 2)
            {
                for (int j = 0; j < sizeA; ++j)
                    (*A)[i*sizeA + j] = oneBasicFuns(p, (int)(knotArr.size() - 1), &(knotArr[0]), j, uVec[i - 1]);
            }
            else if (i == 1)
            {
                (*A)[sizeA] = -1;
                (*A)[sizeA + 1] = 1;
            }
            else if (i == sizeA - 2)
            {
                (*A)[i*sizeA + sizeA - 2] = -1;
                (*A)[i*sizeA + sizeA - 1] = 1;
            }
            else // i == 0 or i == sizeA -1
                (*A)[i*sizeA + i] = 1;
        }
        else
        {
            for (int j = 0; j < sizeA; ++j)
                (*A)[i*sizeA + j] = oneBasicFuns(p, (int)(knotArr.size() - 1), &(knotArr[0]), j, uVec[i]);
        }

    }

    return 0;
}

namespace OUT
{
    ostream& operator <<(ostream& o, const std::vector<double>& a)
    {
        for (auto itr = a.begin(); itr != a.end(); ++itr)
        {
            o << *itr << " ";
        }
        o << endl;
        return o;
    }
}


int NurbsBase::generateCrvControlPoints(InterpolationCurve * crv, bool tri)
{
    int flag = 0;
    if (tri)
    {
        const std::vector<double>& QCoords = (*crv).getInterPointCoords();
        int dim = crv->getDimension();
        if (crv->getDegree() == 1)
        {
            std::vector<double> tmp(QCoords);
            crv->setControlPointCoords(std::move(tmp));
            return 0;
        }
        else if (crv->getDegree() == 2)
        {
            const std::vector<double>& uVec = (*crv).getUparam();
            std::vector<double> N;
            basisFuns(3, 2, uVec[1], crv->getKnots(), &N);

            const std::vector<double> &D = crv->getEndDers();
            int n = 2, p = 2, m = n + p + 3;
            double tmp1 = crv->getKnots()[p + 1] / p;
            double tmp2 = (1 - crv->getKnots()[m - p - 1]) / p;

            std::vector<double> tmp(dim * 5);
            for (int i = 0; i < dim; ++i)
            {
                /* tmp[dim + i] = (QCoords[dim + i] - QCoords[i]* (1 - uVec[1])*(1 - uVec[1]) - QCoords[2 * dim + i]*uVec[1]*uVec[1])
                / (2*uVec[1]*(1 - uVec[1]));*/
                tmp[i] = QCoords[i];
                tmp[dim + i] = QCoords[i] + D[i] * tmp1;
                tmp[4 * dim + i] = QCoords[2 * dim + i];
                tmp[3 * dim + i] = QCoords[2 * dim + i] - D[dim + i] * tmp2;
                tmp[2 * dim + i] = (QCoords[dim + i] - N[0] * tmp[3 * dim + i] - N[2] * tmp[dim + i]) / N[1];
            }
            crv->setControlPointCoords(std::move(tmp));
            return 0;
        }

        int nQ = (*crv).getInterPointNum();
        /*double length = crv->getPolyLineLen();
        double Q_10[2] = { QCoords[2] - QCoords[0], QCoords[3] - QCoords[1] };
        double Q_nn_1[2] = { QCoords[(nQ - 1)*dim] - QCoords[(nQ - 2)*dim], 
                             QCoords[(nQ - 1)*dim + 1] - QCoords[(nQ - 2)*dim + 1] };
        double len1 = Interpolation::norm2(Q_10, 2);
        double len2 = Interpolation::norm2(Q_nn_1, 2);*/

        const std::vector<double>& derivVec = (*crv).getEndDers();
        double Q_10[] = { derivVec[0], derivVec[1] };
        double Q_nn_1[] = { derivVec[2], derivVec[3] };

        std::vector<double> Q(nQ);
        std::vector<double> tmpCtlCoodVec((nQ + 2)*dim);

        const std::vector<double>& knotVec = (*crv).getKnots();
        for (int i = 0; i < dim; ++i)
        {
            for (int k = 0; k < nQ; ++k)
            {
                Q[k] = QCoords[i + k*dim];
            }
            double P1 = Q[0] + knotVec[4] / 3.0 * Q_10[i] /** length / len1*/;
            double Pn1 = Q[nQ - 1] - (1 - knotVec[nQ + 1]) / 3.0 * Q_nn_1[i] /** length / len2*/;
            std::vector<double> P;
            solveTridiagonal(Q, knotVec, P1, Pn1, &P);
            for (int k = 0; k < P.size(); ++k)
                tmpCtlCoodVec[i + k*dim] = P[k];
        }
        (*crv).setControlPointCoords(std::move(tmpCtlCoodVec));
    }
    else
    {
        //using namespace OUT;
        const std::vector<double>& knotVec = (*crv).getKnots();
        const std::vector<double>& uVec = (*crv).getUparam();
        cout << "u parameters are :";
        cout << uVec;
        cout << "knots are:";
        cout << knotVec;
        std::vector<double> A;        
        if (constructMatrix(uVec, knotVec, (*crv).getDegree(), &A, (*crv).getDerivFlag()) == -1)
        {
            cerr << "matrix construct error" << endl;
            flag = -1;
        }

        int sizeA = (int)sqrt(A.size());
        int dim = (*crv).getDimension();
        std::vector<double> B(sizeA * dim);

        std::vector<double> temp1, temp2;
        if ((*crv).getDerivFlag())
        {
            const std::vector<double>& derivVec = (*crv).getEndDers();
            size_t m = knotVec.size() - 1;
            if (m + 1 == uVec.size() + (*crv).getDegree() + 3) // check knotArr
            {
                int p = (*crv).getDegree();
                double tmp1 = knotVec[p + 1] / p;
                double tmp2 = (1 - knotVec[m - p - 1]) / p;

                for (int i = 0; i < dim; ++i)
                {
                    temp1.push_back(derivVec[i] * tmp1);
                    temp2.push_back(derivVec[derivVec.size() - dim + i] * tmp2);
                }
            }
            else
            {
                //...
            }
        }
        const std::vector<double>& interPointCoods = (*crv).getInterPointCoords();

        std::copy(interPointCoods.begin(), interPointCoods.begin() + dim, B.begin());
        if (!temp1.empty())
            std::copy(temp1.begin(), temp1.end(), B.begin() + dim);
        std::copy(interPointCoods.begin() + dim, interPointCoods.end() - dim, B.begin() + dim + temp1.size());
        if (!temp2.empty())
            std::copy(temp2.begin(), temp2.end(), B.end() - dim - temp2.size());
        std::copy(interPointCoods.end() - dim, interPointCoods.end(), B.end() - dim);


        // equation solver method1
        /*EquationSolver solver(&A[0], sizeA, sizeA, &B[0], dim);
        std::vector<double> tmpCtlCoodVec(sizeA * dim);
        if (solver.getSolution(&tmpCtlCoodVec[0]) == false)
        {
            cerr << "solver error" << endl;
            flag = -2;
        }
        (*crv).setControlPointCoords(std::move(tmpCtlCoodVec));*/


        // equation solver method2
        if (crv->getInterPointNum() >= 3)
        {
            doolittleLU(&A[0], sizeA);
            luSolver(&A[0], sizeA, &B[0], dim);
            cout << "new X :\n";
            cout << B;
        }
        (*crv).setControlPointCoords(std::move(B));

    }
    return flag;
}

int NurbsBase::plotNurbs(const InterpolationCurve & crv)
{
    if (!crv.checkKnotNum())
        return -1;
    GLfloat *controlPoints = new GLfloat[crv.getControlPointNum() * 3];
    for (int k = 0; k < crv.getControlPointNum(); ++k)
    {
        controlPoints[3 * k] = (GLfloat)crv.getControlPointCoords()[crv.getDimension()*k + 0];
        controlPoints[3 * k + 1] = (GLfloat)crv.getControlPointCoords()[crv.getDimension()*k + 1];
        controlPoints[3 * k + 2] = 1.0f;
    }
    GLfloat* knots = new GLfloat[crv.getKnots().size()];
    for (int k = 0; k < crv.getKnots().size(); ++k)
        knots[k] = (GLfloat)crv.getKnots()[k];

    gluBeginCurve(ptrNurbs);
    gluNurbsCurve(ptrNurbs, (GLint)crv.getKnots().size(), knots, 3, controlPoints, crv.getDegree() + 1, GL_MAP1_VERTEX_3);
    gluEndCurve(ptrNurbs);

    delete[] knots;
    knots = nullptr;
    delete[] controlPoints;
    controlPoints = nullptr;
    
    return 0;
}

int NurbsBase::doolittleLU(double A[], int rowColnum)
{
    if (A == nullptr)
        return -1;

    for (int t = 0; t < rowColnum; ++t)
    {
        for (int j = t; j < rowColnum; ++j)
        {
            if (j == 0) break;
            double temp = 0;
            for (int k = 0; k < t; ++k)
                temp += A[t*rowColnum + k] * A[k*rowColnum + j];
            A[t*rowColnum + j] -= temp;
        }
        for (int i = t + 1; i < rowColnum; ++i)
        {
            double temp = 0;
            for (int k = 0; k < t; ++k)
                temp += A[i*rowColnum + k] * A[k*rowColnum + t];
            A[i*rowColnum + t] -= temp;
            A[i*rowColnum + t] /= A[t*rowColnum + t];
        }
    }
    return 0;
}

int NurbsBase::luSolver(const double A[], int rowColnum, double B[], int colnum)
{
    assert(A != nullptr && B != nullptr);
    for (int i = 0; i < rowColnum; ++i)
    {
        for (int k = 0; k < colnum; ++k)
        {
            double temp = 0;
            for (int j = 0; j < i; ++j)
                temp += A[i*rowColnum + j] * B[j*colnum + k];
            B[i*colnum + k] -= temp;
        }
    }
    for (int i = rowColnum-1; i >= 0; --i)
    {
        for (int k = 0; k < colnum; ++k)
        {
            double temp = 0;
            for (int j = i+1; j < rowColnum; ++j)
                temp += A[i*rowColnum + j] * B[j*colnum + k];
            B[i*colnum + k] -= temp;
            B[i*colnum + k] /= A[i*rowColnum + i];
        }
    }
    return 0;
}

void NurbsBase::dersBasisFuns(int i, double u, int p, int n, const std::vector<double>& U, double ** ders)
{
    /*���룺i-�ڼ���������������Ƶ��Ӧ��u;p-������n-�����Ĵ�����U-�ڵ����䣬ders��(n+1)x(p+1)ά��*/
    /*�����ders��(n+1)x(p+1)*/

    double **ndu = new double*[p + 1];
    for (int i = 0; i <p + 1; i++)
        ndu[i] = new double[p + 1];
    double **a = new double*[n + 1];
    for (int i = 0; i < n + 1; i++)
        a[i] = new double[n + 1];
    /*double **ders = new double*[n + 1];
    for (int i = 0; i < p + 1; i++)
    ders[i] = new double[p + 1];*/

    int j, r, s1, s2, k, rk, pk, j1, j2;
    double d, saved, temp;
    double *left = new double[p + 1];
    double *right = new double[p + 1];

    ndu[0][0] = 1.0;
    for (j = 1; j <= p; j++)
    {
        left[j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        saved = 0.0;
        for (r = 0; r < j; r++)
        {/*down triangle*/
            ndu[j][r] = right[r + 1] + left[j - r];
            temp = ndu[r][j - 1] / ndu[j][r];
            /*up triangle*/
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    /*for (i = 0; i <= p; i++)
    {
    for (j = 0; j <= p; j++)
    {
    cout << ndu[i][j] << '\t';
    }
    cout << endl;
    }*/
    for (j = 0; j <= p; j++)/*�����������ֵ*/
        ders[0][j] = ndu[j][p];
    for (r = 0; r <= p; r++)//�Ժ����±����ѭ��
    {
        s1 = 0; s2 = 1;//�ı�����a����
        a[0][0] = 1.0;
        /*calculate the kth ders , k=0...n*/
        for (k = 1; k <= n; k++)
        {
            d = 0.0; rk = r - k; pk = p - k;
            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = p - r;
            for (j = j1; j <= j2; j++)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            j = s1; s1 = s2; s2 = j;
        }
    }
    r = p;
    for (k = 1; k <= n; k++)
    {
        for (j = 0; j <= p; j++)
            ders[k][j] *= r;
        r *= (p - k);
    }
    delete[]left;
    delete[]right;


    //double **ders = new double*[n + 1];
    //for (int i = 0; i <= p + 1; i++)
    //	ders[i] = new double[p + 1];

    for (i = 0; i <p + 1; i++)
    {
        delete[] ndu[i];
    }
    delete[] ndu;
    for (i = 0; i <n + 1; i++)
    {
        delete[] a[i];
    }
    delete[]a;
    //return ders;
}

void NurbsBase::curveDer_1(/*int n, */int p, const std::vector<double>& U, 
    const std::vector<double>& PC, double u, int d, int dim, vector<double>* CK)
{
    /*���룺n-���Ƶ�ĽǱ�,p,U,PC-���Ƶ�,u,d-d�׵���;�����CK,CK��һ��dx3�ľ���,��һ�б�ʾ�����ϵĵ㣬�ڶ��б�ʾһ�׵����������б�ʾ���׵�����xyz,������������*/
    //demo:
    //	int main() {
            //int i = 4; 
            //double u = 2.5;
            //int p = 2;
            //int d = 2;
            //vector<double> U = { 0,0,0,1,2,3,4,4,5,5,5 };
            //int m = sizeof(U) / sizeof(double);
            //int n = m - p - 1 - 1;
            //vector<double> P = { 2,3,0  , 4,5,0 , 5,6,0 , 5,6,0 , 8,9,0 , 1,2,0 , 5,4,0  , 5,6,0 };
            //int nn = (((p) > (d)) ? (p) : (d)) + 1;
            //for (auto& m : U) m *= .2;
            //u *= .2;
            //vector<double> CK;
            //NurbsBase::curveDer_1(p, U, P, u, d, &CK);
            //for (i = 0; i <= d; i++)
            //    cout << CK[3 * i] << '\t' << CK[3 * i + 1] << '\t' << CK[3 * i + 2] << endl;
            //return 0;
    //}

    double **ders = new double*[d + 1];
    for (int i = 0; i < d + 1; ++i)
        ders[i] = new double[p + 1];

    int du = d < p ? d : p;    
    int span = findSpan(p, u, U);
    dersBasisFuns(span, u, p, du, U, ders);

    CK->clear();
    CK->resize(dim * (d+1), 0);
    for (int k = 0; k <= du; k++)
    {
        for (int j = 0; j <= p; j++)
        {
            for (int i = 0; i < dim; ++i)
            {
                (*CK)[dim * k + i] += ders[k][j] * PC[dim*(span - p + j)+i];
            }            
        }
    }
    
    for (int i = 0; i <d + 1; ++i)
        delete[] ders[i];
    delete[]ders;

}

void NurbsBase::curveDer_1(const InterpolationCurve & crv, double u, int d, std::vector<double>* der)
{
    assert(u >= 0 && u <= 1);
    assert(der != nullptr);
    if (!const_cast<InterpolationCurve&>(crv).getReadyFlag())
    {
        der = nullptr;
        return;
    }       
    std::vector<double> derTmp;
    curveDer_1(crv.getDegree(), crv.getKnots(), crv.getControlPointCoords(), u, d, crv.getDimension(), &derTmp);
    der->swap(derTmp);
}

int NurbsBase::evaluate(const InterpolationCurve & crv, double u, std::vector<double>* val)
{
    assert(u >= 0 && u <= 1);
    assert(val != nullptr);
    if (!const_cast<InterpolationCurve&>(crv).getReadyFlag())
    {
        return 1;
    }
    val->clear();
    int p = crv.getDegree();
    int idx = findSpan(p, u, crv.getKnots());
    std::vector<double> N(p + 1);
    basisFuns(idx, p, u, crv.getKnots(), &N);
    std::vector<double> result(crv.getDimension());
    for (int i = 0; i < crv.getDimension(); ++i)
    {
        for (int j = 0; j < p + 1; ++j)
        {
            result[i] += N[j] * crv.getControlPointCoords()[2 * (idx - p + j) + i];
        }
    }
    val->swap(result);
    return 0;
}

