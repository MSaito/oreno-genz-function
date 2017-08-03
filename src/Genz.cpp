#include "Genz.hpp"
#include <cmath>

using namespace std;

namespace {
    //*
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int *RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
    int tuple_next( int m1, int m2, int n, int rank, int x[])
    {
        if ( m2 < m1 ) {
            return 0;
        }
        if ( rank <= 0 ) {
            for (int i = 0; i < n; i++) {
                x[i] = m1;
            }
            return 1;
        }
        rank = rank + 1;
        int i = n - 1;

        for ( ; ; ) {
            if ( x[i] < m2 ) {
                x[i] = x[i] + 1;
                break;
            }
            x[i] = m1;
            if ( i == 0 ) {
                rank = 0;
                for (int j = 0; j < n; j++ ) {
                    x[j] = m1;
                }
                break;
            }
            i = i - 1;
        }
        return rank;
    }

    //*************************************************************************
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
    int i4vec_sum ( int n, int a[] )
    {
        int sum = 0;
        for (int i = 0; i < n; i++ ) {
            sum = sum + a[i];
        }
        return sum;
    }

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
    inline double r8_min ( double x, double y )
    {
        if ( y < x ) {
            return y;
        } else {
            return x;
        }
    }

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_PHI estimates the normal cumulative density function.
//
//  Discussion:
//
//    The approximation is accurate to 1.0E-07.
//
//    This routine is based upon algorithm 5666 for the error function,
//    from Hart et al.
//
//  Modified:
//
//    20 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Miller
//    C++ version by John Burkardt
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double Z, a value which can be regarded as the distance,
//    in standard deviations, from the mean.
//
//    Output, double GENZ_PHI, the integral of the normal PDF from negative
//    infinity to Z.
//
//  Local parameters:
//
//    Local, double ROOTPI, despite the name, is actually the
//    square root of TWO * pi.
//
    double genz_phi ( double z )
    {
        double expntl;
        double p;
        const double p0 = 220.2068679123761;
        const double p1 = 221.2135961699311;
        const double p2 = 112.0792914978709;
        const double p3 = 33.91286607838300;
        const double p4 = 6.373962203531650;
        const double p5 = 0.7003830644436881;
        const double p6 = 0.03526249659989109;
        const double q0 = 440.4137358247522;
        const double q1 = 793.8265125199484;
        const double q2 = 637.3336333788311;
        const double q3 = 296.5642487796737;
        const double q4 = 86.78073220294608;
        const double q5 = 16.06417757920695;
        const double q6 = 1.755667163182642;
        const double q7 = 0.08838834764831844;
        const double rootpi = 2.506628274631001;
        double zabs;

        zabs = fabs ( z );
//
//  12 < |Z|.
//
        if ( 12.0 < zabs )
        {
            p = 0.0;
        }
        else
        {
//
//  |Z| <= 12
//
            expntl = exp ( - zabs * zabs / 2.0 );
//
//  |Z| < 7
//
            if ( zabs < 7.0 )
            {
                p = expntl * ((((((
                                      p6
                                      * zabs + p5 )
                                  * zabs + p4 )
                                 * zabs + p3 )
                                * zabs + p2 )
                               * zabs + p1 )
                              * zabs + p0 ) / (((((((
                                                        q7
                                                        * zabs + q6 )
                                                    * zabs + q5 )
                                                   * zabs + q4 )
                                                  * zabs + q3 )
                                                 * zabs + q2 )
                                                * zabs + q1 )
                                               * zabs + q0 );
            }
//
//  CUTOFF <= |Z|
//
            else
            {
                p = expntl / (
                    zabs + 1.0 / (
                        zabs + 2.0 / (
                            zabs + 3.0 / (
                                zabs + 4.0 / (
                                    zabs + 0.65 ))))) / rootpi;
            }
        }

        if ( 0.0 < z )
        {
            p = 1.0 - p;
        }

        return p;
    }

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
    double r8vec_dot ( int n, const double a1[], const double a2[] )
    {
        double value = 0.0;
        for (int i = 0; i < n; i++ ) {
            value = value + a1[i] * a2[i];
        }
        return value;
    }

    void set01(int dim, double a[], double b[])
    {
        for (int i = 0; i < dim; i++) {
            a[i] = 0.0;
            b[i] = 1.0;
        }
    }
}

namespace GenzNS {
    double Oscillatory::operator()(int ndim, const double z[],
                                   const double *,
                                   const double beta[]) const
    {
        double s = 0.0;
        for (int i = 0; i < ndim; i++) {
            s += z[i];
        }
        double total = 2.0 * M_PI * beta[0] + s;
        return cos(total);
    }

    double Oscillatory::integral(int ndim, const double *,
                                 const double *,
                                 const double alpha[],
                                 const double beta[])
        const
    {
        double value = 0.0;
        int rank = 0;
        int ic[ndim];
        for (;;) {
            rank = tuple_next(0, 1, ndim, rank, ic);
            if ( rank == 0 ) {
                break;
            }
            double total = 2.0 * M_PI * beta[0];
            for (int j = 0; j < ndim; j++ ) {
                if ( ic[j] != 1 ) {
                    total = total + alpha[j];
                }
            }
            int isum = i4vec_sum ( ndim, ic );
            double s = 1 + 2 * ( ( isum / 2 ) * 2 - isum );
            if ( ( ndim % 2 ) == 0 ) {
                value = value + s * cos ( total );
            } else {
                value = value + s * sin ( total );
            }
        }
        if ( 1 < ( ndim % 4 ) ) {
            value = - value;
        }
        return value;
    }

    string Oscillatory::name() const
    {
        return "Oscillatory";
    }

    double Oscillatory::difficulty() const
    {
        return 0.9;
    }

    void Oscillatory::setTestParams(int dim,  double a[], double b[],
                                    double alpha[], double *) const
    {
        set01(dim, a, b);
        for (int i = 0; i < dim; i++) {
            b[i] = alpha[i];
        }
    }

// Product Peak
    double ProductPeak::operator()(int ndim, const double z[],
                                   const double alpha[],
                                   const double beta[]) const
    {
        double total = 1.0;
        for (int j = 0; j < ndim; j++ ) {
            total = total * (
                1.0 / pow( alpha[j], 2) + pow ( z[j] - beta[j], 2 ) );
        }
        return 1.0 / total;
    }

    double ProductPeak::integral(int ndim, const double *,
                                 const double *,
                                 const double alpha[],
                                 const double beta[]) const
    {
        double value = 1.0;
        for (int j = 0; j < ndim; j++ ) {
            value = value * alpha[j] * (
                atan ( ( 1.0 - beta[j] ) * alpha[j] )
                + atan (       + beta[j]   * alpha[j] )
                );
        }
        return value;
    }

    string ProductPeak::name() const
    {
        return "Product Peak";
    }

    double ProductPeak::difficulty() const
    {
        return 0.725;
    }
    void ProductPeak::setTestParams(int dim,  double a[], double b[],
                                    double *, double *) const
    {
        set01(dim, a, b);
    }

// Corner Peak
    double CornerPeak::operator()(int ndim, const double z[],
                                  const double alpha[],
                                  const double beta[]) const
    {
        double total = 1.0;
        for (int j = 0; j < ndim; j++ ) {
            if ( beta[j] < 0.5 ) {
                total = total + z[j];
            } else {
                total = total + alpha[j] - z[j];
            }
        }
        return 1.0 / pow ( total, ndim + 1 );
    }

    double CornerPeak::integral(int ndim, const double *,
                                const double *,
                                const double alpha[],
                                const double *) const
    {
        double value = 0.0;
        double sgndm = 1.0;
        for (int j = 1; j <= ndim; j++ ) {
            sgndm = - sgndm / static_cast<double>( j );
        }

        int rank = 0;
        int ic[ndim];
        for ( ; ; ) {
            rank = tuple_next ( 0, 1, ndim, rank, ic );
            if ( rank == 0 ) {
                break;
            }
            double total = 1.0;
            for (int j = 0; j < ndim; j++ ) {
                if ( ic[j] != 1 ) {
                    total = total + alpha[j];
                }
            }
            int isum = i4vec_sum ( ndim, ic );
            double s = 1 + 2 * ( ( isum / 2 ) * 2 - isum );
            value = value + static_cast<double>(s) / total;
        }
        return value * sgndm;
    }

    string CornerPeak::name() const
    {
        return "Corner Peak";
    }

    double CornerPeak::difficulty() const
    {
        return 0.185;
    }

    void CornerPeak::setTestParams(int dim,  double a[], double b[],
                                   double alpha[], double *) const
    {
        set01(dim, a, b);
        for (int i = 0; i < dim; i++) {
            b[i] = alpha[i];
        }
    }

// Gaussian
    double Gaussian::operator()(int ndim, const double z[],
                                const double alpha[],
                                const double beta[]) const
    {
        double total = 0.0;
        for (int j = 0; j < ndim; j++ ) {
            total = total + pow ( alpha[j] * ( z[j] - beta[j] ), 2 );
        }
        total = r8_min( total, 100.0 );
        return exp ( - total );
    }

    double Gaussian::integral(int ndim, const double *,
                              const double *,
                              const double alpha[],
                              const double beta[]) const
    {
        double value = 1.0;
        double ab = sqrt ( 2.0 );
        for (int j = 0; j < ndim; j++ ) {
            value = value * ( sqrt ( M_PI ) / alpha[j] ) *
                (   genz_phi ( ( 1.0 - beta[j] ) * ab * alpha[j] )
                    - genz_phi (       - beta[j]   * ab * alpha[j] ) );
        }
        return value;
    }

    string Gaussian::name() const
    {
        return "Gaussian";
    }

    double Gaussian::difficulty() const
    {
        return 0.703;
    }
    void Gaussian::setTestParams(int dim,  double a[], double b[],
                                 double *, double *) const
    {
        set01(dim, a, b);
    }

// C0 Function
    double C0Function::operator()(int ndim, const double z[],
                                  const double alpha[],
                                  const double beta[]) const
    {
        double total = 0.0;
        for (int j = 0; j < ndim; j++ ) {
            total = total + alpha[j] * fabs ( z[j] - beta[j] );
        }
        return exp ( - total );
    }

    double C0Function::integral(int ndim, const double *,
                                const double *,
                                const double alpha[],
                                const double beta[]) const
    {
        double value = 1.0;
        for (int j = 0; j < ndim; j++ ) {
            double ab = alpha[j] * beta[j];
            value = value *
                ( 2.0 - exp ( - ab ) - exp ( ab - alpha[j] ) ) / alpha[j];
        }
        return value;
    }

    string C0Function::name() const
    {
        return "C0 Function";
    }

    double C0Function::difficulty() const
    {
        return 2.04;
    }

    void C0Function::setTestParams(int dim, double a[], double b[],
                                    double *, double *) const
    {
        set01(dim, a, b);
    }

// Discontinuous
    double Discontinuous::operator()(int ndim, const double z[],
                                     const double alpha[],
                                     const double beta[]) const
    {
        bool test = false;

        for (int j = 0; j < ndim; j++ ) {
            if ( beta[j] < z[j] ) {
                test = true;
                break;
            }
        }
        double value = 0.0;
        if ( test ) {
            value = 0.0;
        } else {
            double total = r8vec_dot ( ndim, alpha, z );
            value = exp ( total );
        }
        return value;
    }

    double Discontinuous::integral(int ndim, const double *, const double *,
                                   const double alpha[],
                                   const double beta[]) const
    {
        double value = 1.0;
        for (int j = 0; j < ndim; j++ ) {
            value = value * ( exp ( alpha[j] * beta[j] ) - 1.0 ) / alpha[j];
        }
        return value;
    }

    string Discontinuous::name() const
    {
        return "Discontinuous";
    }

    double Discontinuous::difficulty() const
    {
        return 0.43;
    }

    void Discontinuous::setTestParams(int dim, double a[], double b[],
                                      double *, double beta[]) const
    {
        set01(dim, a, b);
        for (int i = 0; i < dim; i++) {
            beta[i] = 2 / M_PI;
        }
    }
}
