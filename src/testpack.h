#ifndef TESTPACK_H
#define TESTPACK_H
void adapt ( int ndim, double a[], double b[], int *minpts, int maxpts,
  double functn ( int indx, int ndim, double z[], double alpha[],
  double beta[] ), double rel_tol, int itest, double alpha[], double beta[],
  int lenwrk, double wrkstr[], double *relerr, double *finest, int *ifail );
double genz_function ( int indx, int ndim,
                       const double z[],
                       const double alpha[],
                       const double beta[] );
double genz_integral ( int indx, int ndim, double a[], double b[],
  double alpha[], double beta[] );
char *genz_name ( int indx );
double genz_phi ( double z );
double genz_random ( int *seed );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
int i4vec_sum ( int n, int a[] );
void multst ( int nsamp, int tstlim, int tstfns[], int tstmax, double difclt[],
  double expnts[], int ndiml, int ndims[], char *sbname,
  void subrtn ( int ndim, double a[], double b[], int *minpts, int maxpts,
    double functn ( int indx, int ndim, double z[], double alpha[],
      double beta[] ),
    double rel_tol, int itest, double alpha[], double beta[], int lenwrk,
    double wrkstr[], double *errest, double *finest, int *ifail ),
  double rel_tol, int maxpts );
double r8_abs ( double x );
double r8_epsilon ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8vec_dot ( int n, const double a1[], const double a2[] );
void r8vec_median ( int n, double r[], double rmed[3] );
double r8vec_product ( int n, double a[] );
double r8vec_sum ( int n, const double a[] );
double r8vec_mulsum ( int n, const double a[], const double b[] );
void timestamp ( void );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );
#endif
