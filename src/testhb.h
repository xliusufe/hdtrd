#ifndef RBS_H_INCLUDED
#define RBS_H_INCLUDED

#define EPS_PARAM 1E-6
#define EPS_KERNEL 1E-6
#define MEPS 1e-10
#define FPMIN 1.0e-30
#define MPI 3.1415926
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define IDEX(a, b) (((a)<(b))?1:0)

void QRDecompN(double *E, double *R, double *x, int n, int p);

void _BandMatrix(double *sighalf, double *rho, int p, int T);

#endif // RBS_H_INCLUDED
