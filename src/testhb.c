#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "testhb.h"


void QRDecompN(double *E, double *R, double *x, int n, int p){
	// Input:
	// X is a p*n matrix
	//
	// Output:
	// R is a p*p lower triangular matrix
	// E is a p*n matrix satisfying E*t(E) = I_p
	//
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;

	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p
	for(j=1;j<p;j++){
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); 
	free(znorm);
}

void _BandMatrix(double *sighalf, double *rho, int p, int T){
	int i,j,n=p+T-1;
	for(i=0;i<n*p;i++){
		sighalf[i] = 0.0;
	}
	for(j=0;j<p;j++){
		for(i=0;i<T;i++){
			sighalf[j*n+j+i] = rho[i];
		}
	}
}
