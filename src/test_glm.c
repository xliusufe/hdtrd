#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <float.h>
#include <string.h>
#include "testhb.h"

double _RPtest0(double *y, double *x, int p, int n){
	int i,j;
	double tmp, nume = 0.0, deno=0.0, Tn, *Q, *R;

	Q  = (double*)malloc(sizeof(double)*n*p);
	R  = (double*)malloc(sizeof(double)*p*p);


	QRDecompN(Q, R, x, n, p);

	for(i=0;i<n;i++){
		deno += y[i]*y[i];
	}
	for(j=0;j<p;j++){
		tmp = 0.0;
		for(i=0;i<n;i++){
			tmp += Q[j*n+i]*y[i];
		}
		nume += tmp*tmp;
	}
	deno -= nume;
	Tn = (n-p)*nume/deno/p;

	free(Q);
	free(R);

	return(Tn);

}

double _GCtest0(double *y, double *x, int p, int n){
	int i,j,k;
	double tmp, Tn=0.0, trS2=0.0;

	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			tmp = 0.0;
			for(k=0;k<p;k++){
				tmp	+= x[k*n+i]*x[k*n+j];
			}
			tmp *= y[i]*y[j];
			Tn 	+= tmp;
			trS2 += tmp*tmp;
		}
	}
	Tn *= 2.0/n;
	trS2 *= 2.0/(n*(n-1));
	return(Tn/sqrt(2*trS2));
}

double _GCtest(double *y, double *x, double *psi, int p, int n){
	int i,j,k;
	double tmp, Tn=0.0, trS2=0.0;

	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			tmp = 0.0;
			for(k=0;k<p;k++){
				tmp	+= x[k*n+i]*x[k*n+j];
			}
			tmp *= y[i]*y[j]*psi[i]*psi[j];
			Tn 	+= tmp;
			trS2 += tmp*tmp;
		}
	}
	Tn *= 2.0/n;
	trS2 *= 2.0/(n*(n-1));
	return(Tn/sqrt(2*trS2));
}

double _Test_RD_Max(const double *x, const double *y, double sigma2, double lambda0, const int n, int p){
	int i,j,k;
	double *meanx, *xtx, *normx,*barxtx;
	double tmp, trSigma2=0.0, Tn=0.0, barx2=0.0, deltax, a1, a2, a3;

	meanx	= (double*)malloc(sizeof(double)*p);
	xtx		= (double*)malloc(sizeof(double)*n*n);
	normx	= (double*)malloc(sizeof(double)*n);
	barxtx  = (double*)malloc(sizeof(double)*n);

	for(j=0;j<p;j++){
		tmp = 0.0;
		for(i=0;i<n;i++){
			tmp += x[j*n+i];
		}
		meanx[j]	= tmp/n;
		barx2		+= meanx[j]*meanx[j];
	}
	for(i=0;i<n;i++){
		tmp = 0.0;
		for(k=0;k<p;k++){
			tmp	+= x[k*n+i]*x[k*n+i];
		}
		normx[i] = tmp;
		tmp = 0.0;
		for(k=0;k<p;k++){
			tmp	+= x[k*n+i]*meanx[k];
		}
		barxtx[i] = tmp;
	}
	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			tmp = 0.0;
			for(k=0;k<p;k++){
				tmp	+= x[k*n+i]*x[k*n+j];
			}
			xtx[i*n+j] = tmp;
			deltax	= tmp*(1-1.0/n) - barxtx[i] - barxtx[j] + barx2 + (normx[i]+normx[j])/(2*n);
			Tn 		+= deltax*y[i]*y[j];
		}
	}
	Tn	*= 2.0*n/(n-1)/(n-2)/(n-2);
	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			xtx[j*n+i]	= xtx[i*n+j];
		}
	}

	a1 = a2 = a3 = tmp = 0.0;
	for(i=1;i<n;i++)
		for(j=0;j<i;j++){
			tmp +=	xtx[i*n+j];
			a1	+=	xtx[i*n+j]*xtx[i*n+j];
		}
	a1 *= 2;

	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			for(k=0;k<n;k++){
				if((k!=i)&(k!=j)){
					a2	+=	xtx[i*n+k]*xtx[j*n+k];
				}
			}
		}
	}
	a2	*= 2;
	a3 	= 4*tmp*tmp - 2*a1 - 4*a2;

	trSigma2 = a1  + (a3/(n-3) - 2*a2)/(n-2);
	trSigma2 /= n*(n-1);

	Tn -= lambda0;
	Tn *= n/sqrt(2.0*trSigma2)/sigma2;

	free(meanx);
	free(xtx);
	free(normx);
	free(barxtx);
	return Tn;
}

double _Test_RD_MAX_GC(double *x, double *y, double sigma2, double lambda0, int n, int p, int ishd){
	int i,j,k;
	double tmp, Tn=0.0, trS2=0.0;

	if(!ishd){
		tmp = 1.0*p/n;
		Tn = _RPtest0(y, x, p, n);	
		Tn = (Tn-1-lambda0/(tmp*sigma2))*sqrt(0.5*n*tmp*(1-tmp));
		return(Tn);
	}
	for(i=1;i<n;i++){
		for(j=0;j<i;j++){
			tmp = 0.0;
			for(k=0;k<p;k++){
				tmp	+= x[k*n+i]*x[k*n+j];
			}
			tmp *= y[i]*y[j];
			Tn 	+= tmp;
			trS2 += tmp*tmp;
		}
	}
	Tn *= 2.0/(n*(n-1));
	Tn -= lambda0;
	trS2 *= 2.0/(n*(n-1));
	return(n*Tn/sqrt(2*trS2)/sigma2);
}

void _mpequation(double *xr, double *xi, double *yr, double *yi, double *zr, double *zi, double *d, double *t, int n, int p, int J, int K){
	int i,j;
	double *mr, *mi, *pt, tau=1.0*p/n, znorm2, mnorm, tmp, tmp1, tmp2;

	mr	= (double*)malloc(sizeof(double)*K);
	mi	= (double*)malloc(sizeof(double)*K);

	tmp = d[0] - d[n-1];
	t[0] = d[n-1];
	for(i=1,pt=t+1;i<J-1;pt++,i++){
		*pt = *(pt-1) + tmp/(J-1);
	}
	*pt = d[0];

	for(i=0;i<K;i++){
		znorm2 = 1.0/(zr[i]*zr[i] + zi[i]*zi[i]);
		tmp1 = tmp2 = 0.0;
		for(j=0;j<p;j++){
			tmp = d[j] - zr[i];
			tmp1 += tmp/(tmp*tmp + zi[i]*zi[i]);
			tmp2 += zi[i]/(tmp*tmp + zi[i]*zi[i]);
		}
		mr[i] = -(1-tau)*zr[i]*znorm2 + tmp1/n;
		mi[i] = (1-tau)*zi[i]*znorm2 + tmp2/n;
	}

	for(i=0;i<K;i++){
		znorm2 = 1.0/(mr[i]*mr[i] + mi[i]*mi[i]);
		yr[i] = zr[i] + mr[i]*znorm2;
		yi[i] = zi[i] - mi[i]*znorm2;

		for(j=0;j<J;j++){
			tmp1 = 1.0+t[j]*mr[i];
			tmp2 = t[j]*mi[i];
			mnorm = tau/(tmp1*tmp1 + tmp2*tmp2);

			tmp = t[j]*t[j];
			xr[j*K+i] = (t[j]+tmp*mr[i])*mnorm;
			xi[j*K+i] = -tmp*mi[i]*mnorm;
		}
	}

	free(mr);
	free(mi);
}

void _calMoments(double *V, double *moms, double *ck, double *Sn, double *eigs, double *tJ, int n, int p, int J, int K0, int has_tJ){
	int i,j,k,s,np=MIN(n,p);
	double tmp,choosek=1.0, tr, *G, *Gk, *Gk1, *pt, *ps,*pj;
	G	= (double*)malloc(sizeof(double)*n*n);
	Gk	= (double*)malloc(sizeof(double)*n*n);
	Gk1	= (double*)malloc(sizeof(double)*n*n);

	if(!has_tJ){
		tmp = eigs[0] - eigs[n-1];
		tJ[0] = eigs[n-1];
		for(i=1,pt=tJ+1;i<J-1;pt++,i++){
			*pt = *(pt-1) + tmp/(J-1);
		}
		*pt = eigs[0];
	}

	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			G[j*n+i] = Sn[j*n+i];
		}
	}

	for(pt=Gk,ps=Sn;ps<Sn+n*n;pt++,ps++){
		*pt = *ps;
	}	
	tmp = 0.0;
	for(pt=eigs;pt<eigs+np;pt++){
		tmp += *pt;
	}
	moms[0] = tmp/p;
	ck[0] 	= 4.0/sqrt(1.0*n*p);

	for(j=0;j<J;j++){
		V[j] = tJ[j];
	}

	for(k=1;k<K0;k++){
		choosek *= 1.0*(n-k)/(k+1);
		for(i=0;i<n-1;i++){
			for(j=0;j<n;j++){
				tmp = 0.0;
				for(s=i+1;s<n;s++){
					tmp += G[s*n+i]*Gk[j*n+s];
				}
				Gk1[j*n+i] = tmp;
			}
		}
		for(pj=Gk1,pt=Gk;pt<Gk+n*n;pt++,pj++){
			*pt = n*(*pj);
		}
		tr = 0.0;
		for(i=0;i<n;i++){
			tr += Gk[i*n+i];
		}	
		moms[k] = tr/choosek/p;
		for(j=0;j<J;j++){
			V[k*J+j] = V[(k-1)*J+j]*tJ[j];
		}
		ck[k] = MAX(1.0, pow(1.0*p,0.5*(k+1)-1))*pow(2.0*(k+1),2.0*(k+1))/pow(1.0*n,0.5*(k+1));
	}

	free(G);
	free(Gk);
	free(Gk1);
}

void _mpmom(double *moms, double *E, double *G, double *H, double *f, double *Sn, double *eigs, double *tJ, int n, int p, int J, int K0, int has_tJ){
	int i,j,k,nr=K0+J;
	double *V, *ck, *pt, *ps;
	V	= (double*)malloc(sizeof(double)*K0*J);
	ck	= (double*)malloc(sizeof(double)*K0);

	_calMoments(V, moms, ck, Sn, eigs, tJ, n, p, J, K0, has_tJ);


	for(pt=E; pt<E+J; pt++){
		*pt = 1.0;
	}
	for(; pt<E+J+K0; pt++){
		*pt = 0.0;
	}

	for(pt=H, ps=moms; ps<moms+K0; pt++,ps++){
		*pt = -(*ps);
	}
	for(ps=moms; ps<moms+K0; pt++,ps++){
		*pt = *ps;
	}

	for(pt=f; pt<f+J; pt++){
		*pt = 0.0;
	}
	for(i=0; pt<f+J+K0; pt++, i++){
		*pt = 1.0/(ck[i]*moms[i]);
	}

	for(k=0;k<K0;k++){
		for(j=0;j<J;j++){
			G[k*nr+j] = -V[k*J+j];
		}
		for(j=0;j<K0;j++){
			G[k*nr+J+j] = 0.0;
		}
		G[k*nr+J+k] = 1.0;
	}
	for(k=0;k<K0;k++){
		for(j=0;j<J;j++){
			G[(k+K0)*nr+j] = V[k*J+j];
		}
		for(j=0;j<K0;j++){
			G[(k+K0)*nr+J+j] = 0.0;
		}
		G[(k+K0)*nr+J+k] = 1.0;
	}

	free(V);
	free(ck);
}

void _mplp4(double *E, double *F, double *G, double *H, double *f, double *xr, double *xi, double *yr, double *yi,
			double *zr, double *zi, double *moms, double *eigs, double *tJ, int n, int p, int J, int M, int K0){
	int i,j,k, K=2*K0,nr=K+J,np=MIN(n,p),n2=n*n+n+2;	
	double tmp, *V, *pt, *ps;
	double beta1,beta2,beta3,beta4,cn=1.0*p/n,tau2,tau3,tau4;
	V	= (double*)malloc(sizeof(double)*(M+1)*J);

	beta1=beta2=beta3=beta4=0.0;
	for(ps=eigs;ps<eigs+np;ps++){
		beta1 += (tmp = *ps);
		beta2 += (tmp *= *ps);
		beta3 += (tmp *= *ps);
		beta4 += (tmp *= *ps);
	}
	beta1 /= p;
	beta2 /= p;
	beta3 /= p;
	beta4 /= p;

	tau2	= 1.0*n*n/((n-1)*(n+2));
	tau3	= 1.0*tau2*n*n/((n-2)*(n+4));
	tau4	= 1.0*tau3*n*n2/((n-3)*(n+6)*(n+1));

	moms[0] = beta1;
	moms[1] = tau2*(beta2-cn*beta1*beta1);
	moms[2] = tau3*(beta3-3*cn*beta2*beta1+2*cn*cn*beta1*beta1*beta1);
	moms[3] = tau4*(beta4-4*cn*beta3*beta1-cn*beta2*beta2*(2*n*n+3*n-6)/n2+(2*beta2-cn*beta1*beta1)*n*(5*n+6)*cn*cn*beta1*beta1/n2);


	_mpequation(xr, xi, yr, yi, zr, zi, eigs, tJ, n, p, J, K0);

	// E
	for(pt=V; pt<V+J; pt++){
		*pt = 1.0;
	}
	for(k=1;k<M+1;k++){
		for(j=0;j<J;j++){
			V[k*J+j] = V[(k-1)*J+j]*tJ[j];
		}
	}
	for(k=0,ps=V;k<M+1;k++){
		pt = E + k*(J+K);
		for(i=0; i<J; i++,pt++){
			*pt = *ps++;
		}
		pt = E + k*(J+K)+J;
		for(i=0; i<K; i++,pt++){
			*pt = 0.0;
		}
	}

	// F
	F[0] = 1.0;
	for(pt=F+1,ps=moms; ps<moms+M; pt++,ps++){
		*pt = *ps;
	}	

	// H
	for(pt=H, ps=yr; ps<yr+K0; pt++,ps++){
		*pt = -(*ps);
	}
	for(ps=yi; ps<yi+K0; pt++,ps++){
		*pt = -(*ps);
	}
	for(ps=yr; ps<yr+K0; pt++,ps++){
		*pt = *ps;
	}
	for(ps=yi; ps<yi+K0; pt++,ps++){
		*pt = *ps;
	}

	// f
	for(pt=f; pt<f+J; pt++){
		*pt = 0.0;
	}
	for(i=0; pt<f+J+K; pt++, i++){
		*pt = 1.0;
	}

	// G
	for(k=0;k<K0;k++){
		for(j=0;j<J;j++){
			G[k*nr+j] = -xr[j*K0+k];
		}
		for(j=0;j<K0;j++){
			G[k*nr+J+j] = 0.0;
		}
		G[k*nr+J+k] = 1.0;
	}
	for(k=K0;k<K;k++){
		for(j=0;j<J;j++){
			G[k*nr+j] =  -xi[j*K0+k-K0];
		}
		for(j=K0;j<K;j++){
			G[k*nr+J+j] = 0.0;
		}
		G[k*nr+J+k] = 1.0;
	}

	for(k=0;k<K0;k++){
		for(j=0;j<J;j++){
			G[(k+K)*nr+j] = xr[j*K0+k];
		}
		for(j=0;j<K0;j++){
			G[(k+K)*nr+J+j] = 0.0;
		}
		G[(k+K)*nr+J+k] = 1.0;
	}
	for(k=K0;k<K;k++){
		for(j=0;j<J;j++){
			G[(k+K)*nr+j] = xi[j*K0+k-K0];
		}
		for(j=K0;j<K;j++){
			G[(k+K)*nr+J+j] = 0.0;
		}
		G[(k+K)*nr+J+k] = 1.0;
	}

	free(V);
}

SEXP BANDMATRIX_(SEXP RHO_, SEXP PARA_INT_){
	int p, T;
	p 	= INTEGER(PARA_INT_)[0];
	T 	= INTEGER(PARA_INT_)[1];

	// Outcome
	SEXP _sighalf;
	PROTECT(_sighalf	= allocVector(REALSXP, 	p*(p+T-1)));

	_BandMatrix(REAL(_sighalf), REAL(RHO_), p, T);

	UNPROTECT(1);
	return _sighalf;
}

SEXP GCtest_(SEXP X_, SEXP Y_, SEXP PSI_, SEXP DIM_){
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];
	int ispsi  		= dims[2];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *psi 	= REAL(PSI_);

	// Outcome
	SEXP rTn;
	PROTECT(rTn 		= allocVector(REALSXP, 1));

	if(ispsi==0){
		REAL(rTn)[0] = _GCtest0(y, x, p, n);
	}
	else{
		REAL(rTn)[0] = _GCtest(y, x, psi, p, n);
	}

	UNPROTECT(1);
	return rTn;

}

SEXP RDtest_MAX_(SEXP X_, SEXP Y_, SEXP SIGMA2_, SEXP LAMBDA0, SEXP DIM_){
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double sigma2 	= REAL(SIGMA2_)[0];
	double lambda0 	= REAL(LAMBDA0)[0];

	// Outcome
	SEXP rTn;
	PROTECT(rTn 		= allocVector(REALSXP, 1));

	REAL(rTn)[0] = _Test_RD_Max(x, y, sigma2, lambda0, n, p);

	UNPROTECT(1);
	return rTn;

}

SEXP RDtest_MAX_GC_(SEXP X_, SEXP Y_, SEXP SIGMA2_, SEXP LAMBDA0, SEXP DIM_){
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];
	int ishd 		= dims[2];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double sigma2 	= REAL(SIGMA2_)[0];
	double lambda0 	= REAL(LAMBDA0)[0];

	// Outcome
	SEXP rTn;
	PROTECT(rTn 		= allocVector(REALSXP, 1));

	REAL(rTn)[0] = _Test_RD_MAX_GC(x, y, sigma2, lambda0, n, p, ishd);

	UNPROTECT(1);
	return rTn;

}

SEXP _MP_EQUATION(SEXP ZR, SEXP ZI, SEXP EIGS, SEXP DIMs){
	int n, p, J, K;
	double *zr, *zi, *eigs;
	n 		= INTEGER(DIMs)[0];
	p		= INTEGER(DIMs)[1];
	K		= INTEGER(DIMs)[2];
	J		= INTEGER(DIMs)[3];

	zr 		= REAL(ZR);
	zi  	= REAL(ZI);
	eigs  	= REAL(EIGS);
	// Outcome
	SEXP _output, rtJ, rXr, rXi, rYr, rYi, _r_names;
  	PROTECT(_output		= allocVector(VECSXP, 	5));
  	PROTECT(_r_names   	= allocVector(STRSXP, 	5));
	PROTECT(rtJ 		= allocVector(REALSXP, 	J));
	PROTECT(rXr 		= allocVector(REALSXP, 	K*J));
	PROTECT(rXi 		= allocVector(REALSXP, 	K*J));
	PROTECT(rYr			= allocVector(REALSXP, 	K));
	PROTECT(rYi			= allocVector(REALSXP, 	K));


	_mpequation(REAL(rXr), REAL(rXi), REAL(rYr), REAL(rYi), zr, zi, eigs, REAL(rtJ), n, p, J, K);

  	SET_STRING_ELT(_r_names,	0,	mkChar("xr"));
  	SET_STRING_ELT(_r_names,	1,	mkChar("xi"));
	SET_STRING_ELT(_r_names,	2,	mkChar("yr"));
	SET_STRING_ELT(_r_names,	3,	mkChar("yi"));
	SET_STRING_ELT(_r_names,	4,	mkChar("tJ"));
	SET_VECTOR_ELT(_output, 	0, 	rXr);
	SET_VECTOR_ELT(_output, 	1, 	rXi);
	SET_VECTOR_ELT(_output, 	2, 	rYr);
	SET_VECTOR_ELT(_output, 	3, 	rYi);
	SET_VECTOR_ELT(_output, 	4, 	rtJ);
	setAttrib(_output, 			R_NamesSymbol,	_r_names);

	UNPROTECT(7);
	return _output;
}

SEXP _CAL_MOMENTS(SEXP SN, SEXP EIGS, SEXP DIMs){
	int n, p, J, K, has_tJ;
	double *Sn, *eigs;
	n 		= INTEGER(DIMs)[0];
	p		= INTEGER(DIMs)[1];
	K		= INTEGER(DIMs)[2];
	J		= INTEGER(DIMs)[3];
	has_tJ	= INTEGER(DIMs)[4];

	Sn 		= REAL(SN);
	eigs  	= REAL(EIGS);

	// Outcome
	SEXP _output, rV, rmoms, rck, rtJ, _r_names;
  	PROTECT(_output		= allocVector(VECSXP, 	4));
  	PROTECT(_r_names   	= allocVector(STRSXP, 	4));
	PROTECT(rV 			= allocVector(REALSXP, 	K*J));
	PROTECT(rmoms 		= allocVector(REALSXP, 	K));
	PROTECT(rck			= allocVector(REALSXP, 	K));
	PROTECT(rtJ			= allocVector(REALSXP, 	J));


	_calMoments(REAL(rV), REAL(rmoms), REAL(rck), Sn, eigs, REAL(rtJ), n, p, J, K, has_tJ);

  	SET_STRING_ELT(_r_names,	0,	mkChar("V"));
  	SET_STRING_ELT(_r_names,	1,	mkChar("moms"));
	SET_STRING_ELT(_r_names,	2,	mkChar("ck"));
	SET_STRING_ELT(_r_names,	3,	mkChar("tJ"));
	SET_VECTOR_ELT(_output, 	0, 	rV);
	SET_VECTOR_ELT(_output, 	1, 	rmoms);
	SET_VECTOR_ELT(_output, 	2, 	rck);
	SET_VECTOR_ELT(_output, 	3, 	rtJ);
	setAttrib(_output, 			R_NamesSymbol,	_r_names);

	UNPROTECT(6);
	return _output;
}

SEXP _ARGUEMENTS_MPMOM(SEXP SN, SEXP EIGS, SEXP DIMs){
	int n, p, J, K, has_tJ;
	double *Sn, *eigs;
	n 		= INTEGER(DIMs)[0];
	p		= INTEGER(DIMs)[1];
	K		= INTEGER(DIMs)[2];
	J		= INTEGER(DIMs)[3];
	has_tJ	= INTEGER(DIMs)[4];

	Sn 		= REAL(SN);
	eigs  	= REAL(EIGS);

	// Outcome
	SEXP rmoms, rtJ, rE, rG, rH, rf, _output, _r_names;
  	PROTECT(_output		= allocVector(VECSXP, 	6));
  	PROTECT(_r_names   	= allocVector(STRSXP, 	6));
	PROTECT(rmoms 		= allocVector(REALSXP, 	K));
	PROTECT(rtJ 		= allocVector(REALSXP, 	J));
	PROTECT(rE 			= allocVector(REALSXP, 	J+K));
	PROTECT(rG			= allocVector(REALSXP, 	2*K*(J+K)));
	PROTECT(rH			= allocVector(REALSXP, 	2*K));
	PROTECT(rf			= allocVector(REALSXP, 	J+K));

	_mpmom(REAL(rmoms), REAL(rE), REAL(rG), REAL(rH), REAL(rf), Sn, eigs, REAL(rtJ), n, p, J, K, has_tJ);

	SET_STRING_ELT(_r_names,	0,	mkChar("moms"));
	SET_STRING_ELT(_r_names,	1,	mkChar("tJ"));
  	SET_STRING_ELT(_r_names,	2,	mkChar("E"));
	SET_STRING_ELT(_r_names,	3,	mkChar("G"));
	SET_STRING_ELT(_r_names,	4,	mkChar("H"));
	SET_STRING_ELT(_r_names,	5,	mkChar("f"));
	SET_VECTOR_ELT(_output, 	0, 	rmoms);
	SET_VECTOR_ELT(_output, 	1, 	rtJ);
	SET_VECTOR_ELT(_output, 	2, 	rE);
	SET_VECTOR_ELT(_output, 	3, 	rG);
	SET_VECTOR_ELT(_output, 	4, 	rH);
	SET_VECTOR_ELT(_output, 	5, 	rf);
	
	setAttrib(_output, 			R_NamesSymbol,	_r_names);

	UNPROTECT(8);
	return _output;
}

SEXP _ARGUEMENTS_MPLP4(SEXP ZR, SEXP ZI, SEXP EIGS, SEXP DIMs){
	int n, p, J, M, K, K0;
	double *zr, *zi, *eigs;
	n 		= INTEGER(DIMs)[0];
	p		= INTEGER(DIMs)[1];
	M		= INTEGER(DIMs)[2];
	J		= INTEGER(DIMs)[3];
	K0		= INTEGER(DIMs)[4];

	zr 		= REAL(ZR);
	zi  	= REAL(ZI);
	eigs  	= REAL(EIGS);

	// Outcome
	K = 2*K0;
	SEXP rtJ, rE, rF, rG, rH, rf, rXr, rXi, rYr, rYi, rmoms, _output, _r_names;
  	PROTECT(_output		= allocVector(VECSXP, 	11));
  	PROTECT(_r_names   	= allocVector(STRSXP, 	11));
	PROTECT(rtJ			= allocVector(REALSXP, 	J));
	PROTECT(rE 			= allocVector(REALSXP, 	(M+1)*(J+K)));
	PROTECT(rF			= allocVector(REALSXP, 	1+M));
	PROTECT(rG			= allocVector(REALSXP, 	2*K*(J+K)));
	PROTECT(rH			= allocVector(REALSXP, 	2*K));
	PROTECT(rf			= allocVector(REALSXP, 	J+K));

	PROTECT(rXr 		= allocVector(REALSXP, 	K0*J));
	PROTECT(rXi 		= allocVector(REALSXP, 	K0*J));
	PROTECT(rYr			= allocVector(REALSXP, 	K0));
	PROTECT(rYi			= allocVector(REALSXP, 	K0));
	PROTECT(rmoms		= allocVector(REALSXP, 	4));


	_mplp4(REAL(rE), REAL(rF), REAL(rG), REAL(rH), REAL(rf), REAL(rXr), REAL(rXi), REAL(rYr), REAL(rYi), zr, zi, REAL(rmoms), eigs, REAL(rtJ), n, p, J, M, K0);


	SET_STRING_ELT(_r_names,	0,	mkChar("tJ"));
  	SET_STRING_ELT(_r_names,	1,	mkChar("E"));
	SET_STRING_ELT(_r_names,	2,	mkChar("G"));
	SET_STRING_ELT(_r_names,	3,	mkChar("H"));
	SET_STRING_ELT(_r_names,	4,	mkChar("f"));
	SET_STRING_ELT(_r_names,	5,	mkChar("F"));
  	SET_STRING_ELT(_r_names,	6,	mkChar("xr"));
  	SET_STRING_ELT(_r_names,	7,	mkChar("xi"));
	SET_STRING_ELT(_r_names,	8,	mkChar("yr"));
	SET_STRING_ELT(_r_names,	9,	mkChar("yi"));
	SET_STRING_ELT(_r_names,	10,	mkChar("moms"));

	SET_VECTOR_ELT(_output, 	0, 	rtJ);
	SET_VECTOR_ELT(_output, 	1, 	rE);
	SET_VECTOR_ELT(_output, 	2, 	rG);
	SET_VECTOR_ELT(_output, 	3, 	rH);
	SET_VECTOR_ELT(_output, 	4, 	rf);
	SET_VECTOR_ELT(_output, 	5, 	rF);
	SET_VECTOR_ELT(_output, 	6, 	rXr);
	SET_VECTOR_ELT(_output, 	7, 	rXi);
	SET_VECTOR_ELT(_output, 	8, 	rYr);
	SET_VECTOR_ELT(_output, 	9, 	rYi);
	SET_VECTOR_ELT(_output, 	10, rmoms);

	setAttrib(_output, 			R_NamesSymbol,	_r_names);

	UNPROTECT(13);
	return _output;
}
