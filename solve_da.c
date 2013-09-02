#include <math.h>
#include "fit.h"

void gaussj();

/* This program solves the equation alpha*da = beta for da. */
/* alpha is a matrix and da and beta are vectors */
/* covar is used for temporary space */

void solve_for_da(alpha, covar, beta, da, mfit)
double **alpha, **covar, *beta, *da;
int mfit;
{
int i, j;
double **mbeta;

mbeta = dmatrix(mfit,1);

for(i=0;i<mfit;i++){
	mbeta[i][0] = beta[i];
	for(j=0;j<mfit;j++)	covar[i][j] = alpha[i][j];
}
gaussj(covar,mfit,mbeta,1);
for(i=0; i<mfit; i++) da[i] = mbeta[i][0];
free_dmatrix(mbeta,mfit,1);
}

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

/* gaussj inverts the matrix a and performs */
/* the same operations on the matrix b */
/* This function is very similar to the one in */
/* Numerical Recipes of the same name. */

void gaussj(a,n,b,m)
double **a,**b;
int n,m;
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;

	indxc=ivector(n);
	indxr=ivector(n);
	ipiv=ivector(n);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) myerror("GAUSSJ: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) myerror("GAUSSJ: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free(ipiv);
	free(indxr);
	free(indxc);
}

#undef SWAP

