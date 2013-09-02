#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fit.h"

/* declarations for data declared externally */
extern int grflag;
extern int wiflag;
extern int veflag;
extern int debug;
extern int MYPLOT;

/* linear_fit() does a linear least squares fit to a specified function. */
/* because the functions are defined for non-linear fitting, we have */
/* to jump through some hoops to do the calculation. */


int linear_fit(double **data,struct data_order order,
				int num_indep, int ndata, int itmax, double *a,
				int ma, int *lista, int mfit, double **covar,
				double *chisq, int (*func)(),
				char *filename, char *comment){

int i, j, k, l, m;  /* indices for loops */
int failed=0;
double **alpha;     /* See Chapter 14 of Numerical Recipes */
double *beta;
double *X;
double *x;
double *dyda;       /* Not really needed, but function might return */
                    /* values for dyda[i], so we have to allocate space */

double *aj;         /* The function needs to return the value of each basis */
                    /* function at a given x, X[j](x[i]).  However, function */
                    /* is set up to return value of f(x,a).   */
                    /* To get X[j](x[i]), we send a parameter list where */
                    /* a[i] = 0 if a != j, and a[i] = 1 if a ==j. */
                    /* I call this parameter list aj to distinguish it from a. */

int *fita;          /* fita tells the user function  func() which dyda[i]'s to */
                    /* compute.  We don't need any of them for linear fitting */
                    /* so we set all the fita[i]'s to zero */

int *dydx_flag;     /* dydx_flag tells the user function which dydx[i]'s to compute. */
                    /* We can't use them for linear fitting, because we are not */
                    /* iterating and don't have a guess for the parameters until */
                    /* the fit is finished.  As a result, our linear fitting does */
                    /* not consider errors in the independent variables.  We set */
                    /* all of the dydx_flag[i]'s to zero, so the function does not */
                    /* have to waste time computing them. */

double *dydx;       /* However, since user functions might compute the dydx[i]'s */
                    /* anyway, we need to allocate storage. */

double sig2i;       /* sigmayi * sigmay1 */
double *atry;       /* These are the parameter values multiplying the basis functions */
                    /* used in the fit.  If a basis function is not used in the fit, */
                    /* it is multiplied by zero.  atry has dimension mfit and holds */
                    /* only the parameters multiplying basis functions which are fit */
double chisqr;


alpha=dmatrix(mfit,mfit);
X = dvector(ma);          /* value of basis function at each data point */
beta = dvector(mfit);
dyda=dvector(ma);
fita = ivector(ma);
x = dvector(num_indep);
dydx = dvector(num_indep);
dydx_flag = ivector(num_indep);
aj = dvector(ma);
atry = dvector(mfit);

for(j = 0; j < num_indep; j++) dydx_flag[j] = 0;  /* don't calculate dydx[i]'s */
for(j = 0; j < ma; j++) fita[j] = 0;              /* don't calculate dyda[i]'s */


/* set all alpha[i][j]'s equal to zero */

for(i = 0; i < mfit; i++){
	for(j = 0; j < mfit; j++){
		alpha[i][j] = 0;
	}
	beta[i] = 0;
}

/* loop summing over data points */
for(i = 0; i < ndata; i++){

	/* set up array of independent variable values at i'th data point */
	for(j = 0; j < num_indep; j++){
		x[j] = data[order.x[j]][i];
		if(debug > 1)
			printf("i: %d order.x[%d]: %d x[%d]: %g\n", i, j, order.x[j], j, x[j]);
	}

	/* If windowing is on, only consider data points within window */
	if(wiflag == 0 || window(x, num_indep)){

		/* calculate value of relevant basis functions at current data point */
		for(j = 0; j < mfit; j++){

			/* set up aj[i]'s so function value = basis function value */
			for(l=0; l < ma; l++){
				if(l != lista[j]) aj[l] = 0;
				else aj[l] = 1;
			}
			failed = (*func)(x,aj,&X[lista[j]],dyda,ma,0, fita, dydx_flag, dydx, data[order.y][i]);
			if(debug > 1) printf("X[%d]: %g ",lista[j], X[lista[j]]);
			if(failed) break;
		}
		if(debug > 1) printf("\n");
		if(failed) break;
 
		sig2i = data[order.sig][i]*data[order.sig][i];
		if(debug > 1) printf("order.sig: %d data[order.sig][i]: %g sig2i: %g\n",
									order.sig, data[order.sig][i], sig2i);

		/* add this data point to alpha and beta */
		for(k = 0; k < mfit; k++){
			for(j = 0; j <= k; j++){
/*				if(debug > 1) printf("alpha[%d][%d] += X[%d]*X[%d]/sig2i;\n",
					j, k, lista[j], lista[k]); */
				alpha[j][k] += X[lista[j]]*X[lista[k]]/sig2i;
			}
			beta[k] += data[order.y][i]*X[lista[k]]/sig2i;
		}
	}
}

if(debug) printf("done looping over data points\n");

if(!failed){

	if(debug){
		printf("\nalpha: ");
		for(j=0;j<mfit;j++){
			printf("\n %d   ",j);
			for(k=0;k<mfit;k++) printf(" %g ",alpha[j][k]);
		}
	}

	/* fill in rest of alpha, alpha is symmetric, and we only  */
	/*calculated upper triangular, we need the whole thing */
	for(k = mfit-1; k >= 0; k--){
		for(j = mfit-1; j > k; j--){
			alpha[j][k] = alpha[k][j];
		}
	}

	if(veflag==2){
		printf("\nalpha: ");
		for(j=0;j<mfit;j++){
			printf("\n %d   ",j);
			for(k=0;k<mfit;k++) printf(" %g ",alpha[j][k]);
		}
		printf("\n");
	}

	/* solve the matrix equation alpha*atry = beta, where alpha is */
	/* a matrix, and atry and beta are vectors */
	/* covar is the inverse of alpha */
	if(debug) printf("calling solve_for_da()\n");
	solve_for_da(alpha, covar, beta, atry, mfit);
	if(debug) printf("returned from solve_for_da()\n");

	if(veflag==2){
		printf("\ncovar: ");
		for(j=0;j<mfit;j++){
			printf("\n %d   ",j);
			for(k=0;k<mfit;k++) printf("\n covar[%d][%d]= %g ",j, k, covar[j][k]);
		}
	printf("\n");
	}

if(debug){
for(j=0; j < mfit; j++)printf("\n atry%d = %g", j, atry[j]);
printf("\n");
}

for(j = 0; j < ma; j++) a[j] = 0;
for(j = 0; j < mfit; j++) a[lista[j]] = atry[j];

/* call calc_yfit to compute chisqr */
failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisqr);


if(grflag){
	if(MYPLOT && (num_indep == 1)){
#ifdef OS2
		if(debug)printf("in linear_fit(), calling myplot()\n");
		failed=myplot(func,data, order, num_indep, ndata, filename, comment, a, ma);
		if(debug)printf("in linear_fit(), myplot() returned %d\n", failed);
#endif
	}
	else{
		if(debug)printf("in linear_fit(), calling plot()\n");
		failed=plot(func,data, order, num_indep, ndata, filename, comment, a, ma);
		if(debug)printf("in linear_fit(), plot() returned %d\n", failed);
	}
}

for(j=0; j < ma; j++)printf("\n a%d = %g", j, a[j]);
*chisq = chisqr;
printf("\n chisqr = %g", *chisq);
printf("\n");
}

free_dmatrix(alpha,mfit,mfit);
free(X);
free(beta);
free(dyda);
free(fita);
free(x);
free(dydx);
free(dydx_flag);
free(aj);
free(atry);
return failed;
}
