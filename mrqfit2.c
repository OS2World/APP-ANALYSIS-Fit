#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fit.h"

/* calculates alpha, beta, and chisq */
/* see Numerical Recipes in C for definition of these */
/* the Levenberg Marquardt algorithm requires us to slove */
/* the equation alpha*da = beta (where alpha is a matrix */
/* and da and beta are vectors) for da */
/* this function does this by gauss-jordan */
/* elimination */

/* declarations for data declared externally */
extern int grflag;
extern int gnuopen;
extern FILE *pstream;
extern double *xmin, *xmax;
extern int wiflag;
extern int veflag;
extern int debug;
extern char *GNUPLOT;
extern int MYPLOT;
extern int doflag;

int mrqfit(data,order,num_indep,ndata,itmax,a,ma,lista,mfit,
				covar,chisq,func,filename,comment)
double **data;   
struct data_order order;
int num_indep;
int ndata;
int itmax;
double *a;
int ma, *lista, mfit;
double **covar, *chisq;
int (*func)();
char *filename;
char *comment;
{
int i, j, k;  /* indices for loops */
int failed;
double **alpha;
/* alamda determines how much we change the fitting parameters */
double alamda;
char cmd[20]="";
double *beta;
double *atry; /* parameters to try and see if we reduce chisqr */
double **alpha_try;
double ochisq; /* best chisqr so far */
double *da;
FILE *stream;

/* allocate space for arrays */
atry = dvector(ma);
alpha_try=dmatrix(ma,ma);
beta = dvector(ma);
da = dvector(ma);
alpha=dmatrix(ma,ma);

stream=NULL;

/* calculate alpha, beta, and chisq for current value of parameters */
if(debug) printf("in mrqfit, calling alpha_beta_chisqr()\n");
failed=alpha_beta_chisq(data,order,num_indep,ndata,a,ma,lista,mfit,
							alpha,beta,chisq,func);
if(debug) printf("in mrqfit, alpha_beta_chisqr() returned %d\n", failed);
if(failed){
	free_dmatrix(alpha,ma,ma);
	free(atry);
	free_dmatrix(alpha_try,ma,ma);
	free(beta);
	free(da);
	return 0;
}

/* start alamda small */
alamda = 1e-3;
i = 0;

/* loop until user tells us to quit */
while( strcmp(cmd,"q") != 0){

/* loop specified number of iterations */
while( i <= itmax ){
	ochisq = *chisq;
	for(j=0;j<mfit;j++){
		for(k=0;k<mfit;k++){
			alpha_try[j][k] = alpha[j][k];
		}
		alpha_try[j][j] = alpha[j][j] + alamda;
	}
	if(debug)printf("in mrqfit, calling solve_for_da()\n");
	solve_for_da(alpha_try, covar, beta, da, mfit);
	if(debug)printf("in mrqfit, returned from solve_for_da()\n");
 
	for(j=0;j<ma;j++) atry[j] = a[j];
	for(j = 0; j < mfit; j++) atry[lista[j]] += da[j];

if(veflag==2){
		printf("\nalpha: ");
		for(j=0;j<mfit;j++){
			printf("\n %d   ",lista[j]);
			for(k=0;k<mfit;k++) printf(" %g ",alpha_try[j][k]);
		}
		printf("\ncovar: ");
		for(j=0;j<mfit;j++){
			printf("\n %d   ",lista[j]);
			for(k=0;k<mfit;k++) printf("\n covar[%d][%d]= %g ",lista[j], lista[k], covar[j][k]);
		}
		printf("\natry:");
		for(j=0;j<mfit;j++) printf("\n atry%d = %g",lista[j], atry[lista[j]]);
		printf("\n");
	}

	if(debug)printf("in mrqfit, calling alpha_beta_chisqr()\n");
	failed=alpha_beta_chisq(data,order,num_indep,ndata,atry,ma,
							lista,mfit,alpha,beta,chisq,func);
	if(debug)printf("in mrqfit, alpha_beta_chisqr() returned %d\n", failed);

if(failed){
	free_dmatrix(alpha,ma,ma);
	free(atry);
	free_dmatrix(alpha_try,ma,ma);
	free(beta);
	free(da);
	return 0;
}

	if(*chisq >= ochisq){
		alamda *= 10;
		if(alamda > 1e15) alamda /= 3e14;
		*chisq = ochisq;
	}
	else{
		if(debug){
			if(stream == NULL){
				printf("opening lamda.sts\n");
				stream = fopen("lamda.sts","at");
			}
			if(stream != NULL){
				fprintf(stream,"%g\n",alamda);
				printf("alamda: %g\n",alamda);
			}
		}
		alamda *= 0.1;
		if(alamda < 1e-15) alamda *= 3e14;
		for(j=0;j < mfit; j++)a[lista[j]] = atry[lista[j]];
	}

	if(veflag){
		for(j=0; j < ma; j++)printf("\n a%d = %g", j, a[j]);
		printf("\n chisqr = %g", *chisq);
		printf("\n alamda = %g", alamda);
		printf("\n");
	}
	else printf("iteration: %d chisqr: %g\n", i, *chisq);
	i++;

#ifdef OS2
	if( grflag && MYPLOT && (num_indep == 1) && veflag){
		if(debug)printf("in mrqfit, calling myplot()\n");
		failed=myplot(func,data, order, num_indep, ndata, filename, comment, a, ma);
		if(debug)printf("in mrqfit, myplot() returned %d\n", failed);
	}
#endif
}

if(grflag && ( !MYPLOT || !veflag || (num_indep != 1) )){
		if(debug)printf("in mrqfit, calling plot()\n");
		failed=plot(func,data, order, num_indep, ndata, filename, comment, a, ma);
		if(debug)printf("in mrqfit, plot() returned %d\n", failed);
}

if(doflag) strcpy(cmd, "q");
else{
	printf("q quit, g graphs, n does not graph, anything else continues ");
	fflush(stdout);
	gets(cmd);
}
i = 0;
if(strcmp(cmd,"g") == 0) grflag = 1;
if(strcmp(cmd,"n") == 0) grflag = 0;

}
free_dmatrix(alpha,ma,ma);
free(atry);
free_dmatrix(alpha_try,ma,ma);
free(beta);
free(da);
if(stream != NULL) fclose(stream);

return 0;
}
 

int alpha_beta_chisq(data, order,num_indep,ndata,a,ma,
								lista,mfit,alpha,beta,chisq,funcs)
double **data;
struct data_order order;
int num_indep;
int ndata;
double *a;
int ma, *lista,mfit;
double **alpha, *beta, *chisq;
int (*funcs)(double *,double *,double *, double *,int, int, 
					int *, int *, double *, double);
{
int i,j,k;
int failed;
double *x, *y, *sig,*xsig;
double ymod, *dyda, sig2i, wt,dy;
int *fita, *dydx_flag;
double *dydx;
double sigi;

if(debug)printf("in abc, allocating data\n");
dyda=dvector(ma);
fita = ivector(ma);
x = dvector(num_indep);
dydx_flag = ivector(num_indep);
dydx = dvector(num_indep);

if(debug)printf("in abc, setting pointers for y and sig\n");
y = data[order.y];
sig = data[order.sig];

if(debug)printf("in abc, setting pointers for xsig\n");
for(i = 0; i < num_indep; i++){
	if(order.xsig[i] >= 0){
		xsig = data[order.xsig[i]];
		dydx_flag[i] = 1;
	}
	else dydx_flag[i] = 0;
}
if(debug)printf("in abc, done setting pointers for xsig\n");

for(i = 0; i < ma; i++){
	fita[i] = 0;
	for(j = 0; j < mfit; j++) if(lista[j] == i) fita[i] = 1;
}

*chisq = 0;
for(j=0; j < mfit;j++){
	beta[j] = 0;
	for(k=0;k<mfit;k++) alpha[j][k] = 0;
}

if(debug)printf("in abc, setting pointers for x\n");
for(i=0;i<ndata;i++){
	for(j = 0; j < num_indep; j++){
		x[j] = data[order.x[j]][i];
	}
if(debug > 1)printf("in abc, done setting pointers for x\n");

	if(wiflag == 0 || window(x, num_indep)){
		if(debug > 1)printf("in abc, about to call funcs\n");
		failed = (*funcs)(x,a,&ymod,dyda,ma,1, fita, dydx_flag, dydx, y[i]);
		if(debug > 1)printf("in abc, returned from funcs\n");
		if(failed){
		free(dyda);
		free(fita);
		free(x);
		free(dydx_flag);
		free(dydx);
		return 1;
		}
		data[order.yfit][i] = ymod;
		sig2i = (sig[i]*sig[i]);
		for(j = 0; j < num_indep; j++)
			if(order.xsig[j] >= 0) sig2i += (dydx[j])*(dydx[j])*xsig[i]*xsig[i];
		if(debug==2) printf("i: %d x[0]: %g ymod: %g y: %g sig2i: %g\n", i, x[0], ymod, y[i], sig2i);
		dy = y[i] - ymod;
		for(j=0;j<mfit;j++){
			wt=dyda[lista[j]]/sig2i;
			for(k=0;k<=j;k++)alpha[j][k] += wt*dyda[lista[k]];
			beta[j] += dy*wt;
		}
		*chisq += dy*dy/sig2i;
	}
}

for(j=1;j<mfit;j++)
	for(k=0;k<=j-1;k++) alpha[k][j]=alpha[j][k];

free(dyda);
free(fita);
free(x);
free(dydx_flag);
free(dydx);
return 0;
}

int window(double *x, int num_indep){

int i;

for(i = 0; i < num_indep; i++)
	if(x[i] < xmin[i] || x[i] > xmax[i]) return 0;

return 1;
}

