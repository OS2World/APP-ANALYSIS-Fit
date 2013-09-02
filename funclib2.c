#include <math.h>
#include <string.h>
#include "fit.h"
#define NUM_FCNS 20   /* CHANGE if you define more than 20 functions */

extern int debug;
/* a few function declarations */
int fgauss();
int fgaussc();
int fgaussn();
int florenz();
int florenz2();
int fline();
int fpoly();
int fexpn();
int fxyquad();
int fxygauss();
int fconic();
int fsincos();

/* The structure fcn holds the information about each function */
/* Specifically, it is here that the name used by the command line, */
/* the name of the function in c, the number of parameters, */
/* and a comment are defined.  If gnuplot is used for plotting, */
/* then the comment must be a valid definition of the function */
/* on the gnuplot command line. */
/* to add a function, you must add a line to the fcn initialization */
/* For a function with a variable number of parameters */
/* you must initialize na to be -1.  You have to choose the number */
/* of parameters when you choose the function at the fit> command */
/* prompt. */
/* If a comment begins with "FILE", plotting is done through a */
/* tempory file.  For functions with variable number of parameters, */
/* you must begin your comment with "FILE" */
/* Functions of two independent variables which are not plotted */
/* through a file should be expressed in terms of u and v instead */
/* of x and y.  This is because gnuplot wants it that way. */

/* The elements of the struct routine are as follows: */
/* char name[20]    the name of the function at the fit> command line */
/* int num_indep   the number of independent variables */
/* int linflag   Is the function linear in the parameters? */
/* linflag = 0 for a non-linear function.  linflag = 1 for */
/* a function linear in the parameters */
/* int (*func)()   the name of the function as known to the */
/*                         compiler */
/* int na        the number of parameters, -1 for variable number */
/* char comment[160] a comment about the function */

struct routine{
	char name[20];
	int num_indep;
	int linflag;
	int (*func)(double *x, double *a, double *y, double *dyda,
           int na, int dyda_flag, int *fita, int dydx_flag, 
           double *dydx, double ydat);
	int na;
	char comment[160];
} fcn[NUM_FCNS]={

   {"gauss",1,0, &fgauss, 3, "f(x) = a2*exp(-((x-a0)/a1)**2)" },
   {"gaussc",1,0, &fgaussc, 4,"f(x) = a2*exp(-((x-a0)/a1)**2) + a3" },
   {"ngauss",1,0, &fgaussn, -1,"FILE f(x) = sum( a[i+2]*exp(-((x-ai)/a[i+1])**2) ) + a[n-1]"},
   {"lorenz", 1,0, &florenz,3,"f(x) = a1*a2/(4*(x-a0)**2 + a1)" },
   {"2lorenz",1,0,&florenz2,7,"f(x) = a1*a2/(4*(x-a0)**2 + a1) + a4*a5/(4*(x-a3)**2 + a4) + a6" },
   {"line",1,1,&fline, 2, "f(x) = a0 + a1*x"},
   {"poly", 1,1,&fpoly, -1, "FILE f(x) = sum( ai * pow(x,i) )"},
   {"nexp", 1,0,&fexpn, -1, "FILE f(x) = sum( a[i+2]*exp((x-a[i])/a[i+1]) ) + a[n-1]"},
   {"xyquad", 2,1,&fxyquad, 6, "f(u,v) = a0 + a1*u + a2*v + a3*u*u + a4*v*v + a5*u*v"},
   {"xygauss",2,0,&fxygauss, 6, "f(u,v) = a4*exp(-((u-a0)/a2)**2 - ((v-a1)/a3)**2 ) + a5"},
   {"conic",1,0,&fconic, 6, "FILE a0*x*x + a1*x*y + a2*y*y + a3*x + a4*y +a5 = 0"},
   {"sincos",1,0,&fsincos, -1, "FILE a0 + a[i]*sin(a[i+1]*x) + a[i+2]*cos(a[i+3]*x)"},
};

/* The function getfcnptr() looks through the struct fcn to find a */
/* pointer corresponding to the function specified at the fit> */
/* prompt with the fn command */

int (*getfcnptr(char *name, int *num_indep, int *linflag, int *na, char *comment))
{
int i=0;
while(fcn[i].na != 0 ){
	if(strcmp(name, fcn[i].name) == 0){
		*na = fcn[i].na;
		*linflag = fcn[i].linflag;
		*num_indep = fcn[i].num_indep;
		strcpy(comment,fcn[i].comment);
		return fcn[i].func;
	}
	i++;
}

/* If function name not found return pointer to NULL */
*na = 0;
return (int *)NULL;
}

/* listfcns() lists the functions available */
int listfcns(){
int i=0 ;
while(fcn[i].na != 0 ){
	printf("%s %s\n", fcn[i].name, fcn[i].comment);
	i++;
	if(i == 20){
		i = 0;
		gets();
	}
}
return 0;
}

/* the definition of a function must be of the form: */
/* int function(double *x, double *a, double *y, double *dyda, */
/*                      int na, int dyda_flag, int *fita, int dydx_flag, */
/*                      double *dydx, double ydat); */

/*
/* The x[i]'s are the independent variables passed to the function */
/* *y should equal the value of the function at x. */
/* The a[i]'s are the parameters. */
/* the dyda[i]'s should be set equal to the first derivative */
/* of the fitting function with respect to each parameter */
/* dyda_flag, fita[], and dydx_flag[] are flags which are */
/* passed to the function and used for optimization */
/* purposes only.  These need not be used at all, */
/* but may result in a significant performance boost. */
/* if calculating the derivative is slow and it is not */
/* needed by the calling function. */
/* dydx[i] should be set equal to the derivative of the fitting */
/* function with respect to x[i].  It is used only if */
/* the fitting algorithm is told to consider errors in the x[i] */
/* ydat is the data value that you are fitting to. */
/* It might be useful in multivalued functions */
/* Using it might destroy the statistical relevance of your fit. */

/* The linear fitting routine makes different use of user defined */
/* fitting functions.  the dyda[i]'s and dydx[i]'s are not needed. */
/* Linear fitting functions which are set up properly for non-linear fitting */
/* work fine with the linear fitting routine.  Using dydx_flag, */
/* dyda_flag, and fita[i] will speed up linear fitting. */

/* a conic section, doesn't work very well */
int fconic(x, a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[];
double a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double A, B, C, D, E, F, Q, R, S;
double yp, ym, x2;
int pm;

A = a[0];
B = a[1];
C = a[2];
D = a[3];
E = a[4];
F = a[5];

x2 = x[0]*x[0];
R = B*x[0] + E;
S = A*x2 + D*x[0] + F;
if( (R*R - 4*C*S) < 0 ){
	printf("Function conic failed, wanted to take sqrt of negative #.\n");
	return 1;
}
Q = sqrt(R*R - 4*C*S);

yp = (-R + Q)/2*C;
ym = (-R - Q)/2*C;

if( fabs(yp - ydat) < fabs(ym - ydat) ){
	pm = +1;
	*y = yp;
}
else{
	pm = -1;
	*y = ym;
}

if(dyda_flag){
	if(fita[0]) dyda[0] = -pm*x2/Q;
	if(fita[1]) dyda[1] = -x[0]/(2*C) + pm*R*x[0]/(2*Q*C);
	if(fita[2]) dyda[2] = (R/(2*C) + pm*S/Q)/C;
	if(fita[3]) dyda[3] = -pm*x[0]/Q;
	if(fita[4]) dyda[4] = -1/(2*C) + pm/(2*Q*C);
	if(fita[5]) dyda[5] = -pm/Q;
}
if(dydx_flag[0])
	dydx[0] = -B/(2*C) + pm*(2*B*(B*x[0]+E) - 
									4*C*(2*x[0]*A + D))/(2*C*Q);

return 0;
}

/*f(u,v) = a4*exp(-((u-a0)/a2)**2 - ((v-a1)/a3)**2 ) + a5 */
/* a two dimensional gaussian */
int fxygauss(x, a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[];
double a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double fac0, ex0, arg0,fac1, ex1, arg1;
arg0 = (x[0]-a[0])/a[2];
ex0 = exp(-arg0*arg0);
fac0 = a[4]*ex0*2.0*arg0;
arg1 = (x[1]-a[1])/a[3];
ex1 = exp(-arg1*arg1);
fac1 = a[4]*ex1*2.0*arg1;

*y = a[4]*ex0*ex1 + a[5];

if(dyda_flag){
	if(fita[0]) dyda[0] = ex1*fac0/a[2];
	if(fita[1]) dyda[1] = ex0*fac1/a[3];
	if(fita[2]) dyda[2] = ex1*fac0*arg0/a[2];
	if(fita[3]) dyda[3] = ex0*fac1*arg1/a[3];
	if(fita[4]) dyda[4] = ex0*ex1;
	if(fita[5]) dyda[5] = 1;
}
if(dydx_flag[0]) dydx[0] = -ex1*fac0/a[2];
if(dydx_flag[1]) dydx[1] = -ex0*fac1/a[3];
return 0;
}


/* a sum of sines and cosines */
int fsincos(x, a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[];
double a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{

int i;
double temp;

*y = a[0];
dydx[0] = 0;

dyda[0] = 1;
for(i = 1; i < na-3; i+=4){
	dyda[i] = sin(a[i+1]*x[0]);
	*y += a[i]*dyda[i];
	if(dydx_flag[0] || fita[i+1]) temp = a[i]*cos(a[i+1]*x[0]);
	dydx[0] += a[i+1]*temp;
	dyda[i+1] = x[0]*temp;
}
 
for(i = 3; i < na-1; i+=4){
	dyda[i] = cos(a[i+1]*x[0]);
	*y += a[i]*dyda[i];
	if(dydx_flag[0] || fita[i+1]) temp = -a[i]*sin(a[i+1]*x[0]);
	dydx[0] += a[i+1]*temp;
	dyda[i+1] = x[0]*temp;
}

return 0;
}



/* a quadradic in x and y */
int fxyquad(x, a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[];
double a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double x0, x1, x02, x12;

x0 = x[0];
x1 = x[1];
x02 = x0*x0;
x12 = x1*x1;

*y = a[0] + a[1]*x0 + a[2]*x1 + a[3]*x02 + a[4]*x12 + a[5]*x1*x0;
if(dyda_flag){
	if(fita[0]) dyda[0] = 1;
	if(fita[1]) dyda[1] = x0;
	if(fita[2]) dyda[2] = x1;
	if(fita[3]) dyda[3] = x02;
	if(fita[4]) dyda[4] = x12;
	if(fita[5]) dyda[5] = x1*x0;
}
if(dydx_flag[0]) dydx[0] = a[1] + 2*x0*a[3] + x1*a[5];
if(dydx_flag[1]) dydx[1] = a[2] + 2*x1*a[4] + x0*a[5];
return 0;
}

/* a single lorenzian */
int florenz(x, a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[];
double a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double denom, del, del2;

del = x[0] - a[0];
del2 = del * del;
denom = 4.0*del2 + a[1];
*y = a[2]*a[1]/denom;
if(dyda_flag){
	if(fita[0]) dyda[0] = 8.0*del*(*y)/denom;
	if(fita[1]) dyda[1] = (*y)/a[1] - (*y)/denom;
	if(fita[2]) dyda[2] = (*y)/a[2];
}
if(dydx_flag[0]) dydx[0] = -8.0*del*(*y)/denom;
return 0;
}

/* sum of two lorenzians */
int florenz2(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double denom, del, del2,y1, y2;

del = x[0] - a[0];
del2 = del * del;
denom = 4.0*del2 + a[1];
 y1 = a[2]*a[1]/denom;
if(dyda_flag){
	if(fita[0]) dyda[0] = 8.0*del*y1/denom;
	if(fita[1]) dyda[1] = y1/a[1] - y1/denom;
	if(fita[1]) dyda[2] = y1/a[2];
}
if(dydx_flag[0]) dydx[0] =  -8.0*del*y1/denom;

del = x[0] - a[3];
del2 = del * del;
denom = 4.0*del2 + a[4];
y2 = a[4]*a[5]/denom;
if(dyda_flag){
	if(fita[3]) dyda[3] = 8.0*del*y2/denom;
	if(fita[4]) dyda[4] = y2/a[4] - y2/denom;
	if(fita[5]) dyda[5] = y2/a[5];
	if(fita[6]) dyda[6] = 1;
}
if(dydx_flag[0]) dydx[0] +=  -8.0*del*y2/denom;
*y = y1+y2 + a[6];
return 0;
}

/* a gaussian */
int fgauss(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double fac, ex, arg;
arg = (x[0]-a[0])/a[1];
ex = exp (-arg*arg);
fac = a[2]*ex*2.0*arg;
*y = a[2]*ex;

if(dyda_flag){
	if(fita[0]) dyda[0] = fac/a[1];
	if(fita[1]) dyda[1] = fac*arg/a[1];
	if(fita[2]) dyda[2] = ex;
}
if(dydx_flag[0]) dydx[0] = -fac/a[1];
return 0;
}

/* a gaussian plus a constant */
int fgaussc(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double fac, ex, arg;
arg = (x[0]-a[0])/a[1];
ex = exp (-arg*arg);
fac = a[2]*ex*2.0*arg;
*y = a[2]*ex + a[3];
if(dyda_flag){
	if(fita[0]) dyda[0] = fac/a[1];
	if(fita[1]) dyda[1] = fac*arg/a[1];
	if(fita[2]) dyda[2] = ex;
	if(fita[3]) dyda[3] = 1;
}
if(dydx_flag[0]) dydx[0] = -fac/a[1];
return 0;
}

/* sum of gaussians plus a constant */
int fgaussn(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double fac, ex, arg;
int i;
dydx[0] = 0;
*y = 0;
for(i = 0; i < na-1; i+=3){
	arg = (x[0]-a[i])/a[i+1];
	ex = exp (-arg*arg);
	fac = a[i+2]*ex*2.0*arg;
	*y += a[i+2]*ex;
	if(dyda_flag){
		if(fita[i]) dyda[i] = fac/a[i+1];
		if(fita[i+1]) dyda[i+1] = fac*arg/a[i+1];
		if(fita[i+2]) dyda[i+2] = ex;
	}
	if(dydx_flag[0]) dydx[0] += -fac/a[i+1];
}
dyda[na-1] = 1;
*y += a[na-1];
return 0;
}


int fline(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
dyda[0] = 1;
dyda[1] = x[0];
*y = a[0] + a[1]*x[0];
dydx[0] = a[1];
return 0;
}

int fpoly(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
int i;
double xn, ynum;
xn = 1;
ynum = 0;
dydx[0] = 0;
for( i = 0; i < na; i++){
	ynum += xn*a[i];
	if(dydx_flag[0] && i > 0) dydx[0] += a[i-1] * xn;
	if(fita[i]) dyda[i] = xn;
	xn *= x[0];
}
*y = ynum;
return 0;
}

/* sum of exponentials plus a constant */
int fexpn(x,a,y,dyda,na,dyda_flag, fita, dydx_flag, dydx, ydat)
double x[],a[],*y,dyda[];
int na;
int dyda_flag;
int fita[];
int dydx_flag[];
double dydx[], ydat;
{
double fac, ex, arg;
int i;
double s;

*y = 0;
dydx[0] = 0;
for(i = 0; i < na-1; i+=3){
	 if(x[0] < a[i]) s = 0;
        else s = 1;
	arg = (x[0]-a[i])/a[i+1];
	ex = exp (-arg);
	fac = a[i+2]*ex;
	*y += s*fac;
	if(dyda_flag){
	if(fita[i]) dyda[i] = s*fac/a[i+1];
	if(fita[i+1]) dyda[i+1] = s*fac*arg/a[i+1];
	if(fita[i+2]) dyda[i+2] = s*ex;
	}
	if(dydx_flag[0]) dydx[0] += -s*fac/a[i+1];
}
dyda[na-1] = 1;
*y += a[na-1];
return 0;
}

