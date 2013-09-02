#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fit.h"

extern int debug;

/* function to get data and return number of data points */
int get_data(data, command, num_indep, inbuf, order,filename,maxrows, maxcols)
double **data;
char *command;     /* command which was typed at the fit> prompt */
int num_indep;
char *inbuf;           /* inbuf is really output buffer */
struct data_order *order;
char *filename;
int maxrows;
int maxcols;
{
int i, j, jmax, ndata;
char str[20][30];
FILE *instream;
int iopt1;
char  inbuf1[512];
iopt1 = 20;

sscanf(command,"%s %s %d", inbuf, filename, &iopt1);

if(( instream=fopen(filename,"r") ) == NULL) return 0;
i = 0;
while( fgets(inbuf1,512,instream) != NULL && i < maxrows){
	if(inbuf1[0] == '#') continue;
	jmax=sscanf(inbuf1,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
			str[0],str[1],str[2],str[3],str[4],str[5],str[6],str[7],str[8],str[9],
			str[10],str[11],str[12],str[13],str[14],str[15],str[16],str[17],str[18],str[19]);
	for(j=0; j < jmax && j <= iopt1 && j < maxcols-3; j++){
		data[j][i] = atof(str[j]);
	}
	i++;
}
fclose(instream);
ndata = i;

/* assign first columns to be independent variables */
for(i=0; i < num_indep; i++)order->x[i] = i;

/* assign next column to be dependent variables */
order->y = num_indep;

/* if there is a column num_indep+1, assign it to be error in y */
if(jmax > num_indep) order->isig = num_indep + 1;

/* assign a data column to be all ones: no weighting */
for(i = 0; i < ndata; i++)	data[jmax][i] = 1;
order->nsig = jmax+0;

/* assign a data column to be used for statistical weighting */
order->ssig= jmax+2;

/* assign a data column to be used for yfit */
order->yfit= jmax+1;

if(jmax == 2 && num_indep == 1) order->sig = 2;
if(jmax == 3 && num_indep == 2) order->sig = 3;

/* Assigning data columns for different weightings */
/* takes up extra memory, but is faster and simpler */
/* especially in the case of statistical weighting: */
/* we don't have to take the sqrt of y at every point */
/* for every iteration.  It is done once when statistical */
/* weighting is selected */

return ndata;
}

/* This function initializes the parameters */
void ip(command, inbuf, func, a, ma)
char *command, *inbuf;
int (*func)();
double *a;
int ma;
{
int i;
char cmdargs[NUM_ARGS][30];

if(func == NULL)
	printf("ip failed, you must select a function first\n");
else if(parse(command,cmdargs) == ma){
	for( i = 0; i < ma; i++){
		/* a * in place of a number leaves the parameter unchanged */
		if(strcmp(cmdargs[i],"*") != 0) a[i] = atof(cmdargs[i]);
	}
}
else
	printf("ip failed, there should be %d parameter values entered\n",ma);
}

/* This function changes one of the parameters */
void cp(command, inbuf, func, a, ma)
char *command, *inbuf;
int (*func)();
double *a;
int ma;
{
int i;
char cmdargs[NUM_ARGS][30];

if(func == NULL)
	printf("cp failed, you must select a function first\n");
else if(parse(command,cmdargs) == 2 && atoi(cmdargs[0]) > -1
		&& atoi(cmdargs[0]) < ma)
			a[atoi(cmdargs[0])] = atof(cmdargs[1]);
else
	help("cp",1);
}

/* sp selects the parameters to vary */
int sp(command, inbuf, func, lista, ma)
char *command, *inbuf;
int (*func)();
int *lista, ma;
{
int i, mfit, max=0, min = 0;
char cmdargs[NUM_ARGS][30];

if(func == NULL){
	printf("sp failed, you must select a function first\n");
	return 0;
}
else{
	mfit = parse(command,cmdargs);
	if(debug) printf("mfit: %d ma: %d\n", mfit, ma);
	/* if mfit > ma, default to fitting all parameters */
	if(ma < mfit){
		mfit = ma;
		for(i=0; i < ma; i++) lista[i] = i;
		printf("sp failed, too many parameters selected\n");
		if(debug)printf("mfit: %d ma: %d\n", mfit, ma);
	}
	else{
		/* assign lista[i] to command arguments */
		for( i = 0; i < mfit; i++) lista[i] = atoi(cmdargs[i]);
		/* if any lista[i]'s are too big or too small, Default to all parameters */
		for(i = 0; i < mfit; i++){
			if( lista[i] > max) max = lista[i];
			if( lista[i] < min) min = lista[i];
		}
		if(max > ma-1 || min < 0){
			printf("sp failed, parameter out of range, all parameters fitted\n");
			mfit = ma;
			for(i=0; i < ma; i++) lista[i] = i;
		}
	}
	}
return mfit;
}

int make_data(int (*func)(), int num_indep, double *a, int ma, char *command, char *inbuf){

char cmdargs[NUM_ARGS][30];
char filename[80];
FILE *stream;
double *x, *xmin, *xmax, *xstep, y;
int i,j, *fita;
double  *dydx, *dyda;
int *dydx_flag;
int failed;
double ydat = 0;

/* parse and check number of args entered on command line */
if( (parse(command,cmdargs)) != 1 + 3*num_indep || func == NULL){
	help("md",1);
	return 1;
}

if( (stream = fopen(cmdargs[0],"w")) == NULL){
	printf("md error, cannot open file %s\n", cmdargs[0]);
	return 1;
}

if( num_indep > 2 ){
	printf("md error, md not implemented for more then 2 independent parameters\n");
	return 1;
}

/* need to send all these to the function */
dydx = dvector(num_indep);
dydx_flag = ivector(num_indep);
fita = ivector(ma);
dyda = dvector(ma);

/* we only need function value, not derivatives */
for(j = 0; j < num_indep; j++) dydx_flag[j] = 0;
for(j = 0; j < ma; j++) fita[j] = 0;

x = dvector(num_indep);
xmin = dvector(num_indep);
xmax = dvector(num_indep);
xstep = dvector(num_indep);

/* assign command line arguments */
for( j = 0; j < num_indep; j++){
	x[j] = xmin[j] = atof(cmdargs[3*j+1]);
	xmax[j] = atof(cmdargs[3*j+2]);
	xstep[j] = atof(cmdargs[3*j+3]);
	if(debug)printf("xmin: %g xmax: %g xstep: %g\n", xmin[j], xmax[j], xstep[j]);
}

if(num_indep == 1){
	while( x[0] <= xmax[0] ){
		failed = (*func)(x,a,&y,dyda,ma,0,fita,dydx_flag,dydx,ydat);
		if(failed){
			free(dydx);
			free(dydx_flag);
			free(fita);
			free(dyda);
			free(x);
			free(xmin);
			free(xmax);
			free(xstep);
		return 1;
		}
		fprintf(stream,"%g %g\n", x[0], y);
		x[0] += xstep[0];
	}
}

if(num_indep == 2){
	while( x[0] <= xmax[0] ){
		x[1] = xmin[1];
		while( x[1] <= xmax[1] ){
			failed = (*func)(x,a,&y,dyda,ma,0,fita,dydx_flag,dydx,ydat);
			if(failed){
				free(dydx);
				free(dydx_flag);
				free(fita);
				free(dyda);
				free(x);
				free(xmin);
				free(xmax);
				free(xstep);
				return 1;
			}
			fprintf(stream,"%g %g %g \n", x[0], x[1], y);
			if(debug)printf("%g %g %g\n", x[0], x[1], y);
			x[1] += xstep[1];
		}
		x[0] += xstep[0];
	}
}

fclose(stream);
strcpy(inbuf,command);

free(dydx);
free(dydx_flag);
free(fita);
free(dyda);
free(x);
free(xmin);
free(xmax);
free(xstep);
return 0;

}

void est_errors(factor,func, data, order, num_indep, a, ndata, ma, chisqr)
double factor;
int (*func)();
double **data;
struct data_order order;
int num_indep;
double a[];
int ndata, ma;
double chisqr;
{
int i, j, failed,k;
double tchisqr, tai, x1, x2, y1, y2, x;
double perr, merr;

for(i = 0; i < ma; i++){  /* loop over parameters */
	tai = a[i];              /* temporary a[i] */
	tchisqr = chisqr;        /* temporary chisqr */
	x = factor*chisqr;       /* chisqr we are looking for */

	/* change a[i] and use an iterative linear interpolation */
	a[i] += 1e-7*tai;
	failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &tchisqr);
	x2 = chisqr;
	x1 = tchisqr;
	y1 = a[i];
	y2 = tai;
	a[i] = ( (x - x2)*y1 - (x - x1)*y2 )/(x1-x2);

	k = 0;
	while( fabs(tchisqr - factor*chisqr) > 1e-5*fabs(chisqr) && k < 1000){
	if(debug)printf("tai: %g a[%d]: %g chisqr: %g tchisqr: %g\n", tai, i, a[i], chisqr, tchisqr);
	x2 = x1;
	y2 = y1;
	failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &tchisqr);
	x1 = tchisqr;
	y1 = a[i];
	a[i] = ( (x - x2)*y1 - (x - x1)*y2 )/(x1-x2);
	k++;
	}
	perr = a[i] - tai;
	a[i] = tai;

	/* repeat the process for negative error, except begin with
	assuming negative error is the same as the positive error */

	tchisqr = chisqr;
	a[i] -= perr;
	failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &tchisqr);
	x2 = chisqr;
	x1 = tchisqr;
	y1 = a[i];
	y2 = tai;
	a[i] = ( (x - x2)*y1 - (x - x1)*y2 )/(x1-x2);

	k = 0;
	while( fabs(tchisqr - factor*chisqr) > 1e-5*fabs(chisqr) && k < 1000){
	if(debug)printf("tai: %g a[%d]: %g chisqr: %g tchisqr: %g\n", tai, i, a[i], chisqr, tchisqr);
	x2 = x1;
	y2 = y1;
	failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &tchisqr);
	x1 = tchisqr;
	y1 = a[i];
	a[i] = ( (x - x2)*y1 - (x - x1)*y2 )/(x1-x2);
	k++;
	}

	merr = -a[i] + tai;
	a[i] = tai;

	printf("a[%d] = %g (+ %g, -%g)\n", i, a[i], fabs(perr), fabs(merr));
}

}
