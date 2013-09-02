#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fit.h"

/* This file contains a bunch of function for allocating */
/* vectors and matrices.  */

extern int debug;
extern int wiflag;
void myerror(char *);

double **dmatrix(nr,nc)
int nr,nc;
{
	int i;
	double **m;
	char s[80];

	m=(double **) malloc((unsigned) (nr)*sizeof(double*));
	if (!m) myerror("allocation failure 1 in dmatrix()");

	for(i=0;i<nr;i++) {
		m[i]=(double *) malloc((unsigned) (nc)*sizeof(double));
		sprintf(s,"allocation failure 2 in dmatrix() on column %d",i);
		if (!m[i]) myerror(s);
	}
	return m;
}

int **imatrix(nr,nc)
int nr,nc;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nr)*sizeof(int*));
	if (!m) myerror("allocation failure 1 in imatrix()");

	for(i=0;i<nr;i++) {
		m[i]=(int *)malloc((unsigned) (nc)*sizeof(int));
		if (!m[i]) myerror("allocation failure 2 in imatrix()");
	}
	return m;
}

void free_dmatrix(m,nr,nc)
double **m;
int nr,nc;
{
	int i;

	for(i=nr-1;i>=0;i--) free((char*) (m[i]));
	free((char*) (m));
}

void free_imatrix(m,nr,nc)
int **m;
int nr,nc;
{
	int i;

	for(i=nr-1;i>=0;i--) free((char*) m[i]);
	free((char*) (m));
}

void myerror(s)
char s[80];
{
printf("%s\n", s);
}

double *dvector(int n){
return (double *)malloc((unsigned) (n)*sizeof(double));
}

int *ivector(int n){
return (int *)malloc((unsigned) (n)*sizeof(int));
}


/* parse parses command line arguments */

int parse(command, cmdarg)
char *command, cmdarg[NUM_ARGS][30];
{
int i=0, j=0, k, toomany=0;

while( command[i] != '\x0' ){
	if(j >= NUM_ARGS){
		toomany = 1;
		break;
	}
	if( command[i] == ' ' ){
		i++;
		k = 0;
		while(command[i] == ' ') i++;
		while( command[i] != ' ' && k < 30 && command[i] != '\x0'){
			if(command[i] != ' '){
				cmdarg[j][k] = command[i];
				k++;
			}
		i++;
	   }
	cmdarg[j][k] = '\x0';
	j++;
	}
if(j==0)i++;
}
if(toomany){ 
	printf("Warning: there is a maximum of %d command arguments\n",NUM_ARGS);
	printf("Change NUM_ARGS in fit.h and recompile\n");
}
if(k == 0) j--;
if(debug){
	printf("#args: %d\n",j);
	for( i = 0; i < j; i++) printf("%s\n", cmdarg[i]);
}
return j;
}

/* calc_yfit calculates the value of the */
/* fitting function for all the values of */
/* the independent variable */

int calc_yfit(func, data, order, num_indep, a, ndata, ma, chisqr)
int (*func)();
double **data;
struct data_order order;
int num_indep;
double a[];
int ndata, ma;
double *chisqr;
{
int i, j;
int failed;
double *dyda;
double *dydx;
int *fita;
double *y;
double *x;
int *dydx_flag;
double sig2i;

y = data[order.yfit];

dyda = dvector(ma);
fita = ivector(ma);
x = dvector(num_indep);
dydx = dvector(num_indep);
dydx_flag = ivector(num_indep);

for(j = 0; j < num_indep; j++) dydx_flag[j] = 0;
for(j = 0; j < ma; j++) fita[j] = 0;

*chisqr = 0;
for(i=0; i < ndata; i++){
	for(j = 0; j < num_indep; j++){
		x[j] = data[order.x[j]][i];
		if(debug > 1)
			printf("i: %d order.x[%d]: %d x[%d]: %g\n", i, j, order.x[j], j, x[j]);
	}
	if(debug > 1) printf("about to call (*func)()\n");
	failed = (*func)(x,a,&y[i],dyda,ma,0,fita,dydx_flag,dydx,data[order.y][i]);
	if(failed) break;
	sig2i = data[order.sig][i]*data[order.sig][i];
	if(sig2i > 1e-30 && (wiflag == 0 || window(x, num_indep)))
		*chisqr += (y[i] - data[order.y][i])*(y[i] - data[order.y][i])/sig2i;
}
free(dyda);
free(fita);
free(x);
free(dydx);
free(dydx_flag);
return failed;
}


/* help lists comands */
int help(char topic[12], int maxlines){
char inbuf[80];
FILE *stream;
int i;
char buf[20];
if(strcmp(topic,"") == 0)strcpy(topic,"COMMANDS");

stream = NULL;
if(getenv("FITHELP") == NULL && (stream = fopen("fit.hlp","r")) == NULL){
	printf("help file not found, set environment variable FITHELP to point to fit.hlp\n");
	return 1;
}

if( (stream == NULL) && (stream = fopen(getenv("FITHELP"),"r")) == NULL){
	printf("help file not found, set environment variable FITHELP to point to fit.hlp\n");
	return 1;
}

while(fgets(inbuf,80,stream) != NULL){
	if( strncmp(inbuf, topic, strlen(topic)) == 0){
		i = 0;
		printf("%s",inbuf);
		while( strncmp( fgets(inbuf,80,stream), "END",3) ){
			printf("%s",inbuf);
			i++;
			if(i == 20){ i = 0; gets(buf); }
			if( i == maxlines) return 0;
		}
	}
}
fclose(stream);
return 0;
}

int nonzero(double *array, int ndata){
int i;
for(i = 0; i < ndata; i++) if(array[i] < 1e-60) return 0;
return 1;
}

