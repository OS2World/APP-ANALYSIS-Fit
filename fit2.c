/*** fit2.c last modified 17 AUG 1993 ***/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fit.h"

#ifdef ICC
FILE *popen(char *, char *);
void pclose(FILE *);
void sleep(int);
#endif

/* GLOBAL DATA */

int grflag = 1;       /* do we graph data and fit? */
int gnuopen = 0;      /* is there a stream to gnuplot open */
FILE *pstream;        /* stream to gnuplot */
FILE *pstream2;        /* stream to fitplot */
double *xmin, *xmax;  /* min and max independent variables used in windowing */
int wiflag = 0;       /* windowing flag */
int veflag = 1;       /* verbosity flag */
int debug = 0;        /* debugging flag */
char GNUPLOT[60];     /* gnuplot command line */
int MYPLOT = 0;       /* flag for using fitplot */
int doflag = 0;       /* flag for command file processing */
int noset = 0;        /* flag so we do not set any gnuplot labels or titles */

/* static data */

static double *a;            /* array of parameters for fit */
static double **covar;       /* covariant matrix for fit */
static double chisq;         /* squared error for fit */
static double **data;        /* matrix which holds data */
static int i, j, k;          /* indices for loops */
static int ndata=0;          /* number of data points */
static int ma;               /* number of parameters */
static int *lista;           /* array which has the list of parameters to vary */
static int mfit;             /* number of parameters being varied */
static int  itmax=20;        /* number of iterations */
static int iopt1, iopt2;     /* integers used in command options */
static char file[30];        /* name of file for reading or writing parameters*/
static char inbuf[80]="";    /* buffers used for input and output of strings */
static char buf[25];
static char fname[20];        /* function name */
static char cmdargs[NUM_ARGS][30];  /* command line arguments */
static char gnubuf[COMMAND_SIZE];      /* command which is sent to gnuplot */
static char syscmd[COMMAND_SIZE];      /* command which is sent to system */
static char filename[80];     /* name of data file */
static char topic[20];        /* help topic */
static int *co_lista;         /* lista for writing covariant matrix */
static double **alpha;        /* matrix and vectors used in LM algorithm */
static double *beta, *da;
static int num_indep = 1;     /* number of independent variables */
static int last_ma;           /* ma for last function */
static int last_num_indep;    /* num_indep for last function */
static int linflag = 0;       /* is function linear */
static int failed;            /* did function fail? */
/* each function has a comment, this holds the comment of the */
/* current function.  If function is plotted using gnuplot, comment */
/* should be f(x) = .... , where the comment defines the function */
/* as it would on the gnuplot command line */
static char comment[160];

/* mode in which file is opened for reading and writing parameters */
char mode[5];

/* see description of data_order in fit.h */
/* note that the order is the assignment of the data columns */
/* internally. */
static struct data_order order;
static FILE *stream;         /* stream for reading and writing parameters to a file */
static int datarows = 1024;  /* number of rows in the data matrix */
static int datacols = 6;     /* number of columns in the data matrix */
static int echo = 0;

int (*func)(double *x, double *a, double *y, double *dyda,
				int na, int dyda_flag, int *fita, int dydx_flag, 
				double *dydx, double ydat);

void process_command(char *);

main(int argc, char **argv){
char command[COMMAND_SIZE]="";   /* what was typed at command line */
/* assign these pointers to point to NULL */
a=NULL;
covar=NULL;
lista=NULL;
func=NULL;
pstream = NULL;
pstream2 = NULL;

help("NOTE",100);

/* assign the proper command string to GNUPLOT */

strcpy(GNUPLOT,"gnuplot");

#ifdef OS2 
strcpy(GNUPLOT,"gnuplot 2> gnuout");
if(argc > 1){
	MYPLOT = 1;
	fitcmd("doiky\n");
}
#endif

#ifdef DOS
grflag = 0;
#endif

/* under UNIX, this rather complicated command string is needed. */
/* If modified, you should comment this one out and make another */
/* so you can go back to this one if needed. */
/* Things work pretty well with just a "gnuplot" command, */
/* except that some of the things gnuplot writes to stderr */
/* make it to the user's screen.  Simply using "gnuplot >& gnuout" */
/* does not work, because popen() uses sh commands, and  */
/* "gnuplot >& gnuout" uses csh redirection. */

#ifdef UNIX 
strcpy(GNUPLOT,"/bin/csh -f -c \"gnuplot >& /tmp/gnuout\""); 
/* strcpy(GNUPLOT,"gnuplot"); */
if(debug)printf("*%s*\n", GNUPLOT);
#endif

/* set up for one independent variable by default */
order.x = ivector(num_indep);
order.x[0] = 0;
order.xsig = ivector(num_indep);
for(i = 0; i < num_indep; i++) order.xsig[i] = -1;

xmin = dvector(num_indep);
xmax = dvector(num_indep);

last_num_indep = num_indep;

last_ma = ma = 0;

/* allocate space for data matrix */
data=dmatrix(datacols,datarows);

/* command processing loop */

while(strncmp(command,"quit",4) != 0){
	strcpy(command,"\x0\x0\x0\x0\x0\x0\x0\x0");
	printf("fit2> ");
	fflush(stdout);
	gets(command);
	process_command(command);
}

/* free allocated arrays and close open pipe to gnuplot */
free_dmatrix(data,datacols,datarows);
if(a != NULL) free(a);
if(covar != NULL) free_dmatrix(covar,ma,ma);
if(lista != NULL) free(lista);
free(order.x);
free(order.xsig);
free(xmin);
free(xmax);
#ifndef DOS
if(pstream != NULL){
	sprintf(gnubuf,"quit\n");
	gnucmd(gnubuf);
	pclose(pstream);
}
#endif
#ifdef OS2
if(pstream2 != NULL){
	sprintf(gnubuf,"quit\n");
	fitcmd(gnubuf);
	pclose(pstream2);
}
#endif
#ifdef UNIX
if(!debug){
	system("rm /tmp/fit.tmp");
	system("rm /tmp/gnuout");
}
#endif
printf("\n");
return 0;
}

int gnucmd(char *command){

#ifndef DOS

grflag = 1;
/* open the pipe */
if(pstream == NULL){
	if(debug) printf("attempting to open pipe to gnuplot\n");
	pstream = popen(GNUPLOT,"w");
}       
if(pstream == NULL){
	printf("Failed to open pipe to gnuplot\n");
	return 1;
}
fprintf(pstream, command);
if(debug) printf("command to gnuplot: %s\n", command);
fflush(pstream);
#endif

#ifdef DOS
printf("The DOS version of this program does not support gnuplot\n");
printf("directly.  You need a real operating system for that.  \n");
printf("I suggest OS/2 2.1. You can use wf to write the fit  \n");
printf("to a file and use gnuplot by itself for plotting  \n\n");
grflag = 0;
#endif

return 0;
}

#ifdef OS2

int fitcmd(char *command){

grflag = 1;
/* open the pipe */
if(pstream2 == NULL){
	if(debug) printf("attempting to open pipe to fitplot\n");
	pstream2 = popen("fitplot.cmd > nul","w");
}       
if(pstream2 == NULL){
	printf("Failed to open pipe to fitplot\n");
	return 1;
}
fprintf(pstream2, command);
if(debug) printf("command to fitplot: %s\n", command);
fflush(pstream2);

return 0;
}
#endif

void process_command(char *command){

char cmd[COMMAND_SIZE];
FILE *dostream;

/* if there is no command, reprompt */
if(strcmp(command,"")==0) return;

if(echo) printf("%s\n",command);

/* gd stands for get data */
if(strncmp(command,"gd",2) == 0){
	if(debug) printf("in main, calling get_data()\n");
	ndata=get_data(data,command,num_indep,inbuf,&order,filename,
					datarows, datacols);
	if(debug) printf("in main, get_data() returned %d\n", ndata);
	if(ndata ==0) printf("NO DATA, check filename\n");
	printf("gd: %d data points\n", ndata);
}

/* do executes commands from a file */
else if(strncmp(command,"do",2) == 0){
	doflag = 1;
	if(parse(command,cmdargs) < 1) help("do", 1);
	else{
		if((dostream = fopen(cmdargs[0],"r")) == NULL){
			printf("cannot open file %s\n", cmdargs[0]);
			return;
		}
	while(fgets(cmd, COMMAND_SIZE, dostream) != NULL)
		if(cmd[0] != '#') process_command(cmd);
	}
	fclose(dostream);
	doflag = 0;
}

/* md stands for make data */
else if(strncmp(command,"md",2) == 0){
	if(debug) printf("in main, calling make_data()\n");
	failed=make_data(func, num_indep, a, ma, command,inbuf);
	if(debug) printf("in main, make_data() returned %d\n", failed);
}

/* sh stands for show */
else if(strncmp(command,"sh",2) == 0){
	for(i = 0;  i < num_indep; i++)printf("order.x[%d]: %d\n", i, order.x[i]);
	printf("order.y: %d\n", order.y);
	printf("order.yfit: %d\n", order.yfit);
	printf("order.sig: %d\n", order.sig);
	printf("order.nsig: %d\n", order.nsig);
	printf("order.ssig: %d\n", order.ssig);
	printf("order.isig: %d\n", order.isig);
	printf("order.osig: %d\n", order.osig);
	for(i = 0;  i < num_indep; i++)printf("order.xsig[%d]: %d\n", i, order.xsig[i]);
	if(func != NULL){
		printf("function: %s %s\n", fname, comment);
		printf("varying parameters: ");
		for(i = 0; i < mfit; i++) printf(" %d",lista[i]);
	printf("\n");
	}
	if(ndata) printf("data file: %s # data points: %d\n", filename, ndata);
	if(wiflag){
		printf("\n windowing on. ");
		for(i = 0; i < num_indep; i++)
		printf("xmin%d: %g xmax%d: %g\n", i, xmin[i], i, xmax[i]);
	printf("\n");
	}
	else printf("windowing off\n");
	if(grflag)printf("graphing turned on\n");
	else printf("graphing turned off\n");
}

/* fn command tells us which function to fit to */
else if(strncmp(command,"fn",2) == 0){
	last_ma = ma;
	last_num_indep = num_indep;
	iopt2 = sscanf(command,"%s %s %d ", inbuf, fname, &iopt1);
	if(debug)printf("in main, calling getfcnptr()\n");
	if((func = (int *)getfcnptr(fname,&num_indep, &linflag,&ma, comment)) != NULL){
		if(debug)printf("in main, returned from getfcnptr()\n");
		/* get ma from command line or prompt for it if needed */
		if( ma == -1){
			if(iopt2 > 2) ma = iopt1;
			else{
				printf("\n Enter number of parameters: ");
				fflush(stdout);
				gets(inbuf);
				sscanf(inbuf,"%d",&ma);
			}
		}
		if(debug)("ma: %d last_ma: %d", ma, last_ma);
	
		/* if covar, a, and lista are allocated and need to change size, we free them */
		if(covar != NULL && ma != last_ma){
			if(debug)printf("free_dmatrix(covar,ma,ma)\n");
			free_dmatrix(covar,last_ma,last_ma);
		}
		if(lista !=NULL && ma != last_ma){
			if(debug)printf("free(lista)\n");
			free(lista);
		}
		if(a != NULL && ma != last_ma){
			if(debug)printf("free(a)\n");
			free(a);
		}

		/* if num_indep changed, we reallocate some stuff */
		if(num_indep != last_num_indep){
			free(order.x);
			free(order.xsig);
			free(xmin);
			free(xmax);
			order.x = ivector(num_indep);
			order.xsig = ivector(num_indep);
			xmin = dvector(num_indep);
			xmax = dvector(num_indep);

			/* assign default values */
			for(i = 0; i < num_indep; i++){
			order.xsig[i] = -1;
			order.x[i] = i;
			}
		}
	
		/* allocate space for covar, lista, and a */
		if(ma != last_ma){
			covar=dmatrix(ma,ma);
			lista=ivector(ma);
			a=dvector(ma);
			/* initialize a's to 5, default to fitting all of the parameters */
			for(i=0;i<ma;i++){
				lista[i] = i;
				a[i] = 5;
			}
			mfit=ma;
			printf("%s: %d parameters\n", command, ma);
		}
	}
	else printf("Function %s not found, try lf to list functions\n", fname);
}

/* fit command does the fit */
else if(strncmp(command,"fit",2) == 0){
	if(func == NULL)printf("fi failed, choose function first\n");
	else if(ndata == 0)printf("fi failed, no data\n");
	else{
		sscanf(command,"%s %d ", inbuf, &itmax);
		/* nonzero returns 1 if all sigma's are non-zero */
		if(nonzero(data[order.sig], ndata)){
			if(debug)printf("in main, calling mrqfit(), itmax: %d\n", itmax);
			failed = mrqfit(data,order,num_indep,ndata,itmax,a,ma,
				lista,mfit,covar,&chisq,func,filename,comment);
			if(debug)printf("in main, mrqfit() returned %d\n", failed);
		}
		else{
			failed = 1;
			printf("all sigmay's must be non-zero, check weighting and order\n");
		}
		if(failed == 0)printf("fit done\n");
		if(failed != 0)printf("fit failed\n");
		printf("%s\n",comment);
	}
}


/* ip command initializes the parameters */
else if(strncmp(command,"ip",2) == 0){
	if(debug)printf("in main, calling ip()\n");
	ip(command, inbuf, func, a, ma);
	if(debug)printf("in main, returned from ip()\n");
}

/* cp command changes one of the parameters */
else if(strncmp(command,"cp",2) == 0){
	if(debug)printf("in main, calling ip()\n");
	cp(command, inbuf, func, a, ma);
	if(debug)printf("in main, returned from ip()\n");
}

/* sp command selects the parameters */
else if(strncmp(command,"sp",2) == 0){
	if(debug)printf("in main, calling sp()\n");
	mfit = sp(command, inbuf, func, lista, ma);
	if(debug)printf("in main, sp() returned %d\n", mfit);
}

/* pp command prints the parameters */
else if(strncmp(command,"pp",2) == 0){
	if(func == NULL)printf("pp failed, you must select a function first\n");
	else{
		printf("%s\n",comment);
		for(i = 0; i < ma; i++) printf("a%d= %15.9g\n", i, a[i]);
		printf("chisqr = %g\n", chisq);
		printf("%s\n",command);
	}
}

/* wp command writes the parameters to a file */
else if(strncmp(command,"wp",2) == 0){
	if(func == NULL)printf("wp failed, you must select a function first\n");
	else{
		strcpy(file,"a.dat");  /* filename defaults to a.dat */
		if(sscanf(command,"%s %s %s", inbuf, file,mode) < 3)
		strcpy(mode,"w");  /* mode defaults to overwriting */
		stream = fopen(file,mode);
		if(stream){
			for(i=0; i < ma; i++)
			fprintf(stream,"%19.14e\n", a[i]);
			fclose(stream);
		}
		else printf("cannot open file %s in mode %s\n", file, mode);
		printf("%s\n",command);
	}
}

/* wf writes the fitting function data to a file */
else if(strncmp(command,"wf",2) == 0){
	if(func == NULL)printf("wf failed, you must select a function first\n");
	else{
		if(debug)printf("in main, calling calc_yfit()\n");
		failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisq);
		if(debug)printf("in main, calc_yfit() returned %d\n", failed);
		if(!failed){
			strcpy(file,"fit.dat");  /* filename defaults to fit.dat */
			sscanf(command,"%s %s", inbuf, file);
			stream = fopen(file,"w");
			if(stream){
				for(i=0; i < ndata; i++){
					for(j = 0; j < num_indep; j++)
					fprintf(stream,"%g ", data[order.x[j]][i]);
					fprintf(stream,"%g\n", data[order.yfit][i]);
				}
				fclose(stream);
			}
			else printf("cannot open file %s\n", file);
		}
		else printf("%s failed \n",command);
	}
}

/* rp reads parameters from a file */
else if(strncmp(command,"rp",2) == 0){
	if(func == NULL)printf("rp failed, you must select a function first\n");
	else{
		strcpy(file,"a.dat");  /* filename defaults to a.dat */
		sscanf(command,"%s %s", inbuf, file);
		stream = fopen(file,"r");
		i = 0;
		if(stream){
			while(fgets(inbuf,28,stream) && i < ma){
				a[i] = atof(inbuf);
				i++;
			}
			fclose(stream);
			printf("%s: %d parameters read in\n", command, i);
		}
		else printf("cannot open file %s\n", file);
	}
}

/* plot opens a pipe to gnuplot and tells it to plot the data and fit */
else if(strncmp(command,"plot",2) == 0){
	if(func == NULL) printf("plot failed, you must select a function first\n");
	else{
		grflag = 1;
		if(MYPLOT){
#ifdef OS2
			if(debug) printf("in main, about to call calc_yfit()\n");
			failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisq);
			if(debug) printf("in main, returned from calc_yfit()\n");
			if(debug)printf("in main, calling myplot()\n");
			failed=myplot(func, data, order, num_indep, ndata,
				filename, comment, a, ma);
			if(debug)printf("in main, myplot() returned %d\n", failed);
#endif
		}
		else{
			if(debug)printf("in main, calling plot()\n");
			failed=plot(func, data, order, num_indep, ndata,
			filename, comment, a, ma);
			if(debug)printf("in main, plot() returned %d\n", failed);
		}
		printf("%s\n",command);
		if(failed)printf("%s failed\n",command);
	}
}

else if(strncmp(command,"help",1) == 0){
	strcpy(topic,"");
	sscanf(command,"%*s %s", topic);
	if(debug)printf("in main, calling help()\n");
	help(topic,100);
	if(debug)printf("in main, returned from help()\n");
}

/* pr opens a pipe to gnuplot and tells it to plot the residual error */
/* "pr 1" plots bestfit vs y, "pr 2" plots (bestfit-y) vs x */
else if(strncmp(command,"pr",2) == 0){
	if(func == NULL)printf("pr failed, you must select a function first\n");
	else{
		grflag = 1;
		if(debug)printf("in main, calling pr()\n");
		failed=pr(command,func, data, order, num_indep, ndata,
			filename, comment, a, ma);
		if(debug)printf("in main, pr() returned %d\n", failed);
		if(failed) printf("%s failed\n",command);
	}
}


else if(strncmp(command,"",1) == 0){
	strcpy(topic,"");
	sscanf(command,"%*s %s", topic);
	if(debug)printf("in main, calling help()\n");
	help(topic,100);
	if(debug)printf("in main, returned from help()\n");
}

/* send commands beginning with set or load to gnuplot */
	else if(strncmp(command,"set",2) == 0 ||  strncmp(command,"lo",2) == 0){
	sprintf(gnubuf,"%s\n",command);
	gnucmd(gnubuf);
}

/* wt select weighting for fit */
else if(strncmp(command,"wt",2) == 0){
	if(ndata == 0) printf("wt failed, you must get data first\n");
	else{
		sscanf(command,"%*s %s", inbuf);
		if(strncmp(inbuf,"statistical",1) == 0){
			order.sig = order.ssig;
			for(i=0; i<ndata; i++){
				if(data[order.y][i] > 1e-30)
					data[order.sig][i] = sqrt(fabs(data[order.y][i]));
				else data[order.sig][i] = 1;
			}
		}
		else if(strncmp(inbuf,"instrumental",1) == 0)
			order.sig = order.isig;
		else if(strncmp(inbuf,"none",1) == 0)
			order.sig = order.nsig;
		else if(strncmp(inbuf,"other",1) == 0)
			order.sig = order.y;
		else help("wt",1);
	}
}

/* or redefines order of colums in data */
else if(strncmp(command,"or",2) == 0){
	if(ndata == 0)printf("or failed, you must get data first\n");
	else{
		i = parse(command, cmdargs);
		if(debug) printf("num_indep: %d num args: %d\n", num_indep, i);
		for(j = 0; j < num_indep; j++) if(i > j){
			if(debug) printf("cmdargs[%d]: %s\n",j, cmdargs[j]);
			order.x[j] = atoi(cmdargs[j]);
			if(debug) printf("order.x[%d]: %d\n", j, order.x[j]);
		}
		if(i > num_indep) order.y = atoi(cmdargs[num_indep]);
		if(i > num_indep + 1) order.isig = atoi(cmdargs[num_indep + 1]);
		for(j = 0; j < num_indep; j++)
		if(i > j+num_indep+2) order.xsig[j] = atoi(cmdargs[j+num_indep+2]);
		if(i == 0) help("or",1);
	}
}

/* co prints the covariant matrix from the fit */
/* if there are no arguments, output is to */
/* screen.  An argument is taken as a filename */
/* and matrix is written to a file */
else if(strncmp(command,"co",2) == 0){
	if(func == NULL)printf("co failed, you must select a function first\n");
	else{
		strcpy(file,"");
		sscanf(command,"%s %s", inbuf, file);
		co_lista = ivector(ma);
		beta = dvector(ma);
		da = dvector(ma);
		alpha=dmatrix(ma,ma);
		for(i = 0; i < ma;  i++) co_lista[i] = i;
		if(debug)printf("in main, calling alpha_beta_chisqr()\n");
		alpha_beta_chisq(data, order, num_indep, ndata,a,ma,
			co_lista,ma,alpha,beta,&chisq,func);
		if(debug)printf("in main, returned from alpha_beta_chisqr()\n");
		if(debug)printf("in main, calling solve_for_da()\n");
		solve_for_da(alpha, covar, beta, da, ma);
		if(debug)printf("in main, returned from solve_for_da()\n");
		free(co_lista);
		free(beta);
		free(da);
		free_dmatrix(alpha,ma,ma);
		if( strcmp(file,"") == 0){
			for(i = 0; i < ma; i++){
				strcpy(inbuf,"");
				for(j = 0; j < ma; j++){
					sprintf(buf,"%9.2g ",covar[i][j]);
					strcat(inbuf, buf);
				}
				printf("%s\n",inbuf);
			}
		}
		else{
			stream = fopen(file,"w");
			for(i = 0; i < ma; i++){
				strcpy(inbuf,"");
				for(j = 0; j < ma; j++){
					sprintf(buf,"%10g ",covar[i][j]);
					strcat(inbuf, buf);
				}
				fprintf(stream,"%s\n",inbuf);
			}
			fclose(stream);
		}
	}
}

/* ad allocates data matrix to a specified size */
else if(strncmp(command,"ad",2) == 0){
	iopt1 = 0;
	iopt2 = 0;
	sscanf(command,"%s %d %d", inbuf, &iopt1, &iopt2);
	if(debug)printf("cols: %d rows: %d\n", iopt1, iopt2);
	if(iopt1 > 3 && iopt2 > 1 && iopt1 < 128){
		if(debug) printf("calling free_dmatrix(data,datacols,datarows);\n");
		free_dmatrix(data,datacols,datarows);
		datacols = iopt1;
		datarows = iopt2;
		if(debug) printf("calling data=dmatrix(data,datacols,datarows);\n");
		data=dmatrix(datacols,datarows);
		if(debug) printf("returned from data=dmatrix(data,datacols,datarows);\n");
		ndata = 0;
	}
else printf("ad: columns rows, columns => columns in data file + 3\n");
}

/* gr turns graphing on and off */
else if(strncmp(command,"gr",2) == 0){
	sscanf(command,"%*s %d", &iopt1);
	if(iopt1 == 1|| iopt1 == 0){
		grflag =iopt1;
	}
	else printf("must be \"gr 0\" (graphing off) or \"gr 1\" (graphing on)\n");
}

/* lf lists the functions available for fitting */
else if(strncmp(command,"lf",2) == 0){
	listfcns();
}

/* wi with one parameter turns windowing on and off */
/* with two parameters, it selects xmin and xmax */
else if(strncmp(command,"wi",2) == 0){
	iopt1 = parse(command,cmdargs);
	if(iopt1 < 1) printf("wi failed, try help wi\n");
	else if(iopt1 == 1){
		wiflag = atoi(cmdargs[0]);
	}
	else if(iopt1 > 2*num_indep) printf("wi failed, try help wi\n");
	else{
		if(debug)printf("%d arguments\n", iopt1);
		wiflag = 1;
		i = 0;
		while(i <= iopt1 && (i/2) < num_indep){
			xmin[i/2] = atof(cmdargs[i]);
			xmax[i/2] = atof(cmdargs[i + 1]);
		i +=2;
		}
	}
}

/* ve selects verbosity */
else if(strncmp(command,"ve",2) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && ( iopt1 ==2 || iopt1 == 1|| iopt1 == 0)){
		veflag =iopt1;
	}
	else printf("must be ve 0 1 or 2\n");
}

#ifdef OS2
/* fp selects every iteration plotting */
else if(strncmp(command,"fp",2) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && ( iopt1 ==2 || iopt1 == 1|| iopt1 == 0)){
		MYPLOT =iopt1;
		if(MYPLOT) fitcmd("doiky\n");
		if(!MYPLOT && pstream2 != NULL){
			fitcmd("quit\n");
			pclose(pstream2);
			pstream2 = NULL;
		}
	}
else printf("must be fp 0 1\n");
}
#endif

/* pa pauses for a number of seconds */
else if(strncmp(command,"pause",2) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && ( iopt1 > 0)){
		sleep(iopt1);
	}
	else if(sscanf(command,"%*s %d",&iopt1) > 0 && (iopt1 < 0)){
		printf("Hit Enter to continue\n");
		fflush(stdout);
		getc(stdin);
	}
	else help("pause",1);
}

/* li does a linear least squares fit */
else if(strncmp(command,"li",2) == 0){
	if(func == NULL){
		printf("You must select a function first\n");
		failed = 1;
	}
	else if(ndata < 1){
		printf("No data\n");
		failed = 1;
	}
	else if(linflag ==0){
		printf("Not a linear function\n");
		failed = 1;
	}
	else if(nonzero(data[order.sig], ndata)){
		if(debug)printf("in main, calling linear_fit()\n");
		failed = linear_fit(data,order,num_indep,ndata,itmax,a,ma,
			lista,mfit,covar,&chisq,func,filename,comment);
		if(debug)printf("in main, linear_fit() returned %d\n", failed);
	}
	else{
		failed = 1;
		printf("all sigmay's must be non-zero, check weighting and order\n");
	}
	if(failed == 0){
		printf("fit done\n");
		printf("%s\n",comment);
	}
	if(failed != 0)printf("fit failed\n");
}

/* de selects debugging option */
else if(strncmp(command,"de",2) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && ( iopt1 ==2 || iopt1 == 1|| iopt1 == 0)){
		debug =iopt1;
	}
	else printf("must be de 0 1 or 2\n");
}

/* noset*/
else if(strncmp(command,"noset",5) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && ( iopt1 == 1|| iopt1 == 0)){
		noset =iopt1;
	}
	else printf("must be noset 0 or 1\n");
}

/* echo selects debugging option */
else if(strncmp(command,"echo",2) == 0){
	if(sscanf(command,"%*s %d",&iopt1) > 0 && (iopt1 == 1|| iopt1 == 0)){
		echo =iopt1;
	}
	else printf("must be echo 0 or 1\n");
}

else if(strncmp(command,"run",2) == 0){
	if((i = parse(command, cmdargs)) > 0){
		strcpy(syscmd,"");
		for(j = 0; j < i; j++){
			strcat(syscmd, cmdargs[j]);
			strcat(syscmd," ");
			if(debug) printf("j: %d syscmd *%s* \n",j,syscmd);
		}
		if(debug) printf("Sending command %s to system\n",syscmd);
		system(syscmd);
	}
}

else if(strncmp(command,"!",1) == 0){
		if(debug) printf("Sending command %s to system\n",command + 1);
		system(command + 1);
}

/* gn send commands to gnuplot */
else if(strncmp(command,"gn",2) == 0){
	if((i = parse(command, cmdargs)) > 0){
		strcpy(gnubuf,"");
		for(j = 0; j < i; j++){
			strcat(gnubuf, cmdargs[j]);
			strcat(gnubuf," ");
			if(debug) printf("j: %d gnubuf *%s* \n",j,gnubuf);
		}
		sprintf(gnubuf,"%s\n",gnubuf);
		gnucmd( gnubuf);
	}
	else help("gn",1);
}

/* er estimates errors in the parameters */
else if(strncmp(command,"er",2) == 0){
	if(parse(command, cmdargs) > 0)
		est_errors(atof(cmdargs[0]),func, data, order, num_indep, a, ndata, ma, chisq);
		else help("er",1);
}

else if(strncmp(command,"quit",4)) printf("command %s not recognized\n",command);

fflush(stdout);
}
