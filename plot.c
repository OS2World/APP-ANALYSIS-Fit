#include <stdio.h>
#include <stdlib.h>
#include "fit.h"

extern FILE *pstream;
extern double *xmin, *xmax;
extern int wiflag;
extern int debug;
extern int noset;

int plot(func,data, order, num_indep, ndata, filename, comment, a, ma)
int (*func)();
double **data;
struct data_order order;
int num_indep;
int ndata;
char *filename;
char *comment;
double *a;
int ma;
{
int i, failed;
FILE *stream;
char gnubuf[100];
double umin = 1e200, umax= -1e200, vmin = 1e200, vmax=-1e200;
double chisqr;
int num_samples;
char fitfile[40];

#ifdef DOS
gnucmd(gnubuf);
return 1;
#endif
#ifdef UNIX
strcpy(fitfile,"/tmp/fit.tmp");
#endif
#ifndef UNIX
strcpy(fitfile,"fit.tmp");
#endif
 
if(ndata < 1 ){
	printf("pl failed, no data\n");
	return 0;
}
if( num_indep > 2){
	printf("pl failed, can only plot 1 or 2 independent variables\n");
	return 0;
}

sprintf(gnubuf,"clear\n");
gnucmd( gnubuf);

if(num_indep == 1 && !noset){
	sprintf(gnubuf,"set aut\n");
	gnucmd( gnubuf);
	sprintf(gnubuf,"set xlabel 'x'\n");
	gnucmd( gnubuf);
	if(wiflag){
		sprintf(gnubuf,"set xrange [%g:%g]\n", xmin[0], xmax[0]);
		gnucmd( gnubuf);
	}
	sprintf(gnubuf,"set ylabel 'f(x)'\n");
	gnucmd( gnubuf);
	sprintf(gnubuf,"set title 'data and fit'\n");
	gnucmd( gnubuf);
	}

	if(num_indep == 2 && !noset){
	sprintf(gnubuf,"set aut\n");
	gnucmd( gnubuf);
	sprintf(gnubuf,"set xlabel 'x'\n");
	gnucmd( gnubuf);
	if(wiflag){
		sprintf(gnubuf,"set xrange [%g:%g]\n", xmin[0], xmax[0]);
		gnucmd( gnubuf);
		sprintf(gnubuf,"set yrange [%g:%g]\n", xmin[1], xmax[1]);
		gnucmd( gnubuf);
	}
	sprintf(gnubuf,"set ylabel 'y'\n");
	gnucmd( gnubuf);
	sprintf(gnubuf,"set zlabel 'f(x,y)'\n");
	gnucmd( gnubuf);
	sprintf(gnubuf,"set title 'data and fit'\n");
	gnucmd( gnubuf);
}


/* function is plotted through a defined gnuplot function */
if(strncmp(comment,"FILE",4) != 0){

        /* define the parameters in gnuplot */
	for(i = 0; i < ma; i++){
		sprintf(gnubuf,"a%d = %g\n",i,a[i]);
		gnucmd(gnubuf);
	}

        /* define the function in gnuplot */
	sprintf(gnubuf,"%s\n",comment);
	gnucmd( gnubuf);

        /* tell gnuplot to do the plot */
	if(num_indep == 1){
		sprintf(gnubuf,"plot '%s' using %d:%d, f(x) \n",
			filename,order.x[0] + 1, order.y + 1);
		gnucmd( gnubuf);
		if(!noset){
			if(ndata > 40) num_samples = ndata/2;
			else if(ndata > 20) num_samples = ndata;
			else num_samples = 20;
			if(num_samples > 500) num_samples = 500;
			sprintf(gnubuf,"set samples %d\n",num_samples);
			gnucmd(gnubuf);
		}
	}

/* for 2 independent variables, we need to use the gnuplot */
/* command splot and plot parametrically.  See gnuplot documentation */
/* for the reasons why.  We also need to know tell gnuplot */
/* the ranges of the independent variables, so we calculate these also */

	if(num_indep == 2){
		if(wiflag == 0)
			for(i = 0; i < ndata; i++){
				if(data[order.x[0]][i] < umin) umin = data[order.x[0]][i];
				if(data[order.x[0]][i] > umax) umax = data[order.x[0]][i];
				if(data[order.x[1]][i] < vmin) vmin = data[order.x[1]][i];
				if(data[order.x[1]][i] > vmax) vmax = data[order.x[1]][i];
			}
		else{
			umin = xmin[0];
			umax = xmax[0];
			vmin = xmin[1];
			vmax = xmax[1];
		}
		sprintf(gnubuf,"set parametric\n");
		gnucmd( gnubuf);
		if(1){
			sprintf(gnubuf,"set urange [%g:%g]\n", umin, umax);
			gnucmd( gnubuf);
			sprintf(gnubuf,"set vrange [%g:%g]\n", vmin, vmax);
			gnucmd(gnubuf);
			sprintf(gnubuf,"set samples 40\n");
    		gnucmd(gnubuf);
		}
		sprintf(gnubuf,"splot '%s' using %d:%d:%d, u,v,f(u,v) \n",
			filename,order.x[0] + 1, order.x[1] + 1, order.y + 1);
		gnucmd(gnubuf);
		sprintf(gnubuf,"set noparametric\n");
		gnucmd(gnubuf);
	}

}

/* function is plotted through a temporary file */
else if(strncmp(comment,"FILE",4) == 0){
	if(debug) printf("in plot(), about to call calc_yfit()\n");

/* calculate the values of the fitting function for current parameters */
	failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisqr);
	if(failed) return 1;
	if(debug) printf("in plot(), returned from calc_yfit()\n");
	stream = fopen(fitfile,"w");
	if(num_indep == 1)
		for(i = 0; i < ndata; i++) 
			if(wiflag == 0 || (data[order.x[0]][i] >= xmin[0] && data[order.x[0]][i] <= xmax[0]))
				fprintf(stream,"%g %g\n", data[order.x[0]][i], data[order.yfit][i]);
	if(num_indep == 2)
		for(i = 0; i < ndata; i++)
			if(wiflag == 0 || (data[order.x[0]][i] >= xmin[0] && data[order.x[0]][i] <= xmax[0]) && (data[order.x[1]][i] >= xmin[1] && data[order.x[1]][i] <= xmax[1]))
				fprintf(stream,"%g %g %g\n", data[order.x[0]][i], data[order.x[1]][i], data[order.yfit][i]);
	fclose(stream);
	if(num_indep == 1){
		sprintf(gnubuf,"plot '%s' using %d:%d,'%s'\n",
                                filename,order.x[0] + 1, order.y + 1, fitfile);
		gnucmd(gnubuf);
	}
	if(num_indep == 2){
		sprintf(gnubuf,"set parametric\n");
		gnucmd(gnubuf);
		sprintf(gnubuf,"splot '%s' using %d:%d:%d,'%s' \n",
                                filename,order.x[0] + 1, order.x[1] + 1,order.y + 1, fitfile);
		gnucmd(gnubuf);
		sprintf(gnubuf,"set noparametric\n");
		gnucmd(gnubuf);
	}

}
return 0;
}

/* the pr function plot residual errors */
int pr(command,func,data, order, num_indep, ndata, filename, comment, a, ma)
char *command;
int (*func)();
double **data;
struct data_order order;
int num_indep;
int ndata;
char *filename;
char *comment;
double *a;
int ma;
{
int i, flag, failed;
FILE *stream;
char gnubuf[100];
double chisqr;
char fitfile[40];

#ifdef DOS
gnucmd(gnubuf);
return 1;
#endif
#ifdef UNIX
strcpy(fitfile,"/tmp/fit.tmp");
#endif
#ifndef UNIX
strcpy(fitfile,"fit.tmp");
#endif

if(debug)printf("in pr() command: %s\n", command);
if(ndata < 1 ){
	printf("pr failed, no data\n");
	return 0;
}
if( num_indep > 2){
	printf("pr failed, can only plot 1 or 2 independent variables\n");
	return 0;
}

sprintf(gnubuf,"clear\n");
if( sscanf(command,"%*s %d", &flag ) == 0) flag = 1;
if(debug)printf("flag: %d\n", flag);

failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisqr);
if(failed) return 1;

/* option 1, plot yfit vs. y */
if(flag ==1){
	stream = fopen(fitfile,"w");
	for(i = 0; i < ndata; i++) 
		if(!wiflag || (data[order.x[0]][i] >= xmin[0] && data[order.x[0]][i] <= xmax[0])) fprintf(stream,"%g %g\n", data[order.yfit][i], data[order.y][i]);
	fclose(stream);
	if(!noset){
		sprintf(gnubuf,"set aut\n");
		gnucmd( gnubuf);
		sprintf(gnubuf,"set title 'data vs. fit'\n");
		gnucmd( gnubuf);
		sprintf(gnubuf,"set xlabel 'f(x)'\n");
		gnucmd( gnubuf);
		if(wiflag){
			sprintf(gnubuf,"set xrange [%g:%g]\n", xmin[0], xmax[0]);
			gnucmd( gnubuf);
		}
		sprintf(gnubuf,"set ylabel 'y'\n");
		gnucmd( gnubuf);
	}
	sprintf(gnubuf,"plot '%s' with points \n", fitfile);
	gnucmd(gnubuf);
}

/* option 2, plot y-yfit as a function of independent variables */
else{
	if(debug)printf("stream = fopen();\n");
	stream = fopen(fitfile,"w");
	if(debug)printf("if(num_indep ==1);\n");
	if(!noset){
		sprintf(gnubuf,"set aut\n");
		gnucmd( gnubuf);
		sprintf(gnubuf,"set title 'data-fit'\n");
		gnucmd( gnubuf);
	}
	if(num_indep ==1){
		for(i = 0; i < ndata; i++)
			if(!wiflag || (data[order.x[0]][i] >= xmin[0] && data[order.x[0]][i] <= xmax[0]))
				fprintf(stream,"%g %g\n", data[order.x[0]][i],
						data[order.y][i]-data[order.yfit][i]);
	}
	else if(num_indep ==2){
		for(i = 0; i < ndata; i++)
			if(!wiflag || (data[order.x[0]][i] >= xmin[0] && data[order.x[0]][i] <= xmax[0]))
				fprintf(stream,"%g %g %g\n", data[order.x[0]][i],data[order.x[1]][i],
						data[order.y][i]-data[order.yfit][i]);
	}
	fclose(stream);
	if(num_indep ==1){
		if(!noset){
			sprintf(gnubuf,"set aut\n");
			gnucmd( gnubuf);
			if(wiflag){
				sprintf(gnubuf,"set xrange [%g:%g]\n", xmin[0], xmax[0]);
				gnucmd( gnubuf);
			}
			sprintf(gnubuf,"set title 'y - f(x)'\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set xlabel 'x'\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set ylabel 'y - f(x)'\n");
			gnucmd( gnubuf);
		}
		sprintf(gnubuf,"plot '%s'\n", fitfile);
		gnucmd(gnubuf);
	}	
	if(num_indep ==2){
		if(!noset){
			sprintf(gnubuf,"set aut\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set title 'y - f(x)'\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set xlabel 'x'\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set zlabel 'y - f(x)'\n");
			gnucmd( gnubuf);
			sprintf(gnubuf,"set ylabel 'y'\n");
			gnucmd( gnubuf);
		}
		sprintf(gnubuf,"set parametric\n");
		gnucmd(gnubuf);
		sprintf(gnubuf,"splot '%s'\n", fitfile);
		gnucmd(gnubuf);
		sprintf(gnubuf,"set noparametric\n");
		gnucmd(gnubuf);
	}
}
return 0;
}

#ifdef OS2

int myplot(func,data, order, num_indep, ndata, filename, comment, a, ma)
int (*func)();
double **data;
struct data_order order;
int num_indep;
int ndata;
char *filename;
char *comment;
double *a;
int ma;
{
int i, failed, j;
FILE *stream;
char buf[100];
double chisqr;
char fitfile[40];
double x;

strcpy(fitfile,"fit.tmp");

if(ndata < 1 ){
	printf("pl failed, no data\n");
	return 0;
}
if( num_indep > 1){
	printf("pl failed, can only plot 1 indep variable\n");
	return 0;
}

/* function is plotted through a temporary file */

/* calculate the values of the fitting function for current parameters */
/*
if(debug) printf("in plot(), about to call calc_yfit()\n");
failed = calc_yfit(func, data, order, num_indep, a, ndata, ma, &chisqr);
if(failed) return 1;
if(debug) printf("in myplot(), returned from calc_yfit()\n");
*/

stream = fopen(fitfile,"w");
for(i = 0; i < ndata; i++){
		x = data[order.x[0]][i];
		if( !wiflag || (x > xmin[0])  && (x < xmax[0]))
			fprintf(stream,"%g %g %g\n", 
					x, data[order.y][i], data[order.yfit][i]);
}
fclose(stream);

sprintf(buf,"gd fit.tmp 2\n");
fitcmd(buf);
sprintf(buf,"sy 0 L\n");
fitcmd(buf);
sprintf(buf,"pl\n");
fitcmd(buf);
return 0;
}

#endif

