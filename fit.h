#define NUM_ARGS 20
#define COMMAND_SIZE 80

/* This structure tells us how the data is represented internally in
     the data matrix.  If order.x = 0, then the 0th column is the x's.
/* 	int *x;        independent variables */
 /*	int y;         dependent variable */
 /*	int yfit;      value of function with current parameters */
 /*	int sig;       sigma: an error estimate for y, used in fitting */
 /*	int nsig;     no sigma, column of data full of 1's */
 /*	int ssig;     sigma for statistical weighting */
 /*	int isig;      sigma for instrumental weighting */
 /*	int osig;     sigma for other weighting */

struct data_order{
	int *x;
	int y;
	int yfit;
	int sig;
	int nsig;
	int ssig;
	int isig;
	int osig;
	int *xsig;
};

/*** a few function declarations ***/

double **dmatrix(int,int);/* allocates a two dimensional array of doubles */
double *dvector(int);            /* allocates a one dimensional array of doubles */
int *ivector(int);                   /* allocates a one dimensional array of ints */
void free_dmatrix(double **, int,int);     /* frees a two dimensional array of doubles */
int listfcns(void);                /* lists the functions available for fitting */

/* pointer to the fitting function */
extern int (*func)(double *x, double *a, double *y, double *dyda,
				int na, int dyda_flag, int *fita, int dydx_flag, 
				double *dydx, double ydat);

/* returns pointer to the fitting function */
int (*getfuncptr())(char *function_name, int *num_indep,
			int linflag, int *num_parameters, char *comment);

/* calculates f(x) with current parameter values */
int calc_yfit(int (*func)(), double **data, struct data_order order, 
				int num_indep, double *a, int ndata, int ma, double *chisqr);

int help(char *topic, int maxlines);

/* reads the data from a file */
int get_data(double **data, char *command, int num_indep, char *inbuf,
				struct data_order *order, char *filename,
				int maxrows, int maxcols);

/* does the nonlinear fit */
int mrqfit(double **data,struct data_order order,
				int num_indep, int ndata, int itmax, double *a,
				int ma, int *lista, int mfit, double **covar,
				double *chisq, int (*func)(), 
				char *filename, char *comment);

void solve_for_da();
int alpha_beta_chisq(double **data, struct data_order order,
				int num_indep, int ndata, double *a, int ma, 
				int *lista, int mfit, double **alpha, double *beta,
				double *chisq, int (*funcs)());

int plot( int (*func)(), double **data, struct data_order order,
				int num_indep, int ndata, char *filename, char *comment, 
				double *a, int ma);

#ifdef OS2
int fitcmd(char *command);

int myplot( int (*func)(), double **data, struct data_order order,
				int num_indep, int ndata, char *filename, char *comment,
				double *a, int ma);
#endif

int parse(char *command, char cmdarg[30][30]);

int make_data(int (*func)(), int num_indep, double *a, int ma, 
				char *command, char *inbuf);

/* nonzero returns 1 if all elements of an array are non-zero */
int nonzero(double *, int ndata);

int gnucmd(char *command);

int linear_fit(double **data,struct data_order order,
				int num_indep, int ndata, int itmax, double *a,
				int ma, int *lista, int mfit, double **covar,
				double *chisq, int (*func)(),
				char *filename, char *comment);

int window(double *x, int num_indep);

void est_errors(double factor,int (*func)(), double **data,
				struct data_order order, int num_indep,
				double a[],int ndata, int ma, double chisqr);

