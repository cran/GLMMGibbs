/* myrand.h */

/* Macros */

/* Choose one of the three gamma generators according to shape parameter */
#define gg(alpha) ((alpha)==1.0? unitexp():(alpha)<1.0? gs(alpha): gcf(alpha))

/* Random variates */
void seed(int s);
double u_random();
void asran_seed(int i1, int i2, int i3);
double asran(void);
int *rbern(int length, 
	   double *probabilities);

void snd_set(int dummy);
double snd(void);
void norm2(double *y1, double *y2);
double unitexp(void);
double chi2(int df);
double gs(double alpha);
double gcf(double alpha);
double gcra(double alpha);
double rbeta(double a,double b);
int sample_dd(int m, const double cdf[]);
double sample_cd(int ngrid, double low, double high,
		 const double grid[], double work[]);
double samp2_cd( int ngrid, const double x[], const double dens[],
		 const double d1[], double work[]);
double poissln(int ob, double ex, double priormean, double priorvar);
int r_poisson(double mu);
int r_binomial(int n, double p);
void r_multinom(int n, int k, double p[], int x[]);
int *rcat(int length,
          int categories,
	  double *probabilities,
	  int offset);

int rmvn(int n, double *mean, double *prec, int ix, double *x, double *w1,
         double *w2, int nn);




void generate_ctime_markov_random_chain(int number_of_states,
					double **transition_rates,
					int initial_state,
					int maximum_length,
					double time_period);






















