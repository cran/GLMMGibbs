
typedef double (*PFD1)(double);  /* Pointer to a function returning double */
typedef double (*PFD2)(double, double); /* Ditto ... with 2 arguments */


typedef struct 
{
  double **all_values;
  double **effect_values;
  double **hyperparameter_values;
  char   **all_names;
  char   **effect_names;
  char   **hyperparameter_names;
  int    number_of_effects;
  int    number_of_hyperparameters;
  int    number_of_iterations;
  int    program_length;
  int    *program;
} 
OUTPUT;


typedef SPARSE_MAT (*PFS)(int, char *); 
typedef SPARSE_MAT (*PFSI)(int, int, char *);

typedef  struct block
{
  int number_of_fixed_levels;
  int number_of_random_levels;
  SPARSE_MAT *design_matrix;
  double *effects;
  SPARSE_MAT *s_matrix;
  int of_interest;
  double **hyperparameters;
  SPARSE_MAT prior_estimates;
  int degrees_of_freedom;
  struct block *next;
  char **fixed_column_names;
  char **random_factor_names;
  double shape;
  double scale;
  int number_of_effects;
  int number_of_hyperparameters;
  
}
BLOCK;


typedef struct
{
  int number_of_parameters;
  
  char **names;
  double **results;
}
STATISTICS ;



typedef struct {
  int n;                        /* Number of units */
  double scale;                 /* Scale factor */
  
  /* Pointers to arrays */
  
  double *yvar;                 /* Y-variate vector */
  double *lp;                   /* Current linear predictor (local) */
  double *prior_weight;         /* Prior weights */
  double *offset;               /* Offsets */
  
  /* Pointers to functions */
  
  PFD1 inv_link;                /* Inverse link function */
  PFD1 der_link;                /* First derivative of link function */
  PFD1 var;                     /* Variance function */
  PFD2 dev;                     /* Contribution to deviance */
  
  /* Pointers defining right-hand side of model etc. */
  
  BLOCK *rhs;                   /* List of terms in model formula */
  
  /* Pointer to result buffer */
  
  
} GLM;


double log_1(double);
double log_2(double);

double logit_1(double);
double logit_2(double);

double id_1(double);
double id_2(double);

double poisson_v(double);
double poisson_d(double, double);

double bern_v(double);
double bern_d(double, double);

double norm_v(double);
double norm_d(double, double);


void step_b(GLM *glm, SPARSE_VEC *x, double pr_mean, double pr_prec,
	    double b_old, double b_new, double *pd, double *g, double *h);


int glm_fit_v(GLM *glm, SPARSE_VEC *x, double *b, double pr_mean, 
	      double pr_prec);
     



















