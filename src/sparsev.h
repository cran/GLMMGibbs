/* sparsev.h */

/*
   Sparse vector representation.
   Useful for autoregressive schemes, factor represention etc.

   Stores list of non-zero rows and (when they are not equal to 1) the
   corresponding values.

   Can also be used for conventional vector - list of non-zero rows no
   longer needed.
*/

typedef struct {
   int n;                              /* Effective length */
   int non_zero;                       /* Number of non-zero rows */
   int *rows;                          /* Array of row numbers, 0 ... n-1 */
   double *vals;                       /* Array of values */
} SPARSE_VEC;

SPARSE_VEC *sparse_id(int);            /* Create identity matrix */
SPARSE_VEC *sparse_one(int);           /* Create unit vector */
SPARSE_VEC *sparse_ref(int, double*);  /* SPARSE_VEC reference to vector */
SPARSE_VEC *sparse_vec(int, double*);  /* Standard vector -> SPARSE_VEC form */
SPARSE_VEC *sparse_mat(int, int, double*); /* Standard matrix -> SPARSE_VEC */
SPARSE_VEC *sparse_gl(int, int, int);  /* GLIM gl() function */
SPARSE_VEC *sparse_fl(int, int, int*); /* Expand vector of factor levels */
SPARSE_VEC *sparse_rd(FILE*, int, int);/* Read sparse matrix */
void sparse_wr(FILE*, SPARSE_VEC*, int);/* Write a sparse matrix */
SPARSE_VEC *diff1_smooth(int);         /* First difference smoother */
SPARSE_VEC *diff2_smooth(int);         /* Second difference smoother */
double sparse_el(SPARSE_VEC*, int);    /* Return element of sparse vector */
int is_symmetric(SPARSE_VEC*, int);    /* Check for matrix symmetry */
int islands(SPARSE_VEC *m, int size);  /* # islands for adjacency matrix */
SPARSE_VEC *sparse_tr(SPARSE_VEC*, int);/* Transpose a sparse matrix */
SPARSE_VEC *sparse_cr(int, int, ...);  /* Create sparse vector */
SPARSE_VEC *sparse_c1(int, int, ...);  /* Create sparse vector of 1.0 values*/
void sparse_fr(SPARSE_VEC *, int);     /* Free dynamic storage */
SPARSE_VEC *sparse_kpv(SPARSE_VEC*, SPARSE_VEC*);
                                       /* Direct product of vectors */
SPARSE_VEC *sparse_kp(int, SPARSE_VEC*, int, SPARSE_VEC*);
                                       /* Direct product of matrices */
SPARSE_VEC *sparse_0(int, int);        /* Create zero matrix */
SPARSE_VEC *sparse_ccv(int, SPARSE_VEC*, SPARSE_VEC*);
                                       /* Concatenate matrices vertically */
SPARSE_VEC *sparse_cch(int, SPARSE_VEC*, int, SPARSE_VEC*);
                                       /* Concatenate matrices horizontally */
void sparse_cp(SPARSE_VEC*, SPARSE_VEC*); /* Copy vector */
double sparse_ip(SPARSE_VEC*, SPARSE_VEC*, double*); /* Inner product */
double sparse_ss(SPARSE_VEC*, double*); /* Sum of squares */
int sparse_so(SPARSE_VEC*, SPARSE_VEC*); /* Strictly orthogonal vectors?*/
SPARSE_VEC *sparse_gm_cr(int, SPARSE_VEC*); /* Create Grammian template */
void sparse_gm_cal(int, SPARSE_VEC*, double*, SPARSE_VEC*);
                                       /*Calculate Grammian*/
void elements(SPARSE_VEC*, double *where);






SPARSE_VEC *select_elements(SPARSE_VEC *vector,
		   int first,
		   int last);

double last_value(SPARSE_VEC *vector);
