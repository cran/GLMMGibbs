
typedef struct {
  SPARSE_VEC **columns;
  int number_of_rows;
  int number_of_columns;
  double offset;
  char **column_names;
  char type;
  
  }
SPARSE_MAT;

SPARSE_MAT *sparsemat_sran_design(int, int *, SPARSE_MAT contrast);


SPARSE_MAT sparsemat_id(int);
SPARSE_MAT sparsemat_id_dumchar(int, char *dumchar);
SPARSE_MAT sparsemat_id_fname(int, char *factor_name);


SPARSE_MAT sparsemat_contr_id_fname(int n, int dummy, char *factor_name);

SPARSE_MAT sparsemat_contr_sum(int, int);
SPARSE_MAT sparsemat_contr_sum_fname(int, int, char *factor_name);


SPARSE_MAT sparsemat_contr_treatment(int n,
				     int zlevel);

SPARSE_MAT sparsemat_contr_treatment_fname(int, int, char *factor_name);

SPARSE_MAT sparsemat_contr_helmert(int, int);
SPARSE_MAT sparsemat_contr_helmert_fname(int, int, char *factor_name);


SPARSE_MAT sparsemat_contr_diff2(int,int);
SPARSE_MAT sparsemat_contr_diff2_fname(int, int, char *factor_name);


SPARSE_VEC *sparsevec_elementwise_sum(SPARSE_VEC *,
                                    SPARSE_VEC *,
                                    int);

SPARSE_VEC *sparsevec_elementwise_product(SPARSE_VEC *,
                                    SPARSE_VEC *,
                                    int);
SPARSE_VEC *sran_design(int g,
                        int number_of_contrasts,
                        int number_of_units,
                        int *group,
                        SPARSE_VEC *contrast);

SPARSE_MAT sparsemat_elementwise_sum(SPARSE_MAT, SPARSE_MAT);
SPARSE_MAT sparsemat_elementwise_product(SPARSE_MAT, SPARSE_MAT);
SPARSE_MAT sparsemat_column_by_column_product(SPARSE_MAT mat1,
					      SPARSE_MAT mat2);
void sparsemat_pr(SPARSE_MAT);
void sparsemat_pr2(SPARSE_MAT);
SPARSE_VEC *sparsevec_scalar_multiply(double scalar, SPARSE_VEC *mat);
SPARSE_VEC *sparsevec_scalar_add(SPARSE_VEC *, double);
SPARSE_MAT sparsemat_scalar_add(SPARSE_MAT, double);
SPARSE_VEC *sparsevec_ncopies(SPARSE_VEC *, int);
SPARSE_MAT sparsemat_diff1_Umatrix( int);
SPARSE_MAT sparsemat_diff1_Umatrix_dumchar( int, char *dumchar);
SPARSE_MAT sparsemat_diff2_Umatrix( int);
SPARSE_MAT sparsemat_diff2_Umatrix_dumchar( int, char *dumchar);
int iszeroa(double);

SPARSE_MAT sparsemat_rd(FILE *, int, int);
SPARSE_MAT *sparsemat_matrix(int, int, double *);

SPARSE_MAT *sparsemat_matrix_colnames(int, int, double *,char **);
SPARSE_MAT map_U_matrix(SPARSE_MAT);
SPARSE_MAT sparsemat_mapfile_Umatrix(int nregions,
                                     char *filename);

void sparsemat_dump(SPARSE_MAT);


SPARSE_MAT sparsemat_zoff_direct_product(SPARSE_MAT, SPARSE_MAT);
SPARSE_MAT sparsemat_direct_product(SPARSE_MAT, SPARSE_MAT);
SPARSE_MAT sparsemat_ones(int nrows, int ncols);

int sparsemat_is_identity(SPARSE_MAT);

double sparsemat_qform(SPARSE_VEC *vec1,
		       SPARSE_VEC *vec2,
                       SPARSE_MAT mat);

double sparsemat_dqform(double *vec1,
		       double *vec2,
			SPARSE_MAT mat);

SPARSE_MAT sparsemat_btrans_a_b(SPARSE_MAT matb,
				SPARSE_MAT mata);

SPARSE_MAT s_matrix(SPARSE_MAT  contrast_matrix,
		    SPARSE_MAT u_matrix);


void sparsemat_free(SPARSE_MAT);
void sparsemat_freep(SPARSE_MAT  *mat);
void sparsemat_copy(SPARSE_MAT *from, SPARSE_MAT *to);


void sparsemat_n_product(SPARSE_MAT *terms,
			 int n_terms,
			 SPARSE_MAT *product);

SPARSE_MAT sparsemat_scale(SPARSE_MAT mat,
			   double maximum,
			   double minimum,
			   double tolerance);

      
void matrix_elements(SPARSE_MAT mat, double *res);
SPARSE_MAT smatrix_uidentity(SPARSE_MAT contrast_matrix);
SPARSE_MAT sparsemat_diag(double *data,int n);


SPARSE_MAT *submatrix(SPARSE_MAT *matrix,
		      int first_row,
		      int last_row,
		      int first_column,
		      int last_column);





void  spmatrix(SPARSE_MAT *contrast,
	      SPARSE_MAT *u_matrix,
	      SPARSE_MAT *res);

