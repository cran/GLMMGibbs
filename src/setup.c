
#include "S.h"
#include "sparsev.h"
#include "sparsem.h"
#include "setup.h"
#include "misc.h"
#
GLM *g;


BLOCK *free_block(BLOCK *b);


void onemodel_model_init(int *nunits,double *response,
			 int *is_pw,double *pw,
			 int *is_offset,double *offset,     
			 int *family,double *sf)
     
{
  /* This routine is called by the R function "modeli.init" (see file 'glmm.R')
     and initiallises a model, which is pointed at by g, a pointer to an object
     of type GLM  */

  int i;
  

  
  /* We declare four vectors of pointers to functions representing
     the canonical link, its derivative, the variance functions and 
     the deviance contribution function (ie the function $d(y,mu)$
     such that the currrent deviance is $\sum_{i=1}^{n} d(y_{i},mu_{i})$.
     The first three take one | double| argument and the fourth
     takes two, and both return a | double| result.  The definitions
     PFD1 and PFD2 are given in `` c_glmmv2.h ''.  Notice that all
     this is set up so that, for exmple, | link_pointer[1](x)| will
     evaluate the link function at |x|, for the poisson family */

  PFD1 link_pointer[]                  = {logit_1,log_1,id_1};
  PFD1 diff_link_pointer[]             = {logit_2,log_2,id_2};
  PFD1 variance_function_pointer[]     = {bern_v,poisson_v,id_1};
  PFD2 deviance_contribution_pointer[] = {bern_d,poisson_d,id1};


  g=Calloc(1,GLM);

  g->lp = Calloc(*nunits,double);
  
  g->yvar = Calloc(*nunits,double);
  
  for(i=0;i< (int) *nunits;i++)
    {
      (g->yvar)[i]=response[i];
    }
  
  
  g->n = *nunits;
  
  
  
  /* If the prior weights are not  all 1, 
     copy them into |g->pw|.  Otherwise,
     point |g->pw| to a zero. */	    
  if((int) is_pw[0])  /* are there prior weights? */
    {
      g->prior_weight = (double *) 
	Calloc( (int) *nunits, double);
     
      for(i=0;i< (int) *nunits;i++) 
	{
	  (g->prior_weight)[i]=pw[i];
	}
    }
  else
    {
      g->prior_weight = (double *) 0;
    }
  
  
  
  /* If the offsets are not all 0, copy them into |g->offset|.  Otherwise,
     point |g->offset| to a zero. */	    
  if((int) is_offset[0])  /* Are there offsets? */ 
    {
      g->offset =  Calloc( (int) *nunits, double);
  
      
      for(i=0;i<*nunits;i++) 
	{
	  (g->offset)[i]=offset[i];
	}
    }
  else
    {
      g->offset = (double *) 0;
    }
  
  



     for(i=0;i< *nunits;i++) 
	{
	  (g->lp)[i]=g->offset ? offset[i] : 0.0;
	}


 g->inv_link =  link_pointer[*family];
  g->der_link =  diff_link_pointer[*family];
  g->var      =  variance_function_pointer[*family];
  g->dev      =  deviance_contribution_pointer[*family];
  g->scale    =  *sf; 

  g->rhs = (BLOCK *) 0;
  
    

}

void free_glm(int *dummy)
{
  BLOCK *b;
  
  b=g->rhs;
  
  while(b){
    b=free_block(b);
  }
  

  if(g->offset){
    Free(g->offset);
  }
  if(g->prior_weight){
    Free(g->prior_weight);
  }
  
  Free(g->yvar);
  
  Free(g->lp);

  Free(g);
}

BLOCK *free_block(BLOCK *b)
{
  BLOCK *res;
  
  res = b->next;
    if(b->s_matrix){
    
    sparsemat_freep(b->s_matrix);
  }
  if(b->design_matrix){
    
    sparsemat_freep(b->design_matrix);
  }

  if(b->effects){
    Free(b->effects);
  }

  Free(b);  

  return res;
  
}

onemodel_sample(int *program, int *program_length,  
            int *seed1, int *seed2, int *seed3,
		int *progress_info, int *thin, double *sampled_values,
		double *sampled_hyperparameter_values)
{

  sample(program, program_length, g,
	 seed1, seed2, seed3,
	 progress_info, thin, sampled_values, sampled_hyperparameter_values);

}



void r_to_block(int *n_units_p,
		int *NO_fdesign_matrix_columns_p,
                double *design_matrix_data,
                int *n_random_factors_p,
                int *n_levels_of_random_factors,
                int *levels_of_random_factors,
                int *random_describe,
                char **random_map_files,
                int *contrast_describe,
                char **fixed_column_names,
		char **random_factor_names,
                double *shape,
                double *scale,
		int *of_interest_p,
		int *only_ones)
{
  SPARSE_MAT *fixed_design_matrix;
  SPARSE_MAT *random_design_matrix;
  int n_units;
  int i,j;
  int n_random_factors;
  int number_of_fdesign_matrix_columns;
  SPARSE_MAT *precision_matrices,adj_matrix;
  int n_levels;
  int *product_random_factor;  
  int product_random_factor_nlevels;    
  SPARSE_MAT product_s_matrix;
  SPARSE_MAT contrast_matrix,product_contrast_matrix,p2;
  SPARSE_MAT smatrix;
  SPARSE_MAT *s_matrices;
  
  SPARSE_MAT *contrast_matrices;
  int **random_factor_levels;
  int *product_factor_nlevels; 
  SPARSE_MAT prior_estimates;
  int  degrees_of_freedom = 1;  
  int fixed_part,random_part;
  int of_interest;
  double *n_ones;
  
  SPARSE_MAT precision_matrix;
  BLOCK *this_block,*m;

  double zerop[] = {0.0};
  int oneip[] ={1};
  char *interceptstring[] = {"(Intercept)"};
  
  
  

  /* This routine takes the information about a block from R
     and puts the relevant information into the appropiate
     C Structure BLOCK  */




  

  

  

  int neffects;
  
  
  

  PFS  precision_funcs[]= { sparsemat_id_dumchar, 
                            sparsemat_diff1_Umatrix_dumchar, 
                            sparsemat_diff2_Umatrix_dumchar,
			    sparsemat_mapfile_Umatrix};
  
  PFSI contrast_funcs[] = { sparsemat_contr_id_fname,
 			    		    sparsemat_contr_sum_fname,
			    sparsemat_contr_treatment_fname,
			    sparsemat_contr_helmert_fname,
			    sparsemat_contr_diff2_fname };

  n_units                = (int) *n_units_p;
  
  n_random_factors       = (int) *n_random_factors_p;

  of_interest            = (int) *of_interest_p;
  n_units = g->n;

  /* Fistly, reserve space for the block and work out where it needs to
     go.  If we hhave not called "model block" yet, |g->rhs| will point
     to a |BLOCK * 0| object (becuase of beig set to that in |model_init|)
     and we point |glm->rhs|, the pointer to the first block of the model
     to |this_block|.  Otherwise we work through the linked list of blocks
     by moving from each block, |m|, to the one pointed to by |m->next|
     until we come to one for which |m->next = (BLOCK *) 0|. This is the 
     most recent block added, and we attach the new block to |m->next|  */


  this_block = Calloc (1,BLOCK);
  
  m = g->rhs;



  if (m)
    {
      i=2;
      while (m->next)
    	{
	  i++;
	  m  = m->next;
	}
      m->next = this_block;
    }
  else
    {
      i=1;
      g->rhs = this_block;
    }  

  /*Do the case for a single column of ones seperately */

  if (*only_ones){
    this_block->of_interest = 1;
    this_block->shape = 0.0;
    this_block->scale = 0.0;
    this_block->number_of_fixed_levels = 1;
    this_block->number_of_random_levels = 0;

    n_ones = Calloc(n_units,double);
    for(i=0;i<n_units;i++){
      n_ones[i]=1.0;
    }

    fixed_design_matrix = sparsemat_matrix_colnames(n_units,
						    1,
						    n_ones,
						    interceptstring);
    Free(n_ones);
   

    this_block->design_matrix = Calloc(1,SPARSE_MAT);
    
    sparsemat_copy(fixed_design_matrix, this_block->design_matrix);
    sparsemat_freep(fixed_design_matrix);

    this_block->fixed_column_names =interceptstring;
    this_block->random_factor_names = (char **) 0;
    this_block->hyperparameters = (double **) 0;
    this_block->number_of_effects = 1;
    this_block->effects = Calloc(1,double);
    this_block->s_matrix = (SPARSE_MAT *) 0;
    this_block->number_of_hyperparameters=0;
    
    
  }
  else{
    

  this_block->of_interest = of_interest;
  
  this_block->shape       = *shape;
  this_block->scale       = *scale;


  fixed_part = ( *NO_fdesign_matrix_columns_p > 0  );
  random_part = (  *n_random_factors_p > 0  );

  if (fixed_part){
    fixed_design_matrix = sparsemat_matrix_colnames(n_units,
						    *NO_fdesign_matrix_columns_p,
						    design_matrix_data,
						    fixed_column_names);
  }
  

  if(random_part){
    this_block->random_factor_names =
      (char **) calloc(1, (size_t) sizeof ( char *));

   this_block->random_factor_names[0] =
     (char *) calloc(strlen ( random_factor_names[0]) + 1,
		     (size_t) sizeof(char));
   
    

   strcpy(this_block->random_factor_names[0],
	  random_factor_names[0]);


       

    
    this_block->hyperparameters = 
      Calloc(1, double *);
    this_block->hyperparameters[0] =  
      Calloc(1, double );
    this_block->hyperparameters[0][0] = (double) 1.0;


   n_levels              =   (int) n_levels_of_random_factors[0];
	  precision_matrix             
	    = precision_funcs[(random_describe[0]-1)]
	    (n_levels,random_map_files[0]);

	  

    this_block->number_of_random_levels = n_levels_of_random_factors[0];
    contrast_matrix  
      = contrast_funcs[(contrast_describe[0]-1)](n_levels_of_random_factors[0],
						 1,
						 random_factor_names[0]); 

     

    random_design_matrix = sparsemat_sran_design(n_units,
						 levels_of_random_factors,
						 contrast_matrix);
	
  }
  





  if ((fixed_part) && (!(random_part))){
    this_block->design_matrix = Calloc(1,SPARSE_MAT);
    sparsemat_copy( fixed_design_matrix,    this_block->design_matrix);
    this_block->number_of_effects = *NO_fdesign_matrix_columns_p;
    this_block->number_of_hyperparameters = 0;
    
 
  }




  if ((!(fixed_part)) && (random_part)){
     this_block->design_matrix = Calloc(1,SPARSE_MAT);
    sparsemat_copy( random_design_matrix,    this_block->design_matrix);
    this_block->s_matrix = Calloc(1,SPARSE_MAT);

    
    spmatrix(&contrast_matrix,&precision_matrix,this_block->s_matrix);
  

  
    this_block->number_of_effects = contrast_matrix.number_of_columns;
    this_block->number_of_hyperparameters = 1;
    
    
      


 
  }
  neffects = this_block->design_matrix->number_of_columns;
  
    this_block->effects = Calloc(neffects,double );

  if(fixed_part){
    sparsemat_freep(fixed_design_matrix);
  }

  if(random_part){
    sparsemat_freep(random_design_matrix);
    sparsemat_free(contrast_matrix);
    
  }





    

  
  }
  


   
  
}


void model_show(int *dummy)
{
  int j;
  BLOCK *b;

  printf("Model information\n");
  printf("=================\n");
  printf("\n\n");
  
  printf("Number of units = %d", g->n);
  printf("\n");

  printf("response values: \n");
  printf("\n");

  
  for(j=0;j<(g->n);j++)
    {
      printf("[ %d ] %lf\n",j,(g->yvar)[j]);
    }
  printf("\n");
  
  if(g->prior_weight){
    printf("Prior Weights:\n");
    printf("\n");
    for(j=0;j<(g->n);j++)
      {
	printf("[ %d ] %lf\n",j, (g->prior_weight)[j]);
      }
  }
  else{
    printf("No prior weights\n");
  }
  if(g->offset){
    printf("Offsets:\n");
    printf("\n");
    for(j=0;j<(g->n);j++)
      {
	printf("[ %d ] %lf\n",j,(g->offset)[j]);
      }
  }
  else{
    printf("No offsets\n");
  }
  
   b= g->rhs;
   if(b){
     printf ("Design Matrix\n");
     printf ("-------------\n");
     printf ("\n");
     printf ("\n");
     
   }
   
   while(b){
     printf ("Model Block\n");
     printf ("-----------\n");
     printf ("\n");


       
      if(b->of_interest){
	printf("of interest\n");
      }
      else{
	printf("nuisance\n\n");
      }
      

      printf("shape = %lf\n",b->shape);
      printf("scale = %lf\n",b->scale);
      
      printf("Design Matrix: \n");
      
      sparsemat_pr(*(b->design_matrix));
      if(b->number_of_random_levels){
	
	
	printf("\n S Matrix: \n");
	
	sparsemat_pr(*(b->s_matrix));
      }
       
       
     
   b=b->next;
     
   }


     

}
  




void nparameters(int *number_of_parameters,
		 int *number_of_hyperparameters)
{
  BLOCK *b;
  *number_of_parameters      = 0;
  *number_of_hyperparameters = 0;
  
  b=g->rhs;
  
  while (b){
    if(b->of_interest){
      
   *number_of_parameters      += b->number_of_effects;

   

    }
   *number_of_hyperparameters  += b->number_of_hyperparameters;    
  b=b->next;
  
  }
  


  
  
}


void parameternames(char **parameter_names,
		    char **hyperparameter_names)
{
  BLOCK *b;
  int i,j,ii=0,jj=0;
  int n;
  char *x;
  
   

 

  
  b=g->rhs;

  while (b){

    if(b->of_interest){

    for(i=0;i<b->number_of_effects;i++){
      
      
      x=b->design_matrix->column_names[i];
       parameter_names[ii]=x; 
      
      ii++;
     
      
      
    }
    }


 

    if(b->number_of_random_levels){
    
      hyperparameter_names[0]= b->random_factor_names[0]  ;
      
      jj++;
          
      
    }

    
    b=b->next;
  }
   

	  
}





void bugs_output(int *number_of_parameters,
		 int *number_of_iterations,
		 char **parameter_names,
		 double *parameter_values,
		 char **ind_name,
                 char **out_name)
{
  int i,j;
  FILE *find,*fout;

  find = fopen(ind_name[0],"w");
  fout = fopen(out_name[0],"w");
  
  
  
  for (i=0;i<*number_of_parameters;i++){
    fprintf(find,"%-20.20s %-8d %-8d \n",parameter_names[i], ((i * *number_of_iterations) +1),
                ((i+1) * *number_of_iterations )                       );
    
  }
  for (i=0;i<*number_of_parameters;i++){
    for (j=0;j<*number_of_iterations;j++){
      fprintf(fout,"%10.10d %lf \n",(j+1),parameter_values[j + *number_of_iterations * (i+1) ]);
      
    }
  }
  
  fclose(find);
  fclose(fout);
  
    
  
}


