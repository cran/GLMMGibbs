
 double dev_tol = 0.001; /* Convergence test (deviance) */
static double dev_too_far = 5.0; /* Sampling starting point tolerance */
static int arse_max = 20;      /* Maximum sampling envelope size */
static int mode_op = 0;        /* Mode of operation:
                                  0: Gauss-Siedel (Besag's ICM) 
				  >0: Gibbs sampling
                                  1: Adaptive rejection 
                                  2: Gaussian approx to conditionals
                                  3: Gaussian approx + Metropolis rejection */
static int ncycles = 5;        /* Maximum iterations for conditional mode */
static int nmet = 0;           /* Metropolis steps per sample */
static int nsamp = 0;          /* Count of sampled random variables */
static int nfail = 0;          /* Count of failures */
static int silent = 1;          /* 1=working silently, 0=reporting progress */ 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "sparsev.h"
#include "sparsem.h"
#include "setup.h"
#include "myrand.h"
#include "ars.h"
#include "misc.h"
#define SMALL_REAL 0.00001 
#include "S.h"


void Rprintf(char*, ...);


  
void print_current_parameter_values(GLM *glm);

void prior_predict(int n_levels,
                   double *current_parameter_values,
                   double current_hyperparameter_value,
                   int which,
                   SPARSE_MAT precision,
                   double *mean,
                   double *prec);




void sample(int *program, int *program_length,  GLM *g,
            int *seed1, int *seed2, int *seed3,
            int *progress_info, int *thin,double *sampled_values, 
	    double *sampled_hyperparameter_values)
{
  

int i,j,k,l,ii;
BLOCK *b;
FILE *tmp,*bugind;
double *mean,*prec;
int gibbs_steps;
double shape,scale;
SPARSE_VEC *effects_vector;
double ssq;
int number_of_parameters;
int number_of_hyperparameters;
int nrows,ncols;
int parameter_count;
int hyperparameter_count;
 int saved_iteration=0;

 int wheretoput;
asran_seed(*seed1,*seed2,*seed3); 

 


/* We inittialize pointers to doube meean and prec.
   While the levels of a random factor are being 
   generated mean[0] and prec[0] are the mean and
   precision of the current level*/



fflush(stdout);

mean = Calloc (1,double);
prec = Calloc (1,double);
mean[0] = 0.0;
prec[0]  = 0.0;


/* We sample from the model.*/

b=g->rhs;

l=0;
for(k=0;k<(*program_length/2);k++)
  {
    set_mode(program[k*2]);




    if(*progress_info){
      
    Rprintf("\n");
    Rprintf("Mode of operation and operating parameters\n");
    Rprintf("==========================================\n");
    Rprintf("Estimation mode %d : ", mode_op);
    if (mode_op == 0) {
      Rprintf("Iterative Conditional Mode algorithm\n");
    } else {
      Rprintf("Gibbs sampling ");
      switch (mode_op) {
      case 1: Rprintf(" (adaptive rejection)"); break;
      case 2: Rprintf(" (normal approximation)"); break;
      case 3: Rprintf(" (Metropolis algorithm), "); 
      }
      Rprintf("\n");
    }
    if(mode_op ==1)ncycles=0;
    if(mode_op ==0)ncycles=5;
    if (mode_op != 1) Rprintf("Maximum iterations to find mode = %d\n", ncycles);
    }
    
  





 
    for(i=0;i<program[(k*2)+1];i++)
      {
        parameter_count = 0;
	hyperparameter_count = 0;
      	for(b=g->rhs;b;b=b->next)

	  {
	    for(j=0;j<(b->design_matrix)->number_of_columns;j++)
	      {
		
		    if (b->number_of_random_levels) 
		      {
			prior_predict((b->design_matrix)->number_of_columns,
				      b->effects,
				      b->hyperparameters[0][0],
				      j,
				      *(b->s_matrix),
				      mean,
				      prec);
		    
		      }
                else{
		  mean[0]=0.0;
		  prec[0]=0.0;
		}

		    
		glm_fit_v(g,(b->design_matrix)->columns[j],
			  (b->effects)+j,*mean,*prec);

	 
	    
		 if(k==(*program_length/2 - 1))
		  {
		    if(b->of_interest){
		     
		      if (i % *thin == 0){
                        
wheretoput = parameter_count*((program[(k*2)+1])/(*thin))  + saved_iteration ;

			sampled_values[wheretoput]
			  =*( (b->effects)+j );

			
		      } 
		      
		      parameter_count++;		   	           
		    }
		
		    

	      }

              }

	      

            if((program[2*k]==1) && (b->number_of_random_levels)){
            ssq = sparsemat_dqform(  b->effects,
				     b->effects,
				     *(b->s_matrix));
	    scale = b->shape + ssq/2;
	    shape = b->scale + b->number_of_random_levels/2;
	  	      
            b->hyperparameters[0][0] = gg(shape)/scale;
	    }
	    if(k==(*program_length/2 - 1) && (b->number_of_random_levels)){
	      if (i % *thin == 0 ){



wheretoput = hyperparameter_count*((program[(k*2)+1])/(*thin))  + saved_iteration ;
			sampled_hyperparameter_values[wheretoput]
			  =  b->hyperparameters[0][0];


	
	    
	      }
	    hyperparameter_count++;	      

	    
      	    }
	    
	    
	  }

        if(*progress_info > 0){
	  
	if (i % (*progress_info) == 0)  Rprintf("%d ",i);
        if(i % (*progress_info *10) == 0) Rprintf("\n");
        fflush(stdout); 
	}
	

		 if(k==(*program_length/2 - 1))
		  {
		    if(i % *thin == 0){
		      saved_iteration++;
		    }
		  }
		 

	
      }
    

	    
	 	

    
  }

 Free(prec);
 Free(mean);



 
 
}


void prior_predict(int n_levels,
                   double *current_parameter_values,
                   double current_hyperparameter_value,
                   int which,
                   SPARSE_MAT precision,
                   double *mean,
                   double *prec)
{
  int i,ix;
  double xi,sxx;
 

  double sum=0.0;
  SPARSE_VEC *x;
  precision.column_names = (char **) 0;


  *mean = (double) 0.0;
  *prec = precision.offset;
  
  for(i=0;i<n_levels;i++){
      if (i !=which) *mean -= current_parameter_values[i] * precision.offset;
      }
  
  /*one day, let's have a sum element in a structure like v1 */
  

  x=precision.columns[which];
  sxx = sparse_el(x,(which+1)) + precision.offset;
  for(ix=0;ix<x->non_zero;ix++)
    {
      i =  x->rows ? x->rows[ix] : ix;
      xi = x->vals ? x->vals[ix] : 1.0;
      if(i == which){
	*prec += xi;
	  }
      else{
	*mean -= current_parameter_values[i] * xi;
      }
    }
  
      *mean /= *prec;
      *prec *= current_hyperparameter_value;

      
      
}





 
int glm_fit_v(GLM *glm, SPARSE_VEC *x, double *b, double pr_mean, 
	      double pr_prec)
     /*
	Fit variate vector in Bayesian analysis of the generalized linear model.
	
	The routine implements one step of the Iterated Conditional Mode
	algorithm or of the Gibbs sampler. <x> is a sparse vector representation
	of the term to be fitted. On input, <b> must
	hold the value of coefficient currently used in computation of the
	linear predictor. On output, <b> holds the new coefficient value and
	the linear predictor vector for the model is appropriately modified.
	
	Returns 1 if failure, 0 otherwise. Failure is defined as
	(mode_op == 0 or 2)      failure to converge to posterior mode
	(mode_op == 1)           failure of adaptive rejection sampling
	(mode_op == 3)           Metropolis rejection (not really failure!)
	*/
{
  int cy, fail;
  double bs, b0, b1, g, h, pds, pd0, pd1, prop_m, prop_s, alpha,deviance;
  ARSE *e;


   
  bs = b0 = b1 = *b;

  /* Evaluate initial penalized deviance and derivatives */
   

  step_b(glm, x, pr_mean, pr_prec, b0, b1, &pd1, &g, &h);


 
  pds = pd1;
   
  /* Newton-Raphson search for posterior mode */

  for (fail=cy=1; fail && (cy < ncycles); cy++) {

    b0 = b1;
    pd0 = pd1;
    b1 += (g/h);
    

    step_b(glm, x, pr_mean, pr_prec, b0, b1, &pd1, &g, &h);
    fail = ((g*g/h) > dev_tol);
  }

  if (mode_op == 1) {

    /* Gibbs sampling, using adaptive rejection */


    
    e = alloc_arse(arse_max);

    
    for (init_arse(e, b1, -0.5*pd1, g, h);
	 rejection_sample(e);
	 update_arse(e, -0.5*pd1, g)) {
      
      b0 = b1;
      b1 = e->y;
          
      step_b(glm, x, pr_mean, pr_prec, b0, b1, &pd1, &g, 0);
  
      
    }
    /* If accepted on basis of inner hull, correct the linear predictor */

    if (e->y != b1) {
      b0 = b1;
      b1 = e->y;
      step_b(glm, x, pr_mean, pr_prec, b0, b1, 0, 0, 0);
    }
    free_arse(e);
    fail = 0;
  } else if (mode_op > 1) {

    /* Sample from Normal approximation to the posterior */

    prop_m = b1 + g/h;		/* Mean of proposal distribution */
    prop_s = 1.0/sqrt(h);	/* SD of proposal dostribution */
    for (fail=1, cy=0; cy < nmet; cy++){
      b0 = b1;
      b1 = prop_m + prop_s*snd();
      step_b(glm, x, pr_mean, pr_prec, b0, b1, &pd1, 0, 0);
      if (mode_op == 3) {

	/* Metropolis rejection test */

	alpha = exp(-0.5*(pd1 - pds + 
			  h*(bs-b1)*(bs+b1-2.0*prop_m)));
	if ((alpha > 1.0) || (u_random() < alpha)) {
               
	  /* Accept value */

	  bs = b1;
	  pds = pd1;
	  fail = 0;
	}
      } else {
	fail = 0;
      }
    }
    if ((mode_op == 3) && (b1 != bs)) {

      /* Restore last accepted value */

      b0 = b1;
      b1 = bs;
      step_b(glm, x, pr_mean, pr_prec, b0, b1, 0, 0, 0);
    }        
  }
  *b = b1;
  nsamp ++;
  nfail += fail;
  return fail;
}



void step_b(GLM *glm, SPARSE_VEC *x, double pr_mean, double pr_prec,
	    double b_old, double b_new, double *pd, double *g, double *h)
     /*
	Modify linear predictor according to change in one regression coefficient
	Optionally (if relevant pointers are supplied) : 
	Calculate penalized deviance contribution from affected units
	(Penalized deviance is minus twice the log posterior)
	Calculate gradient and (-) curvature of contribution of
	these units to the log posterior 
	*/ 
{
  int i, ix;
  double xi, yvi, lpi, pwi, fvi, devi, gw, hw, dv, dl, wvi, wti;
  double pdo, go;


  /* Likelihood */
  for (ix=0, gw=hw=dv=0.0; ix < x->non_zero; ix++) {
    i = x->rows ? x->rows[ix] : ix;
    xi = x->vals ? x->vals[ix] : 1.0;
    lpi = (glm->lp[i] += xi*(b_new - b_old));
    if (pd) {
      fvi = (*glm->inv_link)(lpi); 
      yvi = glm->yvar[i];
      pwi = glm->prior_weight ? glm->prior_weight[i] : 1.0;
      devi = pwi*(*glm->dev)(yvi, fvi)/glm->scale;
      dv += devi;
      if (g) {
	dl = (*glm->der_link)(fvi);
	wvi = dl*(yvi - fvi);
	wti = pwi/(glm->scale*dl*dl*(*glm->var)(fvi));
	gw += wti*xi*wvi;
	if (h) {
	  hw += wti*xi*xi;
	}
      }
    }
  }
  /* Prior distribution */

  if (pd) {
    wvi = b_new - pr_mean;
    *pd = dv + pr_prec*wvi*wvi;
    if (g) {
      *g = gw - pr_prec*wvi; 
      if (h) {
	*h = hw + pr_prec;
      }
    }
  }
  return;

}


set_mode(int i)
{
  mode_op=i;
}

void print_current_parameter_values(GLM *glm){

  BLOCK *b;
  int i;


  Rprintf("\n");
  Rprintf("\n");
  Rprintf("Current Parameter Values \n");
  Rprintf("======================== \n");
  Rprintf("\n");
  Rprintf("effect_values\n");
  Rprintf("\n");

  b = glm->rhs;

  while(b->next){
  Rprintf("\n");
    for(i=0;i<b->number_of_effects;i++){
      Rprintf("%lf",*(b->effects+i));
    }
  Rprintf("\n");
  b=b->next;
  }

}
