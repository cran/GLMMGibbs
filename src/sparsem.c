
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "sparsev.h"
#include "sparsem.h"
#include "misc.h"
#define SMALL_REAL 0.001
#include "S.h"

     



SPARSE_VEC *sran_design(int g,
                        int number_of_contrasts,
                        int number_of_units,
                        int *group,
                        SPARSE_VEC *contrast)
{
  SPARSE_VEC *res,*contrast_column;
  SPARSE_VEC *column;
  int i,j,factor_value;
  int matrix_ij_value;
  double *values;
 
  
  res =  Calloc (1, SPARSE_VEC);


  values = Calloc (number_of_units, double);
  if (!values) goto noroom; 

  for(i=0;i<number_of_contrasts;i++){
     contrast_column = contrast+i;
     for(j=0;j<number_of_units;j++){
         factor_value = group[j];
         matrix_ij_value = sparse_el(contrast_column,factor_value);
         values[j] = matrix_ij_value;
         }
      
     column = sparse_vec(number_of_units,values);
     (res[i]).n = column->n;
     (res[i]).non_zero = column->non_zero;
     (res[i]).rows = column->rows;
     (res[i]).vals = column->vals;
  }

  return res;

noroom:
  fprintf(stderr, "*** sran_design *** dynamic storage overflow\n");
  exit(1);
}

SPARSE_MAT *sparsemat_sran_design(int number_of_units,
                                 int *group,
                                 SPARSE_MAT contrast)
{
   SPARSE_MAT *res;
   double *values;
   int i,j;   
   SPARSE_VEC *contrast_column;

/* This subroutine constructs the columns of the design matrix
   correspdonding to the product random factor (that is, the
   factor obtained by  \lq\lq multiplying together \rq\rq the
   individual factors in the product (see subroutine
   |product_random_factor| in file \lq\lq v2b.c \rq\rq)) of
   a particular block.

   On entry,  |number_of_units| is the number of observations
   in the data;

   |*group| is a pointer which points at the beginning of a
   \lq column \rq of integers length  |number_of_units| 
   represening the levels of the product random factor.

   |contrast| is a |SPARSE_MAT| object representing the
   contrast matrix for the product random factor.
    
  
   If the random factor has $i$=|contrast.number_of_columns|
   $\times$ $j$=|contrast.number_of_rows| matrix with elements 
   $c_{ij}$ and if the levels of the product random factor are 
   $l_{1}, \ldtos l_{n}$, then the columns in the design matrix
   representing the product random factor are:

   \[ \left(

  
 
*/   
   res = Calloc(1,SPARSE_MAT);
   
   res->number_of_columns = contrast.number_of_columns;
   res->number_of_rows = number_of_units;
   res->offset = 0.0;
   
   
   res->columns = Calloc((res->number_of_columns) ,SPARSE_VEC *);

   if(!(res->columns)) goto noroom;
   
   values =  Calloc(number_of_units, double);
   if (!(values)) goto noroom;
   
   for(i=0;i<res->number_of_columns;i++)
     {
       for(j=0;j<number_of_units;j++)
	 {
	   contrast_column = contrast.columns[i];
	   values[j] = sparse_el(contrast_column,group[j]);
	 }
       
       res->columns[i] = sparse_vec(number_of_units,values);
     }

   Free(values);
   
   res->column_names = Calloc(res->number_of_columns,
			    char *);
   if(contrast.column_names)
     {
       for(i=0;i<res->number_of_columns;i++)
	 {
	   res->column_names[i] = my_strcopy(contrast.column_names[i]);
	 }
     }
   else
     {
       res->column_names = (char **) 0;
     }
   
  
   return res;
   
 noroom:
   fprintf(stderr, "*** sran_design *** dynamic storage overflow\n");
   exit(1);
 }



   
			      

SPARSE_MAT sparsemat_contr_id_fname(int n, int dummy, char *factor_name)
{
  SPARSE_MAT res;
  int i;
  
  res = (sparsemat_id(n));
  res.column_names = level_names(factor_name,1,n);
  res.type = 'i';
  
  return res;
  
}


SPARSE_MAT sparsemat_contr_sum(int n, int dummy)
{
  SPARSE_MAT res;
  int i;
  

  res.number_of_columns = (n-1);
  res.number_of_rows = n;
  res.offset =0.0;
  
  res.columns = Calloc((n-1), SPARSE_VEC *); 
  if (!(res.columns)) goto noroom;
  

  for(i=0;i<(n-1);i++)
    {
      res.columns[i] =  Calloc (1,SPARSE_VEC);
      (res.columns[i])->n = n;
      (res.columns[i])->non_zero = 2;
      (res.columns[i])->rows = Calloc (2, int);
      if (!((res.columns[i])->rows)) goto noroom;
      ((res.columns[i])->rows)[0] = i;
      ((res.columns[i])->rows)[1] = n-1;      
      (res.columns[i])->vals = (double *) Calloc (2, double);
       if (!((res.columns[i])->vals)) goto noroom;
      ((res.columns[i])->vals)[0] = 1;
      ((res.columns[i])->vals)[1] = -1; 
    }
  
  res.column_names = (char **) 0;
  
  return res;
  

noroom:
  fprintf(stderr, "*** sparsemat_contr_sum *** dynamic storage overflow\n");
  exit(1);

}


SPARSE_MAT sparsemat_contr_sum_fname(int n,
                                     int dummy,
				     char *factor_name)
{
  SPARSE_MAT res;
  
  res = sparsemat_contr_sum(n,dummy);
  res.column_names = level_names(factor_name,1,(n-1));
  res.type = 's';
  
  return res;
}





SPARSE_MAT sparsemat_contr_treatment(int n,
				     int zlevel)
{
  SPARSE_MAT res;
  int i;
  
  res.number_of_columns = n-1;
  res.number_of_rows = n;
  res.offset = 0.0;
  
  res.columns = Calloc( (n-1) , SPARSE_VEC *);
  if (!(res.columns)) goto noroom; 
  
  for(i=0;i<(n-1);i++)
    {
      res.columns[i] = Calloc (1,SPARSE_VEC );
      (res.columns[i])->n = n;
      (res.columns[i])->non_zero = 1;
      (res.columns[i])->rows = Calloc (1, int);
      if (!((res.columns[i])->rows)) goto noroom;
      if(i < (zlevel - 1)){
	((res.columns[i])->rows)[0] = (i);
      }
      else{
	((res.columns[i])->rows)[0] = (i+1);
      }
      
      (res.columns[i])->vals = (double *) 0;
      
    }

  return res;
  



noroom:
  fprintf(stderr, "*** sparsemat_contr_treatment *** dynamic storage overflow\n");
  exit(1);

}


SPARSE_MAT sparsemat_contr_treatment_fname(int n,
                                           int zlevel,
					   char *factor_name)
{
  SPARSE_MAT res;
  res = sparsemat_contr_treatment(n,zlevel);
  res.column_names = level_names_zlevel(factor_name,1,n,zlevel);
  res.type = 't';
  
  return res;
}


SPARSE_MAT sparsemat_contr_helmert(int n,int dummy)
{
  SPARSE_MAT res;
  int i,j,k;
  
  
  res.number_of_columns = (n-1);
  res.number_of_rows = n;
  res.offset = (double) 0.0;
  
  res.columns = Calloc((n-1), SPARSE_VEC *);
  for(i=0;i<(n-1);i++)
    { 
      res.columns[i] = Calloc(1,SPARSE_VEC );
      if(!(res.columns[i])) goto noroom;
      (res.columns[i])->n = n;       
      (res.columns[i])->non_zero = (i+2);
      (res.columns[i])->rows = Calloc((i+2), int);
      k=0;
      for(j=0;j<(i+2);j++)
	{
	  ((res.columns[i])->rows)[k] = j;
          k++;
	}	  
      
      (res.columns[i])->vals =  Calloc((i+2), double);
      for(j=0;j<(i+1);j++)
	{
	  ((res.columns[i])->vals)[j] = (double) -1.0;
	}
      ((res.columns[i])->vals)[(i+1)] = i+1;
      
	   
      
      

    }
  res.type='h';
  
  return res;
  

noroom:
  fprintf(stderr, "*** sparsemat_helmert_treatment *** dynamic storage overflow\n");
  exit(1);

}


SPARSE_MAT sparsemat_contr_helmert_fname(int n, int dummy, char *factor_name)
{
  SPARSE_MAT res;
  res = sparsemat_contr_helmert( n, dummy);
  res.column_names = level_names(factor_name,1,(n-1));
  return res;
}




SPARSE_MAT sparsemat_contr_diff2(int n, int dummy)
{
  SPARSE_MAT res;
  int i;
  res.offset = (double) 0.0;
  res.number_of_columns = (n-2);
  res.number_of_rows = n;
  
  res.columns =  Calloc(n-2, SPARSE_VEC*) ;
  if(!(res.columns)) goto noroom;
  
  for(i=0;i<(n-2);i++)
    {
      res.columns[i] = Calloc( 1,SPARSE_VEC);
      if(!(res.columns[i])) goto noroom;
    }
  
  for(i=0;i<(n-2);i++)
    {
      fflush(stdout);    
      res.columns[i]->n = n;
      res.columns[i]->non_zero = 1;
      res.columns[i]->vals = (double *) 0;
      res.columns[i]->rows = Calloc(1,int);
      (res.columns[i]->rows)[0] = (i+1);
    }

   
  return res;

  

noroom:
   fprintf(stderr, "*** sparsemat_contr_diff2 *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_MAT sparsemat_contr_diff2_fname(int n,
				       int dummy,
				       char *factor_name)
{
  SPARSE_MAT res;
  
  res = sparsemat_contr_diff2(n,dummy);
  res.column_names = level_names(factor_name,2,(n-1));
  return res;
  
}

 
SPARSE_VEC *sparsevec_scalar_multiply(double scalar, SPARSE_VEC *vec)
{

    SPARSE_VEC *res,*col;
    int i,j;


   
   res= Calloc(1,SPARSE_VEC);
   if (!res) goto noroom;


   res->n = vec->n;
    res->non_zero = vec->non_zero;
    

    if(!vec->non_zero){
      res->non_zero=0;
      res->rows=(int *) 0;
      res->vals=(double *) 0;
      return res;
    }

    if((!(vec->vals)) && (!(vec->rows)))
      {
	res->rows = (int *) 0;
	res->vals = Calloc(vec->n,double);
	
			
	for(i=0;i<(vec->non_zero);i++)
	  {
	    res->vals[i] = scalar;
	  }
	return res;
	
      }
    
    if(!(vec->vals))
      {
	res->rows = Calloc(vec->non_zero,int);
	res->vals = Calloc(vec->non_zero,double);
	
	
	for(i=0;i<vec->non_zero;i++)
	  {
	    res->rows[i] = vec->rows[i];
            res->vals[i] = scalar;
	    
	  }
	return res;
	
      }
    
    if(!(vec->rows))
      {
	res->rows = (int *) 0;
	res->vals = Calloc(vec->non_zero,double);
	for(i=0;i<vec->non_zero;i++)
	  {
	    res->vals[i] = scalar * vec->vals[i];
	  }
	
	return res;
	
      }
    
    res->vals = Calloc(vec->non_zero,double);
    
 
    res->rows = Calloc (vec->non_zero,int);
        for(i=0;i<vec->non_zero;i++)
      {
	res->rows[i] = vec->rows[i];
	res->vals[i] = scalar * vec->vals[i];
	
      }
    return res;
	
      







noroom:
   fprintf(stderr, "*** sparse_scalar_multiply *** dynamic storage overflow\n");
   exit(1);

  }



SPARSE_MAT sparsemat_scalar_multiply( double scalar,
				     SPARSE_MAT mat1)
{
  SPARSE_MAT res;
  int i;
  
  res.number_of_rows = mat1.number_of_rows;
  res.number_of_columns = mat1.number_of_columns;
  res.offset = scalar * mat1.offset;
  res.columns = Calloc(  mat1.number_of_columns,SPARSE_VEC *);
  
  
  for(i=0;i<mat1.number_of_columns;i++)
    { 
      res.columns[i] = Calloc(1,SPARSE_VEC);
    }
  
  for(i=0;i<mat1.number_of_columns;i++)
    {
      res.columns[i] = sparsevec_scalar_multiply(scalar,mat1.columns[i]);
    }
  
  if(mat1.column_names)
    {
      res.column_names = Calloc(mat1.number_of_columns,char *);
      for(i=0;i<mat1.number_of_columns;i++)
	{
	  res.column_names[i] = my_strcopy(mat1.column_names[i]);
	  
	}
    }
  else
    {
      res.column_names = (char **) 0;
    }
  
  return res;
  
}


     
				     


SPARSE_VEC *sparsevec_elementwise_sum(SPARSE_VEC *mat1, 
				   SPARSE_VEC *mat2,
				   int columns)
{

  SPARSE_VEC *res,*col_vector;
  int i,j,n;
  double *work;
  
  res = Calloc(1,SPARSE_VEC);
  if (!res) goto noroom;
  n = mat1->n; 
  work = Calloc(n, double);
  if (!work) goto noroom;
    
  for(i=0;i<columns;i++)
    {
      
	    
      if (mat1[i].n != n) goto dimerror;
      if (mat2[i].n != n) goto dimerror;
      for(j=0;j<n;j++)
	{
          
	  work[j] = sparse_el((mat1+i),(j+1)) +
	    sparse_el((mat2+i),(j+1));
	}
      
      col_vector = sparse_vec(n,work);

     
      res[i].n = col_vector->n;
      res[i].non_zero = col_vector->non_zero;      
      res[i].rows = col_vector->rows;
      res[i].vals = col_vector->vals;
    }

return res;
  
noroom:
   fprintf(stderr,"*** sparse_elementwise_sum *** dynamic storage overflow\n");
   exit(1);

dimerror:
   fprintf(stderr, "*** sparse_elementwise_sum *** dimension error\n");
   exit(1);
}

SPARSE_VEC *sparsevec_elementwise_product(SPARSE_VEC *mat1,
					 SPARSE_VEC *mat2,
                                         int columns)
{

  SPARSE_VEC *res,*col_vector;
  int i,j,n;
  double *work;
  
  res = Calloc(columns, SPARSE_VEC);
  if (!res) goto noroom;
  n = mat1->n; 
  work = Calloc(n, double);
  if (!work) goto noroom;

  for(i=0;i<columns;i++)
    {
      
	    
      if (mat1[i].n != n) goto dimerror;
      if (mat2[i].n != n) goto dimerror;
      for(j=0;j<n;j++)
	{
          
	  work[j] = sparse_el((mat1+i),(j+1)) *
	    sparse_el((mat2+i),(j+1));
	}
      
      col_vector = sparse_vec(n,work);

     
      res[i].n = col_vector->n;
      res[i].non_zero = col_vector->non_zero;      
      res[i].rows = col_vector->rows;
      res[i].vals = col_vector->vals;
    }
  Free(work);
  
return res;
  
noroom:
   fprintf(stderr,"*** sparse_elementwise_product *** dynamic storage overflow\n");
   exit(1);

dimerror:
   fprintf(stderr, "*** sparse_elementwise_product *** dimension error\n");
   exit(1);
}
      
  
SPARSE_MAT sparsemat_elementwise_sum(SPARSE_MAT mat1,
                                    SPARSE_MAT mat2)
{
   SPARSE_MAT res;
   SPARSE_VEC *work;
   int i;
   work =  Calloc(1,SPARSE_VEC);  



   if (mat1.number_of_rows != mat2.number_of_rows) goto dimerror;
   if (mat1.number_of_columns != mat2.number_of_columns) goto dimerror;
    
   res.number_of_rows = mat1.number_of_rows;
   res.number_of_columns = mat1.number_of_columns;
   res.offset = mat1.offset + mat2.offset;
   res.columns = Calloc(res.number_of_columns,SPARSE_VEC  *);

   for(i=0;i<(res.number_of_columns);i++)
     {
       if((((mat1.columns)[i])->n)!=(res.number_of_rows)) goto dimerror;
       if((((mat2.columns)[i])->n)!=(res.number_of_rows)) goto dimerror;
       work = sparsevec_elementwise_sum(mat1.columns[i],
                                      mat2.columns[i],(int) 1);
       res.columns[i] = work;  
     }

   res.column_names = (char **) 0;
   

return res;


dimerror:

   fprintf(stderr, "*** sparsemat_elementwise_sum *** dimension error\n");
   exit(1); 
}
  
SPARSE_MAT sparsemat_elementwise_product(SPARSE_MAT mat1,
                                         SPARSE_MAT mat2)
{
   SPARSE_MAT res;
   SPARSE_VEC *work;
   int i;
   work = Calloc(1, SPARSE_VEC);  



   if (mat1.number_of_rows != mat2.number_of_rows) goto dimerror;
   if (mat1.number_of_columns != mat2.number_of_columns) goto dimerror;
  
   res.number_of_rows = mat1.number_of_rows;
   res.number_of_columns = mat1.number_of_columns;
   res.offset = mat1.offset + mat2.offset;
   res.columns = Calloc(res.number_of_columns,SPARSE_VEC *);
   
   for(i=0;i<(res.number_of_columns);i++)
     {
       if((((mat1.columns)[i])->n)!=(res.number_of_rows)) goto dimerror;
       if((((mat2.columns)[i])->n)!=(res.number_of_rows)) goto dimerror;
       work = sparsevec_elementwise_product(mat1.columns[i],
                                      mat2.columns[i],(int) 1);
       res.columns[i] = work;  
     }

return res;


dimerror:

   fprintf(stderr, "*** sparsemat_elementwise_product *** dimension error\n");
   exit(1); 
}


SPARSE_MAT sparsemat_column_by_column_product(SPARSE_MAT mat1,
					      SPARSE_MAT mat2)
{
  SPARSE_MAT res;
  int i,j,k;
  
  if (mat1.number_of_rows != mat2.number_of_rows) goto unequal_rows;
  res.number_of_rows    = mat1.number_of_rows;
  res.number_of_columns = mat1.number_of_columns * mat2.number_of_columns;
  res.offset            = mat1.offset * mat2.offset; 
  res.columns           = Calloc( (mat1.number_of_columns * mat2.number_of_columns),
						SPARSE_VEC *);
  if(!res.columns) goto noroom;
  
  k=0;
  
  for(i=0;i<mat1.number_of_columns;i++)
    {
      for(j=0;j<mat2.number_of_columns;j++)
	{
	  if (mat1.columns[i]->n != mat1.number_of_rows) goto mat1_error;
	  if (mat2.columns[j]->n != mat2.number_of_rows) goto mat2_error; 
          res.columns[k] = sparsevec_elementwise_product(mat1.columns[i],
							 mat2.columns[j],1);
	  k++;
	  
	}
    }

if((mat1.column_names) && (mat2.column_names))
{
   res.column_names = string_product(mat1.column_names,
		                             mat2.column_names,
                                     mat1.number_of_columns,
                                     mat2.number_of_columns);
}

else
{
  res.column_names = (char **) 0;
}


return res;  

unequal_rows:
  fprintf(stderr,"*** sparsemat_column_by_column_product ***\n");
  fprintf(stderr,"mat1 has %d rows, mat2 has %d",mat1.columns,mat2.columns); 
  exit(1);
	  

mat1_error:
  fprintf(stderr,"*** sparsemat_column_by_column_product ***\n");
  fprintf(stderr,"mat1.number_of_rows = %d, mat1.columns[%d] = %d\n",
	  mat1.number_of_rows,i,mat1.columns[i]);
  exit(1);

mat2_error:
  fprintf(stderr,"*** sparsemat_column_by_column_product ***\n");
  fprintf(stderr,"mat1.number_of_rows = %d, mat1.columns[%d] = %d\n",
	  mat2.number_of_rows,j,mat2.columns[j]);
  exit(1);
	
noroom:
  fprintf(stderr,"*** sparsemat_column_by_column *** dynamic storage overflow\n");

}


SPARSE_MAT sparsemat_id(int n)
  {
    SPARSE_MAT res;
    int i;

    res.number_of_rows = n;
    res.number_of_columns = n;
    res.offset = (double) 0.0;
    res.columns = Calloc(n, SPARSE_VEC *);
    if(!(res.columns))  goto noroom;
    for (i=0;i<n;i++)
      {
        res.columns[i] = Calloc(1,SPARSE_VEC);
         
        (res.columns[i])->n =n;
        (res.columns[i])->non_zero = 1;
        (res.columns[i])->rows = Calloc(1,int);
        if (!(res.columns[i])) goto noroom;
        ((res.columns[i])->rows)[0] = i;
        (res.columns[i])->vals = ((double *) 0);
      }

     res.column_names = (char **) 0;
     res.type = 'i';
     
   return res;

noroom:

   fprintf(stderr, "*** sparsemat_id *** dimension error\n");
   exit(1); 
   



  }


SPARSE_MAT sparsemat_id_dumchar(int n, char *dum)
{
  SPARSE_MAT res;
  int i;
  
  res = (sparsemat_id(n));
  return res;
  
}


SPARSE_MAT sparsemat_id_fname(int n, char *factor_name)
{
  SPARSE_MAT res;
  int i;
  
  res = (sparsemat_id(n));
  res.column_names = level_names(factor_name,1,n);
  return res;
  
}


void sparsemat_wr(FILE *out, SPARSE_MAT mat)
{
  int i,logical;
  SPARSE_VEC *sparsevp,*sparsevp2;
  
  logical =  iszeroa(mat.offset);
  

  fprintf(out,"number_of_rows=%d\n",mat.number_of_rows);
  fprintf(out,"number_of_columns=%d\n",mat.number_of_columns);
  
  for(i=0;i<mat.number_of_columns;i++)
    {
      sparsevp = (mat.columns)[i];

      if (mat.column_names)
        {
	  fprintf(out,"%s: ",mat.column_names[i]);
	}
      
      if (sparsevp->n != mat.number_of_rows) goto dimerror;
      if(logical)
	{
	  sparse_wr(out,sparsevp,(int) 1);
	}
      else
	{
	  sparsevp2 = sparsevec_scalar_add(sparsevp,mat.offset);
	  sparse_wr(out,sparsevp2,(int) 1);
	}
    }  
      return;
      
    dimerror:
      
      fprintf(stderr, "*** sparsemat_pr *** dimension error\n");

  
      exit(0); 

    
  
}


void sparsemat_pr(SPARSE_MAT mat)
{
  int i,logical;
  SPARSE_VEC *sparsevp,*sparsevp2;
  
  logical =  iszeroa(mat.offset);
  

  printf("number_of_rows=%d\n",mat.number_of_rows);
  printf("number_of_columns=%d\n",mat.number_of_columns);
  printf("\n");
  printf("offset = %lf",mat.offset);
  printf("\n");
  
  for(i=0;i<mat.number_of_columns;i++)
    {
      sparsevp = (mat.columns)[i];

      if (mat.column_names)
        {
	  printf("%s: ",mat.column_names[i]);
	}
      
      if (sparsevp->n != mat.number_of_rows) goto dimerror;

	{
	  sparse_wr(stdout,sparsevp,(int) 1);
	}

    }  
      return;
      
    dimerror:
      
      fprintf(stderr, "*** sparsemat_pr *** dimension error\n");
  sparsemat_dump(mat);
  
      exit(0); 

    
  
}


void sparsemat_pr2(SPARSE_MAT mat)
{
  int i,logical;
  SPARSE_VEC *sparsevp,*sparsevp2;
  
  logical =  iszeroa(mat.offset);
  

  printf("number_of_rows=%d\n",mat.number_of_rows);
  printf("number_of_columns=%d\n",mat.number_of_columns);
  
  for(i=0;i<mat.number_of_columns;i++)
    {
      sparsevp = (mat.columns)[i];

      if (mat.column_names)
        {
	  printf("%s: ",mat.column_names[i]);
	}
      
      if (sparsevp->n != mat.number_of_rows) goto dimerror;
      if(logical)
	{
	  sparse_wr(stdout,sparsevp,(int) 1);
	}
      else
	{
	  sparsevp2 = sparsevec_scalar_add(sparsevp,mat.offset);
	  sparse_wr(stdout,sparsevp2,(int) 1);
	}
    }  
      return;
      
    dimerror:
      
      fprintf(stderr, "*** sparsemat_pr *** dimension error\n");
  sparsemat_dump(mat);
  
      exit(0); 

    
  
}


SPARSE_VEC *sparsevec_scalar_add(SPARSE_VEC *vec, 
                                 double scalar)
{
  SPARSE_VEC *res;
  double *work;
  int i;
  
  work = Calloc (vec->n, double);
  if (!(work)) goto noroom;
  
  for (i=0;i<vec->n;i++)
    {
      work[i] = sparse_el(vec,(i+1)) + scalar;
    }
  
  res = sparse_vec(vec->n,work);
  
  return res; 
  
 noroom:
  
  fprintf(stderr, "*** sparsevec_scalar_add *** dimension error\n");
  exit(1); 
  
}  

SPARSE_MAT sparsemat_scalar_add(SPARSE_MAT mat,
				double scalar)
{
  int i;
  
  SPARSE_MAT res;
  res.number_of_columns = mat.number_of_columns;
  res.number_of_rows    = mat.number_of_rows;
  res.offset = 0.0;
  res.columns = Calloc (mat.number_of_columns,SPARSE_VEC *);
  


  for(i=0;i<mat.number_of_columns;i++)
    {
      res.columns[i] = sparsevec_scalar_add(mat.columns[i],scalar);
    }
      
  
  if(mat.column_names)
    {
      res.column_names = (char **) Calloc(mat.number_of_columns,char *);
      
      for(i=0;i<mat.number_of_columns;i++)
	{
	  res.column_names[i] = my_strcopy(mat.column_names[i]);
	  
	}
    }
  else
    {
      res.column_names = (char **) 0;
    }
  
  return res;
}


SPARSE_VEC *sparsevec_ncopies(SPARSE_VEC *vec, int n)
{
  SPARSE_VEC *res;
  int i,j,position;
  
  res = Calloc(1, SPARSE_VEC);
  res->n = n*(vec->n);
  res->non_zero = n*(vec->non_zero);
  

  if(vec->rows)
    {
      res->rows = Calloc(n*(vec->non_zero),int);
            position=0;
      for(i=0;i<n;i++)
	{
	  for(j=0;j<(vec->non_zero);j++)
	    {
	      (res->rows)[position] = (vec->rows)[j]
		+ i*(vec->n);
	      position++;
	      
	    }
	}
    }
   

  if(vec->vals)
    {
      res->vals = Calloc(n*(vec->non_zero),double);
      
      position=0;
      for(i=0;i<n;i++)
	{
	  for(j=0;j<(vec->non_zero);j++)
	    {
	      (res->vals)[position] = (vec->vals)[j];
	      position++;
	    }
	}
    }
  return res;
  
}         



SPARSE_MAT sparsemat_diff1_Umatrix(int n)
{
  SPARSE_MAT res;
  int i;
  
  res.number_of_rows = n;
  res.number_of_columns = n;
  res.offset = (double) 0.0;


  
  res.columns = Calloc(n, SPARSE_VEC *);
  if (!(res.columns)) goto noroom;

  for(i=0;i<n;i++)
    {
      (res.columns)[i] = Calloc(1, SPARSE_VEC );
      if (!((res.columns)[i])) goto noroom;
    }
  
  
  /* First column */
  ((res.columns)[0])-> n = n;
  ((res.columns)[0])->non_zero = 2;
  res.columns[0]->rows = Calloc(2, int);
  res.columns[0]->vals = Calloc(2, double);  
  *(((res.columns)[0])->rows) = 0;
  *(((res.columns)[0])->vals) = 1.0;
  *(((res.columns)[0])->rows+1) = 1;
  *(((res.columns)[0])->vals+1) = -1.0;

  /* Columns 2 to (n-1) */
  for(i=1;i<(n-1);i++)
    {
      ((res.columns)[i])->n =n;
      ((res.columns)[i])->non_zero =3;
      res.columns[i]->rows = Calloc(3, int) ;
      res.columns[i]->vals = Calloc(3, double) ;      
      *((res.columns[i])->rows) = i-1;
      *((res.columns[i])->vals) = (double) -1.0;
      *((res.columns[i])->rows+1) = i;
      *((res.columns[i])->vals+1) = (double) 2.0;
      *((res.columns[i])->rows+2) = i+1;
      *((res.columns[i])->vals+2) = (double) -1.0;
    }
  
  /* Final column */
  ((res.columns)[(n-1)])-> n = n;
  ((res.columns)[(n-1)])->non_zero = 2;
  res.columns[(n-1)]->rows = Calloc(2,        int);
  res.columns[(n-1)]->vals = Calloc(2, double);
  *(((res.columns)[(n-1)])->rows) = (n-2);
  *(((res.columns)[(n-1)])->vals) = -1.0;
  *(((res.columns)[(n-1)])->rows+1) = (n-1);
  *(((res.columns)[(n-1)])->vals+1) = 1.0;	

  return res;
  
 
 noroom:
  
  fprintf(stderr, "*** sparsemat_diff1_Umatrix *** dimension error\n");
  exit(1); 

}

SPARSE_MAT sparsemat_diff1_Umatrix_dumchar(int n, char *dumchar)
{
  return sparsemat_diff1_Umatrix(n);
}


SPARSE_MAT sparsemat_diff2_Umatrix(int n)
{
  SPARSE_MAT res;
  int i;
  
  res.number_of_rows = n;
  res.number_of_columns = n;
  res.offset = (double) 0.0;
  
  res.columns = Calloc (n, SPARSE_VEC *);
  if (!(res.columns)) goto noroom;
  

  for(i=0;i<n;i++)
    {
      (res.columns)[i] = Calloc(1,SPARSE_VEC );
      if (!((res.columns)[i])) goto noroom;
    }
  
      
  /* First column */     
  ((res.columns)[0])-> n = n;
  ((res.columns)[0])->non_zero = 3;
  ((res.columns)[0])->rows = Calloc(3, int);
  ((res.columns)[0])->vals = Calloc(3, double);
  *(((res.columns)[0])->rows) = 0;
  *(((res.columns)[0])->vals) = 1.0;
  *(((res.columns)[0])->rows+1) = 1;
  *(((res.columns)[0])->vals+1) = -2.0;  
  *(((res.columns)[0])->rows+2) = 2;
  *(((res.columns)[0])->vals+2) = 1.0; 




  /* Second column */
  ((res.columns)[1])-> n = n;
  ((res.columns)[1])->non_zero = 4;
  ((res.columns)[1])->rows = Calloc(4, int);
  ((res.columns)[1])->vals = Calloc(4, double);
  *(((res.columns)[1])->rows) = 0;
  *(((res.columns)[1])->vals) = -2.0;
  *(((res.columns)[1])->rows+1) = 1;
  *(((res.columns)[1])->vals+1) = 5.0; 
  *(((res.columns)[1])->rows+2) = 2;
  *(((res.columns)[1])->vals+2) = -4.0;
  *(((res.columns)[1])->rows+3) = 3;
  *(((res.columns)[1])->vals+3) = 1.0;




  /* Columns 3 to (n-2) */
  for(i=5; i<(n+1); i++)
    {
      ((res.columns)[(i-3)])->n =n;
      ((res.columns)[(i-3)])->non_zero =5;
      ((res.columns)[(i-3)])->rows =  Calloc(5, int);

      ((res.columns)[(i-3)])->vals = Calloc(5, double);

      *(((res.columns)[(i-3)])->rows) = i-5;
      *(((res.columns)[(i-3)])->vals) = 1.0;
      *(((res.columns)[(i-3)])->rows+1) = i-4;
      *(((res.columns)[(i-3)])->vals+1) = -4.0;
      *(((res.columns)[(i-3)])->rows+2) = i-3;
      *(((res.columns)[(i-3)])->vals+2) = 6.0;
      *(((res.columns)[(i-3)])->rows+3) = i-2;
      *(((res.columns)[(i-3)])->vals+3) = -4.0;
      *(((res.columns)[(i-3)])->rows+4) = i-1;
      *(((res.columns)[(i-3)])->vals+4) = 1.0;    
    }
  
  
  /* Penultimate Column */


  ((res.columns)[(n-2)])-> n = n;
  ((res.columns)[(n-2)])->non_zero = 4;
  ((res.columns)[(n-2)])->rows = Calloc(4, int);
  ((res.columns)[(n-2)])->vals = Calloc(4, double);
  *(((res.columns)[(n-2)])->rows) = (n-4);
  *(((res.columns)[(n-2)])->vals) = 1.0;
  *(((res.columns)[(n-2)])->rows+1) = (n-3);
  *(((res.columns)[(n-2)])->vals+1) = -4.0; 
  *(((res.columns)[(n-2)])->rows+2) = (n-2);
  *(((res.columns)[(n-2)])->vals+2) = 5.0;
  *(((res.columns)[(n-2)])->rows+3) = (n-1);
  *(((res.columns)[(n-2)])->vals+3) = -2.0;

  /* Final Column */ 

  ((res.columns)[(n-1)])-> n = n;
  ((res.columns)[(n-1)])->non_zero = 3;
  ((res.columns)[(n-1)])->rows = Calloc(3, int);
  ((res.columns)[(n-1)])->vals = Calloc(3, double);
  *(((res.columns)[(n-1)])->rows) = (n-3);
  *(((res.columns)[(n-1)])->vals) = 1.0;
  *(((res.columns)[(n-1)])->rows+1) = (n-2);
  *(((res.columns)[(n-1)])->vals+1) = -2.0;  
  *(((res.columns)[(n-1)])->rows+2) = (n-1);
  *(((res.columns)[(n-1)])->vals+2)  = 1.0; 

  

  return res;
  



 
 noroom:
  
  fprintf(stderr, "*** sparsemat_diff2_Umatrix *** dimension error\n");
  exit(1); 
}


SPARSE_MAT sparsemat_diff2_Umatrix_dumchar(int n, char *dumchar)
{
  return sparsemat_diff2_Umatrix(n);
}

 
int iszeroa(double x)
{
  int iszero;
  iszero = ((x < SMALL_REAL) && (x > - SMALL_REAL));
  return iszero;
  
}
 
int isone(double x)
{
  int isone;
  isone = iszeroa(x - (double) 1.0);
  return isone;
}  

SPARSE_MAT sparsemat_rd(FILE *infile,
			int number_of_rows,
			int number_of_columns)
{
  SPARSE_MAT res;
  int i,ifv,d,nz,r,col;
  double v;
  char dc;
  
  res.offset = (double) 0.0;
  res.number_of_rows = number_of_rows;
  res.number_of_columns = number_of_columns;
  res.columns = Calloc(number_of_columns, SPARSE_VEC *);
  res.column_names = (char **) 0;
  
  for(col=0;col<number_of_columns;col++)
    {
      res.columns[col] = Calloc(1,SPARSE_VEC );
    }
  


  for(col=0;col<number_of_columns;col++)
    {

      
      (res.columns[col])->n = number_of_rows;
      (res.columns[col])->vals = (double *) 0;
      ifv = 0;
      if (fscanf(infile, " %d %c", &nz, &dc)!=2) goto syntax_error;
      if (dc!=':') goto syntax_error;
      (res.columns[col])->non_zero = nz;
      for (i=0; i<nz; i++) 
	{
	  if (dc=='/') {
	    res.columns[col]->non_zero = i;
	    fprintf(stderr, "*** sparse_rd *** too few elements ");
            goto error;
	  }
	  if(fscanf(infile, " %d %c", &r, &dc)!=2) goto syntax_error;
	  r--;
	  if (!i) {
            (res.columns[col])->rows = Calloc(nz, int);
            if (!res.columns[col]->rows) goto noroom;
	  }
	  res.columns[col]->rows[i] = r;
	  if (r>=number_of_rows || r<0) 
	    {
	      fprintf(stderr, "*** sparse_rd *** bad row number ");
	      goto error;
	    }
	  if (dc=='=') 
	    {
	      if(fscanf(infile,"%c",&dc)!=1) goto syntax_error;
	      if(dc!='=') goto syntax_error;
	      
	      if(fscanf(infile, " %lf %c", &v, &dc)!=2) goto syntax_error;
	      if (!i)
		{
		  ifv = 1;
                  printf("YYY\n");
		  
		  res.columns[col]->vals = Calloc(nz, double);
		  if (!res.columns[col]->vals) goto noroom;
		} 
	      else 
		{
		  if (!ifv) goto syntax_error;
		}
	      res.columns[col]->vals[i] = v;
	    } 

	  else 
	    {
	      if (ifv) goto syntax_error;
	    }
	  if (dc!=',' && dc!='/') goto syntax_error;
	}
      if (nz && dc!='/') 
	{
	  fprintf(stderr, "*** sparse_rd *** too many elements ");
	  goto error;
	}

    }

return res; 

syntax_error:
   fprintf(stderr, "*** sparse_rd *** syntax error ");

error:
   fprintf(stderr, "(column %d)/n", col+1);
   exit(1);

noroom:
   fprintf(stderr, "*** sparse_rd  *** dynamic storage overflow\n");
   exit(1);
 }


SPARSE_MAT *sparsemat_matrix(int number_of_rows,
			     int number_of_columns,
			     double *matf77)
{
  SPARSE_MAT *res;
  int i;
  res = Calloc(1,SPARSE_MAT);
  
  res->number_of_rows = number_of_rows;
  res->number_of_columns = number_of_columns;
  res->offset = (double) 0.0;  
  res->columns = Calloc(number_of_columns,SPARSE_VEC *);
  
  res->column_names = (char **) 0;
  

  for(i=0;i<number_of_columns;i++) 
    {
      res->columns[i] = Calloc( 1,SPARSE_VEC);
    }

  for(i=0;i<number_of_columns;i++)
    {
      res->columns[i] = sparse_vec(number_of_rows, matf77);
      matf77 += number_of_rows;
    }

  
  return res;
  
}




SPARSE_MAT *sparsemat_matrix_colnames(int number_of_rows,
				      int number_of_columns,
				      double *matf77,
				      char **column_names)
{
  SPARSE_MAT *res;
  int i;
  res = Calloc(1,SPARSE_MAT);
  
  res->number_of_rows = number_of_rows;
  res->number_of_columns = number_of_columns;
  res->offset = (double) 0.0;  
  res->columns = Calloc(number_of_columns,SPARSE_VEC *);
  
  res->column_names = (char **) 0;
  

  for(i=0;i<number_of_columns;i++) 
    {
      res->columns[i] = Calloc( 1,SPARSE_VEC);
    }

  for(i=0;i<number_of_columns;i++)
    {
      res->columns[i] = sparse_vec(number_of_rows, matf77);
      matf77 += number_of_rows;
    }

     res->column_names = 
      Calloc(res->number_of_columns, char *);
    for(i=0;i<res->number_of_columns;i++)
      {
	res->column_names[i] = 
	  Calloc (strlen(column_names[i]), char);
	strcpy(res->column_names[i],column_names[i]);
      }

  return res;
  
}


SPARSE_MAT map_U_matrix(SPARSE_MAT adj_matrix)
{
  int i,j;
  int *t;
  SPARSE_MAT U_matrix;

  
  
  U_matrix.number_of_rows = adj_matrix.number_of_rows;
  U_matrix.number_of_columns = adj_matrix.number_of_columns;  
  U_matrix.offset = (double) 0.0;

  
  U_matrix.columns = Calloc( U_matrix.number_of_columns,SPARSE_VEC *);
  
					   
  
  for(i=0;i<U_matrix.number_of_columns;i++)
    {
      U_matrix.columns[i] = Calloc (1,SPARSE_VEC);
    }

  
  for(i=0;i<U_matrix.number_of_columns;i++)
    {
      U_matrix.columns[i]->n = U_matrix.number_of_rows;
      U_matrix.columns[i]->non_zero = (adj_matrix.columns[i]->non_zero) + 1;
      U_matrix.columns[i]->vals =  
	Calloc( U_matrix.columns[i]->non_zero, double);
       U_matrix.columns[i]->rows    =  
	Calloc( U_matrix.columns[i]->non_zero, int);

      j=0;

      
      while((adj_matrix.columns[i]->rows[j] < (i)   && (j < U_matrix.columns[i]->non_zero - 1)))
	{
	  U_matrix.columns[i]->rows[j] = adj_matrix.columns[i]->rows[j];
	  
	  U_matrix.columns[i]->vals[j] = (double) -1.0;
          j++;   
	}

      U_matrix.columns[i]->rows[j] = i;
      U_matrix.columns[i]->vals[j] = (double) (adj_matrix.columns[i]->non_zero);      
      j++;
      while(j<(adj_matrix.columns[i]->non_zero)+1)
	{
	  U_matrix.columns[i]->rows[j] = adj_matrix.columns[i]->rows[(j-1)];
	  U_matrix.columns[i]->vals[j] = (double) -1.0;
	  j++;
	  
	}


    }
  return U_matrix;
  
}

SPARSE_MAT sparsemat_mapfile_Umatrix(int nregions,
                                     char *filename)
{
  FILE *ifp;
  SPARSE_MAT adj_matrix,res;
  
  ifp = fopen(filename,"r");

  
  if (ifp)
  adj_matrix = sparsemat_rd(ifp,
			    nregions,
			    nregions);


  res = map_U_matrix(adj_matrix);
  return res;
}

  



void sparsemat_dump(SPARSE_MAT mat)
{
  int i,j;
  SPARSE_VEC *col;
  
  printf("number_of_rows=%d\n",mat.number_of_rows);
  printf("number_of_columns=%d\n",mat.number_of_columns);  
  
  for(i=0;i<mat.number_of_columns;i++)
    {
      col=mat.columns[i];
      printf("column %d:  \n",i);
      printf("n = %d\n",col->n);
      printf("non_zero = %d\n",col->non_zero);
     
      if(col->rows)
	{
	  
	  for(j=0;j<col->non_zero;j++)
	    {
	      printf("rows[%d] = %d\n",j,col->rows[j]);
	    }

      printf("\n");
	}
      if(col->vals)
	{
	  
      for(j=0;j<col->non_zero;j++)
	{
	  printf("vals[%d] = %lf\n",j,col->vals[j]);
	}
      printf("\n");
    }
    }
  
}


#ifdef FALSE
SPARSE_VEC sparsevec_elementwise_product(SPARSE_VEC *vec1, SPARSE_VEC *vec2)
{
  SPARSE_VEC *res;
  double *data;  
  int n,i;
  
 
  if (!((n=vec1->n) == vec2->n)) goto unequal;
  

  data = Calloc(n, double);
  
  if(!(data)) goto noroom;



  for(i=0;i<n;i++)
    {
      data[i] = sparse_el(vec1,(i+1)) *
	sparse_el(vec2,(i+1));
    }

  res = sparse_vec(n,data);
  
  Free(data);
  
  


unequal:
    printf("*** sparsevec_elementwise_product *** vectors of unequal length"); 
  exit(0);
  



noroom:
   printf("*** sparsevec_elementwise_product *** dynamic storage overflow\n");
   exit(1);
}

#endif




SPARSE_MAT sparsemat_zoff_direct_product(SPARSE_MAT mat1,
					 SPARSE_MAT mat2)
{
  SPARSE_MAT res;
  int i,j,k;
  
  res.offset = (double) 0.0;
  res.number_of_columns = (mat1.number_of_columns *
			   mat2.number_of_columns);
  res.number_of_rows    = (mat1.number_of_rows *
			   mat2.number_of_rows);
  res.columns = Calloc (res.number_of_columns,SPARSE_VEC *);
  

  k=0;
  
  for(i=0;i<mat1.number_of_columns;i++)
    {
      for(j=0;j<mat2.number_of_columns;j++)
       {
	  
	  res.columns[k] = Calloc(1,SPARSE_VEC);
	  res.columns[k] = sparse_kpv(mat1.columns[i],
				      mat2.columns[j]);
	  
	  k++;
	}
    }
  
  return res;
}

SPARSE_MAT sparsemat_direct_product(SPARSE_MAT mat1,
				    SPARSE_MAT mat2)
{
  SPARSE_MAT res;
  SPARSE_MAT temp1,temp2,temp3,temp4;
  
  temp1 = sparsemat_ones(mat1.number_of_columns,
			 mat1.number_of_rows);
  temp2 = sparsemat_scalar_multiply(mat1.offset,
				    temp1);
  temp3 = sparsemat_zoff_direct_product(temp2,
					mat2);
  sparsemat_free(temp1);
  sparsemat_free(temp2);
  temp1 = sparsemat_ones(mat2.number_of_columns,
			 mat2.number_of_rows);
  temp2 = sparsemat_scalar_multiply(mat2.offset,
				    temp1);
  sparsemat_free(temp1);
  
  temp1 = sparsemat_zoff_direct_product(mat1,
					temp2);
  sparsemat_free(temp2);
  temp2 = sparsemat_elementwise_sum(temp1,temp3);
  sparsemat_free(temp1);
  sparsemat_free(temp3);
  temp1 = sparsemat_zoff_direct_product(mat1,
					mat2);
  res = sparsemat_elementwise_sum(temp1,temp2);
  res.offset = (mat1.offset *
		mat2.offset);
  sparsemat_free(temp1);
  sparsemat_free(temp2);

  
  return res;
}

  

SPARSE_MAT sparsemat_ones(int nrows, int ncols)
{
  SPARSE_MAT res;
  int i;
  
  res.number_of_rows = nrows;
  res.number_of_columns = ncols;
  res.offset = (double) 0.0;
  res.columns = Calloc(ncols,SPARSE_VEC *);
  
			
  for(i=0;i<ncols;i++)
    {
      res.columns[i] = Calloc(1,SPARSE_VEC);
    } 

  for(i=0;i<ncols;i++)
    {
      res.columns[i]->n = nrows;
      res.columns[i]->non_zero = nrows;
      res.columns[i]->rows = (int *) 0;
      res.columns[i]->vals = (double *) 0;
    }
  res.column_names = (char **) 0;
  
  return res;
  
noroom:
   printf("*** sparsevec_ones *** dynamic storage overflow\n");
   exit(1);

}


int sparsemat_is_identity(SPARSE_MAT mat)
{
  int i;
  SPARSE_VEC *vec;

  if(mat.number_of_columns != mat.number_of_rows) return 0;
  for(i=0;i<mat.number_of_columns;i++)
    {
      vec = mat.columns[i];
      if (vec->n != mat.number_of_rows) goto dimerror;
      if (vec->non_zero != (int) 1) return (int) 0;
      if (vec->rows[0] != i) return (int) 0;
      if (vec->vals && (!(isone(vec->vals[0])))) return (int) 0;
    }
  return (int) 1;

 dimerror:
  
  
  fprintf(stderr, "*** sparsemat_is_identity *** dimension error\n");
  exit(1); 

}


double sparsemat_qform(SPARSE_VEC *vec1,
		       SPARSE_VEC *vec2,
                       SPARSE_MAT mat)
{
  double res=0.0;
  double add;
  
  int i,ir,j,jr;
  double iv,jv;

  if(sparsemat_is_identity(mat)) return sparse_ip(vec1,vec2,(double *) 0);
  

  for(i=0;i<vec1->non_zero;i++)
    {
      ir = vec1->rows ? vec1->rows[i] : i ;
      iv = vec1->vals ? vec1->vals[i] : (double) 1.0;
      
      for(j=0;j<vec2->non_zero;j++)
	{
	  jr = vec2->rows ? vec2->rows[j] : j ;
	  jv = vec2->vals ? vec2->vals[j] : (double) 1.0;
	  res += iv*jv*(mat.offset+sparse_el(mat.columns[ir],(jr+1)));

	}
    }
  return res;
  
}



double sparsemat_dqform(double *vec1,
		       double *vec2,
			SPARSE_MAT mat)
{
  double res=0.0;
  int i,j,ir;
  double iv;
  
  SPARSE_VEC *vec;
  double sum1=0.0,sum2=0.0;
  
  if(!(iszeroa(mat.offset)))
     {
       for(i=0;i<mat.number_of_columns;i++) sum1 += vec1[i];
       for(i=0;i<mat.number_of_columns;i++) sum2 += vec1[i];
       res = mat.offset * sum1 * sum2;
     }
  
  

  for(i=0;i<mat.number_of_columns;i++){
    vec=mat.columns[i];
    for(j=0;j<vec->non_zero;j++){
      ir =  vec->rows ? vec->rows[j] : j ;
      iv = vec->vals ? vec->vals[j] : (double) 1.0;
      res += iv*vec1[ir]*vec2[i];
    }
  }
  return res;
}




  
SPARSE_MAT sparsemat_btrans_a_b(SPARSE_MAT matb,
				SPARSE_MAT mata)
{
  SPARSE_MAT res;
  double *data;
  int i,j,n;
  SPARSE_VEC  *v;
  FILE *fp;
 
 
  
  n = matb.number_of_columns;

  
  res.number_of_rows = n;
  res.number_of_columns = n;
  res.offset = (double) 0.0;
  res.columns = Calloc(n,SPARSE_VEC *);
  
  for(i=0;i<n;i++)
    {
      res.columns[i] = (SPARSE_VEC *) Calloc(1,SPARSE_VEC);
    } 

  data = Calloc(n, double);
  
  
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	{
      	  data[j] = sparsemat_qform(matb.columns[i],matb.columns[j],mata);
	}
      v =sparse_vec(n,data);
      sparse_cp(v,res.columns[i]);
      
    }
  res.column_names = (char **) 0;

  return res;

  
}


SPARSE_MAT s_matrix(SPARSE_MAT contrast_matrix,
		    SPARSE_MAT u_matrix)
{
  SPARSE_MAT res,lots_ones,work1,a,work3,work4,work5,work6;
  SPARSE_VEC *onest;
  double b; 
  FILE *outfp;

  outfp = fopen("s_debug","w");
  
  lots_ones = sparsemat_ones(u_matrix.number_of_rows,u_matrix.number_of_rows);

  res.offset = (double) 0.0 ;


  
	  

  
  if(sparsemat_is_identity(contrast_matrix))
    {



     
      res = u_matrix;
      res.column_names = (char **) 0;
      return res;

    }

      work1=sparsemat_btrans_a_b(u_matrix,lots_ones);
      
  a=sparsemat_btrans_a_b(contrast_matrix,work1);
  onest=sparse_one(u_matrix.number_of_rows);



  fflush(stdout);

  
  b=sparsemat_qform(onest,onest,u_matrix);



  fflush(stdout);
 
  
  
  work3 = sparsemat_btrans_a_b(contrast_matrix,u_matrix);
 
  
  work3.column_names = (char **) 0;
  

  fflush(stdout);
  
  if(iszeroa(b)) {return work3;}
  else
   
    { 

  fflush(stdout);
      work4 = sparsemat_scalar_multiply((-1.0/b),a);
  fflush(stdout);
     work5 = sparsemat_elementwise_sum(work3,work4);
  fflush(stdout);

     work6 = sparsemat_scale(work5,5.0,-5.0,0.1);
  fflush(stdout);
     
        return work6;}   

  
}

void sparsemat_free(SPARSE_MAT mat)
{
  int i;
  
  for(i=0;i<mat.number_of_columns;i++)
      {
      sparse_fr(mat.columns[i],1);
    }
      
  Free(mat.columns);

  if(mat.column_names){
  for(i=0;i<mat.number_of_columns;i++)
      {
      Free(mat.column_names[i]);
    }
      
  Free(mat.column_names);
  
  }
}



void sparsemat_freep(SPARSE_MAT *mat)
{
  int i;
  
  if(mat->columns){
    
  for(i=0;i<mat->number_of_columns;i++)
      {
      sparse_fr(mat->columns[i],1);
    }
      
  Free(mat->columns);
  }
  
  if(mat->column_names){
  for(i=0;i<mat->number_of_columns;i++)
      {
      Free(mat->column_names[i]);
    }
      
  Free(mat->column_names);
  }
  
  Free(mat);
  


}



void sparsemat_copy(SPARSE_MAT *from, SPARSE_MAT *to)
{
  int i;

  to->number_of_rows = from->number_of_rows;

  to->number_of_columns = from->number_of_columns;
  

  to->offset = from->offset;

  to->columns = Calloc (from->number_of_columns, SPARSE_VEC *);
  
					 
  
  if(!(to->columns)) goto noroom;

  for(i=0;i<from->number_of_columns;i++)
    {
      to->columns[i] =  Calloc(1, SPARSE_VEC);
  if(!(to->columns[i])) goto noroom;
    }
  
	

  for(i=0;i<from->number_of_columns;i++)
    {
      sparse_cp(from->columns[i],
		to->columns[i]);
    }

if(from->column_names){  

  to->column_names = Calloc (from->number_of_columns,char *) ;
  
			
  for(i=0;i<from->number_of_columns;i++)
    {
      to->column_names[i] = my_strcopy(from->column_names[i]);
    }
 }
  
 else
{
  to->column_names = (char **) 0 ;
  
}  

return;

noroom:   
   fprintf(stderr, "*** sparsemat_copy *** dynamic storage overflow\n");
   exit(1);

}



void sparsemat_n_product(SPARSE_MAT *terms,
			 int n_terms,
			 SPARSE_MAT *product)
{
  SPARSE_MAT running_product,old;
  int i;
  
  

  sparsemat_copy(terms,&running_product);
  for(i=1;i<n_terms;i++)
    {

      running_product = sparsemat_direct_product(running_product,
						      terms[i]);
    }

  sparsemat_copy(&running_product,product);
}


      
SPARSE_MAT sparsemat_scale(SPARSE_MAT mat,
			   double maximum,
			   double minimum,
			   double tolerance)
{
  SPARSE_MAT res;
  
  int n,i,j,k;
  int *times;
  SPARSE_VEC *vec;
  
  double *values;
  double element;
  int maxt,position=0;
  double offset;
  
  n = (int) ((maximum- minimum)/tolerance) + 3;
  values = Calloc(n, double);
  times =  Calloc(n, int);
  for(i=0;i<n;i++){
    values[i] = (minimum + (i-1) * tolerance);
    times[i] = 0;
  }

  for(i=0;i<mat.number_of_columns;i++){
    vec = mat.columns[i];
    for(j=0;j<mat.number_of_rows;j++){
      element = sparse_el(vec,(j+1));
      k=0;
      while(values[(k+1)] < element ) k++;
      (times[k])++;
    }
  }
  
  maxt=0;
  for(i=0;i<n;i++){
    if(times[i] > maxt){
      maxt=times[i];
      position=i;
    }
  }
  


  for(i=0;i<mat.number_of_columns;i++){
    vec = mat.columns[i];
    for(j=0;j<mat.number_of_rows;j++){
      element = sparse_el(vec,(j+1));
      if((element > values[position]) && element < values[(position+1)])
	 {
	   offset=element;
	   goto jumphere;
	 }
    }
  }
  
  


jumphere: 
  res = sparsemat_scalar_add(mat,-(offset));
  res.offset = offset; 
  
return res ;

}


void matrix_elements(SPARSE_MAT mat, double *res)
{
  int i;
  double *x;
  
  
  for(i=0;i<mat.number_of_columns;i++){
    elements(mat.columns[i],res+(mat.number_of_rows*i));
    
  }

  return;
  

  
}


SPARSE_MAT smatrix_uidentity(SPARSE_MAT contrast_matrix)
{
  SPARSE_MAT res;
  int g;
  int i;
  double *data;
  
  
  g=contrast_matrix.number_of_columns;
  
  if (contrast_matrix.type == 'i'){
    res=sparsemat_id(g);
  }
  

  if (contrast_matrix.type == 's'){
    res=sparsemat_id(g-1);
    res.offset = 1.0;
    
  }

  if (contrast_matrix.type == 't'){
    res=sparsemat_id(g-1);
    res.offset = -(1.0/((double) g));
    
  }

  if (contrast_matrix.type == 'h'){
    data = Calloc(g-1,double);
    for(i=0;i<g-1;i++){
      data[i]=(double) (i+1)*(i+2);
      
    }
    res=sparsemat_diag(data,g-1);
    Free(data);    
  }
  
  
  res.column_names == (char **) 0;
  

  

  
  return(res);
  
}


SPARSE_MAT sparsemat_diag(double *data,int n)
  {
    SPARSE_MAT res;
    int i;

    res.number_of_rows = n;
    res.number_of_columns = n;
    res.offset = (double) 0.0;
    res.columns = Calloc(n, SPARSE_VEC *);
    if(!(res.columns))  goto noroom;
    for (i=0;i<n;i++)
      {
        res.columns[i] = Calloc(1,SPARSE_VEC);
         
        (res.columns[i])->n =n;
        (res.columns[i])->non_zero = 1;
        (res.columns[i])->rows = Calloc(1,int);
        if (!(res.columns[i])) goto noroom;
        ((res.columns[i])->rows)[0] = i;
        (res.columns[i])->vals = (Calloc(1,double));
        (res.columns[i])->vals[0] = data[i];
      }

     res.column_names = (char **) 0;
     
   return res;

noroom:

   fprintf(stderr, "*** sparsemat_diag *** dimension error\n");
   exit(1); 
   



  }




SPARSE_MAT *submatrix(SPARSE_MAT *matrix,
		      int first_row,
		      int last_row,
		      int first_column,
		      int last_column)
{
  int i;
  SPARSE_MAT *res;
  SPARSE_VEC *y,*x;
  res=Calloc(1,SPARSE_MAT);
  
  res->number_of_columns=last_column - first_column +1;
  res->number_of_rows=last_row - first_row +1;
  res->columns = Calloc(res->number_of_columns,SPARSE_VEC*);

  for(i=first_column-1;i<last_column;i++){
      y=matrix->columns[i];
      x=select_elements(y,first_row,last_row);
      res->columns[(i-(first_column-1))] = Calloc(1,SPARSE_VEC);
      sparse_cp(x,res->columns[(i-(first_column-1))]);
      sparse_fr(x,1);
    }
   
  return res;
  

    
}



void  spmatrix(SPARSE_MAT *contrast,
	      SPARSE_MAT *u_matrix,
	      SPARSE_MAT *res)
{



  int g;
  char *name;
  int i;
  SPARSE_VEC *x,*y;
  SPARSE_MAT *mat1;
  double l;
  int n;
  double *last_row,*last_column;
  double *all_elements;
  int j,k;
  SPARSE_MAT *res1;
  
  
  
   
  g = contrast->number_of_columns;


  if(contrast->type == 'i'){

  res->number_of_columns=g;
  res->number_of_rows=g;
  res->columns = Calloc(g,SPARSE_VEC*);
  for(i=0;i<g;i++){
    res->columns[i] = Calloc(1,SPARSE_VEC);
    (res->columns[i])->n =g;
    (res->columns[i])->non_zero = 1;
    (res->columns[i])->rows = Calloc(1,int);
    ((res->columns[i])->rows)[0] = i;
    (res->columns[i])->vals = ((double *) 0);
  }
  res->offset = 0.0;
  }


  if(contrast->type == 'h'){

  res->number_of_columns=g;
  res->number_of_rows=g;
  res->columns = Calloc(g,SPARSE_VEC*);
  for(i=0;i<g;i++){
    res->columns[i] = Calloc(1,SPARSE_VEC);
    (res->columns[i])->n =g;
    (res->columns[i])->non_zero = 1;
    (res->columns[i])->rows = Calloc(1,int);
    ((res->columns[i])->rows)[0] = i;
    (res->columns[i])->vals = (Calloc(1,double));
    ((res->columns[i])->vals)[0] = (double) (i+1)*(i+2);
  }
  res->offset = 0.0;
  }
  

  if (contrast->type == 't'){
    g = contrast->number_of_columns;
    
    res->number_of_columns=g;
    res->number_of_rows=g;

    res->columns = Calloc(g,SPARSE_VEC*);

    
    for(i=1;i<(g+1);i++){
      y=u_matrix->columns[i];
   
      
      x=select_elements(y,2,(g+1));
          
      res->columns[(i-1)] = Calloc(1,SPARSE_VEC);
      
      sparse_cp(x,res->columns[(i-1)]);

          
    }
    
    
    if(u_matrix->type == 'i'){
      res->offset  = -(1.0/(g+1));
    }
    else{
      res->offset  = 0.0;
    }
    


   
  }
  
    if(contrast->type == 's'){

      if (u_matrix->type == 'i'){
	res->number_of_columns = g;
	res->number_of_rows = g;
	res->columns = Calloc(g, SPARSE_VEC *);
	for(i=0;i<g;i++){
	  res->columns[i] = Calloc(1,SPARSE_VEC);
	  (res->columns[i])->n =g;
	  (res->columns[i])->non_zero = 1;
	  (res->columns[i])->rows = Calloc(1,int);
	  ((res->columns[i])->rows)[0] = i;
	  (res->columns[i])->vals = ((double *) 0);
	}
	res->offset = 1.0;
	
      }
      else{
	


      n=u_matrix->number_of_columns;
      k=0;      
      mat1 = submatrix(u_matrix,1,(n-1),1,(n-1));
      


      last_row = Calloc(n,double);
      last_column = Calloc(n,double);
      all_elements = Calloc((n-1)*(n-1),double);
      

      for(i=0;i<n;i++){
	last_row[i]=last_value(u_matrix->columns[i]);
	
      }

      elements(u_matrix->columns[(n-1)],last_column);

      matrix_elements(*mat1,all_elements);
      
     

      sparsemat_freep(mat1);
      
      



      
      for(i=0;i<(n-1);i++){
	for(j=0;j<(n-1);j++){
	  all_elements[k] = all_elements[k]
	  -last_row[i]
	  -last_row[j];

	  k++;
	  
	}
      }
      
      res1 = sparsemat_matrix((n-1),(n-1),all_elements);      
      sparsemat_copy(res1,res);
      
      res->offset= last_row[(n-1)];
      
      Free(all_elements);
      Free(last_column);
      Free(last_row);
      }
      
      
    }
    



}







