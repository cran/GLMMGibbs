/* sparsev.c */

/* David Clayton, March 1992 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "misc.h"
#include "sparsev.h"
#include "S.h"

SPARSE_VEC *sparse_id(int n)
/*
   Create n*n identity matrix as array of sparse vectors.
*/
{
   SPARSE_VEC *res;
   int i;

   res = Calloc(n, SPARSE_VEC);
   if (!res) goto noroom;
   for (i=0; i<n; i++) {
      res[i].n = n;
      res[i].non_zero = 1;
      res[i].vals = (double *) 0;
      if (!(res[i].rows = Calloc(1,int))) goto noroom;
      *(res[i].rows) = i;
   }
   return res;

noroom:
   fprintf(stderr, "*** sparse_id *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_VEC *sparse_one(int n)
/*
   Create sparse unit vector, length <n>
*/
{
   SPARSE_VEC *res;

   res = Calloc(1,SPARSE_VEC);
   if (!res) goto noroom;
   res->n = n;
   res->non_zero = n;
   res->rows = (int *) 0;
   res->vals = (double *) 0;
   return res;

noroom:
   fprintf(stderr, "*** sparse_one  *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_VEC *sparse_ref(int n, double *vector)
/*
   Create a SPARSE_VEC reference to a standard vector length n.
*/
{
   SPARSE_VEC *res;

   if(!(res = Calloc(1,SPARSE_VEC))) goto noroom;
   res->rows = (int *) 0;
   res->vals = vector;
   res->n = res->non_zero = n;
   return res;

noroom:
   fprintf(stderr, "*** sparse_ref *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_VEC *sparse_vec(int n, double *vector) {
/*
   Returns a pointer to a SPARSE_VEC copy of a standard vector of length <n>
*/
   int i, nz, not1;
   double v;
   SPARSE_VEC *res;

   res = Calloc(1,SPARSE_VEC);
   if (!res) goto noroom;
   for (i=nz=not1=0; i<n; i++)
      if ((v=vector[i]) != 0.0) {
         nz++;
         if (v != 1.0) not1 = 1;
      }
   res->n = n;
   res->non_zero = nz;
   res->rows = (int *) 0;
   res->vals = (double *) 0;
   if (!nz) return res;
   if (nz<n)
     {
       
 res->rows = Calloc(nz, int);
   if (!res->rows) goto noroom;
}
   
   if (not1) {
      res->vals = Calloc(nz, double);
      if (!res->vals) goto noroom;
   } 
   for (i=nz=0; i<n; i++) 
      if ((v=vector[i]) != 0.0) {
         if (res->rows) res->rows[nz] = i;
         if (not1) res->vals[nz] = v;
         nz++;
      }
   return res;

noroom:
   fprintf(stderr, "*** sparse_vec *** dynamic storage overflow\n");
   exit(1);
}


SPARSE_VEC *sparse_mat(int nrow, int ncol, double *matf77) {
/*
   <matf77> is a (pointer to) a Fortran 77 style two-dimensional array
   with <nrow> rows and <ncol> columns. That is, the array is stored as
   a 1-dimensional array with row subscript varying fastest.

   The function returns a pointer to a copy stored as a sparse matrix,
   more specifically as an array of column vectors stored as SPARSE_VECs.
*/
   SPARSE_VEC *res, *r, *col;
   int i;

   r = res = Calloc(ncol, SPARSE_VEC);
   if (!res) goto noroom;

   for (i=0; i<ncol; i++, r++, matf77+=nrow) {
      col = sparse_vec(nrow, matf77);
      /* Probably the next 4 lines can be replaced by:   r=col  */
      r->n = col->n;
      r->non_zero = col->non_zero;
      r->rows = col->rows;
      r->vals = col->vals;
      free(col);
   }
   return res;

noroom:
   fprintf(stderr, "*** sparse_mat *** dynamic storage overflow\n");
   exit(1);
}
 
SPARSE_VEC *sparse_fl(int n, int nl, int x[])
/*
   Take a vector <x> of integer factor levels coded 1...nl, expand
   to <nl> sparse vectors of dummy variables.
*/
{
   SPARSE_VEC *res, *col;
   int i, j, nz, *row;

   res = col = Calloc(nl, SPARSE_VEC);
   if (!res) goto noroom;

   for (j=0; j<nl; j++, col++) {
      col->n = n;
      col->vals = (double *) 0;
      for (i=nz=0; i<n; i++) if ((x[i]-1)==j) nz++;
      col->non_zero = nz;
      if (nz){
         col->rows = row = Calloc(nz, int);
         if (!row) goto noroom;
         for (i=0; i<n; i++) if ((x[i]-1)==j) *(row++) = i;
      } else {
         fprintf(stderr,"*** sparse_fl *** (warning) empty column %d\n",j);
         col->rows = (int *) 0;
      }
   }
   return res;

noroom:
   fprintf(stderr, "*** sparse_fl  *** dynamic storage overflow\n");
   fprintf(stderr, "(creating %d*%d factor design matrix)\n",n, nl);
   exit(1);
}

SPARSE_VEC *sparse_gl(int n, int nl, int repeat)
/*
   Equivalent to generation of factor vector by GLIM gl() function,
   followed by <sparse_fl> to expand to sparse vectors of dummy
   variables. Thus, the routine creates <nl> sparse vectors.
*/
{
}

SPARSE_VEC *sparse_rd(FILE *infile, int n, int m)
/*
   Read sparse matrix consisting of <m> sparse vectors of nominal
   length <n>. Each vector is input in  the following format:
   <non_zero>: <row> (=<val>), <row> (=<val>), ...   /

   Rows are numbered 1:n
*/

/*
   Changed by jonm 30/6/95 to work in case of  non-zero 
   elements <> 1.0.
*/ 
{
   SPARSE_VEC *res;
   char dc;
   int nz, r, ifv, i, col;
   double v;

   res = Calloc(m, SPARSE_VEC);
   if (!res) goto noroom;
   for (col=0; col<m; col++) {
      res[col].n = n;
      ifv = 0;
      if (fscanf(infile, " %d %c", &nz, &dc)!=2) goto syntax_error;
      if (dc!=':') goto syntax_error;
      res[col].non_zero = nz;
      for (i=0; i<nz; i++) {
         if (dc=='/') {
            res[col].non_zero = i;
            fprintf(stderr, "*** sparse_rd *** too few elements ");
            goto error;
         }
         if(fscanf(infile, " %d %c", &r, &dc)!=2) goto syntax_error;
         r--;
         if (!i) {
            res[col].rows = Calloc(nz, int);
            if (!res[col].rows) goto noroom;
         }
         res[col].rows[i] = r;
         if (r>=n || r<0) {
            fprintf(stderr, "*** sparse_rd *** bad row number ");
            goto error;
         }
         if (dc=='=') {
	   if(fscanf(infile,"%c",&dc)!=1) goto syntax_error;
	   if(dc!='=') goto syntax_error;
	   
            if(fscanf(infile, " %lf %c", &v, &dc)!=2) goto syntax_error;
            if (!i) {
               ifv = 1;
               res[col].vals = Calloc(nz, double);
               if (!res[col].vals) goto noroom;
            } else {
               if (!ifv) goto syntax_error;
            }
            res[col].vals[i] = v;
         } else {
            if (ifv) goto syntax_error;
         }
         if (dc!=',' && dc!='/') goto syntax_error;
      }
      if (nz && dc!='/') {
         fprintf(stderr, "*** sparse_rd *** too many elements ");
         goto error;
      }
   }
   return res;

syntax_error:
   fprintf(stderr, "*** sparse_rd *** syntax error ");

error:
   fprintf(stderr, "(column %d)/n", col+1);
   return 0;

noroom:
   fprintf(stderr, "*** sparse_rd  *** dynamic storage overflow\n");
   exit(1);
}

void sparse_wr(FILE *out, SPARSE_VEC *mat, int cols)   
/*
   Write out a sparse matrix in the above format
*/
{
   int i, rows;

   rows = mat->n;
   for(;cols; cols--, mat++){
      if (mat->n!=rows) fprintf(stderr,
         "*** sparse_wr *** (warning) unequal column lengths\n");
      fprintf(out, "%d:", mat->non_zero);
      for (i=0; i<mat->non_zero; i++){
         if (i) fprintf(out, ",");
         if (mat->rows)
           fprintf(out, "%d", 1+mat->rows[i]);
         if (mat->rows && mat->vals)
           fprintf(out, "=");
         if (mat->vals) 
           fprintf(out, "=%lf", mat->vals[i]);
      }
      fprintf(out, "/\n");
   }
}

SPARSE_VEC *diff1_smooth(int n)
/* Matrix of adjacency weights for 1st difference
   (linear interpolation)smoothing */
{
   SPARSE_VEC *res, *col;
   int i;

   res = col = Calloc(n, SPARSE_VEC);
   if (!res) goto noroom;

   /* First column */
   col->n = n;
   col->non_zero =1;
   if (!(col->rows = Calloc(1, int))) goto noroom;
   col->vals = (double *) 0;
   *(col->rows) = 1;

   /* Columns 2 to (n-1) */
   for (i=2; i<n; i++){
      col++;
      col->n = n;
      col->non_zero =2;
      if (!(col->rows = Calloc(2, int))) goto noroom;
      col->vals = (double *) 0;
      *(col->rows) = i-2;
      *(col->rows+1) = i;
   }

   /* Column (n-1) */
   col++;
   col->n = n;
   col->non_zero =1;
   if (!(col->rows = Calloc(1, int))) goto noroom;
   col->vals = (double *) 0;
   *(col->rows) = n-2;

   return res;

noroom:
   fprintf(stderr, "*** diff1_smooth  *** dynamic storage overflow\n");
   exit(1);

}

SPARSE_VEC *diff2_smooth(int n)
/* Matrix of adjacency weights for 2nd difference
   (cubic interpolation) smoothing */
{
   SPARSE_VEC *res, *col;
   int i;

   res = col = Calloc(n, SPARSE_VEC);
   if (!res) goto noroom;

   /* First column */
   col->n = n;
   col->non_zero =2;
   if (!(col->rows = Calloc(2, int))) goto noroom;
   if (!(col->vals = Calloc(2, double))) goto noroom;
   *(col->rows) = 1;
   *(col->vals) = 2.0;
   *(col->rows+1) = 2;
   *(col->vals+1) = -1.0;

   /* Second column */
   col++;
   col->n = n;
   col->non_zero =3;
   if (!(col->rows = Calloc(3, int))) goto noroom;
   if (!(col->vals = Calloc(3, double))) goto noroom;
   *(col->rows) = 0;
   *(col->vals) = 2.0;
   *(col->rows+1) = 2;
   *(col->vals+1) =  4.0;
   *(col->rows+2) = 3;
   *(col->vals+2) = -1.0;

   /* Columns 3 to (n-2) */
   for (i=4; i<n; i++){
      col++;
      col->n = n;
      col->non_zero =4;
      if (!(col->rows = Calloc(4, int))) goto noroom;
      if (!(col->vals = Calloc(4, double))) goto noroom;
      *(col->rows) = i-4;
      *(col->vals) = -1.0;
      *(col->rows+1) = i-3;
      *(col->vals+1) = 4.0;
      *(col->rows+2) = i-1;
      *(col->vals+2) = 4.0;
      *(col->rows+3) = i;
      *(col->vals+3) = -1.0;
   }

   /* Penultimate column */
   col++;
   col->n = n;
   col->non_zero =3;
   if (!(col->rows = Calloc(3, int))) goto noroom;
   if (!(col->vals = Calloc(3, double))) goto noroom;
   *(col->rows) = n-4;
   *(col->vals) = -1.0;
   *(col->rows+1) = n-3;
   *(col->vals+1) =  4.0;
   *(col->rows+2) = n-1;
   *(col->vals+2) = 2.0;

   /* Last column */
   col++;
   col->n = n;
   col->non_zero =2;
   if (!(col->rows = Calloc(2, int))) goto noroom;
   if (!(col->vals = Calloc(2, double))) goto noroom;
   *(col->rows) = n-3;
   *(col->vals) = -1.0;
   *(col->rows+1) = n-2;
   *(col->vals+1) = 2.0;

   return res;

noroom:
   fprintf(stderr, "*** diff2_smooth  *** dynamic storage overflow\n");
   exit(1);

}

int is_symmetric(SPARSE_VEC *m, int size)
/* Check sparse matrix for symmetry - returns 0 if so and, if not, 1+ the
   column on which the first error was discovered */
{
   int i, j, jj, k, nz, nnz, *mirows, *mjrows;
   double val, valk, *mivals, *mjvals;

   
   for (i=0; i<size; i++)
      {
      if (m[i].n != size) return i+1;
      nz = m[i].non_zero;
      mirows = m[i].rows;
      mivals = m[i].vals;
      
      for (j=0; j<nz; j++)
         {
         jj = mirows[j];
         val = (int) mivals ? mivals[j] : 1.0;
         nnz = m[jj].non_zero;
         mjrows = m[jj].rows;
         mjvals = m[jj].vals;
         for (k=0; k<nnz; k++) {
	     if (mjrows[k]==i)

		{
		valk = (int) mjvals? mjvals[k] : 1.0;
		if (valk==val)
		   k=nnz+1;
		else
		   k=nnz;
		}
	  }

         if (k==nnz)

	   {


		   return i+1;
		 }
	 
         }
      }

   return 0;
}

int islands(SPARSE_VEC *m, int size)
/*
   Compute how many ``islands'' the adjacency matrix <m> implies
*/
{
   int *isl, *mirows, nz, nisl, code, i, j, jj, k, isli, islj;

   isl = Calloc(size, int);
   if (!isl) goto noroom;

   for (isl[0]=code=nisl=i=1; i<size; i++){
      isli = 0;
      nz = m[i].non_zero; 
      mirows = m[i].rows;  
      for (j=0; j<nz; j++){
         jj = mirows[j];
         if (jj<i) {
            islj = isl[jj];
            if (isli) {
               if (isli != islj){
                  /* This area bridges two islands */
                  for (k=0; k<i; k++) {
                     if (isl[k] == islj) isl[k] = isli;
                  }
                  nisl--;
               }
            } else {
               isli = islj;
            }
         }
      }
      if (isli) {
         /* This area is part of an existing island */
         isl[i] = isli;
      } else {
         /* This area starts a new island */
         code++;
         nisl++;
         isl[i] = code;
      }
   }
   free(isl);
   return nisl;

noroom:
   fprintf(stderr, "*** islands()  *** dynamic storage overflow\n");
   exit(1);
}

double sparse_el(SPARSE_VEC *x, int i)
/*
  Return i-th element of a sparse vector, i in (1..n).
*/
{
   int *xrows;
   double *xvals;
   int ix, nz;

   if (i<1 || i>x->n) return 0.0;
   
   xrows = x->rows;
   xvals = x->vals;

   if (xrows) {
      nz = x->non_zero;
      for (ix=0; ix<nz; ix++){ 
        if (xrows[ix]==(i-1)) {
          return (xvals? xvals[ix]: 1.0);
        }
      }
      return 0.0;
   } else {
     return xvals? xvals[i-1]: 1.0;
   }
}


SPARSE_VEC *sparse_tr(SPARSE_VEC *byrow, int nrow)
/*
   Transpose a sparse matrix. Input is a matrix stored as an array of
   sparse row vectors. Returned is an array of sparse column vectors.
*/
{
   SPARSE_VEC *res, *resi;
   int ncol, i, j, k, nz, nzr, not1, *ind, *iw;
   double *w, *val;

   ncol = byrow->n;
   res = Calloc(ncol, SPARSE_VEC); /* result vector */
   iw = Calloc(nrow, int);                /* work space */
   w = Calloc(nrow, double);           /* work space */
   if (!w) return 0;  /* no space */

   for (i=0, resi=res; i<ncol; i++, resi++) {
     for (nz=not1=j=0; j<nrow; j++) {
       nzr = byrow[j].non_zero;
       if (nzr) {
         ind = byrow[j].rows;
         if (ind) {
           for (k=0;(k<nzr) &&  (ind[k]!=i);k++) {}
         } else {
           k = i;
         }
         if (k<nzr) {
           iw[nz] = j;
           w[nz] = byrow[j].vals ? byrow[j].vals[k] : 1.0;
           not1 = not1 + (w[nz] != 1.0);
           nz++;
         }
       }
     }
     resi->n = nrow;
     resi->non_zero = nz;
     if (nz < nrow){
       resi->rows = Calloc(nz, int);
       if (!resi->rows) return 0;
       for (j=0; j<nz; j++) resi->rows[j] = iw[j];
     } else {
       resi->rows = (int *) 0;
     }
     if (not1) {
       resi->vals = Calloc(nz, double);
       if (!resi->vals) return 0;
       for (j=0; j<nz; j++) resi->vals[j] = w[j];
     } else {
       resi->vals = (double *) 0;
     }
   }
   free(w);
   free(iw);
   return res;
}

SPARSE_VEC *sparse_cr(int n, int nz, ...)
/*
   Create sparse vector of nominal length <n>, with <nz> non-zero elements. 
   Remaining arguments are in pairs.
   First member is element number (with 1 indicating the first) and the
   second member is the value.
*/
{
  va_list ap;
  int i;
  SPARSE_VEC *res;
  
  res = Calloc(nz, SPARSE_VEC);
  if (!res) goto no_room;
  res->n = n;
  res->non_zero = nz;
  res->rows = Calloc(nz, int);
  res->vals = Calloc(nz, double);
  if (!res->vals) goto no_room;

  va_start(ap, nz);
  for (i=0; i<nz; i++) {
    res->rows[i] = va_arg(ap, int) - 1;
    res->vals[i] = va_arg(ap, double);
  }
  va_end(ap);
  return res;

no_room:
  return (SPARSE_VEC *) 0;
}

SPARSE_VEC *sparse_c1(int n, int nz, ...)
/*
   As sparse_cr(), but all non-zero elements are 1.0 so that the second
   argument of each pair is omitted.
*/
{
  va_list ap;
  int i;
  SPARSE_VEC *res;
  
  res = Calloc(nz, SPARSE_VEC);
  if (!res) goto no_room;
  res->n = n;
  res->non_zero = nz;
  res->rows = Calloc(nz, int);
  if (!res->rows) goto no_room;
  res->vals = (double *) 0;

  va_start(ap, nz);
  for (i=0; i<nz; i++) {
    res->rows[i] = va_arg(ap, int) - 1;
  }
  va_end(ap);
  return res;

no_room:
  return (SPARSE_VEC *) 0;
}

void sparse_fr(SPARSE_VEC *v, int ncol)
/* Free dynamic storage allocated to a sparse matrix */
{
  int i;

  for (i=0; i<ncol; i++){
    if (v[i].rows) free(v[i].rows);
    if (v[i].vals) free(v[i].vals);
  }
  free(v);
}

SPARSE_VEC *sparse_kpv(SPARSE_VEC *a, SPARSE_VEC *b)
/* Kroenecker product of two vectors */
{
  int n, nz, i, j, ia, jb, ij, *r;
  double *v, ai, bj;
  SPARSE_VEC *res;

  n = a->n * b->n;
  nz = a->non_zero * b->non_zero;
  if (nz < n) {
    r = Calloc(nz, int);
    if (!r) goto noroom;
  } else {
    r = (int *) 0;
  }
  if (a->vals || b->vals){
    v = Calloc(nz, double);
    if (!v) goto noroom;
  } else {
    v = (double *) 0;
  }
  res = Calloc(1,SPARSE_VEC);
  if (!res) goto noroom;
  res->rows = r;
  res->vals = v;
  res->n = n;
  res->non_zero = nz;
  if (!r && !v) return res;
  for (ij=i=0; i<a->non_zero; i++) {
    ia = a->rows ? a->rows[i] : i;
    ai = a->vals ? a->vals[i] : 1.0;
    for (j=0; j<b->non_zero; j++, ij++) {
      jb = b->rows ? b->rows[j] : j;
      bj = b->vals ? b->vals[j] : 1.0;
      if (r) r[ij] = jb + ia * b->n;
      if (v) v[ij] = ai * bj;
    }
  }
  return res;

noroom:
   fprintf(stderr, "*** sparse_kpv *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_VEC *sparse_kp(int ma, SPARSE_VEC *a, int mb, SPARSE_VEC *b)
/* Kroenecker product of two matrices */
{
  int na, nb, mab, i, j, ij;
  SPARSE_VEC *w, *res;

  na = a->n;
  nb = b->n;
  mab = ma*mb;
  res = Calloc(mab, SPARSE_VEC);
  if (!res) goto noroom;
  for (ij=i=0; i<ma; i++) {
    if (a[i].n != na) goto dimerror;
    for (j=0; j<mb; j++, ij++) {
      if (b[j].n != nb) goto dimerror;
      w = sparse_kpv(a+i, b+j);
      sparse_cp(w, res+ij);      
      sparse_fr(w, 1);
    }
  }
  return res;

noroom:
   fprintf(stderr, "*** sparse_kp *** dynamic storage overflow\n");
   exit(1);

dimerror:
   fprintf(stderr, "*** sparse_kp *** input matrix dimension error\n");
   exit(1);
}

SPARSE_VEC *sparse_0(int n, int m)
/* Create n*m matrix of zeros */
{
  SPARSE_VEC *r;
  int i;

  r = Calloc(m, SPARSE_VEC);
  if (!r) goto noroom;
  for (i=0; i<m; i++) {
    r[i].n = n;
    r[i].non_zero = 0;
    r[i].rows = (int *) 0;
    r[i].vals = (double *) 0;
  }
  return r;

noroom:
   fprintf(stderr, "*** sparse_0 *** dynamic storage overflow\n");
   exit(1);
}

SPARSE_VEC *sparse_ccv(int m, SPARSE_VEC *a, SPARSE_VEC *b)
/* Concatenate two matrices vertically (a above b) */
{
  int na, nb, i, j, k, jj, n, nz, *rows;
  double vv, *vals;
  SPARSE_VEC *r;

  na = a->n;
  nb = b->n;
  n = na+nb;
  r = Calloc(m, SPARSE_VEC);
  if (!r) goto noroom;
  for (i=0; i<m; i++) {
    if (a[i].n != na || b[i].n != nb) goto dimerror;
    r[i].n = n;
    nz = a[i].non_zero + b[i].non_zero;
    r[i].non_zero = nz;
    if (a[i].rows || b[i].rows) {
      rows = Calloc(nz, int);
      if (!rows) goto noroom;
    } else {
      rows = (int *) 0;
    }
    if (a[i].vals || b[i].vals) {
      vals = Calloc(nz, double);
      if (vals) goto noroom;
    } else {
      vals = (double *) 0;
    }
    r[i].rows = rows; 
    r[i].vals = vals;
    if (rows || vals) {
      for (k=j=0; j < a[i].non_zero; k++, j++) {
        jj = a[i].rows ? a[i].rows[j] : j;
        vv = a[i].vals ? a[i].vals[j] : 1.0;
        if (rows) rows[k] = jj;
        if (vals) vals[k] = vv;
      }
      for (j=0; j < b[i].non_zero; k++, j++) {
        jj = b[i].rows ? b[i].rows[j] : j;
        vv = b[i].vals ? b[i].vals[j] : 1.0;
        if (rows) rows[k] = jj + na;
        if (vals) vals[k] = vv;
      }
    }
  }
  return r;

noroom:
   fprintf(stderr, "*** sparse_ccv *** dynamic storage overflow\n");
   exit(1);

dimerror:
   fprintf(stderr, "*** sparse_ccv *** input matrix dimension error\n");
   exit(1);
}
  
SPARSE_VEC *sparse_cch(int ma, SPARSE_VEC *a, int mb, SPARSE_VEC *b)
/* Concatenate two matrices horizontally (a, b) */
{
  int n, m, i, j;
  SPARSE_VEC *r;

  n = a->n;
  m = ma+mb;
  r = Calloc(m, SPARSE_VEC);
  if (!r) goto noroom;
  for (i=j=0; j<ma; i++, j++) {
    if (a[j].n != n) goto dimerror;
    sparse_cp(a+j, r+i);
  }
  for (j=0; j<ma; i++, j++) {
    if (b[j].n != n) goto dimerror;
    sparse_cp(b+j, r+i);
  }
  return r;

noroom:
   fprintf(stderr, "*** sparse_cch *** dynamic storage overflow\n");
   exit(1);

dimerror:
   fprintf(stderr, "*** sparse_cch *** input matrix dimension error\n");
   exit(1);
}

void sparse_cp(SPARSE_VEC *fr, SPARSE_VEC *to)
/* Copy vector */
{
  int i, nz, *rows;
  double *vals;

  to->n = fr->n;
  nz = fr->non_zero;
  to->non_zero = nz;
  if (fr->rows){
    rows = Calloc(nz, int);
    if (!rows) goto noroom;
    for (i=0; i<nz; i++) rows[i] = fr->rows[i];
  } else {
    rows = (int *) 0;
  }
  if (fr->vals){
    vals = Calloc(nz, double);
    if (!vals) goto noroom;
    for (i=0; i<nz; i++) vals[i] = fr->vals[i];
  } else {
    vals = (double *) 0;
  }
  to->rows = rows;
  to->vals = vals;
  return;

noroom:
   fprintf(stderr, "*** sparse_cp *** dynamic storage overflow\n");
   exit(1);
}

double sparse_ip(SPARSE_VEC *x, SPARSE_VEC *y, double *w)
/*
   Inner product of two sparse vectors: x-transpose.w.y
*/
{
   int n, nzx, nzy, i, j, ii, jj;
   double xw, yw, wt, ip;

   if (x->n != y->n) {
      fprintf(stderr,"*** sparse_ip *** unequal length vectors\n");
      exit(1);
   }
   n = x->n;
   nzx = x->non_zero;
   nzy = y->non_zero;
   for (i=j=0, ip=0.0; i<nzx && j<nzy;) {
      ii = x->rows? x->rows[i] : i;
      jj = y->rows? y->rows[j] : j;
      if (ii == jj) {
         xw = x->vals? x->vals[i] : 1.0;
         yw = y->vals? y->vals[j] : 1.0;
         wt = w? w[ii] : 1.0;
         ip += (xw*yw*wt);
         i++;
         j++;
      }
      else {
         if (ii < jj)
            i++;
         else
            j++;
      }
   }
   return ip;
}

double sparse_ss(SPARSE_VEC *x, double *w)
/*
   (Weighted) sum of squares of the elements of x
*/
{
   int nz, i, ii;
   double ss, xw, wt;

   nz = x->non_zero;
   for (i=0, ss=0.0; i<nz; i++){
      ii = x->rows? x->rows[i] : i;
      xw = x->vals? x->vals[i] : 1.0;
      wt = w? w[ii] : 1.0;
      ss += (xw*xw*wt);
   }
   return ss;
}

int sparse_so(SPARSE_VEC *x, SPARSE_VEC *y)
/*
   Are two vectors strictly orthogonal (inner product zero for any weights)
*/
{
   int n, nzx, nzy, i, j, *rx, *ry;
   
   if (x->n != y->n) {
      fprintf(stderr,"*** sparse_so *** unequal length vectors\n");
      exit(1);
   }
   n = x->n;
   nzx = x->non_zero;
   nzy = y->non_zero;
   if (nzx == 0 || nzy == 0) return 1; /* Must be strictly orthogonal */
   if ((nzx+nzy) > n) return 0; /* Can't be strictly orthogonal */
   rx = x->rows;
   ry = y->rows;
   for (i=j=0; i<nzx && j<nzy;) {
      if (rx[i] == ry[j])
         return 0;
      else
         if (rx[i] < ry[j])
            i++;
         else
            j++;
   }
   return 1;
}

SPARSE_VEC *sparse_gm_cr(int m, SPARSE_VEC *x)
/*
   Create an empty Grammian matrix to hold x-transpose.w.x
   m is the number of columns in x
   Upper triangle only is held
*/
{
   int i, j, k, n, *wk;
   SPARSE_VEC *res;

   res = Calloc(m, SPARSE_VEC);
   wk = Calloc(m, int);
   n = x->n;
   for (i=0; i<m; i++) {
      if (x[i].n != n) {
         fprintf(stderr, "*** sparse_gm_cr *** Unequal length columns\n");
         exit(1);
      }
      for (j=k=0; j<i; j++) {
         if (sparse_so(x+i, x+j)) wk[k++] = j;
      }
      wk[k++] = i;
      res[i].n = m;
      res[i].non_zero = k;
      res[i].vals = Calloc(k, double);
      if (k<m) {
         res[i].rows = Calloc(k, int);
         for (j=0; j<k; j++) res[i].rows[j] = wk[j];
      }
      else {
         res[i].rows = (int *) 0;
      }
   }
   free(wk);
   return res;
}

void sparse_gm_cal(int m, SPARSE_VEC *x, double *w, SPARSE_VEC *gm)
/*
   Calculate the Grammian x-transpose.w.x and store in the empty
   matrix gm. Note that the structure of gm defines which inner
   products need to be calculated --- it must be set up using
   sparse_gm_cr().
*/
{
   int i, j, jj, nz;
   
   if (gm->n != m) {
      fprintf(stderr, "*** sparse_gm_cal *** Wrong size matrix\n");
      exit(1);
   }
   for (i=0; i<m; i++) {
      nz = gm[i].non_zero;
      for (j=0; j<nz; j++) {
         jj = gm[i].rows? gm[i].rows[j] : j;
         gm[i].vals[j] = (i==jj)? sparse_ss(x+i, w) : sparse_ip(x+i, x+jj, w);
      }
   }
   return;
}






void elements(SPARSE_VEC *x, 
	      double *res)
{
  int i=0;
  int pointer=0;
  int t;
  

 
  if(!(x->rows)){
 
    
    for(i=0;i<x->n;i++){
      res[i] = x->non_zero ? (x->vals ? (x->vals)[i] : 1.0 ) : 0.0;
      
    }
    
  }
  

  if(x->rows){
      while((i<x->n) && (pointer<x->non_zero)){
         
      if ((x->rows)[pointer] == i) {
	t= (x->vals ? (x->vals)[pointer] : 1.0) ;
   	
        res[i]=t;

	
	pointer++;
    }
      i++;
      
  }
  
}

  return;
  
  
  
}

  

SPARSE_VEC *select_elements(SPARSE_VEC *vector,
		   int first,
		   int last)
{
  SPARSE_VEC *res;
  double *data;
 

 
  data = Calloc(vector->n,double);

  elements(vector,data);
  
   
  res  =sparse_vec((last - first + 1), data+first-1);
  
  

  
  
  Free(data);

  
  
  
  return res;


    
  
  
}

double last_value(SPARSE_VEC *vector)
{
  int n,nz;
  n=vector->n;
  nz=vector->non_zero;
  
  if (vector->vals){
    if(vector->rows[(nz-1)] == (n-1)){
      return(vector->vals[(nz-1)]);
      
	}else{
	  return (0.0);
	}

  }else{
    if(vector->rows[(nz-1)] == (n-1)){
      return   1.0;
      
	}else{
	  return 0.0;
	}
  }
  
  
}


