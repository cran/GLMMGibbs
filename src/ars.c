#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myrand.h"
#include "ars.h"
#include "S.h"
/* Constants */

#define EXP_TAIL_TOL 0.1  /* Tolerance for length of exponential tail */
#define GRAD_TOL   0.0001 /* Gradient tolerance for switching exp -> linear */
#define PEAK_TOL   3.0    /* Tolerance for peak of envelope above inner hull */
#define BIG_LD     80.0   /* Largest value for a log density (after offset)*/
#define SMALL_LD  -80.0   /* Smallest value */
#define HEAD_ROOM  20.0   /* Headroom when computing offset */
#define long int
/* Global variables */

static char dfname[80] = "ars.dmp"; 

/* Local function declarations */

void ran_pwe(ARSE *);
void join_seg(ARSE *, EXP_SEG *, EXP_SEG *);
void area_arse(ARSE *, EXP_SEG *);
double area_seg(EXP_SEG *);
void dump_arse(char*, ARSE *);

long rejection_sample(ARSE *e){
/*
   If <e> is a valid envelope, this routine samples it and puts the
   result in e->y. It also tests whether this value can be accepted
   without further ado, by testing the lower density bound.

   If <e> is not a valid envelope, it does nothing.

   Returns 0 if value e->y is to be accepted. Returns 1 if the density and



   gradient at e->y must be computed.
*/
   double dy, ihull;
   EXP_SEG *s;

   if (e->status == -1) {            /*  Not a valid envelope */
      return 1;
   }

   if (e->status == 1) {             /* Accept old value ... */
      e->status = 0;                 /* but not more than once! */
      return 0;
   }

   /* Otherwise, sample a new value from the envelope */

   if (e->area ==  0.0) 
      dump_arse("*** ars (rejection_sample): Envelope error", e);

   ran_pwe(e);

   /* Compute acceptance bound:
      log (envelope*uniform_ran) = log_envelope - unit_exp   */

   s = e->y_seg;
   dy = e->y  -  s->y;
   e->lb = (s->lf + dy*s->g) - unitexp();

   /* Compute inner hull and accept if above boundary */

   if (dy < 0.0) {
      if (s->left) {
	 ihull = s->lf + dy*s->g0;
	 if (ihull > e->lb) return 0;
      }
   } else {
      if (s->right) {
	 ihull = s->lf + dy*s->g1;
	 if (ihull > e->lb) return 0;
      }
   }
   return 1;                          /* Return for computation of lf, g */
}

void update_arse(ARSE *e, double lf, double g) {
/*
   If last point is not accepted, updates the sampling envelope
*/

   EXP_SEG *sl, *sr, *sn;
   
   /* Test to see if the current y is acceptable.
      If so, flag for acceptance and return */

   if (e->status == 0 && lf > e->lb){
      e->status = 1;
      return;
   }

   /* Otherwise, use this point to update the envelope */

   /* Test for envelope full */

   if (e->n_seg == e->max_seg) {
      if (e->area > 0.0 ) {
	 return;
      } else
	 dump_arse("*** ars (update_arse): Envelope overflow", e);
   }

   /* Otherwise, add a new segment */

   if (!e->y_seg)
      dump_arse("*** ars (update_arse) : Bug - no y_seg", e);
   sn = e->env + (e->n_seg++);
   sn->y = e->y;
   sn->lf = lf;
   sn->g = g;
   if (e->y < e->y_seg->y) {
      sr = e->y_seg;
      sl = sr->left;
      if (!sl) e->first = sn;
      join_seg(e, sn, sr);
      join_seg(e, sl, sn);
   } else if (e->y > e->y_seg->y) {
      sl = e->y_seg;
      sr = sl->right;
      join_seg(e, sl, sn);
      join_seg(e, sn, sr);
   } else {
      dump_arse("*** ars (update_arse) : Repeated evaluation", e);
   }

   area_arse(e, sn);

   return;
}

ARSE *alloc_arse(long max_seg)
{
   ARSE *e;
   EXP_SEG *s;

   e = Calloc(1,ARSE);
   s = Calloc(max_seg, EXP_SEG);
   if (!e || !s) 
      dump_arse("*** ars (alloc_arse): No dynamic memory", e);

   e->env = s;
   e->first = e->env;
   e->max_seg = max_seg;
   e->n_seg = 0;
   return e;
}

void init_arse(ARSE *e, double y, double lf, double g, double c)
{
   EXP_SEG *s;
   double step, sr;

   if (c < 0.0) 
      dump_arse("*** ars (init_arse): Wrong sign for curvature", e);

   s = e->y_seg = e->first = e->env;
   e->area = 0.0;
   e->n_seg = 1;
   e->status = -1;
   e->min_lf = e->max_lf = lf;
   e->offset = lf - BIG_LD + HEAD_ROOM;

   /* First segment */

   s->lf = lf;
   s->y = y;
   s->g = g;
   s->left = s->right = (EXP_SEG *) 0;
   s->f0 = s->f1 = 0.0;
   
   step = g/c;
   sr = sqrt(c);
   
   /* Suggested evaluation point to the right */

   if (g > -EXP_TAIL_TOL*sr) {
      s->f1 = -1.0;
      if (step > 0.0) y += step;
      s->y1 = e->y = y + 0.5/sr;
   }

   /* Suggested evaluation point to the left */
   
   if (g < EXP_TAIL_TOL*sr) {
      s->f0 = -1.0;
      if (step < 0.0) y += step;
      s->y0 = e->y = y - 0.5/sr;
   }
   return;
}

void free_arse(ARSE *e)
{
   free(e->env);
   free(e);
}

void arse_dump(char *filename)
{
   char *c;
   for (c=dfname; *(c++) = *(filename++);) {}
}

   
/* Local functions */

void ran_pwe(ARSE *e)
{
   double u, cu, fs;
   EXP_SEG *s;

   /* First sample the segments ... */

   u = u_random()*e->area;
   for (s=e->first, cu=0.0; (cu += s->area) < u; s = s->right){
     if (!s->right) 
        dump_arse("*** ars (ran_pwe): Segment out of range", e);
   }
   e->y_seg = s;

   /* ... then the value within the segment. */

   if (!s->right) {
      e->y = s->y0 - unitexp()/s->g ; /* Sample exponential right tail */
      if (e->y < s->y0) {
         dump_arse("*** arse(ran_pwe<R>) : Result outside segment", e);
      }
   } else {
      if (!s->left) {
         e->y = s->y1 - unitexp()/s->g; /* Sample exponential left tail */
         if (e->y > s->y1) {
            dump_arse("*** arse(ran_pwe<L>) : Result outside segment", e);
	 }
      } else {
         u = s->lf1 - s->lf0;
         if (u < GRAD_TOL && u > -GRAD_TOL) {  /* Sample uniform distrn */
            e->y = s->y0 + u_random()*(s->y1 - s->y0);
            if (e->y < s->y0 || e->y > s->y1) {
               dump_arse("*** arse(ran_pwe<U>) : Result outside segment", e);
	    }
         } else {   /* Sample truncated exponential distrn */
	    u = u_random();
            fs = u*s->f0 + (1.0-u)*s->f1;
            e->y = s->y + (log(fs) - (s->lf - e->offset))/s->g;
            if (e->y < s->y0 || e->y > s->y1) {
               dump_arse("*** arse(ran_pwe<I>) : Result outside segment", e);
	    }
	 }
      }
   }
   return;
}

void join_seg(ARSE *e, EXP_SEG *sl, EXP_SEG *sr) {
/*
   Compute boundary of 2 segments
*/

   double y, f, g, lf, c, s, min_width, cut, dy, dll, dlr, inner;

   if (sl) {
      sl->right = sr;
      sl->f1 = 0.0;
   }
   if (sr) {
      sr->left = sl;
      sr->f0 = 0.0;
   }
   
   /* If no leftward segment, <sr> must be the leftmost */
   
   if (!sl) {
      if (!sr->right) {
	 dump_arse("*** ars (join_seg): Open left seg", e);
      } else {
	 c = (sr->g - sr->right->g)/(sr->right->y - sr->y);
	 if (c<0.0)
	    dump_arse("*** ars (join_seg): Positive 2nd derivative <2>",e);
         s = sqrt(c);
         if (sr->g <= EXP_TAIL_TOL*s){  /* Suggest a new evaluation point */
            sr->f0 = -1.0;
 	    sr->y0 = sr->g < 0.0 ?
               sr->y + sr->g/c - 0.5/s :
               sr->y - 0.5/s;
	 }
      }
      return;
   }

   /* If no rightward segment, <sl> must be the rightmost */

   if (!sr) {
      if (!sl->left) {
	 dump_arse("*** ars (join_seg): Open right seg", e);
      } else {
	 c = (sl->g - sl->left->g)/(sl->left->y - sl->y);
	 if (c<0.0)
	    dump_arse("*** ars (join_seg): Positive 2nd derivative <3>",e);
         s = sqrt(c);
         if (sl->g >= - EXP_TAIL_TOL*s) { /* Suggest a new evaluation point */
            sl->f1 = -1.0;
	    sl->y1 = sl->g > 0.0 ?
               sl->y + sl->g/c + 0.5/s:
               sl->y + 0.5/s;
	 }
      }
      return;
   }

   /* Otherwise, join the two segments */
   
   if (sl->y >= sr->y)
      dump_arse("*** ars (join_seg): Segments out of order", e);
   if (sl->g < sr->g)
      dump_arse("*** ars (join_seg): Positive 2nd derivative", e);
   dy = sr->y - sl->y;
   dll = sr->lf - sl->lf - sr->g*dy;
   dlr = sl->lf - sr->lf + sl->g*dy;
   if (dll < 0.0 || dlr < 0.0){
      printf("dll=%lf,  dlr=%lf\n",dll,dlr) ;
      dump_arse("*** ars (join_seg): Errors in funct/grad calculation", e);
    }
   cut = dll/(dll+dlr);
   sl->lf1 = sr->lf0 = lf = sl->lf + sl->g*cut*dy;
   sl->y1 = sr->y0 = y = sl->y + cut*dy;
   sl->g1 = sr->g0 = g = (sr->lf - sl->lf)/dy;

   /* If "steeple", suggest reevaluation at midpoint */
   
   inner = sl->lf + g*cut*dy;
   if ((lf-inner) > PEAK_TOL && sl->g > 0.0 && sr->g < 0.0) {
      sl->y1 = sr->y0 = (sl->y + sr->y)/2.0;
      sl->f1 = sr->f0 = -1.0;
   } else {
      if (lf > e->max_lf) e->max_lf = lf;
      if (lf < e->min_lf) e->min_lf = lf;
   }
   return;
}

double area_seg(EXP_SEG *s) {
/*
   Compute area of segment and return its value
*/

   double test;

   if (!s) return 0.0;

   /* If gradient is very small, use linear approximation */

   if (s->f0 < 0.0 || s->f1 < 0.0)
      return (s->area = -1.0);

   /* Special algorithm for shallow gradient */

   if (s->left && s->right) {
      test = s->lf1 - s->lf0;
      if (test < GRAD_TOL && test > - GRAD_TOL)
	 return (s->area = (s->y1 - s->y0)*(s->f1 + s->f0)/2.0);
   }

   /* The usual method */

   return (s->area = (s->f1 - s->f0)/s->g);
}

void area_arse(ARSE *e, EXP_SEG *s) {
/*
   Recompute densities and areas after inserting new segment <s>.

   If not a valid envelope, sets e->status=-1
   and value and segment for next function evaluation.
   
*/

   EXP_SEG *sl, *sr;
   double off;

   sl = s->left;
   sr = s->right;
   off = e->offset;

   /* If no overflow, ... */
   
   if ((e->max_lf - off) < BIG_LD) { /* Calculate areas of new segments only */
      if (sl && s->f0 != -1.0) {
         sl->f1 = s->f0 = exp(s->lf0 - off);
      }
      if (sr && s->f1 != -1.0) {
         sr->f0 = s->f1 = exp(s->lf1 - off);
      }
      area_seg(sl);
      area_seg(s);
      area_seg(sr);

   } else {  /* ... otherwise reset offset and recalculate all areas */ 

      e->offset = off = e->max_lf - BIG_LD + HEAD_ROOM;
      for (s = e->first; s; s = s->right) {
         if (s->left && s->f0 != -1.0) {
	    s->f0 = exp(s->lf0 - off);
         }
         if (s->right && s->f1 != -1.0) { 
	    s->f1 = exp(s->lf1 - off);
         }
         area_seg(s);
      }
   }
   /*
      Calculate total area under envelope.
      If not a valid envelope, set e->y to the next evaluation point
      and e->status to -1
      Otherwise set e->status to zero
   */

   for (e->area=0.0, s=e->first; s; s = s->right) {
      if (s->area < 0.0) {
	 e->status = -1;
	 e->area = 0.0;
	 e->y_seg = s;
	 e->y = (s->f0 < 0.0) ? s->y0 : s->y1;
	 return;
      } else {
	 e->area += s->area;
      }
   }
   e->status = 0;
   return;
}


/* Diagnostic dump on errors */


void dump_arse(char *message, ARSE *e)
{
   FILE *dumpfile;  
   EXP_SEG *s;
   long i;
 
   fprintf(stderr, "\n *** ARS error ... Terminating\n");
   dumpfile = fopen(dfname, "w");
   if (!dumpfile) {
      fprintf(stderr, "*** Unable to open dumpfile <%s>- revert to stderr\n",
         dfname);
      dumpfile = stderr;
   }
   fprintf(dumpfile, "%s\n\n", message);
   fprintf(dumpfile, "State of sampling envelope when error occurred:\n\n");
   fprintf(dumpfile,
   "========================================================================\n");
   fprintf(dumpfile, "n_seg = %2d, Offset = %8.1le,          Area = %10.2le\n",
      e->n_seg, e->offset, e->area);
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   fprintf(dumpfile,
   "       y        f0       lf       f1        g          Area      e->y\n");
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   for (s = e->first; s; s = s->right){
      if (s->left)
	 fprintf(dumpfile, " == %8.1le)\n", s->y0);
      fprintf(dumpfile, "%2d=%8.1le %8.1le %9.1le %8.1le %8.1le  %10.2le",
	 1 + (int)(s - e->env), s->y, s->f0, s->lf, s->f1, s->g, s->area);
      if (s == e->y_seg){
	 fprintf(dumpfile, "%10.1le", e->y);
	 if (e->status < 0) fprintf(dumpfile, " *");
      
      }
      if (s->right)
	 fprintf(dumpfile, "\n(%8.1le", s->y1);
      else
         fprintf(dumpfile, "\n");
   }
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   fprintf(dumpfile, "All segments (in order of evaluation)\n");
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   for (i=0, s= e->env; i<e->n_seg; i++, s++) {
      fprintf(dumpfile, "%2d=%8.1le %8.1le %9.1le %8.1le %8.1le  %10.2le\n",
	 1 + i, s->y, s->f0, s->lf, s->f1, s->g, s->area);
   }
   fprintf(dumpfile,
   "========================================================================\n");
#ifdef FALSE
   fflush(dumpfile);
#endif 
   error("ARS file error---if you haven't clearly over-parameterised the model, please email
          mylesj@icrf.icnet.uk\n");
}




void show_arse( ARSE *e)
{
   FILE *dumpfile;  
   EXP_SEG *s;
   long i;
 
   dumpfile = stdout;

   fprintf(dumpfile,
   "========================================================================\n");
   fprintf(dumpfile, "n_seg = %2d, Offset = %8.1le,          Area = %10.2le\n",
      e->n_seg, e->offset, e->area);
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   fprintf(dumpfile,
   "       y        f0       lf       f1        g          Area      e->y\n");
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   for (s = e->first; s; s = s->right){
      if (s->left)
	 fprintf(dumpfile, " == %8.1le)\n", s->y0);
      fprintf(dumpfile, "%2d=%8.1le %8.1le %9.1le %8.1le %8.1le  %10.2le",
	 1 + (int)(s - e->env), s->y, s->f0, s->lf, s->f1, s->g, s->area);
      if (s == e->y_seg){
	 fprintf(dumpfile, "%10.1le", e->y);
	 if (e->status < 0) fprintf(dumpfile, " *");
      
      }
      if (s->right)
	 fprintf(dumpfile, "\n(%8.1le", s->y1);
      else
         fprintf(dumpfile, "\n");
   }
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   fprintf(dumpfile, "All segments (in order of evaluation)\n");
   fprintf(dumpfile,
   "------------------------------------------------------------------------\n");
   for (i=0, s= e->env; i<e->n_seg; i++, s++) {
      fprintf(dumpfile, "%2d=%8.1le %8.1le %9.1le %8.1le %8.1le  %10.2le\n",
	 1 + i, s->y, s->f0, s->lf, s->f1, s->g, s->area);
   }
   fprintf(dumpfile,
   "========================================================================\n");
}







