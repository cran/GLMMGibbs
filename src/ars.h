/* ars.h */

/* Header file for adaptive rejection sampling routines */

typedef struct seg {         /* Segment of piecewise exponential envelope */
   struct seg *left, *right; /* Pointers to adjacent segments */
   double y,  lf, g,         /* y, log_density(y), gradient */
	  y0, f0, lf0, g0,   /* Left boundary, g0 is gradient of lower bound */
	  y1, f1, lf1, g1,   /* Right boundary */
	  area;              /* Area of segment */
} EXP_SEG;

typedef struct {             /* Adaptive Rejection Sampling Envelope */
   int max_seg, n_seg;       /* Number of pieces in envelope */
   EXP_SEG *env;             /* Pointer to array of segments */
   EXP_SEG *first;           /* Pointer to leftmost segment */
   double area;              /* Total area of envelope */
   double offset;            /* Current offset for log density */
   double min_lf, max_lf;    /* Min and max of log density in envelope */
   double y;                 /* Current y, and ... */
   EXP_SEG *y_seg;           /* ... segment where it lies, and */
   int status;               /* ... a status flag  :
				    -1   envelope incomplete
				     0   envelope ready for sampling
				     1   current value accepted */
   double lb;                /* Lower bound of log density for acceptance */
} ARSE;

/*
   These routines are used to do adaptive rejection sampling, as follows ...

   envelope = alloc_arse(Max no segments to be used);

   for (
      init_arse(envelope, var, log_dens, grad, curv);
	 rejection_sample(envelope);
	    update_arse(envelope, log_dens, grad)
      ) {
      var = envelope->y;
      log_dens = ?;
      grad = ?;
   }

   On completion, envelope->y will contain an appropriate value. Note,
   however, that it might differ from that used during the last evaluation
   of the body of the loop.

   free_arse(envelope) frees up the dynamic memory.
   arse_dump(filename) opens a dump file for trapping errors.
*/

int rejection_sample(ARSE*);
ARSE *alloc_arse(int);
void init_arse(ARSE*, double, double, double, double);
void update_arse(ARSE*, double, double);
void free_arse(ARSE*);
void arse_dump(char*);
void show_arse(ARSE*);






