/* random.c */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "myrand.h"
#include <limits.h>
#include "misc.h"

#define TWOPI 6.2831853071795865
#define E 2.71828182



static int snd_call=0;


#define znew  ((z=k[ii]*(z&65535)+(z>>16))<<16)
#define wnew  ((ww=k[jj]*(ww&65535)+(ww>>16))&65535)
#define IUNI  (znew+wnew)
#define UNI   (znew+wnew)*2.328306e-10
static unsigned long z=362436069, ww=521288629;
static int ii=0,jj=1;



static long k[] = {18000,
18030,
18273,
18513,
18879,
19074,
19098,
19164,
19215,
19584,
19599,
19950,
20088,
20508,
20544,
20664,
20814,
20970,
21153,
21243,
21423,
21723,
21954,
22125,
22188,
22293,
22860,
22938,
22965,
22974,
23109,
23124,
23163,
23208,
23508,
23520,
23553,
23658,
23865,
24114,
24219,
24660,
24699,
24864,
24948,
25023,
25308,
25443,
26004,
26088,
26154,
26550,
26679,
26838,
27183,
27258,
27753,
27795,
27810,
27834,
27960,
28320,
28380,
28689,
28710,
28794,
28854,
28959,
28980,
29013,
29379,
29889,
30135,
30345,
30459,
30714,
30903,
30963,
31059,
31083,
31215,
31353,
31488,
31743,
32430,
32718,
33105,
33189,
33249,
33375,
33378,
33663,
33768,
33858,
33894,
34158,
34323,
34383,
34590,
34653,
34890,
35355,
35523,
35643,
36309,
36594,
36804,
36969,
37698,
37935,
37959,
38079,
38223,
38283,
38484,
38568,
38610,
38649,
38733,
38850,
39444,
39618,
39690,
39948,
40833,
40995,
41019,
41064,
41289,
41628,
41793,
41874,
42153,
42444,
42513,
42594,
42633,
42699,
42819,
42903,
42975,
43038,
43155,
43473,
43563,
43995,
44019,
44568,
44574,
44994,
45723,
45729,
45780,
45789,
45915,
45939,
46515,
47088,
47529,
48015,
48033,
48195,
48204,
48393,
49209,
49248,
49299,
49458,
50034,
50223,
50580,
50589,
50694,
50853,
50988,
51198,
51558,
51618,
51729,
51744,
51813,
51873,
51933,
52023,
52215,
52275,
52509,
52743,
52950,
53130,
53199,
53529,
53709,
53898,
53934,
53958,
54144,
54168,
54399,
54474,
54564,
54885,
55044,
55074,
55179,
55254,
55680,
55809,
55848,
55869,
56205,
56538,
56604,
56790,
56859,
57039,
57204,
57225,
57525,
57603,
57774,
57780,
57918,
58149,
58368,
58443,
58758,
59253,
59325,
59775,
60009,
60060,
60489,
60735,
60990,
61140,
61578,
61914,
62505,
62634,
62778,
62790,
62865,
62874,
62904,
63129,
63273,
63444,
63663,
63765,
63885,
64185,
64314,
64455,
64545,
64860,
65184,
65595,
65904,
65943,
66429,
66705,
66909,
67014,
67053,
67158,
67434,
67608,
67734,
67884,
67893,
67959,
68169,
68178,
68523,
68838,
69033,
69180,
69414,
69435,
69489,
69564,
69669,
69894,
69918,
70149,
70173,
70269,
70353,
71148,
71370,
72084,
72309,
72564,
72648,
72840,
72954,
73125,
73185,
73398,
73608,
73635,
73695,
73698,
73803,
73863,
73890,
74034,
74433,
74559,
74724,
74949,
75165,
75348,
75423,
75570,
75819,
75963,
75999,
76188,
76200,
76260,
76344,
76398,
76914,
76929,
77304,
77649,
78273,
78480,
78795,
78855,
78879,
79008,
79014,
79074,
79173,
79224,
79245,
79275,
79488,
79809,
79908,
80103,
80160,
80643,
80790,
80793,
80799,
81000,
81018,
81270,
81408,
81510,
81858,
82014,
82203,
83103,
83289,
83340,
83424,
83565,
83688,
83793,
83928,
84240,
84294,
84495,
84885,
85008,
85074,
85659,
86064,
86070,
86145,
86229,
86235,
86295,
86478,
86964,
87318,
87435,
87453,
87555,
88374,
88590,
88815,
89025,
89148,
89238,
89340,
89364,
89379,
89448,
89889,
90300,
90453,
90555,
90624,
90774,
91089,
91458,
91653,
91803,
92043,
92238,
92445,
92664,
92694,
92700,
92820,
93183,
93303,
93438,
93810,
93855,
93903,
94395,
94509,
94614,
94950,
95010,
95025,
95073,
95280,
95664,
95805,
95889,
95934,
95979,
96045,
96060,
96063,
96165,
96168,
96210,
96624,
96873,
98013,
98184,
98280,
98349,
98853,
99105,
		   99228};


double asran()
{
  double res;
  snd_call=0;
  
 
  res=UNI;
  return res;
  
}

void asran_seed(int i, int j, int k)
{

  ii = (i+j)    % 445;
  jj = (j+k)    % 445;
  if(jj=ii){
    jj++;
  }
  
}


  
  



double u_random()
{
double u;
u=asran();
return u;
}

/* N(0,1) deviates */

void norm2(double *g1, double *g2)
/* Pair of N(0,1) deviates (Alg 3.6 of Ripley) */
   {
   double u1, u2, w, c;

   do
      {
      u1 = 2*asran() - 1;
      u2 = 2*asran() - 1;
      w = u1*u1 + u2*u2;
      }
      while (w>=1.0);
   c = sqrt(-2.0*log(w)/w);
   *g1 = u1*c;
   *g2 = u2*c;
   return;
   }


static double snd_save;

void snd_set(int dummy)
{
  snd_call=0;
  
}


double snd()
/* Generate standard normal deviate */
   {
   double this_one;

   snd_call = !snd_call;
   if (snd_call)
      {
      norm2(&this_one, &snd_save);
      return(this_one);
      }
   else
      return(snd_save);
   }

/* Unit exponential variates by the von Neumann method (Alg 3.7 of Ripley) */

double unitexp()
   {
   double i, t, u, ustar;

   for ( i=0.0; ; i++)
      {
      t = u = u_random();
      
      do
	 {
	 ustar = u_random();
	 if (u <= ustar)
	    return i+t;
	 else
	    u = u_random();
	 }
	 while (u<ustar);
      }
   }

/* Chi-squared variates */

double chi2(int df)
   {
   int i;
   double x2, s;

   if (df<=100) /* Use sum of unit exponential and squared normal variates */
      {
      for (x2=0.0, i=df; i>1; i -= 2)
	 x2 += 2.0*unitexp();
      if (i)
	 x2 += (s=snd())*s;
      return(x2);
      }
   else /* Use gamma generator */
      return 2.0*gg((double)df/2.0);
   }




/* Generate gamma variates.
   Taken from algorithms 3.19 and 3.20 from Ripley's book 

   Use macro gg() to select appropriate function
*/

double gs(double alpha)
/* For alpha < 1.0 */
   {
   double b, p, x;

   for(b=(alpha+E)/E;;)
      {
      p = b * u_random();
      if (p<=1.0)
	 {
	 x = pow(p, 1./alpha);
	 if (x<=-log(u_random()))
	    return(x);
	 }
      else
	 {
	 x = -log((b-p)/alpha);
	 if (pow(x,alpha-1.0)>=u_random())
	    return(x);
	 }
      }
   }

static double aprev=0.0, c1, c2, c3, c4, c5;

double gcf(double alpha)
/* For alpha>1.0 */
   {
   double aa, u1, u2, w, x;

   if (alpha!=aprev)
      {
      aprev = alpha;
      c1 = alpha-1.0;
      aa = 1.0/c1;
      c2 = aa*(alpha-1.0/(6.0*alpha));
      c3 = 2.0*aa;
      c4 = c3 + 2.0;
      if (alpha>2.5)
         c5 = 1.0/sqrt(alpha);
      }
   for(;;)
      {
      do
	 {
	 u1 = u_random();
	 u2 = u_random();
	 if (alpha>2.5)
	    u1 = u2 + c5*(1.0-1.86*u1);
	 }
	 while ( (u1>=1.0) || (u1<=0.0) );
      w = c2*u2/u1;
      if ( (c3*u1 + w + 1.0/w)<=c4 || (c3*log(u1) - log(w) + w)<1.0 )
	 return(c1*w);
      }
   }

double gcra(double alpha)
/* Gamma variate by cube root approximation as given by Abramowitz & Stegun
   section 26.4.14. Perhaps necessary for very large gamma? */
   {
   double a9, w;
   a9 = alpha*9.0;
   w = (a9 - 1.0 + sqrt(a9)*snd())/a9;
   return alpha*w*w*w;
   }


double rbeta(double a,double b)
{
  double res,x,y;
  x=gg(a);
  y=gg(b);
  res = x/(x+y);
  return res;
}



/* Sample from a discrete distribution */

int sample_dd(int m, const double cdf[])
/* Sample from the discrete distribution with m states, labelled 0,...,m-1.
   cdf[] is an array, length m, with elements proportional to the
   cumulative probabilities. Thus, cdf[i]/cdf[m-1] i=0,...,m-1 are the
   cumulative probabilities.
   Uses the method of binary search without frills. Intended for situations
   where only a few numbers are required for each distribution, so that
   the set up time involved in more efficient methods is not worth it!
*/
   {
   double p;
   int l, h, i;

   p = cdf[m-1]*u_random();
   for(l=0, h=m-1; l<h; )
      {
      i = (l+h)/2;
      if (p>cdf[i])
	 l = i+1;
      else
	 h = i;
      }
   return(l);
   }

/* Poisson random variates (Ripley, Alg 3.3 ,page 55 )*/

int r_poisson(double mu){
   
  int n;
  double c, p;

  c=exp(-mu); 
  p=1.0;
  n=0;
  do {
    p *= u_random();
    n++;
  } while (p >= c);
  return n-1;
}

/* Binomial random variates (Ripley, Alg 3.14, page 78) */

int r_binom(int n, double p){

  const int kmax = 50;
  
  int k, i, x;
  double theta, v, g1, g2;

  for (k=n, theta=p, x=0; k>kmax; ) {
    i = (int) (1.0 + theta*(double) k);
    g1 = gg((double) i);
    g2 = gg((double) (k+1-i));
    v = g1/(g1+g2);
    if (theta < v) {
      theta /= v;
      k = i-1;
    } else {
      x += i;
      theta = (theta-v)/(1.0-v);
      k -= i;
    }
  }

  for (i=0; i<k; i++) {
    if (u_random() < p) x++;
  }

  return x;
}
    
/* Multinomial random variates, k categories */

void r_multinom(int n, int k, double p[], int x[]) {

  int i, x0, xi;
  double pt, pi;
  double sum=0.0;
  
  for(i=0;i<k;i++){
    sum += p[i];
    
  }

  for(i=0;i<k;i++){
     p[i] /= sum;
    
  }
  
  for (i=k-1, x0 = n, pt=1.0; i>0; i--, x0 -= xi, pt -= pi) {
    pi = p[i];
    xi = r_binom(x0, pi/pt);
    x[i] = xi;
  }
  x[0] = x0;
}
