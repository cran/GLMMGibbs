#include <stdio.h>

#include "sparsev.h"
#include "sparsem.h"
#include "misc.h"
#include <math.h>

     static  int doublecompare(const void *ii, const void *jj)
          {
	    double *i, *j;
	    i = (double *) ii;
	    j = (double *) jj;
                  if (*i > *j)
                          return (1);
                  if (*i < *j)
                          return (-1);
                  return (0);
          }




     static  int intcompare(const void *ii, const void *jj)
          {
	    int *i, *j;
	    i = (int *) ii;
	    j = (int *) jj;
                  if (*i > *j)
                          return (1);
                  if (*i < *j)
                          return (-1);
                  return (0);
          }




double id_1(double x){
  return 1.0;
}

double id_2(double x){
  return 1.0;
}

double id1(double x, double y){
  return 1.0;
}

double logit(double z)
     /* Logit link (inverse transformation) */
{
  double w;

  
  return log(z/(1-z));
}

double logit_1(double z)
     /* Logit link (inverse transformation) */
{
  double w;

  
  if (z > 30.0) z = 30;
  if (z < -30.0) z = -30;
  w = exp(z);
  return w/(1.0+w);
}

double logit_2(double z)
     /* Logit link (first derivative) */
{
  return 1.0/(z*(1.0-z));
}



double log_1(double z)
     /* Log link (inverse transformation) */
{
  if (z > 30.0) z = 30;
  if (z < -30.0) z = -30;
  return exp(z);
}

double log_2(double z)
     /* Log link (first derivative) */
{
  return 1.0/z;
}

double bern_v(double mu)
     /* Bernoulli errors (variance function) */
{
  return mu*(1.0-mu);
}

double bern_d(double y, double mu)
     /* Bernoulli errors (deviance contribution) */
{
  double dv=0.0;
  if (y!=0.0) dv -= 2*(y*log(mu/y));
  if (y!=1.0) dv -= 2*(1-y)*log((1-mu)/(1-y));
  return dv;
}


double poisson_v(double mu)
     /* Poisson errors (variance function) */
{
  return mu;
}

double poisson_d(double y, double mu)
     /* Poisson errors (deviance contribution) */
{
  double res;
  
  res = mu - y;
  if (y > 0.0  && mu > 0.0) res += y*log(y/mu);
  return 2.0*res;
}



char *myitoa(int n)
{
  int i=1;
  int tenpowi=10;
  char *res;
  
  while (n>=tenpowi)
	 {
	   tenpowi*=10;
	   i++;
	 }
	 
  i++;
  
  res=(char *) calloc(i,(size_t) sizeof(char));
  res[(i-1)] = '\0';
  
  i-=2;
  
  do
    {
      res[i--] = n % 10 + '0';
    }
  while ((n/=10) > 0);
  
  return res;
  
  
      
}


char *myconcat(char *a, char *b){
  int len;
  char *res;
  int i,la,lb;
  
  len = (la = strlen(a)) + (lb = strlen(b));
  res = (char *) calloc ( (len + 1), (size_t) sizeof(char)); 
  for(i=0;i<la;i++)  res[i]=a[i];
  for(i=0;i<=lb;i++) res[(i+la)]=b[i];  

  return res;
}


char **level_names(char *factor_name,int i,int j)
{
  char **res;
  int k;
  
  res = (char **) calloc ((j-i+1), (size_t) sizeof(char *));
  for(k=i;k<=j;k++)
    {
      res[(k-i)]=myconcat(factor_name,myitoa(k));
      
    }
  return res;
  
}

char **level_names_zlevel(char *factor_name,int i,int j,int zlevel)
{
  char **res;
  int k;
  
  res = (char **) calloc ((j-i+1), (size_t) sizeof(char *));
  for(k=i;k<zlevel;k++)
    { 
      res[(k-i)]=myconcat(factor_name,myitoa(k));
      
    }
  for(k=zlevel+1;k<=j;k++)
    { 
      res[(k-i-1)]=myconcat(factor_name,myitoa(k));
    }
  return res;
  
}

char *my_strcopy(char *str)
{
  char *res;
  int i=0;
  
  res = (char *) calloc((strlen(str) + 1), (size_t) sizeof(char));
  do
    {
      res[i]=str[i];
      i++;
      
    }
  while(str[i]!='\0');
  
  return res;
}


char **string_product(char **a,
		      char **b,
		      int la,
		      int lb)
{
  char **res;
  int i,j,k;

  res = (char **) calloc(la*lb, (size_t) sizeof (char *));
  if (!res) goto noroom;
  k=0;
  
  for(i=0;i<la;i++)
    {
      for(j=0;j<lb;j++)
	{
	  res[k] = myconcat(a[i],b[j]);
	}
    }
  return res;
  

 noroom:
  fprintf(stderr,"*** string_product *** dynamic storage overflow\n");
  exit(0);
  
} 

long *longp(long i){
  long *res;
  if (!(res)) goto noroom;
  
  res = (long *) malloc((size_t) sizeof(long));
  res[0] = i;
  return res;

 noroom:
  fprintf(stderr,"*** longp *** dynamic storage overflow\n");
  exit(0);
}


double *ones(int n)
{
  double *res;
  int i;
  
  res=(double *) calloc(n, (size_t) sizeof (double));
  if (!(res)) goto noroom;
  for(i=0;i<n;i++){
    res[i]=(double) 1.0;
      }

  return res;

 noroom:
  fprintf(stderr,"*** ones *** dynamic storage overflow\n");
  exit(0);
  
}

int *iones(int n)
{
  int *res;
  int i;
  
  res=(int *) calloc(n, (size_t) sizeof (int));
  if (!(res)) goto noroom;
  for(i=0;i<n;i++){
    res[i]=(int) 1;
      }

  return res;

 noroom:
  fprintf(stderr,"*** iones *** dynamic storage overflow\n");
  exit(0);  

}
