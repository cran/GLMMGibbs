#include <string.h>



#define MYTOL 1e-08
#define LOG2  0.6931472
double id_1 (double);
double id_2 (double);
double id1 (double, double);
double logit(double);

double logit_1(double );
double logit_2(double );
double log_1(double );
double log_2(double );
double bern_v(double );
double bern_d(double , double );
double poisson_v(double );
double poisson_d(double , double );


char *myitoa(int n);
char *myconcat(char *a, char *b);
char **level_names(char *factor_name,int i,int j);
char **level_names_zlevel(char *factor_name,int i,int j,int zlevel);
char *my_strcopy(char *str);
char **string_product(char **, char **, int,int);


