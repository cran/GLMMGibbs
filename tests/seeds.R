library(GLMMGibbs)


#These commands use the R functions Ra and glmm to
#etsimate the posterior distributions of the parameters
#in models for Crowder's (1979) seeds data.

#load in the seeds data from the glmm library

data(seeds)

#declare the plate factor to be a randoom effect with
#a hyperprior with a gamma distribution with shape and scale
#parameters both euqal to 0.001.
   seeds$plate <- as.factor(1:21)

   seeds$plate <- Ra(seeds$plate,
	shape=0.001,
        scale=0.001,
        type="identity",
        contrast="sum",
	of.interest=F)


#Fit the two models by Gibbs sampling.  In each case, take 500 steps
# of ICM towards the posterior maximum, then do 1000 setps of Gibbs sampling
# as "burn-in" and  20000 iterations, storing every fifth.



g.seeds.1 <- glmm(formula=cbind(r,n-r)~x1+x2+plate,
     data=seeds,
     family=binomial,
     store.results=T,             
     icm=500,
     burnin=1000,
     keep=20000,
     progress.info=500,
     thin=5)


g.seeds.2 <- glmm(formula=cbind(r,n-r)~x1*x2+plate,
     data=seeds,
     family=binomial,
     store.results=T,             
     icm=500,
     burnin=1000,
     keep=20000,
     progress.info=500,
     thin=5)

g.seeds.1

g.seeds.2



