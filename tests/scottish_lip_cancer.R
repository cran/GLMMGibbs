library(GLMMGibbs)

#Fit the data from Breslow and Clayton (1993) on the number 
#of cases of lip ccaner in the 56 (pre-reorganisation)
#counties of scotland, in generalised linear mixerd models
#which either do or do not include the precentage of workers employed
#in Agriculture, Forestry and Fishing, and which contain a county
#effect which represent either heterogeneous extra-poisson variability
#or spatially correlated errors.


#load in the data from the glmm library
data(scottish.lip.cancer)


scottish.lip.cancer$AFF10 <- scottish.lip.cancer$AFF/10
#For compatibility with Brreslow and Clayton.  Notice that
#we have to do it like this because arithmetical operations
#on covariates are not allowed within glmm.


#define two random effects, county.het modelling
#heterogeneous extra-poisson variability and county.spatial
#modelling spatially correlated varriability
scottish.lip.cancer$county.spatial <- Ra(as.factor(1:56)                      ,
                                shape=0.0001,scale=0.0001,
                                type="adj.map",
                                map="scotland.adj",
                                contrast="sum")

scottish.lip.cancer$county.het <- Ra(as.factor(1:56),
                                shape=0.0001,scale=0.0001,
                                contrast="sum")


#
#Fit the models by Gibbs Sampling, taking 2000 steps of
#ICM towards an approximate maximum of the posterior, 4000
#steps of burnin, and then another 40000 steps, saving every 10.
#
g.scotland.1 <- glmm(
    formula = observed ~  county.het + offset(log(expected)), 
    family = poisson,
    data = scottish.lip.cancer,
    icm = 2000, 
    burnin = 4000,
    keep = 40000,
    progress.info = 1000,
    thin = 10)

g.scotland.2 <- glmm(
    formula = observed ~  county.spatial + offset(log(expected)), 
    family = poisson,
    data = scottish.lip.cancer,
    icm = 2000, 
    burnin = 4000,
    keep = 40000,
    progress.info = 1000,
    thin = 10)

g.scotland.3 <- glmm(
    formula = observed ~ AFF10 + county.het + offset(log(expected)), 
    family = poisson,
    data = scottish.lip.cancer,
    icm = 2000, 
    burnin = 4000,
    keep = 40000,
    progress.info = 1000,
    thin = 10)

g.scotland.4 <- glmm(
    formula = observed ~ AFF10 + county.spatial + offset(log(expected)), 
    family = poisson,
    data = scottish.lip.cancer,
    icm = 2000, 
    burnin = 4000,
    keep = 40000,
    progress.info = 1000,
    thin = 10)










