#Fit a Generalise Linear Mixed Model to the datadescribed
#by Thall and Vail (1990) on a clinical trial of progadine.
library(GLMMGibbs)


#load the data
data(epil.bugs)

#create a  data frame having one row for each visit,
#using log(baseline seizures/4) and log(age) as
#covariates, and centring these about their mean.



### k copies of each element of x in turn.
inner.rep <- function(x,k){rep(x,rep(k,(length(x))))}

centre <- function(x){x-mean(x)}

epil <- data.frame(y=epil.bugs$y,
                   base=centre(log(inner.rep(epil.bugs$Base/4,4))),
                   age=centre(log(inner.rep(epil.bugs$Age,4))),
                   trt=inner.rep(epil.bugs$Trt,4),
                   v4=rep((c(0,0,0,1)),59),
                   subject=as.factor(inner.rep(1:59,4))
                   )



#declare the  factor subject to represent a random effect.

epil$subject <- Ra(epil$subject,shape=0.0001,scale=0.0001,
                  contrast="treatment",of.interest=FALSE)


#
#Fit the model by Gibbs Sampling, taking 2000 steps of
#ICM towards an approximate maximum of the posterior, 4000
#steps of burnin, and then another 40000 steps, saving every 10.
#

g.epil  <- glmm(y~base+age+trt+base*trt+v4+subject,
            family=poisson,
            data=epil,
            icm=2000,
            burnin=2000,
            keep=40000,
            store.results=TRUE,
            progress.info=1000,
                thin=10)


g.epil
                   







