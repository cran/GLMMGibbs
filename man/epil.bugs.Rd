\name{epil.bugs}
\alias{epil.bugs}
\title{Results of a Clinical trial of Progadine}
\description{The dataset describes the results of the pre-crossover
  phase of a crossover trial of the effects of
  progadine, an anti-epileptic drug.
  Patients were assigned blindly and
  at random to receive either the drug or a placebo in 
  addition to standard chemotherapy.  Each patient visited
  the trial centre on four subsequent occasions,  on each of
  which the number of seizures in the previous two weeks was
  recorded. Also recorded was the patient`s age, whether they
  were on treatment or placebo and the number of seizures
  suffered in the two-week period before the trial (known as
  the baseline seizures). The data is stored in the same format as in
  the BUGS example book}
\usage{data(epil.bugs)}  
\format{ \tabular{rlll}{
    [,1] \tab y      \tab numeric(236) \tab number of seizures
    (one observation per visit) \cr
    [,2] \tab Trt    \tab numeric(59)  \tab whether the patient
    received Progadine (=1) or placebo (=0)
    (one observation per patient) \cr
    [,3] \tab Base    \tab numeric(59)    \tab number of baseline seizures
    (one observation per patient) \cr
    [,4] \tab Age     \tab numeric(59) \tab Age
    (one observation per patient) \cr
  }
}
\source{Thall, P . F. and Vail, S. C. 1990,
  Some covariance  models for longitudinal count data with
  overdispersion
  \emph{Biometrics}, pages 657-671, volume 46.}
\reference{
  Spiegelhalter, D. J., Thomas, A., Best, N. G., and Gilks,
  W. R. (1995)
  \emph{BUGS: Examples, Version 0.50 Volume 1. MRC Biostatistics Unit, Cambridge.} 
}    
\keywords{datasets}





