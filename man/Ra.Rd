\name{Ra}
\alias{Ra}
\title{Create Random Factor}
\description{ Create a random factor to be used in the \code{glmm}
    function for  fitting generalised linear mixed models by
    Gibbs sampling}
\usage{Ra(data, shape=0.001, scale=0.001, type="identity",
contrast="identity", map="none", zlevel=1, of.interest=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
    \item{data}{ The values of the factor which is to be declared
	random. }
    \item{shape}{ The shape parameter of the gamma prior
	distribution to be placed on the hyperparameter. }
    \item{scale}{ The scale parameter of the gamma prior
	distribution to be placed on the hyperparameter. }
    \item{type}{ The ``type'' of random effect.  Current
	possible options are \code{"identity"} which creates an
	exchangeable random effect,  \code{"adj.map"} which
	creates a spatially smoothed effect, assuming the
	levels to represent the values of the effect. The
	default value is \code{"identity"}.   }
    \item{type}{ A character string representing the way in
	which the factor is to be encoded in the model matrix.
	The possible values are \code{"identity"}, which causes a
	``dummy'' coding to be used, \code{"treatment"} which
	causes a coding which sets the first level of the
	factor to zero, \code{sum} which forces the effects to
	sum to zero, and {\tt helmert} the  S-plus default
	coding mechanism. }
    \item{map}{ If the factor is defined to be of type {\tt adj.map},
        the full path name to a file containing information about
	which levels of the factor represent adjacent regions, in
	a format described below. }    
    \item{of.interest}{ if \code{TRUE}, all values of each level
	of the effect are output during the ``keep'' phase of
	sampling (see the instructions for {\tt glmm} below).
	If \code{FALSE}, only the sampled values of the
	hyperparameter are kept.}
}
\details{
  \code{Ra} returns an a factor, which when included in the model
  passed to the \code{glmm}
  function, is treated as a random effect.

  The effect values have a multivariate normal prior with mean zero
  and a covariance matrix which depends on the `type' of random
  effect and the parameterisation being used.  Its inverse is the
  product of a fixed matrix and a random scalar hyperparameter, which has
  a gamma prior distribution with given shape and scale parameters.

  
  Full details are given in the document
  "\code{GLMMGibbs}: An R Package for Estimating Bayesian Generalised
  Linear Mixed Models by Gibbs Sampling", supplied with this package.
}
\value{
   A value of class \code{factor}, but with additional attributes.
}
\references{
    
    Clayton, D.G. (1996) Generalized Linear Mixed Models in
    \em Ph{Markov chain Monte Carlo in Practice}, ed. Gilks, W. R. and
    Richardson, S. and  and Spiegelhalter, D. J., Chapman & Hall.
    
}
\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\note{ if the factor is set to be of type \code{adj.map}, the argument
  \code{mapfile} must be set to the name of a file which describes which
  areas are adjacent to which.  The format of the file is as follows:
  each area should be assigned a number, and the file should then have
  successive lines representing the consecutive areas, each consisting of
  the number of neighbours of the area, followed by a colon, followed
  by the numbers of the neighbouring  areas, followed by a 

    
  Thus, for example, if area 1 has three neighbours, numbered 7, 11, and
  12, the first line of the file should be:

  \code{3:7,11,12}

  Compare the file \code{scotland.adj} with the data given in the
  document
   \code{GLMMGibbs}: An R Package for Estimating Bayesian Generalised
  Linear Mixed Models by Gibbs Sampling, supplied with this package,
  for a further example.}



\seealso{  \code{\link{glmm}}}
\keyword{models}








