\name{glmm}
\alias{glmm}
\title{GLMMs By Gibbs Sampling}
\description{
    \code{glmm} estimates the posterior distribution of the fixed
    effects of a
    Generalised Linear Mixed Model.  It also
    estimates the hyperparameter related to each random effect, and the
    effect values of each random effect declared to be ``of interest'' by
    function \code{Ra}.}

\usage{
 glmm(formula, family, data, weights, offset, icm = 50,
      burnin = 1000, keep = 1000, model.show = FALSE,
      progress.info = 0, store.results = FALSE,  thin = 1,
      bugsfile, seed1 = 6762, seed2 = 92928, seed3 = 19729)
}
\arguments{
 \item{formula}{ The model formula for the model to be fitted.}
 \item{family}{A description of the error distribution and link
     function to be used. At present the available
     options are  \code{"binomial"} and \code{"poisson"}.  The
     canonical link must always be used, so there is no
     \code{link} argument. }
 \item{data}{The name of a data frame in which the data is stored.
     (at present, the function can only be used with the data
      in a data frame)}
 \item{weights}{ A vector of weights.}
 \item{icm}{The number of steps of deterministic maximisation of
     the posterior distribution  by the Iterative Conditional Mode
     algorithm before sampling begins.  All hyperparameters are
     fixed to 1.0 . }
 \item{burnin}{ The number of steps of Gibbs sampling taken
     before storage of the parameter values begins. }
 \item{keep}{ The number of steps of Gibbs sampling taken once 
     storage of the parameter values has begun. }
 \item{thin}{An integer \eqn{t} specifying that one storage of
     the parameter values has begun every \eqn{(t)} th iteration is
     saved}
 \item{model.show}{ A debugging argument used by he developers
     which will be deleted at the final release.}
 \item{progress.info}{ If an integer, \eqn{n} say, the functions
     reports when every \eqn{n} iterations have been carried out}
 \item{bugsfile}{A character string.  If not missing,
     the function creates two files, with this string
     their primary filenames and \code{.ind} and \code{.out}
     their secondary filenames  (extensions)
     .These contain the output in {\tt BUGS }
     format. }
 \item{store.results}{if \code{TRUE}, returns the sampled values of the
     parameters (see ``Value'' below) }
 \item{seed1}{The first of three seeds supplied to the Random Number
     Generator used by the underlying \code{C} code} 
 \item{seed2}{The second seed.}
 \item{seed3}{The third seed.}
}
\details{
    Generalised Linear Mixed Models (GLMMs) are an extension
    of GLMs with the addition of ``random effects''  given whose
    values the response values are conditionally independent.
    The function \code{glmm} fits these models in a Bayesian
    paradigm by Gibbs sampling.

    Full details are given in the document
    "\code{GLMMGibbs}: An R Package for Estimating Bayesian
    Generalised Linear Mixed Models by Gibbs Sampling", supplied
    with this package.
    }
\value{
    an object  of class \code{glmmfit}, which contains the sample
    statistics of the sampled values and, if the
    \code{store.results} argument is set to \code{TRUE}, the
    sampled values themselves.}
\references{
    
    Clayton, D.G. (1996) Generalized Linear Mixed Models in
    \emph{Markov chain Monte Carlo in Practice}, ed. Gilks, W. R. and
    Richardson, S. and  Spiegelhalter, D. J., Chapman & Hall.

}
\note{ \code{GLMMgibbs} (the package from which glmm comes) is a beta release
    and we strongly recommend the use of \code{save.image()}
    before \code{glmm()} is used}
\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\seealso{ \code{\link{Ra}},
  \code{\link{plot.glmm}},
  \code{\link{print.glmm}}
\keyword{models}


