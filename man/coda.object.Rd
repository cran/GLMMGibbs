\name{coda.object}
\alias{coda.object}
\title{Convert Glmm Output For Use With CODA}
\description{Convert the result of a glmm fit for use
    with the CODA package}
\usage{
rcoda.object(g)
}
\arguments{
    \item{g}{The returned value from a
	call to the \code{glmm} function,
	an object of class \code{glmmfit}.}
}
\value{
    an object of class \code{mcmc}, which can be passed to
    functions in the \code{CODA} library}
\details{
    \code{CODA} is a package which can be used to check the convergence
    of output produced by the \code{BUGS} program.  The
    \code{rcoda.object} function converts an object returned by
    a call to the \code{glmm} into an object which can be used
    with the \code{CODA} package.
  The \code{glmm} function  yielding \code{g}
  must have been called with the argument
  \code{store.results} set to \code{TRUE}.}

    For full details, see the section
    "Using CODA and BOA with the
    output of \code{glmm}" of the document
    "\code{GLMMGibbs}: An R Package for Estimating Bayesian
    Generalised  Linear Mixed Models by Gibbs Sampling",
    supplied with this package.
}
\references{

    \url{http://www-fis.iarc.fr/coda/}

    \url{http://www.mrc-bsu.cam.ac.uk/bugs/welcome.shtml}
    
}
\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\seealso{\code{\link{glmm}}}
\keyword{models}





