\name{print.glmmfit}
\alias{print.glmmfit}
\title{Print Result Of Glmm Fit}
\description{Print out information about a glmm fit.}
\usage{print.glmm(g)}
\arguments{\item{g}{The result of a call to the \code{glmm} function,
    an object of class \code{glmmfit}.}}
\details{The call made, the number of steps, and the sample
    statistics (mean, standard deviation, and 2.5, 25, 50,
    75 and 97.5) of all the parameters of interest in the model
    are printed out. Sample statistics for hyperparameters are
    plotted on both the tau and sigma scales.}
    The \code{glmm} function  yielding \code{g}
  must have been called with the argument
  \code{store.results} set to \code{TRUE}.}
\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\seealso{ \code{\link{glmm}}}
\keyword{models}
