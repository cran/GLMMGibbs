\name{bugs.out}
\alias{bugs.out}
\title{Write Glmm Fit  To Bugs Output Files }
\description{\code{bugs.out} writes out the returned value from a call to the
       \code{glmm} function to two files,
       with extensions \code{.ind} and \code{.out}, in
       \code{BUGS} output format.}
\usage{
bugs.out(g,file)
}
\arguments{
    \item{g}{The result of a call to the \code{glmm} function,
	an object of class \code{glmmfit}.}
    \item{file}{The primary filename of the files to which
	the the object is to be written.}
}
\details{The sampled values and parameter names of the fitted
    model object are written to two files, \code{"filemame.ind"}
    and \code{"filename.out"}, where \code{"filename"} is the
    value of the character object argument
    \code{file} . The format of these files is
    described    in the section ``Using CODA and BOA with the
    output of \code{glmm}'' of the document
    ``\code{GLMMGibbs}: An R Package for Estimating Bayesian
    Generalised  Linear Mixed Models by Gibbs Sampling'',
    supplied with this package. }
    The \code{glmm} function  yielding \code{g}
  must have been called with the argument
  \code{store.results} set to \code{TRUE}.}

\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\seealso{ \code{\link{glmm}}, \code{\link{plot.glmm}},
          \code{\link{print.glmm}}
\keyword{models}




