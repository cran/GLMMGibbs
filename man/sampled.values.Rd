\name{sampled.values}
\alias{sampled.values}
\title{Sampled Values from a glmm fit}
\description{Returns a vector consisting of the sampled values
  from an object returned from a call to the glmm function}
\usage{ sampled.values(g,parameter.name)}

\arguments{{g}{The result of a call to the \code{glmm} function,
    an object of class \code{glmmfit}.}
  {parameter.name} The name of a parameter in the model}
\details{The function takes and object \code{g}, the returned value
  of a call to the \code{glmm} function, and the name of a parameter
  in the model fitted, and returns a vector of the sampled values.
  The \code{glmm} function  yielding \code{g}
  must have been called with the argument
  \code{store.results} set to \code{TRUE}.}
\seealso{ \code{\link{glmm}} }
\keyword{models}  

