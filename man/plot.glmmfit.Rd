\name{plot.glmmfit}
\alias{plot.glmmfit}
\title{Plot Result Of Glmm Fit}
\description{Plot either the trace plots, the histograms,
    or both, of the sampled values for  parameters of
    interest in a model fitted by \code{glmm}}
\usage{
plot.glmmfit(g, parameters, type = "b", hscale = "b",
plots.per.page = 4, postscript)
}
\arguments{
    \item{g}{The result of a call to the \code{glmm} function,
	an object of class \code{glmmfit}.}
    \item{parameters}{A vector of character strings giving names
	parameters in the model.  If not missing, only sampled
	values for these parameters are plotted.  If missing,
	sampled values for all parameters of interest are plotted.}
    \item{type}{A character, determining which type of
	plot is given. If \code{"t"}, \code{"T"}, \code{trace} or
	\code{Trace} a trace plot is made. If \code{"h"}, \code{"H"},
	\code{"hist"} or \code{"Hist"} a histogram is plotted.
	If \code{"b"}, \code{"B"}, \code{"both"} or \code{"Both"} both
	plots are made, side by side, with the trace plot on the
	left and the histogrm on the right. Default is \code{"b"}.}
    \item{hscale}{A character, determining the scale on which
	scale sampled values of hyperparameter values are plotted.
	If \code{"t"}, \code{"T"}, \code{"tau"} or \code{"Tau"} sampled
	values are plotted on the tau scale. If \code{"s"}, \code{"S"}
	\code{"sigma"} or \code{"Sigma"} values are plotted on the
	sigma scale.  If \code{"b"}, \code{"B"}, \code{"both"} or
	\code{"Both"}, sampled values are plotted on both scales.
	Default is \code{"b"}.}
    \item{plots.per.page}{The number of variables plotted per
	page of output. Default is 4.}
    \item{postscript}{A character string. If it is not missing and
	takes the value \code{"file"}, say, the output is written
        to files \code{"file_page1.ps"}, \code{"file_page2.ps"}, etc. }
}
\details{The sampled values of some or all of the parameters of
    interest in a model fitted
    by the \code{glmm} function are plotted out. The parameters
    for which values are plotted, whether trace plots, histograms,
    or both are plotted, the scale on which hyperparameter values
    are plotted, the number of variables plotted per page and
    whether the plots are made to the screen or to postscript files
    are all decided by the arguments to the function}
    The \code{glmm} function  yielding \code{g}
  must have been called with the argument
  \code{store.results} set to \code{TRUE}.}
\author{ Jonathan Myles, Imperial Cancer Research Fund, and David
  Clayton, Wellcome Trust \email{mylesj@icrf.icnet.uk}}
\seealso{  \code{\link{glmm}} }
\keyword{models}
    
	
	



