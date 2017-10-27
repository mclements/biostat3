\name{plot.bshazard}
\alias{plot.bshazard}
\alias{lines.bshazard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for 'bshazard' (minor update of the same function in the
\code{bshazard} package) }
\description{
A plot of hazard rate is produced. The overall option allows to plot an
hazard rate for each covariate value (assuming proportional hazard).
The updated function now passess \code{...} to \code{plot.default}.
}
\usage{
\method{plot}{bshazard}(x,conf.int=T,overall=T,col,lwd,xlab="Time",ylab="Hazard rate", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 	the result of a call to the bshazard function. 
}
  \item{conf.int}{
Determines whether confidence intervals will be plotted. The default is to do so if there is only 1 curve, i.e., no strata. }
  \item{overall}{
Determines whether an overall curve is plotted (default overall=T) or a curve for each covariate value in the data (overall=F). 
It works only if there are covariates. }
  \item{col}{
  a vector of integers specifying colors for each curve. The default value is 1. }
  \item{lwd}{
a vector of integers specifying line types for each curve. The default value is 1. }
  \item{xlab}{
 	label given to the x-axis. }
  \item{ylab}{
 	label given to the y-axis. }
  \item{\dots}{
    other arguments to \code{plot.default} (now correctly passed through)
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  bshazard,summary.bshazard,print.bshazard}
\examples{
library(bshazard)
fit<-bshazard(Surv(time, status==2) ~ 1,data=lung)
plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }

