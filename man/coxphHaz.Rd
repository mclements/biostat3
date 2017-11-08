\name{coxphHaz}
\alias{coxphHaz}
\alias{print.coxphHaz}
\alias{plot.coxphHaz}
\alias{plot.coxphHazList}
\alias{lines.coxphHazList}
\title{
Smoothed hazard estimates for \code{coxph}
}
\description{
Smoothed hazard estimates for \code{coxph}
}
\usage{
coxphHaz(object, newdata, n.grid = 300, kernel = "epanechnikov", from,
to, ...)
\method{print}{coxphHaz}(x, digits=NULL, ...)
\method{plot}{coxphHaz}(x, xlab="Time", ylab="Hazard", type="l", ...)
\method{plot}{coxphHazList}(x, xlab="Time", ylab="Hazard", type="l", col=1:length(x), lty=1, legend.args=list(), ...)
\method{lines}{coxphHazList}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{coxph} object
}
  \item{newdata}{
data-frame with covariates for prediction
}
  \item{n.grid}{
the number of grid values for which the hazard is calculated
}
  \item{kernel}{
the kernel used for smoothing
}
  \item{from}{
argument for \code{density}. Defaults to the minimum time.
}
  \item{to}{
argument for \code{density}. Defaults to the maximum time.
}
  \item{x}{
object
}
  \item{digits}{
argument passed to \code{print.density}
}
  \item{col}{
graphics argument
}
  \item{lty}{
graphics argument
}
  \item{xlab}{
graphics argument
}
  \item{ylab}{
graphics argument
}
  \item{type}{
graphics argument
}
  \item{legend.args}{
    a list of options that are passed to the legend call. Defaults are
    \code{list(x="topright",legend=strata(attr(x,"newdata")),col=col,lty=lty)}.
}
  \item{\dots}{
other arguments. For \code{coxphHaz}, these arguments are passed to
\code{density}. For the \code{plot} and \code{lines} methods, these are
passed to the relevant \code{plot}, \code{matplot} and \code{matlines} functions.
}
}
\details{
Smooth hazard estimates from a Cox model using kernel smoothing of the
Nelson-Aalen estimator.
}
\value{
The \code{coxphHaz} function returns either a class of type
\code{c("coxphHaz","density")} when \code{newdata} has one row or, for multiple rows in
\code{newdata}, a class of type \code{"coxphHazList"}, which is a list of
type \code{c("coxphHaz","density")}.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

\seealso{
 \code{\link{coxph}}, \code{\link{survfit}}, \code{\link{density}}
}
\examples{
fit <- coxph(Surv(surv_mm/12,status=="Dead: cancer")~agegrp, data=colon)
newdata <- data.frame(agegrp=levels(colon$agegrp))
haz <- suppressWarnings(coxphHaz(fit,newdata))
plot(haz, xlab="Time since diagnosis (years)")
legend("topright", legend=newdata$agegrp, lty=1:4, col=1:4, bty="n")
}
\keyword{ survival }

