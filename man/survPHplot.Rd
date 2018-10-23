\name{survPHplot}
\alias{survPHplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot to assess non-proportionality
}
\description{
Plot of log(time) versus -log(-log(survival)) to assess
non-proportionality. A constant distance between curves suggest proportionality.
}
\usage{
survPHplot(formula, data, subset, contrasts, weights, col = 1:5,
           lty = 1:5, pch = 19, xlab = "Time (log scale)",
           ylab = "-log(-log(Survival))", log = "x",
           legend.args = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
either (i) formula with a \code{Surv} object on the left-hand-side and
stratification covariates on the right-hand-side, or (ii) a \code{survfit} object
}
  \item{data}{
data argument passed to \code{survfit}
}
  \item{subset}{
subset argument passed to \code{survfit}

}
  \item{contrasts}{
contrasts argument passed to \code{survfit}
}
  \item{weights}{
weights argument passed to \code{survfit}

}
  \item{col}{
colours of the curves passed to \code{lines}
}
  \item{lty}{
line type of the curves passed to \code{lines}
}
  \item{pch}{
pch for the curves passed to \code{points}
}
  \item{xlab}{
xlab graphics argument passed to \code{plot.default}
}
  \item{ylab}{
ylab graphics argument passed to \code{plot.default}
}
  \item{log}{
log graphics argument passed to \code{plot.default}
}
  \item{legend.args}{
list of arguments passed to \code{legend}. These arguments update the
  base arguments, which are
  \code{list(x="topright",legend=names(survfit$strata),col=col,lty=lty,pch=pch)} 
}
  \item{\dots}{
Other arguments passed to \code{plot.default}
}
}
\details{
The default plot is to use straight lines between the transformed
  survival values for each strata, rather than using steps.
}
\value{
Primary purpose is for plotting (side effect). The return value is
  initial plot.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
survPHplot(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594,
           data=colon, subset=(stage=="Localised"),
           legend.args=list(bty="n"))
}
\keyword{ survival }

