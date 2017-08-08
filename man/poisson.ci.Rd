\name{poisson.ci}
\alias{poisson.ci}
\title{
Exact Poisson confidence intervals.
}
\description{
A wrapper for the \code{\link{poisson.test}} that allows for vector values.
}
\usage{
poisson.ci(x, T = 1, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
number of events.}
  \item{T}{
time base for event count. 
}
  \item{conf.level}{
confidence level for the returned confidence interval.
}
}
\details{
This uses \code{stats::poisson.test} for the calculations.
}
\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{poisson.test}}
}
\examples{
### These are paraphrased from data sets in the ISwR package

## SMR, Welsh Nickel workers
poisson.ci(137, 24.19893)

## eba1977, compare Fredericia to other three cities for ages 55-59
poisson.ci(c(11, 6+8+7), c(800, 1083+1050+878))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%_\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
