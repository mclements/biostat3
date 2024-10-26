\name{survRate}
\alias{survRate}
\title{
Describe rates
}
\description{
Describe rates using the \code{\link[survival]{Surv}} function.
}
\usage{
survRate(formula, data, subset, addvars = TRUE, ci=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
formula with the left-hand-side being a \code{\link[survival]{Surv}} function and
the right-hand-side being any stratification variables.
}
  \item{data}{
source dataset
}
  \item{subset}{
subset conditions for the source dataset
}
  \item{addvars}{
logical for whether to add the stratification variables to the output
(default=TRUE). This is useful for subsequent analysis.
}
  \item{ci}{
logical for whether to calculate the confidence interval (default=TRUE). 
}
  \item{\dots}{
other arguments to the \code{\link{poisson.test}} function for
calculation of the confidence intervals.
}
}
%% \details{
%% }
\value{
data-frame with columns \code{tstop}, \code{event}, \code{rate},
\code{lower} and \code{upper}. Covariates are appended if
\code{addvar=TRUE}.

Confidence intervals use \code{stats::poisson.test}.
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
## incidence rates for CHD for low- or high-energy diets
survRate(Surv(y,chd) ~ hieng, data=diet)
}
\keyword{ survival }

