\name{survrate}
\alias{survrate}
\title{
Describe rates
}
\description{
Describe rates using the \code{\link{Surv}} function.
}
\usage{
survrate(formula, data, subset, addvars = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
formula with the left-hand-side being a \code{\link{Surv}} function and
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
(default=FALSE). This is useful for subsequent analysis.
}
  \item{\dots}{
other arguments to the \code{\link{poisson.test}} function for
calculation of the confidence intervals.
}
}
%% \details{
%% }
\value{
data-frame.
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
survrate(Surv(y,chd) ~ hieng, data=diet)
}
\keyword{ survival }

