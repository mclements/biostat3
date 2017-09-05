\name{eform}
\alias{eform}
\alias{irr}
\alias{or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the exponential form for coefficients and their confidence intervals using either profile likelihood-based or Wald-based confidence intervals.
}
\description{
\code{irr} and \code{or} use eform with a different name for the estimator.
}
\usage{
eform(object, parm, level = 0.95, method = c("Delta","Profile"), name = "exp(beta)")
irr(..., name = "IRR")
or(..., name = "OR")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted model object with \code{coef} and \code{confint} methods
}
  \item{parm}{
a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
}
  \item{level}{
the confidence level required
}
  \item{method}{
string to determine method to use the delta method (\code{stats::confint.default}), which assumes that the parameters are asymptotically normal, or profile likelihood-based confidence intervals (\code{MASS:::confint.gllm}), respectively.
}
  \item{name}{
name of the estimator.
}
\item{...}{
  arguments to pass from \code{irr} or \code{or} to \code{eform}.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
}
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

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
## from example(glm)
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
eform(glm.D93) # needs MASS to be installed
eform(glm.D93, method="Wald")  # based on asymptotic normality
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
