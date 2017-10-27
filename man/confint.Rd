% file MASS/man/confint.Rd
% copyright (C) 1999-2010 W. N. Venables and B. D. Ripley
%
\name{confint-MASS}
\alias{confint.glm}
\alias{confint.profile.glm}
\title{
Confidence Intervals for Model Parameters
}
\description{
  Computes confidence intervals for one or more parameters in a fitted
  model.  
}
\usage{
\method{confint}{glm}(object, parm, level = 0.95, trace = FALSE, \dots)

}
\arguments{
  \item{object}{
    a fitted model object.  Methods currently exist for the classes
    \code{"glm"}, \code{"nls"} and for profile objects from these classes.
  }
  \item{parm}{
    a specification of which parameters are to be given confidence
    intervals, either a vector of numbers or a vector of names.  If
    missing, all parameters are considered.
  }
  \item{level}{
    the confidence level required.
  }
  \item{trace}{
    logical.  Should profiling be traced?
  }
  \item{\dots}{
    additional argument(s) for methods.
  }
}
\value{
  A matrix (or vector) with columns giving lower and upper confidence
  limits for each parameter.  These will be labelled as (1 - level)/2 and
  1 - (1 - level)/2 in \% (by default 2.5\% and 97.5\%).
}
\details{
  \code{\link[stats]{confint}} is a generic function in package \code{stats}.
  
  These \code{confint} methods call the appropriate profile method,
  then find the confidence intervals by interpolation in the profile
  traces.  If the profile object is already available it should be used
  as the main argument rather than the fitted model object itself.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\seealso{
  \code{\link[stats]{confint}} (the generic and \code{"lm"} method),
  \code{\link{profile}}
}
\examples{
confint(glm(chd ~ hieng + offset(log(y)), data=diet, family=poisson))
}
\keyword{models}
