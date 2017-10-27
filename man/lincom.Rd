\name{lincom}
\alias{lincom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear combination of regression parameters.
}
\description{
Using results calculated by the \code{\link{linearHypothesis}} function in
the \code{car} package, calculate a linear combination of regression parameters.
}
\usage{
lincom(model, specification, level = 0.95, eform = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
regression model object (as per the \code{model} argument in \code{\link{linearHypothesis}})
}
  \item{specification}{
    specification of the linear combination. This is the same as a
    single component of the \code{hypothesis.matrix} argument in \code{\link{linearHypothesis}}.
}
  \item{level}{
the confidence level required
}
  \item{eform}{
logical for whether to exponentiate the confidence interval (default=FALSE)
}
  \item{\dots}{
other arguments to the \code{\link{linearHypothesis}} function.
}
}
\details{
Multiple specifications of linear combinations are called individually.
}
\value{
  A matrix with columns including the estimate, a normal-based confidence interval,
  test statistic and p-values.
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

\seealso{
See Also \code{\link{linearHypothesis}}.
}
\examples{
  fit <- glm(chd ~ hieng*job + offset(log(y)), data=diet, family=poisson)
  lincom(fit, c("hienghigh+hienghigh:jobconductor",
                "hienghigh+hienghigh:jobbank"),
         eform=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
