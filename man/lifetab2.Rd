\name{lifetab2}
\alias{lifetab2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula wrapper for \code{\link{lifetab}} from the
\code{\link{KMsurv}} package.
}
\description{
Calculate a life table using the actuarial method using a formula and a
data-frame with optional breaks. 
}
\usage{
lifetab2(surv, data = NULL, breaks = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv}{
    formula for the \code{Surv} object, including a time and event indicator.
}
  \item{data}{
    optional \code{data.frame} for the \code{Surv} object. If this is not provided, then the parent frame
    is used for the \code{Surv} object.
}
  \item{breaks}{
    optional numeric vector of breaks. If this is not provided, then the
    unique time values from the \code{Surv} object are used together
    with \code{Inf}.
}
}
\details{
See \code{lifetab} for details. This wrapper is meant to make
life easier.
}
\value{
  A \code{data.frame} as per \code{\link{lifetab}}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark Clements for the wrapper.
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
## we can use unique transformed times
lifetab2(Surv(floor(surv_yy),status=="Dead: cancer"), colon_sample)

## we can also use the breaks argument
lifetab2(Surv(surv_yy,status=="Dead: cancer"), colon_sample, breaks=0:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")

