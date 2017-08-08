\name{muhaz2}
\alias{muhaz2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula wrapper for the \code{\link{muhaz}} function from the
\code{muhaz} package.}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
muhaz2(formula, data, subset, ...)
plot.muhazList(x, lty=1:length(x), col=1:length(x), ...)
lines.muhazList(x, lty=1:length(x), col=1:length(x), ...)
summary.muhazList(object)
as.data.frame.muhazList(x)
as.data.frame.muhaz(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    formula with the left-hand side being a \code{Surv} object, including a time and event indicator, and the right-hand side indicated stratification.
}
  \item{data}{
    optional \code{data.frame} for the \code{Surv} object. If this is not provided, then the parent frame
    is used for the \code{Surv} object.
}
  \item{subset}{
    optional \code{subset} statement
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For a single strata, this is a \code{\link{muhaz}} object. For multiple strata, this is a \code{muhazList} object, which includes methods for 
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
% 
% %% ~Make other sections like Warning with \section{Warning }{....} ~
% 
% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
  plot(muhaz2(Surv(surv_mm, status == "Dead: cancer")~1, melanoma))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
