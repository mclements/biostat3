\name{addIndicators}
\alias{addIndicators}
\title{
Utility to add indicators from a data-frame based on a formula.
}
\description{
Column-bind a model matrix to the source data-frame
}
\usage{
addIndicators(data, formula, drop.intercept = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
source data-frame or matrix.
}
  \item{formula}{
model formula used to add columns.
}
  \item{drop.intercept}{
logical as to whether to drop the column named '(Intercept)'.
}
}
\details{
This function calls \code{model.matrix}, conditionally checks for and removes '(Intercept)', and binds with the original data-frame (or matrix).
}
\value{
data-frame or matrix.
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
addIndicators(data.frame(f = c("a","a","b")), ~f+0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
