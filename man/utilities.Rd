\name{utilities}
\alias{updateList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility functions for the \code{biostat3} package
}
\description{
Utility functions for the \code{biostat3} package.
}
\usage{
updateList(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
base object (list)
}
  \item{...}{
arguments to update
}
}
\details{
Update the names in the base \code{object} list that are specified in
the arguments to update.
}
\value{
  list
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
  updateList(list(a=1,b=2), a=10, c=30)
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
