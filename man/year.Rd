\name{year}
\alias{year}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a \code{Date} vector to a numeric vector
}
\description{
Convert a \code{Date} vector to a numeric vector (either continuous or truncated).
}
\usage{
year(date, trunc = FALSE, year.length = 365.24)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{
Date vector
}
  \item{trunc}{
logical for whether to truncate the date to a whole year or consider the
date as a double (default).
}
  \item{year.length}{
assumed length of a year
}
}
\details{
For the double calculation, we use (truncated year of Date) + (date - 1 Jan of Year)/year.length.
}
\value{
numeric vector
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

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
  c(year(as.Date("2001-07-01")),year(as.Date("2001-01-01"),trunc=TRUE))
}
% one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
