\name{smoothHaz}
\alias{smoothHaz}
\alias{plot.smoothHaz}
\title{
Simple implementation for kernel density smoothing of the Nelson-Aalen
estimator. 
}
\description{
Simple implementation for kernel density smoothing of the Nelson-Aalen
estimator. Prefer \code{muhaz} for right censored data and
\code{bshazard} for left truncated and right censored data.
}
\usage{
smoothHaz(object, n.grid = 300, kernel = "epanechnikov", from = NULL, to = NULL, min.n.risk = 1, ...)
\method{plot}{smoothHaz}(x, xlab = "Time", ylab = "Hazard", type = "l", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    \code{survfit} object
}
  \item{n.grid}{
    number of grid points; passed to \code{density}
}
  \item{kernel}{
    kernel used; passed to \code{density}
}
  \item{from}{
    left boundary; passed to \code{density}
}
  \item{to}{
    right boundary; passed to \code{density}
}
  \item{min.n.risk}{
    minimum number at risk
}
  \item{x}{
    object of class \code{smoothHaz}
}
  \item{xlab}{
    graphics argument
}
  \item{ylab}{
    graphics argument
}
  \item{type}{
    graphics argument
}
  \item{\dots}{
    Other arguments
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
%% \value{
%% %%  ~Describe the value returned
%% %%  If it is a LIST, use
%% %%  \item{comp1 }{Description of 'comp1'}
%% %%  \item{comp2 }{Description of 'comp2'}
%% %% ...
%% }
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
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% function (object, n.grid = 300, kernel = "epanechnikov", from = NULL, 
%%     to = NULL, min.n.risk = 1, ...) 
%% {
%%     x <- object
%%     index <- x$n.risk >= min.n.risk
%%     time <- x$time[index]
%%     weights <- x$n.event[index]/x$n.risk[index]
%%     if (is.null(from)) 
%%         from <- min(x$time)
%%     if (is.null(to)) 
%%         to <- max(x$time)
%%     newobject <- suppressWarnings(stats::density(time, weight = weights, 
%%         kernel = kernel, n = n.grid, from = from, to = to, ...))
%%     structure(newobject, n.grid = n.grid, call = match.call(), 
%%         class = c("smoothHaz", "density"))
%%   }
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")

