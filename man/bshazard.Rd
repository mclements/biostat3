\name{as.data.frame.bshazard}
\alias{as.data.frame.bshazard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to work with bshazard objects.
}
\description{
Convert a bshazard object to a data-frame.
}
\usage{
\method{as.data.frame}{bshazard}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
bshazard object
}
  \item{\dots}{
other arguments
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
Returns a data-frame with names time, hazard, conf.low and conf.high
(cf. lower.ci and upper.ci provided in the object).
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, ...) 
{
    with(x, data.frame(time, hazard, conf.low = lower.ci, conf.high = upper.ci))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
