\name{biostat3-package}
\alias{biostat3-package}
\alias{biostat3}
\docType{package}
\title{
\packageTitle{biostat3}
}
\description{
\packageDescription{biostat3}
}
%% \details{

%% The DESCRIPTION file:
%% \packageDESCRIPTION{biostat3}
%% \packageIndices{biostat3}
%% ~~ An overview of how to use the package, including the most important ~~
%% ~~ functions ~~
%% }
\author{
\packageAuthor{biostat3}

Maintainer: \packageMaintainer{biostat3}
}
%% \references{
%% ~~ Literature or other references for background information ~~
%% }
\keyword{ package }
\keyword{ survival }
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
  plot(muhaz2(Surv(surv_mm, status == "Dead: cancer")~1, melanoma))
}
