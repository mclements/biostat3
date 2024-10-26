\name{muhaz2}
\alias{muhaz2}
\alias{plot.muhaz2}
\alias{plot.muhazList}
\alias{lines.muhazList}
\alias{lines.muhaz2}
\alias{ggplot.muhazList}
\alias{summary.muhazList}
\alias{as.data.frame.muhaz}
\alias{as.data.frame.muhazList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula wrapper for the \code{\link[muhaz]{muhaz}} function from the
\code{muhaz} package.
}
\description{
Formula wrapper for the \code{\link[muhaz]{muhaz}} function from the
\code{muhaz} package.
}
\usage{
muhaz2(formula, data, subset, max.time, ...)
\method{plot}{muhaz2}(x, haz.scale=1, ylab="Hazard", ylim=NULL, log="", ...)
\method{plot}{muhazList}(x, lty=1:5, col=1:length(x), log="", legend.args=list(), ...)
\method{lines}{muhaz2}(x, ..., haz.scale = 1)
\method{lines}{muhazList}(x, lty=1, col=1:length(x), ...)
\method{summary}{muhazList}(object, ...)
\method{ggplot}{muhazList}(data, mapping=NULL, 
                            xlab="Time", ylab="Hazard", ...,
                            environment = parent.frame())
\method{as.data.frame}{muhazList}(x, row.names, optional, ...)
\method{as.data.frame}{muhaz}(x, row.names, optional, ...)
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
    subset predictate for the dataset
}
  \item{max.time}{
    maximum follow-up time for the hazards
}
  \item{xlab}{
    graphics argument for xlab (x-axis label)
}
  \item{ylab}{
    graphics argument for ylab (y-axis label)
}
  \item{lty}{
    graphics argument for line type
}
  \item{col}{
    graphics argument for line colour
}
  \item{legend.args}{
    a list of options that are passed to the legend call. Defaults are \code{list(x="topright",legend=names(x),col=col,lty=lty)}.
}
  \item{haz.scale}{
    scale for the hazard in the plot
}
  \item{row.names}{
    not currently used
}
  \item{object}{
    muhazList object
}
  \item{ylim}{
    graphics argument for the limits of the y axis
}
  \item{log}{
    graphics argument for a log transformation of the x or y axes
}
  \item{x}{
    muhazList or muhaz object
}
\item{environment}{*[Deprecated]* Used prior to tidy evaluation.}
  \item{optional}{
    not currently used
  }

  \item{mapping}{
  Default list of aesthetic mappings to use for plot. If not
          specified, must be supplied in each layer added to the plot.
  }
  \item{...}{
    other arguments
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
For a single strata, this is a \code{\link[muhaz]{muhaz}} object. For multiple strata, this is a \code{muhazList} object, which includes methods for 
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
