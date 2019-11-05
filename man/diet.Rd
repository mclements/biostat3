\name{diet}
\alias{diet}
\docType{data}
\title{
Diet data set
}
\description{
Diet data set
}
\usage{data("diet")}
\format{
  A data frame with 337 observations on the following 15 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{chd}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{hieng}}{a factor with levels \code{low} \code{high}}
    \item{\code{energy}}{a numeric vector}
    \item{\code{job}}{a factor with levels \code{driver} \code{conductor} \code{bank}}
    \item{\code{month}}{a numeric vector}
    \item{\code{height}}{a numeric vector}
    \item{\code{weight}}{a numeric vector}
    \item{\code{doe}}{a Date for date of study entry}
    \item{\code{dox}}{a Date for date of study exit}
    \item{\code{dob}}{a Date for date of birth}
    \item{\code{yob}}{a numeric vector for continuous year of birth}
    \item{\code{yoe}}{a numeric vector for continuous year of entry}
    \item{\code{yox}}{a numeric vector for continuous year of exit}
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the __description__ above ~~
%% }
%% \source{
%% %%  ~~ reference to a publication or URL from which the data were obtained ~~
%% }
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(diet)
## maybe str(diet) ; plot(diet) ...
}
\keyword{datasets}
