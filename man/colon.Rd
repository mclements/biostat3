\name{colon}
\alias{colon}
\docType{data}
\title{
Colon cancer dataset
}
\description{
Colon cancer patients based on a simulation of a population-based cancer registry.
}
\usage{data("colon")}
\format{
  A data frame with 15564 observations on the following 18 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{Male} \code{Female}}
    \item{\code{age}}{a numeric vector age at diagnosis in years (12-99)}
    \item{\code{stage}}{a factor with levels \code{Unknown} \code{Localised} \code{Regional} \code{Distant}}
    \item{\code{mmdx}}{a numeric vector for the month of diagnosis (1-12)}
    \item{\code{yydx}}{a numeric vector for the year of diagnosis (1975-1994)}
    \item{\code{surv_mm}}{a numeric vector observed time from diagnosis
      to end of follow-up in months (0.5-251.5)}
    \item{\code{surv_yy}}{a numeric vector observed time from diagnosis
      to end of follow-up in years (0.5-20.5)}
    \item{\code{status}}{a factor with levels \code{Alive} \code{Dead: cancer} \code{Dead: other} \code{Lost to follow-up}}
    \item{\code{subsite}}{a factor with levels \code{Coecum and ascending} \code{Transverse} \code{Descending and sigmoid} \code{Other and NOS}}
    \item{\code{year8594}}{a factor with levels \code{Diagnosed 75-84} \code{Diagnosed 85-94}}
    \item{\code{agegrp}}{a factor with levels \code{0-44} \code{45-59} \code{60-74} \code{75+}}
    \item{\code{dx}}{a Date date of diagnosis}
    \item{\code{exit}}{a Date date of exit}
    \item{\code{id}}{a numeric vector individual identifier}
    \item{\code{ydx}}{a numeric vector for continuous year of diagnosis}
    \item{\code{yexit}}{a numeric vector for continuous year of exit}
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
data(colon)
## maybe str(colon) ; plot(colon) ...
}
\keyword{datasets}
