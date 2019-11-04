\name{colon}
\alias{colon}
\docType{data}
\title{
Colon cancer dataset
}
\description{
Colon cancer dataset
}
\usage{data("colon")}
\format{
  A data frame with 15564 observations on the following 14 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{Male} \code{Female}}
    \item{\code{age}}{a numeric vector}
    \item{\code{stage}}{a factor with levels \code{Unknown} \code{Localised} \code{Regional} \code{Distant}}
    \item{\code{mmdx}}{a numeric vector}
    \item{\code{yydx}}{a numeric vector}
    \item{\code{surv_mm}}{a numeric vector}
    \item{\code{surv_yy}}{a numeric vector}
    \item{\code{status}}{a factor with levels \code{Alive} \code{Dead: cancer} \code{Dead: other} \code{Lost to follow-up}}
    \item{\code{subsite}}{a factor with levels \code{Coecum and ascending} \code{Transverse} \code{Descending and sigmoid} \code{Other and NOS}}
    \item{\code{year8594}}{a factor with levels \code{Diagnosed 75-84} \code{Diagnosed 85-94}}
    \item{\code{agegrp}}{a factor with levels \code{0-44} \code{45-59} \code{60-74} \code{75+}}
    \item{\code{dx}}{a Date}
    \item{\code{exit}}{a Date}
    \item{\code{id}}{a numeric vector}
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
