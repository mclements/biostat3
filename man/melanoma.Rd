\name{melanoma}
\alias{melanoma}
\docType{data}
\title{
Melanoma cancer dataset
}
\description{
Melanoma cancer dataset
}
\usage{data("melanoma")}
\format{
  A data frame with 7775 observations on the following 14 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{Male} \code{Female}}
    \item{\code{age}}{a numeric vector}
    \item{\code{stage}}{a factor with levels \code{Unknown} \code{Localised} \code{Regional} \code{Distant}}
    \item{\code{mmdx}}{a numeric vector}
    \item{\code{yydx}}{a numeric vector}
    \item{\code{surv_mm}}{a numeric vector}
    \item{\code{surv_yy}}{a numeric vector}
    \item{\code{status}}{a factor with levels \code{Alive} \code{Dead: cancer} \code{Dead: other} \code{Lost to follow-up}}
    \item{\code{subsite}}{a factor with levels \code{Head and Neck} \code{Trunk} \code{Limbs} \code{Multiple and NOS}}
    \item{\code{year8594}}{a factor with levels \code{Diagnosed 75-84} \code{Diagnosed 85-94}}
    \item{\code{dx}}{a Date}
    \item{\code{exit}}{a Date}
    \item{\code{agegrp}}{a factor with levels \code{0-44} \code{45-59} \code{60-74} \code{75+}}
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
data(melanoma)
## maybe str(melanoma) ; plot(melanoma) ...
}
\keyword{datasets}
