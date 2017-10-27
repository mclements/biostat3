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
  A data frame with 7775 observations on the following 15 variables.
  \describe{
    \item{\code{sex}}{a factor for the sex of the cancer patient with
      levels \code{Male} and \code{Female}}
    \item{\code{age}}{a numeric vector for the age at cancer diagnosis}
    \item{\code{stage}}{a factor for the cancer stage (or extent) as
      cancer diagnosis, with levels \code{Unknown}, \code{Localised},
      \code{Regional} and \code{Distant}}
    \item{\code{mmdx}}{a numeric vector for the month of cancer diagnosis}
    \item{\code{yydx}}{a numeric vector for the (truncated) year of cancer diagnosis}
    \item{\code{surv_mm}}{a numeric vector for the number of months
      since cancer diagnosis}
    \item{\code{surv_yy}}{a numeric vector for the number of years since
    cancer diagnosis}
    \item{\code{status}}{a factor for status at end of follow-up, with
      levels \code{Alive}, \code{Dead: cancer}, \code{Dead: other} and \code{Lost to follow-up}}
    \item{\code{subsite}}{a factor with levels \code{Head and Neck},
      \code{Trunk}, \code{Limbs} and \code{Multiple and NOS}}
    \item{\code{year8594}}{a factor for calendar period of cancer
      diagnosis with levels \code{Diagnosed 75-84} \code{Diagnosed 85-94}}
    \item{\code{dx}}{a Date for the date of cancer diagnosis}
    \item{\code{exit}}{a Date for the date of study exit}
    \item{\code{agegrp}}{a factor for age group at cancer diagnosis with
      levels \code{0-44}, \code{45-59}, \code{60-74} and \code{75+}}
    \item{\code{id}}{a numeric vector for row ID}
    \item{\code{ydx}}{a numeric for the (continuous) year of cancer diagnosis}
    \item{\code{yexit}}{a numeric for the (continuous) year of study exit}
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
%% \examples{
%% data(melanoma)
%% ## maybe str(melanoma) ; plot(melanoma) ...
%% }
\keyword{datasets}
