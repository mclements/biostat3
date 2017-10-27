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
  A data frame with 15564 observations on the following 13 variables.
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
    \item{\code{subsite}}{a factor for the sub-site of the cancer, with
      levels \code{Coecum and ascending}, \code{Transverse},
      \code{Descending and sigmoid}, and \code{Other and NOS}}
    \item{\code{year8594}}{a factor for calendar period of cancer
      diagnosis with levels \code{Diagnosed 75-84} \code{Diagnosed 85-94}}
    \item{\code{agegrp}}{a factor for age group at cancer diagnosis with levels \code{0-44} \code{45-59} \code{60-74} \code{75+}}
    \item{\code{dx}}{a Date for the date of cancer diagnosis}
    \item{\code{exit}}{a Date for the date of study exit}
    \item{\code{id}}{a numeric for a row ID}
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
%% data(colon)
%% ## maybe str(colon) ; plot(colon) ...
%% }
\keyword{datasets}
