\name{colon_sample}
\alias{colon_sample}
\docType{data}
\title{
Sample from the \code{\link{colon}} dataset used for teaching. 
}
%% \description{
%% %%  ~~ A concise (1-5 lines) description of the dataset. ~~
%% }
\usage{data("colon_sample")}
\format{
  A data frame with 35 observations on the following 9 variables.
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
%% data(colon_sample)
%% ## maybe str(colon_sample) ; plot(colon_sample) ...
%% }
\keyword{datasets}
