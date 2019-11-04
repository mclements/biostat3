\name{brv}
\alias{brv}
\docType{data}
\title{
Bereavement dataset
}
\description{
Bereavement dataset
}
\usage{data("brv")}
\format{
  A data frame with 399 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{a numeric vector the id of a subject}
    \item{\code{couple}}{a numeric vector for the id of a couple}
    \item{\code{dob}}{a Date for the date of birth}
    \item{\code{doe}}{a Date for the date of entry into study}
    \item{\code{dox}}{a Date for the date of exit from study}
    \item{\code{dosp}}{a Date for the date of bereavement}
    \item{\code{fail}}{a numeric vector for status at study exit 0=alive 1=died}
    \item{\code{group}}{a numeric vector for Group}
    \item{\code{disab}}{a numeric vector for disability level}
    \item{\code{health}}{a numeric vector for perceived health status}
    \item{\code{sex}}{a numeric vector for sex 1=M 2=F}
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
data(brv)
## maybe str(brv) ; plot(brv) ...
}
\keyword{datasets}
