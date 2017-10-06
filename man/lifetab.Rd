\name{lifetab}
\alias{lifetab}
\title{ Create cohort life table }
\description{
  Create cohort life table.
}
\usage{
lifetab(tis, ninit, nlost, nevent)
}
\arguments{
  \item{tis}{ a vector of end points of time intervals, whose length is
    1 greater than nlost and nevent. }
  \item{ninit}{ the number of subjects initially entering the study. }
  \item{nlost}{ a vector of the number of individuals lost follow or
    withdrawn alive for whatever reason. }
  \item{nevent}{ a vector of the number of individuals who experienced
    the event }
}
\value{
  A data.frame with the following columns:
  \item{nsubs}{ the number of subject entering the intervals
    who have not experienced the event.}
  \item{nlost}{ the number of individuals lost follow or
    withdrawn alive for whatever reason. }
  \item{nrisk}{ the estimated number of individuals at risk of
    experiencing the event. }
  \item{nevent}{ the number of individuals who experienced
    the event. }
  \item{surv}{ the estimated survival function at the start
    of the intervals. }
  \item{pdf}{ the estimated probability density function at
    the midpoint of the intervals. }
  \item{hazard}{ the estimated hazard rate at the midpoint
    of the intervals. }
  \item{se.surv}{ the estimated standard deviation of survival at the
    beginning of the intervals. }
  \item{se.pdf}{ the estimated standard deviation of the prbability
    density function at the midpoint of the intervals. }
  \item{se.hazard}{ the estimated standard deviation of the hazard
    function at the midpoint of the intervals}
  The row.names are the intervals.
}
\author{ Jun Yan \email{jyan@stat.uconn.edu} }
\examples{
tis <- c(0, 2, 3, 5, 7, 11, 17, 25, 37, 53, NA)
nsubs <- c(927, 848, 774, 649, 565, 449, 296, 186, 112, 27)
nlost <- c(2, 3, 6, 9, 7, 5, 3, rep(0, 3))
nevent <- c(77, 71, 119, 75, 109, 148, 107, 74, 85, 27)

lifetab(tis, nsubs[1], nlost, nevent)
}
\keyword{ manip }
