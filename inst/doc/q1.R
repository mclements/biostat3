## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-17, 2016-03-01
## Edited: Mark Clements, 2017-08-02
###############################################################################

###############################################################################
## Exercise 1b
###############################################################################

## @knitr loadDependecies
library(biostat3)

## @knitr lifeTable
print(lifetab2(Surv(floor(surv_yy), status == "Dead: cancer")~1, colon, breaks=0:10), digits=2)

## @knitr KaplanMeier
mfit <- survfit(Surv(surv_mm, status == "Dead: cancer") ~ 1, data = colon_sample) # make Kaplan-Meier estimates
summary(mfit)                                                  # print Kaplan-Meier table
plot(mfit,                                                     # plot Kaplan-Meier curve
     ylab="S(t)",
     xlab="Time since diagnosis in months",
     main = "Kaplan−Meier estimates of cause−specific survival")
