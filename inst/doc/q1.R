## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-17, 2016-03-01
## Edited: Mark Clements, 2017-08-02
###############################################################################

###############################################################################
## Exercise 1b
###############################################################################

## @knitr loadDependencies
library(biostat3)
library(survminer)
library(knitr)

## @knitr printData
kable(biostat3::colon_sample, "html")

## @knitr lifeTable
lifetab2(Surv(floor(surv_yy), status == "Dead: cancer")~1,
         colon_sample, breaks=0:10) |>
    kable("html", digits=2)

## @knitr KaplanMeier
mfit <- survfit(Surv(surv_mm/12, status == "Dead: cancer") ~ 1,
                data = colon_sample) # make Kaplan-Meier estimates
summary(mfit) # print Kaplan-Meier table
plot(mfit,    # plot Kaplan-Meier curve
     ylab="S(t)",
     xlab="Time since diagnosis (years)",
     main = "Kaplan−Meier estimates of cause−specific survival")

ggsurvplot(mfit, # plot Kaplan-Meier curve
     ylab="S(t)",
     xlab="Time since diagnosis (years)",
     main = "Kaplan−Meier estimates of cause−specific survival",
     risk.table = TRUE,
     conf.int = TRUE,
     ggtheme = theme_minimal())
