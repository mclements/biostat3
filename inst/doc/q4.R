## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28
## Edited: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 4
###############################################################################
## @knitr loadDependencies
library(biostat3) 
library(dplyr)    # for data manipulation

## @knitr loadPreprocess
melanoma <- biostat3::melanoma %>%
    filter(stage=="Localised") %>%
    mutate(year = floor(surv_yy),
           month = floor(surv_mm),
           death_cancer = ifelse( status == "Dead: cancer", 1, 0))

## @knitr actuarialYears
lifetab2(Surv(floor(surv_yy),death_cancer)~1, data = melanoma)[,1:7]

## @knitr actuarialMonths
lifetab2(Surv(floor(surv_mm),death_cancer)~1, data = melanoma)[110:130,1:7]

## @knitr kmYears
mfit_years <- survfit(Surv(year, death_cancer) ~ 1, data = melanoma)
summary(mfit_years)

## @knitr kmMonths
mfit_months <- survfit(Surv(month, death_cancer) ~ 1, data = melanoma)
summary(mfit_months,times=110:130)
