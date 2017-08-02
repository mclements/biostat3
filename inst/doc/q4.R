## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28
###############################################################################

###############################################################################
## Exercise 4
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit
require(KMsurv)
require(dplyr)    # for data manipulation

## @knitr loadPreprocess
melanoma_raw<- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma <- melanoma_raw %>%
    filter(stage=="Localised") %>%
    mutate(year = floor(surv_yy),
           month = floor(surv_mm),
           death_cancer = ifelse( status == "Dead: cancer", 1, 0))

## @knitr actuarialYears
melanomaByYear <- melanoma %>%
    group_by(year) %>%
    summarise(nevent = sum(death_cancer), nlost = length(death_cancer)-sum(death_cancer))
with(melanomaByYear, lifetab(c(year,tail(year,1)+1), nrow(melanoma), nlost, nevent))[,1:7]

## @knitr actuarialMonths
melanomaByMonth <- melanoma %>%
    group_by(month) %>%
    summarise(nevent = sum(death_cancer), nlost = length(death_cancer)-sum(death_cancer))
with(melanomaByMonth, lifetab(c(month,tail(month,1)+1), nrow(melanoma), nlost, nevent))[110:130,1:7]

## @knitr kmYears
mfit_years <- survfit(Surv(year, death_cancer) ~ 1, data = melanoma)
summary(mfit_years)

## @knitr kmMonths
mfit_months <- survfit(Surv(month, death_cancer) ~ 1, data = melanoma)
data.frame(summary(mfit_months)[c(2:4,6,8,10,9)])[110:130,]
