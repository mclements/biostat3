## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("foreign")
## install.packages("dplyr")


###############################################################################
## Exercise 12
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit
require(dplyr)    # for data manipulation



## @knitr loadPreprocess
melanoma_raw<- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma <- melanoma_raw %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0))


## @knitr 12.a

# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence
# whereas Stata (and most other software) use the Breslow method
cox1 <- coxph(Surv(surv_mm, death_cancer) ~ sex, method=c("breslow"), data=melanoma)
summary(cox1)


## @knitr 12.b
cox2 <- coxph(Surv(surv_mm, death_cancer) ~ sex + agegrp + stage + subsite + year8594, method=c("breslow"), data=melanoma)
summary(cox2)


## @knitr 12.c
cox3 <- coxph(Surv(surv_mm, death_cancer) ~ agegrp + agegrp * sex, method=c("breslow"), data=melanoma)
summary(cox3)


## @knitr 12.d
cox4 <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp + subsite + stage, method=c("breslow"), data=melanoma)
summary(cox4)

## Test proportional hazards assumption
cox.zph(cox4, transform="identity") #Stata appears to be using 'identity'

## Stratified Cox model; separate baseline hazard functions are fit for each strata.
cox5 <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp + subsite + strata(stage), method=c("breslow"), data=melanoma)
summary(cox5)

## Test proportional hazards assumption
cox.zph(cox5, transform="identity") #Stata appears to be using 'identity'

cox5 <- coxph(Surv(surv_mm, death_cancer) ~ sex * agegrp + year8594 + agegrp + subsite + strata(stage), method=c("breslow"), data=melanoma)
summary(cox5)
