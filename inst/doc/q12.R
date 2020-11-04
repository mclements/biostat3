## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
## Revised: Mark Clements, 2017-11-03
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")


###############################################################################
## Exercise 12
###############################################################################
## @knitr loadDependencies
library(biostat3) 


## @knitr loadPreprocess
melanoma.l <- transform(biostat3::melanoma,
                      death_cancer = ifelse( status == "Dead: cancer", 1, 0))


## @knitr 12.a

# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence of ties,
# whereas Stata (and some other software) uses the Breslow method
cox1 <- coxph(Surv(surv_mm, death_cancer) ~ sex, data=melanoma.l)
summary(cox1)


## @knitr 12.b
cox2 <- coxph(Surv(surv_mm, death_cancer) ~ sex + agegrp + stage + subsite + year8594, data=melanoma.l)
summary(cox2)


## @knitr 12.c
cox3 <- coxph(Surv(surv_mm, death_cancer) ~ agegrp + agegrp * sex, data=melanoma.l)
summary(cox3)


## @knitr 12.d.i
cox4 <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp + subsite + stage, data=melanoma.l)
summary(cox4)

## Test proportional hazards assumption
cox.zph(cox4, transform="log") 
cox.zph(cox4, transform="identity") # Stata default

## @knitr 12.d.ii
## Stratified Cox model; separate baseline hazard functions are fit for each strata.
cox5 <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp + subsite + strata(stage), data=melanoma.l)
summary(cox5)

## Test proportional hazards assumption
cox.zph(cox5, transform="log") 
cox.zph(cox5, transform="identity") 

## @knitr 12.d.iii
cox5 <- coxph(Surv(surv_mm, death_cancer) ~ sex * agegrp + year8594 + agegrp + subsite + strata(stage), data=melanoma.l)
summary(cox5)
anova(cox5)
