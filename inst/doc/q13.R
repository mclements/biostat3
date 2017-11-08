## Purpose: To do the solution for the Biostat III exercise 13 in R
## Author: Johan Zetterqvist, 2015-03-04
## Revised: Mark Clements, 2017-11-03
###############################################################################

## ## Install needed packages only need to be done once
## Install needed packages only need to be done once
## install.packages("foreign")

###############################################################################
## Exercise 13
###############################################################################
## @knitr loadDependencies
library(biostat3) # for Surv and survfit

## @knitr loadPreprocess


## @knitr 13.a
## y is the observed time
## so it already measures time since entry
poisson13a <- glm( chd ~ hieng + offset( log( y ) ), family=poisson, data=diet)
summary(poisson13a)
eform(poisson13a)

cox13a <- coxph(Surv(y, chd) ~ hieng, data=diet)
summary(cox13a)

## @knitr 13.b

## counting process data with coxph
scale <- 365.24
cox13b <- coxph(Surv((doe-dob)/scale, (dox-dob)/scale, chd) ~ hieng, data=diet)
summary(cox13b)
