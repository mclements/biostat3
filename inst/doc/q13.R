## Purpose: To do the solution for the Biostat III exercise 13 in R
## Author: Johan Zetterqvist, 2015-03-04
###############################################################################

## ## Install needed packages only need to be done once
## Install needed packages only need to be done once
## install.packages("foreign")

###############################################################################
## Exercise 13
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit

## @knitr loadPreprocess
diet <- data.frame(read.dta("http://biostat3.net/download/diet.dta"))

## @knitr 13.a
## y is the observed time
## so it already measures time since entry
poisson13a <- glm( chd ~ hieng + offset( log( y ) ), family=poisson, data=diet)
summary(poisson13a)
exp(cbind(coef(poisson13a), confint(poisson13a)))

cox13a <- coxph(Surv(y, chd) ~ hieng, data=diet, ties="breslow")
summary(cox13a)

## @knitr 13.b
poisson13a <- glm( chd ~ hieng + offset( log( y ) ), family=poisson, data=diet)
summary(poisson13a)
exp(cbind(coef(poisson13a), confint(poisson13a)))

## Create two new variables for age
## age at study entry
diet$entry_age <-  as.numeric(diet$doe - diet$dob) / 365.24
## age at study exit
diet$exit_age <- as.numeric(diet$dox - diet$dob) / 365.24

## Use the new age variables to provide
## counting process data to coxph
cox13b <- coxph(Surv(entry_age, exit_age, chd) ~ hieng, data=diet, ties="breslow")
summary(cox13b)
