## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercise 9 in R
## Author: Annika Tillander,  Johan Zetterqvist
###############################################################################

## ## Install needed packages only need to be done once
## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("muhaz")
## install.packages("car")

###############################################################################
## Exercise 9
###############################################################################
## @knitr loadDependecies

library(biostat3)
library(dplyr)    # for data manipulation

## @knitr loadPreprocess
## Read melanoma data, select subcohorts and create a death indicator
melanoma.l <- biostat3::melanoma %>%
    filter(stage=="Localised") %>%
    mutate(death_cancer = as.numeric(status=="Dead: cancer"))

melanoma.l2 <-    mutate(melanoma.l,
                         ## Create a death indicator (only count deaths within 120 months)
                         death_cancer = death_cancer * as.numeric( surv_mm <= 120),
                         ## Create a new time variable
                         surv_mm = pmin(surv_mm, 120))

## @knitr 9.a

summary( coxfit9a <- coxph(Surv(surv_mm, death_cancer) ~ year8594,
                           data = melanoma.l2) )

## @knitr 9.b

summary( coxfit3c <- coxph(Surv(surv_mm, death_cancer) ~ year8594,
                           data = melanoma.l) )

## @knitr 9.c

summary( coxfit9c <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex + agegrp,
                           data = melanoma.l2) )

## Test if the effect of age is significant
## Use a Wald test with the car package
library(car)
linearHypothesis(coxfit9c,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))

## @knitr 9.d

## Test if the effect of age is significant
## Use an LR test with the anova function
## The model without agegrp
summary( coxfit9d <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex,
                           data = melanoma.l2) )
## LR test
anova(coxfit9c, coxfit9d)

## @knitr 9.e

summary( coxfit9e <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex + agegrp,
                           data = melanoma.l2) )

## Split follow up by year
melanoma.spl <- survSplit(melanoma.l2, cut=12*(0:10), end="surv_mm", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl <- mutate(melanoma.spl,
                       pt = surv_mm - start,
                       fu = as.factor(start) )

## Run Poisson regression
summary(poisson9e <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))

eform(poisson9e)
summary(coxfit9e)

## @knitr 9.f

sfit9f <- survfit(Surv(surv_mm, event=death_cancer) ~ 1,
                  data = melanoma.l2 )
## Have a look at the fitted object
str(sfit9f, 1)

head(sfit9f$time)

## Split follow up by year
melanoma.spl2 <- survSplit(melanoma.l2, cut=sfit9f$time, end="surv_mm", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl2 <- mutate(melanoma.spl2,
                        pt = surv_mm - start,
                        fu = as.factor(start) )

## Run Poisson regression (slow)
if (redo <- FALSE) {
  poisson9f <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ),
                    family = poisson,
                    data = melanoma.spl2 )
  save(poisson9f,file="poisson9f.RData")
} else {
  load("poisson9f.RData")
}

## IRR
coef9f <- eform(poisson9f)
## We are not interested in nuisance parameters fu1, fu2, etc
npar <- nrow(coef9f)
pars <- (npar-4):npar
coef9f[pars,]

summary(coxfit9e)
