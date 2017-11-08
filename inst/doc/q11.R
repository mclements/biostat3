## Date: 2015-03-04
## Purpose: To do the solution for Biostat III exercise 11 in R
## Author: Johan Zetterqvist
## Revised: Mark Clements 2017-11-03
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("dplyr")
## install.packages("foreign")

###############################################################################
## Exercise 11
###############################################################################
## @knitr loadDependencies

library(biostat3) # 
library(dplyr)    # for data manipulation

## @knitr loadPreprocess

## Read melanoma data
## and select subcohorts
data(melanoma)
melanoma.l <- melanoma %>%
  filter(stage=="Localised") %>%
  mutate(
    ## Create a death indicator
    death_cancer = as.numeric(status=="Dead: cancer"),
    death_any = as.numeric(status=="Dead: cancer" | status=="Dead: other") )

## Truncate follow-up time

melanoma.l2 <-
  mutate(melanoma.l,
         ## Create new death indicators (only count deaths within 120 months)
         death_cancer = death_cancer * as.numeric( surv_mm <= 120),
         death_any = death_any * as.numeric( surv_mm <= 120),
         ## Create a new time variable
         surv_mm = pmin(surv_mm, 120))

## @knitr 11.a

summary( coxfit11a <- coxph(Surv(surv_mm, death_any) ~ sex + year8594 + agegrp,
                           data = melanoma.l2) )

## @knitr 11.b

summary( coxfit11b <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp,
                           data = melanoma.l2) )
