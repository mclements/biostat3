## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("survival")
## install.packages("dplyr")


###############################################################################
## Exercise 22
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit
require(dplyr)    # for data manipulation

###########################################
### A help function to calculate ###
### and print incidence (hazard) ratios
### from a fitted poisson regression
### from glm
###########################################

IRR <- function(fit){
    summfit <- summary(fit )$coefficients
    IRfit <- exp( cbind( summfit[, 1:2], summfit[, 1] - 1.96*summfit[, 2], summfit[, 1] +
                        1.96*summfit[, 2] ) )
    colnames(IRfit) <- c("IRR", "Std. err", "CI_lower", "CI_upper")
    print(IRfit)
}


## @knitr loadPreprocess
brv_raw <- read.dta("http://biostat3.net/download/brv.dta")


## @knitr 22.a
head(brv_raw)

##Look at the five first couples
brv_raw %>% select(couple, id, sex, doe, dosp, dox, fail) %>% filter(couple<=5) %>% arrange(couple, id)


## @knitr 22.b
brv <- brv_raw %>%
    mutate(age_entry = as.numeric(doe - dob) / 365.24, # Calc age at entry
           att_age = as.numeric(dox - dob) / 365.24,   # Calc attained age
           t_at_risk = att_age - age_entry)            # Calc time at risk

## crude rates
brv %>%
    select(fail, t_at_risk, sex) %>%
    group_by(sex) %>%
    summarise(D = sum(fail), Y = sum(t_at_risk)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))


poisson22b <- glm( fail ~ sex + offset( log( t_at_risk) ), family=poisson, data=brv)
summary( poisson22b )
IRR(poisson22b)

brv %>%
    select(sex, age_entry) %>%
    group_by(sex) %>%
    summarise(meanAgeAtEntry = mean(age_entry))

## @knitr 22.c

## Creating times relativ to spouse death (year=0)
brv2 <- mutate(brv_raw,
               y_before_sp_dth =  as.numeric(doe -dosp) / 365.24,
               y_after_sp_dth = as.numeric(dox - dosp) / 365.24)

## Splitting at spouse death (year=0)
brvSplit <- survSplit(brv2, cut = 0, end="y_after_sp_dth", start="y_before_sp_dth", id="id",event="fail")

## Calculating risk times
brvSplit <- mutate(brvSplit,
                   t_sp_at_risk =   y_after_sp_dth - y_before_sp_dth,
                   brv = ifelse(y_after_sp_dth > 0, 1, 0))

## Look at the five first couples
brvSplit %>% select(couple, id, sex, doe, dosp, dox, fail, y_before_sp_dth, y_after_sp_dth, t_sp_at_risk) %>% filter(couple<=5) %>% arrange(couple, id)

## @knitr 22.d
poisson22d <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=brvSplit)
summary(poisson22d)
IRR(poisson22d)

## @knitr 22.e
## Poisson regression for sex==1
poisson22e1 <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=filter(brvSplit, sex==1))
summary(poisson22e1)
IRR(poisson22e1)

## Poisson regression for sex==2
poisson22e2 <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=filter(brvSplit, sex==2))
summary(poisson22e2)
IRR(poisson22e2)

## Poisson regression, interaction with sex
brvSplit2 <- mutate(brvSplit,
                    sex = as.factor(sex),
                    brv = as.factor(brv))
poisson22e3 <- glm( fail ~ sex + brv:sex + offset( log(t_sp_at_risk) ), family=poisson, data=brvSplit2)
summary(poisson22e3)
IRR(poisson22e3)

## @knitr 22.f
## Translate time scale from years from spouse death to ages
brvSplit3 <- brvSplit2 %>%
    mutate(age_sp_dth =  as.numeric(dosp - dob) / 365.24, # Age at spouse death
           age_start = age_sp_dth + y_before_sp_dth,      # Age at start of timeband
           age_end = age_sp_dth + y_after_sp_dth)         # Age at end of timeband

age_cat <- seq(70,100,5) # Split at these ages
brvSplit4 <- survSplit(brvSplit3, cut=age_cat, start="age_start", end="age_end", event="fail", zero = 0)

brvSplit4 <- mutate(brvSplit4,
                    t_at_risk = age_end- age_start, # Creating new time at risk
                    age = cut(age_end, age_cat))   # Creating age band category

## Calculate crude rates
brvSplit4 %>%
    select(fail, age, t_at_risk) %>%
    group_by(age) %>%
    summarise(D = sum(fail), Y = sum(t_at_risk), Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

poisson22f1 <- glm( fail ~ brv + age + offset( log(t_at_risk) ), family=poisson, data=brvSplit4)
summary(poisson22f1)
IRR(poisson22f1)

poisson22f2 <- glm( fail ~ brv +age + sex + offset( log(t_at_risk) ), family=poisson, data=brvSplit4)
summary(poisson22f2)
IRR(poisson22f2)

## @knitr 22.g
poisson22g <- glm( fail ~ age + sex + brv:sex + offset( log(t_at_risk) ), family=poisson, data=brvSplit4)
summary(poisson22g)
IRR(poisson22g)

## @knitr 22.i
summary(coxph(Surv(age_start, age_end, fail) ~ brv,
              method=c("breslow"), data = brvSplit4))


summary(coxph(Surv(age_start, age_end, fail) ~ brv + sex,
              method=c("breslow"), data = brvSplit4))

## @knitr 22.j
summary(coxph(Surv(age_start, age_end, fail) ~ sex + sex:brv,
              method=c("breslow"), data = brvSplit4))
