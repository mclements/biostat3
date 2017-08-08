## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist
## Modified: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 7
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("muhaz")
## install.packages("car")

## @knitr loadDependecies
library(biostat3)
library(dplyr)    # for data manipulation
library(car)      # for linearHypothesis and deltaMethod


## @knitr loadPreprocess
## Read melanoma data

## Create a new dataset with only localised cancer
melanoma.l <- filter(biostat3::melanoma, stage=="Localised")
head( melanoma.l )
summary(melanoma.l)

## @knitr 7.a.i

## Plot Kaplan-Meier curve using survfit
## Create a new event indicator
melanoma.l <- mutate(melanoma.l,
                   death_cancer = as.numeric(status=="Dead: cancer") )

## Create a fitted object for our subcohort
## using survfit
sfit7a1 <- survfit(Surv(surv_mm, event=death_cancer) ~ year8594,
                   data = melanoma.l )

## Have a look at the fitted object
str(sfit7a1, 1)

## Plot the survival curve (with some bells and whistles)
plot(sfit7a1,
     ## No automatic labelling of the curve (we do that ourselves)
     mark.time=F,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ## Make the plot prettier
     xlab="Years since diagnose",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
## Add legend too
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),lty=c("solid","dashed"))

### TRY IF YOU WANT ###
## install.packages("survMisc")
## require(survMisc)
## autoplot(sfit7a1)


## @knitr 7.a.ii

## To plot smoothed hazards, we use the muhaz package (using the muhaz2 wrapper)

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"), lty=1:2)
legend("topright",legend=levels(melanoma.l$year8594), col=c("blue","red"), lty=1:2)

## @knitr 7.a.iii
## Compare with Kaplan-Meier plot
par(mfrow=c(1,2)) ## Two graphs in the same window

plot(sfit7a1,
     ## No automatic labelling of the curve (we do that ourselves)
     mark.time=F,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ## Make the plot prettier
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),lty=c("solid","dashed"))

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"),lty=c("solid","dashed"))
legend("topright",legend=levels(melanoma.l$year8594), col=c("blue","red"),lty=c("solid","dashed"))

## @knitr 7.b

## Calculate the incidence rate by time of diagnosis
rates_by_diag_yr <- melanoma.l %>%

    ## Stratify on year of diagnosis
    group_by(year8594) %>%

    ## We only need the deaths and the persontime
    select(death_cancer, surv_mm) %>%

    summarise(
        ## Calculate the number of deaths
        D = sum(death_cancer),
        ## Calculate the total person time
        Y = sum(surv_mm)/12/1000,
        ## Calculate the rate
        Rate = D/Y,
        ## Add confidence intervals for the rates
        CI_low = poisson.ci(D,Y)[1],
        CI_high = poisson.ci(D,Y)[2]) %>% print
        


## @knitr 7.c.i

## Calculate the incidence rate by time of diagnosis
## but with new variables
melanoma.l2 <- mutate(melanoma.l,
                       ## Update the death indicator (only count deaths within 120 months)
                       death_cancer = death_cancer * as.numeric(surv_mm<=120),
                       ## Create a new time variable
                       surv_mm = pmin(surv_mm, 120) )

## Calculate the rates on the truncated data
rates_by_diag_yr2 <- melanoma.l2 %>%

    ## Stratify on year of diagnosis
    group_by(year8594) %>%

    ## We only need the deaths and the persontime
    select(death_cancer, surv_mm) %>%

    summarise(
        ## Calculate the number of deaths
        D = sum(death_cancer),
        ## Calculate the total person time
        Y = sum(surv_mm)/12/1000,
        ## Calculate the rate
        Rate = D/Y,
        ## Add confidence intervals for the rates
        CI_low = poisson.ci(D,Y)[1],
        CI_high = poisson.ci(D,Y)[2]) %>% print

## @knitr 7.c.ii

## MRR full data
rates_by_diag_yr2[2, "Rate"] / rates_by_diag_yr2[1, "Rate"]
## ## MRR truncated data
## rates_by_diag_yr2[2, "Rate"] / rates_by_diag_yr2[1, "Rate"]

## @knitr 7.c.iii
## Use glm to estimate the rate ratios
## we scale the offset term to personyears
poisson7c <- glm( death_cancer ~ year8594 + offset( log( surv_mm/12/1000 ) ), family=poisson, data=melanoma.l2 )
summary( poisson7c )

## IRR
eform(poisson7c)


## Note that the scaling of the offset term only has an impact on the intercept
summary( glm( death_cancer ~ year8594 + offset( log( surv_mm ) ),
             family=poisson, data=melanoma.l2 ) )

## @knitr 7.d

## Add a new variable for year
melanoma.l2 <- mutate( melanoma.l2, surv_yy1 = surv_mm/12)

## Split follow up by year
melanoma.spl <- survSplit(melanoma.l2, cut=0:9, end="surv_yy1", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl <- mutate(melanoma.spl,
                       pt = surv_yy1 - start,
                       fu = as.factor(start) )

## @knitr 7.e

## Calculate the incidence rate by observation year
yearly_rates <- melanoma.spl %>%

    ## Stratify on follow-up year
    group_by(fu) %>%

    ## We only need the deaths and the persontime
    select(death_cancer, pt) %>%

    summarise(
        ## Calculate the number of deaths
        D = sum(death_cancer),
        ## Calculate the total person time measured in 1000
        Y = sum(pt)/1000,
        ## Calculate the rate
        Rate = D/Y,
        CI_low = poisson.ci(D,Y)[1],
        CI_high = poisson.ci(D,Y)[2]) %>% print


## Plot by year
with(yearly_rates, matplot(fu,
                           cbind(Rate, CI_low,
                             CI_high),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by years since diagnosis",
                           ylab="Incidence rate per 1000 person-years",
                           xlab="Years since diagnosis") )

## @knitr 7.f
# Plot smoothed hazards

par(mfrow=c(1,2))
with(yearly_rates, matplot(as.numeric(as.character(fu))+0.5,
                           cbind(Rate, CI_low,
                             CI_high),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by time since diagnosis",
                           ylab="Mortality rate per 1000 person-years",
                           xlab="Years since diagnosis") )

hazfit7f <- muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ 1, data = melanoma.l)
## scale hazard by 1000
plot(hazfit7f, xlab="Years since diagnosis",col="blue",lty="solid", haz.scale=1000)


## @knitr 7.g
## Run Poisson regression
summary(poisson7g <- glm( death_cancer ~ fu + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7g, method="Wald")

## @knitr 7.h
summary(poisson7h <- glm( death_cancer ~ fu + year8594 + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7h, method="Wald")

# Add interaction term
summary(poisson7h2 <- glm( death_cancer ~ fu*year8594 + offset( log(pt) ), family=poisson, data=melanoma.spl ))
## IRR
eform(poisson7h2, method="Wald")

## @knitr 7.i

summary(poisson7i <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
eform(poisson7i, method="Wald")

## Test if the effect of age is significant using a likelihood ratio test
drop1(poisson7i, ~agegrp, test="Chisq")
## For this we can also use the car package and a Wald test
linearHypothesis(poisson7i,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))
## ADVANCED:
## Alternative approach for the likelihood ratio test
# poisson7i_2 <- update(poisson7i,. ~ . - agegrp)
# anova(poisson7i_2,poisson7i,test="Chisq")

## @knitr 7.j

summary(poisson7j <- glm( death_cancer ~ fu + agegrp + year8594*sex + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
eform(poisson7j, method="Wald")

## @knitr 7.k.i
# hand calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

## @knitr 7.k.ii
linearHypothesis(poisson7j,c("sexFemale + year8594Diagnosed 85-94:sexFemale = 0"))

## calculate the confidence interval on the linear predictor scale and then transform
exp(deltaMethod(poisson7j, "b14+b15", parameterNames = paste0("b",0:15)))[-2] # NB: drop SE
## We can get the standard error (and drop the confidence interval)
deltaMethod(poisson7j, "exp(b14+b15)", parameterNames = paste0("b",0:15))[-(3:4)]

# ADVANCED:
# hand calculations with confidence intervals
# use estimates with covariance matrix from glm
# to obtain estimates with ci's
linvec <- c(rep(0,14),c(1,1))
coef.Female8594 <- crossprod(coef(poisson7j),linvec)
se.Female8594 <- sqrt(diag(t(linvec)%*%vcov(poisson7j)%*%linvec))
ci.lower <- coef.Female8594 - 1.96*se.Female8594
ci.upper <- coef.Female8594 + 1.96*se.Female8594
exp(c(coef.Female8594,ci.lower,ci.upper))

## @knitr 7.k.iii
## Create dummies and Poisson regression
melanoma.spl <- melanoma.spl %>%

    ## Add confidence intervals for the rates
    mutate(femaleEarly = sex=="Female" & year8594=="Diagnosed 75-84",
           femaleLate = sex=="Female" & year8594=="Diagnosed 85-94")

summary(poisson7k <- glm( death_cancer ~ fu + agegrp + year8594 + femaleEarly +
                         femaleLate + offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))

## IRR
eform(poisson7k, method="Wald")

## @knitr 7.k.iv
## Add interaction term
summary(poisson7k2 <- glm( death_cancer ~ fu + agegrp + year8594 + year8594:sex +
                         offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))

## IRR
eform(poisson7k2, method="Wald")


## @knitr 7.l

summary( poisson7l.early <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 75-84" ) )
eform(poisson7l.early, method="Wald")

summary( poisson7l.late <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 85-94" ) )

eform(poisson7l.late, method="Wald")

# compare with results in i
eform(poisson7i, method="Wald")

# compare with results in j
eform(poisson7j, method="Wald")


# Poisson-regression with effects specific for diagnose period
summary(poisson7l2 <- glm( death_cancer ~ fu + fu:year8594 + agegrp + agegrp:year8594
                          + sex*year8594 + offset( log(pt) ),
                          family=poisson, data=melanoma.spl ))
eform(poisson7l2, method="Wald")
