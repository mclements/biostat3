## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist
###############################################################################

###############################################################################
## Exercise 7
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("muhaz")
## install.packages("car")

## @knitr loadDependecies
require(survival) # for Surv and survfit
require(dplyr)    # for data manipulation
require(foreign)  # for reading data set from Stata


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
## Read melanoma data
melanoma <- tbl_df( read.dta("http://biostat3.net/download/melanoma.dta") )

## Create a new dataset with only localised cancer
melanoma.l <- filter(melanoma, stage=="Localised")
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

## To plot smoothed hazards, we need use the muhaz package
require(muhaz)

## The plot.muhaz cannot rescale so we create a new variable
## Create a new variable base on the number of months
## translated to years (more accurate than surv_yy)
melanoma.l <- mutate( melanoma.l, surv_yy1 = surv_mm/12)

## Since the muhaz interface does not support data frames
## we do this the old school way
hazfit7a.early <- with(melanoma.l,
                       muhaz(times=surv_yy1,
                             delta=status=="Dead: cancer",
                             subset=year8594=="Diagnosed 75-84"))

hazfit7a.late <- with(melanoma.l,
                      muhaz(times=surv_yy1,
                            delta=status=="Dead: cancer",
                            subset=year8594=="Diagnosed 85-94",
                            max.time=max(surv_yy1)))

## With bells and whistles
plot(hazfit7a.early,xlab="Years since diagnose",col="blue",lty="solid")
lines(hazfit7a.late,col="red",lty="dashed")
legend("topright",legend=levels(melanoma$year8594),col=c("blue","red"),lty=c("solid","dashed"))

## @knitr 7.a.iii
## Compare with Kaplan-Meier plot
par(mfrow=c(1,2)) ## Two graphs in the same window

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
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),lty=c("solid","dashed"))

plot(hazfit7a.early,xlab="Years since diagnose",col="blue",lty="solid")
lines(hazfit7a.late,col="red",lty="dashed")
legend("topright",legend=levels(melanoma$year8594),col=c("blue","red"),lty=c("solid","dashed"))

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
        Rate = D/Y) %>%

    ## Add confidence intervals for the rates
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

rates_by_diag_yr

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
        Rate = D/Y) %>%

    ## Add confidence intervals for the rates
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

rates_by_diag_yr2

rates_by_diag_yr

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
IRR(poisson7c)

## with confidence intervals (takes some time to compute)
## exp(cbind(coef(poisson7c),confint(poisson7c)))


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
        Rate = D/Y) %>%

    ## Add confidence intervals for the rates
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

yearly_rates

## Plot by year
with(yearly_rates, matplot(fu,
                           cbind(Rate, CI_low,
                             CI_high),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths per 1000 person-year by years since diagnose",
                           ylab="IR",
                           xlab="Years since diagnose") )

## @knitr 7.f
hazfit7f <- with(melanoma.l2, muhaz(times=surv_yy1, delta=death_cancer, max.time=10.5) )

# Plot smoothed hazards
frame() # Open new graphics device (window)
plot(hazfit7f,xlab="Years since diagnose",col="blue",lty="solid")

par(mfrow=c(1,2))
with(yearly_rates, matplot(fu,
                           cbind(Rate, CI_low,
                             CI_high),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths per 1000 person-year by years since diagnose",
                           ylab="IR",
                           xlab="Years since diagnose") )
plot(hazfit7f,xlab="Years since diagnose",col="blue",lty="solid")


## @knitr 7.g
## Run Poisson regression
summary(poisson7g <- glm( death_cancer ~ fu + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
IRR(poisson7g)

## @knitr 7.h
summary(poisson7h <- glm( death_cancer ~ fu + year8594 + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
IRR(poisson7h)

# Add interaction term
summary(poisson7h2 <- glm( death_cancer ~ fu*year8594 + offset( log(pt) ), family=poisson, data=melanoma.spl ))
## IRR
IRR(poisson7h2)

## @knitr 7.i

summary(poisson7i <- glm( death_cancer ~ fu + agegrp + year8594 + sex + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
IRR(poisson7i)

## Test if the effect of age is significant
## For this we use the car package
require(car)
linearHypothesis(poisson7i,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))
## ADVANCED:
## Alternative by comparing deviances
## poisson7i_2 <- update(poisson7i,. ~ . - agegrp)
## anova(poisson7i_2,poisson7i,test="Chisq")

## @knitr 7.j

summary(poisson7j <- glm( death_cancer ~ fu + agegrp + year8594*sex + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
IRR(poisson7j)

## @knitr 7.k.i
# hand calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

## @knitr 7.k.ii
linearHypothesis(poisson7j,c("sexFemale = 0","year8594Diagnosed 85-94:sexFemale = 0"))

# ADVANCED:
# hand calculations with confidence intervals
# use estimates with covariance matrix from glm
# to obtain estimates with ci's
length(coef(poisson7j))
linvec <- c(rep(0,14),c(1,1))
coef.Female8594 <- crossprod(coef(poisson7j),linvec)
var.Female8594 <- t(linvec)%*%vcov(poisson7j)%*%linvec
ci.d <- 1.96*sqrt(var.Female8594)
ci.lower <- coef.Female8594 - ci.d
ci.upper <- coef.Female8594 + ci.d
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
IRR(poisson7k)

## @knitr 7.k.iv
## Add interaction term
summary(poisson7k2 <- glm( death_cancer ~ fu + agegrp + year8594 + year8594:sex +
                         offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))

## IRR
IRR(poisson7k2)


## @knitr 7.l

melanoma.spl %>%

    group_by(year8594) %>%

    summarize()

summary( poisson7l.early <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = which(year8594 == "Diagnosed 75-84" ) ) )
IRR(poisson7l.early)

summary( poisson7l.late <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = which(year8594 == "Diagnosed 85-94" ) ) )

IRR(poisson7l.late)

# compare with results in i
IRR(poisson7i)

# compare with results in j
IRR(poisson7j)


# Poisson-regression with effects specific for diagnose period
summary(poisson7l2 <- glm( death_cancer ~ fu + fu:year8594 + agegrp + agegrp:year8594
                          + sex*year8594 + offset( log(pt) ),
                          family=poisson, data=melanoma.spl ))
IRR(poisson7l2)
