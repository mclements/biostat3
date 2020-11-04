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
## @knitr loadDependencies

library(biostat3)
library(dplyr)    # for data manipulation
library(splines)   # ns (recommended package)

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

survdiff(Surv(surv_mm, death_cancer) ~ year8594, data = melanoma.l)
summary( coxfit9b <- coxph(Surv(surv_mm, death_cancer) ~ year8594,
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
anova(coxfit9c)
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

## Split follow up by event times
melanoma.spl2 <- survSplit(melanoma.l2, cut=sfit9f$time, end="surv_mm", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl2 <- mutate(melanoma.spl2,
                        pt = surv_mm - start,
                        fu = as.factor(start) )

## Collapse
melanoma.spl3 <- melanoma.spl2 %>%
    group_by(fu,year8594,sex,agegrp) %>%
    summarise(pt=sum(pt), death_cancer=sum(death_cancer)) %>%
    data.frame

## Run Poisson regression
poisson9f <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ),
                 family = poisson,
                 data = melanoma.spl3 )

## IRR
coef9f <- eform(poisson9f)
## We are not interested in nuisance parameters fu1, fu2, etc
npar <- nrow(coef9f)
pars <- (npar-4):npar
coef9f[pars,]

summary(coxfit9e)


## @knitr 9.g
## split and collapse
cuts.splines <- seq(0,max(sfit9f$time),by=3)
mid.splines <- cuts.splines + 1.5
melanoma.spl4 <-
    survSplit(Surv(surv_mm,death_cancer)~., data=melanoma.l2, cut=cuts.splines,
              start="tstart", end="tstop") %>%
    mutate(cut=cut(tstop,cuts.splines),
           mid=mid.splines[unclass(cut)]) %>%
    group_by(mid,year8594,sex,agegrp) %>%
    summarise(pt=sum(tstop-tstart), death_cancer=sum(death_cancer))

poisson9g <- glm( death_cancer ~ ns(mid,df=3) + year8594 + sex + agegrp + offset( log(pt) ),
                 family = poisson,
                 data = melanoma.spl4 )
summary(poisson9g)
eform(poisson9g)

## quick approach: use the effects package
library(effects)
plot(Effect("mid", poisson9g))

## utility function to draw a confidence interval
polygon.ci <- function(time, interval, col="lightgrey") 
    polygon(c(time,rev(time)), c(interval[,1],rev(interval[,2])), col=col, border=col)

## define exposures
times <- seq(0,max(cuts.splines),length=1000)
delta <- diff(times)[1]
newdata <- data.frame(mid=times+delta/2, year8594="Diagnosed 85-94",
                      sex="Male", agegrp="45-59",
                      pt=1)
## predict rates and 95% CI
pred <- predict(poisson9g, newdata=newdata, se.fit=TRUE)
predrate <- exp(pred$fit)
ci <- with(pred, exp(cbind(fit-1.96*se.fit, fit+1.96*se.fit)))
## plot
matplot(newdata$mid, ci, type="n", xlab="Time since diagnosis (months)",
        ylab="Rate", main="Males aged 45-59 years diagnosed 1985-94")
polygon.ci(newdata$mid, ci) 
lines(newdata$mid, predrate)

## predict survival and 95% CI
library(rstpm2)
logcumhaz <- rstpm2::predictnl(poisson9g,
          fun=function(fit,newdata) log(cumsum(delta*predict(fit, newdata, type="response"))),
          newdata=newdata)
surv <- exp(-exp(logcumhaz$Estimate))
ci <- with(logcumhaz, exp(-exp(cbind(Estimate-1.96*SE, Estimate+1.96*SE))))

## plot
matplot(newdata$mid, ci, type="n", xlab="Time since diagnosis (months)",
        ylab="Survival", main="Males aged 45-59 years diagnosed 1985-94")
polygon.ci(newdata$mid, ci) 
lines(newdata$mid, surv)
