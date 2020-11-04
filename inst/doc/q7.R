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

## @knitr loadDependencies
library(biostat3)
library(dplyr)    # for data manipulation
library(car)      # for linearHypothesis


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
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
## Add legend too
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),lty=c("solid","dashed"), bty="n")

### TRY IF YOU WANT ###
if (FALSE) {
    library(survMisc)
    ## Note: `autoplot(sfit7a1)` was broken; I have submitted a pull request to fix this
    ## autoplot(sfit7a1)
    ## alternatively:
    autoplot(sfit7a1, timeTicks = "custom", times= seq(0, 20*12, 5*12))
}

## @knitr 7.a.ii

## To plot smoothed hazards, we use the muhaz package (using the muhaz2 wrapper)

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"), lty=1:2)


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

plot(muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ year8594, data=melanoma.l),
     xlab="Years since diagnosis", col=c("blue","red"),lty=c("solid","dashed"))

## @knitr 7.b

survRate(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma.l)


## @knitr 7.c.i

## Calculate the incidence rate by time of diagnosis
## but with new variables
melanoma.l2 <- mutate(melanoma.l,
                      ## Update the death indicator (only count deaths within 120 months)
                      ## death_cancer = death_cancer * as.numeric(surv_mm<=120),
                      death_cancer = ifelse(surv_mm<=120 & status == "Dead: cancer",1,0),
                      ## Create a new time variable
                      ## surv_mm = pmin(surv_mm, 120)
                      surv_mm = ifelse(surv_mm<=120, surv_mm, 120)
                      )

## Calculate the rates on the truncated data
rates_by_diag_yr2 <- survRate(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma.l2)
rates_by_diag_yr2

## @knitr 7.c.ii

## MRR full data
rates_by_diag_yr2[2, "rate"] / rates_by_diag_yr2[1, "rate"]
with(rates_by_diag_yr2[2:1,], poisson.test(event, tstop))

## @knitr 7.c.iii
## Use glm to estimate the rate ratios
## we scale the offset term to 1000 person-years
poisson7c <- glm( death_cancer ~ year8594 + offset( log( surv_mm/12/1000 ) ), family=poisson, data=melanoma.l2 )
summary( poisson7c )

## also for collapsed data
summary(glm( event ~ year8594 + offset( log( tstop/12/1000 ) ), family=poisson, data=rates_by_diag_yr2))


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
yearly_rates <- survRate(Surv(pt/1000,death_cancer)~fu, data=melanoma.spl)

## Plot by year
with(yearly_rates, matplot(fu,
                           cbind(rate, lower,
                             upper),
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
                           cbind(rate, lower,
                                 upper),
                           ylim=c(0,max(upper)),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by time since diagnosis",
                           ylab="Mortality rate per 1000 person-years",
                           xlab="Years since diagnosis") )

hazfit7f <- muhaz2(Surv(surv_mm/12, status == "Dead: cancer") ~ 1, data = melanoma.l)
## scale hazard by 1000
plot(hazfit7f, xlab="Years since diagnosis",col="blue",lty="solid", haz.scale=1000, xlim=c(0,10))


## @knitr 7.g
## Run Poisson regression
summary(poisson7g <- glm( death_cancer ~ fu + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7g)

## @knitr 7.h
summary(poisson7h <- glm( death_cancer ~ fu + year8594 + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
eform(poisson7h)

# Add interaction term
summary(poisson7h2 <- glm( death_cancer ~ fu*year8594 + offset( log(pt) ), family=poisson, data=melanoma.spl ))
## IRR
eform(poisson7h2)

## @knitr 7.i

summary(poisson7i <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ), family=poisson, data=melanoma.spl ))

## IRR
eform(poisson7i)

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
eform(poisson7j)

## @knitr 7.k.i
# hand calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

## @knitr 7.k.ii
biostat3::lincom(poisson7j,c("sexFemale + year8594Diagnosed 85-94:sexFemale"),eform=TRUE)


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
eform(poisson7k)

## @knitr 7.k.iv
## Add interaction term
summary(poisson7k2 <- glm( death_cancer ~ fu + agegrp + year8594 + year8594:sex +
                         offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))
eform(poisson7k2)


## @knitr 7.l

summary( poisson7l.early <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 75-84" ) )
eform(poisson7l.early)

summary( poisson7l.late <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 85-94" ) )
eform(poisson7l.late)

# compare with results in i
eform(poisson7i)

# compare with results in j
eform(poisson7j)


# Poisson-regression with effects specific for diagnose period
summary(poisson7l2 <- glm( death_cancer ~ fu + fu:year8594 + agegrp + agegrp:year8594
                          + sex*year8594 + offset( log(pt) ),
                          family=poisson, data=melanoma.spl ))
eform(poisson7l2)

## @knitr 7.m
## Split follow up by month
library(splines)
time.cut <- seq(0,10,by=1/12)
nrow(biostat3::melanoma)
melanoma.spl <- survSplit(Surv(surv_mm/12,status=="Dead: cancer")~., data=biostat3::melanoma,
                          cut=time.cut,
                          subset=stage=="Localised")
nrow(melanoma.spl)
melanoma.spl <- transform(melanoma.spl, mid=(tstop+tstart)/2, risk_time=tstop-tstart)
poisson7m <- glm(event ~ ns(mid,df=6) + agegrp + year8594 +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl)
df <- data.frame(agegrp="0-44", year8594="Diagnosed 75-84",
                 mid=time.cut[-1], risk_time=1)
## plot the rate at the baseline values
plot(df$mid, predict(poisson7m, newdata=df, type="response"),
     type="l", ylab="Rate", xlab="Time since diagnosis (years)",
     ylim=c(0,0.05))

## @knitr 7.n
## using melanoma.spl and df from previous chunk
poisson7n <- glm(event ~ ns(mid,df=4) + agegrp + year8594 +
                     ifelse(year8594=="Diagnosed 85-94",1,0):ns(mid,df=3) +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl)
library(rstpm2)
## get log(RR) confidence interval using predictnl (delta method)
pred <- predictnl(poisson7n, function(object)
    log(predict(object, newdata=transform(df, year8594="Diagnosed 85-94"),
                type="response") /
        predict(object, newdata=df, type="response")))
gpred2 <- transform(pred, time = df$mid, rr = exp(fit), ci = exp(confint(pred)))
ggplot(pred2, aes(x=time, y=rr, ymin=ci.2.5.., ymax=ci.97.5..)) +
    ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=0.6) +
    xlab("Time since diagnosis (years)") +
    ylab("Rate ratio")
    
## Calculate the rate difference
band <- function(x,yy,col="grey")
    polygon(c(x,rev(x)), c(yy[,1], rev(yy[,2])), col=col, border=col)
pred <- predictnl(poisson7n, function(object)
    predict(object, newdata=transform(df, year8594="Diagnosed 85-94"),
                type="response") -
    predict(object, newdata=df, type="response"))
rd <- pred$fit
ci <- confint(pred)
matplot(df$mid,
        ci,
        type="n", # blank plot area
        xlab="Time since diagnosis (years)",
        ylab="Rate difference") 
band(df$mid,ci) # add confidence band
lines(df$mid, rd) # add rate difference

## @knitr 7.o
## Calculate survival from a smooth Poisson regression model 
twoState <- function(object, ...) {
    out <- as.data.frame(markov_msm(list(object),matrix(c(NA,1,NA,NA),2,byrow=TRUE), ...))
    transform(subset(out, state==1),
              S=P,
              S.lower=P.lower,
              S.upper=P.upper)
}
df2 <- expand.grid(agegrp=levels(biostat3::melanoma$agegrp),
                   year8594=levels(biostat3::melanoma$year8594))
df2 <- transform(df2,risk_time=1)
pred <- twoState(poisson7n, t=c(0,df$mid), newdata = df2, tmvar = "mid")
ggplot(pred, aes(x=time,y=S,ymin=S.lower,ymax=S.upper,fill=year8594)) +
    ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=0.6) +
    facet_grid(~agegrp) +
    xlab("Time since diagnosis (years)") +
    ylab("Survival")
