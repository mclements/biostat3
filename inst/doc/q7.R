## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist
## Modified: Mark Clements, 2017-08-07, 2024-10-27
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
library(car)      # for car::linearHypothesis -> biostat3::lincom
library(knitr)
library(broom)
options(pillar.sigfig = 5) # increase the number of digits printed in tibbles

## @knitr loadPreprocess
## Read melanoma data

## Create a new dataset with only localised cancer
melanoma.l <- subset(biostat3::melanoma, stage=="Localised")
head(melanoma.l)
summary(melanoma.l)

## @knitr 7.a.i

## Plot Kaplan-Meier curve using survfit
## Create a new event indicator
melanoma.l <- transform(melanoma.l,
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
     mark.time=FALSE,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ## Make the plot prettier
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
## Add legend too
legend("bottomleft",legend=levels(melanoma.l$year8594),col=c("blue","red"),
       lty=c("solid","dashed"), bty="n")

### TRY IF YOU WANT ###
if (FALSE) {
    library(survMisc)
    ## Note: `autoplot(sfit7a1)` was broken; I have submitted a pull request to fix this
    ## autoplot(sfit7a1)
    ## alternatively:
    autoplot(sfit7a1, timeTicks = "custom", times= seq(0, 20*12, 5*12))
}

## @knitr 7.a.ii

## To plot smoothed hazards, we use the bshazard package
library(bshazard)
par(mfrow=1:2)
for(level in levels(melanoma.l$year8594))
    plot(bshazard(Surv(surv_mm/12, status == "Dead: cancer") ~ 1,
                  data=subset(melanoma.l, year8594==level)),
         xlab="Years since diagnosis", main=level, xlim=c(0,20), ylim=c(0,0.08))

## @knitr 7.a.iii
## Compare with Kaplan-Meier plot
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sfit7a1,
     ## No automatic labelling of the curve (we do that ourselves)
     mark.time=FALSE,
     ## Time is measured in months,  but we want to see it in years
     xscale=12,
     ylim=c(0.6,1),
     ## Make the plot prettier
     xlab="Years since diagnosis",
     ylab="S(t)",
     col=c("blue","red"),
     lty=c("solid","dashed"))
legend("bottomleft",legend=levels(melanoma.l$year8594),
       col=c("blue","red"),lty=c("solid","dashed"), bty="n")
cols = c("Diagnosed 75-84"="blue", "Diagnosed 85-94"="red")
ltys = c("Diagnosed 75-84"="solid", "Diagnosed 85-94"="dashed")
for(level in levels(melanoma.l$year8594))
    plot(bshazard(Surv(surv_mm/12, status == "Dead: cancer") ~ 1,
                  data=subset(melanoma.l, year8594==level),
                  verbose=FALSE),
         xlab="Years since diagnosis", main=level,
         xlim=c(0,22), ylim=c(0,0.08), col=cols[level],
         lty=ltys[level])

## @knitr 7.b
survRate(Surv(surv_mm/12, death_cancer) ~ year8594, data=melanoma.l) |>
    kable("html")


## @knitr 7.c.i

## Calculate the incidence rate by time of diagnosis
## but with new variables
melanoma.l2 <-
    transform(melanoma.l,
              ## Update the death indicator (only count deaths within 120 months)
              ## death_cancer = death_cancer * as.numeric(surv_mm<=120),
              death_cancer = ifelse(surv_mm<=120 & status == "Dead: cancer",1,0),
              ## Create a new time variable
              surv_mm = pmin(surv_mm,120))

## Calculate the rates on the truncated data
rates_by_diag_yr2 <- survRate(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma.l2)
rates_by_diag_yr2 |> kable("html")

## @knitr 7.c.ii

## MRR full data
rates_by_diag_yr2[2, "rate"] / rates_by_diag_yr2[1, "rate"]
with(rates_by_diag_yr2[2:1,], poisson.test(event, tstop)) # reverse the row order:)

## @knitr 7.c.iii
## Use glm to estimate the rate ratios
## we scale the offset term to 1000 person-years
poisson7c <- glm(death_cancer ~ year8594 + offset(log(surv_mm/12/1000)),
                 family=poisson, data=melanoma.l2)
summary(poisson7c)

## also for collapsed data
summary(glm(event ~ year8594 + offset(log(tstop/12/1000)), family=poisson,
            data=rates_by_diag_yr2))

## IRR - using broom::tidy
tidy(poisson7c, conf.int=TRUE, exponentiate=TRUE)
## or (less generally)
biostat3::eform(poisson7c)

## Note that the scaling of the offset term only has an impact on the intercept
summary(glm(death_cancer ~ year8594 + offset(log(surv_mm)),
            family=poisson, data=melanoma.l2))

## @knitr 7.d

## Add a new variable for year
melanoma.l2 <- transform(melanoma.l2, surv_yy1 = surv_mm/12)

## Split follow up by year
melanoma.spl <- survSplit(melanoma.l2, cut=0:9, end="surv_yy1", start="start",
                           event="death_cancer")

## Calculate persontime and
## recode start time as a factor
melanoma.spl <- transform(melanoma.spl,
                          pt = surv_yy1 - start,
                          fu = as.factor(start) )

## @knitr 7.e

## Calculate the incidence rate by observation year
yearly_rates <- survRate(Surv(pt/1000,death_cancer)~fu, data=melanoma.spl) |>
    transform(start=as.numeric(levels(fu))[fu]) |>
    transform(mid=start+0.5)

library(tinyplot) # lightweight base graphics extension
with(yearly_rates, {
    plt(rate~start, ymin=lower, ymax=upper, type="ribbon",
        main="Cancer deaths by years since diagnosis",
        ylab="Incidence rate per 1000 person-years",
        xlab="Years since diagnosis")
})


## @knitr 7.f
# Plot smoothed hazards

library(bshazard)
library(tinyplot)
par(mfrow=1:2)
survRate(Surv(pt,death_cancer)~fu, data=melanoma.spl) |>
    transform(start=as.numeric(levels(fu))[fu]) |>
    transform(mid=start+0.5) |>
    with({
        plt(rate~mid, ymin=lower, ymax=upper,
            type="ribbon",
            main="Rates",
            ylab="Mortality rate per person-year",
            xlab="Years since diagnosis",
            xlim=c(0,10))
    })
bshazard(Surv(surv_mm/12, status == "Dead: cancer") ~ 1, data = melanoma.l) |>
    plot(xlab="Years since diagnosis", xlim=c(0,10),
         main="Smoothed hazard")
box() # redo the box around the plot

## @knitr 7.g
## Run Poisson regression
summary(poisson7g <- glm(death_cancer ~ fu + offset(log(pt)),
                         family = poisson,
                         data = melanoma.spl))
## IRR
## tidy(poisson7g, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7g) # faster:)

## @knitr 7.h
summary(poisson7h <- glm( death_cancer ~ fu + year8594 + offset( log(pt) ),
                         family = poisson,
                         data = melanoma.spl ))
## IRR
## tidy(poisson7h, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7h)

# Add interaction term
summary(poisson7h2 <- glm(death_cancer ~ fu*year8594 + offset(log(pt)),
                          family=poisson, data=melanoma.spl))
## IRR
## tidy(poisson7h2, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7h2)

## @knitr 7.i

summary(poisson7i <- glm(death_cancer ~ fu + year8594 + sex + agegrp + offset(log(pt)),
                         family=poisson, data=melanoma.spl))

## IRR
## tidy(poisson7i, conf.int=TRUE, exponentiate=TRUE)
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

summary(poisson7j <- glm(death_cancer ~ fu + agegrp + year8594*sex + offset(log(pt)),
                         family=poisson, data=melanoma.spl))

## IRR
## tidy(poisson7j, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7j)

## @knitr 7.k.i
# hand calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

## @knitr 7.k.ii
## You will need the "car" package to use lincom. If it is not already installed:
## install.packages("car")
lincom(poisson7j,c("sexFemale + year8594Diagnosed 85-94:sexFemale"),eform=TRUE)


## @knitr 7.k.iii
## Create dummies and Poisson regression
melanoma.spl <- melanoma.spl |>
    ## Add confidence intervals for the rates
    transform(femaleEarly = sex=="Female" & year8594=="Diagnosed 75-84",
              femaleLate = sex=="Female" & year8594=="Diagnosed 85-94")

summary(poisson7k <- glm( death_cancer ~ fu + agegrp + year8594 + femaleEarly +
                         femaleLate + offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))

## IRR
## tidy(poisson7k, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7k)

## @knitr 7.k.iv
## Add interaction term
summary(poisson7k2 <- glm( death_cancer ~ fu + agegrp + year8594 + year8594:sex +
                         offset( log(pt) ), family=poisson,
                         data=melanoma.spl ))
## tidy(poisson7k2, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7k2)

## @knitr 7.l

summary( poisson7l.early <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 75-84" ) )
## tidy(poisson7l.early, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7l.early)

summary( poisson7l.late <- glm( death_cancer ~ fu + agegrp + sex + offset( log(pt) ),
                       family = poisson, data = melanoma.spl,
                       subset = year8594 == "Diagnosed 85-94" ) )
## tidy(poisson7l.late, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7l.late)

# compare with results in i
## tidy(poisson7i, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7i)

# compare with results in j
## tidy(poisson7j, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7j)


# Poisson-regression with effects specific for diagnose period
summary(poisson7l2 <- glm( death_cancer ~ fu + fu:year8594 + agegrp + agegrp:year8594
                          + sex*year8594 + offset( log(pt) ),
                          family=poisson, data=melanoma.spl ))
## tidy(poisson7l2, conf.int=TRUE, exponentiate=TRUE)
eform(poisson7l2)

## @knitr 7.m
## Split follow up by month
library(splines)
library(tinyplot)
time.cut <- seq(0,10,by=1/12)
nrow(biostat3::melanoma)
melanoma.spl2 <- survSplit(Surv(surv_mm/12,status=="Dead: cancer")~.,
                           data=biostat3::melanoma,
                           cut=time.cut,
                           subset=stage=="Localised") |>
    transform(mid=(tstop+tstart)/2, risk_time=tstop-tstart)
nrow(melanoma.spl2)
poisson7m <- glm(event ~ ns(mid,df=6) + agegrp + year8594 +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl2)
df <- data.frame(agegrp="0-44", year8594="Diagnosed 75-84",
                 mid=time.cut[-1], risk_time=1)
pred <- predict(poisson7m, newdata=df, se.fit=TRUE)
qq <- qnorm(0.975)
df <- with(pred, transform(df,
                           fit=exp(fit), # parallel (not sequential:)
                           conf.low=exp(fit-qq*se.fit),
                           conf.high=exp(fit+qq*se.fit)))
## plot the rate at the baseline values
with(df, plt(fit~mid, ymin=conf.low, ymax=conf.high, type="ribbon",
     ylab="Rate", xlab="Time since diagnosis (years)",
     ylim=c(0,0.05)))

## @knitr 7.n
## using melanoma.spl2 and df from previous chunk
poisson7n <- glm(event ~ ns(mid,df=4) + agegrp + year8594 +
                     ifelse(year8594=="Diagnosed 85-94",1,0):ns(mid,df=3) +
                     offset(log(risk_time)),
                 family=poisson,
                 data=melanoma.spl2)
library(rstpm2)
library(tinyplot)
## get log(RR) confidence interval using predictnl (delta method)
df <- data.frame(agegrp="0-44", year8594="Diagnosed 75-84",
                 mid=time.cut[-1], risk_time=1)
predictnl(poisson7n, function(object)
    predict(object, newdata=transform(df, year8594="Diagnosed 85-94"), type="link") -
    predict(object, newdata=df, type="link")) |>
    cbind(df) |>
    as.data.frame() |>
    transform(fit = exp(fit),
              conf.low=exp(fit-1.96*se.fit),
              conf.high=exp(fit+1.96*se.fit)) |>
    with(plt(fit~mid, ymin=conf.low, ymax=conf.high, type="ribbon",
             xlab="Time since diagnosis (years)",
             ylab="Rate ratio"))

predictnl(poisson7n,
          function(object)
              predict(object, newdata=transform(df, year8594="Diagnosed 85-94"),
                      type="response") -
              predict(object, newdata=df, type="response")) |>
    cbind(df) |>
    transform(conf.low=fit-1.96*se.fit,
              conf.high=fit+1.96*se.fit) |>
    with(plt(fit~mid, ymin=conf.low, ymax=conf.high, type="ribbon",
             xlab="Time since diagnosis (years)",
             ylab="Rate difference"))

## @knitr 7.o
## Calculate survival from a smooth Poisson regression model 
if (requireNamespace("deSolve")) {
    twoState <- function(object, ...) {
        markov_msm(list(object),matrix(c(NA,1,NA,NA),2,byrow=TRUE), ...) |>
            as.data.frame() |>
            subset(state==1)
    }
    df3 <- expand.grid(agegrp=levels(biostat3::melanoma$agegrp),
                       year8594=levels(biostat3::melanoma$year8594)) |>
        transform(risk_time=1)
    library(tinyplot)
    twoState(poisson7n, t=c(0,time.cut), newdata = df3, tmvar = "mid") |>
        with(plt(P~time|year8594, ymin=P.lower,ymax=P.upper,
                 type="ribbon",facet=agegrp,
                 xlab="Time since diagnosis (years)",
                 ylab="Survival"))
} else cat("To run this example, please install the deSolve package\n")

## @knitr 7.o.b
if (requireNamespace("pracma")) {
library(tinyplot)
library(rstpm2)
library(pracma)
df3 <- with(biostat3::melanoma,
            expand.grid(agegrp=levels(agegrp),
                        year8594=levels(year8594))) |>
    transform(risk_time=1)
oneState <- function(object,newdata=df3) {
    f <- function(t,y)
        -y*predict(object, newdata=transform(newdata,mid=t), type="response")
    pracma::ode45(f, t0=0, y0=rep(1,nrow(newdata)), tfinal=10, hmax=0.1)
}
out0 = oneState(poisson7n, newdata=df3)
out = predictnl(poisson7n,
                function(object, newdata)
                    log(-log(as.vector(oneState(object, newdata)$y))),
                newdata=df3, conf.int=TRUE) |> # slow
    transform(fit=exp(-exp(fit)),
              conf.low=exp(-exp(fit-1.96*se.fit)),
              conf.high=exp(-exp(fit+1.96*se.fit))) |>
    cbind(with(biostat3::melanoma,
               expand.grid(mid=out0$t,
                           agegrp=levels(agegrp),
                           year8594=levels(year8594))))
with(out,plt(fit~mid|year8594,type="ribbon",facet=agegrp,
             ymin=conf.low, ymax=conf.high,
             xlab="Time since diagnosis (years)",
             ylab="Survival"))
} else cat("To run this example, please install the pracma package\n")



