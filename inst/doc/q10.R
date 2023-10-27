## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
## Revised: Mark Clements, 2017-11-03
###############################################################################

###############################################################################
## Exercise 10
###############################################################################
## @knitr loadDependencies
library(biostat3)
library(dplyr)    # for data manipulation
library(ggplot2)
library(car)      # car::linearHypothesis -> biostat3::lincom

## @knitr loadPreprocess
localised <- dplyr::filter(biostat3::melanoma, stage == "Localised") %>%
    dplyr::mutate(death_cancer = ifelse( status == "Dead: cancer" & surv_mm <= 120, 1, 0), #censoring for > 120 months
           trunc_yy = pmin(surv_mm/12,10))  #scale to years and truncate to 10 years

## @knitr 10.a
# Using muhaz2 to smooth the Kaplan-Meier hazards by strata
hazDiaDate <- muhaz2(Surv(trunc_yy,death_cancer)~year8594, data=localised)
hazDiaDateDf <- as.data.frame(hazDiaDate)

## Comparing hazards
plot(hazDiaDate, haz.scale=1000,
     xlab="Time since diagnosis (years)", 
     ylab="Hazard per 1000 person-years")
# or using ggplot2
ggplot(hazDiaDateDf, aes(x=x, y=y*1000, colour= strata)) + geom_line() +
    xlab("Time since diagnosis (years)") +
    ylab("Hazard per 1000 person-years")


## @knitr 10.b
## Comparing hazards on a log scales
plot(hazDiaDate, log="y")


## @knitr 10.c
## Calculating -log cumulative hazards per strata
survfit1 <- survfit(Surv(trunc_yy,death_cancer)~year8594, data=localised)
plot(survfit1, col=1:2, fun=function(S) -log(-log(S)), log="x",
     xlab="log(time)", ylab="-log(H)")
legend("topright",legend=levels(localised$year8594),col=1:2,lty=1)

## or we can use
biostat3::survPHplot(Surv(trunc_yy,death_cancer)~year8594, data=localised)

## @knitr 10.d
# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence of ties
# whereas Stata (and some other software) use the Breslow method
cox1 <- coxph(Surv(trunc_yy, death_cancer==1) ~ year8594, data=localised)
summary(cox1)


## @knitr 10.e

cox2 <- coxph(Surv(trunc_yy, death_cancer==1) ~ sex + year8594 + agegrp, data=localised)
summary(cox2)

## Plot of the scaled Schoenfeld residuals for calendar period 1985–94.
## The smooth line shows the estimated log hazard ratio as a function of time.
cox2.phtest <- cox.zph(cox2, transform="identity") #Stata appears to be using 'identity'
plot(cox2.phtest[2],resid=TRUE, se=TRUE, main="Schoenfeld residuals", ylim=c(-4,4))

## @knitr 10.g
## The results from the previous proportional hazards assumption test
print(cox2.phtest)

## @knitr 10.h
melanoma2p8Split <- survSplit(localised, cut=c(2), end="trunc_yy", start="start",
                              event="death_cancer", episode="fu") %>%
    mutate(fu = as.factor(fu))

##Tabulate ageband including risk_time
melanoma2p8Split %>% select(id, start, trunc_yy) %>% filter(id<=3) %>% arrange(id, trunc_yy)

head(melanoma2p8Split)

cox2p8Split1 <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + agegrp*fu,
                      data=melanoma2p8Split)
summary(cox2p8Split1)

cox2p8Split1b <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + agegrp +
                           I(agegrp=="45-59" & fu=="2") + I(agegrp=="60-74" & fu=="2") +
                           I(agegrp=="75+" & fu=="2"), data=melanoma2p8Split)
summary(cox2p8Split1b)


## @knitr 10.i
cox2p8Split2 <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + fu + fu:agegrp, data=melanoma2p8Split)
summary(cox2p8Split2)

## @knitr 10.ib

## Alternative approach using tt():
## http://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
cox2p8tvc2 <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp +
                        tt(agegrp=="45-59") + tt(agegrp=="60-74") + tt(agegrp=="75+"),
                    data=localised,
                    tt = function(x, t, ...) x*(t>=2))
summary(cox2p8tvc2)
## The tt labels do not play nicely with lincom:(

cox2p8tvct <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp + tt(agegrp),
                    data=localised,
                    tt = function(x, t, ...) cbind(`45-59`=(x=="45-59")*t,
                                                   `60-74`=(x=="60-74")*t,
                                                   `75+`=(x=="75+")*t))
summary(cox2p8tvct)
lincom(cox2p8tvct, "agegrp75+",eform=TRUE)                   # t=0
lincom(cox2p8tvct, "agegrp75+ + tt(agegrp)75+",eform=TRUE)   # t=1
lincom(cox2p8tvct, "agegrp75+ + 2*tt(agegrp)75+",eform=TRUE) # t=2

cox2p8tvclogt <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp + 
                           tt(agegrp),
                       data=localised,
                    tt = function(x, t, ...) cbind(`45-59`=(x=="45-59")*log(t),
                                                   `60-74`=(x=="60-74")*log(t),
                                                   `75+`=(x=="75+")*log(t)))
summary(cox2p8tvclogt)
lincom(cox2p8tvclogt, "agegrp75+ - 0.6931472*tt(agegrp)75+",eform=TRUE) # t=0.5 => log(t)=-0.6931472
lincom(cox2p8tvclogt, "agegrp75+", eform=TRUE) # t=1 => log(t)=0
lincom(cox2p8tvclogt, "agegrp75+ + 0.6931472*tt(agegrp)75+",eform=TRUE) # t=2 => log(t)=0.6931472



## @knitr 10.j
library(splines)
time.cuts <- seq(0,10,length=100)
delta <- diff(time.cuts)[1]
## split and collapse
melanoma2p8Split2 <- survSplit(Surv(trunc_yy,death_cancer)~., data=localised,
                               cut=time.cuts, end="tstop", start="tstart",
                               event="death_cancer") %>%
    mutate(fu=cut(tstop,time.cuts),
           mid=time.cuts[unclass(fu)]+delta/2) %>%
    group_by(mid,sex,year8594,agegrp) %>%
    summarise(pt=sum(tstop-tstart), death_cancer=sum(death_cancer)) %>%
    mutate(age75 = (agegrp=="75+")+0)

poisson2p8tvc <- glm(death_cancer ~ sex + year8594 + agegrp + ns(mid,df=3) +
                         age75:ns(mid,df=3) + offset(log(pt)),
    data=melanoma2p8Split2, family=poisson)

## utility function to draw a confidence interval
polygon.ci <- function(time, interval, col="lightgrey") 
    polygon(c(time,rev(time)), c(interval[,1],rev(interval[,2])), col=col, border=col)
## define exposures
newdata <- data.frame(mid=seq(0,max(time.cuts),length=100), year8594="Diagnosed 85-94",
                      sex="Male", agegrp="75+", age75=1, pt=1)

library(rstpm2)
logirr <- rstpm2::predictnl(poisson2p8tvc,
          fun=function(fit,newdata) predict(fit, newdata) -
                        predict(fit, transform(newdata, agegrp='0-44', age75=0)),
          newdata=newdata)
pred <- exp(logirr$fit)
ci <- exp(confint(logirr))
## plot
matplot(newdata$mid, ci, type="n", xlab="Time since diagnosis (months)",
        ylab="Rate ratio", main="Ages 75+ compared with ages 0-44 years")
polygon.ci(newdata$mid, ci) 
lines(newdata$mid, pred)
