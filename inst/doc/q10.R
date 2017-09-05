## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("foreign")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("muhaz")

###############################################################################
## Exercise 10
###############################################################################
## @knitr loadDependecies
library(biostat3)
library(dplyr)    # for data manipulation
library(ggplot2)


## @knitr loadPreprocess
melanoma_raw <- biostat3::melanoma
melanoma <- melanoma_raw %>%
    filter(stage == "Localised") %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer" & surv_mm <= 120, 1, 0), #censuring for > 120 monts
           trunc_yy = ifelse(surv_mm <=  120, surv_mm/12, 10))  #scale to years and truncate to 10 years

## @knitr 10.a
# Using muhaz2 to smooth the Kaplan-Meier hazards by strata; the kernel and bandwidth options were selected based on smoother performance.
hazDiaDate <- muhaz2(Surv(trunc_yy,death_cancer)~year8594, data=melanoma, bw.method="g", bw.grid=5)
hazDiaDateDf <- as.data.frame(hazDiaDate)

## Max hazard ratio
maxHaz <- hazDiaDateDf %>% group_by(strata) %>%
    summarise(stratMax=max(haz.est))
print(maxHaz)
maxHaz$stratMax[2]/maxHaz$stratMax[1]

## Comparing hazards
plot(hazDiaDate)
legend("topright", legend=levels(melanoma$year8594), col=1:2, lty=1)
# or using ggplot2
ggplot(hazDiaDateDf, aes(x=est.grid, y=haz.est, colour= strata)) + geom_line()

## @knitr 10.b
## Comparing hazards on a log scales
ggplot(hazDiaDateDf, aes(x=est.grid, y=haz.est, colour= strata)) + geom_line() +
  scale_y_continuous(trans='log')

## @knitr 10.c
## Calculating -log cumulative hazards per strata
survfit1 <- survfit(Surv(trunc_yy,death_cancer)~year8594, data=melanoma)
plot(survfit1, col=1:2, fun=function(S) -log(-log(S)), log="x",
     xlab="log(time)", ylab="-log(H)")
legend("topright",legend=levels(melanoma$year8594),col=1:2,lty=1)

## @knitr 10.d
# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence
# whereas Stata (and most other software) use the Breslow method
cox1 <- coxph(Surv(trunc_yy, death_cancer==1) ~ year8594, method=c("breslow"), data=melanoma)
summary(cox1)


## @knitr 10.e

cox2 <- coxph(Surv(trunc_yy, death_cancer==1) ~ sex + year8594 + agegrp, method=c("breslow"), data=melanoma)
summary(cox2)

## Plot of the scaled Schoenfeld residuals for calendar period 1985–94.
## The smooth line shows the estimated log hazard ratio as a function of time.
cox2.phtest <- cox.zph(cox2, transform="identity") #Stata appears to be using 'identity'
plot(cox2.phtest[2],resid=TRUE, se=TRUE, main="Schoenfeld residuals", ylim=c(-4,4))

## @knitr 10.g
## The results from the previous proportional hazards assumption test
print(cox2.phtest)

## @knitr 10.h
melanoma2p8Split <- survSplit(melanoma, cut=c(2), end="trunc_yy", start="start", event="death_cancer")
melanoma2p8Split <- mutate(melanoma2p8Split, fu = as.factor(start))

##Tabulate ageband including risk_time
melanoma2p8Split %>% select(id, start, trunc_yy) %>% filter(id<=3) %>% arrange(id, trunc_yy)

head(melanoma2p8Split)

cox2p8Split1 <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + agegrp*fu, method=c("breslow"), data=melanoma2p8Split)
summary(cox2p8Split1)

cox2p8Split1b <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + agegrp + I(agegrp=="75+" & fu=="2"), method=c("breslow"), data=melanoma2p8Split)
summary(cox2p8Split1b)


## @knitr 10.i
cox2p8Split2 <- coxph(Surv(start, trunc_yy, death_cancer) ~ sex + year8594 + fu + fu:agegrp, method=c("breslow"), data=melanoma2p8Split)
summary(cox2p8Split2)

## Alternative approach using tt():
## http://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
cox2p8tvc2 <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp + tt(agegrp), data=melanoma,
                   tt = function(x, t, ...) (x=="75+")*(t>=2))
summary(cox2p8tvc2)

cox2p8tvct <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp + tt(agegrp), data=melanoma,
                   tt = function(x, t, ...) (x=="75+")*t)
summary(cox2p8tvct)

cox2p8tvclogt <- coxph(Surv(trunc_yy, death_cancer) ~ sex + year8594 + agegrp + tt(agegrp), data=melanoma,
                    tt = function(x, t, ...) (x=="75+")*log(t))
summary(cox2p8tvclogt)
