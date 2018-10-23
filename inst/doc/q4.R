## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28
## Edited: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 4
###############################################################################
## @knitr loadDependencies
library(biostat3) 

## @knitr loadPreprocess
data(melanoma)
localised <- subset(melanoma, stage=="Localised")
localised <- transform(localised,
                       death_cancer = ifelse(status == "Dead: cancer", 1, 0))

## @knitr actuarialYears
lifetab2(Surv(floor(surv_yy),death_cancer)~1, data = localised)[,1:7]

## @knitr actuarialMonths
lifetab2(Surv(floor(surv_mm),death_cancer)~1, data = localised)[110:130,1:7]

## @knitr kmYears
mfit_years <- survfit(Surv(surv_yy, death_cancer) ~ 1, data = localised)
summary(mfit_years)

## @knitr kmMonths
mfit_months <- survfit(Surv(surv_mm, death_cancer) ~ 1, data = localised)
summary(mfit_months,times=110:130)

## @knitr comparisonPlot
plot(mfit_months, conf.int=FALSE,
     ylim=c(0.7,1), cex=5,
     xlab="Time from cancer diagnosis (months)",
     ylab="Survival")
lty <- lifetab2(Surv(floor(surv_yy),death_cancer)~1, data = localised)
ltm <- lifetab2(Surv(floor(surv_mm),death_cancer)~1, data = localised)
kmy <- summary(mfit_years, times=0:20)
lines((0:20)*12, lty$surv, col="green")
lines(0:251, ltm$surv, col="orange")
lines((0:20)*12, kmy$surv, col="red", type="s")
legend("topright",
       legend=c("KM (months)", "KM (years)", "Actuarial (months)", "Actuarial (years)"),
       lty=1,
       col=c("black", "red", "orange", "green"),
       bty="n")
