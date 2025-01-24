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
library(knitr) # kable

## @knitr loadPreprocess
localised <- subset(biostat3::melanoma, stage=="Localised") |>
    transform(death_cancer = ifelse(status == "Dead: cancer", 1, 0))

## @knitr actuarialYears
lifetab2(Surv(floor(surv_yy),death_cancer)~1, data = localised)[,1:7] |>
    kable("html")

## @knitr actuarialMonths
lifetab2(Surv(floor(surv_mm),death_cancer)~1, data = localised)[110:130,1:7] |>
    kable("html")

## @knitr kmYears
mfit_years <- survfit(Surv(surv_yy, death_cancer) ~ 1, data = localised)
summary(mfit_years)

## @knitr kmMonths
mfit_months <- survfit(Surv(surv_mm/12, death_cancer) ~ 1, data = localised)
summary(mfit_months,times=(110:130)/12)

## @knitr comparisonPlot
plot(mfit_months, conf.int=FALSE,
     ylim=c(0.7,1), cex=5,
     xlab="Time from cancer diagnosis (years)",
     ylab="Survival")
survfit(Surv(surv_yy, death_cancer) ~ 1, data = localised) |>
    lines(col="red", conf.int=FALSE)
lifetab2(Surv(surv_mm/12,death_cancer)~1, data = localised) |>
    lines(col="orange")
lifetab2(Surv(floor(surv_yy),death_cancer)~1, data = localised) |>
    lines(col="green")
legend("topright",
       legend=c("KM (months)", "KM (years)", "Actuarial (months)", "Actuarial (years)"),
       lty=1,
       col=c("black", "red", "orange", "green"),
       bty="n")
