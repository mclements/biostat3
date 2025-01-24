## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28, 2016-03-08
## Edited: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 3
###############################################################################
## @knitr loadDependencies
library(biostat3)
library(survival)
library(dplyr)    # for data manipulation
library(survminer) # for ggsurvplot
library(bshazard) # smooth hazards

## @knitr loadPreprocess
melanoma <- biostat3::melanoma |>
    subset(stage=="Localised") |>
    transform(death_cancer = ifelse(status == "Dead: cancer", 1, 0),
              death_all = ifelse(status == "Dead: cancer" |
                                 status == "Dead: other", 1, 0))

## @knitr a_survDiaDate
survfit(Surv(surv_mm/12, death_cancer==1) ~ year8594, data = melanoma) |>
    plot(col = 1:2,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$year8594), col=1:2, lty = 1)

## @knitr b_hazDiaDate
par(mfrow=1:2)
for(level in levels(melanoma$year8594))
    plot(bshazard(Surv(surv_mm/12,death_cancer)~1, subset(melanoma,year8594==level),
                  verbose=FALSE), main=level, xlim=c(0,20), ylim=c(0,0.1),
         xlab="Time since diagnosis\n(years)")

## @knitr b_hazDiaDate_ggplot
group_by(melanoma, year8594) |> 
    do(as.data.frame(bshazard(Surv(surv_mm/12, death_cancer)~1, data=., verbose=FALSE))) |>
    ungroup() |>
    ggplot(aes(x=time,y=hazard,group=year8594)) + geom_line(aes(col=year8594)) +
    geom_ribbon(aes(ymin=lower.ci, ymax=upper.ci, fill=year8594), alpha=0.3) +
    xlim(0,20) + ylim(0,0.1) + 
    xlab('Time since diagnosis (years)') + ylab('Hazard')

## @knitr c_testDiaDate
## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma)
## Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma, rho=1)

## @knitr d_crudeRates1000_agegrp
survRate(Surv(surv_mm/1000,death_cancer)~agegrp, data=melanoma)
## melanoma |>
##     select(death_cancer, surv_mm, year8594) |>
##     group_by(year8594) |>
##     summarise(D = sum(death_cancer), Y = sum(surv_mm)/1000, Rate = D/Y,
##               lower.ci = stats::poisson.test(D,Y)$conf.int[1],
##               upper.ci = stats::poisson.test(D,Y)$conf.int[2]) 

survfit(Surv(surv_mm, death_cancer) ~ agegrp, data = melanoma) |>
    plot(col = 1:4,
         xlab = "Months since diagnosis",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$agegrp), col=1:4, lty = 1)

## @knitr e_crudeRates1000_agegrp
survfit(Surv(surv_mm/12, death_cancer) ~ agegrp, data = melanoma) |>
    plot(col = 1:4,
         xlab = "Years since diagnosis",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$agegrp), col=1:4, lty = 1)


survRate(Surv(surv_mm/12/1000,death_cancer)~year8594, data=melanoma)

## @knitr f_sexDiff
aplot = function() {
    survfit(Surv(surv_mm/12, death_cancer) ~ sex, data = melanoma) |>
    plot(col = 1:2,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates")
}
if (requireNamespace("muhaz", quietly=TRUE)) {
    par(mfrow=c(1, 2))
    aplot()
    plot(muhaz2(Surv(surv_mm/12,death_cancer)~sex, data=melanoma), lty=1,
         xlab="Time since diagnosis (years)")
} else {
    par(mfrow=c(1, 1))
    aplot()
    legend("bottomleft", legend=levels(melanoma$sex), lty=1, col=1:2)
}

## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma)

