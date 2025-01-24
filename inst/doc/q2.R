## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-24, 2016-03-07
## Edited: Mark Clements, 2017-08-02
###############################################################################

###############################################################################
## Exercise 2
###############################################################################
## @knitr loadDependencies
library(biostat3) # loads the survival and muhaz packages
library(dplyr)    # for data manipulation
library(bshazard) # for smoothed hazards
library(knitr)    # for html tables
library(tinyplot) # plt() for nicer base graphics
## utility function
as.data.frame.bshazard <- function(x, ...)
    with(x, data.frame(time,hazard,lower.ci,upper.ci))

## @knitr loadPreprocess
## Create 0/1 outcome variables and then show first six rows of the dataset
melanoma <- 
    transform(biostat3::melanoma,
              death_cancer = ifelse(status == "Dead: cancer", 1, 0),
              death_all = ifelse(status %in% c("Dead: cancer", "Dead: other"), 1, 0))
head(melanoma) |> knitr::kable("html")

## @knitr a_tabulate
## Tabulate by stage
Freq <- xtabs(~stage, data=melanoma)
cbind(Freq, Prop=proportions(Freq)) |> kable("html")

## @knitr a_plotSurv
survfit(Surv(surv_mm/12, death_cancer) ~ stage, data = melanoma) |>
    plot(col=1:4,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr a_plotHazard
library(dplyr)
hazards <- group_by(melanoma, stage) |> 
    do(as.data.frame(bshazard(Surv(surv_mm/12, death_cancer)~1, data=., verbose=FALSE))) |>
    ungroup()
with(hazards, plt(hazard~time|stage, ymin=lower.ci, ymax=upper.ci,
                  type="ribbon",
                  xlab="Time since diagnosis (years)",
                  ylab="Hazard",
                  col=1:4, lty=1, xlim=c(0,20), ylim=0:1))
ggplot(hazards,aes(x=time,y=hazard,group=stage)) + geom_line(aes(col=stage)) +
    geom_ribbon(aes(ymin=lower.ci, ymax=upper.ci, fill=stage), alpha=0.3) +
    xlim(0,20) + ylim(0,1) + 
    xlab('Time since diagnosis (years)') + ylab('Hazard')

## @knitr b_crudeRates
survRate(Surv(surv_mm/12, death_cancer) ~ stage, data=melanoma) |> kable("html")

## @knitr b2_crudeRates
library(dplyr)
melanoma |>
    group_by(stage) |>
    summarise(D = sum(death_cancer), M = sum(surv_mm/12), Rate = D/M,
              CI_low = stats::poisson.test(D,M)$conf.int[1],
              CI_high = stats::poisson.test(D,M)$conf.int[2]) |>
    kable("html")

## @knitr c_crudeRates1000
survRate(Surv(surv_mm/12/1000, death_cancer) ~ stage, data=melanoma) |>
    kable("html")

## @knitr d_crudeRates1000_sex
survRate(Surv(surv_mm/12/1000, death_cancer) ~ sex, data=melanoma) |>
    kable("html")

## @knitr d_plotSurv_sex
survfit(Surv(surv_mm/12, death_cancer) ~ sex, data = melanoma) |>
    plot(col=1:2,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival")
legend("topright", legend=levels(melanoma$sex), lty=1, col=1:2)

library(dplyr)
hazards <- group_by(melanoma, sex) |> 
    do(as.data.frame(bshazard(Surv(surv_mm/12, death_cancer)~1, data=., verbose=FALSE))) |>
    ungroup()
with(hazards, plt(hazard~time|sex, ymin=lower.ci, ymax=upper.ci,
                  type="ribbon",
                  xlab="Time since diagnosis (years)",
                  ylab="Hazard",
                  col=1:2, lty=1, xlim=c(0,20)))

## @knitr e_tabByAge
Freq <- xtabs(~status+agegrp, melanoma)
Freq
round(proportions(Freq,"agegrp")*100,1) # proportions for status by agegrp
chisq.test(Freq[-4,])

## @knitr f_survStage
par(mfrow=c(1, 1))
survfit(Surv(surv_mm/12, death_all) ~ stage, data = melanoma) |>
    plot(col=1:4,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates\nAll-cause")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr g_allCa75p
par(mfrow=c(1, 2))
survfit(Surv(surv_mm/12, death_cancer) ~ stage, data = subset(melanoma,agegrp=="75+")) |>
    plot(col=1:4,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates\nCancer | Age 75+")
survfit(Surv(surv_mm/12, death_all) ~ stage, data = subset(melanoma,agegrp=="75+")) |>
    plot(col=1:4,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier survival estimates\nAll-cause | Age 75+")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr h_allCaAgeGrp
par(mfrow=c(1, 2))
survfit(Surv(surv_mm/12, death_cancer) ~ agegrp, data = melanoma) |>
    plot(col=1:4,
         xlab = "Time since diagnosis (years)",
         ylab = "Survival",
         main = "Kaplan-Meier estimates of\ncancer survival by age group")
survfit(Surv(surv_mm/12, death_all) ~ agegrp, data = melanoma) |>
    plot(col=1:4,
         xlab = "Time since diagnois",
         ylab = "Survival",
         main = "Kaplan-Meier estimates of\nall-cause survival by age group")
legend("topright", levels(melanoma$agegrp), col=1:4, lty = 1)
