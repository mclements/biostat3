# Purpose: To do the solution for Biostat III exercises in R
# Author: Annika Tillander, 2014-01-30
# Edited: Andreas Karlsson, 2015-02-24, 2016-03-07
###############################################################################

###############################################################################
## Exercise 2
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit
require(muhaz)    # for hazard estimates
require(dplyr)    # for data manipulation

## @knitr loadPreprocess
## Get the data for exercise 2
melanoma_raw <- read.dta("http://biostat3.net/download/melanoma.dta")

## Examine the data
head(melanoma_raw)

## Create 0/1 outcome variable
melanoma <- melanoma_raw %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0),
           death_all = ifelse( status == "Dead: cancer" |
                               status == "Dead: other", 1, 0))

## @knitr a_tabulate
## Tabulate by stage
melanoma %>%
    group_by(stage) %>%
    summarise(Freq = n(), Percent = n() / nrow(.)) %>%
    mutate(Cum = cumsum(Percent))

## @knitr a_plotSurv
par(mfrow=c(1, 2))
mfit <- survfit(Surv(surv_mm, death_cancer) ~ stage, data = melanoma)

plot(mfit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

hazByGroup <- function(group){
    with(subset(melanoma, group),
         muhaz(times=surv_mm,
               delta=death_cancer,
               min.time = min(surv_mm),
               max.time = max(surv_mm)))
}

plot(hazByGroup(melanoma$stage=="Unknown"), col=1, xlim=c(0,250), ylim=c(0,0.08))
lines(hazByGroup(melanoma$stage=="Localised"), col=2)
lines(hazByGroup(melanoma$stage=="Regional"), col=3)
lines(hazByGroup(melanoma$stage=="Distant"), col=4)
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr b_crudeRates
melanoma %>%
    select(death_cancer, surv_mm, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), M = sum(surv_mm), Rate = D/M) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

melanoma %>%
    select(death_cancer, surv_yy, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy), Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr c_crudeRates1000
melanoma %>%
    select(death_cancer, surv_yy, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr d_crudeRates1000_sex
melanoma %>%
    select(death_cancer, surv_yy, sex) %>%
    group_by(sex) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr d_plotSurv_sex
par(mfrow=c(1, 2))
sfit <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)

plot(sfit, col=1:2,
     xlab = "Follow-up Time",
     ylab = "Survival")
legend("bottomleft", levels(melanoma$sex), col=1:2, lty = 1)

plot(hazByGroup(melanoma$sex=="Male"), col=1, xlim=c(0,250))
lines(hazByGroup(melanoma$sex=="Female"), col=2)
legend("topright", levels(melanoma$sex), col=1:2, lty = 1)

## @knitr e_tabByAge
with(melanoma,table(status, agegrp))

## @knitr f_survStage
par(mfrow=c(1, 1))
afit <- survfit(Surv(surv_mm, death_all) ~ stage, data = melanoma)
plot(afit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nAll-cause")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr g_allCa75p
par(mfrow=c(1, 2))
mfit75 <- survfit(Surv(surv_mm, death_cancer) ~ stage, data = subset(melanoma,agegrp=="75+"))
plot(mfit75, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nCancer | Age 75+")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

afit75 <- survfit(Surv(surv_mm, death_all) ~ stage, data = subset(melanoma,agegrp=="75+"))
plot(afit75, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nAll-cause | Age 75+")
legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

## @knitr h_allCaAgeGrp
par(mfrow=c(1, 2))
mfitage <- survfit(Surv(surv_mm, death_cancer) ~ agegrp, data = melanoma)
plot(mfitage, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier estimates of\ncancer survival by age group")
legend("topright", levels(melanoma$agegrp), col=1:4, lty = 1)

afitage <- survfit(Surv(surv_mm, death_all) ~ agegrp, data = melanoma)
afit75 <- survfit(Surv(surv_mm, death_all) ~ stage, data = subset(melanoma,agegrp=="75+"))
plot(afitage, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier estimates of\nall-cause survival by age group")
legend("topright", levels(melanoma$agegrp), col=1:4, lty = 1)
