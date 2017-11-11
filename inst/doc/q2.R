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

## @knitr loadPreprocess
## Examine the data
data(melanoma)
head(melanoma)

## Create 0/1 outcome variables
melanoma <- 
    transform(melanoma,
              death_cancer = ifelse( status == "Dead: cancer", 1, 0),
              death_all = ifelse( status == "Dead: cancer" |
                                  status == "Dead: other", 1, 0))

## @knitr a_tabulate
## Tabulate by stage
Freq <- xtabs(~stage, data=melanoma)
cbind(Freq, Prop=prop.table(Freq))

## @knitr a_plotSurv
par(mfrow=c(1, 2))
mfit <- survfit(Surv(surv_mm, death_cancer) ~ stage, data = melanoma)

plot(mfit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival")
## legend("topright", levels(melanoma$stage), col=1:4, lty = 1)

hazards <- muhaz2(Surv(surv_mm, death_cancer)~stage, melanoma)
plot(hazards,
     col=1:4, lty=1, xlim=c(0,250), ylim=c(0,0.08),
     legend.args=list(bty="n"))

## @knitr b_crudeRates
survRate(Surv(surv_mm/12, death_cancer) ~ stage, data=melanoma)
melanoma %>%
    select(death_cancer, surv_mm, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), M = sum(surv_mm/12), Rate = D/M,
              CI_low = stats::poisson.test(D,M)$conf.int[1],
              CI_high = stats::poisson.test(D,M)$conf.int[2]) 

## @knitr c_crudeRates1000
survRate(Surv(surv_mm/12/1000, death_cancer) ~ stage, data=melanoma)

## @knitr d_crudeRates1000_sex
survRate(Surv(surv_mm/12/1000, death_cancer) ~ sex, data=melanoma)

## @knitr d_plotSurv_sex
par(mfrow=c(1, 2))
sfit <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)

plot(sfit, col=1:2,
     xlab = "Follow-up Time",
     ylab = "Survival")

plot(muhaz2(Surv(surv_mm,death_cancer)~sex, data = melanoma), 
     col = 1:2, lty = 1)

## @knitr e_tabByAge
xtabs(~status+agegrp, melanoma)

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
plot(afitage, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier estimates of\nall-cause survival by age group")
legend("topright", levels(melanoma$agegrp), col=1:4, lty = 1)
