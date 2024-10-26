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
library(dplyr)    # for data manipulation
library(survminer) # for ggsurvplot

## @knitr loadPreprocess
melanoma <- biostat3::melanoma %>%
    filter(stage=="Localised") %>% # subset those with localised cancer
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0),
           death_all = ifelse( status == "Dead: cancer" |
                               status == "Dead: other", 1, 0))

## @knitr a_survDiaDate
mfityear8594 <- survfit(Surv(surv_mm/12, death_cancer==1) ~ year8594, data = melanoma)

plot(mfityear8594, col = 1:2,
     xlab = "Time since diagnosis (years)",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$year8594), col=1:2, lty = 1)

## @knitr b_hazDiaDate
plot(muhaz2(Surv(surv_mm,death_cancer)~year8594, data=melanoma))

ggplot.muhazList(muhaz2(Surv(surv_mm/12,death_cancer)~year8594, data=melanoma)) +
    xlab("Time since diagnosis (years)")


## @knitr c_testDiaDate
## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma)
## Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer) ~ year8594, data=melanoma, rho=1)

## @knitr d_crudeRates1000_agegrp
survRate(Surv(surv_mm/1000,death_cancer)~agegrp, data=melanoma)
## melanoma %>%
##     select(death_cancer, surv_mm, year8594) %>%
##     group_by(year8594) %>%
##     summarise(D = sum(death_cancer), Y = sum(surv_mm)/1000, Rate = D/Y,
##               CI_low = stats::poisson.test(D,Y)$conf.int[1],
##               CI_high = stats::poisson.test(D,Y)$conf.int[2]) 

mfit_agegrp <- survfit(Surv(surv_mm, death_cancer) ~ agegrp, data = melanoma)
plot(mfit_agegrp, col = 1:4,
     xlab = "Months since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$agegrp), col=1:4, lty = 1)

## @knitr e_crudeRates1000_agegrp
mfit_agegrp_year <- survfit(Surv(surv_mm/12, death_cancer) ~ agegrp, data = melanoma)
plot(mfit_agegrp_year, col = 1:4,
     xlab = "Years since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", levels(melanoma$agegrp), col=1:4, lty = 1)


survRate(Surv(surv_mm/12/1000,death_cancer)~year8594, data=melanoma)

## @knitr f_sexDiff
par(mfrow=c(1, 2))
mfit_sex <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)
plot(mfit_sex, col = 1:2,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
plot(muhaz2(Surv(surv_mm,death_cancer)~sex, data=melanoma))

## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma)

