## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28, 2016-03-08
###############################################################################

###############################################################################
## Exercise 3
###############################################################################
## @knitr loadDependecies
require(foreign)  # needed to read data set from Stata
require(survival) # for Surv and survfit
require(muhaz)    # for hazard estimates
require(dplyr)    # for data manipulation

## @knitr loadPreprocess
melanoma_raw <- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma <- melanoma_raw %>%
    filter(stage=="Localised") %>% # subset those with localised cancer
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0),
           death_all = ifelse( status == "Dead: cancer" |
                               status == "Dead: other", 1, 0))

## @knitr a_survDiaDate
mfityear8594 <- survfit(Surv(surv_mm, death_cancer==1) ~ year8594, data = melanoma)

plot(mfityear8594, col = 1:2,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("topright", c("Diagnosed 75-84", "Diagnosed 85-94"), col=1:2, lty = 1)

## @knitr b_hazDiaDate
hazByGroup <- function(group){
    with(subset(melanoma, group),
         muhaz(times=surv_mm,
               delta=death_cancer,
               min.time = min(surv_mm),
               max.time = max(surv_mm)))
}

plot(hazByGroup(melanoma$year8594 == "Diagnosed 75-84"), col=1, main="Smoothed hazard estimates")
lines(hazByGroup(melanoma$year8594 == "Diagnosed 85-94"), col=2)
legend("topright", c("Diagnosed 75-84", "Diagnosed 85-94"), col=1:2, lty = 1)

## @knitr c_testDiaDate
## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma)
## Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma, rho=1)

## @knitr d_crudeRates1000_agegrp
melanoma %>%
    select(death_cancer, surv_mm, agegrp) %>%
    group_by(agegrp) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_mm)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

mfit_agegrp <- survfit(Surv(surv_mm, death_cancer==1) ~ agegrp, data = melanoma)
plot(mfit_agegrp, col = 1:4,
     xlab = "Months since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", c("0-44", "45-59", "60-74", "75+"), col=1:4, lty = 1)

## @knitr e_crudeRates1000_agegrp
mfit_agegrp_year <- survfit(Surv(surv_mm/12, death_cancer==1) ~ agegrp, data = melanoma)
plot(mfit_agegrp_year, col = 1:4,
     xlab = "Months since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", c("0-44", "45-59", "60-74", "75+"), col=1:4, lty = 1)

melanoma %>%
    select(death_cancer, surv_mm, agegrp) %>%
    group_by(agegrp) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_mm)/12/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr f_sexDiff
par(mfrow=c(1, 2))
mfit_sex <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)

plot(mfit_sex, col = 1:2,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("topright", c("Male", "Female"), col=1:2, lty = 1)

plot(hazByGroup(melanoma$sex == "Male"), col=1, main="Smoothed hazard estimates")
lines(hazByGroup(melanoma$sex == "Female"), col=2)
legend("topright", c("Male", "Female"), col=1:2, lty = 1)

## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma)
