## Purpose: Solution for Biostat III exercises in R
## Author: Mark Clements, 2017-11-03
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("survival")
## install.packages("dplyr")


###############################################################################
## Exercise 23
###############################################################################
## @knitr loadDependencies
library(biostat3) # for Surv and survfit
library(dplyr)    # for data manipulation
library(tinyplot) # for some nice plots

## @knitr loadPreprocess
calculate_smr = function(data)
    summarise(data,
              Observed = sum(observed),
              Expected = sum(expected)) |>
        mutate(SMR=Observed/Expected,
               poisson.ci(Observed,Expected))

## @knitr 23.a1
mel <- filter(biostat3::melanoma, stage == "Localised") |> 
        mutate( dead = (status %in% c("Dead: cancer","Dead: other") & surv_mm <= 120)+0, 
                surv_mm = pmin(120, surv_mm))
head(mel) 

## @knitr 23.a2
## Define the age at start and end of follow-up 
mel <- mutate(mel, adx = age+0.5,   # age at diagnosis  (mid-point approximation) 
              astart = adx, 
              astop  = adx+surv_mm/12)
## Split by age 
mel.split <- survSplit(mel, cut = 1:105, event = "dead", 
                       start = "astart", end = "astop")
## Quick check: the first two ids 
subset(mel.split, id<=2, select = c(id, astart, astop, dead)) 

## @knitr 23.b1
# For each age time band from (a), we calculate the start and stop in calendar time 
# We calculate the time since diagnosis as difference between age at start/stop and 
# age at diagnosis, and add that interval to year at diagnosis
mel.split2 <- mutate(mel.split, 
                    ystart = ydx + astart - adx, 
                    ystop  = ydx + astop - adx)
subset(mel.split2, id<=2, select = c(id, adx, astart, astop, dead, ydx, ystart, ystop))

## @knitr 23.b2
## Now we can split along the calendar time 
## For each of the new age-calendar time bands, we now have to adjust the age at 
## start and end in the same way as above 
mel.split2 <- survSplit( mel.split2, cut = 1970:2000, event = "dead", 
                         start = "ystart", end = "ystop" ) |>
              mutate( astart = adx + ystart - ydx, 
                      astop  = adx + ystop - ydx)
## Quick check: this seems ok 
subset(mel.split2, id<=2, select = c(id, ystart, ystop, astart, astop, dead)) 

## @knitr 23.c1
## We calculate the total person time at risk for each time band
mel.split2 <- mutate( mel.split2, 
                      age  = floor(astart),  # Age at which person time was observed 
                      year = floor(ystart),  # Calendar year during which person time was observed 
                      pt  =  ystop - ystart  # ... or astop - astart, works the same
                    ) 
subset(mel.split2, id<=2, select = c(id, ystart, ystop, astart, astop, dead, age, year, pt))

## @knitr 23.c2
## Now tabulate: sum of person time across all combinations of age & year 
## (for some years, ages) 
xtabs(pt ~ age + year, data=mel.split2, subset = age>=50 & age<60 & year>=1980 & year<1990)

## @knitr 23.c3
## Same: sum up 0/1 alive/dead for total count of deaths 
xtabs(dead ~ age + year, data=mel.split2, subset = age>=50 & age<60 & year>=1980 & year<1990) 

## @knitr 23.d
mel.split2 <- mutate(mel.split2, 
                     age10  = cut(age, seq(0, 110 ,by=10), right=FALSE), 
                     year10 = cut(year, seq(1970, 2000, by=5), right=FALSE) 
                    ) 
sr <- survRate(Surv(pt, dead) ~ sex + age10 + year10, data=mel.split2) 
rownames(sr) <-1:nrow(sr) ## Simple rownames for display 
head(sr, n = 20) 

## @knitr 23.e1
pt <- mutate(mel.split2, sex = unclass(sex)) |>    # make sex integer to be in line with popmort 
    group_by(sex, age, year)                 |>    # aggregate by sex, age, year 
    summarise(pt = sum(pt), observed = sum(dead), .groups="keep") |> # sum the person time, deaths 
    ungroup()  # For convenience
head(pt) 

## @knitr 23.e2
head(rstpm2::popmort)
summary(rstpm2::popmort) 

## @knitr 23.e3
joint <- left_join(pt, rstpm2::popmort)
head(joint) 

## @knitr 23.e4
joint <- mutate(joint, expected = pt * rate) 
head(joint) 


## @knitr 23.f1
calculate_smr(joint) |> kable("html")

## @knitr 23.f2
group_by(joint, sex) |> calculate_smr() |> kable("html")

## @knitr 23.f3
SMR_byYear <- group_by(joint, year) |> calculate_smr()
SMR_byYear |> kable("html")
plot(SMR~year, data=SMR_byYear, type = "o")  # quick & dirty plot 
with(SMR_byYear,
     plt(SMR~year, type = "ribbon", ymin=`2.5 %`, ymax=`97.5 %`))

## @knitr 23.f4
joint <- mutate(joint, age_group = cut(age, seq(0, 110, by=10), right = FALSE))
SMR_byAge <- group_by(joint, age_group) |> calculate_smr()
SMR_byAge |> kable("html")
plot(SMR~age_group, SMR_byAge, xlab="Age group (years)")
abline( h = 1:2, lty = 2)  # two reference lines at 1 & 2

## @knitr 23.f5
by(joint, joint$sex, function(data) poisson.test(sum(data$observed), sum(data$expected)))


## @knitr 23.g1
joint2 <- transform(joint,
                    sex  = factor(sex, levels = 1:2, labels = c("m", "f")),
                    year =  factor(year) |> relevel(ref = "1985"),  # mid-study
                    age_group = relevel(age_group, ref = "[70,80)")  # Close to one already
                   )
## Model & parameters
summary(fit <- glm(observed ~ sex + year + age_group + offset(log(expected)), data=joint2, family=poisson)) 
eform(fit)

## @knitr 23.g2
drop1(fit, test = "Chisq")
