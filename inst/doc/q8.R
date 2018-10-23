## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01
##         Benedicte Delcoigne, 2015-03-03
##         Mark Clements, 2017-11-05
###############################################################################


###############################################################################
## Exercise 8
###############################################################################
## @knitr loadDependencies
library(biostat3) 
library(dplyr)     # for data manipulation
library(bshazard)  # for bshazard

## @knitr loadPreprocess
diet <- biostat3::diet
head(diet)
summary(diet)

## @knitr 8a1_Haz_att_age
scale <- 365.24
plot.bshazard(bshazard(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ 1, data=subset(diet,hieng=="low")),
                         ylim=c(0,0.03), conf.int=FALSE, xlab="Attained age (years)")
lines(bshazard(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ 1, data=subset(diet,hieng=="high")),
      col="red", conf.int=FALSE)
legend("topleft", legend=c('hieng=="low"','hieng=="high"'), col=1:2, lty=1, bty="n")

## @knitr 8a2_Haz_time_entry
scale <- 365.24
plot(muhaz2(Surv((dox - doe)/scale, chd) ~ hieng, data=diet), lty=2,
     xlab="Time since study entry (years)", ylim=c(0,0.025),
     legend=FALSE)
lines(bshazard(Surv((dox - doe)/scale, chd) ~ 1, data=subset(diet,hieng=="low")),
      conf.int=FALSE)
lines(bshazard(Surv((dox - doe)/scale, chd) ~ 1, data=subset(diet,hieng=="high")),
      col="red", conf.int=FALSE)
legend("topleft", legend=c('hieng=="low" (bshazard)','hieng=="high" (bshazard)',
                           'hieng=="low" (muhaz)','hieng=="high" (muhaz)'), col=1:2, lty=c(1,1,2,2), bty="n")


## @knitr 8b_ir
diet <- mutate(diet, y1k = y / 1000)
poisson8b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8b)
eform(poisson8b)


## @knitr 8c_ir
## Create BMI variable
diet$bmi <- diet$weight/((diet$height/100)^2)

## Create orderly varable instead of categorical, start at zero
levels(diet$job)
diet <- mutate(diet, jobNumber = unclass(job) - 1)

poisson8c <- glm( chd ~ hieng + jobNumber + bmi + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8c)
eform(poisson8c)

##########################################################

## @knitr 8d_agebands

## Split time at 30,50,60 and 72 with time scale age at entry to attained age
scale <- 365.24
age.cuts <- c(30,50,60,72)
diet.spl.dob <- survSplit(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ ., data=diet, cut=age.cuts,start="tstart",end="tstop")


## Tabulate ageband
diet.spl.dob %>% select(id, tstart, tstop, y) %>% filter(id<=3) %>% arrange(id, tstart)

## Create an agespan variable
diet.spl.dob <- mutate(diet.spl.dob,
                       agespan = cut(tstop, age.cuts))

## Make the numeric variables factors since we want to model them with dummie variables and calculate time at risk
diet.spl.dob <- mutate(diet.spl.dob,
                       jobNumber = as.factor(jobNumber),
                       risk_time = (tstop-tstart))

## Tabulate ageband including risk_time
diet.spl.dob %>% select(id,  tstart, tstop, y,risk_time) %>% filter(id<=3) %>% arrange(id, tstart)


## Tabulate number of events per agespan
xtabs(~agespan+chd,diet.spl.dob)

## @knitr 8d_model1
poisson8d <- glm( chd ~ hieng + agespan + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.dob)

summary(poisson8d)
eform(poisson8d)

## @knitr 8d_model2
poisson8d <- glm( chd ~ hieng + agespan + jobNumber + bmi + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.dob)

summary(poisson8d)
eform(poisson8d)

## @knitr 8.e.i
time.cuts <- c(0, 5, 10, 15, 22)
diet.spl.t_entry <- survSplit(Surv((dox-doe)/365.24, chd) ~ ., data=diet, cut=time.cuts, end="tstop", start="tstart", event="chd")

##Tabulate ageband
diet.spl.t_entry %>% select(id, tstart, tstop, y) %>% filter(id<=3) %>% arrange(id, tstart)

diet.spl.t_entry <- mutate(diet.spl.t_entry,
                           fu = cut(tstop, time.cuts),
                           risk_time = (tstop-tstart))

##Tabulate ageband including risk_time
diet.spl.t_entry %>% select(id, fu, tstart, y, risk_time) %>% filter(id<=3) %>% arrange(id, tstart)

poisson8e1 <- glm( chd ~ fu + hieng + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary(poisson8e1)
eform(poisson8e1)

## @knitr 8.e.ii
poisson8e2 <- glm( chd ~ fu + hieng + job + bmi + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary(poisson8e2)
eform(poisson8e2)
