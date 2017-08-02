## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01
##         Benedicte Delcoigne, 2015-03-03
###############################################################################


###############################################################################
## Exercise 8
###############################################################################
## @knitr loadDependecies
require(survival) # for Surv and survfit
require(foreign)  # for reading data set from Stata
require(ggplot2)
require(dplyr)    # for data manipulation
require(muhaz)    # for hazard estimates


###########################################
### A help function to calculate ###
### and print incidence (hazard) ratios
### from a fitted poisson regression
### from glm
###########################################

IRR <- function(fit){
    summfit <- summary(fit )$coefficients
    IRfit <- exp( cbind( summfit[, 1:2], summfit[, 1] - 1.96*summfit[, 2], summfit[, 1] +
                        1.96*summfit[, 2] ) )
    colnames(IRfit) <- c("IRR", "Std. err", "CI_lower", "CI_upper")
    print(IRfit)
}


## @knitr loadPreprocess
diet <- read.dta("http://biostat3.net/download/diet.dta")
head(diet)
summary(diet)

## @knitr 8a1_first_step_att_age
## Creating a variable for attained age
diet <- diet %>% mutate(att_age = as.numeric(dox - dob) / 365.24)
summary(diet$att_age)


## @knitr 8a1_Haz_att_age
## Run muhaz smoother for Kaplan-Meier hazards ones per strata, note that this is left censured data and the muhaz smoother is not made for this.
dietHiengAge <- diet %>%  group_by(hieng) %>%
    do( h = muhaz(times = .$att_age, delta = .$chd, min.time = min(.$att_age[as.logical(.$chd)]), max.time = max(.$att_age[as.logical(.$chd)])-2)) %>%
    do( data.frame(Hazard = .$h$haz.est, Time = .$h$est.grid, Strata = .$hieng))
ggplot(dietHiengAge, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")

## Have a look at the raw hazard if you dare...
## hazard <- data.frame(with( diet, kphaz.fit(time=att_age, status=chd, strata=hieng)))
## ggplot(hazard, aes(time, haz, group=strata, color=factor(strata))) + geom_point()

## @knitr 8a2_Haz_time_entry
diet <- diet %>% mutate(t_entry = as.numeric(dox - doe) / 365.24)
summary(diet$t_entry)

dietHiengEntry <- diet %>%  group_by(hieng) %>%
    do( h = muhaz(times = .$t_entry, delta = .$chd,
            min.time = min(.$t_entry[as.logical(.$chd)]), max.time = max(.$t_entry[as.logical(.$chd)]))) %>%
    do( data.frame(Hazard = .$h$haz.est, Time = .$h$est.grid, Strata = .$hieng))
ggplot(dietHiengEntry, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")


## @knitr 8b_ir
diet <- mutate(diet, y1k = y / 1000)
poisson8b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8b)
exp(cbind(coef(poisson8b),confint(poisson8b)))


## @knitr 8c_ir
## Create BMI variable
diet$bmi <- diet$weight/((diet$height/100)^2)

## Create orderly varable instead of categorical, start at zero
diet <- diet %>% mutate(jobNumber = match(job, c("driver","conductor","bank")) -1)

poisson8c <- glm( chd ~ hieng + jobNumber + bmi + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8c)
exp(cbind(coef(poisson8c),confint(poisson8c)))

##########################################################

## @knitr 8d_agebands
## Create a variable for age at entry
diet <- mutate(diet, entry_age = as.numeric( doe - dob ) / 365.24)

## Split time at 30,50,60 and 72 with time scale age at entry to attained age, the zero=0 makes output be in absolute age
diet.spl.dob <- survSplit(diet, cut=c(30,50,60,72), end="att_age", start="entry_age", event="chd", zero = 0)

## Tabulate ageband
diet.spl.dob %>% select(id, entry_age, att_age, y) %>% filter(id<=3) %>% arrange(id, entry_age)

## Create an agespan variable
diet.spl.dob <- mutate(diet.spl.dob,
                       agespan = NA,
                       agespan = ifelse( entry_age>=30 & entry_age <50 , 30, agespan),
                       agespan = ifelse( entry_age>=50 & entry_age <60 , 50, agespan),
                       agespan = ifelse( entry_age>=60 & entry_age <72, 60, agespan))

## Make the numeric variables factors since we want to model them with dummie variables and calculate time at risk
diet.spl.dob <- mutate(diet.spl.dob,
                       agespan = as.factor(agespan),
                       jobNumber = as.factor(jobNumber),
                       risk_time = (att_age - entry_age))

## Tabulate ageband including risk_time
diet.spl.dob %>% select(id,  entry_age, att_age, y,risk_time) %>% filter(id<=3) %>% arrange(id, entry_age)


## Tabulate number of events per agespan
diet.spl.dob %>%
    group_by(agespan) %>%
    summarise(noFailure=sum(chd==0), Failure=sum(chd), Total=n())

## @knitr 8d_model1
poisson8d <- glm( chd ~ hieng + agespan + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.dob)

summary(poisson8d)
IRR(poisson8d)

## @knitr 8d_model2
poisson8d <- glm( chd ~ hieng + agespan + jobNumber + bmi + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.dob)

summary(poisson8d)
IRR(poisson8d)

## @knitr 8.e.i

diet.spl.t_entry <- survSplit(diet, cut=c(0, 5, 10, 15, 22), end="t_entry", start="start", event="chd")

##Tabulate ageband
diet.spl.t_entry %>% select(id, start, t_entry, y) %>% filter(id<=3) %>% arrange(id, t_entry)

diet.spl.t_entry <- mutate(diet.spl.t_entry,
                           fu = as.factor(start) ,
                           risk_time = (t_entry-start))

##Tabulate ageband including risk_time
diet.spl.t_entry %>% select(id, start, t_entry, y, risk_time) %>% filter(id<=3) %>% arrange(id, t_entry)

poisson8e1 <- glm( chd ~ fu + hieng + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary(poisson8e1)
IRR(poisson8e1)

## @knitr 8.e.ii
poisson8e2 <- glm( chd ~ fu + hieng + job + bmi + offset( log( t_entry) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary(poisson8e2)
IRR(poisson8e2)
