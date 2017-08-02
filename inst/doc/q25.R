## Purpose: To do the solution for Biostat III exercises in R
## Author: Benedicte Delcoigne, 2014-03-05
###############################################################################

###############################################################################
## Exercise 25
###############################################################################

## In order to load data from Stata you need package "foreign"
## In order to read data from Stata 13 you need package "foreign"

## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # cox and conditional logistic analyses
require(Epi)      # sample a nested case-control study from a cohort

## @knitr loadPreprocess
## Get the data for exercise 25 and have a look at it
mel <- read.dta("http://biostat3.net/download/melanoma.dta")
head(mel)
mel <- subset(mel , stage=="Localised")               # restrict the cohort to stage==1
num_dupli  <- length(mel$id) - length(unique(mel$id)) # check if one line is one individual
num_dupli                                             # will be 0  if true
str(mel$status)                                       # check the structure of "status"
table (mel$status, useNA="always")                    # check the categories of "status"
mel$dc <- 0                                           # create an indicator for the event "Dead: cancer"
mel$dc[mel$status=="Dead: cancer"] <- 1               # dc == 1 if "Dead: Cancer" ; otherwise 0
table (mel$dc, useNA="always")                        # check
mel$surv_10y  <- ifelse( mel$surv_mm < 120 , mel$surv_mm, 120)   # restrict follow-up to 120 months
dim(mel[mel$surv_mm >= 120 & mel$dc==1,])             # how many cases after 120 months ?
mel[mel$surv_mm >= 120 & mel$dc==1,]$dc  <- 0         # events occuring after 120 months are changed to status = 0
table(mel$dc , mel$status)                            # check



## @knitr ex_25_coxph
str(mel$sex) ; str(mel$year8594) ; str(mel$agegrp)    #Check structure of risk factors/confounders
out_coh <- coxph(Surv(surv_10y,dc) ~ sex + year8594 + agegrp, data = mel, ties="breslow" )
summary(out_coh)



## @knitr n_ind
n_ind <- length(mel$id)
n_ind

## @knitr n_event
table(mel$dc, useNA="always")
ncase <-  table (mel$dc, useNA="always")[2]
ncase


## @knitr gen_ncc
nccdata <-ccwc( entry=0, exit=surv_10y , fail=dc, origin=0, controls=1, include=list(sex,year8594,agegrp,dc,id), data=mel )
tail(nccdata, 8)


## @knitr clogit
out_ncc <- clogit(Fail ~ sex + year8594 + agegrp + strata(Set), data=nccdata)
summary(out_ncc)


## @knitr n_unique_ncc
n_uni <- length(unique(nccdata$id))
n_uni

## @knitr compare_coh_ncc
comp <- data.frame(summary(out_coh)$coef[c(1:5),2],summary(out_ncc)$coef[c(1:5),2])
colnames(comp) <- c("cohort HR", "NCC HR")
comp                           # print the HR estimates from the cohort and the NCC
var <- data.frame((summary(out_coh)$coef[c(1:5),3])^2,(summary(out_ncc)$coef[c(1:5),3])^2,(summary(out_coh)$coef[c(1:5),3])^2 / (summary(out_ncc)$coef[c(1:5),3])^2 )
colnames(var) <- c("cohort var", "NCC var", "ratio coh/ncc" )
var                            # print the variances of estimates from the cohort and the NCC and their ratio


## @knitr loop_ncc
M <- 10                   # Number of loops: change the M value to change the number of loops
param   <- matrix(0,M,5)  # Define the matrice of the coefficients
for (i in 1:M)  {         # Start of the loop, create NCC data and analyse it
    nccdata <-ccwc( entry=0, exit=surv_10y , fail=dc, origin=0, controls=1, include=list(sex,year8594,agegrp), data=mel )
    out_ncc <- clogit(Fail ~ sex + year8594 + agegrp + strata(Set), data=nccdata)
    param [ i , 1:5] <- coef(out_ncc)   # store the 5 coefficients M times
}                  # End of the loop

## @knitr output
param <- as.data.frame(param)        # data frame with the coefficients
param <- exp(param)                  # data frame with the HR
param <- rbind (summary(out_coh)$coef[c(1:5),2], param)
colnames(param) <- c("sex Female","year_8594","agegrp45-59","agegrp60-74","agegrp75+")
rownames(param) <- c("cohort", c(1:M))
param
mean_param <- apply(param[c(2:M),], 2, mean)  # compute the mean of the HR for the M loops
sd_param <-  apply(param[c(2:M),], 2, sd)     # compute the sd of the HR for the M loops
par_sum <- rbind (summary(out_coh)$coef[c(1:5),2], mean_param, sd_param)
rownames(par_sum) <- c("cohort HR","mean HR NCC","sd HR NCC")
colnames(par_sum) <- c("sex Female","year_8594","agegrp45-59","agegrp60-74","agegrp75+")
par_sum

par(mfrow=c(2,3) )                                   # allow 2*3 graphs on the same page
hist(param[,1], main="sex Female", xlab="HR value")  # histogram of HR for variable sex (female vs male)
abline(v=summary(out_coh)$coef[1,2], col="green")
abline(v=mean_param[1], col="red")
hist(param[,2], main="year_8594", xlab="HR value")   # histogram of HR for variable period (later vs earlier)
abline(v=summary(out_coh)$coef[2,2], col="green")
abline(v=mean_param[2], col="red")
hist(param[,3], main="agegrp45-59", xlab="HR value") # histogram of HR for variable agegrp (45-59 vs reference)
abline(v=summary(out_coh)$coef[3,2], col="green")
abline(v=mean_param[3], col="red")
hist(param[,4], main="agegrp60-74", xlab="HR value") # histogram of HR for variable agegrp (60-74 vs reference)
abline(v=summary(out_coh)$coef[4,2], col="green")
abline(v=mean_param[4], col="red")
hist(param[,5], main="agegrp75+", xlab="HR value")   # histogram of HR for variable agegrp (75+ vs reference)
abline(v=summary(out_coh)$coef[5,2], col="green")
abline(v=mean_param[5], col="red")
                                        # obs: the vertical green line is at the cohort's HR value, the red line at the mean of the M loops
