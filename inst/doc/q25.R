## Purpose: To do the solution for Biostat III exercises in R
## Author: Benedicte Delcoigne, 2014-03-05
## Revised: Mark Clements, 2017-11-03
###############################################################################

###############################################################################
## Exercise 25
###############################################################################

## @knitr loadDependencies
library(biostat3) # cox and conditional logistic analyses
library(Epi)      # sample a nested case-control study from a cohort

## @knitr loadPreprocess
## Get the data for exercise 25 and have a look at it
data(melanoma)
mel <- subset(melanoma, stage=="Localised")           # restrict the cohort to stage==1
mel <- transform(mel,
                 dc = (mel$status=="Dead: cancer" & surv_mm<120)+0,
                 surv_10y = pmin(120, surv_mm))
table(mel$dc , mel$status)                            # check


## @knitr ex_25_coxph
str(mel$sex) ; str(mel$year8594) ; str(mel$agegrp)    #Check structure of risk factors/confounders
out_coh <- coxph(Surv(surv_10y,dc) ~ sex + year8594 + agegrp, data = mel)
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
comp <- exp(data.frame(coef(out_coh),coef(out_ncc)))
colnames(comp) <- c("cohort HR", "NCC HR")
comp                           # print the HR estimates from the cohort and the NCC
var <- data.frame(diag(vcov(out_coh)), diag(vcov(out_ncc)), diag(vcov(out_coh))/diag(vcov(out_ncc)))
colnames(var) <- c("cohort var", "NCC var", "ratio coh/ncc" )
var                            # print the variances of estimates from the cohort and the NCC and their ratio


## @knitr loop_ncc
M <- 20                   # Number of loops: change the M value to change the number of loops
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
colnames(param) <- names(coef(out_coh))
rownames(param) <- c("cohort", c(1:M))
param
mean_param <- apply(param[c(2:M),], 2, mean)  # compute the mean of the HR for the M loops
sd_param <-  apply(param[c(2:M),], 2, sd)     # compute the sd of the HR for the M loops
par_sum <- rbind (summary(out_coh)$coef[c(1:5),2], mean_param, sd_param)
rownames(par_sum) <- c("cohort HR","mean HR NCC","sd HR NCC")
colnames(par_sum) <- names(coef(out_coh))
par_sum

par(mfrow=c(2,3) )                                   # allow 2*3 graphs on the same page
for (i in 1:5) {
    hist(param[,i], main=colnames(par_sum)[i], xlab="HR value")  # histogram of HR for variable sex (female vs male)
    abline(v=summary(out_coh)$coef[i,2], col="green")
    abline(v=mean_param[i], col="red")
}
plot.new()
legend(0, 0.5, c("Cox HR value","Mean NCC"),lty=1,col=c("green","red"),bty="n")
