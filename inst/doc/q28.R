###############################################################################
## Exercise 28 (Modelling cause-specific mortality)
## Flexible parametric models (fpm) implemented in Stata (stpm2) and R (rstpm2)
## One difference: Restricted cubic splines
##                 using truncated power basis (stpm2, Stata)
##                 using B-spline basis (rstpm2, R)
## R codes suggested by Xing-Rong Liu and Mark Clements (2016-03-06)
## For more detail interpretation, please see solutions_biostat3_2016.pdf file
###############################################################################

## @knitr loadDependencies
library(biostat3)  
library(rstpm2)  # for the flexible parametric model
library(dplyr)   # for data manipulation

## @knitr loadPreprocess
data(melanoma)

## Extract subset and create 0/1 outcome variables
melanoma0 <- melanoma %>% filter(stage=="Localised") %>%
             mutate(event = ifelse(status=="Dead: cancer" & surv_mm<120, 1, 0),
                    time = pmin(120, surv_mm)/12,
                    agegrp1 = (agegrp=="0-44")+0,  # used by time-dependent effect
                    agegrp2 = (agegrp=="45-59")+0, # used by time-dependent effect
                    agegrp3 = (agegrp=="60-74")+0, # used by time-dependent effect
                    agegrp4 = (agegrp=="75+")+0)   # used by time-dependent effect

## @knitr a_flex
## (a) Flexible parametric model with df=4

fpma <- stpm2(Surv(time,event) ~ year8594, data=melanoma0, df=4)
exp(cbind(IRR=coef(fpma),confint(fpma)))
## updated package: eform(fpma)
summary(fpma)

## @knitr a_cox
cox <- coxph(Surv(time, event) ~ year8594,
             data=melanoma0) # year8594 is a categorical variable
summary(cox)

## @knitr b_surv
## (b) Prediction and plots of survival and hazard by calendar period
years <- levels(melanoma0$year8594)
alegend <- function() legend("topright", legend=years, col=1:2, lty=1, bty="n")

plot(fpma,newdata=data.frame(year8594=years[1]),
     xlab="Time since diagnosis (years)")
lines(fpma,newdata=data.frame(year8594=years[1]), ci=FALSE) # repeated
alegend()

## @knitr b_haz
plot(fpma,newdata=data.frame(year8594=years[1]), type="haz",
     xlab="Time since diagnosis (years)", ylab="Hazard")
lines(fpma,newdata=data.frame(year8594=years[2]), type="haz", col=2)
alegend()

## @knitr c_haz_log
## (c) hazards on log scale, adding log="y"
plot(fpma,newdata=data.frame(year8594=years[1]), type="haz",
     xlab="Time since diagnosis (years)", log="y", ci=FALSE,
     ylab="Hazard (log scale)")
lines(fpma,newdata=data.frame(year8594=years[2]), type="haz", col=2)
alegend()

## @knitr d_AIC_BIC
summary(fpma)
aicc <- bicc <- beta <- se <- rep(NULL,6)
BIC <- function(object, nknots){
    -2 * as.numeric(bbmle::logLik(object)) + log(sum(melanoma0$event)) * (1 + nknots)
}

for (i in 1:6 ) {
  fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
  aicc[i] <- AIC(fitaic)
  bicc[i] <- BIC(fitaic, i)
  beta[i] <- as.numeric(coef(fitaic)[2])
  se[i] <- as.numeric(sqrt(fitaic@vcov[2,2]))
}
rbind(beta=beta, se=se, aicc=aicc, bicc=bicc)

## @knitr e_base_surv
## Baseline survival
fitaic0 <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=6)
plot(fitaic0,newdata=data.frame(year8594=years[1]), lty=6, ci=F,
     xlab="Time since diagnosis (years)")

dfs <- c("s_df1","s_df2","s_df3","s_df4","s_df5","s_df6")
for (i in 1:5 ) {
  fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
  plot(fitaic,newdata=data.frame(year8594=years[1]), add=TRUE,lty=i,
       xlab="Time since diagnosis (years)")
}
legend("topright", legend=dfs[1:6], lty=1:6)

## Baseline hazard
fitaic1 <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=6)
plot(fitaic1,newdata=data.frame(year8594=years[1]), lty=6, type="haz",
     ci=F, xlab="Time since diagnosis (years)", ylab="Hazard")

for (i in 1:5 ) {
  fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
  plot(fitaic,newdata=data.frame(year8594=years[1]), add=TRUE, lty=i,
       type="haz", xlab="Time since diagnosis (years)", ylab="Hazard")
}
dfs2 <- c("h_df1","h_df2","h_df3","h_df4","h_df5","h_df6")
legend("topright", legend=dfs2[1:6], lty=1:6)

## @knitr f_sex_age
fpmf <- stpm2(Surv(time, event) ~ sex + year8594 + agegrp,
              data=melanoma0, df=4)

beta <- coef(fpmf)[2:6] ## log(HR)
se <- sqrt(diag(fpmf@vcov[2:6,2:6]))            # standard errors of beta
z <- beta/se                                    # z-scores
confin <- cbind(beta - 1.96*se, beta + 1.96*se) # 95% confidence interval of beta=log(HR)
p_value <- 1-pchisq((beta/se)^2, 1)             # Wald-type test of beta

## Summaries of Beta and HR
cbind(beta = beta, se = se, z = beta/se, p_value = p_value) # summaries of beta
cbind(HR = exp(beta), lower.95 = exp(confin)[,1], upper.95 =  exp(confin)[,2])

## To test the joint significance of categorized age with Wald-type test
## Joint H0:
##          beta_(agegrp45-59) = 0
##          beta_(agegrp60-74) = 0
##          beta_(agegrp75+) = 0
statistic <- t(coef(fpmf)[4:6]) %*% solve(vcov(fpmf)[4:6,4:6]) %*% coef(fpmf)[4:6]
p_value <- 1 - pchisq(statistic, 3)
cbind(chi2_statistics =statistic, P_value = p_value)

## To test the overall effect of age with LR test
fpmf2 <- stpm2(Surv(time,event) ~ sex + year8594 , data=melanoma0, df=4)
anova(fpmf, fpmf2)

## @knitr h_time_varying
## (h) Change to time-varying effect of agegrp2:4
## NB: including main effect of agerp
fpmh <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
              data=melanoma0, tvc=list(agegrp2 = 2, agegrp3 = 2, agegrp4 = 2),
              smooth.formula=~ nsx(log(time),df=4)  )
## NB: no main effect of agegrp
fpmh1 <- stpm2(Surv(time,event) ~ sex + year8594 ,
               data=melanoma0, tvc=list(agegrp2 = 2, agegrp3 = 2, agegrp4 = 2),
               smooth.formula=~ nsx(log(time),df=4)  )

beta <- coef(fpmh)[2:6] ## log(HR)
se <- sqrt(diag(fpmh@vcov[2:6,2:6])) ## standard errors of beta
z <- beta/se ## z-scores
confin <- cbind(beta - 1.96*se, beta + 1.96*se) ## 95% confidence interval of beta=log(HR)
p_value <- 1-pchisq((beta/se)^2, 1) ## Wald-type test of beta

## Summaries of Beta and HR
cbind(beta = beta, se = se, z = beta/se, p_value = p_value) ## summaries of beta
cbind(HR = exp(beta), lower.95 = exp(confin)[,1], upper.95 =  exp(confin)[,2])

## LR test comparing fpmh (non-PH for agegrp2:4) with fpmf(PH for agegrp2:4)
anova(fpmh, fpmf)
anova(fpmh1, fpmf)

## I investigated the non-proportional effect of age with penalized models
## here, sp is the optimal smoothing parameters estimated from models without sp argument
pfit0 <- pstpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                smooth.formula=~s(log(time)), data=melanoma0, sp=0.1359685)

## The time-dependent effects including linear forms of age groups
pfit1 <- pstpm2(Surv(time,event) ~ sex + year8594,
                smooth.formula=~s(log(time)) + s(log(time),by=agegrp2) +
                                s(log(time),by=agegrp3) + s(log(time),by=agegrp4),
                data=melanoma0, sp=c( 0.1429949, 1.6133966, 1.3183117, 1.9958815))
anova(pfit1, pfit0)## the results also suggest there is strong evidence

## @knitr i_plot_base_haz
## (i) Plot of baseline hazard with fpmh
sexs <- levels(melanoma$sex)
years <- levels(melanoma$year8594)
agegrps <- levels(melanoma$agegrp)
par(mfrow=c(1,1))
newdata1 <- data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0)
plot(fpmh,newdata=newdata1, xlab="Time since diagnosis (years)", ylab="Hazard",
     type="haz", ci=FALSE, ylim=c(0,0.051), rug=FALSE)
legend("topright", legend=c(sexs[1], years[1], agegrps[1]), lty=1, bty="n")

## @knitr j_age_HR
plot(fpmh, newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="hr", rug=FALSE, line.col=1, ci=FALSE,log="y", ylim=c(1,500), lty=1,
     exposed=function(data) transform(data,sex=sexs[1],year8594=years[1],agegrp2=1, agegrp3=0, agegrp4=0),
     xlab="Time since diagnosis (years)",ylab="Hazards ratio")
plot(fpmh, newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="hr", rug=FALSE, line.col=1, ci=FALSE, add=T,log="y", lty=2,
     exposed=function(data) transform(data,sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=1, agegrp4=0),
     xlab="Time since diagnosis (years)",ylab="Hazards ratio")
plot(fpmh, newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="hr", rug=FALSE, line.col=1, ci=FALSE,add=T,log="y", lty=3,
     exposed=function(data) transform(data,sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=1),
     xlab="Time since diagnosis (years)",ylab="Hazards ratio")
legend("topright", legend=agegrps[2:4], lty=1:3)

## @knitr j_oldest
plot(fpmh, newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="hr", rug=FALSE, line.col=1, ci=T, log="y",
     exposed=function(data) transform(data,sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=1),
     xlab="Time since diagnosis (years)",ylab="Hazards ratio")

## @knitr k_haz_diff
plot(fpmh,newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="hdiff", rug=FALSE, line.col=1, ci=T,
     exposed=function(data) transform(data,sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=1),
     xlab="Time since diagnosis (years)",ylab="Hazards difference")

## @knitr l_surv_diff
plot(fpmh,newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
     type="sdiff", rug=FALSE, line.col=1, ci=T,
     exposed=function(data) transform(data, sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=1),
     xlab="Time since diagnosis (years)",ylab="Survivals difference")

## @knitr m_time_dep_eff
aicc <- bicc <- beta <- se <- rep(1,3)
nevent <- sum(melanoma0$event)
BIC <- function(object, nknots){
  -2 * as.numeric(bbmle::logLik(object)) + log(nevent) * (1 + nknots)
}
for (i in 1:3 ) {
  fitdf <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                  data=melanoma0, tvc=list(agegrp2 = i, agegrp3 = i, agegrp4 = i),
                  smooth.formula=~ nsx(log(time),df=4)  )
  aicc[i] <- AIC(fitaic)
  bicc[i] <- BIC(fitaic, i)
}
cbind(aicc=aicc, bicc=bicc)

## PLots with different df
fitdf1 <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4, data=melanoma0,
                tvc = list(agegrp2 = 1, agegrp3 = 1, agegrp4 = 1), smooth.formula=~ nsx(log(time),df=4))

plot(fitdf1, newdata = data.frame(sex=sexs[1], year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=0),
     type="hr", rug=FALSE, line.col=1, ci=F, log="y", lty=1,
     exposed=function(data) transform(data, sex = sexs[1],year8594=years[1],
                                        agegrp2=0, agegrp3 = 0, agegrp4=1),
     xlab="Time since diagnosis (years)", ylab="Hazards ratio")

for (i in 2:3 ) {
    fitdf <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                   data=melanoma0, tvc=list(agegrp2 = i, agegrp3 = i, agegrp4 = i),
                   smooth.formula=~ nsx(log(time),df=4))

    plot(fitdf, newdata=data.frame(sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0),
         type="hr", rug=FALSE, line.col=1, ci=F, log="y", add=T, lty=i,
         exposed=function(data) transform(data, sex=sexs[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=1),
         xlab="Time since diagnosis (years)",ylab="Hazards ratio")
}
legend("topright", legend=c("1 df", "2 df", "3 df"), lty=1:3)
