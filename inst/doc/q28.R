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
library(ggplot2)

## @knitr loadPreprocess
data(melanoma)

## Extract subset and create 0/1 outcome variables
melanoma0 <- biostat3::melanoma %>% filter(stage=="Localised") %>%
             mutate(event = ifelse(status=="Dead: cancer" & surv_mm<120, 1, 0),
                    time = pmin(120, surv_mm)/12,
                    agegrp1 = (agegrp=="0-44")+0,  # used by time-dependent effect
                    agegrp2 = (agegrp=="45-59")+0, # used by time-dependent effect
                    agegrp3 = (agegrp=="60-74")+0, # used by time-dependent effect
                    agegrp4 = (agegrp=="75+")+0)   # used by time-dependent effect

## @knitr a_flex
## (a) Flexible parametric model with df=4
fpma <- stpm2(Surv(time,event) ~ year8594, data=melanoma0, df=4)
summary(fpma)
eform(fpma)["year8594Diagnosed 85-94",]

## @knitr a_cox
cox <- coxph(Surv(time, event) ~ year8594,
             data=melanoma0) # year8594 is a categorical variable
summary(cox)

## @knitr b_surv
## (b) Prediction and plots of survival and hazard by calendar period
years <- levels(melanoma0$year8594)

s <- predict(fpma,newdata=data.frame(year8594=years),grid=TRUE,full=TRUE,se.fit=TRUE,
             type="surv")
head(s)
ggplot(s, aes(x=time,y=Estimate,fill=year8594,ymin=lower,ymax=upper)) +
    xlab("Time since diagnosis (years)") +
    ylab("Survival") +
    geom_ribbon(alpha=0.6) +
    geom_line()

## @knitr b_haz
plot(fpma,newdata=data.frame(year8594=years[1]), type="haz",
     xlab="Time since diagnosis (years)", ylab="Hazard")
lines(fpma,newdata=data.frame(year8594=years[2]), type="haz", col=2)
legend("topright", legend=years, col=1:2, lty=1, bty="n")


## @knitr c_haz_log
## (c) hazards on log scale, adding log="y"
plot(fpma,newdata=data.frame(year8594=years[1]), type="haz",
     xlab="Time since diagnosis (years)",
     ci=FALSE,
     ylab="Hazard (log scale)",main=years[1], log="y", ylim=c(0.01,0.07))
lines(fpma,newdata=data.frame(year8594=years[2]), type="haz", col=2)
legend("topright", legend=years, col=1:2, lty=1, bty="n")

## @knitr d_AIC_BIC
## utility function to row bind from a list
Rbind <- function(object) do.call(rbind,object)
out <- lapply(1:6, function(i) {
    fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
    data.frame(
        i,
        AIC=AIC(fitaic),
        BIC=BIC(fitaic),
        beta=as.numeric(coef(fitaic)[2]),
        se=coef(summary(fitaic))[2,2])
})
out %>% Rbind

## @knitr e_base_surv
## Baseline survival
fitaic0 <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=6)
plot(fitaic0,newdata=data.frame(year8594=years[1]), lty=6, ci=FALSE,
     xlab="Time since diagnosis (years)")
for (i in 1:5 ) {
  fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
  lines(fitaic,newdata=data.frame(year8594=years[1]), lty=i)
}
legend("topright", legend=paste0("df=",1:6), lty=1:6)

## @knitr e_base_haz
## Baseline hazard
fitaic1 <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=6)
plot(fitaic1,newdata=data.frame(year8594=years[1]), lty=6, type="haz",
     ci=FALSE, xlab="Time since diagnosis (years)", ylab="Hazard")
for (i in 1:5 ) {
  fitaic <- stpm2(Surv(time, event) ~ year8594, data=melanoma0, df=i)
  lines(fitaic,type="haz",newdata=data.frame(year8594=years[1]), lty=i)
}
legend("topright", legend=paste0("df=",1:6), lty=1:6)

## @knitr f_sex_age
fpmf <- stpm2(Surv(time, event) ~ sex + year8594 + agegrp,
              data=melanoma0, df=4)
summary(fpmf)
eform(fpmf)[2:6,]

## To test the overall effect of age with LR test
fpmf2 <- stpm2(Surv(time,event) ~ sex + year8594, data=melanoma0, df=4)
anova(fpmf, fpmf2)


## @knitr g_sex_age
summary(fit <- coxph(Surv(time, event) ~ sex + year8594 + agegrp,
                     data=melanoma0))
anova(fit)

## @knitr h_time_varying
## (h) Change to time-varying effect of agegrp2:4
## NB: including main effect of agegrp
fpmh <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
              data=melanoma0, tvc=list(agegrp2 = 2, agegrp3 = 2, agegrp4 = 2),
              df=4)
summary(fpmh)

## LR test comparing fpmh (non-PH for agegrp2:4) with fpmf(PH for agegrp2:4)
anova(fpmh, fpmf)

## @knitr h_time_varying_penalised
## I investigated the non-proportional effect of age with penalized models
## here, sp is the optimal smoothing parameters estimated from models without sp argument
pfit0 <- pstpm2(Surv(time,event) ~ sex + year8594 + agegrp,
                data=melanoma0, sp=0.1359685)

## The time-dependent effects including linear forms of age groups
pfit1 <- pstpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                tvc=list(agegrp2=7,agegrp3=7,agegrp4=7),
                data=melanoma0, sp=c( 0.1429949, 1.6133966, 1.3183117, 1.9958815))
anova(pfit1, pfit0)## the results also suggest there is strong evidence

## @knitr i_plot_base_haz
## Plot of baseline hazard with fpmh
sexes <- levels(melanoma$sex)
years <- levels(melanoma$year8594)
agegrps <- levels(melanoma$agegrp)
newdata1 <- data.frame(sex=sexes[1],year8594=years[1],agegrp2=0, agegrp3=0, agegrp4=0)
plot(fpmh, newdata=newdata1, xlab="Time since diagnosis (years)", 
     type="haz")

## @knitr j_age_HR
plot(fpmh, newdata=newdata1, xlab="Time since diagnosis (years)",
     type="hr",var="agegrp2", ci=FALSE, ylim=c(0,6))
lines(fpmh, newdata=newdata1, type="hr", var="agegrp3", lty=2)
lines(fpmh, newdata=newdata1, type="hr", var="agegrp4", lty=3)
legend("topright", legend=paste0(agegrps[-1]," vs 0-44"), lty=1:3)

## @knitr j_oldest
plot(fpmh, newdata=newdata1,
     type="hr", log="y",
     exposed=function(data) transform(data,agegrp4=agegrp4+1), # same as var="agegrp4"
     xlab="Time since diagnosis (years)")

## @knitr j_age_HR_ggplot
pred <- lapply(2:4, function(i)
    predict(fpmh, newdata=newdata1, type="hr",var=paste0("agegrp",i),
            grid=TRUE, se.fit=TRUE, full=TRUE) %>%
    mutate(ageGroup=paste0(agegrps[i]," vs 0-44")))
pred <- Rbind(pred)
ggplot(pred, aes(x=time,y=Estimate,fill=ageGroup,ymin=lower,ymax=upper)) +
    geom_ribbon(alpha=0.3) +
    geom_line() +
    ylim(0,8) +
    xlab("Time since diagnosis (years)") +
    ylab("Hazard ratio") +
    labs(fill="Age group")

## @knitr k_haz_diff
plot(fpmh,newdata=newdata1,
     type="hdiff", var="agegrp4",
     xlab="Time since diagnosis (years)")

## @knitr l_surv_diff
plot(fpmh,newdata=newdata1,
     type="sdiff", var="agegrp4",
     xlab="Time since diagnosis (years)")

## @knitr m_time_dep_eff
out <- lapply(1:3, function(i) {
    fitdf <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                   data=melanoma0, tvc=list(agegrp2 = i, agegrp3 = i, agegrp4 = i),
                   df=4)
    data.frame(i,
               aic=AIC(fitdf),
               bic=BIC(fitdf))
}) 
out %>% Rbind

## Plots with different df
fitdf1 <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4, data=melanoma0,
                tvc = list(agegrp2 = 1, agegrp3 = 1, agegrp4 = 1), df=4)

plot(fitdf1, newdata = newdata1,
     type="hr", ci=FALSE, log="y", var="agegrp4",
     xlab="Time since diagnosis (years)")
for (i in 2:3 ) {
    fitdf <- stpm2(Surv(time,event) ~ sex + year8594 + agegrp2 + agegrp3 + agegrp4,
                   data=melanoma0, tvc=list(agegrp2 = i, agegrp3 = i, agegrp4 = i),
                   df=4)
    lines(fitdf, newdata=newdata1,
          type="hr", lty=i,
          var="agegrp4")
}
legend("topright", legend=c("1 df", "2 df", "3 df"), lty=1:3)
