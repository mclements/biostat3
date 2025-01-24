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
mel <- subset(biostat3::melanoma, stage=="Localised") |>
    transform(dc = (status=="Dead: cancer" & surv_mm<120)+0,
              surv_10y = pmin(120, surv_mm))
with(mel,table(dc, status))

## @knitr ex_25_coxph
str(subset(mel, select=c(sex,year8594,agegrp))) # Check structure of risk factors/confounders
out_coh <- coxph(Surv(surv_10y,dc) ~ sex + year8594 + agegrp, data = mel)
summary(out_coh)


## @knitr n_ind
print(n_ind <- length(mel$id))

## @knitr n_event
table(mel$dc, useNA="always")
print(ncase <-  table (mel$dc, useNA="always")[2])

## @knitr gen_ncc
set.seed(12345)
nccdata <- ccwc(entry=0, exit=surv_10y , fail=dc, origin=0, controls=1,
                include=list(sex,year8594,agegrp,dc,id), data=mel)
tail(nccdata, 8) |> kable("html")


## @knitr clogit
out_ncc <- clogit(Fail ~ sex + year8594 + agegrp + strata(Set), data=nccdata)
summary(out_ncc)


## @knitr n_unique_ncc
print(n_uni <- length(unique(nccdata$id)))
n_uni

## @knitr compare_coh_ncc
library(broom)
cat("Full cohort:\n")
print(tidy_coh <- tidy(out_coh, conf.int=TRUE, exponentiate=TRUE))
cat("Nested case-control study:\n")
print(tidy_ncc <- tidy(out_ncc, conf.int=TRUE, exponentiate=TRUE))
tibble(term = tidy_coh$term, variance.ratio = (tidy_ncc$std.error/tidy_coh$std.error)^2)

tidy_coh <- tidy(out_coh)
tidy_ncc <- tidy(out_ncc)
f <- function(x,mu,sigma) {
    mu=as.numeric(mu)
    sigma=as.numeric(sigma)
    xp=exp(mu+x*sigma)
    cbind(x=xp,y=dlnorm(xp,mu,sigma))
}
x <- seq(-4,4,length=301)
par(mfrow=c(2,3))
for (i in 1:5) {
    f_coh = f(x,tidy_coh[i,"estimate"], tidy_coh[i,"std.error"])
    f_ncc = f(x,tidy_ncc[i,"estimate"], tidy_ncc[i,"std.error"])
    plot(rbind(f_coh,f_ncc), type="n", xlab="Hazard ratio", ylab="Density",
         main=as.character(tidy_coh$term[i]))
    polygon(f_coh,
            col=alpha("green",0.2), border=alpha("green",0.2))
    polygon(f_ncc,
            col=alpha("blue",0.2), border=alpha("blue",0.2))
}
legend("topright", c("Cohort","NCC"), col=c(alpha("green",0.2),alpha("blue",0.2)), lwd=5,
       bty="n")

## @knitr loop_ncc
set.seed(54321)
M <- 20
tidys_ncc <- lapply(1:M, function(i)
    ccwc(entry=0, exit=surv_10y , fail=dc, origin=0, controls=1,
         include=list(sex,year8594,agegrp), data=mel, silent=TRUE) |>
    clogit(formula = Fail ~ sex + year8594 + agegrp + strata(Set)) |> 
    suppressWarnings() |> tidy())
tidy_coh <- tidy(out_coh)
f <- function(x,mu,sigma) {
    mu=as.numeric(mu)
    sigma=as.numeric(sigma)
    xp=exp(mu+x*sigma)
    cbind(x=xp,y=dlnorm(xp,mu,sigma))
}
x <- seq(-4,4,length=301)
par(mfrow=c(2,3))
for (i in 1:5) {
    f_coh = f(x,tidy_coh[i,"estimate"], tidy_coh[i,"std.error"])
    f_ncc = lapply(tidys_ncc, function(object) 
        f(x,object[i,"estimate"], object[i,"std.error"]))
    plot(do.call(what=rbind,c(list(f_coh),f_ncc)), type="n", xlab="Hazard ratio", ylab="Density",
         main=as.character(tidy_coh$term[i]))
    for (j in 1:length(tidys_ncc))
        polygon(f_ncc[[j]],
                col=alpha("blue",0.02), border=alpha("blue",0.02))
    polygon(f_coh,
            col=alpha("green",0.2), border=alpha("green",0.2))
}
plot(0:1,0:1,type="n",axes=FALSE, xlab="",ylab="")
legend("center", c("Cohort","NCC"), col=c(alpha("green",0.2),alpha("blue",0.2)), lwd=5,
       bty="n")
