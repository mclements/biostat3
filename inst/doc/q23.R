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


## @knitr loadPreprocess

## @knitr 23.a
data(melanoma)
scale <- 365.24
mel <- mutate(melanoma,
              ydx=biostat3::year(dx),
              adx=age+0.5, # mid-point approximation
              dead=(status %in% c("Dead: cancer","Dead: other") & surv_mm<110)+0,
              surv_mm=pmin(110,surv_mm),
              astart=adx, 
              astop=adx+surv_mm/12)
mel.split <- survSplit(mel,
                       cut=1:110,
                       event="dead",start="astart", end="astop")
subset(mel.split, id<=2, select=c(id,astart,astop,dead))

## @knitr 23.b
mel.split <- mutate(mel.split,
                    ystart=year(dx)+astart-adx,
                    ystop=year(dx)+astop-adx)
mel.split2 <- survSplit(mel.split,
                       cut=1970:2000,event="dead",
                       start="ystart", end="ystop") %>%
    mutate(astart=adx+ystart-ydx,
           astop=adx+ystop-ydx,
           age=floor(astop),
           year=floor(ystop),
           pt = ystop - ystart)
subset(mel.split2, id<=2, select=c(id,ystart,ystop,astart,astop,dead))


## @knitr 23.c
xtabs(pt ~ age+year, data=mel.split2, subset = age>=50 & age<60)
xtabs(dead ~ age+year, data=mel.split2, subset = age>=50 & age<60)


## @knitr 23.d
mel.split2 <- mutate(mel.split2,
                     age10=cut(age,seq(0,110,by=10),right=FALSE),
                     year10=cut(year,seq(1970,2000,by=5),right=FALSE))
head(survRate(Surv(pt,dead)~sex+age10+year10, data=mel.split2))


## @knitr 23.e

pt <- mutate(mel.split2,sex=unclass(sex)) %>%
    group_by(sex, age, year) %>%
    summarise(pt=sum(pt))
expected <- inner_join(popmort, pt) %>%
    mutate(pt=ifelse(is.na(pt),0,pt)) %>%
    group_by(sex,year) %>%
    summarise(E=sum(rate*pt)) %>% ungroup
observed <- mutate(mel.split2, sex=as.numeric(unclass(sex))) %>%
    group_by(sex, year) %>%
    summarise(O=sum(dead)) %>% ungroup
joint <- inner_join(observed,expected) %>%
    mutate(SMR = O/E)

## @knitr 23.f

## overall SMRs
by(joint, joint$sex, function(data) poisson.test(sum(data$O), sum(data$E)))

## utility function to draw a confidence interval
polygon.ci <- function(time, interval, col="lightgrey") 
    polygon(c(time,rev(time)), c(interval[,1],rev(interval[,2])), col=col, border=col)

## modelling by calendar period
summary(fit <- glm(O ~ sex*ns(year,df=3)+offset(log(E)), data=joint, family=poisson))
##
pred <- predict(fit,type="response",newdata=mutate(joint,E=1),se.fit=TRUE)
full <- cbind(mutate(joint,fit=pred$fit), confint.predictnl(pred))
ci.cols <- c("lightgrey", "grey")
matplot(full$year, full[,c("2.5 %", "97.5 %")], type="n", ylab="SMR", xlab="Calendar year")
for (i in 1:2) {
    with(subset(full, sex==i), {
        polygon.ci(year, cbind(`2.5 %`, `97.5 %`), col=ci.cols[i])
    })
}
for (i in 1:2) {
    with(subset(full, sex==i), {
        lines(year,fit,col=i)
    })
}
legend("topright", legend=levels(mel.split2$sex), lty=1, col=1:2, bty="n")
