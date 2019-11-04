## Purpose: Collapsibility
## Author: Mark Clements, 2019-10-23
###############################################################################

###############################################################################
## Exercise 14
###############################################################################
## @knitr loadDependencies
library(biostat3) 
library(rstpm2)
library(ggplot2)

## @knitr loadPreprocess

## @knitr 14.simulate
set.seed(12345)
d <- local({
    n <- 1e4
    x <- rbinom(n, 1, 0.5)
    u <- rnorm(n, 0, 3)
    t <- rexp(n, exp(-5+x+u))
    c <- runif(n, 0, 10)
    y <- pmin(t, c)
    delta <- (t < c)
    data.frame(y,x,u,delta)
})
head(d)


## @knitr 14.dag
library(dagitty)
g1 <- dagitty( "dag {
    X -> T -> \"(Y,Delta)\"
    U -> T
    C -> \"(Y,Delta)\"
}")
plot(graphLayout(g1))

## @knitr 14.a
summary(fit1 <- glm(delta~x+u+offset(log(y)),data=d,family=poisson))

summary(fit2 <- coxph(Surv(y,delta)~x+u,data=d))

summary(fit3 <- stpm2(Surv(y,delta)~x+u,data=d,df=4))

## summary table for the coefficients for X
rbind(Poisson=coef(summary(fit1))["x",c("Estimate","Std. Error")],
      Cox=coef(summary(fit2))["x",c("coef","se(coef)")],
      Stpm2=coef(summary(fit3))["x",c("Estimate","Std. Error")])

## @knitr 14.tvc.xu
fit <- stpm2(Surv(y,delta)~x+u,data=d,df=4, tvc=list(x=2))
plot(fit, type="hr", newdata=data.frame(x=0,u=0), var="x", ylim=c(1,4))
s <- predict(fit, type="surv", newdata=data.frame(x=0:1,u=3), grid=TRUE, full=TRUE,
             se.fit=TRUE)

ggplot(s, aes(x=y,y=Estimate,fill=factor(x),ymin=lower,ymax=upper)) +
    ylab("Survival") +
    geom_ribbon(alpha=0.6) +
    geom_line()

plot(fit, type="meanhr", newdata=transform(d,x=0), var="x", seqLength=31, ylim=c(1,4))
plot(fit, type="meansurvdiff", newdata=transform(d,x=0), var="x", seqLength=31)
plot(fit, type="meansurv", newdata=transform(d,x=0), seqLength=101)


## @knitr 14.b
summary(fit1 <- glm(delta~x+offset(log(y)),data=d,family=poisson))

summary(fit2 <- coxph(Surv(y,delta)~x,data=d))

summary(fit3 <- stpm2(Surv(y,delta)~x,data=d,df=4))

## @knitr 14.b.2
## summary table for the coefficients for X
rbind(Poisson=coef(summary(fit1))["x",c("Estimate","Std. Error")],
      Cox=coef(summary(fit2))["x",c("coef","se(coef)")],
      Stpm2=coef(summary(fit3))["x",c("Estimate","Std. Error")])

## @knitr 14.tvc.x
fit <- stpm2(Surv(y,delta)~x,data=d,df=4, tvc=list(x=2))
plot(fit, type="hr", newdata=data.frame(x=0), var="x", ylim=c(1,4))

## @knitr 14.simulate.2
set.seed(12345)
d <- local({
    n <- 1e4*10 # CHANGED N
    x <- rbinom(n, 1, 0.5)
    u <- rnorm(n, 0, 3)
    t <- rexp(n, exp(-5+x+u))
    c <- runif(n, 0, 10/1000) # CHANGED FROM 10 TO 0.01
    y <- pmin(t, c)
    delta <- (t < c)
    data.frame(y,x,u,delta)
})

## @knitr 14.simulate.3
set.seed(12345)
d <- local( {
    n <- 1e4
    x <- rbinom(n, 1, 0.5)
    u <- rnorm(n, 0, 1) # CHANGED SD FROM 3 TO 1
    t <- rexp(n, exp(-5+x+u))
    c <- runif(n, 0, 10) 
    y <- pmin(t, c)
    delta <- (t < c)
    data.frame(y,x,u,delta)
})

## @knitr 14.aft
fit <- aft(Surv(y,delta)~x+u,data=d,df=4)
summary(fit)
fit <- aft(Surv(y,delta)~x,data=d,df=4)
summary(fit)
