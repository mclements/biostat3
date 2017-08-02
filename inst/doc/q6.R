## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01, 2016-03-08
###############################################################################

###############################################################################
## Exercise 6
###############################################################################
## @knitr loadDependecies
require(foreign)  # for reading data set from Stata
require(survival) # for Surv and survfit
require(epiR)     # for epi.conf
require(dplyr)    # for data manipulation

## @knitr loadPreprocess
diet <- data.frame(read.dta("http://biostat3.net/download/diet.dta"))
head(diet)
summary(diet)


## @knitr 6a_ir
diet$y1k <- diet$y/1000

diet.ir6a <- diet %>%
    group_by(hieng) %>%
    summarise(Event = sum(chd), Time = sum(y1k)) %>%      # group sums
    cbind(select(., Event, Time) %>%                      # formating for epi.conf
          as.matrix() %>%                                 # formating for epi.conf
          epi.conf(ctype="inc.rate", method="exact")) %>% # calc rate + CI
    print()

## IRR
diet.ir6a[diet.ir6a$hieng=="high","est"] / diet.ir6a[diet.ir6a$hieng=="low","est"]

## @knitr 6b_ir
poisson6b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson6b)
exp(cbind(coef(poisson6b),confint(poisson6b)))


## @knitr 6c_energyDist
hist6c <- hist(diet$energy, breaks=25, probability=TRUE)
curve(dnorm(x, mean=mean(diet$energy), sd=sd(diet$energy)), col = "red", add=TRUE)
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
# For kurtosis and skewness, see package e1071

## @knitr 6d_engCat
diet$eng3 <- cut(diet$energy, breaks=c(1500,2500,3000,4500),labels=c("low","medium","high"))
diet %>% group_by(eng3) %>%
    summarise(Freq = n(), Percent = n()/dim(diet)[1]) %>%
    mutate(Cum = cumsum(Percent))

## @knitr 6e_irEng
## IRR
diet.ir6e <- diet %>%
    group_by(eng3) %>%
    summarise(Event = sum(chd), Time = sum(y1k)) %>%      # group sums
    cbind(select(., Event, Time) %>%                      # formating for epi.conf
          as.matrix() %>%                                 # formating for epi.conf
          epi.conf(ctype="inc.rate", method="exact")) %>% # calc rate + CI
    print()

diet.ir6e[diet.ir6e$eng3=="medium","est"]/diet.ir6e[diet.ir6e$eng3=="low","est"]
diet.ir6e[diet.ir6e$eng3=="high","est"]/diet.ir6e[diet.ir6e$eng3=="low","est"]

## @knitr 6f_irEng
dummies6f <- model.matrix(~eng3,data=diet)
diet$X2 <- dummies6f[,"eng3medium"]
diet$X3 <- dummies6f[,"eng3high"]
diet$X1 <- (1-diet$X2)*(1-diet$X3)
colSums(diet[c("X1","X2","X3")])

## @knitr 6g_irEng
head(diet[which(diet$eng3=="low"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="medium"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="high"),c("energy","eng3","X1","X2","X3")])

## @knitr 6h
poisson6h <- glm( chd ~ X2 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6h)
exp(cbind(coef(poisson6h),confint(poisson6h)))

## @knitr 6i
poisson6i <- glm( chd ~ X1 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6i )
exp(cbind(coef(poisson6i),confint(poisson6i)))

## @knitr 6j
poisson6j <- glm( chd ~ eng3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6j )
exp(cbind(coef(poisson6j),confint(poisson6j)))

## @knitr 6k
sums6k <- colSums(diet[c("chd","y")])
sums6k[1]/sums6k[2]
