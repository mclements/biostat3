## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01, 2016-03-08
## Edited: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 6
###############################################################################
## @knitr loadDependencies
library(biostat3)
library(dplyr)    # for data manipulation

## @knitr loadPreprocess
head(diet)
summary(diet)

## @knitr 6a_ir
diet <- biostat3::diet
diet$y1k <- diet$y/1000

diet.ir6a <- survRate(Surv(y/1000,chd) ~ hieng, data=diet)
## or
diet %>%
    group_by(hieng) %>%
    summarise(Event = sum(chd), Time = sum(y1k), Rate = Event/Time,      # group sums
              CI_low = poisson.test(Event,Time)$conf.int[1],
              CI_high = poisson.test(Event,Time)$conf.int[2]) 

## IRR
with(diet.ir6a, poisson.test(event,tstop)) 
with(diet.ir6a, poisson.test(rev(event),rev(tstop)))

## @knitr 6b_ir
poisson6b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson6b)
eform(poisson6b)
eform(poisson6b, method="Profile")

## @knitr 6c_energyDist
hist6c <- hist(diet$energy, breaks=25, probability=TRUE, xlab="Energy (units)")
curve(dnorm(x, mean=mean(diet$energy), sd=sd(diet$energy)), col = "red", add=TRUE)
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
# For kurtosis and skewness, see package e1071

## @knitr 6d_engCat
diet$eng3 <- cut(diet$energy, breaks=c(1500,2500,3000,4500),labels=c("low","medium","high"), 
                 right = FALSE)
cbind(Freq=table(diet$eng3),
      Prop=table(diet$eng3)/nrow(diet))


## @knitr 6e_irEng
## rates and IRRs
diet.ir6e <- survRate(Surv(y/1000,chd) ~ eng3, data=diet)
print(diet.ir6e)

# calculate IRR and confidence intervals
with(diet.ir6e, rate[eng3=="medium"] / rate[eng3=="low"])
with(diet.ir6e[c(2,1),], { # compare second row with first row
  poisson.test(event, tstop)
})
with(diet.ir6e, rate[eng3=="high"] / rate[eng3=="low"])
with(diet.ir6e[c(3,1),], { # compare third row with first row
  poisson.test(event, tstop)
})


## @knitr 6f_irEng
diet <- mutate(diet, 
               X1 = as.numeric(eng3 == "low"),
               X2 = as.numeric(eng3 == "medium"),
               X3 = as.numeric(eng3 == "high"))
# or
diet <- biostat3::addIndicators(diet, ~eng3+0) %>%
    mutate(X1 = eng3low, X2 = eng3medium, X3 = eng3high)
colSums(diet[c("X1","X2","X3")])

## @knitr 6g_irEng
filter(diet, eng3=="low")    %>% select(c(energy,eng3,X1,X2,X3)) %>% head
filter(diet, eng3=="medium") %>% select(c(energy,eng3,X1,X2,X3)) %>% head
filter(diet, eng3=="high")   %>% select(c(energy,eng3,X1,X2,X3)) %>% head

## @knitr 6h
poisson6h <- glm( chd ~ X2 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary(poisson6h)
eform(poisson6h)

## @knitr 6i
poisson6i <- glm( chd ~ X1 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
# or 
poisson6i <- glm( chd ~ I(eng3=="low") + I(eng3=="high") + offset( log( y1k ) ), family=poisson, data=diet )

summary( poisson6i )
eform( poisson6i )

## @knitr 6j
poisson6j <- glm( chd ~ eng3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6j )
eform( poisson6j )

## @knitr 6k
summarise(diet, rate = sum(chd) / sum(y))
