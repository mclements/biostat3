## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01, 2016-03-08
## Edited: Mark Clements, 2017-08-07
###############################################################################

###############################################################################
## Exercise 6
###############################################################################
## @knitr loadDependencies
library(biostat3) # diet dataset
library(dplyr)    # for data manipulation
library(knitr)    # kable()
library(broom)    # tidy()

## @knitr loadPreprocess
head(diet) |> kable("html")
summary(diet) |> kable("html")

## @knitr 6a_ir
diet <- biostat3::diet

diet.ir6a <- survRate(Surv(y/1000,chd) ~ hieng, data=diet)
kable(diet.ir6a, "html")
## or
diet |>
    group_by(hieng) |>
    summarise(Event = sum(chd), Time = sum(y/1000), Rate = Event/Time,      # group sums
              lower.ci = poisson.test(Event,Time)$conf.int[1],
              upper.ci = poisson.test(Event,Time)$conf.int[2]) |> kable("html")

## @knitr 6a_ir_b
## IRR
with(diet.ir6a, poisson.test(event,tstop)) 
with(diet.ir6a, poisson.test(rev(event),rev(tstop)))

## @knitr 6b_ir
poisson6b <- glm(chd ~ hieng + offset( log(y/1000)), family=poisson, data=diet)
summary(poisson6b)
eform(poisson6b)
eform(poisson6b, method="Profile")
## using tidyverse...
tidy(poisson6b, conf.int=TRUE, exponentiate=TRUE)

## @knitr 6c_energyDist
hist(diet$energy, breaks=25, probability=TRUE, xlab="Energy (units)")
curve(dnorm(x, mean=mean(diet$energy), sd=sd(diet$energy)), col = "red", add=TRUE)
lines(density(diet$energy), col = "blue")
legend("topright", legend=c("Normal density","Smoothed density"),
       lty=1, col=c("red", "blue"), bty="n")
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
# For kurtosis and skewness, see package e1071

## @knitr 6d_engCat
diet <- transform(diet, eng3=cut(energy, breaks=c(1500,2500,3000,4500),
                                 labels=c("low","medium","high"), 
                                 right = FALSE))
cbind(Freq=table(diet$eng3),
      Prop=proportions(table(diet$eng3))) |> kable("html")


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
diet <- transform(diet, 
                  X1 = as.numeric(eng3 == "low"),
                  X2 = as.numeric(eng3 == "medium"),
                  X3 = as.numeric(eng3 == "high"))
## or (names eng3low, eng3medium and eng3high)
diet = cbind(diet, model.matrix(~eng3+0, diet))
colSums(diet[c("X1","X2","X3","eng3low","eng3medium","eng3high")])

## @knitr 6g_irEng
subset(diet, eng3=="low", select=c(energy,eng3,X1,X2,X3)) |> head() |> kable("html")
subset(diet, eng3=="medium", select=c(energy,eng3,X1,X2,X3)) |> head() |> kable("html")
subset(diet, eng3=="high", select=c(energy,eng3,X1,X2,X3)) |> head() |> kable("html")

## @knitr 6h
poisson6h <- glm(chd ~ X2 + X3 + offset(log(y/1000)), family=poisson, data=diet)
summary(poisson6h)
tidy(poisson6h, conf.int=TRUE, exponentiate=TRUE)

## @knitr 6i
poisson6i <- glm(chd ~ X1 + X3 + offset(log(y/1000)), family=poisson, data=diet)
# or 
poisson6i <- glm(chd ~ I(eng3=="low") + I(eng3=="high") + offset(log(y/1000)),
                 family=poisson, data=diet)
summary(poisson6i)
tidy(poisson6i, conf.int=TRUE, exponentiate=TRUE)

## @knitr 6j
poisson6j <- glm(chd ~ eng3 + offset(log(y/1000)), family=poisson, data=diet)
summary(poisson6j)
tidy(poisson6j, conf.int=TRUE, exponentiate=TRUE)

## @knitr 6k
summarise(diet, rate = sum(chd) / sum(y))
## or
survRate(Surv(y/1000,chd)~1, data=diet)
