
library(foreign)
getData <- function(name) {
    assign(name,read.dta(sprintf("http://biostat3.net/download/%s.dta",name)))
    save(list=name,file=sprintf("~/src/R/biostat3/data/%s.rda",name))
    prompt(name=name,file=sprintf("~/src/R/biostat3/man/%s.Rd",name))
}
getData("melanoma")
getData("colon")
getData("brv")
getData("colon_sample")
getData("diet")
getData("popmort")

load("~/src/R/biostat3/data/colon.rda")
load("~/src/R/biostat3/data/melanoma.rda")
load("~/src/R/biostat3/data/diet.rda")
year <- function(date) as.numeric(format(as.Date(date),"%Y"))
dmy <- function(day,month,year)
  as.Date(paste(year,month,day,sep="-"), "%Y-%m-%d")
Date2year <- function(date) {
    y <- year(date)
    len <- as.numeric(dmy(31,12,y)-dmy(1,1,y))+1
    y + as.numeric(date-dmy(1,1,y))/len
}
Date2year(as.Date("2015-07-01"))
year2Date <- function(year) {
    y <- floor(year)
    len <- as.numeric(dmy(31,12,y)-dmy(1,1,y))+1
    dmy(1,1,y) + (year-y)*len
}
## year2Date(Date2year(as.Date("2011-12-31")))
## year2Date(Date2year(as.Date("2012-12-31"))) # leap-year

set.seed(12345)
colon2 <-
    transform(colon,
              id=1:nrow(colon),
              ydx=Date2year(dx),
              yexit=Date2year(exit),
              bdate=year2Date(Date2year(dx)-(age+runif(nrow(colon)))))
colon2 <- transform(colon2, ybdate = Date2year(bdate))
melanoma2 <-
    transform(melanoma,
              id=1:nrow(melanoma),
              ydx=Date2year(date = dx),
              yexit=Date2year(exit),
              bdate=year2Date(Date2year(dx)-(age+runif(nrow(melanoma)))))
melanoma2 <- transform(melanoma2, ybdate = Date2year(bdate))
diet <- transform(diet, yob=Date2year(dob), yoe=Date2year(doe), yox=Date2year(dox))

colon <- colon2
melanoma <- melanoma2
save(colon, file="~/src/R/biostat3/data/colon.rda", version=2)
save(melanoma, file="~/src/R/biostat3/data/melanoma.rda", version=2)
save(diet, file="~/src/R/biostat3/data/diet.rda", version=2)
