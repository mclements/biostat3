
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



