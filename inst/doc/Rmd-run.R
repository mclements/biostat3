library(rmarkdown)
library(prettydoc)
setwd("~/src/R/biostat3/inst/doc")

files <- paste0("q", c(1:4,6:13,22,23,25,28),".Rmd")
## files <- paste0("q", c(8),".Rmd")

rmarkdown::render("index.Rmd", html_pretty())

setwd("labs")
for (file in files) {
    rmarkdown::render(file, html_pretty())
}
setwd("..")

setwd("solutions")
for (file in files) {
    rmarkdown::render(file, html_pretty())
}
setwd("..")


## setwd("~/src/R/biostat3/inst/doc/solutions")
## rmarkdown::render("q2.Rmd", "html_document")
