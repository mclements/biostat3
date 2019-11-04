year <- function(date, trunc=FALSE, year.length=365.24) {
    stopifnot(class(date)=="Date")
    year <- as.numeric(format(date,"%Y"))
    if (trunc) return(year)
    start <- as.Date(sprintf("%i-01-01",year))
    year + as.numeric(date-start)/year.length
}

## confint.anova <- function(x, level = 0.95) {
##       if(!all(c("value","vcov") %in% names(attributes(x))))
##         stop("confint.anova is only intended for linearHypothesis objects")
##       cf <- as.vector(attr(x, "value"))
##       vcov <- attr(x, "vcov")
##       ses <- sqrt(diag(vcov))
##       a <- (1 - level)/2
##       a <- c(a, 1 - a)
##       pct <- format.perc(a, 3)
##       fac <- qnorm(a)
##       ci <- cbind(cf, cf + ses %o% fac)
##       dimnames(ci) <- list(rownames(vcov), c("Estimate",pct))
##       ci
##   }

updateList <- function(object, ...) {
    revised <- list(...)
    for (name in setdiff(names(object),names(revised)))
        revised[[name]] <- object[[name]]
    revised
}

format_perc <- 
    function (probs, digits) 
        paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
              "%")
