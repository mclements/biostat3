lincom <- function(model, specification, level = 0.95, eform=FALSE, ...) {
    if(!requireNamespace("car"))
        stop("Install the 'car' package to use the lincom function")
    if (length(specification)>1)
        return(t(sapply(specification, function(spec) lincom(model, spec, level, eform, ...))))
    x <- car::linearHypothesis(model, specification, ...)
    cf <- as.vector(attr(x, "value"))
    ses <- sqrt(as.vector(attr(x,"vcov")))
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format_perc(a, 3)
    fac <- qnorm(a)
    ci <- cbind(cf, cf + ses %o% fac)
    if (eform) ci <- exp(ci)
    dimnames(ci) <- list(specification, c("Estimate",pct))
    cbind(ci, x[-1,3:4])
}
