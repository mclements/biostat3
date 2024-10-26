as.data.frame.bshazard <- function(x, ...) {
    with(x, data.frame(time,hazard,conf.low=lower.ci,conf.high=upper.ci))
}

## plot.bshazard <-
## function (x, conf.int = T, overall = T, col = 1, lwd = 1, xlab = "Time", 
##     ylab = "Hazard rate", ...) 
## {
##     if (overall == T) {
##         plot(x$time, x$hazard, xlab = xlab, type = "l", ylab = ylab, 
##             lwd = lwd, lty = 1, col = col, ...)
##         if (conf.int == T) {
##             lines(x$time, x$low, lty = 2, col = col, lwd = lwd)
##             lines(x$time, x$up, lty = 2, col = col, lwd = lwd)
##         }
##     }
##     if (overall == F & !is.null(x$cov.value)) {
##         covs <- unique(x$raw.data[, attr(x$cov.value, "names")])
##         lin <- (covs - matrix(as.numeric(x$cov.value), nrow(covs), 
##             ncol(covs), byrow = T)) * matrix(x$coefficients, 
##             nrow(covs), ncol(covs), byrow = T)
##         HR <- exp(rowSums(lin))
##         plot(x$time, x$hazard, xlab = xlab, type = "n", ylab = ylab, 
##             lwd = lwd, lty = 1, col = col, ...)
##         for (i in 1:nrow(covs)) {
##             h <- x$hazard * HR[i]
##             lines(x$time, h, xlab = xlab, type = "l", ylab = ylab, 
##                 lwd = lwd, lty = 1, col = col)
##         }
##     }
## }
