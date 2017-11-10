smoothHaz <-
    function(object, n.grid = 300, kernel = "epanechnikov", from=NULL, to=NULL, min.n.risk=1, ...){
        x <- object
        index <- x$n.risk >= min.n.risk
        time <- x$time[index]
        weights <- x$n.event[index]/x$n.risk[index]
        if (is.null(from)) from <- min(x$time)
        if (is.null(to)) to <- max(x$time)
        newobject <- suppressWarnings(stats::density(time, weight = weights, kernel = kernel,
                                                     n=n.grid, from = from, to = to, ...))
        structure(newobject, n.grid=n.grid, call=match.call(),
                  class=c("smoothHaz","density"))
    }
plot.smoothHaz <- function(x, xlab="Time", ylab="Hazard", type="l", ...) {
    plot.default(x, xlab=xlab, ylab=ylab, type=type, ...)
}
