## Using ideas from https://www.biostat.washington.edu/sites/default/files/modules//SISCR_2016_12_all.v2.pdf and http://intlpress.com/site/pub/files/_fulltext/journals/sii/2009/0002/0001/SII-2009-0002-0001-a004.pdf
coxphHaz <-
    function(object, newdata, n.grid = 300, kernel = "epanechnikov", from=NULL, to=NULL, ...){
        stopifnot(inherits(object,"coxph"))
        stopifnot(!missing(newdata))
        if(nrow(newdata)>1)
            return(structure(lapply(1:nrow(newdata),
                                    function(i) coxphHaz(object, newdata[i,,drop=FALSE], n.grid, kernel, from, to, ...)),
                             newdata=newdata,
                             class="coxphHazList"))
        x <- survival::survfit(object, conf.type = "none")
        index <- x$n.risk > 0 & x$n.event > 0
        time <- x$time[index]
        weights <- diff(c(0,x$cumhaz[index]))
        Sh0 <- x$n.event[index]/weights
        weights.var <- weights/Sh0
        if (is.null(from)) from <- min(x$time)
        if (is.null(to)) to <- max(x$time)
        newobject <- suppressWarnings(stats::density(time, weight = weights, kernel = kernel,
                                                     n=n.grid, from = from, to = to, ...))
        newobject$y = newobject$y*exp(predict(object, newdata))
        ## variance only correct for X=0 and for kernel="epanechnikov"
        newobject.var <- suppressWarnings(stats::density(time, weight = weights.var, kernel = "biweight",
                                                         n=n.grid, from = from, to = to, ...))
        newobject$var = 3/5/newobject$bw*newobject.var$y
        structure(newobject, newdata=newdata, n.grid=n.grid, call=match.call(),
                  class=c("coxphHaz","density"))
    }
print.coxphHaz <- function(x, digits=NULL, ...) {
    cat("\nnewdata:\n")
    print(attr(x,"newdata"))
    NextMethod()
}
plot.coxphHaz <- function(x, xlab="Time", ylab="Hazard", type="l", ...) {
    plot.default(x, xlab=xlab, ylab=ylab, type=type, ...)
}
plot.coxphHazList <- function(x, xlab="Time", ylab="Hazard", type="l", col=1:length(x), lty=1, legend.args=list(), ...) {
    plot <- matplot(x[[1]]$x, do.call("cbind", lapply(x, function(item) item$y)),
                    xlab=xlab, ylab=ylab, type=type, col=col, lty=lty, ...)
    base.legend.args <- list(x="topright",legend=strata(attr(x,"newdata")),col=col,lty=lty)
    legend.args <- modifyList(base.legend.args, legend.args)
    do.call("legend", legend.args)
    invisible(plot)
}
lines.coxphHazList <- function(x, ...)
    matlines(x[[1]]$x, do.call("cbind", lapply(x, function(item) item$y)), 
             ...)

as.data.frame.coxphHaz <- function(x, row.names=NULL, optional = FALSE, level=0.95, ...) {
    newdata = attr(x,"newdata")
    ## To avoid "row names were found from a short variable and have been discarded":
    rownames(newdata) = NULL
    alpha <- (1-level)/2
    transform(data.frame(newdata, x=x$x, y=x$y, var=x$var),
              y.lower = y*exp(qnorm(alpha)*sqrt(var)/y),
              y.upper = y*exp(qnorm(1-alpha)*sqrt(var)/y))
}
as.data.frame.coxphHazList <- function(x, row.names=NULL, optional = FALSE, ...) {
    do.call(rbind, lapply(x, as.data.frame, rownames, optional, ...))
}


## ## adjusted muhaz (not used)
## library(dplyr)
## haz <- as.data.frame(muhaz2(Surv(surv_mm,dead)~1, data=colon2))
## hazfun <- with(haz, approxfun(x = est.grid, y = haz.est))
## expbeta <- exp(predict(fit))
## colon3 <- dplyr::arrange(colon2,surv_mm)
## n <- nrow(colon3)
## nstar <- n-100
## predhaz0 <- sapply(1:nstar, function(i) hazfun(colon3[i,"surv_mm"])*(n-i+1)/sum(expbeta[i:n]))
## index <- tapply(1:nstar, colon3$surv_mm[1:nstar], function(i) floor(median(i)))
## predhaz0 <- predhaz0[index]
## times <- colon3$surv_mm[index]
## ## issue: many ties
## plot(times, predhaz0,type="l",col="blue",ylim=c(0,0.1))
## lines(times, predhaz0*exp(coef(fit)),col="blue",lty=2)
