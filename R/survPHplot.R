survPHplot <- function(formula, data, subset, contrasts, weights,
                       col=1:5, lty=1:5, pch=19,
                       xlab="Time (log scale)",
                       ylab="-log(-log(Survival))",
                       log="x",
                       legend.args=list(),
                       ...) {
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"), 
        names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1]] <- as.name("survfit")
    fit <- eval(mf, parent.frame())
    n <- length(fit$strata)
    index <- rep(names(fit$strata), fit$strata)
    time <- fit$time
    trans <- -log(-log(fit$surv))
    plot <- plot(trans~time, type="n", log=log, xlab=xlab, ylab=ylab, ...)
    col <- rep(col, length.out=n)
    lty <- rep(lty, length.out=n)
    for (i in 1:n) {
        j <- names(fit$strata)[i]==index
        lines(time[j],trans[j],col=col[i],lty=lty[i])
        points(time[j],trans[j],pch=19,col=col[i])
    }
    base.legend.args <- list(x="topright",legend=names(fit$strata),col=col,lty=lty,pch=pch)
    legend.args <- do.call("updateList",c(list(base.legend.args), legend.args))
    do.call("legend", legend.args)
    invisible(plot)
}
