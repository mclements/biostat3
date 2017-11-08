muhaz2 <-
  function (formula, data, subset, max.time, ...) 
  {
    Call <- match.call()
    Call[[1]] <- as.name("muhaz2")
    indx <- match(c("formula", "data", "subset"), 
                  names(Call), nomatch = 0)
    if (indx[1] == 0) 
      stop("a formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(temp)
    Terms <- terms(formula, c("strata", "cluster"))
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
      stop("Interaction terms are not valid for this function")
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Surv(Y)) 
      stop("Response must be a survival object")
    if (!is.null(attr(Terms, "offset"))) 
      warning("Offset term ignored")
    ll <- attr(Terms, "term.labels")
    if (length(ll) == 0) 
      X <- factor(rep(1, n))
    else X <- strata(m[ll])
    if (!is.Surv(Y)) 
      stop("y must be a Surv object")
    ## newY <- aeqSurv(Y)
    if (missing(max.time)) {
      max.time <- NULL
      formula2 <- formula
      formula2[[3]] <- quote(1)
    }
    if (attr(Y, "type") == "right" || attr(Y, "type") == 
        "counting")  {
      temp <- tapply(1:nrow(Y), X, 
                     function(index) {
                       if (is.null(max.time)) {
                         sfit <- survfit(formula2, data[index,,drop=FALSE])
                         max.time <- approx(sfit$n.risk, sfit$time, xout = 10)$y
                         if (is.na(max.time) || max.time>max(Y[index,1]))
                           max.time <- max(Y[index,1])
                       }
                       val <- muhaz::muhaz(Y[index,1,drop=FALSE], 
                                                  Y[index,2,drop=FALSE], 
                                                  max.time = max.time, ...)
                       class(val) <- c("muhaz2","muhaz")
                       val
                       }, simplify=FALSE)
    } else {
      stop("survival type not supported")
    }
    if (length(temp)==1) {
      temp <- temp[[1]]
    }
    else {
        class(temp) <- "muhazList"
        names(temp) <- levels(strata(m[ll]))
        attr(temp,"vars") <- lapply(m[ll], levels)
    }
    structure(temp, call=Call)
  }
plot.muhaz2 <- function(x, haz.scale = 1, ylab="Hazard", ylim=NULL, log="", ...) {
    x$haz.est <- x$haz.est * haz.scale
    if (log %in% c("y","xy","yx")) x$haz.est <- ifelse(x$haz.est==0,NA,x$haz.est)
    if (is.null(ylim)) {
        ylim <- if (log %in% c("y","xy","yx"))
                    c(min(x$haz.est,na.rm=TRUE), max(x$haz.est,na.rm=TRUE)) else c(0, max(x$haz.est,na.rm=TRUE))
        }
    muhaz::plot.muhaz(x, ylab=ylab, ylim=ylim, log=log, ...)
}
lines.muhaz2 <- function(x, ..., haz.scale = 1) {
    x$haz.est <- x$haz.est * haz.scale
    muhaz::lines.muhaz(x, ...)
}
plot.muhazList <- function(x, lty=1:5, col=1:length(x), log="", legend.args=list(), ...) {
  lty <- rep(lty, length=length(x))
  col <- rep(col, length=length(x))
  est.grid <- unlist(lapply(x, "[[", "est.grid"))
  haz.est <- unlist(lapply(x, "[[", "haz.est"))
  if (log %in% c("y","xy","yx")) haz.est <- ifelse(haz.est==0,NA,haz.est)
  plot <- plot.muhaz2(list(est.grid=est.grid,haz.est=haz.est), type="n", log=log, ...)
  for(i in 1:length(x)) {
    lines.muhaz2(x[[i]], lty=lty[i], col=col[i], type="l", ...)
  }
  if (!is.logical(legend.args) || legend.args) {
      if (is.logical(legend.args)) legend.args <- list()
      base.legend.args <- list(x="topright",legend=names(x),col=col,lty=lty)
      legend.args <- do.call("updateList",c(list(base.legend.args), legend.args))
      do.call("legend", legend.args)
  }
  invisible(plot)
}
lines.muhazList <- function(x, lty=1, col=1:length(x), ...) {
  lty <- rep(lty, length=length(x))
  col <- rep(col, length=length(x))
  for(i in 1:length(x)) {
    lines(x[[i]], lty=lty[i], col=col[i], type="l", ...)
  }
}
summary.muhazList <- function(object, ...)
  lapply(object, muhaz::summary.muhaz)
as.data.frame.muhaz <- function(x, row.names, optional, ...) {
    if ("est.grid" %in% names(x)) {
        data.frame(x=x$est.grid, y=x$haz.est)
        } else {
            est.grid <- unlist(lapply(x, "[[", "est.grid"))
            haz.est <- unlist(lapply(x, "[[", "haz.est"))
            data.frame(x=est.grid,y=haz.est)
        }
}
as.data.frame.muhazList <- function(x, row.names, optional, ...) {
  est.grid <- unlist(lapply(x, "[[", "est.grid"))
  haz.est <- unlist(lapply(x, "[[", "haz.est"))
  strata <- unlist(lapply(1:length(x), function(i) rep(names(x)[i],length(x[[i]]$est.grid))))
  out <- data.frame(x=est.grid, y=haz.est, strata, row.names=1:length(strata))
  values <- do.call("expand.grid",rev(attr(x,"vars")))
  newdata <- values[rep(1:length(x), sapply(x, function(xi) length(xi$est.grid))),,drop=FALSE]
  rownames(newdata) <- 1:nrow(newdata)
  out <- cbind(out,newdata)
  out
}

plot.bshazard <-
function (x, conf.int = T, overall = T, col = 1, lwd = 1, xlab = "Time", 
    ylab = "Hazard rate", ...) 
{
    if (overall == T) {
        plot(x$time, x$hazard, xlab = xlab, type = "l", ylab = ylab, 
            lwd = lwd, lty = 1, col = col, ...)
        if (conf.int == T) {
            lines(x$time, x$low, lty = 2, col = col, lwd = lwd)
            lines(x$time, x$up, lty = 2, col = col, lwd = lwd)
        }
    }
    if (overall == F & !is.null(x$cov.value)) {
        covs <- unique(x$raw.data[, attr(x$cov.value, "names")])
        lin <- (covs - matrix(as.numeric(x$cov.value), nrow(covs), 
            ncol(covs), byrow = T)) * matrix(x$coefficients, 
            nrow(covs), ncol(covs), byrow = T)
        HR <- exp(rowSums(lin))
        plot(x$time, x$hazard, xlab = xlab, type = "n", ylab = ylab, 
            lwd = lwd, lty = 1, col = col, ...)
        for (i in 1:nrow(covs)) {
            h <- x$hazard * HR[i]
            lines(x$time, h, xlab = xlab, type = "l", ylab = ylab, 
                lwd = lwd, lty = 1, col = col)
        }
    }
}
