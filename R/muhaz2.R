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
    else class(temp) <- "muhazList"
    structure(temp, call=Call)
  }
plot.muhaz2 <- function(x, ..., haz.scale = 1) {
  x$haz.est <- x$haz.est * haz.scale
  muhaz::plot.muhaz(x, ...)
}
lines.muhaz2 <- function(x, ..., haz.scale = 1) {
  x$haz.est <- x$haz.est * haz.scale
  muhaz::lines.muhaz(x, ...)
}
plot.muhazList <- function(x, lty=1, col=1:length(x), ...) {
  lty <- rep(lty, length=length(x))
  col <- rep(col, length=length(x))
  est.grid <- unlist(lapply(x, "[[", "est.grid"))
  haz.est <- unlist(lapply(x, "[[", "haz.est"))
  plot.muhaz2(list(est.grid=est.grid,haz.est=haz.est), type="n", ...)
  for(i in 1:length(x)) {
    lines(x[[i]], lty=lty[i], col=col[i], type="l", ...)
  }
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
as.data.frame.muhaz <- function(x) {
  est.grid <- unlist(lapply(x, "[[", "est.grid"))
  haz.est <- unlist(lapply(x, "[[", "haz.est"))
  data.frame(est.grid,haz.est)
}
as.data.frame.muhazList <- function(x) {
  est.grid <- unlist(lapply(x, "[[", "est.grid"))
  haz.est <- unlist(lapply(x, "[[", "haz.est"))
  strata <- unlist(lapply(1:length(x), function(i) rep(names(x)[i],length(x[[i]]$est.grid))))
  data.frame(est.grid, haz.est, strata, row.names=1:length(strata))
}
