poisson.ci <- function(x, T=1, conf.level = 0.95) {
  stopifnot(length(x)==length(T) || length(T)==1)
  T <- rep(T,length(x))
  val <- t(mapply(function(xi,Ti) 
    stats::poisson.test(xi, Ti, conf.level)$conf.int,
    x,
    T))
  val <- data.frame(val,row.names=names(x))
  a <- (1 - conf.level)/2
  a <- c(a, 1 - a)
  names(val) <- format_perc(a, 3)
  if (nrow(val)>1) val else as.vector(unlist(val))
}
