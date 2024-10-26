poisson.ci <-
    function (x, T = 1, alternative = c("two.sided", "less", "greater"),
              conf.level = 0.95) {
        if ((l <- length(x)) != length(T)) 
            if (length(T) == 1L) 
                T <- rep(T, l)
            else stop("'x' and 'T' have incompatible length")
        xr <- round(x)
        if (any(!is.finite(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
            stop("'x' must be finite, nonnegative, and integer")
        x <- xr
        if (any(is.na(T) | (T < 0))) 
            stop("'T' must be nonnegative")
        if ((k <- length(x)) < 1L) 
            stop("not enough data")
        alternative <- match.arg(alternative)
        p.L <- function(x, alpha) {
            ifelse(x == 0, 0, qgamma(alpha, x))
        }
        p.U <- function(x, alpha) qgamma(1 - alpha, x + 1)
        CINT <- switch(alternative,
                       less = cbind(0, p.U(x, 1 - conf.level)), 
                       greater = cbind(p.L(x, 1 - conf.level), Inf),
                       two.sided = {
                           alpha <- (1 - conf.level)/2
                           cbind(p.L(x, alpha), p.U(x, alpha))
                       }) / cbind(T,T) |> as.data.frame()
        names(CINT) <- format_perc(switch(alternative,
                                          less = c(0, conf.level),
                                          greater = c(1-conf.level, 1),
                                          two.sided = c((1 - conf.level)/2,
                                                        1-(1 - conf.level)/2)),
                                   3)
        CINT
    }

