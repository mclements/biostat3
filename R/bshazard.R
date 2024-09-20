as.data.frame.bshazard <- function(x, ...) {
    with(x, data.frame(time,hazard,conf.low=lower.ci,conf.high=upper.ci))
}
