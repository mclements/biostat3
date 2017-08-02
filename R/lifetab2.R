lifetab2 <-
    function(surv, data=NULL, breaks=NULL) {
        stopifnot(require(KMsurv))
        y <- eval(substitute(surv),data,parent.frame())
        time <- y[,1]
        event <- y[,2]
        if (is.null(breaks))
            breaks <- c(sort(unique(time)), Inf)
        cut_time <- cut(time,breaks,include.lowest=TRUE,right=FALSE)
        NA2zero <- function(x) {if (any(is.na(x))) x[is.na(x)] <- 0; x}
        nevent <- NA2zero(tapply(event,cut_time,sum))
        nlost <- NA2zero(tapply(event,cut_time,length)) - nevent
        KMsurv::lifetab(tis = breaks, # should be one element longer for the intervals
                        ninit = nrow(data),           # number of individuals at the start
                        nlost = nlost,                 # number lost for each interval
                        nevent = nevent)           # number of events for each interval
    }
