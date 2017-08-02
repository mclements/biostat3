require(epiR)
require(dplyr)
## wrapper function for epi.conf
ir.est <- function(df, grp, event, p.time, method="exact"){
    grpSums <- df %>%
        select_(.dots = c(grp, event, p.time)) %>%
            group_by_(grp) %>%
                summarise_each(funs(sum))
    cbind(grpSums,epi.conf(data.matrix(select_(grpSums,.dots=c(event, p.time))), ctype="inc.rate", method)) %>%
        rename(IR=est, 'CI lower'=lower, 'CI upper'=upper)
}

## wrapper function for muhaz hazard smoother by strata (todo: rewrite with dplyr)
smoothHazard <- function(df, strat){
    tmp <- lapply(as.list(sort(levels(df[,strat]))),
                  function(x) c(strat=x, with(df[df[strat]==x,],
                                    muhaz(times=surv_mm, delta=death_cancer, max.time = max(surv_mm)))))
    out <- do.call("rbind", lapply(tmp,function(obj)
                                   data.frame(obj$haz.est, obj$est.grid, obj$strat)))
    colnames(out) <- c("Hazard","Time", strat)
    return(out)
}

