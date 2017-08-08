lifetab2 <-
  function (formula, data, subset, breaks=NULL, ...) 
  {
    Call <- match.call()
    Call[[1]] <- as.name("lifetab2")
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
    if (is.null(breaks))
      breaks <- c(sort(unique(Y[,1,drop=FALSE])), Inf)
    if (attr(Y, "type") == "right" || attr(Y, "type") == 
        "counting")  {
      NA2zero <- function(x) {if (any(is.na(x))) x[is.na(x)] <- 0; x}
      temp <- tapply(1:nrow(Y), X, 
                     function(index) {
                       time <- Y[index,1,drop=FALSE]
                       cut_time <- cut(time,breaks,include.lowest=TRUE,right=FALSE)
                       event <- Y[index,2,drop=FALSE]
                       nevent <- NA2zero(tapply(event,cut_time,sum))
                       nlost <- NA2zero(tapply(event,cut_time,length)) - nevent
                       KMsurv::lifetab(tis = breaks, # should be one element longer for the intervals
                                       ninit = nrow(data),           # number of individuals at the start
                                       nlost = nlost,                 # number lost for each interval
                                       nevent = nevent)
                       },           # number of events for each interval                       
                     simplify=FALSE)
    } else {
      stop("survival type not supported")
    }
    if (length(temp)==1) {
      temp <- temp[[1]]
    }
    structure(temp, call=Call)
  }
