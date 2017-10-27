survrate <- function (formula, data, subset, addvars = FALSE, ...) 
  {
    Call <- match.call()
    Call[[1]] <- as.name("strate")
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
    if (attr(Y, "type") == "right" || attr(Y, "type") == 
        "counting")  {
      NA2zero <- function(x) {if (any(is.na(x))) x[is.na(x)] <- 0; x}
      temp <- tapply(1:nrow(Y), X, 
                     function(index) {
                         if(attr(Y, "type") ==  "counting") {
                             T <- sum(Y[index,2]-Y[index,1], na.rm=TRUE)
                             x <- sum(Y[index,3], na.rm=TRUE)
                         } else {
                             T <- sum(Y[index,1], na.rm=TRUE)
                             x <- sum(Y[index,2], na.rm=TRUE)
                         }
                         test <- stats::poisson.test(NA2zero(x), NA2zero(T), ...)
                         out <- data.frame(x=x, T=T, rate=test$estimate, lower=test$conf.int[1], upper=test$conf.int[2])
                         category <- lapply(m[ll][index,,drop=FALSE], unique)
                         if (length(category)>0 && addvars)
                             out <- cbind(category,out)
                         out
},
                     simplify=FALSE)
      if (length(temp)>1) {
          rownames <- names(temp)
          temp <- do.call("rbind", temp)
          row.names(temp) <- rownames
      } else {
          temp <- temp[[1]]
          row.names(temp) <- 1
      }
    } else {
      stop("survival type not supported")
    }
    structure(temp, call=Call)
  }
