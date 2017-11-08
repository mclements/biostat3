addIndicators <- function(data, formula, drop.intercept = TRUE) {
  mm <- stats::model.matrix(formula,data)
  if (drop.intercept && any(intercept <- "(Intercept)" %in% names(mm)))
      mm <- mm[, !intercept, drop = FALSE]
  for (name in colnames(mm))
      data[[name]] <- mm[, name]
  data
}
