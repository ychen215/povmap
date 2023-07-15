wtd.quantile <- function (x, weights = NULL, probs = NULL) 
{
  n <- length(x)
  order <- order(x)
  x <- x[order]
  weights <- weights[order]
  if (is.null(weights)) {
    rw <- seq_len(n)/n
  }
  else {
    rw <- cumsum(weights)/sum(weights)
  }
  q <- vapply(probs, function(p) {
    if (p == 0) {
      return(x[1])
    }
    else if (p == 1) {
      return(x[n])
    }
    select <- min(which(rw >= p))
    if (rw[select] == p) {
      mean(x[select:(select + 1)])
    }
    else {
      x[select]
    }
  }, numeric(1))
  return(unname(q))
}
emdi_direct <- direct(y = "eqIncome", smp_data = eusilcA_smp,
                        +                       smp_domains = "district", weights = "weight", threshold =
                          +                           function(y, weights){0.6 * wtd.quantile(y, weights, 0.5)}, na.rm = TRUE)