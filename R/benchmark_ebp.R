# Internal documentation -------------------------------------------------------

# Benchmark function for the EBP

# This function is called within the EBP-function and the agruments benchmark
# and benchmark_type are documented there.

benchmark_ebp <- function (point_estim, framework, benchmark, benchmark_type) {

  if (is.list(point_estim)) {
    estim <- as.list(point_estim$ind[names(benchmark)])
  } else {
    estim <- as.list(as.data.frame(point_estim)[names(benchmark)])
  }

  EBP_bench <- as.list(as.data.frame(
    matrix(NA, nrow = length(estim[[1]]), ncol = length(benchmark))
  ))
  names(EBP_bench) <- names(benchmark)

  share <- framework$n_pop / framework$N_pop

  for(i in names(benchmark)) {

    if (benchmark_type == "raking") {
      EBP_bench[[i]] <- estim[[i]] + benchmark[[i]] - sum(share * estim[[i]])
    } else if (benchmark_type == "ratio") {
      phi <- share / estim[[i]]
      EBP_bench[[i]] <- estim[[i]] + (1 / (sum(share^2 / phi))) *
        (benchmark[[i]] - sum(share * estim[[i]])) * (share / phi)
    }
  }

  names(EBP_bench) <- c(paste0(names(benchmark),"_bench"))

  if (is.list(point_estim)) {
    point_estim_bench <- data.frame(point_estim$ind, EBP_bench)
  } else {
    point_estim_bench <- as.matrix(data.frame(point_estim, EBP_bench))
  }

  return(point_estim_bench)
}
