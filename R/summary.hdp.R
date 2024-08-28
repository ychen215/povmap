# Summarizes a povmap hdp Object
#' @export
#' @rdname povmap_summaries

summary.hdp <- function(object, ...) {
  #throw_class_error(object, "hdp") #I commented this line since it shows error
  #that "hdp object has to be created by the emdi". However, hdp function is not
  #part of emdi.

  call_povmap <- object$call

  N_dom_unobs <- object$framework$N_dom_unobs
  N_dom_smp <- object$framework$N_dom_smp

  smp_size <- object$framework$N_smp
  pop_size <- object$framework$N_pop

  smp_size_dom <- summary(as.data.frame(
    table(object$framework$smp_domains_vec)
  )[, "Freq"])
  pop_size_dom <- summary(as.data.frame(
    table(object$framework$pop_domains_vec)
  )[, "Freq"])
  sizedom_smp_pop <- rbind(
    Sample_domains = smp_size_dom,
    Population_domains = pop_size_dom
  )

  if (object$transformation == "box.cox" || object$transformation == "dual") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Method = object$method,
      Optimal_lambda =
        object$transform_param$optimal_lambda,
      Shift_parameter =
        round(object$transform_param$shift_par, 3),
      row.names = ""
    )
  } else if (object$transformation == "log.shift") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Method = object$method,
      Optimal_lambda =
        object$transform_param$optimal_lambda,
      row.names = ""
    )
  } else if (object$transformation == "log") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Shift_parameter =
        round(object$transform_param$shift_par, 3),
      row.names = ""
    )
  } else if (object$transformation == "ordernorm") {
    transform_method <- data.frame(Transformation  = object$transformation,
                                   Shift_parameter = 0,
                                   row.names       = ""
    )
  } else if (object$transformation == "arcsin") {
    transform_method <- data.frame(Transformation  = object$transformation,
                                   Shift_parameter = 0,
                                   row.names       = ""
    )
  } else if (object$transformation == "no") {
    transform_method <- NULL
  }


  icc_mixed <- summary(icc(object$model))

  sum_povmap <- list(
    out_of_smp = N_dom_unobs,
    in_smp = N_dom_smp,
    size_smp = smp_size,
    size_pop = pop_size,
    size_dom = sizedom_smp_pop,
    smp_size_tab = NULL,
    transform = transform_method,
    #normality = norm,
    icc = icc_mixed,
    # coeff_determ = coeff_det,
    call = call_povmap
  )

  class(sum_povmap) <- c("summary.hdp", "povmap")
  sum_povmap
}


#' @export
print.summary.hdp <- function(x, ...) {
  #throw_class_error(x, "hdp")
  cat("Empirical Best Prediction with high-dimensional parameter\n")
  cat("\n")
  cat("Call:\n ")
  print(x$call)
  cat("\n")
  cat("Out-of-sample domains: ", x$out_of_smp, "\n")
  cat("In-sample domains: ", x$in_smp, "\n")
  cat("\n")
  cat("Sample sizes:\n")
  cat("Units in sample: ", x$size_smp, "\n")
  cat("Units in population: ", x$size_pop, "\n")
  print(x$size_dom)
  cat("\n")
  # cat("ICC: ", x$icc, "\n")
  cat("ICC:\n")
  print(x$icc)
  cat("\n")
  if (is.null(x$transform)) {
    cat("Transformation: No transformation \n")
  } else {
    cat("Transformation:\n")
    print(x$transform)
  }
}


#  ICC

icc <- function(model) {
  u <- as.numeric(model$sigma2u)
  e <- model$sigma2e
  return(u / (u + e))
}
