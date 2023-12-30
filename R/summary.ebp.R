# Summarizes an emdi ebp Object

#' @export
#' @importFrom moments skewness kurtosis
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom nlme predict 
#' @rdname emdi_summaries

summary.ebp <- function(object, ...) {
  throw_class_error(object, "ebp")

  call_emdi <- object$call

  N_dom_unobs <- object$framework$N_dom_unobs
  N_dom_smp <- object$framework$N_dom_smp

  N_subdom_unobs <-   object$framework$N_subdom_unobs
  N_subdom_smp <- object$framework$N_subdom_smp
  
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

  skewness_res <- skewness(residuals(object$model,
    level = 0,
    type = "pearson"
  ))
  kurtosis_res <- kurtosis(residuals(object$model,
    level = 0,
    type = "pearson"))
  variance_res <- object$model_par$sigmae2est 
  variance_ran <- object$model_par$sigmau2est 
  
  if (!is.null(object$framework$smp_subdomains) && !is.null(object$framework$pop_subdomains))
  { # two fold model
   skewness_ran <- skewness(ranef(object$model)[[object$framework$smp_domains]]) 
   kurtosis_ran <- kurtosis(ranef(object$model)[[object$framework$smp_domains]])
   skewness_ran_sub <- skewness(ranef(object$model)[[object$framework$smp_subdomains]])
   kurtosis_ran_sub <- kurtosis(ranef(object$model)[[object$framework$smp_subdomains]])
   variance_ran_sub <- object$model_par$sigmah2est 
   norm <- data.frame(
   Skewness = c(skewness_res, skewness_ran,skewness_ran_sub), 
   Kurtosis = c(kurtosis_res, kurtosis_ran,kurtosis_ran_sub),
   Shapiro_W = c(NA, NA, NA),
   Shapiro_p = c(NA, NA, NA),
   row.names = c("Error", paste0(object$framework$smp_domains," random effect"),paste0(object$framework$smp_subdomains," random effect"))
   ) 
   var <- data.frame(Variance = c(variance_res,variance_ran,variance_ran_sub),
                     row.names = c("Error", paste0(object$framework$smp_domains," random effect"),paste0(object$framework$smp_subdomains," random effect"))         
  } # close two fold model 
  else {
  skewness_ran <- skewness(ranef(object$model)$"(Intercept)")
  kurtosis_ran <- kurtosis(ranef(object$model)$"(Intercept)")
  if (length(residuals(object$model, level = 0, type = "pearson")) > 3 &&
    length(residuals(object$model, level = 0, type = "pearson")) < 5000) {
    shapiro_p_res <-
      shapiro.test(residuals(object$model, level = 0, type = "pearson"))[[2]]
    shapiro_W_res <-
      shapiro.test(residuals(object$model, level = 0, type = "pearson"))[[1]]
  } else {
    warning(strwrap(prefix = " ", initial = "",
                    "Number of observations exceeds 5000 or is lower then 3 and
                    thus the Shapiro-Wilk test is not applicable for the
                    residuals."))
    shapiro_p_res <- NA
    shapiro_W_res <- NA
  }

  if (length(ranef(object$model)$"(Intercept)") > 3 &&
    length(ranef(object$model)$"(Intercept)") < 5000) {
    shapiro_p_ran <- shapiro.test(ranef(object$model)$"(Intercept)")[[2]]
    shapiro_W_ran <- shapiro.test(ranef(object$model)$"(Intercept)")[[1]]
  } else {
    warning(strwrap(prefix = " ", initial = "",
                    "Number of domains exceeds 5000 or is lower then 3 and thus
                    the Shapiro-Wilk test is not applicable for the random
                    effects."))
    shapiro_p_ran <- NA
    shapiro_W_ran <- NA
  }

  norm <- data.frame(
    Skewness = c(skewness_res, skewness_ran),
    Kurtosis = c(kurtosis_res, kurtosis_ran),
    Shapiro_W = c(shapiro_W_res, shapiro_W_ran),
    Shapiro_p = c(shapiro_p_res, shapiro_p_ran),
    row.names = c("Error", "Random_effect")
  )
  var <- data.frame(Variance = c(variance_res,variance_ran),
                    row.names = c("Error", "Random_effect")
                    )
  } # close one fold model  
  tempMod <- object$model
  tempMod$call$fixed <- object$fixed
  r_squared <- suppressWarnings(r.squaredGLMM(tempMod))
  if (is.matrix(r_squared)) {
    r_marginal <- r_squared[1, 1]
    r_conditional <- r_squared[1, 2]
  } else {
    r_marginal <- r_squared[1]
    r_conditional <- r_squared[2]
  }
  icc_mixed <- icc(object$model)

  groups=tempMod$data[,object$framework$smp_domains]
  
  y_yhat <- data.frame("Y" = tempMod$data[,1],"marginal" = tempMod$fitted[,1], "conditional"=tempMod$fitted[,2],weights=object$model$data[,object$framework$weights])
  y_yhat_bar <- aggregate_weighted_mean(y_yhat,by=list(groups),w=object$model$data[,object$framework$weights])[,-c(1,5)]
  
  r2_area <- cor(y_yhat_bar[])[2:3,1]^2
  r_marginal_area <- r2_area[1]
  r_conditional_area <- r2_area[2]
  
  
  
  
  coeff_det <- data.frame(
    Marginal_R2    = r_marginal,
    Conditional_R2 = r_conditional,
    Marginal_Area_R2 = r_marginal_area,
    Conditional_Area_R2 = r_conditional_area, 
    row.names      = ""
  )
  
  
  

  sum_emdi <- list(
    out_of_smp = N_dom_unobs,
    in_smp = N_dom_smp,
    out_of_smp_sub = N_subdom_unobs, 
    in_smp_sub = N_subdom_smp, 
    size_smp = smp_size,
    size_pop = pop_size,
    size_dom = sizedom_smp_pop,
    smp_size_tab = NULL,
    transform = transform_method,
    normality = norm,
    variance = var, 
    icc = icc_mixed,
    coeff_determ = coeff_det,
    call = call_emdi
  )

  class(sum_emdi) <- c("summary.ebp", "emdi")
  sum_emdi
}


#' @export
print.summary.ebp <- function(x, ...) {
  throw_class_error(x, "ebp")
  cat("Empirical Best Prediction\n")
  cat("\n")
  cat("Call:\n ")
  print(x$call)
  cat("\n")
  cat("Out-of-sample domains: ", x$out_of_smp, "\n")
  cat("In-sample domains: ", x$in_smp, "\n")
  if (!is.null(x$in_smp_sub)) {
  cat("Out-of-sample subdomains: ", x$out_of_smp_sub, "\n")
  cat("In-sample subdomains: ", x$in_smp_sub, "\n")  
  }
  cat("\n")
  cat("Sample sizes:\n")
  cat("Units in sample: ", x$size_smp, "\n")
  cat("Units in population: ", x$size_pop, "\n")
  print(x$size_dom)
  cat("\n")
  if (is.null(x$call$weights)) {
    cat("Explanatory measures:\n")
  } else {
    cat("Explanatory measures for the mixed model:\n")
  }
  print(x$coeff_determ)
  cat("\n")
  if (is.null(x$call$weights)) {
    cat("Residual diagnostics:\n")
  } else {
    cat("Residual diagnostics for the mixed model:\n")
  }
  print(x$normality)
  cat("\n")
  cat("Estimated variance of random effects")
  print(x$variance)
  cat("\n")
  cat("ICC: ", x$icc, "\n")
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
  u <- as.numeric(VarCorr(model)[1, 1])
  e <- model$sigma^2
  u / (u + e)
}
