# Point estimation function

# This function implements the transformation of data, estimation of the unconditional random effects model 
# and the monte-carlo approximation to predict
# the desired indicators. If the weighted version of the approach is used, then
# additional estimation steps are taken in order to calculate weighted
# regression coefficients before the monte-carlo approximation. See
# corresponding functions below.


point_estim_ell <- function(framework,
                        fixed,
                        transformation,
                        interval,
                        L,
                        keep_data = FALSE,
                        Ydump = NULL) {
  
  # Transformation of data -----------------------------------------------------
  
  # Estimating the optimal parameter by optimization
  # Optimal parameter function returns the minimum of the optimization
  # functions from generic_opt; the minimum is the optimal lambda.
  # The function can be found in the script optimal_parameter.R
  optimal_lambda <- optimal_parameter(
    generic_opt = generic_opt,
    fixed = fixed,
    smp_data = framework$smp_data,
    smp_domains = framework$smp_domains,
    transformation = transformation,
    interval = interval,
    framework = framework
  )
  
  # Data_transformation function returns transformed data and shift parameter.
  # The function can be found in the script transformation_functions.R
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = optimal_lambda
  )
  shift_par <- transformation_par$shift
  
  # Model estimation, model parameter and parameter of generating model --------
  
  # Estimation of the nested error linear regression model
  # See Molina and Rao (2010) p. 374
  # lme function is included in the nlme package which is imported.

  # Do an unconditional random effect model 
  random_arg <- NULL 
  random_arg[framework$smp_domains] <- list(as.formula(~1))
  names(random_arg) <- c(framework$smp_domains)
  weights_vec <- transformation_par$transformed_data[,framework$weights]
  browser()
  
  re_model <- plm(
    formula = fixed,
    #formula = pchinc ~ wall_improved, 
    data = transformation_par$transformed_data,
    weights = "FACTOR", 
    model = "random",
    index = framework$smp_domains
  )
  
    
    
  # Function model_par extracts the needed parameters theta from the nested
  # error linear regression model. It returns the beta coefficients (betas),
  # sigmae2est, sigmau2est and the random effect (rand_eff).
  est_par <- model_par_ell (
    mixed_model = mixed_model,
    framework = framework,
    fixed = fixed,
    transformation_par = transformation_par
  )
  
  
  
  if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
    # Do two fold model 
    random_arg <- list(as.formula(~1),as.formula(~1))
    names(random_arg) <- c(framework$smp_domains,framework$smp_subdomains)
    mixed_model2f <- nlme::lme(
      fixed = fixed,
      data = transformation_par$transformed_data,
      random = random_arg, 
      method = framework$nlme_method,
      control = nlme::lmeControl(maxIter = framework$nlme_maxiter,
                                 tolerance = framework$nlme_tolerance,
                                 opt = framework$nlme_opt,
                                 optimMethod = framework$nlme_optimmethod, 
                                 msMaxIter=framework$nlme_msmaxiter,
                                 msTol=framework$nlme_mstol,
                                 returnObject = framework$nlme_returnobject 
      ),
      keep.data = keep_data,
      weights = quiet(cat(weights_arg))
    )
    est_par$sigma2u2f <- as.numeric(VarCorr(mixed_model2f)[2])
    est_par$sigma2h2f <- as.numeric(VarCorr(mixed_model2f)[4])
    est_par$sigma2e2f <- mixed_model2f$sigma^2 
  } 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Function gen_model calculates the parameters in the generating model.
  # See Molina and Rao (2010) p. 375 (20)
  # The function returns sigmav2est and the constant part mu.
  gen_par <- gen_model(
    model_par = est_par,
    fixed = fixed,
    framework = framework
  )
  
  
  
  
  
  
  # Monte-Carlo approximation --------------------------------------------------
  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(
        y =
          as.numeric(framework$smp_data[[paste0(fixed[2])]])
      )
  }
  
  # The monte-carlo function returns a data frame of desired indicators.
  if (L>0) {
    indicator_prediction <- monte_carlo(
      transformation = transformation,
      L = L,
      framework = framework,
      lambda = optimal_lambda,
      shift = shift_par,
      model_par = est_par,
      gen_model = gen_par,
      fixed = fixed,
      Ydump = Ydump 
    )
  }
  else {
    indicator_prediction <- analytic(
      transformation=transformation,
      framework = framework,
      lambda = optimal_lambda,
      shift = shift_par,
      model_par = est_par,
      gen_model = gen_par, 
      fixed = fixed, 
      Ydump=Ydump
    )
  }
  
  
  mixed_model$coefficients_weighted <- if (!is.null(framework$weights)) {
    as.numeric(est_par$betas)
  } else {
    NULL
  }
  names(mixed_model$coefficients_weighted) <- if (!is.null(framework$weights)) {
    rownames(est_par$betas)
  } else {
    NULL
  }
  return(list(
    ind = indicator_prediction,
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model_par = est_par,
    gen_model = gen_par,
    model = mixed_model
  ))
} # End point estimation function