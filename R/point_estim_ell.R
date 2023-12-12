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
 
  # Using do.call passes the name of the weight vector from framework to the plm function 
  # if weight is NULL, that is appropriately passed to PLM 
 args <- list(formula=fixed, 
           data = transformation_par$transformed_data, 
           weights = framework$smp_data[,framework$weights],
           model="random",
           index = framework$smp_domains)
 
 re_model <- do.call(plm, args)
 
  # Function model_par extracts the needed parameters theta from the random
  # effects linear regression model. It returns the beta coefficients (betas),
  # sigmae2est, and sigmau2est.
  est_par <- model_par_ell (
    re_model = re_model,
    framework = framework,
    fixed = fixed,
    transformation_par = transformation_par
  )
  
  
  
  
  # Monte-Carlo approximation --------------------------------------------------
  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(
        y =
          as.numeric(framework$smp_data[[paste0(fixed[2])]])
      )
  }
  
  # Monte Carlo ELL 
  
  # The monte-carlo function returns a data frame of desired indicators.
  if (L>0) {
    indicator_prediction <- monte_carlo_ell(
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
      stop(strwrap(prefix = " ", initial = "",
                   "L must be positive when using ELL estimation"))
    }
        
  
  return(list(
    ind = indicator_prediction,
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model_par = est_par,
    gen_model = gen_par,
    model = re_model
  ))
} # End point estimation function





# All following functions are only internal ------------------------------------

# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).

model_par_ell <- function(framework,
                      re_model,
                      fixed,
                      transformation_par) {
  
  # fixed parameters
  betas <- re_model$coefficients 
  # Estimated error variance
  sigmae2est <- re_model$ercomp$sigma2[1]
  # VarCorr(fit2) is the estimated random error variance
  sigmau2est <- re_model$ercomp$sigma2[2]
  vcov <- re_model$vcov
  
  

    
    return(list(
      betas = betas,
      sigmae2est = sigmae2est,
      sigmau2est = sigmau2est,
      vcov = vcov
    ))
  } 

# Monte-Carlo approximation ----------------------------------------------------

# The function approximates the expected value (Molina and Rao (2010)
# p.372 (6)). For description of monte-carlo simulation see Molina and
# Rao (2010) p. 373 (13) and p. 374-375
monte_carlo_ell <- function(transformation,
                        L,
                        framework,
                        lambda,
                        shift,
                        model_par,
                        gen_model,
                        fixed,
                        Ydump) {
  
  # Preparing matrices for indicators for the Monte-Carlo simulation
  
  if(!is.null(framework$aggregate_to_vec)){
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  } else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }
  
  
  if (!is.null(Ydump)) {
    Ydumpdf <- data.frame(matrix(ncol = 8+length(model_par$betas), nrow = 0))
    colnames(Ydumpdf) <- c("L","Domain","Simulated_Y","XBetahat","eta","epsilon","sigma2eta","sigma2eps",colnames(model_par$betas))
    write.csv(Ydumpdf,Ydump,row.names = FALSE)
  }
  
  
  ests_mcmc <- array(dim = c(
    N_dom_pop_tmp,
    L,
    length(framework$indicator_names)
  ))
  
  for (l in seq_len(L)) {
    
    # Errors in generating model: individual error term and random effect
    # See below for function errors_gen.
    parameters <- parameters_gen(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model
    )
    
    browser()
    # generate gen_model$mu here based parameters$beta, if you have the data conveniently around 
    
    # Prediction of population vector y
    # See below for function prediction_y.
    population_vector <- prediction_y(
      transformation = transformation,
      lambda = lambda,
      shift = shift,
      gen_model = gen_model,
      parameters_gen = parameters,
      framework = framework,
      fixed = fixed
    )
    
    if(!is.null(framework$pop_weights)){
      pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
    }else{
      pop_weights_vec <- rep(1, nrow(framework$pop_data))
    }
    if (!is.null(Ydump)){
      Ydumpdf <- data.frame(rep(l,nrow(framework$pop_data)), framework$pop_domains_vec,population_vector,gen_model$mu,errors$vu,errors$epsilon)
      #write.csv(Ydumpdf,Ydump,row.names = FALSE,append=TRUE)
      write.table(Ydumpdf,file=Ydump,row.names = FALSE,append=TRUE,col.names=F, sep=",") 
    }
    
    # Calculation of indicators for each Monte Carlo population
    ests_mcmc[, l, ] <-
      matrix(
        nrow = N_dom_pop_tmp,
        data = unlist(lapply(framework$indicator_list,
                             function(f, threshold) {
                               matrix(
                                 nrow = N_dom_pop_tmp,
                                 data = unlist(mapply(
                                   y = split(population_vector, pop_domains_vec_tmp),
                                   pop_weights = split(pop_weights_vec, pop_domains_vec_tmp),
                                   f,
                                   threshold = framework$threshold
                                 )), byrow = TRUE
                               )
                             },
                             threshold = framework$threshold
        ))
      )
  } # End for loop
  
  
  # Point estimations of indicators by taking the mean
  
  point_estimates <- data.frame(
    Domain = unique(pop_domains_vec_tmp),
    apply(ests_mcmc, c(3), rowMeans)
  )
  colnames(point_estimates) <- c("Domain", framework$indicator_names)
  return(point_estimates)
} # End Monte-Carlo


# The function parameters_gen returns error terms of the generating model.
# and draws betas 
# See Molina and Rao (2010) p. 375 (20)

parameters_gen <- function(framework, model_par, gen_model) {
  epsilon <- rnorm(framework$N_pop, 0, sqrt(model_par$sigmae2est))
  
  # empty vector for new random effect in generating model
  vu <- vector(length = framework$N_pop)
  # draw random effect
  vu <- rep(rnorm(framework$N_dom_pop,0,sqrt(model_par$sigmau2est)
    ),framework$n_pop)
  

    #betas <- mvrnorm(n=1, mu=model_par$betas,Sigma=as.matrix(model_par$vcov))
    # This is equivalent and needs no package   
    R <- chol(model_par$vcov)   
    betas <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + model_par$betas 
    
  return(list(epsilon = epsilon, vu = vu, betas=betas))
} # End parameters_gen

# The function prediction_y returns a predicted income vector which can be used
# to calculate indicators. Note that a whole income vector is predicted without
# distinction between in- and out-of-sample domains.
prediction_y <- function(transformation,
                         lambda,
                         shift,
                         gen_model,
                         parameters_gen,
                         framework,
                         fixed) {
  # construct vector of Xes by starting with intercept and binding it to the population
  # dataframe constructed from all the coefficient names except the 
  Xes <- cbind(rep(1,framework$N_pop),framework$pop_data[rownames(parameters_gen$betas)[-1]]) 
  
  # predicted population income vector
  y_pred <- as.matrix(Xes) %*% (parameters_gen$betas) + parameters_gen$epsilon + parameters_gen$vu
  
  # back-transformation of predicted population income vector
  y_pred <- back_transformation(
    y = y_pred,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  y_pred[!is.finite(y_pred)] <- 0
  
  return(y_pred)
} # End prediction_y
