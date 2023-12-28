# Internal documentation -------------------------------------------------------

#quiet <- function(x) { 
#  sink(tempfile()) 
#  on.exit(sink()) 
#  invisible(force(x)) 
#} 

# Point estimation function

# This function implements the transformation of data, estimation of the nested
# error linear regression model and the monte-carlo approximation to predict
# the desired indicators. If the weighted version of the approach is used, then
# additional estimation steps are taken in order to calculate weighted
# regression coefficients before the monte-carlo approximation. See
# corresponding functions below.


point_estim <- function(framework,
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
  weights_arg <- NULL 
  if(!is.null(framework$weights) &&
     any(framework$weights_type %in% c("nlme", "nlme_lambda"))) {
    #transformation_par$transformed_data$weights_scaled <-
    #  framework$smp_data[,framework$weights] /
    #    mean(framework$smp_data[,framework$weights], na.rm = TRUE)
    transformation_par$transformed_data$weights_temp <- framework$smp_data[,framework$weights]
    #weights_arg <- nlme:::varComb(nlme:::varIdent(form = ~1),nlme:::varFixed(~1/weights_temp))
    weights_arg <- ~1/weights_temp
  }   
  
  random_arg <- NULL 
  if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
    # Do two fold model 
    random_arg <- list(as.formula(~1),as.formula(~1))
    names(random_arg) <- c(framework$smp_domains,framework$smp_subdomains)  
    
      }
  else {
  # Do one-fold model 
  random_arg[framework$smp_domains] <- list(as.formula(~1))
  } 
  
  
  args <- list(fixed=fixed,
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
               weights = weights_arg)
  
  #ptm <- proc.time()
  mixed_model <- do.call(nlme:::lme,args)
  #proc.time() - ptm     
 
    
    # Function model_par extracts the needed parameters theta from the nested
    # error linear regression model. It returns the beta coefficients (betas),
    # sigmae2est, sigmau2est, sigmah2est (which is set to 0 for one-fold models) and the random effect (rand_eff)
  
    est_par <- model_par(
      mixed_model = mixed_model,
      framework = framework,
      fixed = fixed,
      transformation_par = transformation_par
    )
    
  
  

     
    
        
    
   
   





    dep_var <- transformation_par$transformed_data[[
      as.character(mixed_model$terms[[2]])]]
       # this is needed for guadarrama weights calculation in gen_model 


  # Function gen_model calculates the parameters in the generating model.
  # See Molina and Rao (2010) p. 375 (20)
  # The function returns sigmav2est and the constant part mu.
  gen_par <- gen_model(
    model_par = est_par,
    fixed = fixed,
    framework = framework,
    dep_var = dep_var 
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


# All following functions are only internal ------------------------------------

# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).

model_par <- function(framework,
                      mixed_model,
                      fixed,
                      transformation_par) {

  # fixed parameters
  betas <- nlme::fixed.effects(mixed_model)
  # Estimated error variance
  sigmae2est <- mixed_model$sigma^2
  # VarCorr(fit2) is the estimated random error variance
  # if subdomains specified, save the results in sigmae2est, sigmah2est, and sigmau2est 
  # Random effect: vector with zeros for all domains, filled with 0
  rand_eff <- rep(0, framework$N_dom_pop)  
  rand_eff_h <- rep(0, framework$N_subdom_pop)  
  if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
    sigmau2est <- as.numeric(VarCorr(mixed_model)[2,1])
    sigmah2est <- as.numeric(VarCorr(mixed_model)[4,1])
    # Random sub-area effect for sample subdomains
    rand_eff[framework$dist_obs_dom] <- random.effects(mixed_model,level=1)$"(Intercept)"[framework$dist_obs_smp_dom]
    # Random area-effect 
    rand_eff_h[framework$dist_obs_subdom] <- random.effects(mixed_model,level=2)$"(Intercept)"[framework$dist_obs_smp_subdom]
  } 
  else {
    sigmau2est <- as.numeric(nlme::VarCorr(mixed_model)[1, 1])  
    sigmah2est <- 0 
    # random area effect 
    rand_eff[framework$dist_obs_dom] <- (random.effects(mixed_model)[[1]])
  }
  
  
  
  
 
  #Variance of cluster components 
    # random effect for in-sample domains (dist_obs_dom)
  
  
    
    
    return(list(
      betas = betas,
      sigmae2est = sigmae2est,
      sigmah2est = sigmah2est, 
      sigmau2est = sigmau2est,
      rand_eff = rand_eff,
      rand_eff_h = rand_eff_h
    ))

} # End model_par



  
  
  
# Function gen_model calculates the parameters in the generating model with a one-fold model.
# See Molina and Rao (2010) p. 375 (20)
gen_model <- function(fixed,
                      framework,
                      model_par,
                      dep_var) {

  if (is.null(framework$weights)) {
    weight_smp <- rep(1,framework$N_smp)  
  } 
  else {
    weight_smp <- framework$smp_data[[as.character(framework$weights)]]
  }
    
  
  # calculate gamma 
  if (any(framework$weights_type %in% c("nlme", "nlme_lambda")) | is.null(framework$weights)) {
    rand_eff <- model_par$rand_eff
    rand_eff_h <- model_par$rand_eff_h 
    
  weight_sum <- rep(0, framework$N_dom_smp)
  sums <- aggregate(data.frame(weight_smp,weight_smp^2), by=list(framework$smp_domains_vec),FUN=sum)
  delta2 <- sums[,3] / sums[,2]^2 # sum of the squares divided by the square of the sum 
  gamma <- model_par$sigmau2est / (model_par$sigmau2est + ((model_par$sigmae2est + model_par$sigmah2est) * delta2))
  if (model_par$sigmah2est>0) {
    sums_sub <- aggregate(data.frame(weight_smp, weight_smp^2)[framework$smp_subdomains_vec,], by = list(framework$smp_subdomains_vec), FUN = sum)
    sums_sub <- sums_sub[framework$dist_obs_smp_subdom,]
    delta2_sub <- sums_sub[,3] / sums_sub[,2]^2
    gamma_sub <- model_par$sigmah2est / (model_par$sigmah2est + model_par$sigmae2est * delta2_sub)
  }
  
  
  } 
  else {
      # Calculations needed for pseudo EB for Guadarrama option 
    betas <- model_par$betas
      weight_sum <- rep(0, framework$N_dom_smp)
      mean_dep <- rep(0, framework$N_dom_smp)
      mean_indep <- matrix(0, nrow = framework$N_dom_smp, ncol = length(betas))
      delta2 <- rep(0, framework$N_dom_smp)
      
      gamma <- rep(0, framework$N_dom_smp)
      num <- matrix(0, nrow = length(betas), ncol = 1)
      den <- matrix(0, nrow = length(betas), ncol = length(betas))
      
      for (d in 1:framework$N_dom_smp) {
        #cat(d)
        domain <- names(table(framework$smp_domains_vec)[d])
        dep_smp <- dep_var[
            framework$smp_domains_vec == domain
          ]
        # Domain means of of the dependent variable
        dep_smp <- dep_var[framework$smp_domains_vec == domain]
        
        
        
        
        weight_smp <- framework$smp_data[[
          as.character(framework$weights)]][framework$smp_domains_vec == domain]
        weight_sum[d] <- sum(weight_smp)
        
         if(length(weight_smp) == 1) {
           indep_smp <- matrix(model.matrix(fixed, framework$smp_data)[framework$smp_domains_vec == domain,]
                               , ncol = length(betas), nrow = 1)
        } else {
          indep_smp <- model.matrix(fixed, framework$smp_data)[framework$smp_domains_vec == domain,]
        }
        
        # weighted mean of the dependent variable
        mean_dep[d] <- sum(weight_smp * dep_smp) / weight_sum[d]
        
        # weighted means of the auxiliary information
        for (k in 1:length(betas)) {
          mean_indep[d, k] <- sum(weight_smp * indep_smp[, k]) / weight_sum[d]
        }
        
        delta2[d] <- sum(weight_smp^2) / (weight_sum[d]^2)
        gamma[d] <- model_par$sigmau2est / (model_par$sigmau2est + model_par$sigmae2est * delta2[d])
        weight_smp_diag <- diag(x=weight_smp,nrow=length(weight_smp), ncol=length(weight_smp))
        dep_var_ast <- dep_smp - gamma[d] * mean_dep[d]
        indep_weight <- t(indep_smp) %*% weight_smp_diag
        indep_var_ast <- indep_smp - matrix(rep(
          gamma[d] *
            mean_indep[d, ],
          framework$n_smp[d]
        ),
        nrow = framework$n_smp[d],
        byrow = TRUE
        )
        
        
        
        num <- num + (indep_weight %*% dep_var_ast)
        den <- den + (indep_weight %*% indep_var_ast)
      } # close loop over domains 
      
      
      betas <- solve(den) %*% num
      
      # random effect for in-sample domains (dist_obs_dom)
      rand_eff <- rep(0, framework$N_dom_pop)  
      rand_eff[framework$dist_obs_dom] <- gamma * (mean_dep -
                                                            mean_indep %*% betas)
      
    
    
    
  }
 

  
  framework$pop_data[[paste0(fixed[2])]] <- seq_len(nrow(framework$pop_data))
  X_pop <- model.matrix(fixed, framework$pop_data)
  
  # Constant part of predicted y
  mu_fixed <- X_pop %*% model_par$betas
  sigmav2est <- model_par$sigmau2est * (1 - gamma)
  rand_eff_pop <- rep(rand_eff, framework$n_pop)
  
   
   
  if (model_par$sigmah2est==0) {
    mu <- mu_fixed + rand_eff_pop
    sigmai2est <- 0 
  }
  else {
    sigmai2est <- model_par$sigmah2est * (1-gamma_sub)
    rand_eff_h_pop <- rep(rand_eff_h,framework$n_pop_subdom)
    mu <- mu_fixed + rand_eff_pop + rand_eff_h_pop 
  }
  
    
    
    return(list(sigmav2est = sigmav2est, sigmai2est = sigmai2est, mu = mu, mu_fixed = mu_fixed,rand_eff=rand_eff))
} # End gen_model


#calculate weighted mean using R's aggregate function 
aggregate_weighted_mean  <-function(df,by,w) {
  wdf <- cbind(df*w,w)  
  wdf_sum <- aggregate(wdf,by=by, FUN=sum)
  aggregate_weighted_mean <- cbind(Domain = wdf_sum[,1],wdf_sum[,-1]/wdf_sum$w)
  aggregate_weighted_mean <-within(aggregate_weighted_mean,rm(w))
  return(aggregate_weighted_mean)
}

# calculated weighted quantile 
aggregate_weighted_quantile  <-function(df,by,w,q) {
   aggregate_weighted_quantile <- mapply(FUN=wtd.quantile,x=split(df,by),weight=split(w,pop_domains_vec_tmp),probs=q)
}
  


#analytic functions 
# This function calculates the expected value, intially in the case of no transformation 
analytic <- function(transformation,
                     framework,
                     lambda,  
                     shift, 
                     model_par, 
                     gen_model, 
                     fixed, 
                     Ydump
                     ) {

  if(!is.null(framework$aggregate_to_vec)){
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  } else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }

  if(!is.null(framework$pop_weights)) {
    pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
  }else{
    pop_weights_vec <- rep(1, framework$N_pop)
  }
  
  
  


  
  

# construct vector for variance of random effect, copied from errors_gen line 568
sigma2vu <- vector(length = framework$N_pop)
# variance of random effect for out-of-sample domains
sigma2vu[!framework$obs_dom] <- model_par$sigmau2est
# variance of random effect for in-sample domains
sigma2vu[framework$obs_dom] <- rep(gen_model$sigmav2est,framework$n_pop[framework$dist_obs_dom])
#sigma2vu[framework$obs_dom] <- rep(gen_model$sigmav2est,framework$n_pop[framework$dist_obs_dom])

if (model_par$sigmah2est==0) {
  sigma2eta <-0 
} 
  else {
    sigma2eta <- vector(length=framework$N_pop)
    # variance of subarea random effect for out-of-sample domains    
    sigma2eta[!framework$obs_subdom] <- model_par$sigmah2est 
    sigma2eta[framework$obs_subdom] <- rep(gen_model$sigmai2est,framework$n_pop_subdom[framework$dist_obs_subdom])
  }


var <- sigma2vu+sigma2eta+model_par$sigmae2est



# do mean with function 
indicators <- matrix(ncol=length(framework$indicator_names),nrow=framework$N_pop)
colnames(indicators) <- framework$indicator_names

indicators[,"Mean"] <- expected_transformed_mean(gen_model$mu,var=var, transformation=transformation,lambda=lambda) 
indicators[,"Head_Count"] <- expected_head_count(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,threshold=framework$threshold)
#indicators[,"Median"] <- transformed_percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.5)
#indicators[,"Quantile_10"] <- transformed_percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.1) 
#indicators[,"Quantile_25"] <- transformed_percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.25) 
#indicators[,"Quantile_75"] <- transformed_percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.75) 
#indicators[,"Quantile_90"] <- transformed_percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.9) 



point_estimates <- aggregate_weighted_mean(indicators[,c("Mean","Head_Count")],by=list("Domain" = pop_domains_vec_tmp),w=pop_weights_vec) 
point_estimates <- cbind(point_estimates, data.frame(matrix(ncol=length(framework$indicator_names)-2,nrow=N_dom_pop_tmp)))
colnames(point_estimates) <- c("Domain",framework$indicator_names)
return(point_estimates)
} # end analytic 

expected_transformed_mean <- function(mu=mu,var=var,transformation=transformation,lambda=lambda) {
  if (transformation=="no") {
    expected_mean <- mu   
  }
  else if (transformation=="arcsin") {
    term1 <- arcsin_transform_back(mu) # y = sin(x)^2,  
    # dy/dx = 2(sin(x)*cos(x))  
    dy2dx <- -2*sin(mu)^2+2*cos(mu)^2 # product rule for differentiation, f=2sin(x),g=cos(x),f'=2cos(x),g'=-sin(x),dy2dx=f*g'+g*f'=-2sin(x)^2+2*cos(X)^2   
    #dy3dx = -4*sin(mu)*cos(m)-4*cos(x)*sin(x)=-8*sin(x)*cos(x)=-4*dy/dyx
    #dy4/dx <- -4*dy2dx, dy5/dx=-4*dy3dx = 32*sin(x)*cos(x)=16*dydx, dy6/dx <- 16*dy2dx, dy8/dx <- -64*dy2dx
    expected_mean <- term1+0.5*dy2dx*(var)-(1/24)*-4*dy2dx*3*(var^2)+(1/720)*-16*dy2dx*(var^3)+(1/40320)*-64*dy2dx*(var^4)  # 9th order Taylor series
  }
  else if (transformation=="log") {
    expected_mean <- exp(mu+(0.5*var))
  }
  else if (transformation=="log.shift") {
  expected_mean <- exp(mu+(0.5*var))-lambda
  }
  
    return(expected_mean)
}

expected_head_count <- function(mu=mu,threshold=threshold,var=var,transformation=transformation,lambda=lambda,shift=shift) {
  # do poverty by transforming threshold, then applying normal CDF. 
  transformed_threshold <- transformation(y=threshold,transformation=transformation,lambda=lambda,shift=shift)
  expected_head_count <- pnorm(transformed_threshold$y - mu,sd=var^0.5) # formula for head count 
  return(expected_head_count)
}


transformed_percentile <- function(mu=mu,threshold=threshold,var=var,transformation=transformation,lambda=lambda,shift=shift,p=p) {
   mu_percentile <- mu+qnorm(p,0,var)
   transformed_percentile <- back_transformation(y=mu_percentile,transformation=transformation,lambda=lambda,shift=shift) 
   return(transformed_percentile) 
}

# Monte-Carlo approximation ----------------------------------------------------

# The function approximates the expected value (Molina and Rao (2010)
# p.372 (6)). For description of monte-carlo simulation see Molina and
# Rao (2010) p. 373 (13) and p. 374-375
monte_carlo <- function(transformation,
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
    Ydumpdf <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(Ydumpdf) <- c("L","Domain","Simulated_Y","XBetahat","eta","epsilon")
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
    errors <- errors_gen(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model
    )

    # Prediction of population vector y
    # See below for function prediction_y.
    population_vector <- prediction_y(
      transformation = transformation,
      lambda = lambda,
      shift = shift,
      gen_model = gen_model,
      errors_gen = errors,
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


# The function errors_gen returns error terms of the generating model.
# See Molina and Rao (2010) p. 375 (20)

errors_gen <- function(framework, model_par, gen_model) {
  epsilon <- rnorm(framework$N_pop, 0, sqrt(model_par$sigmae2est))

  # empty vector for new random effect in generating model
  vu <- vector(length = framework$N_pop)
  # new random effect for out-of-sample domains
  vu[!framework$obs_dom] <- rep(
    rnorm(
      framework$N_dom_unobs,
      0,
      sqrt(model_par$sigmau2est)
    ),
    framework$n_pop[!framework$dist_obs_dom]
  )
  # new random effect for in-sample-domains
  vu[framework$obs_dom] <- rep(
    rnorm(
      rep(1, framework$N_dom_smp),
      0,
      sqrt(gen_model$sigmav2est)
    ),
    framework$n_pop[framework$dist_obs_dom]
  )
  # individual error term in generating model epsilon

  return(list(epsilon = epsilon, vu = vu))
} # End errors_gen

# The function prediction_y returns a predicted income vector which can be used
# to calculate indicators. Note that a whole income vector is predicted without
# distinction between in- and out-of-sample domains.
prediction_y <- function(transformation,
                         lambda,
                         shift,
                         gen_model,
                         errors_gen,
                         framework,
                         fixed) {

  # predicted population income vector
  y_pred <- gen_model$mu + errors_gen$epsilon + errors_gen$vu

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
