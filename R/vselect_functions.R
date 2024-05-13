
lassoebp_vselect <- function(dt,
                             yvar,
                             candidate_vars,
                             lambda_min = 0,
                             lambda_max = 500,
                             by = 5,
                             domain,
                             scale = FALSE,
                             seed = 1909,
                             folds = 5,
                             family = poisson(link = "log"),
                             return_onlyvars = TRUE){

  set.seed(seed)

  N <- dim(dt)[1]

  ind <- sample(N, N)

  lambda_vector <- seq(lambda_min, lambda_max, by = by)

  nk <- floor(N / folds)

  dev_matrix <- matrix(Inf, ncol = folds, nrow = length(lambda_vector))

  ## first fit good starting model

  bic_vector <- rep(Inf, length(lambda_vector))

  ## first, fit good starting model

  pql_formula <- as.formula(paste(yvar, 1, sep = " ~ "))

  dt$domain <- dt[[domain]]

  if (scale == TRUE){

    dt[, candidate_vars] <- scale(dt[, candidate_vars])

  }


  pql <-
    glmmPQL(fixed = pql_formula,
            random = ~1 | domain,
            family = family,
            data = dt)

  delta_start <- as.matrix(t(c(as.numeric(pql$coeff$fixed),
                           rep(0, length(candidate_vars)),
                           as.numeric(t(pql$coef$random$domain)))))

  q_start <- as.numeric(VarCorr(pql)[1, 1]) ## get variance of the random effect

  ### formula to be used by glmmLasso function
  glmlasso_formula <- as.formula(paste(paste0(yvar, " ~ "),
                                       paste(candidate_vars,
                                             collapse = " + ")))


  ## loop over the folds
  for (i in 1:folds) {

    print(paste("CV Loop ", i, sep = ""))

    if (i < folds) {

      indi <- ind[(i-1) * nk + (1 : nk)]

    } else {

      indi <- ind[((i-1) * nk + 1) : N]

    }

    train_dt <- dt[-indi,]
    test_dt <- dt[indi,]

    delta_temp <- delta_start
    q_temp <- q_start

    ## loop over lambda grid
    for(j in 1:length(lambda_vector))
    {
      #print(paste("Lambda Iteration ", j,sep=""))

      glm_model <- try(glmmLasso(glmlasso_formula,
                                 rnd = list(domain = ~1),
                                 family = family,
                                 data = train_dt,
                                 lambda = lambda_vector[j],
                                 switch.NR = FALSE,
                                 final.re = FALSE,
                                 control = list(start = delta_temp[j,],
                                                q_start = q_temp[j])),
                       silent=TRUE)

      if(!inherits(glm_model, "try-error")) {

        yhat <- predict(glm_model, test_dt)
        delta_temp <- rbind(delta_temp,
                            glm_model$Deltamatrix[glm_model$conv.step,])

        q_temp <- c(q_temp, glm_model$Q_long[[glm_model$conv.step + 1]])

        dev_matrix[j,i] <- sum(family$dev.resids(test_dt[[yvar]],
                                                 yhat,
                                                 wt = rep(1,
                                                          length(yhat))))
      }
    }
  }

  dev_vector <- apply(dev_matrix, 1, sum, na.rm = TRUE)

  optlambda_index <- which.min(dev_vector)

  ## now fit full model until optimnal lambda (which is at opt4)
  for(j in 1:optlambda_index)
  {
    glm_final <- glmmLasso(fix = glmlasso_formula,
                           rnd = list(domain = ~1),
                           family = family,
                           data = dt,
                           lambda = lambda_vector[j],
                           switch.NR = FALSE,
                           final.re = FALSE,
                           control = list(start = delta_start[j,],
                                          q_start = q_start[j]))

    delta_start <- rbind(delta_start, glm_final$Deltamatrix[glm_final$conv.step,])
    q_start <- c(q_start, glm_final$Q_long[[glm_final$conv.step + 1]])
  }

  if (return_onlyvars == TRUE){

    selvars_list <- names(glm_final$coefficients[glm_final$coefficients != 0])

    selvars_list <- selvars_list[selvars_list != "(Intercept)"]

    return(selvars_list)

  }

  return(glm_final)

}




















