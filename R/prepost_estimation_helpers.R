#' Create Descriptive Statistics for Small Area Estimation Report
#'
#' This function estimates the coefficient of variation at level specified, basic
#' statistics such number of units, regions and target areas as well as the
#' threshold on which SAE is applied and the outcome indicator of interest
#' (i.e. poverty line and poverty rate). These indicators are all expressed for
#' the census and survey
#'
#' @param ebp_object the EBP object produced from by EMDI from unit model estimation
#' the object is of class "ebp emdi"
#' @param smp_weights the sample weight variable in the household dataset
#' (i.e. the training data), include column of 1s (DO NOT LEAVE UNSPECIFIED)
#' @param pop_weights the population weight variable in the census (or synthetic)
#' census dataset (i.e. the test data), include column of 1s (NO NOT LEAVE
#' UNSPECIFIED)
#' @param repvar the variable level at which Coefficient of Variation should be
#' computed
#' @param welfare the welfare aggregate variable or outcome variable of interest
#' @param smp_data the survey/training data
#' @param pop_data the population/census/training data
#' @param threshold the poverty line or threshold specified
#' @param pop_domains the target area variable within `pop_data`
#' @param smp_domains the target area variable within `smp_data`
#' @param indicator the poverty headcount indicator in the ebp object. 
#' Defaults to "Head_count"   

#' @examples 
#' data("eusilcA_pop2")
#' data("eusilcA_smp2")
#' #### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'
#'### model estimation
#'ebp_reportdescriptives(ebp_object = emdi_model,
#'                       smp_weights = "weight",
#'                       pop_weights = "popweights",
#'                       repvar = "state",
#'                       welfare = "eqIncome",
#'                       smp_data = eusilcA_smp2, 
#'                       pop_data = eusilcA_pop2, 
#'                       threshold = 11000, 
#'                       pop_domains = "district", 
#'                       smp_domains = "district")
#'                       
#'
#'
#' @export
#'

ebp_reportdescriptives <- function(ebp_object,
                                   smp_weights,
                                   pop_weights,
                                   repvar,
                                   welfare,
                                   smp_data,
                                   pop_data,
                                   threshold,
                                   pop_domains,
                                   smp_domains,
                                   design = NULL,
                                   HT = TRUE,
                                   indicator = "Head_Count") {
  
  

  ###get list of variables
  #hh_varlist <- colnames(ebp_object$framework$smp_data)
  hh_varlist <- as.character(as.list(attributes(ebp_object$model$terms)$variables)[-1])
  pop_varlist <- hh_varlist[!(hh_varlist %in% c(welfare, smp_weights))]

  ### subset the survey and census data
  smp_df <- smp_data[complete.cases(smp_data[,hh_varlist]), 
                     c(hh_varlist, repvar, smp_weights,smp_domains)]
  pop_df <- pop_data[complete.cases(pop_data[,pop_varlist]),
                     c(pop_varlist, repvar, pop_weights,pop_domains)]
  
  
  
  ### ---------- Estimate CV for census and survey at repvar level --------- ###
  
  ##### create a dataset with headcounts, MSEs for survey and census information
  smp_doms <- unique(smp_df[[smp_domains]])
  
  ##### quickly rename column names in the MSE section of the EBP object
  colnames(ebp_object$MSE)[!grepl("Domain",
                                  colnames(ebp_object$MSE))] <-
    paste0("MSE_", colnames(ebp_object$MSE)[!grepl("Domain",
                                                   colnames(ebp_object$MSE))])
  
  df <- merge(x = ebp_object$MSE[, c("Domain", paste("MSE_",indicator,sep=""))],
              y = ebp_object$ind[, c("Domain", indicator)],
              by = "Domain")
  
  
  
  df$in_sample <- ifelse(df$Domain %in% smp_doms, 1, 0)
  
  # df$Domain <- as.integer(as.character(df$Domain))
  add_df <-
    data.frame(Domain = names(tapply(smp_df[[smp_weights]],
                                     smp_df[[smp_domains]],
                                     sum,
                                     na.rm = TRUE)),
               smp_weights = tapply(smp_df[[smp_weights]],
                                    smp_df[[smp_domains]],
                                    sum,
                                    na.rm = TRUE))
  
  df <- merge(x = df, y = add_df, by = "Domain")
  
  add_df <-
    data.frame(Domain = names(tapply(pop_df[[pop_weights]],
                                             pop_df[[pop_domains]],
                                             sum,
                                             na.rm = TRUE)),
               pop_weights = tapply(pop_df[[pop_weights]],
                                    pop_df[[pop_domains]],
                                    sum,
                                    na.rm = TRUE))
  
  df <- merge(x = df, y = add_df, by = "Domain")
  
  ### add the repvar variable to df as well
  pop_df$Domain <- pop_df[[pop_domains]]
  
  add_df <- unique(pop_df[, c("Domain", repvar)])
  
  df <- merge(x = df,
              y = add_df[, c("Domain", repvar)],
              by = "Domain")

  df$CV <- df[,paste("MSE_",indicator,sep="")] / df[,indicator]
  
  ### compute the cvs for census and survey at repvar level
  naivevar_dt <- povmap:::direct(y = welfare,
                        smp_data = ebp_object$framework$smp_data,
                        smp_domains = smp_domains,
                        weights = smp_weights,
                        threshold = threshold,
                        design = design, 
                        var = TRUE,
                        HT = HT)

  naivevar_dt$MSE$Head_Count_bench <- naivevar_dt$MSE$Head_Count
  naivevar_dt$ind$Head_Count_bench <- naivevar_dt$ind$Head_Count
  
  
  naivevar_dt$ind[,paste("Direct_",indicator,"_CV",sep="")] <- sqrt(naivevar_dt$MSE[,indicator]) / naivevar_dt$ind[,indicator]
  
  
  add_df <- data.frame(unique(df[[repvar]]))
  
  colnames(add_df) <- repvar
  
  add_df$sum_smp_weights <- tapply(X = df$smp_weights,
                                   INDEX = df[[repvar]],
                                   FUN = sum,
                                   na.rm = TRUE)
  add_df$sum_pop_weights <- tapply(X = df$pop_weights,
                                   INDEX = df[[repvar]],
                                   FUN = sum,
                                   na.rm = TRUE)
  df <- merge(x = df, y = add_df, by = repvar)
  
  povrate <- weighted.mean(x = df[,indicator],
                           w = df$smp_weights,
                           na.rm = TRUE)
  
  
  df$smp_weights <- df$smp_weights / df$sum_smp_weights
  df$pop_weights <- df$pop_weights / df$sum_pop_weights
  
  df <- merge(x = df,
              y = naivevar_dt$ind[, c("Domain", paste("Direct_",indicator,"_CV",sep=""))],
              by = "Domain")
  
  
  cv_df_region <-
    data.frame(indicator = paste0("CV for Area: ", unique(df[[repvar]])),
               ebp_cv = tapply(X = df$CV * df$pop_weights,
                               INDEX = df[[repvar]],
                               FUN = sum,
                               na.rm = TRUE),
               direct_cv = tapply(X = df[,paste("Direct_",indicator,"_CV",sep="")] * df$smp_weights,
                               INDEX = df[[repvar]],
                               FUN = sum,
                               na.rm = TRUE))
  
  #### ----------------- add other elements of the table ----------------- ####
  ##### compute number of households in census and survey
  basic_df <-
    data.frame(indicator = c("Number of Units", "Number of Regions",
                             "Number of Target Areas"),
               census = c(round(sum(pop_df[[pop_weights]], na.rm = TRUE)),
                          length(unique(pop_df[[repvar]][is.na(pop_df[[repvar]]) == FALSE])),
                          length(unique(pop_df[[pop_domains]][is.na(pop_df[[smp_domains]]) == FALSE]))),
               survey = c(ebp_object$framework$N_smp,
                          length(unique(smp_df[[repvar]][is.na(smp_df[[repvar]]) == FALSE])),
                          length(unique(smp_df[[smp_domains]][is.na(smp_df[[smp_domains]]) == FALSE]))))
  
  basic_df$census <- as.integer(basic_df$census)
  basic_df$survey <- as.integer(basic_df$survey)
  
  ##### compute poverty numbers
  smp_data$poor <- ifelse(smp_data[[welfare]] < threshold, 1, 0)
  
  smp_data[[smp_weights]] <-
    smp_data[[smp_weights]] / sum(smp_data[[smp_weights]],
                                  na.rm = TRUE)
  
  pov_df <-
    data.frame(indicator = c("National Poverty Rate", "National Poverty Line"),
               model = c(povrate, threshold),
               survey = c(sum(smp_data$poor * smp_data[[smp_weights]]),
                          threshold))
  
  row.names(cv_df_region) <- NULL
  
  return(list(cv_table = cv_df_region,
              basicinfo_df = basic_df,
              poverty_df = format(pov_df, scientific = FALSE)))
  
  
}

#' Perform test for difference between survey and census means
#'
#' This function computes weighted means of the same set of variables within the
#' census and the survey. A test for difference of the means are performed for
#' each variable with two-tailed p-values returned.
#'
#' @param smp_data the survey/training data
#' @param pop_data the population/census/training data
#' @param varlist character vector, the set of variables of interest 
#' (ensure all factor variables are converted into integers to avoid error messages)
#' @param smp_weights the sample weight variable in the household dataset
#' (i.e. the training data), include column of 1s (DO NOT LEAVE UNSPECIFIED)
#' @param pop_weights the population weight variable in the census (or synthetic)
#' census dataset (i.e. the test data), include column of 1s (NO NOT LEAVE
#' UNSPECIFIED)
#' 
#' @examples 
#' data("eusilcA_pop2")
#' data("eusilcA_smp2")
#' 
#' #### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'
#'### model estimation
#'ebp_test_means(varlist = variables,
#'               smp_weights = "weight",
#'               pop_weights = "popweights",
#'               smp_data = eusilcA_smp2, 
#'               pop_data = eusilcA_pop2)
#' 
#'
#' @export

ebp_test_means <- function(smp_data,
                           pop_data,
                           varlist,
                           smp_weights,
                           pop_weights){
  
  ### get the set of complete cases of the variables as
  ### would have been used in model estimation
  smp_df <- smp_data[complete.cases(smp_data[,c(varlist, smp_weights)]),
                     c(varlist, smp_weights)]
  pop_df <- pop_data[complete.cases(pop_data[,c(varlist, pop_weights)]),
                     c(varlist, pop_weights)]
  
  weighted.sd <- function(x, w){
    
    delta_sq <- (x - mean(x))^2 ##square deviation of the xs
    
    nzero_w <- (length(w[w > 0]) - 1) / length(w[w > 0])
    
    result <- sqrt(sum(w * (delta_sq)) / (nzero_w * sum(w)))
    
    return(result)
  }
  
  smp_means_df <- data.frame(smp_means = apply(X = smp_df[,varlist],
                                               MARGIN = 2,
                                               FUN = weighted.mean,
                                               w = smp_df[[smp_weights]]),
                             smp_sd = apply(X = smp_df[,varlist],
                                            MARGIN = 2,
                                            FUN = weighted.sd,
                                            w = smp_df[[smp_weights]]),
                             variable = varlist)
  
  pop_means_df <- data.frame(pop_means = apply(X = pop_df[,varlist],
                                               MARGIN = 2,
                                               FUN = weighted.mean,
                                               w = pop_df[[pop_weights]]),
                             pop_sd = apply(X = pop_df[,varlist],
                                            MARGIN = 2,
                                            FUN = weighted.sd,
                                            w = pop_df[[pop_weights]]),
                             variable = varlist)
  
  means_df <- merge(smp_means_df, pop_means_df, by = "variable")
  
  means_df$diff_sd <- sqrt((means_df$smp_sd)^2 + (means_df$pop_sd)^2)
  
  means_df$diff <- means_df$pop_means - means_df$smp_means
  
  means_df$zscore <- means_df$diff / means_df$diff_sd
  
  means_df$pvalue <- 2 * (1 - pnorm(abs(means_df$zscore)))
  
  return(means_df[, c("variable", "smp_means", "pop_means", "diff", "pvalue")])
  
}

#' Produce coefficient table for reporting
#'
#' This function takes the object of class 'ebp' to present the regression
#' model results having specified the number of decimal places.
#'
#' @param ebp_object the EBP object produced from by EMDI from unit model estimation
#' the object is of class "ebp emdi"
#' @param decimals the number of decimals to report on coefficient estimates
#'
#' @examples 
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#' 
#' #### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'
#'ebp_reportcoef_table(emdi_model, 4)
#'
#' @export


ebp_reportcoef_table <- function(ebp_object,
                                 decimals = 3) {
  
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) ### round to decimal place
  
  options(scipen = 999) ##drop the scientific notation
  
  varname_dt <- as.data.frame(rownames(coef(summary(ebp_object$model))))
  
  colnames(varname_dt) <- "Variable"
  
  coef_dt <- as.data.frame(coef(summary(ebp_object$model)))
  
  coef_dt <- cbind(varname_dt, coef_dt)
  
  coef_dt$Sig <- ifelse(coef_dt$`p-value` < 0.001, "***",
                        ifelse(coef_dt$`p-value` < 0.05 &
                                 coef_dt$`p-value` >= 0.001, "**",
                               ifelse(coef_dt$`p-value` < 0.01 &
                                        coef_dt$`p-value` >= 0.05, "*", "")))
  
  # coef_dt[,Sig := ifelse(`p-value` < 0.001, "***",
  #                        ifelse(`p-value` < 0.05 & `p-value` >= 0.001, "**",
  #                               ifelse(`p-value` < 0.01 & `p-value` >= 0.05, "*", "")))]
  
  coef_dt$Value <- ifelse(coef_dt$Value < abs(0.0004999999),
                          signif(coef_dt$Value, 2),
                          specify_decimal(coef_dt$Value, decimals))
  
  # coef_dt[,Value := ifelse(Value < abs(0.0004999999),
  #                          signif(Value, 2),
  #                          specify_decimal(Value, decimals))]
  
  coef_dt$StdError <- ifelse(coef_dt$Std.Error < abs(0.0004999999),
                             signif(coef_dt$Std.Error, 2),
                             specify_decimal(coef_dt$Std.Error, decimals))
  
  # coef_dt[,StdError := ifelse(Std.Error < abs(0.0004999999),
  #                              signif(Std.Error, 2),
  #                              specify_decimal(Std.Error, decimals))]
  
  coef_dt$Value <- paste0(coef_dt$Value, coef_dt$Sig)
  
  ### quick relabellings
  colnames(coef_dt)[colnames(coef_dt) %in% c("Value", "StdError")] <-
    c("coeff", "std_error")
  
  rownames(coef_dt) <- seq(nrow(coef_dt))
  
  return(coef_dt[, c("Variable", "coeff", "std_error")])
  
}

#' Produce EBP Head Count Population/Rate by Rank
#'
#' This function combines the ebp object with the census data to produce report
#' tables that rank head count estimates either by population of poor or the
#' head count rates themselves in descending order. The function allows the user
#' to select the first/last "x" number of areas by name as well.
#'
#' @param ebp_object the EBP object produced from by EMDI from unit model estimation
#' the object is of class "ebp emdi"
#' @param pop_data the population/census/training data
#' @param pop_domains the target area variable within `pop_data`
#' @param pop_domnames the population domain names
#' @param pop_weights the population weight variable in the census
#' @param byrank_indicator if argument is "count", the function ranks the product
#' of Head_Count (from object of class `ebp`) and `pop_weights`, otherwise it
#' the function simply ranks Head_Count output within `ebp` object
#' @param number_to_list an integer, the first `number_to_list` number of
#' target areas to produce from `byrank_indicator` ordering.
#' @param head a logical, if `TRUE` the top `number_to_list` results will be returned
#' and if `FALSE` the bottom `number_to_list` will be returned
#' @param indicator the indicator to rank by. Defaults to "Head_Count" 
#' 
#' @examples 
#' data("eusilcA_pop2")
#' data("eusilcA_smp2")
#' 
#' #### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'                            
#'### full data of highest population below threshold by rank (descending order)
#'ebp_report_byrank(ebp_object = emdi_model,
#'                  pop_data = eusilcA_pop2,
#'                  pop_domnames = "district",
#'                  pop_weights = "popweights")
#'
#'### full data of highest rate below threshold by rank (descending order)                   
#'ebp_report_byrank(ebp_object = emdi_model,
#'                  pop_data = eusilcA_pop2, 
#'                  pop_domains = "district",
#'                  pop_weights = "popweights",
#'                  byrank_indicator = "rate")
#'
#'### bottom 10 poverty count below threshold by rank (in ascending order)                  
#'ebp_report_byrank(ebp_object = emdi_model,
#'                  pop_data = eusilcA_pop2,
#'                  pop_domains = "district",
#'                  pop_weights = "popweights",
#'                  number_to_list = 10,
#'                  head = FALSE)
#'                    
#'                                    
#'
#' @export


ebp_report_byrank <- function(ebp_object,
                              pop_data,
                              pop_domains,
                              pop_weights,
                              byrank_indicator = "count",
                              number_to_list = NULL,
                              head = TRUE,
                              indicator = "Head_Count"){
  
  ### compute population totals
  pop_data <- pop_data[, c(pop_domains, pop_weights)]
  
  result_dt <- tapply(X = pop_data[[pop_weights]],
                      INDEX = pop_data[[pop_domains]],
                      FUN = sum,
                      na.rm = TRUE)
  
  result_dt <- as.data.frame(result_dt)
  
  result_dt <- data.frame(domain = rownames(result_dt),
                          population = result_dt[[1]])
  
  pop_data[[pop_domains]] <- as.character(pop_data[[pop_domains]])
  
  ### include the EBP Head_Count
  
  
  result_dt <- merge(x = result_dt,
                     y = ebp_object$ind[, c("Domain", indicator)],
                     by.x = "domain",
                     by.y = "Domain")

  result_dt$poor_count <- result_dt[,indicator] * result_dt$population
  
  ### rank order the table as requested
  if (byrank_indicator == "count") {
    
    result_dt <- result_dt[order(-result_dt$poor_count),]
    
  } else {
    
    result_dt <- result_dt[order(-result_dt[,indicator]),]
    
  }
  
  if (is.null(number_to_list)){
    
    number_to_list <- nrow(result_dt)
    
  }
  
  if (head == TRUE) {
    
    result_dt <- head(result_dt, number_to_list)
    
  } else if (head == FALSE) {
    
    result_dt <- tail(result_dt, number_to_list)
  }
  
  
  return(result_dt)
  
}

#' Coefficient of Variation (CV) estimations for Unit EBP Model Headcount Estimates
#' 
#' Function \code{ebp_compute_cv} estimates CVs for the headcount of the unit model
#' EBP functions using three different methods. CV, by definition, is the ratio of
#' mean square error of the head count to the  head count estimates. Therefore, the 
#' CV types are distinguished by the method of estimating the mean square. 
#' 
#' Method 1 uses the calibrated/naive bootstrapping of the MSE which allows to calibrate 
#' each bootstrap sample on auxiliary information using the \code{direct} function. 
#' Calibrated bootstrap improves on the bias of the naive bootstrap when used in the 
#' complex survey context (see \cite{Rao and Wu (1988)}) for more details.
#' 
#' Method 2 employs the Horowitz Thompson variance estimation technique to compute 
#' MSE i.e. each household is assigned the probability selection within the sample
#' under a given sampling scheme. The computation employs \code{sae::direct} function.
#'
#' Method 3 finally uses the design effect adjusted naive calibrated MSE. The design
#' effect is estimated using the \code{survey::svydesign} function.  
#' 
#' @param ebp_object the EBP object produced from by EMDI from unit model estimation
#' the object is of class "ebp emdi"
#' @param smp_data the survey/training data
#' @param welfare the welfare aggregate variable or outcome variable of interest
#' @param calibvar the calibration variable to be used in method 1
#' @param domainvar the target area variable
#' @param boot_type the bootstrap type "calibrated" or "naive" to be used in method 1
#' @param designvar the survey design variable to be used in estimating the design 
#' effect for method 3.
#' @param smp_weights the weight variable
#' @param threshold the poverty line or threshold specified
#' 
#' @examples 
#' data("eusilcA_pop2")
#' data("eusilcA_smp2")
#' 
#'#### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'                            
#'### full data of highest population below threshold by rank (descending order)
#'ebp_report_byrank(ebp_object = emdi_model,
#'                  smp_data = eusilcA_smp2,
#'                  welfare = "eqIncome",
#'                  calibvar = "state",
#'                  domainvar = "district",
#'                  threshold = 11000,
#'                  smp_weights = "weight")           
#'                  
#'                                    
#'@export
#'

ebp_compute_cv <- function(ebp_object,
                           smp_data,
                           welfare,
                           calibvar,
                           domainvar,
                           boot_type = "calibrate",
                           designvar = NULL,
                           threshold,
                           smp_weights,
                           HT= TRUE){
  
  ## ******************** Direct Estimate : Mean and CV ************************
  
  
  ## computing direct estimate using calibrated bootstrapping (EMDI + LAEKEN) - direct CV1
  #### first prepare the calibration matrix
  
  calibmatrix <- create_calibmatrix(smp_data[[calibvar]])
  
  
  direct_obj <- emdi::direct(y = welfare,
                             smp_data = ebp_object$framework$smp_data,
                             smp_domains = domainvar,
                             weights = smp_weights,
                             design = designvar,
                             threshold = threshold,
                             var = TRUE,
                             HT = HT, 
                             boot_type = boot_type,
                             X_calib = calibmatrix,
                             totals = NULL,
                             na.rm = TRUE)
  
  direct_obj$ind$Direct_Head_Count_CV <- sqrt(direct_obj$MSE$Head_Count) / direct_obj$ind$Head_Count
  
  ## computing direct estimate using the Horowitz Thompson (HT) indicator - direct CV2
  ### first compute poverty rates
  poor <- as.integer(ebp_object$framework$smp_data[[welfare]] < threshold)
  
  domsize_dt <- as.data.frame(tapply(ebp_object$model$data[[smp_weights]],
                                     ebp_object$model$data[[domainvar]],
                                     sum,
                                     na.rm = TRUE))
  
  colnames(domsize_dt) <- "popsize"
  domsize_dt$Domain <- rownames(domsize_dt)
  
  domsize_dt <- domsize_dt[is.na(domsize_dt$popsize) == FALSE,]
  
  domsize_dt <- domsize_dt[,c("Domain", "popsize")]
  
  ### HT estimator CV for direct estimate
  directht_dt <- sae::direct(y = poor,
                             dom = ebp_object$model$data[[domainvar]],
                             sweight = ebp_object$model$data[[smp_weights]],
                             domsize = domsize_dt)
  
  directht_dt$Domain <- direct_obj$ind$Domain
  ## Compute design effect controlled direct estimates and CVs. (direct CV3)
  #### first estimate naive bootstrap
  #### compute design effect
  #### include psu list into the ebp data object
  ebp_object$model$data$domainvar <- ebp_object$model$data[[domainvar]]
  
  ebp_object$model$data$poor <- as.integer(ebp_object$model$data[[welfare]] > log(threshold))
  
  ebp_object$model$data$weights <- smp_data[[smp_weights]]
  
  if(is.null(designvar)){
    
    ebpobj_svy <- survey::svydesign(ids = ~1,
                                    weights = ~weights,
                                    strata = NULL,
                                    survey.lonely.psu = "adjust",
                                    data = ebp_object$model$data)
    
  } else {
    
    ebp_object$model$data$designvar <- ebp_object$model$data[[designvar]]
    
    ebpobj_svy <- survey::svydesign(ids = ~1,
                                    weights = ~weights,
                                    strata = ~designvar,
                                    survey.lonely.psu = "adjust",
                                    data = ebp_object$model$data)
    
  }
  
  
  deff_adjust <- survey::svymean(x = ~poor, ebpobj_svy, na = TRUE, deff = TRUE)
  deff_adjust <- attr(deff_adjust, "deff")[1,1]
  
  ### multiple design effect with naive calibration
  naivevar_dt <- direct(y = welfare,
                        smp_data = ebp_object$framework$smp_data,
                        smp_domains = domainvar,
                        design = designvar,
                        weights = smp_weights,
                        threshold = threshold,
                        var = TRUE)
  
  naivevar_dt$ind$deff_CV <- 
    sqrt(deff_adjust) * (sqrt(naivevar_dt$MSE$Head_Count) / naivevar_dt$ind$Head_Count)
  
  ## ************************ SAE Model Estimates and CV Estimation ****************************
  
  ## compute standard CV using the EMDI package estimator function
  emdi_dt <- emdi::estimators(object = ebp_object,
                              indicator = "Head_Count",
                              MSE = FALSE,
                              CV = TRUE)
  
  result_dt <- emdi_dt$ind
  
  colnames(result_dt)[colnames(result_dt) %in% c("Head_Count", 
                                                 "Head_Count_CV")] <- 
    c("EBP_Head_Count", "EBP_Head_Count_CV")
  
  direct_dt <- direct_obj$ind[, c("Domain", "Head_Count", "Direct_Head_Count_CV")]
  
  colnames(direct_dt)[colnames(direct_dt) %in% c("Head_Count", 
                                                 "Direct_Head_Count_CV")] <- 
    c("Direct_Head_Count", "CB_Head_Count_CV")
  
  result_dt <- merge(result_dt, 
                     direct_dt,
                     on = "Domain")
  
  direct_dt <- directht_dt[, c("Domain", "CV")]
  
  direct_dt$CV <- direct_dt$CV / 100
  
  direct_dt$Domain <- as.factor(direct_dt$Domain)
  
  colnames(direct_dt)[colnames(direct_dt) %in% "CV"] <- "HT_Head_Count_CV"
  
  result_dt <- merge(result_dt,
                     direct_dt,
                     on = "Domain")
  
  direct_dt <- naivevar_dt$ind[, c("Domain", "deff_CV")]
  
  colnames(direct_dt)[colnames(direct_dt) %in% "deff_CV"] <- "DesignEffect_CV"
  
  result_dt <- merge(result_dt,
                     direct_dt,
                     on = "Domain")
  
  result_dt <- result_dt[,c("Domain", "Direct_Head_Count", "EBP_Head_Count", 
                            "HT_Head_Count_CV", "CB_Head_Count_CV", 
                            "DesignEffect_CV", "EBP_Head_Count_CV")]
  
  
  return(result_dt)
  
}


create_calibmatrix <- function(x){
  
  
  unique_obs <- unique(x)
  
  result <- 
  lapply(unique_obs,
         function(y) {
           
           z <- as.integer(y == x)
           
           return(z)
           
         })
  
  result <- do.call(cbind, result)
  
  colnames(result) <- unique_obs
  
  return(result)

}


#' Compare estimated and direct parameters in aggregate area levels
#' 
#' The function \code{aggregate_saedirect} computes the head count, poverty gap
#' and gini estimates at higher levels of aggregation than that estimated by
#' \code{ebp}. This serves to compare at a more representative level the 
#' model estimates as well as direct survey estimates for each of the three
#' parameters. 
#' 
#' @param ebp_object the EBP object produced from by EMDI from unit model estimation
#' the object is of class "ebp emdi"
#' @param smp_data a data.frame; the survey/training dataframe
#' @param pop_data data.frame; the population/census dataframe
#' @param welfare character; the welfare aggregate variable or outcome variable of interest
#' @param pop_domains character; the domain variable within the population/census dataframe
#' @param smp_domains character; the domain variable within the survey dataframe
#' @param smp_weights character; the sample weight variable
#' @param pop_weights character; the population weight variable
#' @param threshold numeric; the poverty line or threshold specified
#' @param indicator character; an outcome indicator of interest. Options are "Head_Count",
#' "Gini" and "Poverty_Gap"
#' @param in_sample logical, if TRUE only in-sample areas will be used in estimating the 
#' aggregated EBP model estimates. Otherwise, all areas are included in estimating the means
#' 
#' @examples 
#' data("eusilcA_pop2")
#' data("eusilcA_smp2")
#' 
#'#### set of variables used in model estimation
#'variables <- c("gender", "eqsize", "cash", "self_empl",
#'               "unempl_ben", "age_ben", "surv_ben",
#'               "sick_ben", "dis_ben", "rent", "fam_allow",
#'               "house_allow", "cap_inv", "tax_adj")
#'
#'### estimate a unit model
#'emdi_model <- emdiplus::ebp(fixed = as.formula(paste("eqIncome ~ ", paste(variables, 
#'                                                                          collapse= "+"))),
#'                            pop_data = eusilcA_pop2, 
#'                            pop_domains = "district", 
#'                            smp_data = eusilcA_smp2, 
#'                            smp_domains = "district",
#'                            na.rm = TRUE,
#'                            weights = "weight",
#'                            pop_weights = "popweights",
#'                            MSE = TRUE,
#'                            threshold = 11000,
#'                            B = 2,
#'                            L = 2)
#'
#'  ### estimate aggregated poverty headcounts (EBP and Direct) for regions                            
#' aggregate_saedirect(ebp_object = emdi_model,
#'                     smp_data = eusilcA_smp2,
#'                     pop_data = eusilcA_pop2, 
#'                     welfare = "eqIncome",
#'                     smp_domains = "district",
#'                     pop_domains = "district",
#'                     pop_weights = "popweights",
#'                     pop_regions = "state",
#'                     smp_weights = "weight",
#'                     threshold = 11000,
#'                     indicator = "Head_Count",
#'                     in_sample = FALSE,
#'                     smp_regions = "state")
#'                     
#'  ### estimate aggregated gini (EBP and Direct) for regions                            
#' aggregate_saedirect(ebp_object = emdi_model,
#'                     smp_data = eusilcA_smp2,
#'                     pop_data = eusilcA_pop2, 
#'                     welfare = "eqIncome",
#'                     smp_domains = "district",
#'                     pop_domains = "district",
#'                     pop_weights = "popweights",
#'                     pop_regions = "state",
#'                     smp_weights = "weight",
#'                     threshold = 11000,
#'                     indicator = "Gini",
#'                     in_sample = FALSE,
#'                     smp_regions = "state")
#'                     
#'  ### estimate aggregated Poverty_Gap (EBP and Direct) for regions                            
#' aggregate_saedirect(ebp_object = emdi_model,
#'                     smp_data = eusilcA_smp2,
#'                     pop_data = eusilcA_pop2, 
#'                     welfare = "eqIncome",
#'                     smp_domains = "district",
#'                     pop_domains = "district",
#'                     pop_weights = "popweights",
#'                     pop_regions = "state",
#'                     smp_weights = "weight",
#'                     threshold = 11000,
#'                     indicator = "Head_Count",
#'                     in_sample = FALSE,
#'                     smp_regions = "state")
#'                     
#'  @export


aggregate_saedirect <- function(ebp_object,
                                smp_data,
                                pop_data,
                                welfare,
                                smp_domains,
                                pop_domains,
                                pop_weights,
                                pop_regions,
                                smp_regions,
                                smp_weights,
                                threshold,
                                indicator,
                                in_sample = TRUE) {
  
  
  #### compute the small area estimates at regionvar level
  
  ### include pop_weights and pop_regions in the ebp indicator dataset
  pop_df <- pop_data[, c(pop_domains, pop_weights)]
  
  pop_df <- as.data.frame(tapply(X = pop_df[[pop_weights]],
                                 INDEX = pop_df[[pop_domains]],
                                 FUN = sum,
                                 na.rm = TRUE))
  
  colnames(pop_df) <- "pop_weights"
  
  pop_df$Domain <- rownames(pop_df)
  
  ### add indicator for in and out of sample
  ebp_object$ind$in_sample <- ifelse(ebp_object$ind$Domain %in% 
                                       unique(smp_data[,smp_domains]), 
                                     1, 
                                     0)
  
  ebp_object$ind <- merge(x = ebp_object$ind,
                          y = pop_df,
                          by = "Domain")
  
  ebp_object$ind <- merge(x = ebp_object$ind,
                          y = unique(pop_data[, c(pop_domains, 
                                                  pop_regions)]),
                          by.x = "Domain",
                          by.y = pop_domains)
  
  
  ### quick function to compute weighted mean since weighted.mean 
  ### is failing for whatever reason
  wgt.mean <- function(x, na.rm) {
    
    x <- x[1:length(x)/2]
    w <- x[(length(x)/2) + 1 : length(x)]
    
    df <- data.frame(x = x,
                     w = w)
    
    df <- df[complete.cases(df),]
    
    w <- df$w
    x <- df$x
    
    y <- sum(w * x, na.rm = na.rm) / sum(w, na.rm = na.rm)
    
    return(y)
    
  }
  
  ### regional SAE computation
  if (in_sample == TRUE){
    
    sae_df <- as.data.frame(tapply(X = c(ebp_object$ind[ebp_object$ind$in_sample == 1,
                                                        indicator],
                                         ebp_object$ind[ebp_object$ind$in_sample == 1,
                                                        "pop_weights"]),
                                   INDEX = rep(ebp_object$ind[ebp_object$ind$in_sample == 1,
                                                          pop_regions], 2),
                                   FUN = wgt.mean,
                                   na.rm = TRUE))
    
  } else {
    
    sae_df <- as.data.frame(tapply(X = c(ebp_object$ind[[indicator]], 
                                         ebp_object$ind[["pop_weights"]]),
                                   INDEX = rep(ebp_object$ind[[pop_regions]], 
                                               2),
                                   FUN = wgt.mean,
                                   na.rm = TRUE))
    
  }
  
  colnames(sae_df) <- paste0("EBP_", indicator)
  
  sae_df[[pop_regions]] <- rownames(sae_df)
  
  rownames(sae_df) <- NULL
  
  ### compute direct estimates
  call_df <- data.frame(argument = c("Head_Count", 
                                     "Gini", 
                                     "Poverty_Gap"),
                        call = c("compute_headcount", 
                                 "compute_gini", 
                                 "compute_gap"))
  
  function_call <- call_df[call_df$argument == indicator, "call"]
  
  
  if (indicator %in% c("Head_Count", "Poverty_Gap")) {
    
    direct_df <- as.data.frame(tapply(X = c(smp_data[[welfare]],
                                            smp_data[[smp_weights]]),
                                      INDEX = rep(smp_data[[smp_regions]], 2),
                                      FUN = get(function_call),
                                      threshold = threshold))
    
    
  } else {
    
    direct_df <- as.data.frame(tapply(X = c(smp_data[[welfare]], 
                                            smp_data[[smp_weights]]),
                                      INDEX = rep(smp_data[[smp_regions]], 2),
                                      FUN = get(function_call)))
    
  }
  
  direct_df <- cbind(data.frame(rownames(direct_df)), direct_df)
    
  colnames(direct_df) <- c(smp_regions, paste0("Direct_", indicator))
  
  row.names(direct_df) <- NULL
  
  ## combine results
  results_df <- merge(direct_df, sae_df, by = smp_regions)
  
  return(results_df)
  
}



#### --------  internals

compute_gini <- function(welfare) {

  
  welfare <- welfare[1:length(welfare)/2]
  weight <- welfare[(length(welfare)/2) + 1 : length(welfare)]
  
  if (is.null(weight)){
    weight <- rep(1, length(welfare))
    }

  df <- data.frame(welfare = welfare,
                   weight = weight)

  df <- df[complete.cases(df),]

  welfare <- df$welfare
  weight <- df$weight
  
  weight <- weight/sum(weight)
  
  order <- order(welfare)
  welfare <- welfare[order]
  weight <- weight[order]
  p <- cumsum(weight)
  nu <- cumsum(weight * welfare)
  n <- length(nu)
  nu <- nu/nu[n]
  gini <- sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])

  return(gini)  
}

compute_gap <- function(welfare, threshold){
  
  welfare <- welfare[1:length(welfare)/2]
  weight <- welfare[(length(welfare)/2) + 1 : length(welfare)]
  
  df <- data.frame(welfare = welfare,
                   weight = weight)
  
  df <- df[complete.cases(df),]
  
  welfare <- df$welfare
  weight <- df$weight
  
  
  pov_status <- (welfare < threshold)
  
  relative_distance <- (1 - (welfare[pov_status] / threshold))
  
  weight_pov <- weight[pov_status]
  
  weight_total <- sum(weight)
  
  fgt1 <- sum(relative_distance * weight_pov) / weight_total
  
  return(fgt1)
  
}

compute_headcount <- function(welfare, threshold){
  
  welfare <- welfare[1:length(welfare)/2]
  weight <- welfare[(length(welfare)/2) + 1 : length(welfare)]
  
  df <- data.frame(welfare = welfare,
                   weight = weight)
  
  df <- df[complete.cases(df),]
  
  welfare <- df$welfare
  weight <- df$weight
  
  pov_status <- as.integer(welfare < threshold)
  
  fgt0 <- sum(pov_status * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
  
  return(fgt0)
  
}



