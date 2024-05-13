#### A quick script to show how to develop poverty maps using package data and the
#### the post estimation diagnostics of interest

#### start with some model selection
pacman::p_load(glmmLasso, nlme, MASS)

eusilcA_smp$logeqIncome <- log(eusilcA_smp$eqIncome)


eusilcA_smp$fake1 <- rnorm(nrow(eusilcA_smp)) + eusilcA_smp$logeqIncome*2
eusilcA_smp$fake2 <- rnorm(nrow(eusilcA_smp), mean = 2, sd = 3) - eusilcA_smp$eqIncome

eusilcA_pop$fake1 <- rnorm(nrow(eusilcA_pop)) + mean(eusilcA_smp$logeqIncome*2, na.rm = TRUE)
eusilcA_pop$fake2 <- rnorm(nrow(eusilcA_pop), mean = 2, sd = 3)- mean(eusilcA_smp$eqIncome,
                                                                      na.rm = TRUE)

candidate_vars <- colnames(eusilcA_smp)[!grepl("Income|state|district|weight|gender",
                                               colnames(eusilcA_smp))]

eusilcA_smp$male <- ifelse(eusilcA_smp$gender == "male", 1, 0)

### create more variables that obviously work better


### try soccer data to see why the function might be failing
vselect_obj <-
lassoebp_vselect(dt = eusilcA_smp,
                 yvar = "logeqIncome",
                 domain = "district",
                 candidate_vars = candidate_vars,
                 scale = TRUE,
                 family = gaussian(link = "identity"),
                 return_onlyvars = TRUE)


#### now lets try to create a poverty map and some post estimation diagnostics

model_formula <- as.formula(paste("eqIncome ~ ",
                                  paste(vselect_obj,
                                        collapse = " + ")))

unit_model <- povmap::ebp(fixed = model_formula,
                          pop_data = eusilcA_pop,
                          smp_data = eusilcA_smp,
                          smp_domains = "district",
                          pop_domains = "district",
                          L = 2,
                          threshold = 0.5*mean(eusilcA_pop$eqIncome, na.rm = TRUE),
                          transformation = "log",
                          MSE = TRUE,
                          B = 3,
                          weights = "weight",
                          pop_weights = "hhsize")




#### quick post estimation diagnostics and tables that will go into your paper
direct_est <- direct(y = "eqIncome",
                     smp_data = eusilcA_smp,
                     smp_domains = "district",
                     weights = "weight",
                     var = TRUE,
                     B = 2)

descriptives_dt <-
ebp_reportdescriptives(model = unit_model,
                       direct = direct_est,
                       smp_data = eusilcA_smp,
                       weights = "weight",
                       pop_weights = "hhsize",
                       CV_level = "state",
                       pop_data = eusilcA_pop,
                       pop_domains = "district")

### test that we have the same distribution of variables
means_dt <-
  ebp_test_means(varlist = candidate_vars,
                 pop_data = eusilcA_pop,
                 smp_data = eusilcA_smp,
                 weights = "weight")

results_dt <-
ebp_reportcoef_table(unit_model, 4)

# full data of highest population below threshold by rank (descending order)
ebp_report_byrank(model = unit_model,
                  pop_data = eusilcA_pop,
                  pop_domains = "district",
                  pop_weights = "hhsize")

# full data of highest rate below threshold by rank (descending order)
ebp_report_byrank(model = unit_model,
                  pop_data = eusilcA_pop,
                  pop_domains = "district",
                  pop_weights = "hhsize",
                  byrank_indicator = "rate")

# bottom 10 poverty count below threshold by rank (in ascending order)
ebp_report_byrank(model = unit_model,
                  pop_data = eusilcA_pop,
                  pop_domains = "district",
                  pop_weights = "hhsize",
                  number_to_list = 10,
                  head = FALSE)


#### compute the different kinds of CV
cv_dt <-
  ebp_compute_cv(model = unit_model,
                 threshold = 0.5*mean(eusilcA_pop$eqIncome, na.rm = TRUE))

normality_fit <-
  ebp_normalityfit(model = unit_model)





