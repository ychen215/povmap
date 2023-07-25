library(povmap)
data("eusilcA_smp")
data('eusilcA_pop')
emdi_model_Guadarrama <- ebp(
   fixed = eqIncome ~ gender + eqsize + cash + self_empl +
     unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
     fam_allow + house_allow + cap_inv + tax_adj,
  pop_data = eusilcA_pop, pop_domains = "district",
   smp_data = eusilcA_smp, smp_domains = "district", weights = "weight",
   weights_type = "Guadarrama", transformation = "log", na.rm = TRUE)

emdi_model_nlme_log <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
    fam_allow + house_allow + cap_inv + tax_adj,
   pop_data = eusilcA_pop, pop_domains = "district",
  smp_data = eusilcA_smp, smp_domains = "district", weights = "weight",
 weights_type = "nlme", transformation = "log", na.rm = TRUE
  )

emdi_model_nlme_bc <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
   fam_allow + house_allow + cap_inv + tax_adj,
   pop_data = eusilcA_pop, pop_domains = "district",
   smp_data = eusilcA_smp, smp_domains = "district",
  weights = "weight", weights_type = "nlme", na.rm = TRUE
  )

emdi_model_nlme_lambda <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash + self_empl +
     unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
     fam_allow + house_allow + cap_inv + tax_adj,
   pop_data = eusilcA_pop, pop_domains = "district",
   smp_data = eusilcA_smp, smp_domains = "district",
   weights = "weight", weights_type = "nlme_lambda", na.rm = TRUE
   )

head(estimators(emdi_model_Guadarrama, indicator = "Mean"))

benchmark <- mean(eusilcA_smp$eqIncome)
names(benchmark) <- c("Mean")
ebp_bench_external <- ebp(
       fixed = eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben +
             age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow +
             house_allow + cap_inv + tax_adj,
       pop_data = eusilcA_pop, pop_domains = "district",
       smp_data = eusilcA_smp, smp_domains = "district",
       na.rm = TRUE, benchmark = benchmark, benchmark_type = "ratio")
median_state <- tapply(eusilcA_smp$eqIncome, eusilcA_smp$state, median)
benchmark_table <- data.frame(state = names(median_state), Mean = median_state)

ebp_bench_external_state <- ebp(
   fixed = eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben +
     age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow +
     house_allow + cap_inv + tax_adj,
   pop_data = eusilcA_pop, pop_domains = "district",
   smp_data = eusilcA_smp, smp_domains = "district",
   na.rm = TRUE, benchmark = benchmark_table, benchmark_type = "ratio",
   benchmark_level = "state")

#Takes longer because of MSE estimation
ebp_bench_internal_state <- ebp(
   fixed = eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben +
     age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow +
     house_allow + cap_inv + tax_adj,
   pop_data = eusilcA_pop, pop_domains = "district",
   smp_data = eusilcA_smp, smp_domains = "district",
   weights = "weight", weights_type = "nlme",
   na.rm = TRUE, benchmark = c("Mean"), benchmark_type = "ratio",
   benchmark_level = "state", MSE = TRUE)

head(estimators(ebp_bench_external, indicator = "Mean_bench"))

sum(ebp_bench_external$ind$Mean_bench *
      table(ebp_bench_external$framework$pop_domains_vec)/
      ebp_bench_external$framework$N_pop)

mean(eusilcA_smp$eqIncome)

head(estimators(ebp_bench_external, indicator = "Mean_bench"))

# Takes longer due to use of bootstrap 
emdi_direct <- direct(
  y = "eqIncome", smp_data = eusilcA_smp, smp_domains = "district",
  weights = "weight", var = TRUE, boot_type = "naive", B = 50, na.rm = TRUE)
emdi_direct$ind$Mean_bench <- emdi_direct$ind$Mean
emdi_direct$MSE$Mean_bench <- emdi_direct$MSE$Mean
compare_plot(ebp_bench_internal_state, direct = emdi_direct,
             CV = TRUE, indicator = "Mean_bench")


ebp_no <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
    fam_allow + house_allow + cap_inv + tax_adj,
  pop_data = eusilcA_pop, pop_domains = "district",
  smp_data = eusilcA_smp, smp_domains = "district",
  na.rm = TRUE, transformation = "no"
  )

ebp_ordernorm <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
    fam_allow + house_allow + cap_inv + tax_adj,
  pop_data = eusilcA_pop, pop_domains = "district",
  smp_data = eusilcA_smp, smp_domains = "district",
  na.rm = TRUE, transformation = "ordernorm"
  )

summary(ebp_no)$normality
summary(ebp_ordernorm)$normality
qqnorm(ebp_no)
qqnorm(ebp_ordernorm)
eusilcA_smp$eqIncome_prop <- eusilcA_smp$eqIncome / max(eusilcA_smp$eqIncome)

ebp_no <- ebp(
  fixed = eqIncome_prop ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
    fam_allow + house_allow + cap_inv + tax_adj,
  pop_data = eusilcA_pop, pop_domains = "district",
  smp_data = eusilcA_smp, smp_domains = "district",
  na.rm = TRUE, transformation = "no")

ebp_arcsin <- ebp(
  fixed = eqIncome_prop ~ gender + eqsize + cash + self_empl +
    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
    fam_allow + house_allow + cap_inv + tax_adj,
  pop_data = eusilcA_pop, pop_domains = "district",
  smp_data = eusilcA_smp, smp_domains = "district",
  transformation = "arcsin", na.rm = TRUE)

summary(ebp_no)$normality
summary(ebp_arcsin)$normality

qqnorm(ebp_no)
qqnorm(ebp_arcsin)
