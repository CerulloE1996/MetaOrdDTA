

## This model below is NOT fine! need to re-sample - use less chains
Model_D_real_data <- readRDS(file.path(
                            "application_results_real_data",
                            "seed_1_dummy_data_0intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))

## This model below is fine (all main r-hat's < 1.01 w/ high 400+ min_ESS!)
## This is the model used for the data simulation / DGM for the dummy data!
Model_FULL_ALL_4_COVARIATES_UN_RAND_CUTPOINTS_real_data <- readRDS(file.path(
  "application_results_real_data",
  "seed_1_dummy_data_0intercept_only_0N_covs_7param_Xurand_thr_TRUECS_0het_sigma_0_MR_model_5_applied_results.RDS"))


summary_main <- Model_D_real_data$model_summary_and_trace_obj$get_summary_main()
summary_main





Model_D_real_data$model_summary_and_trace_obj$get_HMC_info()


model_summary_and_trace_obj <- Model_FULL_ALL_4_COVARIATES_UN_RAND_CUTPOINTS_real_data$model_summary_and_trace_obj



##
## ---- For 1st table:
##
{
  # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "AUC"))) %>% print(n = 100)
  ##
  tibble_AUC <- model_summary_and_trace_obj$extract_params(params = c("AUC"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_AUC_pred <- model_summary_and_trace_obj$extract_params(params = c("AUC_pred"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_AUC_diff <- model_summary_and_trace_obj$extract_params(params = c("AUC_diff"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  ##
  tibble_Se_baseline  <- model_summary_and_trace_obj$extract_params(params = c("Se_baseline"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_Sp_baseline  <- model_summary_and_trace_obj$extract_params(params = c("Sp_baseline"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}
##
## ---- For 2nd table:
##
{
  tibble_beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu"))  %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_total_SD_inc_C <- model_summary_and_trace_obj$extract_params(params = c("total_SD_inc_C")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_sigma <- model_summary_and_trace_obj$extract_params(params = c("beta_sigma")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_tau <- model_summary_and_trace_obj$extract_params(params = c("beta_tau")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_corr <- model_summary_and_trace_obj$extract_params(params = c("beta_corr")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_rho12 <- model_summary_and_trace_obj$extract_params(params = c("rho12")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}
