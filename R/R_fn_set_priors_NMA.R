






#' R_fn_compute_n_thr_random_NMA
#' @keywords internal
#' @export
R_fn_compute_n_thr_random_NMA <- function(n_thr, 
                                          n_studies_per_test) { 
      
      n_thr_random <- 0
      for (t in 1:n_index_tests) {
        n_thr_random <- n_thr_random + n_thr[t]*n_studies_per_test[t]
      }
      ##
      return(n_thr_random)
  
}



#' R_fn_set_priors_NMA
#' @keywords internal
#' @export
R_fn_set_priors_NMA <- function(  priors,
                                  ##
                                  intercept_only,
                                  cov_data,
                                  ##
                                  n_studies,
                                  n_index_tests,
                                  ##
                                  n_studies_per_test,
                                  ##
                                  cts,
                                  ##
                                  model_parameterisation,
                                  random_thresholds,
                                  Dirichlet_random_effects_type,
                                  ##
                                  softplus,
                                  ##
                                  n_thr,
                                  n_cat,
                                  ##
                                  indicator_index_test_in_study = NULL
)  {
    
    X <- cov_data$X
    baseline_case_nd <- cov_data$baseline_case_nd
    baseline_case_d  <- cov_data$baseline_case_d
    baseline_case    <- cov_data$baseline_case
  
    cov_info_list <- R_fn_get_covariate_info_NMA( intercept_only = intercept_only,
                                                  cov_data = cov_data, 
                                                  model_parameterisation = model_parameterisation, 
                                                  n_index_tests = n_index_tests,
                                                  n_studies = n_studies)
    ##
    n_covariates_max <- cov_info_list$n_covariates_max
    print(paste("n_covariates_max = ", cov_info_list$n_covariates_max))
    ##
    priors$use_probit_link <- if_null_then_set_to(priors$use_probit_link, 1)
    if (priors$use_probit_link == 1) { 
      mult <- 1.0
    } else { 
      mult <- 1.702
    }
    ##
    priors$K_fold_CV_indicator <- if_null_then_set_to(priors$K_fold_CV_indicator, rep(1, n_studies))
    ##
    priors$compound_symmetry <- if_null_then_set_to(priors$compound_symmetry, 0)
    priors$hetero_sigma <- if_null_then_set_to(priors$hetero_sigma, 0)
    ##
    n_thr_max <- max(n_thr)
    n_cat_max <- max(n_cat)
    ##
    if (cts == TRUE) { 
            
                    # ##
                    # ## ---- Set priors for the locations ("beta"):
                    # ##
                    # priors$prior_beta_mu_mean <- if_null_then_set_to( priors$prior_beta_mu_mean, 
                    #                                                   rep(list( matrix(0.0, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    # priors$prior_beta_mu_SD   <- if_null_then_set_to( priors$prior_beta_mu_SD, 
                    #                                                   rep(list( matrix(5.0*mult, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    # ##
                    # priors$prior_beta_tau_SD <- if_null_then_set_to(priors$prior_beta_tau_SD, 
                    #                                                 array(1.0*mult, dim = c(n_index_tests, 2)))
                    # ##
                    # priors$prior_beta_sigma_SD <- if_null_then_set_to(priors$prior_beta_sigma_SD, 
                    #                                                   rep(1.0*mult, 2))
                    # ##
                    # ## ---- Set priors for raw scales ("gamma"):
                    # ##
                    # if (softplus == TRUE) {
                    #       priors$prior_raw_scale_mu_mean <- if_null_then_set_to( priors$prior_raw_scale_mu_mean, 
                    #                                                              rep(list( matrix(0.0, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    #       priors$prior_raw_scale_mu_SD   <- if_null_then_set_to( priors$prior_raw_scale_mu_SD,
                    #                                                              rep(list( matrix(1.0*mult, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    #       ##
                    #       priors$prior_raw_scale_tau_SD <-  if_null_then_set_to( priors$prior_raw_scale_tau_SD, 
                    #                                                              array(1.0*mult, dim = c(n_index_tests, 2)))
                    #       ##alpha_lb
                    #       priors$prior_raw_scale_sigma_SD <- if_null_then_set_to( priors$prior_raw_scale_sigma_SD, 
                    #                                                               rep(1.0*mult, 2))
                    # } else if (softplus == FALSE) {
                    #       priors$prior_raw_scale_mu_mean <-  if_null_then_set_to( priors$prior_raw_scale_mu_mean, 
                    #                                                               rep(list( matrix(0.0, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    #       priors$prior_raw_scale_mu_SD   <-  if_null_then_set_to( priors$prior_raw_scale_mu_SD,
                    #                                                               rep(list( matrix(1.0*mult, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                    #       ##
                    #       priors$prior_raw_scale_tau_SD <- if_null_then_set_to( priors$prior_raw_scale_tau_SD, 
                    #                                                             array(1.0*mult, dim = c(n_index_tests, 2)))
                    #       ##
                    #       priors$prior_raw_scale_sigma_SD <- if_null_then_set_to( priors$prior_raw_scale_sigma_SD,
                    #                                                               rep(1.0*mult, 2))
                    # }
                    # ##
                    # ## ---- Set priors for box-cox:
                    # ##
                    # priors$prior_boxcox_lambda_mean <- if_null_then_set_to( priors$prior_boxcox_lambda_mean, 
                    #                                                         rep(0.0, n_index_tests))
                    # priors$prior_boxcox_lambda_SD   <- if_null_then_set_to( priors$prior_boxcox_lambda_SD,   
                    #                                                         rep(1.0, n_index_tests))
                    # ##
                    # ## ---- Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                    # ##
                    # priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                    # priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                    # priors$prior_beta_corr_LKJ <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                    # ##
                    # ## ---- Set default priors for the between-study corr matricex for raw_scale and "raw_scale_L_Omega"):
                    # ##
                    # priors$raw_scale_corr_lb <- if_null_then_set_to(priors$raw_scale_corr_lb, -1.0)
                    # priors$raw_scale_corr_ub <- if_null_then_set_to(priors$raw_scale_corr_ub, +1.0)
                    # priors$prior_raw_scale_corr_LKJ <- if_null_then_set_to(priors$prior_raw_scale_corr_LKJ, 2.0)
      
    } else if (cts == FALSE) { ## ordinal
      
                    ####
                    if (model_parameterisation %in% c("Xu", "bivariate")) {
                                ##
                                ## ---- Set priors for the locations ("beta"):
                                ##
                                priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean, 
                                                                                 rep(list( matrix(0.0, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                                priors$prior_beta_mu_SD   <- if_null_then_set_to(priors$prior_beta_mu_SD, 
                                                                                 rep(list( matrix(1.0*mult, nrow = 2, ncol = n_covariates_max) ), n_index_tests))
                                for (t in 1:n_index_tests) {
                                  try({ 
                                    priors$prior_beta_mu_SD[[t]][1:2, 2:n_covariates_max] <- 0.5*mult
                                  })
                                }
                                ##
                                priors$prior_beta_tau_SD <- if_null_then_set_to(priors$prior_beta_tau_SD, 
                                                                                array(dim = c(n_index_tests, 2), 0.5*mult))
                                ##
                                priors$prior_beta_sigma_SD <- if_null_then_set_to(priors$prior_beta_sigma_SD, 
                                                                                  rep(0.5*mult, 2))
                                
                                ##
                                ## ---- Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                                ##
                                priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                                priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                                priors$prior_beta_corr_LKJ <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                                ##
                                ## ---- Induced-Dirichlet priors:
                                ##
                                if (random_thresholds == TRUE) { 
                                        
                                        n_thr_random <- R_fn_compute_n_thr_random_NMA( n_thr = n_thr,
                                                                                       n_studies_per_test = n_studies_per_test)
                                        priors$n_total_C_if_random <- n_thr_random
                                        ##
                                        ## array[n_index_tests, 2] vector[n_cat_max] prior_dirichlet_phi;
                                        mat <- array(1.0, dim = c(n_index_tests, n_thr_max + 1))
                                        priors$prior_dirichlet_phi <- if_null_then_set_to( priors$prior_dirichlet_phi, 
                                                                                           list(mat, mat))
                                        ##
                                        ##
                                        prior_kappa_mean_nd <- rep(log(200), n_index_tests)
                                        prior_kappa_mean_d  <- rep(log(50), n_index_tests)
                                        priors$prior_kappa_mean   <- if_null_then_set_to(  priors$prior_kappa_mean, 
                                                                                           list(prior_kappa_mean_nd, prior_kappa_mean_d))
                                        ##
                                        prior_kappa_SD <- rep(0.75, n_index_tests)
                                        priors$prior_kappa_SD   <- if_null_then_set_to(  priors$prior_kappa_SD, 
                                                                                         list(prior_kappa_SD, prior_kappa_SD))
                                        ##
                                        priors$kappa_lb <- if_null_then_set_to( priors$kappa_lb, 
                                                                                1.0)
                                        ##
                                        priors$kappa_ub <- if_null_then_set_to( priors$kappa_ub, 
                                                                                exp(10))
                                        ##
                                        priors_kappa_df <- 15
                                        priors$prior_kappa_df <- if_null_then_set_to( priors$prior_kappa_df, 
                                                                                      list(rep(priors_kappa_df, n_index_tests), 
                                                                                           rep(priors_kappa_df, n_index_tests)))
                                } else { 
                                        priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, 
                                                                                            array(1.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                }
                                
                    } else if (model_parameterisation %in% c("R&G", "HSROC", "Gatsonis")) { 
                                # ##
                                # ## ---- Set priors for the locations ("beta"):
                                # ##
                                # vec <- rep(0.0, n_covariates_max)
                                # priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean, 
                                #                                                  rep(list(vec), n_index_tests))
                                # ##
                                # vec <- rep(1.0*mult, n_covariates_max)
                                # priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD, 
                                #                                                rep(list(vec), n_index_tests))
                                # ##
                                # priors$prior_beta_tau_SD   <- c(if_null_then_set_to(priors$prior_beta_tau_SD, 
                                #                                                     array(dim = c(n_index_tests, 1), 1.0*mult)))
                                # ##
                                # priors$prior_beta_sigma_SD <- c(if_null_then_set_to(priors$prior_beta_sigma_SD, 
                                #                                                     rep(1.0*mult, 1)))
                                # ##
                                # ## ---- Priors for raw_scale ("gamma"):
                                # ##
                                # if (softplus == TRUE) {
                                #       vec <- rep(0.0, n_covariates_max)
                                #       priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                #                                                             rep(list(vec), n_index_tests))
                                #       ##
                                #       vec <- rep(1.0*mult, n_covariates_max)
                                #       priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                #                                                           rep(list(vec), n_index_tests))
                                #       ##
                                #       priors$prior_raw_scale_tau_SD   <-  if_null_then_set_to(priors$prior_raw_scale_tau_SD,  
                                #                                                                rep(1.0*mult, n_index_tests))
                                #       ##
                                #       priors$prior_raw_scale_sigma_SD <-  if_null_then_set_to(priors$prior_raw_scale_sigma_SD, 
                                #                                                                1.0*mult)
                                # } else { 
                                #       vec <- rep(0.0, n_covariates_max)
                                #       priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                #                                                             rep(list(vec), n_index_tests))
                                #       ##
                                #       vec <- rep(0.5, n_covariates_max)
                                #       priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                #                                                           rep(list(vec), n_index_tests))
                                #       ##
                                #       priors$prior_raw_scale_tau_SD   <- if_null_then_set_to(priors$prior_raw_scale_tau_SD,  
                                #                                                              rep(0.5*mult, n_index_tests))
                                #       ##
                                #       priors$prior_raw_scale_sigma_SD <- if_null_then_set_to(priors$prior_raw_scale_sigma_SD, 
                                #                                                              0.5*mult)
                                # }
                                # ##
                                # ## ---- Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                                # ##
                                # priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                                # priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                                # priors$prior_beta_corr_LKJ  <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                                # 
                                # if (random_thresholds == TRUE) { 
                                #         n_thr_random <- R_fn_compute_n_thr_random_NMA( n_thr = n_thr,
                                #                                                        n_studies_per_test = n_studies_per_test)
                                #         priors$n_total_C_if_random <- n_thr_random
                                #         ##
                                #         
                                #         ##
                                #         ##
                                #         prior_dirichlet_phi <- array(1.0, dim = c(n_index_tests, n_thr_max + 1))
                                #         priors$prior_dirichlet_phi <- if_null_then_set_to( priors$prior_dirichlet_phi, 
                                #                                                            prior_dirichlet_phi)
                                #         ##
                                #         prior_kappa_mean <- rep(log(50), n_index_tests)
                                #         priors$prior_kappa_mean   <- if_null_then_set_to(  priors$prior_kappa_mean, 
                                #                                                            prior_kappa_mean)
                                #         ##
                                #         prior_kappa_SD <- rep(1.5, n_index_tests)
                                #         priors$prior_kappa_SD   <- if_null_then_set_to(  priors$prior_kappa_SD, 
                                #                                                          prior_kappa_SD)
                                #         ##
                                #         priors$kappa_lb <- if_null_then_set_to( priors$kappa_lb, 
                                #                                                 0.001)
                                # }
                                # 
                    }
      
    }
    ##
    ## output list:
    ##
    return(priors)
    
}





 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

