




 


#' R_fn_set_priors_NMA
#' @keywords internal
#' @export
R_fn_set_priors_NMA <- function(  priors,
                                  ##
                                  cts,
                                  ##
                                  model_parameterisation,
                                  random_thresholds,
                                  Dirichlet_random_effects_type,
                                  ##
                                  softplus,
                                  ##
                                  n_index_tests,
                                  ##
                                  n_thr,
                                  n_cat
)  {
  
  
    # ##
    # ## Update priors (if necessary):
    # ##
    # n_tests <- n_index_tests
    ##
    if (cts == TRUE) { 
            
                    ##
                    ## Set priors for the locations ("beta"):
                    ##
                    priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean, 
                                                                     0.0, array(dim = c(n_index_tests, 2)))
                    priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD, 
                                                                     5.0, array(dim = c(n_index_tests, 2)))
                    ##
                    priors$prior_beta_tau_SD <- if_null_then_set_to(priors$prior_beta_tau_SD, 
                                                                    1.0, array(dim = c(n_index_tests, 2)))
                    ##
                    priors$prior_beta_sigma_SD <- if_null_then_set_to(priors$prior_beta_sigma_SD, 
                                                                     rep(1.0, 2))
                    ##
                    ## Set priors for raw scales ("gamma"):
                    ##
                    if (softplus == TRUE) {
                          priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                                                           array(dim = c(n_index_tests, 2), 0.0))
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                                                         array(dim = c(n_index_tests, 2), 5.0))
                          ##
                          priors$prior_raw_scale_tau_SD <- if_null_then_set_to(priors$prior_raw_scale_tau_SD, 
                                                                          array(dim = c(n_index_tests, 2), 1.0))
                          ##
                          priors$prior_raw_scale_sigma_SD <- if_null_then_set_to(priors$prior_raw_scale_sigma_SD, 
                                                                            rep(1.0, 2))
                    } else if (softplus == FALSE) {
                          priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                                                                array(dim = c(n_index_tests, 2), 0.0))
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                                                              array(dim = c(n_index_tests, 2), 2.5))
                          ##
                          priors$prior_raw_scale_tau_SD <- if_null_then_set_to(priors$prior_raw_scale_tau_SD, 
                                                                               array(dim = c(n_index_tests, 2), 1.0))
                          ##
                          priors$prior_raw_scale_sigma_SD <- if_null_then_set_to(priors$prior_raw_scale_sigma_SD, 
                                                                                 rep(1.0, 2))
                    }
                    ## Set priors for box-cox:
                    priors$prior_boxcox_lambda_mean <- if_null_then_set_to(priors$prior_boxcox_lambda_mean, rep(0.0, n_index_tests))
                    priors$prior_boxcox_lambda_SD   <- if_null_then_set_to(priors$prior_boxcox_lambda_SD,   rep(1.0, n_index_tests))
                    ##
                    ## Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                    ##
                    priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                    check_vec_length(priors, "beta_corr_lb", 1)
                    ##
                    priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                    check_vec_length(priors, "beta_corr_ub", 1)
                    ##
                    priors$prior_beta_corr_LKJ      <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                    check_vec_length(priors, "prior_beta_corr_LKJ", 1)
                    ##
                    ## Set default priors for the between-study corr matricex for raw_scale and "raw_scale_L_Omega"):
                    ##
                    priors$raw_scale_corr_lb <- if_null_then_set_to(priors$raw_scale_corr_lb, -1.0)
                    check_vec_length(priors, "raw_scale_corr_lb", 1)
                    ##
                    priors$raw_scale_corr_ub <- if_null_then_set_to(priors$raw_scale_corr_ub, +1.0)
                    check_vec_length(priors, "raw_scale_corr_ub", 1)
                    ##
                    priors$prior_raw_scale_corr_LKJ <- if_null_then_set_to(priors$prior_raw_scale_corr_LKJ, 2.0)
                    check_vec_length(priors, "prior_raw_scale_corr_LKJ", 1)
      
    } else if (cts == FALSE) { ## ordinal
      
                    ####
                    if (model_parameterisation == "Xu") {
                                ##
                                ## Set priors for the locations ("beta"):
                                ##
                                priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean, 
                                                                                 array(dim = c(n_index_tests, 2), 0.0))
                                priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD, 
                                                                               array(dim = c(n_index_tests, 2), 1.0))
                                ##
                                priors$prior_beta_tau_SD <- if_null_then_set_to(priors$prior_beta_tau_SD, 
                                                                                array(dim = c(n_index_tests, 2), 1.0))
                                ##
                                priors$prior_beta_sigma_SD <- if_null_then_set_to(priors$prior_beta_sigma_SD, 
                                                                                  rep(1.0, 2))
                                ##
                                ## Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                                ##
                                priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                                check_vec_length(priors, "beta_corr_lb", 1)
                                ##
                                priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                                check_vec_length(priors, "beta_corr_ub", 1)
                                ##
                                priors$prior_beta_corr_LKJ      <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                                check_vec_length(priors, "prior_beta_corr_LKJ", 1)
                                ##
                                ## Induced-Dirichlet priors:
                                ##
                                priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, 
                                                                                    array(1.0, dim = c(n_index_tests, max(n_thr))))
                     }
                    
      
    }
    ##
    ## output list:
    ##
    return(priors)
    
 
  
 
  
}





 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

