




 


#' R_fn_set_priors_NMA
#' @keywords internal
#' @export
R_fn_set_priors_NMA <- function(  priors,
                                  ##
                                  n_studies,
                                  n_index_tests,
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
                    priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean, array(0.0, dim = c(n_index_tests, 2)))
                    priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD, array(5.0, dim = c(n_index_tests, 2)))
                    ##
                    priors$prior_beta_tau_SD <- if_null_then_set_to(priors$prior_beta_tau_SD, array(1.0, dim = c(n_index_tests, 2)))
                    ##
                    priors$prior_beta_sigma_SD <- if_null_then_set_to(priors$prior_beta_sigma_SD, rep(1.0, 2))
                    ##
                    ## Set priors for raw scales ("gamma"):
                    ##
                    if (softplus == TRUE) {
                          priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, array(0.0, dim = c(n_index_tests, 2)))
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, array(5.0, dim = c(n_index_tests, 2)))
                          ##
                          priors$prior_raw_scale_tau_SD <- if_null_then_set_to(priors$prior_raw_scale_tau_SD, array(1.0, dim = c(n_index_tests, 2)))
                          ##alpha_lb
                          priors$prior_raw_scale_sigma_SD <- if_null_then_set_to(priors$prior_raw_scale_sigma_SD, rep(1.0, 2))
                    } else if (softplus == FALSE) {
                          priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, array(0.0, dim = c(n_index_tests, 2)))
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, array(2.5, dim = c(n_index_tests, 2)))
                          ##
                          priors$prior_raw_scale_tau_SD <- if_null_then_set_to(priors$prior_raw_scale_tau_SD, array(1.0, dim = c(n_index_tests, 2)))
                          ##
                          priors$prior_raw_scale_sigma_SD <- if_null_then_set_to(priors$prior_raw_scale_sigma_SD, rep(1.0, 2))
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
                    if (model_parameterisation %in% c("Xu", "bivariate")) {
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
                                priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, array(1.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                ##
                                n_total_pooled_cat  <- sum(n_cat)
                                n_thr_random = n_thr * n_studies
                                n_total_C_if_random <- sum(n_thr_random)
                                ##
                                priors$n_total_C_if_random <- n_total_C_if_random
                                priors$n_total_pooled_cat  <- n_total_pooled_cat
                                ##
                                priors$alpha_lb         <- if_null_then_set_to(priors$alpha_lb, 1.0)
                                priors$prior_alpha_mean <- if_null_then_set_to(priors$prior_alpha_mean, array(0.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                priors$prior_alpha_SD   <- if_null_then_set_to(priors$prior_alpha_SD,   array(10.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                
                    } else if (model_parameterisation %in% c("R&G", "HSROC", "Gatsonis")) { 
                      
                                ##
                                ## Set priors for the locations ("beta"):
                                ##
                                priors$prior_beta_mu_mean  <- c(if_null_then_set_to(priors$prior_beta_mu_mean, array(dim = c(n_index_tests, 1), 0.0)))
                                priors$prior_beta_mu_SD    <- c(if_null_then_set_to(priors$prior_beta_mu_SD, array(dim = c(n_index_tests, 1), 1.0)))
                                ##
                                priors$prior_beta_tau_SD   <- c(if_null_then_set_to(priors$prior_beta_tau_SD, array(dim = c(n_index_tests, 1), 1.0)))
                                ##
                                priors$prior_beta_sigma_SD <- c(if_null_then_set_to(priors$prior_beta_sigma_SD, rep(1.0, 1)))
                                ##
                                ## Priors for raw_scale ("gamma"):
                                ##
                                if (softplus == TRUE) {
                                      priors$prior_raw_scale_mu_mean  <- c(if_null_then_set_to(priors$prior_raw_scale_mu_mean, array(dim = c(n_index_tests, 1), 0.0)))
                                      priors$prior_raw_scale_mu_SD    <- c(if_null_then_set_to(priors$prior_raw_scale_mu_SD,   array(dim = c(n_index_tests, 1), 1.0)))
                                      priors$prior_raw_scale_tau_SD   <- c(if_null_then_set_to(priors$prior_raw_scale_tau_SD,  array(dim = c(n_index_tests, 1), 1.0)))
                                      priors$prior_raw_scale_sigma_SD <- c(if_null_then_set_to(priors$prior_raw_scale_sigma_SD, rep(1.0, 1)))
                                } else { 
                                      priors$prior_raw_scale_mu_mean  <- c(if_null_then_set_to(priors$prior_raw_scale_mu_mean, array(dim = c(n_index_tests, 1), 0.0)))
                                      priors$prior_raw_scale_mu_SD    <- c(if_null_then_set_to(priors$prior_raw_scale_mu_SD,   array(dim = c(n_index_tests, 1), 0.5)))
                                      priors$prior_raw_scale_tau_SD   <- c(if_null_then_set_to(priors$prior_raw_scale_tau_SD,  array(dim = c(n_index_tests, 1), 0.5)))
                                      priors$prior_raw_scale_sigma_SD <- c(if_null_then_set_to(priors$prior_raw_scale_sigma_SD, rep(0.5, 1)))
                                }
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
                                priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, array(1.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                ##
                                n_total_pooled_cat  <- sum(n_cat)
                                n_thr_random = n_thr * n_studies
                                n_total_C_if_random <- sum(n_thr_random)
                                ##
                                priors$n_total_C_if_random <- n_total_C_if_random
                                priors$n_total_pooled_cat  <- n_total_pooled_cat
                                ##
                                priors$alpha_lb         <- if_null_then_set_to(priors$alpha_lb, 1.0)
                                priors$prior_alpha_mean <- if_null_then_set_to(priors$prior_alpha_mean, array(0.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                priors$prior_alpha_SD   <- if_null_then_set_to(priors$prior_alpha_SD,   array(10.0, dim = c(n_index_tests, max(n_thr) + 1)))
                                
                              
                                  #### add stuff for ordinal-HSROC-NMA here
                              
                              # put data for the following data variables:
                              #   ##
                              #   prior_beta_mu_mean, 
                              # prior_beta_mu_SD, 
                              # prior_beta_tau_SD, 
                              # prior_beta_sigma_SD,
                              # ###
                              # prior_raw_scale_mu_mean, 
                              # prior_raw_scale_mu_SD,
                              # prior_raw_scale_tau_SD, 
                              # prior_raw_scale_sigma_SD,
                              # ##
                              # prior_dirichlet_alpha
                      
                    }
      
      
      
      
                    
      
    }
    ##
    ## output list:
    ##
    return(priors)
    
 
  
 
  
}





 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

