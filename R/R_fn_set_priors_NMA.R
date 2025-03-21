




 


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
            ## Set priors for raw scales ("gamma"):
            ##
            if (softplus == TRUE) {
                  priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                                                   array(dim = c(n_index_tests, 2), 0.0))
                  priors$prior_raw_scale_mu_SD <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                                                 array(dim = c(n_index_tests, 2), 1.0))
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
                                                                      array(dim = c(n_index_tests, 2), 1.0))
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
      
    } else if (cts == FALSE) { ## ordinal
      
    }
    ##
    ## output list:
    ##
    return(priors)
    
 
  
 
  
}





 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

