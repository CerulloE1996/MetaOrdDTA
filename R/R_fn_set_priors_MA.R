



 


#' R_fn_set_priors_MA
#' @keywords internal
#' @export
R_fn_set_priors_MA <- function(   priors,
                                  ##
                                  intercept_only,
                                  cov_data,
                                  ##
                                  n_studies = NULL,     ## keep as dummy input!
                                  n_index_tests = NULL, ## keep as dummy input!
                                  n_studies_per_test = NULL, ## keep as dummy input!
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
  
        n_index_tests <- 1 ## since MA (not NMA)
        ##
        cov_info_list <- R_fn_get_covariate_info_MA( intercept_only = intercept_only,
                                                     cov_data = cov_data, 
                                                     model_parameterisation = model_parameterisation,
                                                     n_studies = n_studies)
        ##
        n_covariates_max <- cov_info_list$n_covariates_max
        print(paste("n_covariates_max = ", cov_info_list$n_covariates_max))
        ##
        priors$use_probit_link <- if_null_then_set_to(priors$use_probit_link, 1)
        ##
        if (priors$use_probit_link == 1) { 
          mult <- 1.0
        } else { 
          mult <- 1.702
        }
        ##
        priors$K_fold_CV_indicator <- if_null_then_set_to(priors$K_fold_CV_indicator, rep(1, n_studies))
        ##
        if (cts == TRUE) { 
          
                  message("Configuring priors for continuous (Jones et al) -based model")
                  ##
                  ## Set default priors for the locations ("beta") and/or check that priors are in correct format:
                  ##
                  priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean,
                                                                   matrix(0.0, nrow = 2, ncol = n_covariates_max))
                  ##
                  ## probit scale so this is ~ unif. (but may rule-out "extreme" values (e.g.,0.999) depending on cutpoint value, 
                  ## so still need to think about the priors!!):
                  priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD,
                                                                   matrix(5.0*mult, nrow = 2, ncol = n_covariates_max))
                  ##
                  priors$prior_beta_SD_mean <- if_null_then_set_to(priors$prior_beta_SD_mean, rep(0.0, 2))
                  ##
                  ## this allows a lot of heterogeneity on the probit scale, e.g.: pnorm(0.0 +/- 2*0.5) gives 15.9%/84.1%.:
                  priors$prior_beta_SD_SD   <- if_null_then_set_to(priors$prior_beta_SD_SD,  
                                                                   rep(1.0*mult, 2)) 
                  ##
                  ## Set default priors for the raw scales ("gamma"):
                  ##
                  ## since: exp(0.0) = 1.0 and sp_scaled(0)  = 1.0.:
                  priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean,
                                                                   matrix(0.0, nrow = 2, ncol = n_covariates_max))
                  ##
                  priors$prior_raw_scale_SD_mean <- if_null_then_set_to( priors$prior_raw_scale_SD_mean,
                                                                         rep(0.0, 2))
                  ##
                  ## prior SD's depend on whether using scaled-softplus fn or exp():
                  if (softplus == TRUE) {
                          ##
                          ## since: e.g.: sp_scaled(0.0 +2*2.5) ~ 7.22, sp_scaled(0.0 - 2*2.5) ~ 0.0097:
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(  priors$prior_raw_scale_mu_SD,
                                                                                matrix(2.5*mult, nrow = 2, ncol = n_covariates_max))
                          ##
                          ## sp_scaled(0.0 + 2*1.0) ~ 3.07 (-> ~ 3x the mean (if the mean pooled scale is 1)):
                          priors$prior_raw_scale_SD_SD   <- if_null_then_set_to(priors$prior_raw_scale_SD_SD, 
                                                                                rep(1.0*mult, 2)) 
                          ## if using log-normal / exp() for scales (priors will need to be of smaller magnitude as exp() amplifies
                          ## much more that softplus(x)):
                  } else { 
                          ##
                          ## since: e.g., exp(0.0 + 2*1.0) = 7.39:
                          priors$prior_raw_scale_mu_SD <- if_null_then_set_to(  priors$prior_raw_scale_mu_SD,
                                                                                matrix(1.0*mult, nrow = 2, ncol = n_covariates_max))
                          ##
                          ## since: e.g., exp(0.0 + 2*0.50) ~  2.72 (-> ~ 3x the mean (if the mean pooled scale is 1)):
                          priors$prior_raw_scale_SD_SD   <- if_null_then_set_to(priors$prior_raw_scale_SD_SD, rep(1.0*mult, 2))  
                  }
                  ##
                  ## Set default priors for box-cox ("lamnda"):
                  ##
                  priors$prior_boxcox_lambda_mean <- if_null_then_set_to(priors$prior_boxcox_lambda_mean, 0.0)
                  ##
                  priors$prior_boxcox_lambda_SD   <- if_null_then_set_to(priors$prior_boxcox_lambda_SD,   1.0)
                  ##
                  ## Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                  ##
                  priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                  priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                  ##
                  priors$prior_beta_corr_LKJ      <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                  ##
                  ## Set default priors for the between-study corr matricex for raw_scale and "raw_scale_L_Omega"):
                  ##
                  priors$raw_scale_corr_lb <- if_null_then_set_to(priors$raw_scale_corr_lb, -1.0)
                  priors$raw_scale_corr_ub <- if_null_then_set_to(priors$raw_scale_corr_ub, +1.0)
                  ##
                  priors$prior_raw_scale_corr_LKJ <- if_null_then_set_to(priors$prior_raw_scale_corr_LKJ, 2.0)
            
        } else if (cts == FALSE) { ## ordinal
          
                  message("Configuring priors for ordinal model")
                  ##
                  if (model_parameterisation %in% c("HSROC", "R&G", "Gatsonis")) { 
                    
                          message("Configuring priors for HSROC-based (Rutter & Gatsonis) -based ordinal model")
                          ##
                          ## Set default priors for the locations ("beta"):
                          ##
                          priors$prior_beta_mu_mean <- if_null_then_set_to( priors$prior_beta_mu_mean, 
                                                                            rep(0.0, n_covariates_max))
                          ## probit scale so this is ~ unif. (but may rule-out "extreme" values (e.g.,0.999) depending on cutpoint value, 
                          ## so still need to think about the priors!!):
                          priors$prior_beta_mu_SD   <- if_null_then_set_to( priors$prior_beta_mu_SD, 
                                                                            rep(1.0*mult, n_covariates_max))
                          ##
                          priors$prior_beta_SD_mean <- if_null_then_set_to(priors$prior_beta_SD_mean, 0.0)
                          ##
                          ## this allows a lot of heterogeneity on the probit scale, e.g.: pnorm(0.0 +/- 2*0.5) gives 15.9%/84.1%.:
                          priors$prior_beta_SD_SD   <- if_null_then_set_to(priors$prior_beta_SD_SD,   
                                                                           1.0*mult)
                          ##
                          ## Set default priors for the raw scales ("gamma"):
                          ##
                          ## since: exp(0.0) = 1.0 and sp_scaled(0)  = 1.0:
                          priors$prior_raw_scale_mu_mean <- if_null_then_set_to(priors$prior_raw_scale_mu_mean, 
                                                                                0.0) 
                          priors$prior_raw_scale_SD_mean <- if_null_then_set_to(priors$prior_raw_scale_SD_mean, 
                                                                                0.0)
                          ## prior SD's depend on whether using scaled-softplus fn or exp():
                          if (softplus == TRUE) {
                            
                                      priors$prior_raw_scale_mu_SD   <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                                                                            1.0*mult)
                                      ## allows a lot of heterogeneity on the "probit-sp_scaled" scale: sp_scaled(0.0 - 2*0.75) ~ 0.291, sp_scaled(0.0 + 2*0.75) ~ 2.46:
                                      priors$prior_raw_scale_SD_SD   <- if_null_then_set_to(priors$prior_raw_scale_SD_SD, 
                                                                                            1.0*mult)
                                      
                          } else {    ## if using log-normal / exp() for scales (priors will need to be of smaller magnitude as exp() amplifies much more 
                                      ## that softplus(x)).
                                      priors$prior_raw_scale_mu_SD   <- if_null_then_set_to(priors$prior_raw_scale_mu_SD, 
                                                                                            1.0*mult)
                                      ## allows a lot of heterogeneity on the "probit-exp" scale: exp(0.0 - 2*0.5) ~0.37, exp(0.0 + 2*0.5) ~ 7.39:
                                      priors$prior_raw_scale_SD_SD   <- if_null_then_set_to(priors$prior_raw_scale_SD_SD, 
                                                                                            1.0*mult) 
                                      
                          }
                          ##
                          
                  } else if (model_parameterisation %in% c("bivariate", "Xu")) {
                    
                          message("Configuring priors for Xu et al. -based ordinal model")
                          ##
                          ## Set default priors for the locations ("beta"):
                          ##
                          priors$prior_beta_mu_mean <- if_null_then_set_to(priors$prior_beta_mu_mean,
                                                                           matrix(0.0, nrow = 2, ncol = n_covariates_max))
                          ##
                          priors$prior_beta_mu_SD <- if_null_then_set_to(priors$prior_beta_mu_SD,
                                                                           matrix(1.0*mult, nrow = 2, ncol = n_covariates_max))
                          ## probit scale so this is ~ unif. (but may rule-out "extreme" values (e.g.,0.999) depending on cutpoint value,
                          ## so still need to think about the priors!!):
                          priors$prior_beta_SD_mean <- if_null_then_set_to(priors$prior_beta_SD_mean, 
                                                                           c(0.0, 0.0))
                          ## this allows a lot of heterogeneity on the probit scale, e.g.: pnorm(0.0 +/- 2*0.5) gives 15.9%/84.1%:
                          priors$prior_beta_SD_SD   <- if_null_then_set_to(priors$prior_beta_SD_SD,   
                                                                           rep(1.0*mult, 2))
                          ##
                          ## Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
                          ##
                          priors$beta_corr_lb <- if_null_then_set_to(priors$beta_corr_lb, -1.0)
                          priors$beta_corr_ub <- if_null_then_set_to(priors$beta_corr_ub, +1.0)
                          priors$prior_beta_corr_LKJ      <- if_null_then_set_to(priors$prior_beta_corr_LKJ, 2.0)
                          
                  }
                  ##
                  ## Default priors for the cutpoints and/or induced-Dirichlet parameters:
                  ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                  ##
                  if (random_thresholds == FALSE) { 
                    
                          message("Configuring fixed-thresholds-based ordinal probability (using Induced-Dirichlet as a prior) priors")
                          ##
                          ## Set default priors for the cutpoints (using Dirichlet ** priors ** (not model) on the "Induced-Dirichlet" ordinal probs):
                          ##
                          if (model_parameterisation %in% c("HSROC", "R&G", "Gatsonis")) {
                            priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, 
                                                                                rep(1.0, n_cat))
                          } else { 
                            priors$prior_dirichlet_alpha <- if_null_then_set_to(priors$prior_dirichlet_alpha, 
                                                                                rep(list(rep(1.0, n_cat)), 2))
                          }
                          
                  } else if (random_thresholds == TRUE) {
                    
                              message("Configuring random-thresholds-based ordinal probability (using Induced-Dirichlet between-study model) priors")
                              ##
                              ## Set default priors for the cutpoints  (using Dirichlet ** between-study model ** (not priors) on the "Induced-Dirichlet" ordinal probs):
                              ##
                          if (model_parameterisation %in% c("HSROC", "R&G", "Gatsonis")) {
                            
                                priors$prior_dirichlet_alpha <- if_null_then_set_to( priors$prior_dirichlet_alpha, 
                                                                                     rep(1.0, n_thr + 1))
                                ##
                                priors$prior_kappa_mean   <-  if_null_then_set_to( priors$prior_kappa_mean, 
                                                                                           log(100))
                                ##
                                priors$prior_kappa_SD   <-  if_null_then_set_to( priors$prior_kappa_SD, 
                                                                                 1.0)
                                ##
                                priors$kappa_lb <- if_null_then_set_to( priors$kappa_lb, 
                                                                        0.001)
                                
                          } else { 
                            
                                prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
                                priors$prior_dirichlet_alpha <- if_null_then_set_to( priors$prior_dirichlet_alpha, 
                                                                                   list(prior_dirichlet_alpha, prior_dirichlet_alpha))
                                ##
                                priors$prior_kappa_mean   <- if_null_then_set_to(  priors$prior_kappa_mean, 
                                                                                           c(log(200), log(50)))
                                ##
                                priors$prior_kappa_SD   <- if_null_then_set_to(  priors$prior_kappa_SD, 
                                                                                   c(1.0, 1.0))
                                ##
                                priors$kappa_lb <- if_null_then_set_to( priors$kappa_lb, 
                                                                        0.001)
                                ##
                                
                                
                          }
                  }
          
        }
        ##
        ## output:
        ##
        return(priors)
        
}








 
   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

