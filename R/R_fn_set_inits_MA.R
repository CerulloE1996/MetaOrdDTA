





#' R_fn_set_inits_MA
#' @keywords internal
#' @export
R_fn_set_inits_MA <- function(    inits, 
                                  ##
                                  priors,
                                  ##
                                  X = X,
                                  ##
                                  n_studies,
                                  n_tests = NULL, ##  may be dummy argument (needed for pkg)
                                  n_index_tests = NULL,   ##  may be dummy argument (needed for pkg)
                                  n_thr,                                 
                                  ##
                                  cts,
                                  ##
                                  model_parameterisation,
                                  random_thresholds,
                                  Dirichlet_random_effects_type,
                                  ##
                                  softplus
)  {
  
        
        cov_info_list <- R_fn_get_covariate_info_MA(X = X, 
                                                    model_parameterisation = model_parameterisation)
        ##
        n_covariates_max <- cov_info_list$n_covariates_max
        print(paste("n_covariates_max = ", cov_info_list$n_covariates_max))
        ##
        n_cat <- n_thr + 1
        ##
        { 
                if (cts == TRUE) { ## Jones-based model

                                ##
                                ## Set default inits for the locations ("beta"):
                                ##
                                inits$beta_mu <- if_null_then_set_to( inits$beta_mu,
                                                                      matrix(c(-1, +1), nrow = 2, ncol = n_covariates_max))
                                inits$beta_SD <- if_null_then_set_to( inits$beta_SD,
                                                                      rep(0.001, 2))
                                inits$beta_z  <- if_null_then_set_to( inits$beta_z,  
                                                                      array(0.001, dim = c(n_studies, 2)))
                                ##
                                ##
                                ## Set default inits for the raw scales ("gamma"):
                                ## NOTE: Default values depend of whether using log-normal / exp() 
                                ## for scales (as values will need to be of smaller 
                                ## magnitude as exp() amplifies much more that softplus(x)).
                                ##
                                ## since: softplus_scaled(0) = 1.0 and exp(0) = 1.0:
                                inits$raw_scale_mu <- if_null_then_set_to( inits$raw_scale_mu,
                                                                           matrix(0.001, nrow = 2, ncol = n_covariates_max))
                                inits$raw_scale_SD <- if_null_then_set_to( inits$raw_scale_SD, 
                                                                           rep(0.001, 2))
                                inits$raw_scale_z  <- if_null_then_set_to( inits$raw_scale_z, 
                                                                           array(0.001, dim = c(n_studies, 2)))
                                ##
                                ##
                                ## Set default inits for box-cox ("lambda"):
                                ##
                                inits$lambda <- if_null_then_set_to( inits$lambda, 
                                                                     rep(0.001, 2))
                                ##
                                ## Default inits for the between-study corr matrix ("beta_corr"):
                                ##
                                inits$beta_corr <- if_null_then_set_to( inits$beta_corr, 
                                                                        0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
                                ##
                                ## Default inits for the between-study corr matrix ("raw_scale_corr"):
                                ##
                                inits$raw_scale_corr <- if_null_then_set_to( inits$raw_scale_corr, 
                                                                             0.5*(priors$raw_scale_corr_lb + priors$raw_scale_corr_ub))
                  
                } else if (cts == FALSE) { ## ordinal (Xu-based or R&G-based)
                  
                             if (model_parameterisation %in% c("HSROC", "R&G", "Gatsonis")) {  

                                        ##
                                        ## Set default inits for the locations ("beta"):
                                        ##
                                        inits$beta_mu <- if_null_then_set_to(inits$beta_mu, 0.001)
                                        inits$beta_SD <- if_null_then_set_to(inits$beta_SD, 0.001)
                                        inits$beta_z  <- if_null_then_set_to(inits$beta_z,  rep(0.001, n_studies))
                                        ##
                                        ## Default inits for raw scale parameters ("gamma"):
                                        ##
                                        ## since: softplus_scaled(0) = 1.0 and exp(0) = 1.0:
                                        inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, 0.001)   
                                        inits$raw_scale_SD <- if_null_then_set_to(inits$raw_scale_SD, 0.001)
                                        inits$raw_scale_z  <- if_null_then_set_to(inits$raw_scale_z, rep(0.001, n_studies))
                                        ##
                                        ## Default inits for the cutpoints and/or induced-Dirichlet parameters:
                                        ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                                        ##
                                        if (random_thresholds == FALSE) { 
                                                ##
                                                ## Default inits for the cutpoints:
                                                ##
                                                inits$C_raw_vec <- if_null_then_set_to(inits$C, rep(-2.0, n_thr))
                                                
                                        } else if (random_thresholds == TRUE) { 
                                                ##
                                                ## Default inits for the cutpoints:
                                                ##
                                                mat <- matrix(-2.0, nrow = n_studies, ncol = n_thr)
                                                inits$C_raw <-  if_null_then_set_to(inits$C_raw, mat)
                                                ##
                                                target_sds <- rep(0.001, n_thr + 1)
                                                inits$target_sds <- if_null_then_set_to( inits$target_sds, target_sds)
                                                ##
                                                inits$dirichlet_phi <- if_null_then_set_to( inits$dirichlet_phi, rep(1/n_cat, n_cat))
                                                
                                        }
                                     
                             } else if (model_parameterisation %in% c("bivariate", "Xu")) { 

                                        ##
                                        ## Set default inits for the locations ("beta"):
                                        ##
                                        inits$beta_mu <- if_null_then_set_to( inits$beta_mu,
                                                                              matrix(c(-1, +1), nrow = 2, ncol = n_covariates_max))
                                        inits$beta_SD <- if_null_then_set_to(inits$beta_SD, rep(0.001, 2))
                                        inits$beta_z  <- if_null_then_set_to(inits$beta_z,  array(0.001, dim = c(n_studies, 2)))
                                        ##
                                        ## Default inits for the between-study corr matrix ("beta_corr"):
                                        ##
                                        inits$beta_corr <- if_null_then_set_to(inits$beta_corr, 
                                                                            0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
                                        ##
                                        ## Default inits for the cutpoints and/or induced-Dirichlet parameters:
                                        ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                                        ##
                                        if (random_thresholds == FALSE) { 
                                                ##
                                                ## Default inits for the cutpoints:
                                                ##
                                                C_raw_vec <- rep(-2.0, n_thr) ## seq(from = -2.0, to = 2.0, length = n_thr)
                                                inits$C_raw_vec <- if_null_then_set_to(inits$C_raw_vec,  list(C_raw_vec, C_raw_vec))
                                                # ##

                                        } else if (random_thresholds == TRUE) { 
                                                ##
                                                ## Default inits for the cutpoints:
                                                ##
                                                # cutpoint_vec <- seq(from = -2.0, to = 2.0, length = n_thr)
                                                # inits$C_array <- if_null_then_set_to(inits$C_array, list(cutpoint_vec, cutpoint_vec))
                                                ##
                                                mat <- matrix(-2.0, nrow = n_studies, ncol = n_thr)
                                                inits$C_raw <-  if_null_then_set_to(inits$C_raw, list(mat, mat))
                                                ##
                                                target_sds <- rep(0.001, n_thr + 1)
                                                inits$target_sds <- if_null_then_set_to(inits$target_sds, list(target_sds, target_sds))
                                                ##
                                                dirichlet_phi <- rep(1/n_cat, n_cat)
                                                inits$dirichlet_phi <- if_null_then_set_to( inits$dirichlet_phi, list(dirichlet_phi, dirichlet_phi))
                                                
                                        }
                                     
                             }

                    
                }
          
        }
        ##
        ## output list:
        ##
        return(inits)
        
}



  





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

