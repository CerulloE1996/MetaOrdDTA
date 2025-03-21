





#' R_fn_set_inits_MA
#' @keywords internal
#' @export
R_fn_set_inits_MA <- function(    inits, 
                                  priors,
                                  ##
                                  n_studies,
                                  n_tests,
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
  
        n_cat <- n_thr + 1
        ##
        { 
                
                if (cts == TRUE) { ## Jones-based model
                              ##
                              ## Set default inits for the locations ("beta"):
                              ##
                              inits$beta_mu <- if_null_then_set_to(inits$beta_mu, c(-1, +1))
                              inits$beta_SD <- if_null_then_set_to(inits$beta_SD, rep(0.001, 2))
                              inits$beta_z  <- if_null_then_set_to(inits$beta_z,  array(0.001, dim = c(2, n_studies)))
                              ##
                              check_vec_length(inits, "beta_mu", 2)
                              check_vec_length(inits, "beta_SD", 2)
                              check_array_dims(inits, "beta_z", c(2, n_studies))
                              ##
                              ## Set default inits for the raw scales ("gamma"):
                              ## NOTE: Default values depend of whether using log-normal / exp() 
                              ## for scales (as values will need to be of smaller 
                              ## magnitude as exp() amplifies much more that softplus(x)).
                              ##
                              ## since: softplus_scaled(0) = 1.0 and exp(0) = 1.0:
                              inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, c(0.001, 0.001))  
                              inits$raw_scale_SD <- if_null_then_set_to(inits$raw_scale_SD, rep(0.001, 2))
                              inits$raw_scale_z  <- if_null_then_set_to(inits$raw_scale_z, 
                                                                        array(0.001, dim = c(2, n_studies)))
                              ##
                              check_vec_length(inits, "beta_mu", 2)
                              check_vec_length(inits, "beta_SD", 2)
                              check_array_dims(inits, "beta_z", c(2, n_studies))
                              ##
                              ## Set default inits for box-cox ("lambda"):
                              ##
                              inits$lambda <- if_null_then_set_to(inits$lambda, 0.001)
                              check_vec_length(inits, "lambda", 1)
                              ##
                              ## Default inits for the between-study corr matrix ("beta_corr"):
                              ##
                              inits$beta_corr <- if_null_then_set_to(inits$beta_corr, 
                                                                     0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
                              check_vec_length(inits, "beta_corr", 1)
                              ##
                              ## Default inits for the between-study corr matrix ("raw_scale_corr"):
                              ##
                              inits$raw_scale_corr <- if_null_then_set_to(inits$raw_scale_corr, 
                                                                          0.5*(priors$raw_scale_corr_lb + priors$raw_scale_corr_ub))
                              check_vec_length(inits, "raw_scale_corr", 1)
                  
                } else if (cts == FALSE) { ## ordinal (Xu-based or R&G-based)
                  
                             if (model_parameterisation == "R&G") { 
                                     ##
                                     ## Set default inits for the locations ("beta"):
                                     ##
                                     inits$beta_mu <- if_null_then_set_to(inits$beta_mu, 0.001)
                                     inits$beta_SD <- if_null_then_set_to(inits$beta_SD, 0.001)
                                     inits$beta_z  <- if_null_then_set_to(inits$beta_z,  rep(0.001, n_studies))
                                     ##
                                     check_vec_length(inits, "beta_mu", 1)
                                     check_vec_length(inits, "beta_SD", 1)
                                     check_vec_length(inits, "beta_z", n_studies)
                                     ##
                                     ## Default inits for raw scale parameters ("gamma"):
                                     ##
                                     ## since: softplus_scaled(0) = 1.0 and exp(0) = 1.0:
                                     inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, 0.001)   
                                     inits$raw_scale_SD <- if_null_then_set_to(inits$raw_scale_SD, 0.001)
                                     inits$raw_scale_z  <- if_null_then_set_to(inits$raw_scale_z, rep(0.001, n_studies))
                                     ##
                                     check_vec_length(inits, "raw_scale_mu", 1)
                                     check_vec_length(inits, "raw_scale_SD", 1)
                                     check_vec_length(inits, "raw_scale_z", n_studies)
                                     
                             } else if (model_parameterisation == "Xu") { 
                                     ##
                                     ## Set default inits for the locations ("beta"):
                                     ##
                                     inits$beta_mu <- if_null_then_set_to(inits$beta_mu, c(-1, +1))
                                     inits$beta_SD <- if_null_then_set_to(inits$beta_SD, rep(0.001, 2))
                                     inits$beta_z  <- if_null_then_set_to(inits$beta_z,  array(0.001, dim = c(2, n_studies)))
                                     ##
                                     check_vec_length(inits, "beta_mu", 2)
                                     check_vec_length(inits, "beta_SD", 2)
                                     check_array_dims(inits, "beta_z", c(2, n_studies))
                                     ##
                                     ## Default inits for the between-study corr matrix ("beta_corr"):
                                     ##
                                     inits$beta_corr <- if_null_then_set_to(inits$beta_corr, 
                                                                            0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
                                     ##
                                     check_vec_length(inits, "beta_corr", 1)
                                     
                             }
                             ##
                             ## Default inits for the cutpoints and/or induced-Dirichlet parameters:
                             ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                             ##
                             if (random_thresholds == FALSE) { 
                                     ##
                                     ## Default inits for the cutpoints:
                                     ##
                                     inits$C <- if_null_then_set_to(inits$C, seq(from = -2.0, to = 2.0, length = n_thr))
                                     # ##
                                     # check_vec_length(inits, "C", n_thr)
                                
                             } else if (random_thresholds == TRUE) { 
                                     ##
                                     ## Default inits for the cutpoints:
                                     ##
                                     cutpoint_vec <- seq(from = -2.0, to = 2.0, length = n_thr)
                                     inits$C_array <- if_null_then_set_to(inits$C_array, list(cutpoint_vec, cutpoint_vec))
                                     ##
                                     # check_list_of_vectors_length(inits, "C_array", n_thr)
                                     ##
                                     ## Default inits for "induced-Dirichlet" between-study model (for the cutpoints):
                                     # ##
                                     # inits$dirichlet_cat_means_phi <- if_null_then_set_to(inits$dirichlet_cat_means_phi, rep(1/n_cat, n_cat))
                                     # ##
                                     # inits$alpha <- if_null_then_set_to(inits$alpha, rep(10, n_cat))
                                     # ##
                                     # inits$kappa <- if_null_then_set_to(inits$kappa, 100)
                                     # # ##
                                     # # check_vec_length(inits, "dirichlet_cat_means_phi", n_cat)
                                     # # check_vec_length(inits, "kappa", 1)
                                
                             }
                    
                }
          
        }
        ##
        ## output list:
        ##
        return(inits)
        
}



  





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

