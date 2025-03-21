 

 


#' R_fn_set_inits_NMA
#' @keywords internal
#' @export
R_fn_set_inits_NMA <- function(    inits, 
                                   priors,
                                   ##
                                   n_studies,
                                   n_index_tests,
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
  
          
          #       matrix[n_studies, 2] beta_eta_z;       //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
  #               array[n_index_tests] matrix[n_studies, 2] beta_delta_z;       //// Standard normal RVs for test-specific effects
  

          ##
          n_cat <- n_thr + 1
          n_tests <- n_index_tests
          ##
          if (cts == TRUE) { ## Jones-based model
                        ##
                        ## Set default inits for the locations ("beta"):
                        ##
                        inits$beta_mu <- if_null_then_set_to(inits$beta_mu, array(dim = c(n_tests, 2), c(rep(-1.01, n_tests), rep(+1.01, n_tests))))
                        ##
                        ## NMA params for beta:
                        ##
                        inits$beta_sigma   <- if_null_then_set_to(inits$beta_sigma, rep(0.001, 2))
                        inits$beta_tau     <- if_null_then_set_to(inits$beta_tau,      array(dim = c(n_tests, 2),   0.001))
                        inits$beta_eta_z   <- if_null_then_set_to(inits$beta_eta_z,    array(dim = c(n_studies, 2), 0.001))
                        ##
                        init_beta_delta_z <- list()    
                        for (t in 1:n_tests) { 
                          init_beta_delta_z[[t]] <- array(dim = c(n_studies, 2), 0.001)
                        }
                        inits$beta_delta_z <- if_null_then_set_to(inits$beta_delta_z, init_beta_delta_z)
                        ##
                        ## Set default inits for the raw scales ("gamma"):
                        ## NOTE: Default values depend of whether using log-normal / exp() for scales (as values will need to be of smaller 
                        ## magnitude as exp() amplifies much more that softplus(x)).
                        ##
                        if (softplus == TRUE) {
                          inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, array(dim = c(n_tests, 2), c(rep(0.55, n_tests), rep(0.55, n_tests))))
                        } else { 
                          inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, array(dim = c(n_tests, 2), c(rep(0.001, n_tests), rep(0.001, n_tests))))
                        }
                        ##
                        ## NMA params for raw scale ("gamma")::
                        ##
                        inits$raw_scale_sigma   <- if_null_then_set_to(inits$raw_scale_sigma, rep(0.001, 2))
                        inits$raw_scale_tau     <- if_null_then_set_to(inits$raw_scale_tau,   array(dim = c(n_tests, 2),   0.001))
                        inits$raw_scale_eta_z   <- if_null_then_set_to(inits$raw_scale_eta_z, array(dim = c(n_studies, 2), 0.001))
                        ##
                        init_raw_scale_delta_z <- list()    
                        for (t in 1:n_tests) { 
                          init_raw_scale_delta_z[[t]] <- array(dim = c(n_studies, 2), 0.001)
                        }
                        inits$raw_scale_delta_z <- if_null_then_set_to(inits$raw_scale_delta_z, init_raw_scale_delta_z)
                        ##
                        ## Set default inits for box-cox ("lambda"):
                        ##
                        inits$lambda <- if_null_then_set_to(inits$lambda, rep(0.001, n_tests))
                        
          } else if (cts == FALSE) { ## ordinal (Xu-based or R&G-based)
            
                        if (ord_model_parameterisation == "R&G") {
                                    ##
                                    ## Set default inits for the locations ("beta"):
                                    ##
                                    inits$beta_mu <- if_null_then_set_to(inits$beta_mu, 0.0)
                                    inits$beta_SD <- if_null_then_set_to(inits$beta_SD, 0.001)
                                    inits$beta_z  <- if_null_then_set_to(inits$beta_z,  rep(0.001, n_studies))
                                    ##
                                    ## Default inits for raw scale parameters ("gamma"):
                                    ##
                                    if (softplus == TRUE) {
                                      inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, 0.55)  ## since: softplus(+0.5*0.55) ~ 1,  softplus(-0.5*0.55) ~ 1.
                                    } else { 
                                      inits$raw_scale_mu <- if_null_then_set_to(inits$raw_scale_mu, 0.001) ## since: exp(+0.5*0.001) ~ 1,  exp(-0.5*0.001) ~ 1.
                                    }
                                    inits$raw_scale_SD <- if_null_then_set_to(inits$raw_scale_mu, 0.001)
                                    inits$raw_scale_z  <- if_null_then_set_to(inits$raw_scale_mu, rep(0.001, n_studies))
                          
                          
                        } else if (model_parameterisation == "Xu") {
                                    ##
                                    ## Set default inits for the locations ("beta"):
                                    ##
                                    inits$beta_mu <- if_null_then_set_to(inits$beta_mu, c(-1, +1))
                                    inits$beta_SD <- if_null_then_set_to(inits$beta_SD, rep(0.001, 2))
                                    inits$beta_z  <- if_null_then_set_to(inits$beta_z,  array(0.001, dim = c(2, n_studies)))
                                    ##
                                    ## Default inits for the between-study corr matrix ("beta_L_Omega"):
                                    ##
                                    inits_Omega <- array(0.001, dim = c(2, 2))
                                    diag(inits_Omega) <- c(1.0, 1.0)
                                    inits$beta_L_Omega <- t(chol(inits_Omega))
                          
                        }
                        ##
                        ## Default inits for the cutpoints and/or induced-Dirichlet parameters:
                        ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                        ##
                        if (random_thresholds == FALSE) { 
                                    ##
                                    ## Default inits for the cutpoints:
                                    ##
                                    inits$C <- seq(from = -2.0, to = 2.0, length = n_thr)
                          
                        } else if (random_thresholds == TRUE) { 
                                    ##
                                    ## Default inits for the cutpoints:
                                    ##
                                    cutpoint_vec <- seq(from = -2.0, to = 2.0, length = n_thr)
                                    inits$C_array <- list(cutpoint_vec, cutpoint_vec)
                                    ##
                                    ## Default inits for "induced-Dirichlet" between-study model (for the cutpoints):
                                    ##
                                   
                                    inits$dirichlet_cat_means_phi <- rep(1/n_cat, n_cat)
                                    inits$kappa <- 100
                          
                        }
                        
            
            
            
          }
 
          ##
          ## output list:
          ##
          return(inits)
  
  
  
}





 

  
  
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

