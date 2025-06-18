 

 


#' R_fn_set_inits_NMA
#' @keywords internal
#' @export
R_fn_set_inits_NMA <- function(    inits, 
                                   ##
                                   priors,
                                   ##
                                   X,
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
  
      cov_info_list <- R_fn_get_covariate_info_NMA(X = X, 
                                                   model_parameterisation = model_parameterisation, 
                                                   n_index_tests = n_index_tests,
                                                   n_studies = n_studies)
      ##
      n_covariates_max <- cov_info_list$n_covariates_max
      print(paste("n_covariates_max = ", cov_info_list$n_covariates_max))
      ##
      n_thr_max <- max(n_thr)
      n_cat <- n_thr + 1
      n_cat_max <- max(n_cat)
      ##
      n_tests <- n_index_tests
      ##
      if (cts == TRUE) { ## Jones-based model
        
            ##
            ## Set default inits for the locations ("beta"):
            ##
            ## array[n_index_tests] matrix[2, n_covariates_max] beta_mu; 
            mat <- matrix(NA, nrow = 2, ncol = n_covariates_max)
            mat[1, ] <- -1.01 ;   mat[2, ] <- +1.01
            inits$beta_mu <- if_null_then_set_to(inits$beta_mu, 
                                                 rep(list(mat), n_index_tests))
            ##
            ## NMA params for beta:
            ##
            inits$beta_sigma <- if_null_then_set_to(inits$beta_sigma, rep(0.001, 2))
            inits$beta_tau   <- if_null_then_set_to(inits$beta_tau,   array(dim = c(n_index_tests, 2),   0.001))
            inits$beta_eta_z <- if_null_then_set_to(inits$beta_eta_z, array(dim = c(n_studies, 2), 0.001))
            ##
            init_beta_delta_z <- list()    
            for (t in 1:n_index_tests) { 
              init_beta_delta_z[[t]] <- array(dim = c(n_studies, 2), 0.001)
            }
            inits$beta_delta_z <- if_null_then_set_to(inits$beta_delta_z, init_beta_delta_z)
            ##
            ## Set default inits for the raw scales ("gamma"):
            ## NOTE: Default values depend of whether using log-normal / exp() for scales (as values will need to be of smaller 
            ## magnitude as exp() amplifies much more that softplus(x)).
            ##
            ##  array[n_index_tests] matrix[2, n_covariates_max] raw_scale_mu; 
            ##
            mat <- matrix(NA, nrow = 2, ncol = n_covariates_max)
            mat[,] <- 0.001
            inits$raw_scale_mu <- if_null_then_set_to( inits$raw_scale_mu, 
                                                       rep(list(mat), n_index_tests))
            ##
            ## NMA params for raw scale ("gamma")::
            ##
            inits$raw_scale_sigma   <- if_null_then_set_to(inits$raw_scale_sigma, rep(0.001, 2))
            inits$raw_scale_tau     <- if_null_then_set_to(inits$raw_scale_tau,   array(dim = c(n_index_tests, 2),   0.001))
            inits$raw_scale_eta_z   <- if_null_then_set_to(inits$raw_scale_eta_z, array(dim = c(n_studies, 2), 0.001))
            inits$raw_scale_delta_z   <- if_null_then_set_to(inits$raw_scale_delta_z, array(dim = c(n_studies, 2), 0.001))
            ##
            ## Set default inits for box-cox ("lambda"):
            ##
            lambda <- rep(0.001, n_index_tests)
            inits$lambda <- if_null_then_set_to(inits$lambda, list(lambda, lambda))
            ##
            ## Default inits for the between-study corr matrix ("beta_corr" and "raw_scale_corr):
            ##
            print(paste("beta_corr_lb = ", priors$beta_corr_lb))
            print(paste("beta_corr_ub = ", priors$beta_corr_ub))
            inits$beta_corr      <- if_null_then_set_to(inits$beta_corr, 0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
            ##
            print(paste("raw_scale_corr_lb = ", priors$raw_scale_corr_lb))
            print(paste("raw_scale_corr_ub = ", priors$raw_scale_corr_ub))
            inits$raw_scale_corr <- if_null_then_set_to(inits$raw_scale_corr, 0.5*(priors$raw_scale_corr_lb + priors$raw_scale_corr_ub))
            
      } else if (cts == FALSE) { ## ordinal (Xu-based or R&G-based)
            
              if (model_parameterisation %in% c("R&G", "HSROC", "Gatsonis")) { 
                        
                    ##
                    ## Set default inits for the locations ("beta"):
                    ##
                    vec <- rep(0.0, n_covariates_max)
                    inits$beta_mu <- if_null_then_set_to( inits$beta_mu, 
                                                          rep(list(vec), n_index_tests))
                    ##
                    vec <- rep(0.0, n_covariates_max)
                    inits$raw_scale_mu <- if_null_then_set_to( inits$raw_scale_mu, 
                                                               rep(list(vec), n_index_tests))
                    ##
                    inits$beta_sigma <- if_null_then_set_to( inits$beta_sigma, 
                                                             0.001)
                    inits$beta_tau   <- if_null_then_set_to( inits$beta_tau, 
                                                             rep(0.001, n_index_tests))
                    ##
                    inits$beta_eta_z <- if_null_then_set_to( inits$beta_eta_z,
                                                             rep(0.001, n_studies))
                    ##
                    init_beta_delta_z <- list()    
                    for (t in 1:n_index_tests) { 
                      init_beta_delta_z[[t]] <- rep(0.001, n_studies)
                    }
                    inits$beta_delta_z <- if_null_then_set_to(inits$beta_delta_z, 
                                                              init_beta_delta_z)
                    ##
                    ##
                    ##
                    inits$raw_scale_sigma   <- if_null_then_set_to(inits$raw_scale_sigma, 
                                                                   0.001)
                    ##
                    inits$raw_scale_tau     <- if_null_then_set_to(inits$raw_scale_tau,   
                                                                   rep(0.001, n_index_tests))
                    ##
                    inits$raw_scale_eta_z   <- if_null_then_set_to(inits$raw_scale_eta_z, 
                                                                   rep(0.001, n_studies))
                    ##
                    inits$raw_scale_delta_z   <- if_null_then_set_to(inits$raw_scale_delta_z, 
                                                                          array(0.001, dim = c(n_index_tests, n_studies)))
                    ##
                    if (random_thresholds == FALSE) { 
                      
                            ## array[2] vector<lower=-15.0, upper=15.0>[n_total_C_per_class_if_fixed] C_raw_vec;
                            n_total_C_if_fixed <- sum(n_thr)
                            vec <- rep(-3.0, n_total_C_if_fixed)
                            inits$C_raw_vec   <- if_null_then_set_to( inits$C_raw_vec, 
                                                                      vec)

                    } else if (random_thresholds == TRUE) { 
                            
                            ## array[2] vector[n_total_C_if_random] C_raw_vec;   
                            ##
                            n_total_C_if_random <- sum(n_thr * n_studies)
                            vec <- rep(-3.0, n_total_C_if_random)
                            inits$C_raw_vec   <- if_null_then_set_to( inits$C_raw_vec, 
                                                                      vec)
                            ##
                            ## array[2] simplex[n_index_tests, n_cat_max] raw_dirichlet_phi;
                            ##
                            raw_dirichlet_phi <- array(NA, dim = c(n_index_tests, n_cat_max - 1))
                            for (t in 1:n_index_tests) {
                                raw_dirichlet_phi[t, 1:n_thr[t]] <- MetaOrdDTA:::generate_inits_for_raw_simplex_vec(n_thr = n_thr[t],
                                                                                                                    seed = 123)
                            }
                            raw_dirichlet_phi_vec <- c(t(raw_dirichlet_phi))
                            raw_dirichlet_phi_vec <- raw_dirichlet_phi_vec[!is.na(raw_dirichlet_phi_vec)]
                            inits$raw_dirichlet_phi_vec   <- if_null_then_set_to( inits$raw_dirichlet_phi_vec, 
                                                                                  raw_dirichlet_phi_vec)
                            ##
                            kappa <- rep(50.0, n_index_tests)
                            inits$kappa   <- if_null_then_set_to( inits$kappa, 
                                                                  kappa)
                            
                            
                    }
                          
                  } else if (model_parameterisation %in% c("Xu", "bivariate")) {
                    
                        ##
                        ## Set default inits for the locations ("beta"):
                        ##
                        mat <- matrix(NA, nrow = 2, ncol = n_covariates_max)
                        mat[1, ] <- -1.01 ;   mat[2, ] <- +1.01
                        inits$beta_mu <- if_null_then_set_to(inits$beta_mu, 
                                                             rep(list(mat), n_index_tests))
                        # ##
                        # ## Default inits for the between-study corr matrix ("beta_L_Omega"):
                        # ##
                        # inits_Omega <- array(0.001, dim = c(2, 2))
                        # diag(inits_Omega) <- c(1.0, 1.0)
                        # inits$beta_L_Omega <- t(chol(inits_Omega))
                        ##
                        inits$beta_sigma <- if_null_then_set_to(inits$beta_sigma, rep(0.001, 2))
                        inits$beta_tau   <- if_null_then_set_to(inits$beta_tau, array(0.001, dim = c(n_index_tests, 2)))
                        ##
                        inits$beta_eta_z <- if_null_then_set_to(inits$beta_eta_z, array(0.001, dim = c(n_studies, 2)))
                        ##
                        init_beta_delta_z <- list()
                        for (t in 1:n_index_tests) {
                          init_beta_delta_z[[t]] <- array(0.001, dim = c(n_studies, 2))
                        }
                        inits$beta_delta_z <- if_null_then_set_to(inits$beta_delta_z, init_beta_delta_z)
                        ##
                        ## Default inits for the between-study corr matrix ("beta_corr"):
                        ##
                        inits$beta_corr <- if_null_then_set_to(inits$beta_corr, 0.5*(priors$beta_corr_lb + priors$beta_corr_ub))
                        # check_vec_length(inits, "beta_corr", 1)
                        # ##
                        
                        ##
                        ## Default inits for the cutpoints and/or induced-Dirichlet parameters:
                        ## NOTE: these are the same whether using "R&G" or "Xu" parameterisation
                        ##
                        if (random_thresholds == FALSE) { 

                              ## array[2] vector<lower=-15.0, upper=15.0>[n_total_C_per_class_if_fixed] C_raw_vec;
                              n_total_C_if_fixed <- sum(n_thr)
                              vec <- rep(-3.0, n_total_C_if_fixed)
                              inits$C_raw_vec   <- if_null_then_set_to( inits$C_raw_vec, 
                                                                        list(vec, vec))
                              
                        } else if (random_thresholds == TRUE) { 
                              ##)
                              n_total_C_if_random <- sum(n_thr * n_studies)
                              ##
                              ## array[2] vector[n_total_C_if_random] C_raw_vec;   
                              vec <- rep(-3.0, n_total_C_if_random)
                              inits$C_raw_vec   <- if_null_then_set_to( inits$C_raw_vec, 
                                                                        list(vec, vec))
                              ##
                              ## array[2] simplex[n_index_tests, n_cat_max] raw_dirichlet_phi;
                              ##
                              raw_dirichlet_phi <- array(NA, dim = c(n_index_tests, n_cat_max - 1))
                              for (t in 1:n_index_tests) {
                                raw_dirichlet_phi[t, 1:n_thr[t]] <- MetaOrdDTA:::generate_inits_for_raw_simplex_vec(n_thr = n_thr[t],
                                                                                                                    seed = 123)
                              }
                              raw_dirichlet_phi_vec <- c(t(raw_dirichlet_phi))
                              raw_dirichlet_phi_vec <- raw_dirichlet_phi_vec[!is.na(raw_dirichlet_phi_vec)]
                              inits$raw_dirichlet_phi_vec   <- if_null_then_set_to( inits$raw_dirichlet_phi_vec, 
                                                                                    list(raw_dirichlet_phi_vec, raw_dirichlet_phi_vec))
                              ##
                              kappa <- rep(50.0, n_index_tests)
                              inits$kappa   <- if_null_then_set_to( inits$kappa, 
                                                                    list(kappa, kappa))
                              
                        }
            }
            
          }
          ##
          ## output list:
          ##
          return(inits)
  
}





 

  
  
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

