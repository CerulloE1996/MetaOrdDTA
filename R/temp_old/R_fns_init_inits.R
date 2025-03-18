


#' init_inits
#' @keywords internal
#' @export
init_inits    <- function(init_model_outs,
                          init_lists_per_chain,
                          compile,
                          force_recompile,
                          cmdstanr_model_fit_obj,
                          n_chains_burnin,
                          n_params_main,
                          n_nuisance, 
                          N,
                          sample_nuisance,
                          Stan_model_file_path,
                          Stan_data_list,
                          Stan_cpp_user_header,
                          Stan_cpp_flags,
                          ...) {
 
  
  if (is.null(init_model_outs)) {
    warning("init_model_outs not specified - please create init_model_outs using init_model() and then pass as an argument to init_burnin()")
  }
  
  Model_type <- init_model_outs$Model_type
 
  n_params <- n_params_main + n_nuisance
  
  y <- init_model_outs$y
  
  if (is.null(init_lists_per_chain)) { 
    warning("initial values per chain (init_lists_per_chain) not supplied - using defaults")
  }
  
  
  n_covariates_per_outcome_mat <- init_model_outs$model_args_list$model_args_list$n_covariates_per_outcome_mat
   
  
  if (Model_type == "latent_trait") { 
    n_tests <- ncol(y)
    n_params_main <- 1 + sum(n_covariates_per_outcome_mat) + 2 * n_tests 
    n_nuisance <- n_tests * N
  } else if (Model_type == "LC_MVP") { 
    n_tests <- ncol(y)
    n_params_main <- 1 + sum(n_covariates_per_outcome_mat) + 2 * choose(n_tests, 2)
    n_nuisance <- n_tests * N
  } else if (Model_type == "MVP") { 
    n_tests <- ncol(y)
    n_params_main <-   sum(n_covariates_per_outcome_mat) + choose(n_tests, 2)
    n_nuisance <- n_tests * N
  } else { 
    
    if (is.null(n_params_main)) {
      warning("n_params_main not specified - will compute from Stan model")
    }

    if (is.null(n_nuisance)) {
      stop("n_nuisance not specified - please specify")
    }
    
  }
  
  
  if (Model_type != "Stan") { 
    if (Model_type == "MVP") {
      n_class <- 1
    } else {
      n_class <- 2
    }
  }

  mod <- NULL
 
 ####  n_chains_burnin <- 8 
  
  ##  ----------------------------   starting values - set defaults if user does not supply --- this is only for NON-Stan models - for Stan models inits are specified the same as they are for Stan
  
  mod <- NULL
  if (Model_type == "Stan") { 
    
                   dummy_json_file_path <- NULL
                   dummy_model_so_file <- NULL
                   
                    # //////   ----- if Stan model, then all 3 of:  Stan_model_file_path, Stan_data_list, and init_lists_per_chain need to be user-supplied
                    
                    bs_model <- init_model_outs$bs_model
                    print(paste("bs_model = "))
                    print(bs_model)
                    ##
                    json_file_path <- init_model_outs$json_file_path
                    print(paste("json_file_path = "))
                    print(json_file_path)
 
                    if (sample_nuisance == TRUE) { 
                      n_params_main <- n_params - n_nuisance
                    } else { 
                      n_params_main <- n_params
                    }
                    
                    ##
                    bs_names  <-  (bs_model$param_names())
                    ## if using bs names for param names:
                    param_names <- bs_names
                    ## Inits:
                    inits_unconstrained_vec_per_chain <- list()
                    for (kk in 1:n_chains_burnin) {
                      ## NOTE: param_unconstrain_json() Returns a vector of unconstrained parameters * given the constrained parameters *
                      json_string_for_inits_chain_kk <- BayesMVP:::convert_Stan_data_list_to_JSON(init_lists_per_chain[[kk]])
                      validated_json_string <- paste(readLines(json_string_for_inits_chain_kk), collapse="")
                      ##
                      inits_unconstrained_vec_per_chain[[kk]] <- bs_model$param_unconstrain_json(validated_json_string) ## error here
                      ####  inits_unconstrained_vec_per_chain[[kk]] <-  convert_JSON_string_to_R_vector(json_string_for_inits_chain_kk)
                    }
                    
 
  } else {
    
            if (Model_type == "LC_MVP") {
          
                  # prior_for_corr_a <- init_model_outs$prior_for_corr_a
                  # prior_for_corr_b <- init_model_outs$prior_for_corr_b
                  
                  X <- init_model_outs$X
                  n_covariates_per_outcome_mat <- init_model_outs$model_args_list$model_args_list$n_covariates_per_outcome_mat
                  n_covariates_max_nd <- max(n_covariates_per_outcome_mat[1, ])
                  n_covariates_max_d <-  max(n_covariates_per_outcome_mat[2, ])  
                  n_covariates_max <- max(n_covariates_max_nd, n_covariates_max_d)
                  X_nd <- X[[1]]
                  X_d <-  X[[2]]
                  
                  print(str(X_d))
              
                  prior_coeffs_mean_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_mean_mat
                  prior_coeffs_sd_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_sd_mat
          
                  ub_corr <- init_model_outs$model_args_list$model_args_list$ub_corr
                  lb_corr <- init_model_outs$model_args_list$model_args_list$lb_corr
                  lkj_cholesky_eta <- init_model_outs$model_args_list$model_args_list$lkj_cholesky_eta
                  
                  prev_prior_a <- init_model_outs$model_args_list$model_args_list$prev_prior_a
                  prev_prior_b <- init_model_outs$model_args_list$model_args_list$prev_prior_b
                  
                  known_num <- 0
                  known_values_indicator_list <- init_model_outs$model_args_list$model_args_list$known_values_indicator_list
                  known_values_list <- init_model_outs$model_args_list$model_args_list$known_values_list
                  
                  ## these vars need conversion to different types compatible w/ Stan
                  corr_force_positive <- as.integer(init_model_outs$model_args_list$model_args_list$corr_force_positive)
                  prior_only <-  as.integer(init_model_outs$model_args_list$model_args_list$prior_only)
                  
                  overflow_threshold <-  init_model_outs$model_args_list$model_args_list$overflow_threshold
                  underflow_threshold <- init_model_outs$model_args_list$model_args_list$underflow_threshold
                  
                  for (c in 1:n_class) {
                    for (t in 1:n_tests) {
                      X[[c]][[t]][1:N, 1:n_covariates_max] <- 1
                    }
                  }
                  
                  n_pops <- 1 
                  
                  print(prior_coeffs_mean_mat)
                  print(prior_coeffs_sd_mat)
                  
                  Phi_type <- init_model_outs$model_args_list$model_args_list$Phi_type
                  print(paste("Phi_type = ", Phi_type))
                  ##
                  Phi_type_int <- ifelse(Phi_type == "Phi", 1, 2)
                  
                  ## Make sure vecs of length 1 are R arrays (1d):
                  if (!(is.matrix(prev_prior_a))) {
                    prev_prior_a <-  matrix(prev_prior_a, ncol = 1)
                    prev_prior_b <-  matrix(prev_prior_b, ncol = 1)
                  }
                  
                  if (!(is.matrix(lkj_cholesky_eta))) {
                    lkj_cholesky_eta <-  matrix(lkj_cholesky_eta, ncol = 1)
                  }
                  
                  Stan_data_list <- list(N = N,
                                         n_tests = n_tests,
                                         y = y,
                                         n_class = n_class,
                                         n_pops =  n_pops,  ## multi-pop not supported yet (currently only in Stan version)
                                         pop =  c(rep(1, N)),
                                         ##
                                         n_covariates_max_nd = n_covariates_max_nd,
                                         n_covariates_max_d = n_covariates_max_d,
                                         n_covariates_max = n_covariates_max,
                                         X_nd = X_nd,
                                         X_d = X_d,
                                         n_covs_per_outcome = n_covariates_per_outcome_mat,
                                         ##
                                         corr_force_positive = corr_force_positive,
                                         known_num = known_num,
                                         ##
                                         # lb_corr = lb_corr,
                                         # ub_corr = ub_corr,
                                         overflow_threshold =  overflow_threshold,
                                         underflow_threshold = underflow_threshold,
                                         #### priors
                                         prior_only = prior_only,
                                         prior_beta_mean = prior_coeffs_mean_mat,
                                         prior_beta_sd = prior_coeffs_sd_mat,
                                         prior_LKJ = lkj_cholesky_eta,
                                         prior_p_alpha = prev_prior_a,
                                         prior_p_beta = prev_prior_b,
                                         ##
                                         Phi_type =  Phi_type_int,
                                         handle_numerical_issues = 1,
                                         fully_vectorised = 1)
                  
                  print(paste("Stan_data_list = "))
                  print(str(Stan_data_list))
                  
                  print(paste("prev_prior_a = "))
                  print(str(prev_prior_a))
                  print(paste("prev_prior_b = "))
                  print(str(prev_prior_b))
                  
                  print(paste("Stan_data_list$prior_p_alpha = "))
                  print(str(Stan_data_list$prior_p_alpha))
                  print(paste("Stan_data_list$prior_p_beta = "))
                  print(str(Stan_data_list$prior_p_beta))
                  
                  ## Compile using BridgeStan:
                  outs_init_bs_model <- BayesMVP:::init_bs_model_internal(   Stan_data_list = Stan_data_list,
                                                                             Stan_model_name = "LC_MVP_bin_PartialLog_v5.stan")
                  ##   LC_MVP_bin_PartialLog_v5.stan
                  ##   PO_LC_MVP_bin.stan
                  ##   LC_MVP_bin_w_mnl_cpp_grad_v1.stan
                  bs_model <- outs_init_bs_model$bs_model
                  json_file_path <- outs_init_bs_model$json_file_path         
                  Stan_model_file_path <- outs_init_bs_model$Stan_model_file_path  
                  
                  dummy_json_file_path <- json_file_path
                  dummy_model_so_file <- BayesMVP:::transform_stan_path(Stan_model_file_path)
                  
        } else if (Model_type == "MVP") {
                  
                  X <- init_model_outs$X
                  n_covariates_per_outcome_mat <- init_model_outs$model_args_list$model_args_list$n_covariates_per_outcome_mat
                  n_covariates_max_nd <- 999999 #max(n_covariates_per_outcome_mat[[1]])
                  n_covariates_max_d <-  999999 #max(n_covariates_per_outcome_mat[[2]])  
                  n_covariates_max <- max(unlist(n_covariates_per_outcome_mat))
                  
                  prior_coeffs_mean_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_mean_mat
                  prior_coeffs_sd_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_sd_mat
                  
                  ub_corr <- init_model_outs$model_args_list$model_args_list$ub_corr
                  lb_corr <- init_model_outs$model_args_list$model_args_list$lb_corr
                  lkj_cholesky_eta <- init_model_outs$model_args_list$model_args_list$lkj_cholesky_eta
                  
                  prev_prior_a <- init_model_outs$model_args_list$model_args_list$prev_prior_a
                  prev_prior_b <- init_model_outs$model_args_list$model_args_list$prev_prior_b
                  
                  known_num <- 0
                  known_values_indicator_list <- init_model_outs$model_args_list$model_args_list$known_values_indicator_list
                  known_values_list <- init_model_outs$model_args_list$model_args_list$known_values_list
                  
                  ## these vars need conversion to different types compatible w/ Stan
                  corr_force_positive <- as.integer(init_model_outs$model_args_list$model_args_list$corr_force_positive)
                  prior_only <-  as.integer(init_model_outs$model_args_list$model_args_list$prior_only)
                  
                  overflow_threshold <- init_model_outs$model_args_list$model_args_list$overflow_threshold
                  underflow_threshold <- init_model_outs$model_args_list$model_args_list$underflow_threshold
 
                  # for (c in 1:n_class) {
                  #   for (t in 1:n_tests) {
                  #     X[[c]][[t]][1:N, 1:n_covariates_max] <- 1
                  #   }
                  # }
                  
                  n_pops <- 1 
                  
                  print(prior_coeffs_mean_mat)
                  print(prior_coeffs_sd_mat)
 
                  Stan_data_list <- list(N = N,
                                         n_tests = n_tests,
                                         ## y = y,
                                         #####
                                         n_covariates_max = n_covariates_max,
                                         #X = list(X_nd, X_d),
                                         n_covs_per_outcome = n_covariates_per_outcome_mat[1, ],
                                         #####
                                         corr_force_positive = corr_force_positive,
                                         known_num = known_num,
                                         # lb_corr = lb_corr,
                                         # ub_corr = ub_corr,
                                         overflow_threshold =  overflow_threshold,
                                         underflow_threshold = underflow_threshold,
                                         #### priors
                                         prior_only = prior_only,
                                         prior_beta_mean = prior_coeffs_mean_mat[[1]],
                                         prior_beta_sd = prior_coeffs_sd_mat[[1]],
                                         prior_LKJ =  lkj_cholesky_eta)
                  
                  
                  
                  ## Compile using BridgeStan:
                  outs_init_bs_model <- BayesMVP:::init_bs_model_internal( Stan_data_list = Stan_data_list,
                                                                           Stan_model_name = "LC_MVP_bin_PartialLog_v5.stan")
                  ## PO_MVP_bin.stan
                  ## LC_MVP_bin_w_mnl_cpp_grad_v1.stan
                  bs_model <- outs_init_bs_model$bs_model
                  json_file_path <- outs_init_bs_model$json_file_path         
                  Stan_model_file_path <- outs_init_bs_model$Stan_model_file_path  
                  ##
                  dummy_json_file_path <- json_file_path
                  dummy_model_so_file <- BayesMVP:::transform_stan_path(Stan_model_file_path)
          
        } else if (Model_type == "latent_trait") {
          
                  overflow_threshold <- init_model_outs$model_args_list$model_args_list$overflow_threshold
                  underflow_threshold <- init_model_outs$model_args_list$model_args_list$underflow_threshold
                  
                  LT_b_priors_shape <- init_model_outs$model_args_list$model_args_list$LT_b_priors_shape
                  LT_b_priors_scale <- init_model_outs$model_args_list$model_args_list$LT_b_priors_scale
                  LT_known_bs_values <- init_model_outs$model_args_list$model_args_list$LT_known_bs_values
                  LT_known_bs_indicator <- init_model_outs$model_args_list$model_args_list$LT_known_bs_indicator
                  
                  prior_coeffs_mean_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_mean_mat
                  prior_coeffs_sd_mat <- init_model_outs$model_args_list$model_args_list$prior_coeffs_sd_mat
 
                  prev_prior_a <- init_model_outs$model_args_list$model_args_list$prev_prior_a
                  prev_prior_b <- init_model_outs$model_args_list$model_args_list$prev_prior_b
 
                  ## these vars need conversion to different types compatible w/ Stan
                  corr_force_positive <- as.integer(init_model_outs$model_args_list$model_args_list$corr_force_positive)
                  prior_only <-  as.integer(init_model_outs$model_args_list$model_args_list$prior_only)
  
                  ##  print(Stan_data_list)
                  
                  n_pops <- 1
                  
                  print(LT_b_priors_shape) ; print(LT_b_priors_scale) ; 
                  print(LT_known_bs_values) ;  print(LT_known_bs_indicator) ; 
                  
                  Phi_type <- init_model_outs$model_args_list$model_args_list$Phi_type
                  print(paste("Phi_type = ", Phi_type))
                  ##
                  Phi_type_int <- ifelse(Phi_type == "Phi", 1, 2)
                  
                  ## Make sure vecs of length 1 are R arrays (1d):
                  if (length(prev_prior_a) == 1) {
                    prior_p_alpha <- ifelse(is.array(prev_prior_a), prev_prior_a, array(prev_prior_a))
                    prior_p_beta <-  ifelse(is.array(prev_prior_b), prev_prior_b, array(prev_prior_b))
                  }
                  
                  Stan_data_list <- list(N = N,
                                         n_tests = n_tests,
                                         y = y,
                                         n_class = 2,
                                         n_pops =  1,  ## multi-pop not supported yet (currently only in Stan version)
                                         pop =   (rep(1, N)),
                                         corr_force_positive = corr_force_positive,
                                         overflow_threshold =  overflow_threshold,
                                         underflow_threshold = underflow_threshold,
                                         #### priors
                                         prior_only = prior_only,
                                         prior_beta_mean = prior_coeffs_mean_mat,
                                         prior_beta_sd = prior_coeffs_sd_mat,
                                         ####
                                         LT_b_priors_shape = LT_b_priors_shape,
                                         LT_b_priors_scale = LT_b_priors_scale,
                                         LT_known_bs_values = LT_known_bs_values,
                                         LT_known_bs_indicator = LT_known_bs_indicator, 
                                         ####
                                         prior_p_alpha =  prior_p_alpha,
                                         prior_p_beta = prior_p_beta)
         
                  ## Compile using BridgeStan
                  outs_init_bs_model <- BayesMVP:::init_bs_model_internal( Stan_data_list = Stan_data_list,
                                                                           Stan_model_name = "PO_latent_trait_bin.stan")
                  ## PO_latent_trait_bin.stan
                  ## latent_trait_w_mnl_cpp_grad_v1.stan
                  bs_model <- outs_init_bs_model$bs_model
                  json_file_path <- outs_init_bs_model$json_file_path         
                  Stan_model_file_path <- outs_init_bs_model$Stan_model_file_path  
                  ##
                  dummy_json_file_path <- json_file_path
                  dummy_model_so_file <- BayesMVP:::transform_stan_path(Stan_model_file_path)
       
          
        }
          
          
                # re-compile the user-supplied Stan model to extract model methods and make init's vector and JSON data file 
          
      
          param_names <- NULL
   
 
          
                #### bs_model <- NULL
                #### bs_model <- init_model_outs$bs_model
                ##
                bs_names  <-  (bs_model$param_names())
                ## if using bs names for param names:
                param_names <- bs_names
                ## Inits:
                inits_unconstrained_vec_per_chain <- list()
                for (kk in 1:n_chains_burnin) {
                  ## NOTE: param_unconstrain_json() Returns a vector of unconstrained parameters * given the constrained parameters *
                  json_string_for_inits_chain_kk <- BayesMVP:::convert_Stan_data_list_to_JSON(init_lists_per_chain[[kk]])
                  validated_json_string <- paste(readLines(json_string_for_inits_chain_kk), collapse="")
                  ##
                  inits_unconstrained_vec_per_chain[[kk]] <- bs_model$param_unconstrain_json(validated_json_string) ## error here
                  ####  inits_unconstrained_vec_per_chain[[kk]] <-  convert_JSON_string_to_R_vector(json_string_for_inits_chain_kk)
                }
                

    
    
  } ### /// end of "if not Stan model [else]" block
  
  
   model_so_file <- BayesMVP:::transform_stan_path(Stan_model_file_path)
  
  
 
    ##  --------------------------   End of starting values 
    n_params <- n_nuisance + n_params_main
    index_nuisance <- 1:n_nuisance
    index_main <- (1 + n_nuisance):n_params
    ##
    print(paste("n_chains_burnin = ", n_chains_burnin))
    ##
    print(paste("n_params_main = ", n_params_main))
    print(paste("n_nuisance = ", n_nuisance))
    print(paste("n_params = ", n_params))
    ##
    theta_main_vectors_all_chains_input_from_R  <- array(0, dim = c(n_params_main, n_chains_burnin))
    theta_us_vectors_all_chains_input_from_R    <- array(0, dim = c(n_nuisance, n_chains_burnin))
  
    for (kk in 1:n_chains_burnin) {
      theta_main_vectors_all_chains_input_from_R[, kk] <-     inits_unconstrained_vec_per_chain[[kk]][index_main]
      theta_us_vectors_all_chains_input_from_R[, kk] <-     inits_unconstrained_vec_per_chain[[kk]][index_nuisance]
    }
  
    json_file_path <- normalizePath(json_file_path)
    model_so_file <- normalizePath(model_so_file)
    
    Model_args_as_Rcpp_List <- init_model_outs$Model_args_as_Rcpp_List
    ##
    Model_args_as_Rcpp_List$model_so_file <- model_so_file
    Model_args_as_Rcpp_List$json_file_path <- json_file_path
  
  return(list(  cmdstanr_model_fit_obj = mod,
                bs_model = bs_model,
                json_file_path = json_file_path,
                model_so_file = model_so_file,
                dummy_json_file_path = dummy_json_file_path,
                dummy_model_so_file = dummy_model_so_file,
                param_names = param_names,
                Stan_data_list = Stan_data_list,
                Stan_model_file_path = Stan_model_file_path,
                Model_args_as_Rcpp_List = Model_args_as_Rcpp_List,
                inits_unconstrained_vec_per_chain = inits_unconstrained_vec_per_chain,
                init_lists_per_chain = init_lists_per_chain,
                theta_main_vectors_all_chains_input_from_R = theta_main_vectors_all_chains_input_from_R, 
                theta_us_vectors_all_chains_input_from_R = theta_us_vectors_all_chains_input_from_R))
  
}




