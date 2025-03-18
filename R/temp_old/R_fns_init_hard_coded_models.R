

#' init_hard_coded_model
#' @keywords internal
#' @export
init_hard_coded_model <- function(Model_type, 
                                  y,
                                  N = N,
                                  model_args_list = list(),
                                  ...) {
  
  hard_coded_models_vec <- c("LC_MVP", "MVP", "latent_trait")
  
  if (Model_type == "MVP") { 
    n_class <- 1
  } else { 
    n_class <- 2
  }
  
  print(paste("n_class = ", n_class))
  
  if (!(Model_type %in% hard_coded_models_vec)) { 
    stop("If not using a Stan model (i.e., setting 'Model_type' to 'Stan'), \n 
         then please use one of the following: 'LC_MVP', 'MVP', or 'latent_trait'. \n
         However, if using a Stan model, then please do not use 'init_hard_coded_model()' - use 'init_Stan_model()' instead")
  }
  
  # if (Model_type == "latent_trait") { 
  #   Model_type <- "latent_trait"
  # }
  
  if (is.null(y)) { 
    stop("no data (y) supplied.")
  }
  
  if (!is.matrix(y)) { 
    stop("y must be a matrix where #cols = #outcomes and #rows = #individuals")
  }
  
 
 ## N <- nrow(y)
  n_tests <- ncol(y)
  n_obs <- N * n_tests 
   
  ## load fn args 
  ### args relevant for all 3 models (i.e. MVP, LC_MVP and latent_trait)
  prior_only <- FALSE ## prior_only <- model_args_list$prior_only ## currently not available
  prior_coeffs_mean_mat <- model_args_list$prior_coeffs_mean_mat
  prior_coeffs_sd_mat <- model_args_list$prior_coeffs_sd_mat
  
  vect_type <- model_args_list$vect_type
  nuisance_transformation <- "Phi" #  not modifiable 
  
  num_chunks <-  model_args_list$num_chunks
  Phi_type <- model_args_list$Phi_type
  
  overflow_threshold <-   model_args_list$overflow_threshold  # not modifiable
  underflow_threshold <-  model_args_list$underflow_threshold # not modifiable
  
  ### Correlation matrix (Omega) args - for LC_MVP and MVP only. 
  corr_force_positive <- model_args_list$corr_force_positive
  corr_param <- "Sean" #  model_args_list$corr_param #  not currently modifiable 
  
  lkj_cholesky_eta <- model_args_list$lkj_cholesky_eta
  corr_prior_norm <- FALSE # model_args_list$corr_prior_norm #  not currently modifiable 
  corr_prior_beta <- FALSE # model_args_list$corr_prior_beta #  not currently modifiable 
  
  ub_corr <- model_args_list$ub_corr
  lb_corr <- model_args_list$lb_corr
  
  known_values_indicator_list <- model_args_list$known_values_indicator_list
  known_values_list <- model_args_list$known_values_list
  
  ### latent-class specific args (i.e. for LC_MVP and latent_trait only)
  prev_prior_a <- model_args_list$prev_prior_a
  prev_prior_b <- model_args_list$prev_prior_b
  
  prior_for_skewed_LKJ_a <- NULL ## prior_for_skewed_LKJ_a <- model_args_list$prior_for_skewed_LKJ_a ## currently not available
  prior_for_skewed_LKJ_b <- NULL ## prior_for_skewed_LKJ_b <- model_args_list$prior_for_skewed_LKJ_b ## currently not available
  prior_for_corr_a <- NULL ## prior_for_corr_a <- model_args_list$prior_for_corr_a ## currently not available
  prior_for_corr_b <- NULL ## prior_for_corr_b <- model_args_list$prior_for_corr_b ## currently not available
  
  ### latent_trait-specific args
  LT_b_priors_shape <- model_args_list$LT_b_priors_shape
  LT_b_priors_scale <- model_args_list$LT_b_priors_scale
  LT_known_bs_values <- model_args_list$LT_known_bs_values
  LT_known_bs_indicator <- model_args_list$LT_known_bs_indicator
  
  ### covariate args (at the moment MVP and LC_MVP only)
  n_covariates_per_outcome_mat <- model_args_list$n_covariates_per_outcome_mat
  X <- model_args_list$X
  
  if ((is.null(X)) || (is.null(n_covariates_per_outcome_mat))) { 
          
          warning("No covariates (i.e., X) supplied - will assume intercept-only")
    
          n_covariates_per_outcome_mat <- array(1, dim = c(n_class, n_tests))
          
          ## then make dummy (intercept-only) X
          X_per_class <- array(1, dim = c(n_tests, 1, N))
          X_list <- list()
          
          for (c in 1:n_class) {
            X_list[[c]] <- list()
            for (t in 1:n_tests) {
              X_list[[c]][[t]] <- array(999999, dim = c(N, n_covariates_per_outcome_mat[c, t] ))
              for (k in 1:n_covariates_per_outcome_mat[c, t] ) {
                for (n in 1:N) {
                  X_list[[c]][[t]][n,  k] <- X_per_class[t, k, n]
                }
              }
            }
          }
          
          
          X <- X_list
          
       
          
          
  }
  

  if (is.null(overflow_threshold)) { 
    overflow_threshold <- +5
    model_args_list$overflow_threshold <- +5
  }
  if (is.null(underflow_threshold)) { 
    underflow_threshold <- -5
    model_args_list$underflow_threshold <- -5
  }
  
  if (is.null(Phi_type)) { 
    ##
    Phi_type <- "Phi"
    model_args_list$Phi_type <- "Phi"
    ##
    inv_Phi_type <- "inv_Phi"
    model_args_list$inv_Phi_type <- "inv_Phi"
    ##
  
  }
  if (is.null(nuisance_transformation)) { 
    nuisance_transformation <- "Phi"
    model_args_list$nuisance_transformation <- "inv_Phi"
  }
  
  
  ### call C++ function to detect is user has AVX-2 or AVX-512 vectorisation support 
  if (is.null(vect_type)) { 
    vect_type <-  BayesMVP:::detect_vectorization_support()
    model_args_list$vect_type <-  BayesMVP:::detect_vectorization_support()
  }
  
  
  ### find "optimal" number of chunks (at the moment only looks at number of cores, but should ideally also look at L3 cache and memory bandwidth)
  if (is.null(num_chunks)) { 
    num_chunks <- BayesMVP:::find_num_chunks_MVP(N, n_tests)
    model_args_list$num_chunks <- num_chunks
  }
  
  
  if ( (Model_type == "LC_MVP") || (Model_type == "latent_trait") ) { 
    n_class = 2
    model_args_list$n_class <- 2
  } else if (Model_type == "MVP") { 
    n_class = 1
    model_args_list$n_class <- 1
  }
   
  if (n_class == 2) { 
    if (is.null(prev_prior_a)) { 
      prev_prior_a <- 1
      model_args_list$prev_prior_a <- 1
    }
    if (is.null(prev_prior_b)) { 
      prev_prior_b <- 1
      model_args_list$prev_prior_b <- 1
    }
  }
  
  if (is.null(corr_force_positive)) {
    corr_force_positive <- FALSE
    model_args_list$corr_force_positive <- FALSE
  }
  
  if (is.null(lb_corr)) {
    if (corr_force_positive == TRUE) lb <- 0 
    else lb <- -1
    lb_corr <- list()
    for (c in 1:n_class) {
      lb_corr[[c]] <- array(lb, dim = c(n_tests, n_tests))
    }
    model_args_list$lb_corr <- lb_corr
  }
  
  if (is.null(ub_corr)) {
    ub_corr <- list()
    for (c in 1:n_class) {
      ub_corr[[c]]  = array( +1, dim = c(n_tests, n_tests))
    }
    model_args_list$ub_corr <- ub_corr
  }
  
 
  if (is.null(corr_param)) {
    corr_param <- "Sean"
    model_args_list$corr_param <- "Sean"
  }
  
  if (Model_type %in% c("MVP", "LC_MVP")) {
    if (is.null(lkj_cholesky_eta)) {
      warning("lkj_cholesky_eta not supplied - if using LKJ priors (the default for LC-MVP and MVP) then will assume LKJ(2) in all classes")
      lkj_cholesky_eta <- matrix(rep(2, n_class), ncol = 1)
      model_args_list$lkj_cholesky_eta <- matrix(rep(2, n_class), ncol = 1)
    }
  }
  
  

    
  if ( (is.null(n_covariates_per_outcome_mat)) && (Model_type == "Stan") ) {  # if Stan model used, then make dummy variable
  
      n_covariates_per_outcome_mat <- c(array(1, dim = c(n_tests, 1)))
      model_args_list$n_covariates_per_outcome_mat <- c(array(1, dim = c(n_tests, 1)))
  }

  
  
  if (is.null(prior_only)) { 
    prior_only <- FALSE
    model_args_list$prior_only <- FALSE
  }
  
 
  
  # set lp_and_grad function arguments needed and put them in a list, so that they can be passed on to the generic EHMC-ADAM functions
  
  {
        
        n_covariates_max <- max(unlist(n_covariates_per_outcome_mat))
        n_covariates_total <- sum(unlist(n_covariates_per_outcome_mat))
        
        if (Model_type == "LC_MVP") {
          
              n_corrs <- n_class * 0.5 * n_tests * (n_tests - 1)
              n_params_main <- (n_class - 1)  +  n_corrs +  n_covariates_total 
          
        } else if (Model_type == "MVP") {
          
              n_corrs <- 0.5 * n_tests * (n_tests - 1)
              n_params_main <- n_corrs +  n_covariates_total 
          
        } else if (Model_type == "latent_trait") {
          
              LT_n_bs <- n_tests * n_class
              n_corrs <- LT_n_bs
              n_params_main <- (n_class - 1)  +   n_corrs   +  n_covariates_total 
          
        }
          
          
        n_nuisance <-  N * n_tests 
        n_params <- n_params_main + n_nuisance 
        index_us <- 1:n_nuisance
        index_main <- (n_nuisance+1):n_params
        index_corrs <- (n_nuisance+1):(n_nuisance+n_corrs)
        n_us <- n_nuisance
    
  }
  
 
  
  
  if (is.null(prior_coeffs_mean_mat)) { 
        prior_coeffs_mean_mat <-  list() 
        mat_1 <- array(-1, dim = c(n_covariates_max, n_tests))
        mat_2 <- array(+1, dim = c(n_covariates_max, n_tests))
        if (n_class == 2) { 
          for (k in 1:n_covariates_per_outcome_mat[1, t]) {
            prior_coeffs_mean_mat[[1]]  <- mat_1
          }
          for (k in 1:n_covariates_per_outcome_mat[2, t]) {
            prior_coeffs_mean_mat[[2]]  <- mat_2
          }
        } else { 
          prior_coeffs_mean_mat[[1]] <- array(0, dim = c(n_covariates_max, n_tests))
        }
  }
  

                                           
  if (is.null(prior_coeffs_sd_mat)) { 
        prior_coeffs_sd_mat <- list() #  array(1, dim = c(n_class, n_tests, n_covariates_max))
        for (c in 1:n_class) { 
          prior_coeffs_sd_mat[[c]] <- array(1, dim = c(n_covariates_max, n_tests))
        }
  }
  
  
  
  
  # homog_corr = homog_corr # add in future 
 
  
  {
    
 
    
    list_prior_for_corr_a <-  list_prior_for_corr_b <- list()
    
    if (is.null(prior_for_corr_a)) { 
      prior_for_corr_a <- array(1, dim = c(n_tests, n_tests, n_class))
      prior_for_corr_b <- array(1, dim = c(n_tests, n_tests, n_class))
    }
    
    
    for (c in 1:n_class) {
      list_prior_for_corr_a[[c]] <- prior_for_corr_a[,,c]
      list_prior_for_corr_b[[c]] <- prior_for_corr_b[,,c]
    }
    

    
    
    if (is.null(lb_corr)) { 
      lb_corr <- list()
      for (c in 1:n_class) { 
        lb_corr[[c]] <- array(-1.0, dim = c(n_tests, n_tests))
      }
    }
    
    if (is.null(lb_corr)) {
      for (c in 1:n_class) {                                                                                                                           
        if (corr_force_positive == TRUE)    {                                                                                                         
          lb_corr[[c]][, ] <- 0    
        }  else  { 
          lb_corr[[c]][, ] <- -1  
        }
      }
    }
    
    
    if (is.null(ub_corr)) { 
      ub_corr <- list()
      for (c in 1:n_class) { 
        ub_corr[[c]] <- array(+1.0, dim = c(n_tests, n_tests))
      }
    }
    
 
    
    
    if (Model_type == "latent_trait")   corr_param <- "latent_trait"
    
    
    if (is.null(known_values_indicator_list)) {                                                                                                      
      known_values_indicator_list <- list()                                                                                                        
      known_values_list <- list()                                                                                                                  
      for (c in 1:n_class) {                                                                                                                       
        known_values_indicator_list[[c]] <- diag(n_tests)                                                                                        
        known_values_list[[c]] <- diag(n_tests)                                                                                                  
      }      
    }
    
    
    if (is.null(LT_b_priors_shape)) {  
      LT_b_priors_shape <- array(1, dim = c(n_class, n_tests))
    }
    if (is.null(LT_b_priors_scale)) {  
      LT_b_priors_scale <- array(1, dim = c(n_class, n_tests))
    }
    
    
    if (is.null(LT_known_bs_indicator)) {  
      LT_known_bs_indicator <- array(0, dim = c(n_class, n_tests))
    }
    if (is.null(LT_known_bs_values)) {  
      LT_known_bs_values    <- array(0.00001, dim = c(n_class, n_tests))
    }
    
    
    
    known_corr_index_vec <- rep(0, n_corrs)
    
    if (corr_param != "latent_trait") {
      at_least_one_corr_known_indicator <- 0 
      counter <- 1
      corr_index = n_us + 1
      for (c in 1:n_class) {
        for (i in 2:n_tests) {
          for (j in 1:(i-1)) {
            if (known_values_indicator_list[[c]][i, j] == 1)   {
              known_corr_index_vec[counter] <- corr_index
              at_least_one_corr_known_indicator <- 1
            }
            corr_index <- corr_index + 1
            counter <- counter + 1
          }
        }
      }
    } else { 
      at_least_one_corr_known_indicator <- 0 
      counter <- 1
      corr_index = n_us + 1
      for (c in 1:n_class) {
        for (i in 1:n_tests) {
          if (LT_known_bs_indicator[c, i] == 1)   {
            known_corr_index_vec[counter] <- corr_index
            at_least_one_corr_known_indicator <- 1
          }
          corr_index <- corr_index + 1
          counter <- counter + 1
        }
      }
      
      known_corr_index_vec <- (n_us + 1 + n_class*n_tests):(n_us + n_corrs) # [(counter:n_corrs)] 
      
    }
    
  
    
    index_subset_excl_known <- index_main
    for (i_known in 1:length(known_corr_index_vec)) {
      index_subset_excl_known <-   index_subset_excl_known[index_subset_excl_known != known_corr_index_vec[i_known]]
    }
    
    
    
    if (n_class == 1) {

      try({
        if (length(lkj_cholesky_eta) > 1)  lkj_cholesky_eta = array(lkj_cholesky_eta[1])
        if (length(known_values_indicator_list) > 1)  known_values_indicator_list  = list(known_values_indicator_list[[1]])
        if (length(known_values_list) > 1)  known_values_list  = list(known_values_list[[1]])
        prev_prior_a = 999999
        prev_prior_b = 999999

        if (dim(LT_b_priors_shape)[1] > 1) LT_b_priors_shape = LT_b_priors_shape[1,]
        if (dim(LT_b_priors_scale)[1] > 1) LT_b_priors_scale = LT_b_priors_scale[1,]
        if (dim(LT_known_bs_indicator)[1] > 1) LT_known_bs_indicator = LT_known_bs_indicator[1,]
        if (dim(LT_known_bs_values)[1] > 1) LT_known_bs_values = LT_known_bs_values[1,]
      })

    } else {

    }


    
    
    model_args_list$n_covariates_per_outcome_mat <- n_covariates_per_outcome_mat
    model_args_list$prior_only <- prior_only ### not exposed/modifiable
    model_args_list$prior_coeffs_mean_mat <- prior_coeffs_mean_mat
    model_args_list$prior_coeffs_sd_mat <- prior_coeffs_sd_mat
    model_args_list$corr_force_positive <- corr_force_positive
    model_args_list$corr_param <- corr_param
    model_args_list$lkj_cholesky_eta <- lkj_cholesky_eta
    model_args_list$corr_prior_norm <- corr_prior_norm
    model_args_list$corr_prior_beta <- corr_prior_beta
    model_args_list$ub_corr <- ub_corr
    model_args_list$lb_corr <- lb_corr
    model_args_list$prev_prior_a <- prev_prior_a
    model_args_list$prev_prior_b <- prev_prior_b
    model_args_list$prior_for_skewed_LKJ_a <- prior_for_skewed_LKJ_a
    model_args_list$prior_for_skewed_LKJ_b <- prior_for_skewed_LKJ_b
    model_args_list$prior_for_corr_a <- prior_for_corr_a
    model_args_list$prior_for_corr_b <- prior_for_corr_b
    model_args_list$known_values_indicator_list <- known_values_indicator_list
    model_args_list$known_values_list <- known_values_list
    model_args_list$LT_b_priors_shape <- LT_b_priors_shape
    model_args_list$LT_b_priors_scale <- LT_b_priors_scale
    model_args_list$LT_known_bs_values <- LT_known_bs_values
    model_args_list$LT_known_bs_indicator <- LT_known_bs_indicator
    # model_args_list$init_coeffs <- init_coeffs
    # model_args_list$init_raw_corrs <- init_raw_corrs
    # model_args_list$init_prev <- init_prev
    # model_args_list$same_init_per_chain <- same_init_per_chain
    model_args_list$overflow_threshold <- overflow_threshold ### not exposed/modifiable
    model_args_list$underflow_threshold <- underflow_threshold ### not exposed/modifiable
    
  }
  
  
  { ### dummy Stan stuff for C++ (loads dummy Stan model whenever user isn't using a Stan model) 
    

  
  }
  
 
  
  
  {
    
    
    N_obs <- N*n_tests
    
    
    handle_numerical_issues <- T
    
    vect_type <-  vect_type
    
    vect_type_exp =        vect_type
    vect_type_tanh =       vect_type
    vect_type_log =        vect_type
    vect_type_lse =        vect_type
    vect_type_Phi =        vect_type
    vect_type_inv_Phi =    vect_type
    vect_type_log_Phi =     vect_type
    vect_type_inv_Phi_approx_from_logit_prob =     vect_type
    
    skip_checks  = FALSE
    
    Model_args_bools <-   list()
    Model_args_ints <-    list()
    Model_args_doubles <- list()
    Model_args_strings <- list()
    
    Model_args_col_vecs_double <- list()
    Model_args_mats_double <-     list()
    
    Model_args_vecs_of_mats_double <- list()
    Model_args_vecs_of_mats_int <- list()
    Model_args_vecs_of_col_vecs_int <- list()
    
    Model_args_mats_int <- list()
    
    corr_prior_beta <- FALSE
    exclude_priors <-  FALSE
    
    ## bools
    Model_args_bools[[1]] <- exclude_priors
    Model_args_bools[[2]] <- FALSE # CI
    Model_args_bools[[3]] <- corr_force_positive
    Model_args_bools[[4]] <- FALSE # corr_prior_beta
    Model_args_bools[[5]] <- FALSE #  corr_prior_norm
    Model_args_bools[[6]] <- TRUE  #  handle_numerical_issues
    Model_args_bools[[7]] <-  skip_checks
    Model_args_bools[[8]] <-  skip_checks
    Model_args_bools[[9]] <-  skip_checks
    Model_args_bools[[10]] <- skip_checks
    Model_args_bools[[11]] <- skip_checks
    Model_args_bools[[12]] <- skip_checks
    Model_args_bools[[13]] <- skip_checks
    Model_args_bools[[14]] <- skip_checks
    Model_args_bools[[15]] <- FALSE ### debug
    
    Model_args_bools <- c(unlist(Model_args_bools))
    
    ## ints
    Model_args_ints[[1]] <- 1 # n_cores (dummy)
    Model_args_ints[[2]] <- n_class
    Model_args_ints[[3]] <- 5 # dummy 
    Model_args_ints[[4]] <- num_chunks
    
    Model_args_ints <- c(unlist(Model_args_ints))
    Model_args_ints <- as.integer(Model_args_ints)
    
    ## doubles
    Model_args_doubles[[1]] <- prev_prior_a
    Model_args_doubles[[2]] <- prev_prior_b
    Model_args_doubles[[3]] <- overflow_threshold
    Model_args_doubles[[4]] <- underflow_threshold # underflow_threshold
    
    Model_args_doubles <- c(unlist(Model_args_doubles))
    
    ## strings
    Model_args_strings[[1]] <- vect_type
    
    Model_args_strings[[2]] <- Phi_type
    Model_args_strings[[3]] <- inv_Phi_type
    
    Model_args_strings[[4]] <- vect_type_exp
    Model_args_strings[[5]] <- vect_type_log
    Model_args_strings[[6]] <- vect_type_lse
    Model_args_strings[[7]] <- vect_type_tanh
    Model_args_strings[[8]] <- vect_type_Phi
    Model_args_strings[[9]] <- vect_type_log_Phi
    Model_args_strings[[10]] <- vect_type_inv_Phi
    Model_args_strings[[11]] <- vect_type_inv_Phi_approx_from_logit_prob
    Model_args_strings[[12]] <- "all" # sampling_mode
    Model_args_strings[[13]] <- nuisance_transformation
    
    Model_args_strings <- c(unlist(Model_args_strings))
    
    
    ## other
    if (Model_type  %in% c("MVP", "LC_MVP")) {
          Model_args_col_vecs_double[[1]] <-  (lkj_cholesky_eta)
    } else { 
          Model_args_col_vecs_double[[1]] <-   matrix(rep(1, 10), ncol = 1) # dummy var
    }
    
    
    ## other
    Model_args_mats_double[[1]] <- LT_b_priors_shape
    Model_args_mats_double[[2]] <- LT_b_priors_scale
    Model_args_mats_double[[3]] <- LT_known_bs_indicator
    Model_args_mats_double[[4]] <- LT_known_bs_values
    
    
    ## other
    Model_args_vecs_of_mats_double[[1]] <- prior_coeffs_mean_mat
    Model_args_vecs_of_mats_double[[2]] <- prior_coeffs_sd_mat
    Model_args_vecs_of_mats_double[[3]] <- list_prior_for_corr_a
    Model_args_vecs_of_mats_double[[4]] <- list_prior_for_corr_b
    Model_args_vecs_of_mats_double[[5]] <- lb_corr
    Model_args_vecs_of_mats_double[[6]] <- ub_corr
    Model_args_vecs_of_mats_double[[7]] <- known_values_list
    
    
    ## other
    Model_args_vecs_of_mats_int[[1]] <- known_values_indicator_list
    
    # ## other
    # n_covariates_per_outcome_mat[[1]] <- matrix( n_covariates_per_outcome_mat[[1]])
    # if (n_class > 1) {
    #    n_covariates_per_outcome_mat[[2]] <- matrix( n_covariates_per_outcome_mat[[2]])
    # }
    
   ###  Model_args_vecs_of_col_vecs_int[[1]] <- n_covariates_per_outcome_mat
    
    Model_args_2_later_vecs_of_mats_double <- list(X)
    
    Model_args_mats_int[[1]] <- n_covariates_per_outcome_mat
    
    if (is.matrix(n_covariates_per_outcome_mat) == FALSE) {
       n_covariates_per_outcome_mat <- matrix(n_covariates_per_outcome_mat)
       Model_args_mats_int[[1]] <- n_covariates_per_outcome_mat
    }
    
    if (Model_type  %in% c("MVP", "LC_MVP")) {
        if (is.matrix(lkj_cholesky_eta) == FALSE) {
          lkj_cholesky_eta <- matrix(lkj_cholesky_eta)
        }
    }
    
    ### ---------------------------------------------- Put all in a big list
    Model_args_bools <- matrix(Model_args_bools)
    Model_args_ints <- matrix(Model_args_ints)
    Model_args_doubles <- matrix(Model_args_doubles)
    Model_args_strings <- matrix(Model_args_strings)
    
    Model_args_bools[[1]] <- exclude_priors
    Model_args_bools[[2]] <- FALSE # CI
    Model_args_bools[[3]] <- corr_force_positive
    Model_args_bools[[4]] <- FALSE # corr_prior_beta
    Model_args_bools[[5]] <- FALSE #  corr_prior_norm
    Model_args_bools[[6]] <- TRUE  #  handle_numerical_issues
    Model_args_bools[[7]] <-  skip_checks
    Model_args_bools[[8]] <-  skip_checks
    Model_args_bools[[9]] <-  skip_checks
    Model_args_bools[[10]] <- skip_checks
    Model_args_bools[[11]] <- skip_checks
    Model_args_bools[[12]] <- skip_checks
    Model_args_bools[[13]] <- skip_checks
    Model_args_bools[[14]] <- skip_checks
    Model_args_bools[[15]] <- FALSE ### debug
    
    colnames(Model_args_bools) <- c("Model arguments - boolean")
    rownames(Model_args_bools) <- c("exclude_priors", "Cond. indep.", "Force +'ve corr's",
                                    "corr_prior_beta", "corr_prior_norm",
                                    "handle num. issues",
                                    "skip_checks_exp", "skip_checks_log", "skip_checks_lse", "skip_checks_tanh",
                                    "skip_checks_Phi", "skip_checks_log_Phi", "skip_checks_inv_Phi",
                                    "skip_checks_inv_Phi_approx_from_logit_prob",
                                    "debug")

    colnames(Model_args_ints) <- c("Model arguments - integers")
    rownames(Model_args_ints) <- c("n_cores", "n_class", "ub_threshold_phi_approx", "number of chunks")

    colnames(Model_args_doubles) <- c("Model arguments - doubles")
    rownames(Model_args_doubles) <- c("prev_prior_a", "prev_prior_b", "overflow_threshold", "underflow_threshold")

    colnames(Model_args_strings) <- c("Model arguments - character strings")
    rownames(Model_args_strings) <- c("vect_type",
                                      "Phi_type", "inv_Phi_type",
                                      "vect_type_exp", "vect_type_log", "vect_type_lse", "vect_type_tanh",
                                      "vect_type_Phi", "vect_type_log_Phi", "vect_type_inv_Phi",
                                      "vect_type_inv_Phi_approx_from_logit_prob",
                                      "NA", "nuisance_transformation")

    # 
      # names(Model_args_col_vecs_double)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_mats_double)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_mats_int)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_vecs_of_col_vecs_int)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_vecs_of_mats_double[1]) <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_vecs_of_mats_int)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
      # names(Model_args_mats_double)[1] <- c("LKJ prior for corr matrices (MVP/LC_MVP only")
    
    print(paste("n_params_main = ", n_params_main))
    
  #    print(Model_args_mats_int)
    
    Model_args_as_Rcpp_List <- list(
      N = N,
      n_nuisance = n_nuisance,
      n_params_main = n_params_main,
      #  rstan_model = rstan_model,
      #  rstan_model_fixed_nuisance = rstan_model_fixed_nuisance,
      #  rstan_model_fixed_main = rstan_model_fixed_main,
      Model_args_bools               = Model_args_bools,                 
      Model_args_ints                = Model_args_ints,                  
      Model_args_doubles             = Model_args_doubles,              
      Model_args_strings             = Model_args_strings,   
      Model_args_col_vecs_double      = Model_args_col_vecs_double,   
      # Model_args_col_vecs_int
      Model_args_mats_double          = Model_args_mats_double,
      Model_args_mats_int = Model_args_mats_int,
      #Model_args_vecs_of_col_vecs_double
      Model_args_vecs_of_col_vecs_int = Model_args_vecs_of_col_vecs_int,
      Model_args_vecs_of_mats_double = Model_args_vecs_of_mats_double,
      Model_args_vecs_of_mats_int = Model_args_vecs_of_mats_int,
      # Model_args_2_later_vecs_of_col_vecs_double
      # Model_args_2_later_vecs_of_col_vecs_int
      # Model_args_2_later_vecs_of_mats_double = Model_args_2_later_vecs_of_mats_double,
      # Model_args_2_later_vecs_of_mats_int    
      Model_args_2_later_vecs_of_mats_double  = Model_args_2_later_vecs_of_mats_double
    )
    
    
    Model_args_as_Rcpp_List$model_so_file <- "none"
    Model_args_as_Rcpp_List$json_file_path <- "none"
    
    
    ### MVP_fn_args_for_cpp <- MVP_fn_args_convert_R_List_to_MVP_fn_args_struct(FALSE, "MVP_LC", TRUE, FALSE, theta_main_array, theta_us_array, y, X, MVP_fn_args, 10, 0.01, M_dense_main, M_inv_dense_main, M_inv_dense_main_chol)
    
    
    
    
  }
  
  
  
  
  
  
  return(list(Model_args_as_Rcpp_List = Model_args_as_Rcpp_List, 
              model_args_list = model_args_list,
              Model_type = Model_type, 
              y = y, 
              X = X))
 
  
  
}







