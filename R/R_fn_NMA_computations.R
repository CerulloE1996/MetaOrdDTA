




#' generate_performance_metadata
#' @keywords internal
#' @export
generate_performance_metadata <- function(n_thr, 
                                          test_names = NULL) {
  
          n_tests <- length(n_thr)
          
          if (is.null(test_names)) {
            test_names <- paste0("Test", 1:n_tests)
          }
          
          metadata <- tibble()
          
          for (t in 1:n_tests) {
            for (k in 1:n_thr[t]) {
              metadata <- bind_rows(
                metadata,
                tibble(
                  test = t,
                  test_name = test_names[t],
                  threshold = k
                )
              )
            }
          }
          
          return(metadata)
  
}





#' R_fn_extract_NMA_AUC_summary
#' @keywords internal
#' @export
R_fn_extract_NMA_AUC_summary <- function( tibble_gq,
                                          network,
                                          n_thr,
                                          test_names = NULL) {
  
        require(dplyr)
        require(tidyr)
        
        # Extract AUC summaries if they exist
        auc_summary <- tibble_gq %>%
          filter(grepl("^AUC\\[", parameter)) %>%
          mutate(
            test = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))  # Fixed: added .*
          ) %>%
          mutate(
            test_name = if(!is.null(test_names)) test_names[test] else paste0("Test", test)
          ) %>%
          select(test, test_name,
                 AUC_mean = mean,
                 AUC_median = `50%`,
                 AUC_lower = `2.5%`,
                 AUC_upper = `97.5%`,
                 AUC_sd = sd)
        
        # Extract AUC pred summaries if they exist
        auc_pred_summary <- tibble_gq %>%
          filter(grepl("^AUC_pred\\[", parameter)) %>%
          mutate(
            test = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))  # Fixed: added .*
          ) %>%
          mutate(
            test_name = if(!is.null(test_names)) test_names[test] else paste0("Test", test)
          ) %>%
          select(test, test_name,
                 AUC_mean = mean,
                 AUC_median = `50%`,
                 AUC_lower = `2.5%`,
                 AUC_upper = `97.5%`,
                 AUC_sd = sd)
        
        auc_diff_summary <- NULL
        if (network) {
          # Extract AUC differences if they exist
          # General solution for any number of tests
          n_tests <- length(test_names)  # Use test_names to get number of tests
          
          auc_diff_summary <-  tibble_gq %>%
            filter(grepl("^AUC_diff\\[", parameter)) %>%
            mutate(
              index = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))  # Fixed: added .*
            ) %>%
            # Generate metadata for AUC differences
            left_join(
              expand.grid(
                test1 = 1:(n_tests-1),
                test2 = 1:n_tests
              ) %>%
                filter(test2 > test1) %>%
                mutate(
                  index = row_number(),
                  test1_name = test_names[test1],
                  test2_name = test_names[test2],
                  comparison = paste0(test1_name, " vs ", test2_name)
                ),
              by = "index"
            ) %>%
            select(test1, test2, test1_name, test2_name,
                   AUC_diff_mean = mean,
                   AUC_diff_median = `50%`,
                   AUC_diff_lower = `2.5%`,
                   AUC_diff_upper = `97.5%`,
                   AUC_diff_sd = sd) %>%
            mutate(
              prob_test1_better = NA_real_  # Fixed: changed NAreal to NA_real_
            )
        }
        
        return(list(
          auc = auc_summary,
          auc_pred = auc_pred_summary,
          auc_diff = auc_diff_summary
        ))
  
}







#' Run standalone generated quantities for user-specified covariates
#' @export
R_fn_GQ_covariates_NMA <- function(baseline_case_nd, 
                                   baseline_case_d,
                                   model_summary_and_trace_obj,
                                   model_prep_obj,
                                   chains = min(8, parallel::detectCores())
) {
          
          require(cmdstanr)
          
          ##
          stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
          ##
          stan_model_file_path     <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_model_file_path
          stan_functions_directory <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_functions_directory
          ##
          stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
  
          # 
          # stan_model_file_path <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/stan_models/stan_models_NMA/DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan"
          # ## DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan
          # ## DTA_NMA_Nyaga_Xu_FIXEDthr_GQ.stan
          
          # Prepare data for GQ block
          gq_data <- stan_data_list
          
          # Add/replace user-specified covariate values
          gq_data$baseline_case_nd <- baseline_case_nd
          gq_data$baseline_case_d  <- baseline_case_d
          
          
          # mat <- array(1.0, dim = c(n_index_tests, n_thr_max + 1))
          # gq_data$prior_dirichlet_phi <- if_null_then_set_to( gq_data$prior_dirichlet_phi, 
          #                                                    list(mat, mat))
          # ##
          # ##
          # prior_kappa_mean_nd <- rep(log(200), n_index_tests)
          # prior_kappa_mean_d  <- rep(log(50), n_index_tests)
          # gq_data$prior_kappa_mean   <- if_null_then_set_to(  gq_data$prior_kappa_mean, 
          #                                                    list(prior_kappa_mean_nd, prior_kappa_mean_d))
          # ##
          # prior_kappa_SD <- rep(0.75, n_index_tests)
          # gq_data$prior_kappa_SD   <- if_null_then_set_to(  gq_data$prior_kappa_SD, 
          #                                                  list(prior_kappa_SD, prior_kappa_SD))
          # ##
          # gq_data$kappa_lb <- if_null_then_set_to( gq_data$kappa_lb, 
          #                                         1.0)
          # ##
          # gq_data$kappa_ub <- if_null_then_set_to( gq_data$kappa_ub, 
          #                                         exp(10))
          ##
          priors_kappa_df <- 15
          gq_data$prior_kappa_df <- if_null_then_set_to( gq_data$prior_kappa_df, 
                                                        list(rep(priors_kappa_df, n_index_tests), 
                                                             rep(priors_kappa_df, n_index_tests)))
          
          
          
          
          # Compile the standalone GQ model
          gq_model <- cmdstanr::cmdstan_model( stan_file = stan_model_file_path,
                                               include_paths = stan_functions_directory)
          
          # Generate from the fitted values
          fitted_draws <- stan_mod_samples$draws(format = "draws_matrix")
          
          
          
          gq_results <- gq_model$generate_quantities(
            fitted_params = fitted_draws,  # Your fitted cmdstanr model
            data = gq_data,
            seed = 123
            # parallel_chains = chains
          )
          
          return(gq_results)
  
}







# Convert BridgeStan notation to Stan notation
#' @export
convert_bridgestan_par_names_to_stan <- function(names) {
  
          # Handle the conversion more carefully
          result <- names
          
          # Find first dot and replace with [
          result <- sub("\\.", "[", result)
          
          # Replace all remaining dots with commas
          result <- gsub("\\.", ",", result)
          
          # Add closing bracket if we added an opening one
          needs_bracket <- grepl("\\[", result) & !grepl("\\]", result)
          result[needs_bracket] <- paste0(result[needs_bracket], "]")
          
          return(result)
  
}




#' @export
permute_3d_draws_array <- function(array,
                                   iter_index = 1,
                                   chain_index = 2,
                                   param_index = 3,
                                   new_iter_index = 2,
                                   new_chain_index = 3,
                                   new_param_index = 1) {
  
          # Create the permutation vector
          # We need to map old positions to new positions
          perm_vector <- numeric(3)
          
          # Find which old index goes to position 1
          if (new_iter_index == 1) perm_vector[1] <- iter_index
          else if (new_chain_index == 1) perm_vector[1] <- chain_index
          else if (new_param_index == 1) perm_vector[1] <- param_index
          
          # Find which old index goes to position 2
          if (new_iter_index == 2) perm_vector[2] <- iter_index
          else if (new_chain_index == 2) perm_vector[2] <- chain_index
          else if (new_param_index == 2) perm_vector[2] <- param_index
          
          # Find which old index goes to position 3
          if (new_iter_index == 3) perm_vector[3] <- iter_index
          else if (new_chain_index == 3) perm_vector[3] <- chain_index
          else if (new_param_index == 3) perm_vector[3] <- param_index
          
          # Apply the permutation
          new_array <- aperm(array, perm = perm_vector)
          
          # Update dimension names if they exist
          if (!is.null(dimnames(array))) {
            old_names <- dimnames(array)
            new_names <- list(
              old_names[[perm_vector[1]]],
              old_names[[perm_vector[2]]],
              old_names[[perm_vector[3]]]
            )
            
            # Rename based on new positions
            names(new_names) <- c(
              if (new_iter_index == 1) "iteration" else if (new_chain_index == 1) "chain" else "parameter",
              if (new_iter_index == 2) "iteration" else if (new_chain_index == 2) "chain" else "parameter",
              if (new_iter_index == 3) "iteration" else if (new_chain_index == 3) "chain" else "parameter"
            )
            
            dimnames(new_array) <- new_names
          }
          
          return(new_array)
  
}





## needs debugging - bookmark
R_fn_BayesMVP_BridgeStan_covariates_NMA <- function( new_cov_data,
                                                     model_prep_obj,
                                                     model_summary_and_trace_obj
) {


          # model_prep_obj <- MR_Model_dummy$model_prep_obj ## testing
          # model_summary_and_trace_obj <- MR_Model_dummy$model_summary_and_trace_obj  ## testing
          ##
          stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
          ##
          stan_model_file_path     <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_model_file_path
          stan_functions_directory <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_functions_directory
          ##
          stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
          # str(stan_mod_samples)
          ##
          trace_main <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_main
          trace_main_dims <- dim(trace_main) ; trace_main_dims
          ##
          trace_tp <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_tp
          trace_tp_dims <- dim(trace_tp) ; trace_tp_dims
          ##
          trace_gq <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_gq
          trace_gq_dims <- dim(trace_gq) ; trace_gq_dims
          ##
          n_params_main <- trace_main_dims[3] ; n_params_main
          n_chains <- trace_main_dims[2] ; n_chains
          n_iter <- trace_main_dims[1] ; n_iter
          ##
          init_lists_per_chain <- model_prep_obj$init_lists_per_chain
          init_lists_per_chain <- resize_init_list(init_lists_per_chain, n_chains)
          ##
          Stan_data_list <- stan_data_list
          ##
          ## ---- Update cov data:
          ##
          Stan_data_list$baseline_case_nd <- new_cov_data$baseline_case_nd
          Stan_data_list$baseline_case_d <- new_cov_data$baseline_case_d
          ##
          # ###  -----------  Compile + initialise the model using "MVP_model$new(...)"
          BayesMVP_model_prep_obj <- BayesMVPtest:::MVP_model$new(   Model_type =  "Stan",
                                                                     y = NULL,
                                                                     N = NULL,
                                                                     Stan_data_list = Stan_data_list,
                                                                     Stan_model_file_path = stan_model_file_path,
                                                                     init_lists_per_chain = init_lists_per_chain,
                                                                     sample_nuisance = FALSE,
                                                                     n_chains_burnin = n_chains,
                                                                     ##
                                                                     n_params_main = n_params_main,
                                                                     n_nuisance = 0,
                                                                     ##
                                                                     stanc_args = stanc_args)
          #
          bs_model <- BayesMVP_model_prep_obj$init_object$bs_model
          bs_names  <-  (bs_model$param_names())
          bs_names_inc_tp <-  (bs_model$param_names(include_tp = TRUE))
          bs_names_inc_tp_and_gq <-  (bs_model$param_names(include_tp = TRUE, include_gq = TRUE))
          ##
          stan_names_main <- convert_bridgestan_par_names_to_stan(bs_names)
          stan_names_inc_tp <- convert_bridgestan_par_names_to_stan(bs_names_inc_tp)
          stan_names_inc_tp_and_gq <- convert_bridgestan_par_names_to_stan(bs_names_inc_tp_and_gq)
          stan_names_inc_tp_and_gq
          ##
          n_params  <- length(stan_names_main)
          n_params_inc_tp <- length(stan_names_inc_tp)
          n_params_inc_tp_and_gq <- length(stan_names_inc_tp_and_gq)
          ##
          bs_index <-  1:n_params
          bs_index_inc_tp <-  1:n_params_inc_tp
          bs_index_inc_tp_and_gq <-  1:n_params_inc_tp_and_gq
          ##
          # Mow find names, indexes and N's of tp and gq ONLY
          index_tp <- setdiff(bs_index_inc_tp, bs_index)
          names_tp <- stan_names_inc_tp_and_gq[index_tp]
          n_tp <- length(names_tp) ; n_tp
          ##
          index_gq <- setdiff(bs_index_inc_tp_and_gq, bs_index_inc_tp)
          names_gq <- stan_names_inc_tp_and_gq[index_gq]
          n_gq <- length(names_gq) ; n_gq
          ##
          model_so_file <- BayesMVP_model_prep_obj$init_object$model_so_file ; model_so_file
          json_file_path <- BayesMVP_model_prep_obj$init_object$json_file_path ; json_file_path
          ##
          trace_nuisance <- array(1, dim = c(n_iter, n_chains, 1)) ## dummy array
          trace_main_as_list <- trace_nuisance_as_list <- list()
          ##
          for (kk in 1:n_chains) {
            trace_main_as_list[[kk]] <- t(trace_main[1:n_iter, kk, 1:n_chains])
            trace_nuisance_as_list[[kk]] <- t(trace_nuisance[, kk, ])
          }
          ##
          str(trace_main)
          str(trace_main_as_list)
          ##
          pars_indicies_to_track <- 1:n_par_inc_tp_and_gq
          n_params_full <- n_par_inc_tp_and_gq
          ##
          all_param_outs_trace <-    (BayesMVPtest:::fn_compute_param_constrain_from_trace_parallel(
                                                       unc_params_trace_input_main = trace_main_as_list,
                                                       unc_params_trace_input_nuisance = trace_nuisance_as_list,
                                                       pars_indicies_to_track = pars_indicies_to_track,
                                                       n_params_full = n_params_full,
                                                       n_nuisance = 0,
                                                       n_params_main = n_params_main,
                                                       include_nuisance = FALSE,
                                                       model_so_file = model_so_file,
                                                       json_file_path = json_file_path))
          ##
          draws_full <- array(dim = c(n_iter, n_chains, n_params_full))
          ##
          for (kk in 1:n_chains) {
              draws_full[1:n_iter, kk, 1:n_params_full]  <- t(all_param_outs_trace[[kk]][1:n_params_full, 1:n_iter])
          }
          ##
          dimnames(draws_full) <- list(
            iteration = as.character(1:n_iter),
            chain = as.character(1:n_chains),
            variable = stan_names_inc_tp_and_gq
          )
          str(draws_full)
          ##
          # length(bs_names_inc_tp_and_gq)
          # str(bs_names_inc_tp_and_gq)
          ##
          n_cores <- parallel::detectCores()
          ##
          means <- apply(draws_full, FUN = mean, MARGIN = c(3))
          ##
          str(all_param_outs_trace)
          means[10496]
          ##
          tibble_gq <- BayesMVPtest:::generate_summary_tibble(
                                          n_threads = n_cores,
                                          trace = permute_3d_draws_array(draws_full[,,index_gq]),
                                          param_names = names_gq,
                                          n_to_compute = n_gq,
                                          compute_nested_rhat = FALSE,
                                          n_chains = n_chains,
                                          n_superchains = n_chains)
          ## 
          tibble_gq %>% print(n=100)
          extract_params_from_tibble_batch(tibble = tibble_gq, param_strings_vec = c("AUC"), condition = "containing")
          ##
          return(list(tibble_gq = tibble_gq,
                      trace_main = trace_main,
                      trace_tp = trace_tp,
                      trace_gq = trace_gq))



}








#' Fixed create_summary_tibble function to mimic NMA Stan model's k,t ordering
#' Fixed create_summary_tibble to match Stan's output exactly
#' Fixed create_summary_tibble to match Stan's array[n_tests] vector[n_thr_max] structure
#' @export
create_summary_tibble_for_Se_Sp_baseline <- function( arr, 
                                                      n_thr_vec, 
                                                      var_prefix = "Se_baseline"
) {
        
        require(dplyr)
        require(tidyr)
        
        # Get dimensions
        n_iter <- dim(arr)[1]
        n_chains <- dim(arr)[2]
        n_tests <- dim(arr)[3]
        n_thr_max <- dim(arr)[4]
        
        # Initialize list to store results
        summary_list <- list()
        
        # Process in Stan's actual linearization order: k outer, t inner
        for (k in 1:n_thr_max) {        # threshold outer loop
          for (t in 1:n_tests) {        # test inner loop
            # Extract values across all iterations and chains
            vals <- as.vector(arr[, , t, k])
            
            # Check if all values are -1
            if (all(vals < -0.5)) {
              # All -1, this threshold doesn't exist for this test
              summary_row <- data.frame(
                variable = paste0(var_prefix, "[", t, ",", k, "]"),
                mean = -1,
                `2.5%` = -1,
                `50%` = -1,
                `97.5%` = -1,
                check.names = FALSE
              )
            } else {
              # Valid threshold with real data
              vals <- vals[vals > -0.5]  # Remove any -1s (shouldn't be any)
              
              summary_row <- data.frame(
                variable = paste0(var_prefix, "[", t, ",", k, "]"),
                mean = mean(vals),
                `2.5%` = quantile(vals, 0.025, names = FALSE),
                `50%` = quantile(vals, 0.5, names = FALSE),
                `97.5%` = quantile(vals, 0.975, names = FALSE),
                check.names = FALSE
              )
            }
            
            summary_list[[length(summary_list) + 1]] <- summary_row
          }
        }
        
        # Combine all rows into a tibble
        summary_tibble <- bind_rows(summary_list)
        
        tibble(summary_tibble)
  
}






#' Compute Se/Sp baseline for custom covariate values
#' @export
R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline <- function(debugging = FALSE,
                                        model_prep_obj,
                                        model_summary_and_trace_obj,
                                        baseline_case_nd,  # List of vectors, one per test
                                        baseline_case_d,   # List of vectors, one per test
                                        use_probit_link = TRUE,
                                        test_names = NULL
) {
  
          # Extract traces
          trace_main <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_main
          trace_main_dims <- dim(trace_main)
          
          trace_gq <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_gq
          trace_gq_dims <- dim(trace_gq)
          
          stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
          ##
          stan_model_file_path     <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_model_file_path
          stan_functions_directory <- model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_functions_directory
          ##
          stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
          # str(stan_mod_samples)
          ##
          trace_main <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_main
          trace_main_dims <- dim(trace_main) ; trace_main_dims
          ##
          trace_tp <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_tp
          trace_tp_dims <- dim(trace_tp) ; trace_tp_dims
          ##
          trace_gq <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_gq
          trace_gq_dims <- dim(trace_gq) ; trace_gq_dims
          ##
          n_params_main <- trace_main_dims[3] ; n_params_main
          n_chains <- trace_main_dims[2] ; n_chains
          n_iter <- trace_main_dims[1] ; n_iter
          ##
          init_lists_per_chain <- model_prep_obj$init_lists_per_chain
          init_lists_per_chain <- resize_init_list(init_lists_per_chain, n_chains)
          ##
          Stan_data_list <- stan_data_list
          ##
          ## ---- Update cov data:
          ##
          Stan_data_list$baseline_case_nd <- baseline_case_nd
          Stan_data_list$baseline_case_d <- baseline_case_d
          ##
          # Get model data
          n_index_tests <- stan_data_list$n_index_tests ; n_index_tests
          n_thr <- stan_data_list$n_thr ; n_thr
          n_thr_max <- stan_data_list$n_thr_max ; n_thr_max
          n_covariates_nd <- stan_data_list$n_covariates_nd ; n_covariates_nd
          n_covariates_d <- stan_data_list$n_covariates_d ; n_covariates_d
          
          if (debugging) {
            cat("n_index_tests:", n_index_tests, "\n")
            cat("n_thr:", n_thr, "\n")
            cat("n_thr_max:", n_thr_max, "\n")
          }
          
          # Extract beta_mu
          trace_beta_mu <- MetaOrdDTA:::extract_params_from_array_batch(
            debugging = FALSE,
            array_3d = trace_main,
            param_strings_vec = c("beta_mu"),
            condition = "exact_match"
          )
          ##
          trace_beta_tau <- MetaOrdDTA:::extract_params_from_array_batch(
            debugging = FALSE,
            array_3d = trace_tp,
            param_strings_vec = c("beta_tau"),
            condition = "exact_match"
          )
          ##
          trace_beta_L_Sigma <- MetaOrdDTA:::extract_params_from_array_batch(
            debugging = FALSE,
            array_3d = trace_tp,
            param_strings_vec = c("beta_L_Sigma"),
            condition = "exact_match"
          )
          
          # Extract C_vec (or C_array/C_mu depending on model)
          # First check which type of C we have
          gq_param_names <- dimnames(trace_gq)[[3]]
          
          if (any(grepl("^C_vec\\[", gq_param_names))) {
            C_param_name <- "C_vec"
          } else if (any(grepl("^C_array\\[", gq_param_names))) {
            C_param_name <- "C_array"
          } else if (any(grepl("^C_mu\\[", gq_param_names))) {
            C_param_name <- "C_mu"
          } else {
            stop("Could not find C parameters (C_vec, C_array, or C_mu) in generated quantities")
          }
          
          trace_C <- MetaOrdDTA:::extract_params_from_array_batch(
            debugging = FALSE,
            array_3d = trace_gq,
            param_strings_vec = c(C_param_name),
            condition = "exact_match"
          )
          ##
          # Check for C_pred first (random effects model)
          use_C_pred <- FALSE
          C_param_name_for_pred <- NULL
          
          if (any(grepl("^C_pred\\[", gq_param_names))) {
            C_param_name_for_pred <- "C_pred"
            use_C_pred <- TRUE
          } else {
            # Fall back to regular C parameters
            if (any(grepl("^C_vec\\[", gq_param_names))) {
              C_param_name_for_pred <- "C_vec"
            } else if (any(grepl("^C_array\\[", gq_param_names))) {
              C_param_name_for_pred <- "C_array"
            } else if (any(grepl("^C_mu\\[", gq_param_names))) {
              C_param_name_for_pred <- "C_mu"
            } else {
              stop("Could not find C parameters in generated quantities")
            }
          }
          ##
          if (debugging) {
            cat("Using C parameter:", C_param_name_for_pred, "\n")
            cat("use_C_pred:", use_C_pred, "\n")
          }
          ##
          # Extract the C trace
          trace_C_for_pred <- MetaOrdDTA:::extract_params_from_array_batch(
            debugging = FALSE,
            array_3d = trace_gq,
            param_strings_vec = c(C_param_name_for_pred),
            condition = "exact_match"
          )
          
          # Get dimensions
          n_iter <- trace_main_dims[1] ; n_iter
          n_chains <- trace_main_dims[2] ; n_chains
          n_samples <- n_iter * n_chains ; n_samples
          
          # Get parameter names for indexing
          beta_mu_names <- dimnames(trace_beta_mu)[[3]]
          C_names <- dimnames(trace_C)[[3]]
          
          # Flatten the arrays for Rcpp
          trace_beta_mu_flat <- as.vector(trace_beta_mu)
          trace_C_flat <- as.vector(trace_C)
          ##
          trace_beta_tau_flat <- as.vector(trace_beta_tau)
          trace_beta_L_Sigma_flat <- as.vector(trace_beta_L_Sigma)
          trace_C_for_pred_flat <- as.vector(trace_C_for_pred)
          
          # Call Rcpp function
          rcpp_results <- Rcpp_compute_Se_Sp_baseline_loop(
            trace_beta_mu_flat = trace_beta_mu_flat,
            beta_mu_dims = dim(trace_beta_mu),
            trace_C_flat = trace_C_flat,
            C_dims = dim(trace_C),
            beta_mu_names = dimnames(trace_beta_mu)[[3]],
            C_names = dimnames(trace_C)[[3]],
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            n_covariates_nd = n_covariates_nd,
            n_covariates_d = n_covariates_d,
            n_thr = n_thr,
            n_index_tests = n_index_tests,
            n_thr_max = n_thr_max,
            C_param_name = C_param_name,
            use_probit_link = use_probit_link
          )
          ##
          # Extract results with proper dimensions
          Se_baseline <- array(rcpp_results$Se_baseline, 
                               dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          Sp_baseline <- array(rcpp_results$Sp_baseline,
                               dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          Fp_baseline <- array(rcpp_results$Fp_baseline,
                               dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          ##
          # Call Rcpp function
          pred_results <- Rcpp_compute_Se_Sp_baseline_pred_separated(
            trace_beta_mu_flat = trace_beta_mu_flat,
            beta_mu_dims = dim(trace_beta_mu),
            beta_mu_names = dimnames(trace_beta_mu)[[3]],
            trace_beta_tau_flat = trace_beta_tau_flat,
            beta_tau_dims = dim(trace_beta_tau),
            beta_tau_names = dimnames(trace_beta_tau)[[3]],
            trace_beta_L_Sigma_flat = trace_beta_L_Sigma_flat,
            beta_L_Sigma_dims = dim(trace_beta_L_Sigma),
            beta_L_Sigma_names = dimnames(trace_beta_L_Sigma)[[3]],
            trace_C_flat = trace_C_flat,
            C_dims = dim(trace_C_for_pred),
            C_names = dimnames(trace_C_for_pred)[[3]],
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            n_covariates_nd = n_covariates_nd,
            n_covariates_d = n_covariates_d,
            n_thr = n_thr,
            n_index_tests = n_index_tests,
            n_thr_max = n_thr_max,
            C_param_name = C_param_name_for_pred,
            use_probit_link = use_probit_link
          )
          ##
          # Extract results with proper dimensions
          Se_baseline_pred <- array(pred_results$Se_baseline_pred, 
                                    dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          Sp_baseline_pred <- array(pred_results$Sp_baseline_pred,
                                    dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          Fp_baseline_pred <- array(pred_results$Fp_baseline_pred,
                                    dim = c(n_iter, n_chains, n_index_tests, n_thr_max))
          
          # sum(pred_results$Se_baseline_pred == -1)
          
          # Add dimension names
          dimnames(Se_baseline) <- list(
            iteration = as.character(1:n_iter),
            chain = as.character(1:n_chains),
            test = if (!is.null(test_names)) test_names else paste0("Test", 1:n_index_tests),
            threshold = paste0("Thr", 1:n_thr_max)
          )
          
          dimnames(Sp_baseline) <- dimnames(Se_baseline)
          dimnames(Fp_baseline) <- dimnames(Se_baseline)
          ##
          dimnames(Se_baseline_pred) <- dimnames(Se_baseline)
          dimnames(Sp_baseline_pred) <- dimnames(Se_baseline)
          dimnames(Fp_baseline_pred) <- dimnames(Se_baseline)
          ##
          ## ---- Compute AUC:
          ##
          # Add this to the end of your R_fn_compute_Se_Sp_baseline function, before the return statement:
          
          # Compute AUC
          AUC_array <- Rcpp_compute_AUC_all_samples(
            Se_array = as.vector(Se_baseline),
            Sp_array = as.vector(Sp_baseline),
            array_dims = c(n_iter, n_chains, n_index_tests, n_thr_max),
            n_thr = n_thr,
            missing_value_marker = -1.0
          )
          
          # Reshape to proper dimensions
          AUC_baseline <- array(AUC_array, dim = c(n_iter, n_chains, n_index_tests))
          
          # Compute AUC_pred using predictive values
          AUC_pred_array <- Rcpp_compute_AUC_all_samples(
            Se_array = as.vector(Se_baseline_pred),
            Sp_array = as.vector(Sp_baseline_pred),
            array_dims = c(n_iter, n_chains, n_index_tests, n_thr_max),
            n_thr = n_thr,
            missing_value_marker = -1.0
          )
          
          AUC_baseline_pred <- array(AUC_pred_array, dim = c(n_iter, n_chains, n_index_tests))
          
          # Compute AUC differences
          n_test_pairs <- n_index_tests * (n_index_tests - 1) / 2
          AUC_diff <- array(NA, dim = c(n_iter, n_chains, n_test_pairs))
          
          pair_idx <- 1
          for (t1 in 1:(n_index_tests - 1)) {
            for (t2 in (t1 + 1):n_index_tests) {
              AUC_diff[, , pair_idx] <- AUC_baseline[, , t1] - AUC_baseline[, , t2]
              pair_idx <- pair_idx + 1
            }
          }
          
          # Create summaries for AUC
          create_AUC_summary_tibble <- function(auc_array,
                                                test_names = NULL) {
            
                n_iter <- dim(auc_array)[1]
                n_chains <- dim(auc_array)[2]
                n_tests <- dim(auc_array)[3]
                
                summary_list <- list()
                
                for (t in 1:n_tests) {
                  vals <- as.vector(auc_array[, , t])
                  
                  summary_row <- data.frame(
                    variable = paste0("AUC[", t, "]"),
                    test = if (!is.null(test_names)) test_names[t] else paste0("Test", t),
                    mean = mean(vals),
                    `2.5%` = quantile(vals, 0.025, names = FALSE),
                    `50%` = quantile(vals, 0.5, names = FALSE),
                    `97.5%` = quantile(vals, 0.975, names = FALSE),
                    check.names = FALSE
                  )
                  
                  summary_list[[t]] <- summary_row
                }
                
                tibble(bind_rows(summary_list))
                
          }
          
          
          # Return results
          return(list(
            Se_baseline = Se_baseline,
            Sp_baseline = Sp_baseline,
            Fp_baseline = Fp_baseline,
            ##
            Se_baseline_pred = Se_baseline_pred,
            Sp_baseline_pred = Sp_baseline_pred,
            Fp_baseline_pred = Fp_baseline_pred,
            ##
            Se_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Se_baseline, 
                                                 n_thr_vec = n_thr, 
                                                 var_prefix = "Se_baseline"),
            Sp_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Sp_baseline, 
                                                 n_thr_vec = n_thr, 
                                                 var_prefix = "Sp_baseline"),
            Fp_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Fp_baseline, 
                                                 n_thr_vec = n_thr, 
                                                 var_prefix = "Fp_baseline"),
            ##
            ## ---- Prediction intervals:
            ##
            Se_pred_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Se_baseline_pred,
                                                 n_thr_vec = n_thr,
                                                 var_prefix = "Se_baseline_pred"),
            Sp_pred_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Sp_baseline_pred,
                                                 n_thr_vec = n_thr,
                                                 var_prefix = "Sp_baseline_pred"),
            Fp_pred_summaries = create_summary_tibble_for_Se_Sp_baseline(arr = Fp_baseline_pred,
                                                 n_thr_vec = n_thr,
                                                 var_prefix = "Fp_baseline_pred"),
            # ##
            # Xbeta_baseline_nd = Xbeta_baseline_nd,  
            # Xbeta_baseline_d = Xbeta_baseline_d
            ##
            ## ---- AUC stuff:
            ##
            AUC_baseline = AUC_baseline,
            AUC_baseline_pred = AUC_baseline_pred,
            AUC_diff = AUC_diff,
            ##
            AUC_summary = create_AUC_summary_tibble(AUC_baseline, test_names),
            AUC_pred_summary = create_AUC_summary_tibble(AUC_baseline_pred, test_names)#,
            # AUC_diff_summary = create_AUC_summary_tibble(AUC_diff, test_names)
          ))
  
}








#' Update trace_gq with new Se/Sp baseline values
#' @keywords internal
update_trace_gq_with_new_baseline <- function(trace_gq, 
                                              Se_baseline_new, 
                                              Sp_baseline_new,
                                              n_thr) {
  
          # Get dimensions
          n_iter <- dim(Se_baseline_new)[1]
          n_chains <- dim(Se_baseline_new)[2]
          n_tests <- dim(Se_baseline_new)[3]
          n_thr_max <- dim(Se_baseline_new)[4]
          
          # Get parameter names from trace_gq
          param_names <- dimnames(trace_gq)[[3]]
          
          # Create a copy of trace_gq to modify
          trace_gq_updated <- trace_gq
          
          # Update Se_baseline values
          for (t in 1:n_tests) {
            for (k in 1:n_thr_max) {
              if (k <= n_thr[t]) {
                # Find the index for this parameter
                param_name <- paste0("Se_baseline[", t, ",", k, "]")
                param_idx <- which(param_names == param_name)
                
                if (length(param_idx) > 0) {
                  # Update all iterations and chains for this parameter
                  trace_gq_updated[, , param_idx] <- Se_baseline_new[, , t, k]
                }
              }
            }
          }
          
          # Update Sp_baseline values
          for (t in 1:n_tests) {
            for (k in 1:n_thr_max) {
              if (k <= n_thr[t]) {
                # Find the index for this parameter
                param_name <- paste0("Sp_baseline[", t, ",", k, "]")
                param_idx <- which(param_names == param_name)
                
                if (length(param_idx) > 0) {
                  # Update all iterations and chains for this parameter
                  trace_gq_updated[, , param_idx] <- Sp_baseline_new[, , t, k]
                }
              }
            }
          }
          
          return(trace_gq_updated)
  
}











#' Wrapper function that computes everything for new baseline
#' @export
R_fn_compute_MR_complete_new_baseline_analysis_v1 <- function( debugging = FALSE,
                                                            model_prep_obj,
                                                            model_summary_and_trace_obj,
                                                            baseline_case_nd,
                                                            baseline_case_d,
                                                            use_probit_link = TRUE,
                                                            test_names = NULL,
                                                            compute_comparisons = TRUE
) {
  
          # Step 1: Compute Se/Sp/AUC for new baseline
          baseline_results <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(
            debugging = debugging,
            model_prep_obj = model_prep_obj,
            model_summary_and_trace_obj = model_summary_and_trace_obj,
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            use_probit_link = use_probit_link,
            test_names = test_names
          )
          
          stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
          
          if (compute_comparisons) {
            # Step 2: Get the original trace_gq
            trace_gq <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_gq
            
            # Step 3: Update trace_gq with new Se/Sp values
            trace_gq_updated <- update_trace_gq_with_new_baseline(
              trace_gq = trace_gq,
              Se_baseline_new = baseline_results$Se_baseline,
              Sp_baseline_new = baseline_results$Sp_baseline,
              n_thr = stan_data_list$n_thr
            )
            
            # Step 4: Compute comparisons using existing function
            comparison_results <- R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(
              trace_gq = trace_gq_updated,
              n_thr = stan_data_list$n_thr,
              test_names = test_names
            )
            
            # Add comparison results to the output
            baseline_results$comparison_results <- comparison_results
          }
          
          return(baseline_results)
  
}



 






#' Compute NMA comparisons directly from Se/Sp arrays
#' @export
R_fn_compute_NMA_comparisons_from_arrays <- function(Se_baseline, 
                                                     Sp_baseline,
                                                     n_thr,
                                                     test_names = NULL) {
  
          # Get dimensions
          array_dims <- dim(Se_baseline)
          
          # Flatten arrays for C++
          Se_flat <- as.vector(Se_baseline)
          Sp_flat <- as.vector(Sp_baseline)
          
          # Call C++ function
          result <- Rcpp_compute_NMA_comparisons_from_arrays(
            Se_array = Se_flat,
            Sp_array = Sp_flat,
            array_dims = array_dims,
            n_thr = n_thr
          )
          
          # Add metadata
          comp_metadata <- generate_comparison_metadata(n_thr, test_names)
          
          # Combine with summaries
          final_result <- tibble(cbind(comp_metadata, result$summaries))
          
          return(final_result)
  
}








#' Complete analysis for new baseline including comparisons
#' @export
R_fn_compute_MR_complete_new_baseline_analysis_v2 <- function( debugging = FALSE,
                                                               model_prep_obj,
                                                               model_summary_and_trace_obj,
                                                               baseline_case_nd,
                                                               baseline_case_d,
                                                               use_probit_link = TRUE,
                                                               test_names = NULL,
                                                               compute_comparisons = TRUE
) {
          
          # Step 1: Compute Se/Sp/AUC for new baseline
          baseline_results <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(
            debugging = debugging,
            model_prep_obj = model_prep_obj,
            model_summary_and_trace_obj = model_summary_and_trace_obj,
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            use_probit_link = use_probit_link,
            test_names = test_names
          )
          
          stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
          
          if (compute_comparisons) {
            
                # Step 2: Compute comparisons directly from arrays
                comparison_results <- R_fn_compute_NMA_comparisons_from_arrays(
                  Se_baseline = baseline_results$Se_baseline,
                  Sp_baseline = baseline_results$Sp_baseline,
                  n_thr = stan_data_list$n_thr,
                  test_names = test_names
                )
                
                # For predictive comparisons
                comparison_results_pred <- R_fn_compute_NMA_comparisons_from_arrays(
                  Se_baseline = baseline_results$Se_baseline_pred,
                  Sp_baseline = baseline_results$Sp_baseline_pred,
                  n_thr = stan_data_list$n_thr,
                  test_names = test_names
                )
                
                baseline_results$comparison_results <- comparison_results
                baseline_results$comparison_results_pred <- comparison_results_pred
                
          }
          
          return(baseline_results)
  
}




# R_fn_using_Rcpp_compute_NMA_AUC_custom_baseline <- function( model_obj,
#                                                               baseline_case_nd,  # List of vectors, one per test
#                                                               baseline_case_d,   # List of vectors, one per test
#                                                               test_names = NULL
# ) {
#   
#           # Get trace
#           trace_gq <- model_obj$get_trace_generated_quantities()
#           trace_dims <- dim(trace_gq)
#           trace_vec <- as.vector(trace_gq)
#           
#           # Get parameter info from model
#           param_names <- dimnames(trace_gq)[[3]]
#           
#           # Extract necessary indices
#           n_tests <- length(test_names)
#           
#           # Get C_array indices
#           C_array_indices <- list()
#           for (t in 1:n_tests) {
#             # Non-diseased
#             C_nd_idx <- grep(paste0("^C_array\\[1,", t, ",\\d+\\]$"), param_names)
#             C_array_indices[[2*(t-1) + 1]] <- as.list(C_nd_idx)
#             
#             # Diseased  
#             C_d_idx <- grep(paste0("^C_array\\[2,", t, ",\\d+\\]$"), param_names)
#             C_array_indices[[2*(t-1) + 2]] <- as.list(C_d_idx)
#           }
#           
#           # str(C_array_indices)
#           
#           # Get beta_mu indices
#           beta_mu_indices <- list()
#           for (t in 1:n_tests) {
#             # Non-diseased
#             beta_nd_idx <- grep(paste0("^beta_mu\\[", t, ",1,\\d+\\]$"), param_names)
#             beta_mu_indices[[2*(t-1) + 1]] <- as.list(beta_nd_idx)
#             
#             # Diseased
#             beta_d_idx <- grep(paste0("^beta_mu\\[", t, ",2,\\d+\\]$"), param_names)
#             beta_mu_indices[[2*(t-1) + 2]] <- as.list(beta_d_idx)
#           }
#           
#           # Get n_thr and n_covariates
#           n_thr <- sapply(1:n_tests, function(t) {
#             length(grep(paste0("^C_array\\[1,", t, ",\\d+\\]$"), param_names))
#           })
#           
#           n_covariates_nd <- sapply(baseline_case_nd, length)
#           n_covariates_d <- sapply(baseline_case_d, length)
#           
#           # Call Rcpp function
#           result <- Rcpp_fn_compute_NMA_AUC_custom_baseline(
#             trace_gq_vec = trace_vec,
#             trace_dims = trace_dims,
#             C_array_indices = C_array_indices,
#             baseline_case_nd = baseline_case_nd,
#             baseline_case_d = baseline_case_d,
#             beta_mu_indices = beta_mu_indices,
#             n_thr = n_thr,
#             n_covariates_nd = n_covariates_nd,
#             n_covariates_d = n_covariates_d,
#             use_probit_link = TRUE  # Adjust based on your model
#           )
#           
#           # Format output
#           auc_summary <- data.frame(
#             test = 1:n_tests,
#             test_name = if(!is.null(test_names)) test_names else paste0("Test", 1:n_tests),
#             result$auc_summaries
#           )
#           colnames(auc_summary)[3:7] <- c("AUC_mean", "AUC_sd", "AUC_lower", "AUC_median", "AUC_upper")
#           
#           # Format pairwise differences
#           auc_diff_summary <- data.frame(
#             test1 = result$pair_indices[,1],
#             test2 = result$pair_indices[,2],
#             result$auc_diff_summaries
#           )
#           
#           if (!is.null(test_names)) {
#             auc_diff_summary$test1_name <- test_names[auc_diff_summary$test1]
#             auc_diff_summary$test2_name <- test_names[auc_diff_summary$test2]
#           }
#           
#           colnames(auc_diff_summary)[3:7] <- c("AUC_diff_mean", "AUC_diff_sd", 
#                                                "AUC_diff_lower", "AUC_diff_median", "AUC_diff_upper")
#           
#           return(list(
#             auc = auc_summary,
#             auc_diff = auc_diff_summary
#           ))
#   
# }
# 






 




#' generate_comparison_metadata
#' @keywords internal
#' @export
generate_comparison_metadata <- function(n_thr, 
                                         test_names = NULL) {
  
  n_tests <- length(n_thr)
  
  if (is.null(test_names)) {
    test_names <- paste0("Test", 1:n_tests)
  }
  
  # Generate in EXACT same order as Stan
  metadata <- tibble()
  
  for (t1 in 1:(n_tests-1)) {
    for (t2 in (t1+1):n_tests) {
      for (k1 in 1:n_thr[t1]) {
        for (k2 in 1:n_thr[t2]) {
          metadata <- bind_rows(
            metadata,
            tibble(
              test1 = t1,
              test2 = t2,
              test1_name = test_names[t1],
              test2_name = test_names[t2],
              threshold1 = k1,
              threshold2 = k2,
              comparison = paste0(test_names[t1], "_vs_", test_names[t2]),
              threshold_pair = paste0("thr", k1, "_vs_thr", k2)
            )
          )
        }
      }
    }
  }
  
  return(metadata)
  
}








#' R_fn_using_Rcpp_compute_NMA_comparisons_posthoc
#' @keywords internal
#' @export
R_fn_using_Rcpp_compute_NMA_comparisons_posthoc <- function( trace_gq, 
                                                       n_thr,
                                                       test_names = NULL
) {
  
        require(Rcpp)
        require(dplyr)
        
        # Source the C++ code (or have it in your package src/)
        # sourceCpp("nma_comparisons.cpp")
        
        # Pre-compute indices
        n_tests <- length(n_thr)
        param_names <- dimnames(trace_gq)[[3]]
        
        Se_indices_list <- list()
        Sp_indices_list <- list()
        
        for (t in 1:n_tests) {
          Se_indices_list[[t]] <- list()
          Sp_indices_list[[t]] <- list()
          
          for (k in 1:n_thr[t]) {
            Se_param <- paste0("Se_baseline[", t, ",", k, "]")
            Sp_param <- paste0("Sp_baseline[", t, ",", k, "]")
            
            Se_indices_list[[t]][[k]] <- which(param_names == Se_param)[1]
            Sp_indices_list[[t]][[k]] <- which(param_names == Sp_param)[1]
          }
        }
        
        # Flatten array for C++
        trace_dims <- dim(trace_gq)
        trace_gq_vec <- as.vector(trace_gq)
        
        # Call C++ function
        result <- Rcpp_fn_compute_NMA_comparisons(
          trace_gq_vec = trace_gq_vec,
          trace_dims = trace_dims,
          Se_indices_list = Se_indices_list,
          Sp_indices_list = Sp_indices_list,
          n_thr = n_thr
        )
        
        # Add metadata
        comp_metadata <- generate_comparison_metadata(n_thr, test_names)
        
        # Combine with summaries
        final_result <- tibble(cbind(comp_metadata, result$summaries))
        
        return(final_result)
  
}











#' R_fn_using_Rcpp_compute_NMA_performance_posthoc
#' @keywords internal
#' @export
R_fn_using_Rcpp_compute_NMA_performance_posthoc <- function(trace_gq,
                                                      n_thr,
                                                      test_names = NULL) {
          
          require(Rcpp)
          require(dplyr)
          
          # Source the C++ code (or have it in your package src/)
          # sourceCpp("nma_comparisons.cpp")
          
          # Pre-compute indices
          n_tests <- length(n_thr)
          param_names <- dimnames(trace_gq)[[3]]
          
          Se_indices_list <- list()
          Sp_indices_list <- list()
          
          for (t in 1:n_tests) {
            Se_indices_list[[t]] <- list()
            Sp_indices_list[[t]] <- list()
            
            for (k in 1:n_thr[t]) {
              Se_param <- paste0("Se_baseline[", t, ",", k, "]")
              Sp_param <- paste0("Sp_baseline[", t, ",", k, "]")
              
              Se_indices_list[[t]][[k]] <- which(param_names == Se_param)[1]
              Sp_indices_list[[t]][[k]] <- which(param_names == Sp_param)[1]
            }
          }
          
          # Flatten array for C++
          trace_dims <- dim(trace_gq)
          trace_gq_vec <- as.vector(trace_gq)
          
          # Call C++ function
          result <- Rcpp_fn_compute_NMA_performance(
            trace_gq_vec = trace_gq_vec,
            trace_dims = trace_dims,
            Se_indices_list = Se_indices_list,
            Sp_indices_list = Sp_indices_list,
            n_thr = n_thr
          )
          
          # Add metadata
          perf_metadata <- generate_performance_metadata(n_thr, test_names)
          
          # Combine with summaries
          final_result <- cbind(perf_metadata, result$summaries)
          
          return(final_result)
  
}



 


 





#' R_fn_generate_study_accuracy_metadata
#' @keywords internal
#' @export
R_fn_generate_study_accuracy_metadata <- function(n_studies, 
                                             n_thr, 
                                             test_names = NULL) {
  
            n_tests <- length(n_thr)
            
            if (is.null(test_names)) {
              test_names <- paste0("Test", 1:n_tests)
            }
            
            metadata <- tibble()
            
            for (t in 1:n_tests) {
              for (s in 1:n_studies) {
                for (k in 1:n_thr[t]) {
                  metadata <- bind_rows(
                    metadata,
                    tibble(
                      test = t,
                      test_name = test_names[t],
                      study = s,
                      threshold = k
                    )
                  )
                }
              }
            }
            
            return(metadata)
  
}







#' R_fn_extract_study_accuracy_summary
#' @keywords internal
#' @export
R_fn_extract_study_accuracy_summary <- function(  summary_tibble, 
                                                  n_studies, 
                                                  n_thr, 
                                                  test_names = NULL) {
            
            require(dplyr)
            require(tidyr)
            
            # Generate metadata for mapping
            study_metadata <- R_fn_generate_study_accuracy_metadata(n_studies, 
                                                                    n_thr, 
                                                                    test_names)
            
            # Extract study-specific accuracy measures
            study_accuracy_summary <- summary_tibble %>%
              filter(grepl("^(sp_flat|se_flat)\\[", parameter)) %>%
              mutate(
                metric = gsub("_flat\\[.*", "", parameter),
                index = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter))
              ) %>%
              # Remove duplicates if any
              distinct(metric, index, .keep_all = TRUE) %>%
              # Join with metadata
              left_join(
                study_metadata %>% mutate(index = row_number()),
                by = "index"
              ) %>%
              # Remove any rows where the join failed
              filter(!is.na(test)) %>%
              # Select columns before pivoting
              select(test, test_name, study, threshold, metric, mean, sd, `2.5%`, `50%`, `97.5%`) %>%
              # Pivot wider to get all metrics in columns
              pivot_wider(
                names_from = metric,
                values_from = c(mean, sd, `2.5%`, `50%`, `97.5%`),
                names_glue = "{metric}_{.value}"
              ) %>%
              # Clean up column names
              rename_with(~gsub("2\\.5%", "lower", .x)) %>%
              rename_with(~gsub("97\\.5%", "upper", .x)) %>%
              rename_with(~gsub("50%", "median", .x)) %>%
              # Reorder columns for clarity
              select(test, test_name, study, threshold,
                     starts_with("se_"),
                     starts_with("sp_")
                     # starts_with("fp_")
                     )
            
            return(study_accuracy_summary)
  
  
}
  
  
  



  
#' R_fn_deviance_summary
#' @keywords internal
#' @export
R_fn_deviance_summary <- function( summary_tibble, 
                                   test_names = NULL) {
  
            require(dplyr)
            require(tidyr)
            
            # Extract deviance summaries
            deviance_summary <- summary_tibble %>%
              filter(grepl("^(deviance_nd|deviance_d|deviance)\\[", parameter)) %>%
              mutate(
                # Extract test and study indices
                metric = gsub("\\[.*", "", parameter),
                indices = gsub("^[^\\[]+\\[(.*)\\]", "\\1", parameter),
                test = as.numeric(sub(",.*", "", indices)),
                study = as.numeric(sub(".*,", "", indices))
              ) %>%
              mutate(
                test_name = if(!is.null(test_names)) test_names[test] else paste0("Test", test)
              ) %>%
              select(-parameter, -indices, -n_eff, -Rhat, -n_Rhat) %>%
              pivot_wider(
                names_from = metric,
                values_from = c(mean, sd, `2.5%`, `50%`, `97.5%`),
                names_glue = "{metric}_{.value}"
              ) %>%
              # Clean up column names
              rename_with(~gsub("2\\.5%", "lower", .x)) %>%
              rename_with(~gsub("97\\.5%", "upper", .x)) %>%
              rename_with(~gsub("50%", "median", .x))
            
            return(deviance_summary)
            
}







 







