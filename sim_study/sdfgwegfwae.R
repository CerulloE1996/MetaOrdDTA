
# ============================================================================================================
# Run Multiple MetaOrdDTA Simulations
# ============================================================================================================


run_multiple_sims <- function(MCMC_seed,
                              sim_ids_chunk,
                              n_chains,
                              scenario_row,
                              sim_study_dir,
                              output_dir,
                              local_pkg_dir,
                              file_name_string,
                              max_seed_to_check,
                              min_runs_needed = 30,
                              target_MCSE_RMSE,
                              check_MCSE_interval,
                              stop_early_if_target_met = TRUE
) {
  
        require(MetaOrdDTA)
        require(cmdstanr)
        require(dplyr)
        
        # Extract from scenario_row
        DGM <- scenario_row$DGM
        test <- scenario_row$index_test_MA
        n_studies <- scenario_row$n_studies
        model_parameterisation <- scenario_row$model_parameterisation
        alpha_pi <- scenario_row$alpha_pi
        # box_cox <- scenario_row$box_cox
        softplus <- scenario_row$softplus
        random_thresholds <- scenario_row$random_thresholds
        Dirichlet_random_effects_type <- scenario_row$Dirichlet_random_effects_type
        
        # Set threshold ranges
        if (test == "GAD_2") {
          vec_index_inner <- 1:6
        } else if (test == "HADS") {
          vec_index_inner <- 3:17
        } else if (test == "BAI") {
          vec_index_inner <- 3:43
        }
        
        missing_thr <- TRUE
        ##
        N_per_study_mean <- scenario_row$N_per_study_mean
        ##
        N_per_study_SD <- N_per_study_mean
        ##
        assume_perfect_GS <- TRUE
        ##
        true_Mean_prev <- 0.20 ; true_SD_probit_prev <- 0.25
        
        iii <- 1
        
  for (iii in 1:length(sim_ids_chunk)) {
          
          sim_index <- sim_ids_chunk[iii]
          seed <- sim_index
          
          cat(sprintf("\n--- Running simulation %d with seed = %d ---\n", sim_index, seed))
          set.seed(seed)
          
          # Source files - these create true_DGM_Se and true_DGM_Sp
          source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA_v2.R"), local = TRUE)
          source(file.path(local_pkg_dir, "inst", "examples", "MA_missing_thr_prep_data.R"), local = TRUE)
          source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"), local = TRUE)
          
          true_vals <- get_true_values(test = test,
                                       DGM = DGM)
          true_DGM_Se <- true_vals$Se
          true_DGM_Sp <- true_vals$Sp
          true_DGM_Fp <- true_vals$Fp
          
          # Select data
          if (test == "GAD_2") {
            if (missing_thr == TRUE) x <- x_GAD2_w_missing_thr
            else x <- x_GAD2
          } else if (test == "HADS") {
            if (missing_thr == TRUE) x <- x_HADS_w_missing_thr
            else x <- x_HADS
          } else if (test == "BAI") {
            if (missing_thr == TRUE) x <- x_BAI_w_missing_thr
            else x <- x_BAI
          }
          
          n_thr <- ncol(x[[1]]) - 1
          n_studies_actual <- nrow(x[[1]])
          
          use_custom_file <- FALSE
          custom_file_name <- NULL
          ##
          if (model_parameterisation == "Bivariate") {
              use_custom_file  <- TRUE
              custom_file_name <- "DTA_MA_strat_bivariate.stan"
          }
          
          
          use_custom_file  <- TRUE
          custom_file_name <- "DTA_MA_Gat_RANDthr_kappa_2.stan"
          
          
         
          # library(cmdstanr)
          # cmdstan_version()  # Check current version
          # rebuild_cmdstan()  # This rebuilds the toolchain
          
          
          
          quick_kappa_prior(log_mean = log(50), log_sd = 1, df = 5)
          
          quick_kappa_prior(log_mean = log(100), log_sd = 1, df = 5)
          
          quick_kappa_prior(log_mean = log(200), log_sd = 1, df = 5)
          
          
          # Initialize model
          model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
                                        debugging = TRUE,
                                        ##
                                        x = x,
                                        indicator_index_test_in_study = NULL,
                                        intercept_only = TRUE,
                                        cov_data = NULL,
                                        ##
                                        n_chains = n_chains,
                                        ##
                                        network = FALSE,
                                        compound_symmetry = NULL,
                                        hetero_sigma = NULL,
                                        ##
                                        prior_only = FALSE,
                                        ##
                                        model_parameterisation = model_parameterisation,
                                        random_thresholds = random_thresholds,
                                        Dirichlet_random_effects_type = Dirichlet_random_effects_type,
                                        ##
                                        custom_file_name = custom_file_name)
          
          # Set priors
          priors <- model_prep_obj$priors
          
          if (alpha_pi == "one") {
            
                prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
            
          } else if (alpha_pi == "flat") {
            
                prior_dirichlet_alpha <- rep(2.0, n_thr + 1)
            
          }
          
          if (model_parameterisation == "Xu") {
            
                priors$prior_dirichlet_alpha <- list(prior_dirichlet_alpha,
                                                     prior_dirichlet_alpha)
            
          } else {
            
                 priors$prior_dirichlet_alpha <- prior_dirichlet_alpha
            
          }
          
          priors$use_probit_link <- 1
          priors$K_fold_CV_indicator <- rep(1, n_studies_actual)
          
          if (model_parameterisation == "R&G") {
              priors$prior_raw_scale_mu_SD <- 1
              priors$prior_raw_scale_SD_SD <- 1
          }
          
          # Set initial values
          if (model_parameterisation == "Bivariate") {
            
                init_lists_per_chain <- list()
                mu <- list()
                sigma <- list()
                for (k in 1:n_thr) {
                  mu[[k]] <- c(-1, 1)
                  sigma[[k]] <- c(0.01, 0.01)
                }
                for (kk in 1:n_chains) {
                  init_lists_per_chain[[kk]] <- list(mu = mu, sigma = sigma)
                }
            
          } else {
            
               init_lists_per_chain <- NULL
               
          }
          
          # Sample model
          if (n_studies == 50) { 
            n_iter_for_4_chains <- ifelse(model_parameterisation == "Bivariate", 2000, 1000)
          } else if (n_studies == 10) { 
            n_iter_for_4_chains <- ifelse(model_parameterisation == "Bivariate", 2000, 1000)
          }
          ##
          n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
          n_burnin <- 500
          
          model_samples_obj <- model_prep_obj$sample(
            n_burnin = n_burnin,
            n_iter = n_iter,
            adapt_delta = 0.65,
            max_treedepth = 10,
            metric_shape = "diag_e",
            priors = priors,
            init_lists_per_chain = init_lists_per_chain,
            n_chains = n_chains
          )
          
          # Get summaries
          model_summary_obj <- model_samples_obj$summary(
            compute_main_params = TRUE,
            compute_transformed_parameters = TRUE,
            compute_generated_quantities = TRUE,
            save_log_lik_trace = FALSE,
            use_BayesMVP_for_faster_summaries = TRUE
          )
          
          # internal_obj <- model_samples_obj$internal_obj
          # MCMC_params <- model_samples_obj$MCMC_params
          # model_results
          # init_object
          # n_nuisance
          
          # model_summary_obj <- create_summary_and_traces(package = "MetaOrdDTA",
          #                           use_bridgestan = FALSE,
          #                           use_BayesMVP_for_faster_summaries = TRUE,
          #                           internal_obj = internal_obj,
          #                           MCMC_params = MCMC_params,
          #                           ##
          #                           model_results  = NULL,
          #                           init_object = NULL,
          #                           n_nuisance = 0,
          #                           ##
          #                           compute_main_params = TRUE,
          #                           compute_transformed_parameters = TRUE,
          #                           compute_generated_quantities = TRUE,
          #                           save_log_lik_trace = FALSE,
          #                           ##
          #                           compute_nested_rhat = FALSE)
          
          
          # package = "MetaOrdDTA"
          # use_bridgestan = FALSE
          # use_BayesMVP_for_faster_summaries = TRUE
          # internal_obj = internal_obj
          # MCMC_params = MCMC_params
          # ##
          # model_results  = NULL
          # init_object = NULL
          # n_nuisance = 0
          # ##
          # compute_main_params = TRUE
          # compute_transformed_parameters = TRUE
          # compute_generated_quantities = TRUE
          # save_log_lik_trace = FALSE
          # ##
          # compute_nested_rhat = FALSE
          
          # model_summary_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
          
          # model_summary_obj$get_summary_transformed()
          # model_summary_obj$get_summary_generated_quantities()
          ##
          tibble_main <- model_summary_obj$get_summary_main()
          tibble_gq   <- model_summary_obj$get_summary_generated_quantities()
          tibble_all  <- rbind(tibble_main, tibble_gq)
          
          # Extract Se/Sp
          tibble_Se <- dplyr::filter(tibble_gq, stringr::str_detect(parameter, "Se_baseline\\["))
          tibble_Sp <- dplyr::filter(tibble_gq, stringr::str_detect(parameter, "Sp_baseline\\["))
          
          Se_mean_vec <- Se_median_vec <- Se_lower_vec <- Se_upper_vec <- c()
          Sp_mean_vec <- Sp_median_vec <- Sp_lower_vec <- Sp_upper_vec <- c()
          
          counter <- 1
          for (k in 1:n_thr) {
            Se_median_vec[k] <- 100 * tibble_Se$`50%`[counter]
            Sp_median_vec[k] <- 100 * tibble_Sp$`50%`[counter]
            Se_mean_vec[k] <- 100 * tibble_Se$mean[counter]
            Sp_mean_vec[k] <- 100 * tibble_Sp$mean[counter]
            Se_lower_vec[k] <- 100 * tibble_Se$`2.5%`[counter]
            Sp_lower_vec[k] <- 100 * tibble_Sp$`2.5%`[counter]
            Se_upper_vec[k] <- 100 * tibble_Se$`97.5%`[counter]
            Sp_upper_vec[k] <- 100 * tibble_Sp$`97.5%`[counter]
            counter <- counter + 1
          }
          
          # Check ESS
          ESS_vec <- c(tibble_Se$n_eff, tibble_Sp$n_eff)
          min_ESS <- min(ESS_vec, na.rm = TRUE) ; min_ESS
          min_ESS_threshold <- ifelse(n_studies == 10,
                                      500, 
                                      100)
          
          if (min_ESS > min_ESS_threshold) { # Store results
         
            results <- list(
              seed = seed,
              ##
              x = x,
              ##
              n_studies = n_studies,
              ##
              model_parameterisation = model_parameterisation,
              box_cox = box_cox,
              softplus = softplus,
              random_thresholds = random_thresholds,
              ##
              tibble_all = tibble_all,
              ##
              true_DGM_Se = true_DGM_Se,
              true_DGM_Sp = true_DGM_Sp,
              ##
              Se_median_vec = Se_median_vec,
              Sp_median_vec = Sp_median_vec,
              Se_mean_vec = Se_mean_vec,
              Sp_mean_vec = Sp_mean_vec,
              Se_lower_vec = Se_lower_vec,
              Sp_lower_vec = Sp_lower_vec,
              Se_upper_vec = Se_upper_vec,
              Sp_upper_vec = Sp_upper_vec
            )
            
       
            
            # Save
            results_file <- file.path(output_dir, 
                                      paste0("seed_", seed, "_", file_name_string, "_results_list.RDS"))
            saveRDS(results, results_file)
            cat(sprintf("File saved: %s\n", basename(results_file)))
            
            # Check MCSE periodically (like your LC-MVP)
            if (iii %% check_MCSE_interval == 0 || iii == length(sim_ids_chunk)) {
              
                    cat("\n=== Checking MCSE status ===\n")
                    
                    # Load all results to check MCSE
                    complete_results_list <- list()
                    for (s in 1:max_seed_to_check) {
                      check_file <- file.path(output_dir, 
                                              paste0("seed_", s, "_", file_name_string, "_results_list.RDS"))
                      if (file.exists(check_file)) {
                        try({
                          complete_results_list[[s]] <- readRDS(check_file)
                        }, silent = TRUE)
                      }
                    }
                    
                    n_loaded <- sum(!sapply(complete_results_list, is.null))
                    
                    if (n_loaded >= min_runs_needed) {
                      # Compute metrics
                      metrics <- compute_simulation_metrics(
                        complete_results_list = complete_results_list,
                        ##
                        true_DGM_Se = true_DGM_Se,
                        true_DGM_Sp = true_DGM_Sp,
                        min_k = vec_index_inner[1],
                        max_k = tail(vec_index_inner, 1)
                      )
                      
                      tables <- create_metrics_summary_table(
                        metrics = metrics,
                        min_k = vec_index_inner[1],
                        max_k = tail(vec_index_inner, 1)
                      )
                      
                      # Get worst MCSE
                      MCSE_RMSE_Se_Sp <- c(mean(tables$by_threshold$MCSE_RMSE_Se),
                                           mean(tables$by_threshold$MCSE_RMSE_Sp))
                      max_MCSE <- max(MCSE_RMSE_Se_Sp, na.rm = TRUE)
                      
                      cat(sprintf("Completed sims: %d, Current MCSE: %.4f%%\n", n_loaded, max_MCSE * 100))
                      
                      if (!is.na(max_MCSE) && max_MCSE <= target_MCSE_RMSE) {
                        cat("✓ TARGET MCSE ACHIEVED! Stopping early.\n")
                        if (stop_early_if_target_met) {
                          return(list(stopped_early = TRUE, final_n_sims = n_loaded))
                        }
                      }
                    }
                    
            }
          } else {
            cat(sprintf("Skipping - ESS too low: %d\n", min_ESS))
          }
        }
        
        return(list(completed = TRUE))
        
}





