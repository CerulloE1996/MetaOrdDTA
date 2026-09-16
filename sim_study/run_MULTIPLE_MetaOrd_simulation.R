
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
        
        {
                # Extract from scenario_row
                DGM <- scenario_row$DGM
                test <- scenario_row$index_test_MA
                n_studies <- scenario_row$n_studies
                model_parameterisation <- scenario_row$model_parameterisation
                alpha_pi <- scenario_row$alpha_pi
                # box_cox <- scenario_row$box_cox
                softplus <- FALSE ## scenario_row$softplus
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
        }
  
        
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
          
          # Initialize model
          model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
                                        debugging = FALSE,
                                        ##
                                        # seed = MCMC_seed, ## MCMC seed
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
                priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
          } else if (alpha_pi == "flat") {
                priors$prior_dirichlet_alpha <- list( rep(2.0, n_thr + 1), rep(2.0, n_thr + 1))
          } else if (alpha_pi == "informative") { 
            
                generate_neutral_ordinal_dirichlet_prior <- function( n_cat,
                                                                      peak = 0.50,
                                                                      max_alpha = 5,
                                                                      min_alpha = 1,
                                                                      smoothness = 0.25) {

                        stopifnot(n_cat >= 3, peak > 0, peak < 1, max_alpha >= min_alpha, min_alpha > 0, smoothness > 0)
                        ##
                        positions <- seq(0, 1, length.out = n_cat)
                        ##
                        alpha <- min_alpha + (max_alpha - min_alpha) * exp(-((positions - peak)^2) / (2 * smoothness^2))
                        alpha <- round(alpha, 1)
                        ##
                        list( alpha        = alpha,          ## the shared vector
                              non_diseased = alpha,          ## same object, named for drop-in compatibility with the old function
                              diseased     = alpha,
                              sum_alpha    = sum(alpha),
                              prior_mean   = alpha / sum(alpha) )

                }

                outs <- generate_neutral_ordinal_dirichlet_prior( n_cat = n_thr + 1,
                                                                  peak = 0.50,
                                                                  max_alpha = 3,
                                                                  min_alpha = 1,
                                                                  smoothness = 0.25)
                prior_dirichlet_alpha_nd <- outs$alpha
                ##
                outs <- generate_neutral_ordinal_dirichlet_prior( n_cat = n_thr + 1,
                                                                  peak = 0.50,
                                                                  max_alpha = 3,
                                                                  min_alpha = 1,
                                                                  smoothness = 0.25)
                prior_dirichlet_alpha_d <- outs$alpha
                ##
                priors$prior_dirichlet_alpha <- list( prior_dirichlet_alpha_nd, prior_dirichlet_alpha_d)
                priors$prior_dirichlet_alpha
                
          }
          ##
           priors$prior_dirichlet_alpha 
          
          # if (model_parameterisation == "Xu") {
          # 
          #       # priors$prior_dirichlet_alpha <- list(prior_dirichlet_alpha,
          #       #                                      prior_dirichlet_alpha)
          #       priors$prior_dirichlet_alpha <- generate_dual_ordinal_dirichlet_priors(n_cat = n_thr + 1)
          # 
          # } else {
          # 
          #        # priors$prior_dirichlet_alpha <- prior_dirichlet_alpha
          #        priors$prior_dirichlet_alpha <- generate_ordinal_dirichlet_prior(n_cat = n_thr + 1)
          # 
          # }
          
          if (alpha_pi == "one") {
            
              if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                
                      priors$prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
                      priors$prior_kappa_mean <- log(100)
                      priors$prior_kappa_SD <- 1.0
                      priors$kappa_lb <- 0.001
                
              } else { 
                      
                      prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
                      priors$prior_dirichlet_alpha <-  list(prior_dirichlet_alpha, prior_dirichlet_alpha)
                      priors$prior_kappa_mean   <-  c(log(100), log(100))
                      priors$prior_kappa_SD   <-  c(1.0, 1.0)
                      priors$kappa_lb <-  0.001
                
              }
            
          }
          ##
          priors$use_probit_link <- 1
          priors$K_fold_CV_indicator <- rep(1, n_studies_actual)
          ##
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
          if (n_studies == 10) { 
              
              if (model_parameterisation == "Bivariate") { 
                       n_iter_for_4_chains <- 10000
              } else if (model_parameterisation == "Jones") {
                       n_iter_for_4_chains <- 5000
              } else {
                  if (random_thresholds == TRUE) {
                       n_iter_for_4_chains <- 5000
                  } else { 
                       n_iter_for_4_chains <- 5000
                  }
              }
            
          } else { 
            
              if (model_parameterisation == "Bivariate") { 
                     n_iter_for_4_chains <- 10000
              } else if (model_parameterisation == "Jones") {
                     n_iter_for_4_chains <- 5000
              } else {
                if (random_thresholds == TRUE) {
                     n_iter_for_4_chains <- 5000
                } else { 
                     n_iter_for_4_chains <- 5000
                }
              }
            
          }
          ##
          # n_iter_for_4_chains <- 2500 ## --------------------------------------
          ##
          n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
          n_burnin <- 500
          ##
          require(cmdstanr)
          
          # ## n_iter <- 5000
          # n_iter <- 1000
          ##
          # priors$prior_beta_SD_SD <- c(0.5, 0.5)
          
          init_lists_per_chain <- model_prep_obj$init_lists_per_chain
          ##
          # n_chains <- 4 ## -----------------------------
          ##
          init_lists_per_chain <- resize_init_list(  init_lists_per_chain = init_lists_per_chain, 
                                                     n_chains_new = n_chains)
          ##
          if (random_thresholds == TRUE) {
              if (model_parameterisation == "Xu") { 
                  priors$prior_beta_SD_SD <- c(0.5, 0.5)
                  priors$prior_kappa_mean <- c(log(100), log(100))
                  priors$prior_kappa_SD <- c(1.0, 1.0)
                  priors$prior_kappa_dof <- c(5, 3) ## default = c(5, 5)#
                  ##
                  for (kk in 1:n_chains) {
                    init_lists_per_chain[[kk]]$kappa <- c(2000, 500)
                  }
              } else { 
                  priors$prior_beta_SD_SD <- c(0.5, 0.5)
                  priors$prior_kappa_mean <- log(100)
                  priors$prior_kappa_SD <- 1.0
                  priors$prior_kappa_dof <- 5
                  ##
                  for (kk in 1:n_chains) {
                    init_lists_per_chain[[kk]]$kappa <- 1000
                  }
              }
          }  
          ##
          if (model_parameterisation == "Xu") { 
                priors$prior_beta_SD_SD <- c(0.5, 0.5)
                ##
                if (alpha_pi == "informative") {
                  priors$prior_beta_mu_mean[1, 1] <- -1.0 ## ## 1.25
                  priors$prior_beta_mu_mean[2, 1] <- +1.0 ## ## 1.25
                  priors$prior_beta_mu_SD[1, 1] <- 0.5
                  priors$prior_beta_mu_SD[2, 1] <- 0.5
                }
          } else if (model_parameterisation == "R&G") { 
                priors$prior_beta_SD_SD <- 0.5
                priors$prior_raw_scale_SD_SD <- 0.5
          } else if (model_parameterisation == "Jones") { 
                priors$prior_beta_SD_SD <- c(0.5, 0.5)
                priors$prior_raw_scale_SD_SD <- c(0.5, 0.5)
          }
         
          
          ##
          # tibbles_list_testing <- list()
          # ess_vec <- c()
          ##
          options(scipen = 999999)
          options(max.print = 1000) 
          ##
          priors$prior_dirichlet_alpha 
          priors$prior_beta_mu_mean
          priors$prior_beta_mu_SD
          # for (MCMC_seed in 1:5) {
      
            model_samples_obj <- model_prep_obj$sample( seed = MCMC_seed, ## MCMC seed
                                                        n_burnin = n_burnin,
                                                        n_iter = n_iter,
                                                        adapt_delta = 0.65,
                                                        max_treedepth = 10,
                                                        metric_shape = "diag_e",
                                                        priors = priors,
                                                        init_lists_per_chain = init_lists_per_chain,
                                                        n_chains = n_chains)
            ## ---- Get summaries:
            model_summary_obj <- model_samples_obj$summary( compute_main_params = TRUE,
                                                            compute_transformed_parameters = FALSE,
                                                            compute_generated_quantities = TRUE,
                                                            save_log_lik_trace = FALSE,
                                                            use_BayesMVP_for_faster_summaries = TRUE)
            ##
            tibble_main <- model_summary_obj$get_summary_main()
            tibble_gq   <- model_summary_obj$get_summary_generated_quantities()
            tibble_all  <- rbind(tibble_main, tibble_gq)
            ##
            # tibbles_list_testing[[MCMC_seed]] <- tibble_gq
      
              
              try({  
                 model_summary_obj$extract_params(params = "kappa")
                 model_summary_obj$extract_params(params = "C_raw_vec")
                 model_summary_obj$extract_params(params = "beta_mu")
                 model_summary_obj$extract_params(params = "beta_SD")
              })
            
            model_summary_obj$get_efficiency_metrics()$time_burnin
            model_summary_obj$get_efficiency_metrics()$time_sampling
            
              
              # Extract Se/Sp
              tibble_Se <- dplyr::filter(tibble_gq, stringr::str_detect(parameter, "Se_baseline\\["))
              tibble_Sp <- dplyr::filter(tibble_gq, stringr::str_detect(parameter, "Sp_baseline\\["))
              
              # tibble_Sp %>% print(n = 10000)
              
              Se_mean_vec <- Se_median_vec <- Se_lower_vec <- Se_upper_vec <- c()
              Sp_mean_vec <- Sp_median_vec <- Sp_lower_vec <- Sp_upper_vec <- c()
              
              counter <- 1
              for (k in 1:n_thr) {
                Se_median_vec[k] <- 100 * tibble_Se$`50%`[counter]
                Sp_median_vec[k] <- 100 * tibble_Sp$`50%`[counter]
                Se_mean_vec[k]  <- 100 * tibble_Se$mean[counter]
                Sp_mean_vec[k]  <- 100 * tibble_Sp$mean[counter]
                Se_lower_vec[k] <- 100 * tibble_Se$`2.5%`[counter]
                Sp_lower_vec[k] <- 100 * tibble_Sp$`2.5%`[counter]
                Se_upper_vec[k] <- 100 * tibble_Se$`97.5%`[counter]
                Sp_upper_vec[k] <- 100 * tibble_Sp$`97.5%`[counter]
                counter <- counter + 1
              }
              
              # Check ESS (filter out using "vec_index_inner" though!!):
              ESS_vec_Se <- tibble_Se$n_eff
              ESS_vec_Sp <- tibble_Sp$n_eff
              ##
              ESS_vec_Se_main_thr <- ESS_vec_Se[vec_index_inner]
              ESS_vec_Sp_main_thr <- ESS_vec_Sp[vec_index_inner]
              ##
              ESS_vec <- c(ESS_vec_Se_main_thr, ESS_vec_Sp_main_thr) ; ESS_vec
              min_ESS <- min(ESS_vec, na.rm = TRUE) ; min_ESS
              min_ESS_threshold <- ifelse(n_studies == 10,
                                          5000, 
                                          1000)
              ##
              # ess_vec[MCMC_seed] <- min_ESS
              
          # }
          
         ## havent run w/ " sd_thr_jitter <- 0.001"
          
 
          # 
          # min(ess_vec) ; mean(ess_vec) ; max(ess_vec)
          # 
          # ess_thr <- 1000
          # 
          # which(ess_vec > ess_thr)
          # 
          # 
          # library(dplyr)
          # library(tidyr)
          # library(stringr)
          # 
          # Se_wide <- bind_rows(
          #   lapply(which(ess_vec > ess_thr), function(seed) {
          #     filter(tibbles_list_testing[[seed]], str_detect(parameter, "Se_baseline\\[")) %>%
          #       select(parameter, mean) %>%
          #       mutate(seed = seed)
          #   })
          # ) %>%
          #   pivot_wider(names_from = seed, values_from = mean, names_prefix = "seed_") %>%
          #   mutate(
          #     min_Se = apply(select(., starts_with("seed_")), 1, min, na.rm = TRUE),
          #     max_Se = apply(select(., starts_with("seed_")), 1, max, na.rm = TRUE),
          #     range_Se = max_Se - min_Se
          #   )
          # Se_wide
          # Se_wide <- Se_wide %>%
          #   mutate_if(is.numeric, ~. * 100)
          # 
          # 
          # Sp_wide <- bind_rows(
          #   lapply(which(ess_vec > ess_thr), function(seed) {
          #     filter(tibbles_list_testing[[seed]], str_detect(parameter, "Sp_baseline\\[")) %>%
          #       select(parameter, mean) %>%
          #       mutate(seed = seed)
          #   })
          # ) %>%
          #   pivot_wider(names_from = seed, values_from = mean, names_prefix = "seed_") %>%
          #   mutate(
          #     min_Se = apply(select(., starts_with("seed_")), 1, min, na.rm = TRUE),
          #     max_Se = apply(select(., starts_with("seed_")), 1, max, na.rm = TRUE),
          #     range_Se = max_Se - min_Se
          #   )
          # Sp_wide <- Sp_wide %>%
          #   mutate_if(is.numeric, ~. * 100)
          # 
          # Se_and_Sp <- rbind(Se_wide, Sp_wide)
          # 
          # round(mean(Se_and_Sp$range_Se), 3)
          # round(max(Se_and_Sp$range_Se), 3)
      
          ## ---- Xu model, fixed-cutpoints:
          ## 0.172 ## w/ jitte = 0.01 (or 0.1?!) (8 chains)
          ## 0.114 w/ jitter = 0.001 (8 chains) - min_ESS ~ 13,000
          ## 0.194 w/ jitter = 0.001 (4 chains) - min_ESS ~ 6500
          ## ---- Xu model, RANDOM-cutpoint (all 4 chains):
          ## 0.172 w/ df = 5, log(100), sd = 1 kappa priors
          ## 0.13 w/ df = c(5, 3)
          
          # priors$prior_kappa_dof
          
          
          # 
          # try({  
          #   model_summary_obj$extract_params(params = "kappa")
          # })
          
          
          if (min_ESS > min_ESS_threshold) { # Store results
         
            results <- list(
              seed = seed,
              ##
              x = x,
              ##
              n_studies = n_studies,
              ##
              model_parameterisation = model_parameterisation,
              # box_cox = box_cox,
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
            # Check MCSE periodically
            if (iii %% check_MCSE_interval == 0 || iii == length(sim_ids_chunk)) {
              
                    cat("\n=== Checking MCSE status ===\n")
                    
                    # Use the check function
                    current_mcse_status <- check_MCSE_status( output_dir = output_dir,
                                                              file_name_string = file_name_string,
                                                              max_seed_to_check = max_seed_to_check,
                                                              ##
                                                              min_runs_needed = min_runs_needed,
                                                              target_MCSE_RMSE = target_MCSE_RMSE,
                                                              ##
                                                              test = test,
                                                              DGM = DGM,
                                                              n_studies = n_studies,
                                                              ##
                                                              min_ESS_thr = NULL,
                                                              vec_index_inner = vec_index_inner,
                                                              min_studies_per_threshold = NULL)
                                                            
                                                         
                    
                    cat(sprintf("Completed sims: %d, Current MCSE: %.4f%%\n", 
                                current_mcse_status$n_sims, 
                                current_mcse_status$mcse_value * 100))
                    
                    # Check if MCSE is valid AND target is met
                    if (is.finite(current_mcse_status$mcse_value) && 
                        current_mcse_status$mcse_met) {
                      cat("✓ TARGET MCSE ACHIEVED! Stopping early.\n")
                      if (stop_early_if_target_met) {
                        return(list(stopped_early = TRUE, 
                                    final_n_sims = current_mcse_status$n_sims))
                      }
                    }
              
            }
          } else {
            cat(sprintf("Skipping - ESS too low: %d\n", min_ESS))
          }
        }
        
        return(list(completed = TRUE))
        
}





