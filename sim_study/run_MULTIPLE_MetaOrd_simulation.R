
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
                              stop_early_if_target_met = TRUE,
                              sampling_algorithm,
                              min_ESS_of_10_studies,
                              min_ESS_of_50_studies,
                              ordinal_prior_type
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


        ## ("Dirichlet" = each ordinal model's own Dirichlet prior on the cutpoints, "Se_Sp" = the direct Se / Sp prior)
        stopifnot(ordinal_prior_type %in% c("Dirichlet", "Se_Sp"))
        candidate_fields <- c("custom_file_name", "custom_cpp_user_header_file_path", "prior_fields")
        candidate_route <- candidate_fields %in% names(scenario_row)
        if (any(candidate_route)) {
            ## (the pop-median models carry their Se / Sp prior fields; the centred median model has none)
            is_a_pop_median_model <- grepl("popmed", scenario_row$custom_file_name[[1]], fixed = TRUE)
            is_a_centred_C_model <- grepl("_centred_C_", scenario_row$custom_file_name[[1]], fixed = TRUE)
            stopifnot(all(candidate_route[1:2]), candidate_route[3] == is_a_pop_median_model,
                      model_parameterisation %in% c("ord_bivariate", "R&G"),
                      isTRUE(random_thresholds), Dirichlet_random_effects_type == "kappa",
                      is.character(file_name_string), length(file_name_string) == 1,
                      !is.na(file_name_string),
                      grepl("popmed", file_name_string, fixed = TRUE) == is_a_pop_median_model,
                      grepl("_centred_C", file_name_string, fixed = TRUE) == is_a_centred_C_model)
            stopifnot(scenario_row$custom_file_name[[1]] %in%
                            c("DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan",
                              "DTA_MA_ord_biv_RANDthr_kappa_popmed_centred_C_Se_Sp_prior_NicoStan_fused_AVX.stan",
                              "DTA_MA_ord_biv_RANDthr_kappa_centred_C_NicoStan_fused_AVX.stan",
                              "DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan"))
            stopifnot(grepl("ord_HSROC", scenario_row$custom_file_name[[1]], fixed = TRUE) ==
                            (model_parameterisation == "R&G"))
            stopifnot(is.character(scenario_row$custom_cpp_user_header_file_path[[1]]),
                      length(scenario_row$custom_cpp_user_header_file_path[[1]]) == 1,
                      file.exists(scenario_row$custom_cpp_user_header_file_path[[1]]))
            if (is_a_pop_median_model) {
                  prior_fields <- scenario_row$prior_fields[[1]]
                  stopifnot(is.list(prior_fields), !is.null(names(prior_fields)),
                            !anyDuplicated(names(prior_fields)),
                            setequal(names(prior_fields),
                                     c("use_Se_Sp_prior", "Se_Sp_prior_type", "n_chosen_max", "n_chosen",
                                       "chosen_cutoffs", "prior_Se_Sp_dirichlet_A", "prior_Se_Sp_beta",
                                       "prior_Se_Sp_link_mean_SD")))
            }
        }

  for (iii in 1:length(sim_ids_chunk)) {

          sim_index <- sim_ids_chunk[iii]
          seed <- sim_index

          ## (a seed whose results are saved already, by whichever model, is done and is never fitted again)
          file_of_the_results_of_this_seed <- file.path(output_dir,
                                                        paste0("seed_", seed, "_", file_name_string,
                                                               "_results_list.RDS"))
          if (file.exists(file_of_the_results_of_this_seed)) {
            next
          }

          cat(sprintf("\n--- Running simulation %d with seed = %d ---\n", sim_index, seed))
          set.seed(seed)
          time_of_the_start_of_this_fit <- Sys.time()

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
               custom_file_name <- "DTA_MA_strat_bivariate_NicoStan_fused_AVX.stan"
          }
          ## (fused AVX ord_biv / ord_HSROC models: NicoSIMD header as the user header, run on cmdstanr)
          if (model_parameterisation %in% c("ord_bivariate", "R&G")) {
              use_custom_file  <- TRUE
              name_of_the_model_family <- if (model_parameterisation == "ord_bivariate") "DTA_MA_ord_biv_" else
                                          "DTA_MA_ord_HSROC_"
              custom_file_name <- paste0(name_of_the_model_family,
                                         if (random_thresholds) "RANDthr_kappa" else "FIXEDthr",
                                         "_NicoStan_fused_AVX.stan")
          }
          if (any(candidate_route)) {
              custom_file_name <- scenario_row$custom_file_name[[1]]
              use_custom_file <- TRUE
          }
          CPU_has_AVX512 <- any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
          ## (the compiler flags: AVX-512 on the local-HPC, AVX2 on the laptop)
          AVX_FLAGS_of_this_machine <- if (CPU_has_AVX512) {
              "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
          } else {
              "-march=native -mtune=native -mavx2"
          }
          # NicoSIMD_user_header <- file.path(dirname(dirname(local_pkg_dir)), "R_packages", "NicoStan", "inst",
          #                                   "NicoSIMD", "inst", "include", "NicoSIMD",
          #                                   "stan_external_functions_of_NicoSIMD.hpp")
          ## (the frozen NicoSIMD of the sim, first in .libPaths())
          NicoSIMD_user_header <- file.path(find.package("NicoSIMD"), "include", "NicoSIMD",
                                            "stan_external_functions_of_NicoSIMD.hpp")

          if (any(candidate_route)) {
              NicoSIMD_user_header <- scenario_row$custom_cpp_user_header_file_path[[1]]
          }

          time_after_the_simulated_data <- Sys.time()
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
                                        custom_file_name = custom_file_name,
                                        ##
                                        advanced_compile = use_custom_file,
                                        force_recompile = FALSE,
                                        quiet = TRUE,
                                        set_custom_CXX_CPP_flags = use_custom_file,
                                        custom_cpp_user_header_file_path = if (use_custom_file)
                                                                  NicoSIMD_user_header else NULL,
                                        FMA_FLAGS = "-mfma",
                                        AVX_FLAGS = AVX_FLAGS_of_this_machine)

          # Set priors
          priors <- model_prep_obj$priors

          if (alpha_pi == "one") {
                priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
          } else if (alpha_pi == "flat") {
                priors$prior_dirichlet_alpha <- list( rep(2.0, n_thr + 1), rep(2.0, n_thr + 1))
          } else if (alpha_pi == "informative_from_the_literature") {
                stop("No informative prior from the literature is derived by Stan sampling yet")
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

          # if (model_parameterisation == "ord_bivariate") {
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
                       n_iter_for_4_chains <- 1500 ## 1500 # #5000
                  } else {
                       n_iter_for_4_chains <- 3500 ## 1500 ## 5000
                  }
              }

          } else {

              if (model_parameterisation == "Bivariate") {
                     n_iter_for_4_chains <- 10000
              } else if (model_parameterisation == "Jones") {
                     n_iter_for_4_chains <- 5000
              } else {
                if (random_thresholds == TRUE) {
                     n_iter_for_4_chains <-  1000 ## 5000
                } else {
                     n_iter_for_4_chains <-  2000 ## 1000 ## 5000
                }
              }

          }
          ##
          # n_iter_for_4_chains <- 2500 ## --------------------------------------
          ##
          ## ---- BAI, 10 studies, ord_biv_random_C: 2x the iterations of the saved fits (3000 -> 6000) ------------------------------------
          ##
          ## (~93-95% of fits reach min ESS 1000, from ~88-90%)
          if (test == "BAI" &&
              DGM == "ord_biv_random_C" &&
              n_studies == 10 &&
              N_per_study_mean == 500 &&
              model_parameterisation == "ord_bivariate" &&
              random_thresholds == TRUE &&
              alpha_pi == "one") {

                      # n_iter_for_4_chains <- 3000
                      n_iter_for_4_chains <- 6000

          }
          ##
          ## ---- GAD-2, 10 studies, R&G_random_C, ord_biv: 2x the iterations of the saved fits (1500 -> 3000) ------------------------
          ##
          ## (~98-99% of fits reach min ESS 1000, from ~70-76%)
          if (test == "GAD_2" &&
              DGM == "R&G_random_C" &&
              n_studies == 10 &&
              N_per_study_mean == 500 &&
              model_parameterisation == "ord_bivariate" &&
              random_thresholds == TRUE &&
              alpha_pi == "one") {

                      n_iter_for_4_chains <- 3000

          }
          ##
          ## ---- BAI, 50 studies, ord_biv_random_C: 3.5x the iterations -------------------------------------------------------------------
          ##
          ## (the min ESS of Se / Sp over all 63 cutoffs reaches 200)
          if (test == "BAI" &&
              DGM == "ord_biv_random_C" &&
              n_studies == 50 &&
              N_per_study_mean == 500 &&
              model_parameterisation == "ord_bivariate" &&
              random_thresholds == TRUE &&
              alpha_pi == "one") {

                      n_iter_for_4_chains <- 3500

          }
          ##
          n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
          ##
          ## (NicoStan's burn-in is 125-250 iterations; Stan's is never fewer than 200)
          n_burnin <- if (sampling_algorithm == "NicoStan") 500 else 500
          # n_burnin <- 500
          # n_burnin <- 250
          # n_burnin <- 200
          # n_burnin <- 125
          ##
          require(cmdstanr)
          ##
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
                
                    if (model_parameterisation == "ord_bivariate") {
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
          if (model_parameterisation == "ord_bivariate") {
            
                    priors$prior_beta_SD_SD <- c(0.5, 0.5)
                    ##
                    if (alpha_pi == "informative_from_the_literature") {
                          priors$prior_beta_mu_mean[1, 1] <-  0
                          priors$prior_beta_mu_mean[2, 1] <-  0
                          priors$prior_beta_mu_SD[1, 1] <-  c(HADS = 0.20, GAD_2 = 0.25, BAI = 0.35)[[test]]
                          priors$prior_beta_mu_SD[2, 1] <-  c(HADS = 0.20, GAD_2 = 0.25, BAI = 0.35)[[test]]
                    }
                    if (alpha_pi == "informative") {
                          priors$prior_beta_mu_mean[1, 1] <- -1.0 ## ## 1.25
                          priors$prior_beta_mu_mean[2, 1] <- +1.0 ## ## 1.25
                          priors$prior_beta_mu_SD[1, 1]   <- 0.5
                          priors$prior_beta_mu_SD[2, 1]   <- 0.5
                    }
                
          } else if (model_parameterisation == "R&G") {
                
                    priors$prior_beta_SD_SD      <- 0.5
                    priors$prior_raw_scale_SD_SD <- 0.5
                    ##
                    ## (the RC ord_HSROC median vs pop-median sim study: both models on the priors set above)
                    if (alpha_pi %in% c("one", "informative_from_the_literature") &&
                        !("RC_ord_HSROC_model" %in% names(scenario_row))) {
                      stop(paste0("No ordinal-HSROC prior is derived by Stan sampling yet for alpha_pi = ", alpha_pi))
                    }
                
          } else if (model_parameterisation == "Jones") {
            
                    priors$prior_beta_SD_SD      <- c(0.5, 0.5)
                    priors$prior_raw_scale_SD_SD <- c(0.5, 0.5)
                
          } else if (model_parameterisation == "Bivariate") {
                
                    ##
                    ## (the probit Se and Sp of each threshold: [1] Se, [2] Sp; "uniform_Se_and_Sp" is N(0, 1), i.e.
                    ##  Uniform(0, 1) on Se and Sp at every threshold)
                    if (alpha_pi == "one") {
                      file_of_the_prior <-  file.path(sim_study_dir,
                                                      "prior_derivation_tools_checks",
                                                      paste0("Bivariate_flat_1_prior_by_Stan_sampling_",
                                                             test,
                                                             "_at_the_lower_bound_minus_15.rds"))
                      derived_Bivariate <-  readRDS(file_of_the_prior)$fit_of_each_quantity
                      priors$prior_mu_mean <-  list(derived_Bivariate$Se$mu_mean,
                                                    derived_Bivariate$Sp$mu_mean)
                      priors$prior_mu_SD <-  list(derived_Bivariate$Se$mu_SD,
                                                  derived_Bivariate$Sp$mu_SD)
                    } else if (alpha_pi == "uniform_Se_and_Sp") {
                      priors$prior_mu_mean <-  list(rep(0, n_thr), rep(0, n_thr))
                      priors$prior_mu_SD <-  list(rep(1, n_thr), rep(1, n_thr))
                    } else {
                      stop(paste0("No Bivariate prior is derived for alpha_pi = ", alpha_pi))
                    }
            
          }
          ##
          ## ---- The direct Se/Sp prior: type 4, N(0, 1) on probit Sp and probit Se at every cutoff ----------------------------------
          ##
          ## (the pop-median models too: with "Dirichlet" their own fields, the prior off, replace these below)
          if (model_parameterisation %in% c("ord_bivariate", "R&G") && 
              alpha_pi == "one" && 
              # !any(candidate_route)) {
              (!any(candidate_route) || is_a_pop_median_model)) {
    
                    # priors$use_Se_Sp_prior <- 1 ## switches the direct Se/Sp prior on. With 0, the model uses its own Dirichlet prior on the cutpoints, unchanged.
                    # priors$use_Se_Sp_prior <- 0
                    priors$use_Se_Sp_prior <- if (ordinal_prior_type == "Se_Sp") 1 else 0
                    ##
                    priors$Se_Sp_prior_type <- 4 ## type 4 is a Normal on the probit scale of Sp and Se at the chosen cutoffs. Type 1 is a Dirichlet, types 2-3 are Betas.
                    ##
                    priors$n_chosen_max <- n_thr ## the row count of the padded arrays below; here every cutoff.
                    priors$n_chosen <- c(n_thr, n_thr) ## how many cutoffs get a prior in each class: all of them for Sp and all of them for Se.
                    ##
                    priors$chosen_cutoffs <- rbind(1:n_thr, 1:n_thr) ## which cutoffs, one row per class (row 1 Sp, row 2 Se), in rising order.
                    ##
                    ## the type 1 numbers. Type 4 never reads them, but the Stan data block declares them, so they must be passed. Zeros pass its <lower=0.0> bound:
                    priors$prior_Se_Sp_dirichlet_A <- list(rep(1, n_thr + 1), rep(1, n_thr + 1)) 
                    ##
                    priors$prior_Se_Sp_beta <- list(matrix(0, n_thr, 2), matrix(0, n_thr, 2)) ## the type 2-3 Beta numbers. Same reason: declared but unused under type 4.
                    ##
                    ## the type 4 numbers. Column 1 is the mean and column 2 the SD of probit(Sp_k) for class 1 and probit(Se_k) for class 2. Here mean 0 and SD 1 at every cutoff:
                    priors$prior_Se_Sp_link_mean_SD <- list(cbind(rep(0, n_thr), rep(1, n_thr)),
                                                            cbind(rep(0, n_thr), rep(1, n_thr)))

                    ## (fixed cutpoints only): centres the cutpoints at the baseline location, so the prior sits on Phi(C_k − beta), the summary Sp and Se. 
                    ## The random-cutpoint files have no such field, hence the "if":
                    if (random_thresholds == FALSE) {
                      priors$anchor_at_location <- 1
                    }

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

            # model_samples_obj <- model_prep_obj$sample( seed = MCMC_seed, ## MCMC seed
            #                                             n_burnin = n_burnin,
            #                                             n_iter = n_iter,
            #                                             # adapt_delta = 0.65,
            #                                             ## (adapt_delta 0.80 in every fit)
            #                                             adapt_delta = 0.80,
            #                                             max_treedepth = 10,
            #                                             metric_shape = "diag_e",
            #                                             priors = priors,
            #                                             init_lists_per_chain = init_lists_per_chain,
            #                                             n_chains = n_chains)



            ## (the nuisance blocks come from MetaOrdDTA itself, never written here: one source of truth)
            nuisance_params <- MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                  name_of_the_model      = sub("[.]stan$", "", custom_file_name),
                  model_parameterisation = model_parameterisation,
                  random_thresholds      = random_thresholds)
            ##
            ##
            interest_only_params <- if (model_parameterisation == "ord_bivariate") {
              
                                              if (random_thresholds) c("beta_mu", "dirichlet_phi")
                                              else c("beta_mu", "C_raw_vec")
              
                                    } else if (model_parameterisation == "R&G") {
                                      
                                              if (random_thresholds) c("beta_mu", "raw_scale_mu", "dirichlet_phi")
                                              else c("beta_mu", "raw_scale_mu", "C_raw_vec")
                                      
                                    } else {
                                          
                                              NULL
                                      
                                    }
            ##
            # if (any(candidate_route) && is_a_pop_median_model) {
            if (any(candidate_route) && is_a_pop_median_model && ordinal_prior_type == "Dirichlet") {
                priors[names(prior_fields)] <- prior_fields
            }
            ##
            ## ---- The centred cutpoints: each C_v starts at the cutpoints of the off-centred start (C_raw) --------------------------
            ##
            ## (C from C_raw as in the off-centred model, softplus off; the category probabilities of C at the
            ##  anchor, the prior mean of the intercept; C_v = alpha * log(probabilities), alpha = 0.01 +
            ##  dirichlet_phi * kappa of the same start, so softmax(C_v ./ alpha) gives the same C back)
            if (isTRUE(grepl("_centred_C_", custom_file_name, fixed = TRUE))) {
                  for (kk in 1:n_chains) {
                        C_v_of_this_chain <- list()
                        for (c in 1:2) {
                              C_raw_of_this_class <- init_lists_per_chain[[kk]]$C_raw[[c]]
                              C_of_this_class <- C_raw_of_this_class
                              if (ncol(C_raw_of_this_class) > 1) {
                                    for (k in 2:ncol(C_raw_of_this_class)) {
                                          C_of_this_class[, k] <- C_of_this_class[, k - 1] +
                                                                  exp(C_raw_of_this_class[, k])
                                    }
                              }
                              P_of_this_class <- pnorm(C_of_this_class - priors$prior_beta_mu_mean[c, 1])
                              probabilities_of_this_class <- t(apply(cbind(0, P_of_this_class, 1), 1, diff))
                              alpha_of_this_class <- 0.01 + init_lists_per_chain[[kk]]$dirichlet_phi[[c]] *
                                                            init_lists_per_chain[[kk]]$kappa[c]
                              C_v_of_this_chain[[c]] <- sweep(log(probabilities_of_this_class), 2,
                                                              alpha_of_this_class, "*")
                              stopifnot(all(is.finite(C_v_of_this_chain[[c]])))
                        }
                        init_lists_per_chain[[kk]]$C_v <- C_v_of_this_chain
                        init_lists_per_chain[[kk]]$C_raw <- NULL
                  }
            }
            time_after_the_model_prep <- Sys.time()
            model_samples_obj <- model_prep_obj$sample( seed = MCMC_seed, ## MCMC seed
                                            n_burnin = n_burnin,
                                            n_iter = n_iter,
                                            adapt_delta = 0.80,
                                            # adapt_delta = 0.655,
                                            max_treedepth = 10,
                                            # max_treedepth = 9,
                                            metric_shape = "diag_e",
                                            priors = priors,
                                            init_lists_per_chain = init_lists_per_chain,
                                            n_chains = n_chains,
                                            ##
                                            algorithm = sampling_algorithm,
                                            # algorithm = "NicoStan",
                                            # algorithm = "Stan",
                                            ##
                                                  NicoStan_args = list(
                                            ## (every $sample setting that changes sampling, written out)
                                            adapt_delta                          = 0.80,  ## the PS7 value 
                                            ###
                                            # metric_shape_main                    = "dense", ## --------
                                            metric_shape_main                    = "diag",    ## --------
                                            ##
                                            nuisance_params                      = nuisance_params,
                                            sample_nuisance                      = !is.null(nuisance_params),
                                            diffusion_HMC                        = TRUE,
                                            diffusion_HMC_integrator             = "kick_flow_kick",  ## signature
                                            partitioned_HMC                      = FALSE,  ## R_fn_sample.R:1606
                                            ##
                                            burnin_algorithm                     = "ESJD", ## --------
                                            ##
                                            n_adapt                              = n_burnin - round(n_burnin / 10),
                                            clip_iter                            = NULL, ## 20,
                                            clip_iter_tau                        = NULL, ## 50,
                                            n_refresh                            = round(n_burnin / 20),
                                            parallel_method                      = "RcppParallel",
                                            ##
                                            n_superchains                        = 4,
                                            ##
                                            learning_rate                        = 0.05,
                                            # learning_rate                        = 0.025,
                                            ##
                                            tau_mult                             = 1.60,  ## R_fn_sample.R:1355
                                            manual_tau                           = FALSE,  ## R_fn_sample.R:1367
                                            tau_if_manual                        = 3.0,  ## R_fn_sample.R:1368
                                            ##
                                            tau_weight_by_p_jump                 = TRUE,  ## R_fn_sample.R:1421
                                            tau_ramp                             = "staged",
                                            eps_reinit_at_ChEES_handover         = FALSE, ## ----
                                            ##
                                            eps_acceptance_mean                  = "harmonic",
                                            ##
                                            eps_reinit_after_pre_burnin          = TRUE,  ## R_fn_sample.R:1425 ## --------
                                            ##
                                            tau_shrink_on_divergence             = FALSE,  ## R_fn_sample.R:1426
                                            tau_shrink_factor                    = 0.95,  ## R_fn_sample.R:1427
                                            tau_shrink_min_divergent_chains      = 2,  ## R_fn_sample.R:1428
                                            ##
                                            share_tau_ii_across_chains_in_burnin = FALSE,  ## R_fn_sample.R:1429
                                            ##
                                            tau_initial                          = "adaptive", ## pi ## --------
                                            tau_adaptation_block                 = "interest_only",  ## --------
                                            interest_only                        = interest_only_params,
                                            ##
                                            tau_sampling_scale                   = "none",  ## R_fn_sample.R:1446
                                            tau_jitter_burnin                    = "uniform",  ## signature
                                            tau_jitter_sampling_type             = "uniform",
                                            ##
                                            randomize_tau_burnin                 = TRUE,  ## signature
                                            randomize_tau_sampling               = TRUE,  ## signature
                                            ##
                                            theta_hat_us_rule                    = "running_mean_frozen",
                                            ##
                                            metric_type_main                     = "Empirical",
                                            ratio_M_main                         = 0.90,  ## R_fn_sample.R:1686
                                            interval_width_main                  = 1,  ## R_fn_sample.R:1689
                                            ##
                                            metric_type_nuisance                 = "Empirical",
                                            metric_shape_nuisance                = "diag",  ## R_fn_sample.R:1697
                                            ratio_M_nuisance                     = 0.90,  ## R_fn_sample.R:1700
                                            interval_width_nuisance              = 1,  ## R_fn_sample.R:1703
                                            ##
                                            max_eps_main                         = 1.0,  ## R_fn_sample.R:1708
                                            max_eps_nuisance                     = 1.0,  ## R_fn_sample.R:1709
                                            max_tau_main                         = 50,
                                            max_tau_nuisance                     = 50,
                                            max_L                                = 1024,  ## R_fn_sample.R:1764
                                            ##
                                            ## (M_decay_*: helper fns 433 type, 436 scale, 440 power)
                                            M_decay_type                         = "inverse",
                                            M_decay_power                        = 0.5,
                                            M_decay_scale                        = (n_burnin - round(n_burnin / 10)) / 5,
                                            ##
                                            beta1_adam                           = 0,  ## R_fn_sample.R:1620
                                            beta2_adam                           = 0.95,  ## R_fn_sample.R:1621
                                            eps_adam                             = 1e-8,  ## R_fn_sample.R:1622
                                            ##
                                            force_autodiff                       = FALSE,  ## R_fn_sample.R:1624
                                            force_PartialLog                     = FALSE,  ## R_fn_sample.R:1625
                                            multi_attempts                       = TRUE,  ## R_fn_sample.R:1626
                                            ##
                                            metric_estimator                     = "pooled",  ## signature
                                            metric_pooled_offdiagonal_shrinkage  = "adaptive",  ## signature
                                            pre_burnin_n_iter                    = NULL,  ## signature
                                            pre_burnin_L                         = NULL,  ## signature
                                            ##
                                            burnin_schedule                      = "automatic",  ## signature
                                            esjd_jump_power                      = 2,  ## signature
                                            tau_cost_exponent                    = 1,  ## signature
                                            ##
                                            metric_update_rescales_eps_and_tau   = TRUE, ## signature ## ----
                                            ##
                                            pooled_metric_warm_up                = 8, ## says how many draws it waits for before its first estimate replaces the starting identity metric.
                                            metric_start_iter                    = round(n_burnin/5) ## says WHEN the metric starts collecting draws.
                                            ))
            ##SO
            ## ---- Get summaries:
            ##
            time_after_the_sample_call <- Sys.time()
            model_summary_obj <- model_samples_obj$summary( compute_main_params = TRUE,
                                                            compute_transformed_parameters = FALSE,
                                                            compute_generated_quantities = TRUE,
                                                            save_log_lik_trace = FALSE,
                                                            use_BayesMVP_for_faster_summaries = TRUE)
            time_after_the_summary <- Sys.time()
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
              # min_ESS_threshold <- ifelse(n_studies == 10,
              #                             5000,
              #                             1000)
              ## (the min ESS set once in the runner: one value for 10 studies, one for 50 studies)
              min_ESS_threshold <- if (n_studies == 10) min_ESS_of_10_studies else min_ESS_of_50_studies
              # min_ESS_threshold <- ifelse(n_studies == 10,
              #                             1000, ## 500, ## 250,
              #                             200 ## 100 ## 50
              #                             )
              ##
              ## ---- % of divergent sampling transitions (fits above max_percent_divergent are not kept):
              ##
              max_percent_divergent <- 1
              ##
              n_divergent_per_chain <- model_summary_obj$get_divergences()$n_divs_per_chain
              if (is.null(n_divergent_per_chain)) {
                stop("the number of divergent transitions of this fit was not found")
              }
              n_divergent <- sum(n_divergent_per_chain)
              percent_divergent <- 100 * n_divergent / (n_chains * n_iter)
              cat(sprintf("Seed %d: min ESS = %.0f (threshold %.0f), divergent = %d of %d (%.3f%%, max %.1f%%)\n",
                          seed,
                          min_ESS,
                          min_ESS_threshold,
                          n_divergent,
                          n_chains * n_iter,
                          percent_divergent,
                          max_percent_divergent))
              ##
              ## ---- one row per fit, kept or dropped, in the gate log of this cell:
              ##
              file_of_the_gate_log <- file.path(output_dir,
                                                paste0("gate_of_every_fit_", file_name_string, ".csv"))
              efficiency_of_this_fit <- tryCatch(model_summary_obj$get_efficiency_metrics(),
                                                 error = function(e) list())
              seconds_burnin         <- c(efficiency_of_this_fit$time_burnin, NA)[1]
              seconds_sampling       <- c(efficiency_of_this_fit$time_sampling, NA)[1]
              n_grad_evals_sampling  <- c(efficiency_of_this_fit$Min_ESS_main /
                                          efficiency_of_this_fit$Min_ess_per_grad_main_samp, NA)[1]
              min_ESS_per_1000_grads <- 1000 * min_ESS / n_grad_evals_sampling
              seconds_of_each_step   <- tryCatch(diff(as.numeric(c(time_of_the_start_of_this_fit,
                                                     time_after_the_simulated_data,
                                                     time_after_the_model_prep,
                                                     time_after_the_sample_call,
                                                     time_after_the_summary,
                                                     Sys.time()))),
                                                 error = function(e) rep(NA_real_, 5))
              suppressWarnings(utils::write.table(data.frame(seed = seed,
                                                             min_ESS = min_ESS,
                                                             min_ESS_threshold = min_ESS_threshold,
                                                             n_divergent = n_divergent,
                                                             n_transitions = n_chains * n_iter,
                                                             percent_divergent = percent_divergent,
                                                             kept = (min_ESS > min_ESS_threshold) &&
                                                                    (percent_divergent <= max_percent_divergent),
                                                             seconds_burnin = seconds_burnin,
                                                             seconds_sampling = seconds_sampling,
                                                             n_grad_evals_sampling = n_grad_evals_sampling,
                                                             min_ESS_per_1000_grads = min_ESS_per_1000_grads,
                                                             seconds_data = seconds_of_each_step[1],
                                                             seconds_model_prep = seconds_of_each_step[2],
                                                             seconds_sample_call = seconds_of_each_step[3],
                                                             seconds_summary = seconds_of_each_step[4],
                                                             seconds_summary_to_gate = seconds_of_each_step[5]),
                                                  file = file_of_the_gate_log,
                                                  sep = ",",
                                                  row.names = FALSE,
                                                  col.names = !file.exists(file_of_the_gate_log),
                                                  append = TRUE))
              ##
              ## ---- the step size and trajectory length of this fit, in the HMC log of this cell:
              ##
              HMC_info_of_this_fit <- tryCatch(model_summary_obj$get_HMC_info(),
                                               error = function(e) list())
              HMC_info_of_the_sampling <- HMC_info_of_this_fit$stan_HMC_info_list$out_list_sampling
              L_max_per_chain <- HMC_info_of_the_sampling$L_max_per_chain_sampling
              L_mean_of_this_fit <- c(HMC_info_of_this_fit$L_main_during_sampling, NA)[1]
              file_of_the_HMC_log <- file.path(output_dir,
                                               paste0("HMC_of_every_fit_", file_name_string, ".csv"))
              suppressWarnings(utils::write.table(data.frame(seed = seed,
                                                             sampling_algorithm = sampling_algorithm,
                                                             n_burnin = n_burnin,
                                                             eps_sampling = c(HMC_info_of_this_fit$eps, NA)[1],
                                                             tau_sampling = c(HMC_info_of_this_fit$tau, NA)[1],
                                                             L_mean_sampling = L_mean_of_this_fit,
                                                             L_max_sampling = if (length(L_max_per_chain) > 0)
                                                                                max(L_max_per_chain) else NA),
                                                  file = file_of_the_HMC_log,
                                                  sep = ",",
                                                  row.names = FALSE,
                                                  col.names = !file.exists(file_of_the_HMC_log),
                                                  append = TRUE))
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

          ## ---- ord_bivariate model, fixed-cutpoints:
          ## 0.172 ## w/ jitte = 0.01 (or 0.1?!) (8 chains)
          ## 0.114 w/ jitter = 0.001 (8 chains) - min_ESS ~ 13,000
          ## 0.194 w/ jitter = 0.001 (4 chains) - min_ESS ~ 6500
          ## ---- ord_bivariate model, RANDOM-cutpoint (all 4 chains):
          ## 0.172 w/ df = 5, log(100), sd = 1 kappa priors
          ## 0.13 w/ df = c(5, 3)

          # priors$prior_kappa_dof


          #
          # try({
          #   model_summary_obj$extract_params(params = "kappa")
          # })


          if (min_ESS > min_ESS_threshold && percent_divergent <= max_percent_divergent) { # Store results

            results <- list(
              seed = seed,
              ##
              percent_divergent = percent_divergent,
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

            # Check MCSE periodically (like the LC-MVP)
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
                                                              # min_ESS_thr = NULL,
                                                                min_ESS_thr = min_ESS_threshold,
                                                                # min_ESS_thr = if (n_studies == 10) 1000 else 200,
                                                                # min_ESS_thr = if (n_studies == 10) 500 else 100,
                                                                # min_ESS_thr = if (n_studies == 10) 250 else 50,
                                                              vec_index_inner = vec_index_inner,
                                                              min_studies_per_threshold = NULL)



                    ## (check_MCSE_status returns MCSE_val once it has min_runs_needed results, and
                    ##  mcse_value = NA before then; reading mcse_value alone stopped every worker at its first
                    ##  full check)
                    MCSE_value_of_this_check <- c(current_mcse_status$MCSE_val, current_mcse_status$mcse_value)[1]
                    cat(sprintf("Completed sims: %d, Current MCSE: %.4f%%\n",
                                current_mcse_status$n_sims,
                                MCSE_value_of_this_check * 100))

                    # Check if MCSE is valid AND target is met
                    if (is.finite(MCSE_value_of_this_check) &&
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





