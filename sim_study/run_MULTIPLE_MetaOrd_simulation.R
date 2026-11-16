
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
                              min_runs_needed,
                              target_MCSE_RMSE,
                              check_MCSE_interval,
                              stop_early_if_target_met,
                              sampling_algorithm,
                              min_ESS_of_10_studies,
                              min_ESS_of_50_studies,
                              ordinal_prior_type,
                              ## (real_data_dir: NULL, or the folder of the real data to fit in place of simulated data;
                              ##  n_iter_for_4_chains_of_this_run: NULL = the table of this file, or one number for every cell of a pilot;
                              ##  nuisance_params_of_this_run: NULL = the package's nuisance blocks of the model, or their names)
                              real_data_dir,
                              Stan_math_only_model,
                              log_lik_indicator,
                              n_burnin,
                              n_iter_for_4_chains_of_this_run,
                              nuisance_params_of_this_run,
                              ## (keep_the_data_file_in, keep_the_divergence_trace_in: NULL, or a folder to save the Stan data
                              ##  of every cmdstanr fit / the divergence trace of every NicoStan fit in;
                              ##  sampler_settings_to_test_file: NULL, or the .rds file of a list of sampler settings to test;
                              ##  target_proportion_missing_of_this_run: NULL = each test gets its real reporting pattern of
                              ##  missing thresholds, or the proportion of the thresholds to make missing (e.g. 0.85))
                              keep_the_data_file_in,
                              keep_the_divergence_trace_in,
                              sampler_settings_to_test_file,
                              target_proportion_missing_of_this_run
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
                prior_set <- if (is.null(scenario_row$prior_set)) "set_0" else scenario_row$prior_set
                # box_cox <- scenario_row$box_cox
                softplus <- FALSE ## scenario_row$softplus
                ##
                ## ---- THE CUTPOINT BOUNDS OF THE ORDINAL MODELS: set HERE, ONCE, and nowhere else in the sim. The step between two
                ##      neighbouring cutpoints is kept between exp(C_raw_lower) and exp(C_raw_upper). They reach every Stan model as its
                ##      data C_raw_lower and C_raw_upper: no model file holds a number.
                C_raw_lower <- -10
                C_raw_upper <- 2.5
                random_thresholds <- scenario_row$random_thresholds
                Dirichlet_random_effects_type <- scenario_row$Dirichlet_random_effects_type
                ##

                # Set threshold ranges
                if (test == "GAD_2") {
                  vec_index_inner <- 1:6
                } else if (test == "HADS") {
                  vec_index_inner <- 3:17
                } else if (test == "BAI") {
                  vec_index_inner <- 5:43
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
            ## (the pop-median models carry their Se / Sp prior fields)
            is_a_pop_median_model <- grepl("popmed", scenario_row$custom_file_name[[1]], fixed = TRUE)
            stopifnot(all(candidate_route[1:2]), candidate_route[3] == is_a_pop_median_model,
                      model_parameterisation %in% c("ord_bivariate", "ord_HSROC"),
                      isTRUE(random_thresholds), Dirichlet_random_effects_type == "kappa",
                      is.character(file_name_string), length(file_name_string) == 1,
                      !is.na(file_name_string),
                      grepl("popmed", file_name_string, fixed = TRUE) == is_a_pop_median_model)
            stopifnot(scenario_row$custom_file_name[[1]] %in%
                            c("DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan",
                              "DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan"))
            stopifnot(grepl("ord_HSROC", scenario_row$custom_file_name[[1]], fixed = TRUE) ==
                            (model_parameterisation == "ord_HSROC"))
            stopifnot(is.character(scenario_row$custom_cpp_user_header_file_path[[1]]),
                      length(scenario_row$custom_cpp_user_header_file_path[[1]]) == 1,
                      file.exists(scenario_row$custom_cpp_user_header_file_path[[1]]))
            if (is_a_pop_median_model) {
                  prior_fields <- scenario_row$prior_fields[[1]]
                  stopifnot(is.list(prior_fields), !is.null(names(prior_fields)),
                            !anyDuplicated(names(prior_fields)),
                            setequal(names(prior_fields),
                                     c("use_Se_Sp_prior", "Se_Sp_prior_type", "n_chosen_cutoffs_max",
                                       "n_chosen_cutoffs", "chosen_cutoffs", "prior_Se_Sp_dirichlet_A",
                                       "prior_Se_Sp_beta", "prior_Se_Sp_link_mean_SD")))
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

          ## (a seed is not started when the check of a lane found this cell's MCSE target met in the last 10 minutes; the
          ##  time now is read as the time of a file written now in the same folder: one clock, the folder's)
          if (stop_early_if_target_met && file.exists(file.path(output_dir, paste0("MCSE_target_met_", file_name_string, ".txt")))) {
                file_written_now <- file.path(output_dir, paste0("time_now_of_", Sys.info()[["nodename"]], "_process_", Sys.getpid(), ".txt"))
                writeLines("", file_written_now)
                seconds_since_the_target_was_found_met <- as.numeric(file.mtime(file_written_now)) -
                      as.numeric(file.mtime(file.path(output_dir, paste0("MCSE_target_met_", file_name_string, ".txt"))))
                unlink(file_written_now)
                if (isTRUE(seconds_since_the_target_was_found_met < 600)) {
                      cat("TARGET MCSE ACHIEVED (the check of another lane found it)! This seed is not started.\n")
                      return(list(stopped_early = TRUE,
                                  final_n_sims = NA_integer_))
                }
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

          if (!is.null(real_data_dir)) {
                real_x <- readRDS(file.path(real_data_dir,
                                           paste0("real_x_", c(GAD_2 = "GAD2", HADS = "HADS", BAI = "BAI")[[test]], ".RDS")))
                x <- list(real_x$nondiseased, real_x$diseased) ## the model takes the non-diseased class first
          }
          n_thr <- ncol(x[[1]]) - 1
          n_studies_actual <- nrow(x[[1]])

          use_custom_file <- FALSE
          custom_file_name <- NULL
          ##
          ## (one file per model: its fused AVX file; ordinal and Jones models: use_Se_Sp_prior = 0 is the model's own prior;
          ##  the stratified bivariate has one prior only, the Normal on probit Se and Sp, so no such switch)
          if (model_parameterisation == "Bivariate") {
              use_custom_file  <- TRUE
              custom_file_name <- "DTA_MA_strat_bivariate_NicoStan_fused_AVX.stan"
          }
          ## (the stratified bivariate is sampled from its fused file)
          Stan_file_of_the_fused_version <- NULL
          if (model_parameterisation == "Bivariate") {
              Stan_file_of_the_fused_version <- system.file("stan_models", "stan_models_MA",
                                                            "DTA_MA_strat_bivariate_NicoStan_fused_AVX.stan",
                                                            package = "MetaOrdDTA")
              stopifnot(nzchar(Stan_file_of_the_fused_version))
          }
          if (model_parameterisation == "ord_bivariate") {
              use_custom_file  <- TRUE
              custom_file_name <- paste0("DTA_MA_ord_biv_",
                                         if (random_thresholds) "RANDthr_kappa_popmed" else "FIXEDthr",
                                         "_Se_Sp_prior_NicoStan_fused_AVX.stan")
          }
          ## (the ordinal-HSROC model has random cutpoints only)
          if (model_parameterisation == "ord_HSROC") {
              stopifnot(random_thresholds)
              use_custom_file  <- TRUE
              custom_file_name <- "DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan"
          }
          if (model_parameterisation == "Jones") {
              use_custom_file  <- TRUE
              custom_file_name <- "DTA_MA_Jones_Se_Sp_prior_NicoStan_fused_AVX.stan"
          }
          if (any(candidate_route)) {
              custom_file_name <- scenario_row$custom_file_name[[1]]
              use_custom_file <- TRUE
          }
          ## (Stan_math_only_model: the model's own Stan file, Stan::math only, no user header)
          stopifnot(is.logical(Stan_math_only_model), length(Stan_math_only_model) == 1)
          custom_file_name_to_compile <- custom_file_name
          if (Stan_math_only_model && use_custom_file) {
              custom_file_name_to_compile <- sub("_NicoStan_fused_AVX[.]stan$", ".stan", custom_file_name)
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
                                        ## (the cutpoint bounds, set once at the top of this file; the class puts them in the Stan data)
                                        C_raw_lower = C_raw_lower,
                                        C_raw_upper = C_raw_upper,
                                        ##
                                        custom_file_name = custom_file_name_to_compile,
                                        ##
                                        advanced_compile = use_custom_file,
                                        force_recompile = FALSE,
                                        quiet = TRUE,
                                        set_custom_CXX_CPP_flags = use_custom_file,
                                        custom_cpp_user_header_file_path = if (use_custom_file && !Stan_math_only_model)
                                                                  NicoSIMD_user_header else NULL,
                                        FMA_FLAGS = "-mfma",
                                        AVX_FLAGS = AVX_FLAGS_of_this_machine)

          ## (log_lik_indicator = 0: log_lik_study is not computed in generated quantities)
          stopifnot(log_lik_indicator %in% c(0, 1))
          model_prep_obj$advanced_model_options$log_lik_indicator <- as.integer(log_lik_indicator)
          model_prep_obj$internal_obj$outs_data$stan_data_list$log_lik_indicator <- as.integer(log_lik_indicator)

          ## (the cutpoint bounds go to the Stan model as data: the two values set once at the top of this function)
          ## (check: the Stan data of the model hold exactly the bounds set at the top of this file)
          if (!isTRUE(all.equal(as.numeric(model_prep_obj$internal_obj$outs_data$stan_data_list$C_raw_lower), as.numeric(C_raw_lower))) ||
              !isTRUE(all.equal(as.numeric(model_prep_obj$internal_obj$outs_data$stan_data_list$C_raw_upper), as.numeric(C_raw_upper)))) {
                stop("the cutpoint bounds in the Stan data are not the bounds set at the top of the sim runner")
          }
          print(paste("cutpoint bounds in the Stan data of the model: C_raw_lower =", model_prep_obj$internal_obj$outs_data$stan_data_list$C_raw_lower,
                      ", C_raw_upper =", model_prep_obj$internal_obj$outs_data$stan_data_list$C_raw_upper))
          
          # Set priors
          priors <- model_prep_obj$priors

          ## (every prior of every model is set in ONE block further down: "THE PRIORS OF EACH MODEL")

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

          ##
          ## ---- THE SAMPLING ITERATIONS PER CHAIN (4 chains) OF EACH CELL: the ONE place where they are set ------------------------------
          ##
          ## (from the pilot: 10 fits of 1000 iterations per chain of each of the 216 cells - prior sets 0, 1, 2 x
          ##  O-biv-FC, O-biv-RC, ord-HSROC, stratified bivariate x the 3 ordinal DGMs x GAD-2, HADS, BAI x 10 and 50 studies - under
          ##  the priors, the cutpoint bounds (-10, 2.5) and the sampler settings of this file. The rule: the min ESS (Se, Sp) asked of
          ##  the cell (2500 at 10 studies, 500 at 50) over the 10th percentile of its pilot fits' min ESS per iteration, rounded up
          ##  to 50, at least 250, at most 15000. A cell that is not in this table has no number: the runner stops.)
          n_iter_of_each_cell <- as.data.frame(matrix(byrow = TRUE, ncol = 7, data = list(
                ## prior_set, model, random_thresholds, n_studies, test, DGM, n_iter
                ## ---- O-biv-FC, O-biv-RC and ord-HSROC at prior sets 1 and 2 with the Se / Sp prior of these files: every cell of these
                ##      two sets has its number from a pilot of 10 fits of 1000 iterations per chain of the cell (the rows of these
                ##      cells further down the table included); the rule: the fewest iterations (a multiple of 50) at which at least
                ##      90% of the pilot fits reach the min ESS asked of the cell, so that between 85% and 95% of the fits are kept.
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "Hoyer_2018",       4000,
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "diagmeta",         3350,
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "logit_GLMM",       3200,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "Hoyer_2018",       1050,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "diagmeta",         1050,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "logit_GLMM",       1250,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "Hoyer_2018",       900,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "diagmeta",         1500,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "logit_GLMM",       1200,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "Hoyer_2018",       2750,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "diagmeta",         3300,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "logit_GLMM",       3150,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "Hoyer_2018",       2950,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "diagmeta",         3050,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "logit_GLMM",       3150,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "Hoyer_2018",       2350,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "diagmeta",         2250,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "logit_GLMM",       2650,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "Hoyer_2018",       1500,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "diagmeta",         2050,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "logit_GLMM",       1550,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "Hoyer_2018",       1500,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "diagmeta",         1650,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "logit_GLMM",       1750,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "Hoyer_2018",       1850,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "diagmeta",         1500,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "logit_GLMM",       1850,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "Hoyer_2018",       1850,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "diagmeta",         2750,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "logit_GLMM",       2450,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "Hoyer_2018",       6950,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "diagmeta",         5750,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "logit_GLMM",       5850,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "Hoyer_2018",       3750,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "diagmeta",         3800,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "logit_GLMM",       4850,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "Hoyer_2018",       3950,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "diagmeta",         4950,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "logit_GLMM",       4600,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "Hoyer_2018",       1550,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "diagmeta",         2300,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "logit_GLMM",       1700,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "Hoyer_2018",       1500,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "diagmeta",         1650,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "logit_GLMM",       2000,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "Hoyer_2018",       2500,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "diagmeta",         1650,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "logit_GLMM",       1750,
                "set_0", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  4700,
                "set_0", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 7750,
                "set_0", "ord_bivariate", FALSE, 10, "BAI",   "ord_HSROC_random_C",     5700,
                "set_0", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3950,
                "set_0", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 3850,
                "set_0", "ord_bivariate", FALSE, 10, "GAD_2", "ord_HSROC_random_C",     3800,
                "set_0", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  5050,
                "set_0", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 4700,
                "set_0", "ord_bivariate", FALSE, 10, "HADS",  "ord_HSROC_random_C",     3750,
                "set_0", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2550,
                "set_0", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2700,
                "set_0", "ord_bivariate", FALSE, 50, "BAI",   "ord_HSROC_random_C",     2600,
                "set_0", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  2150,
                "set_0", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2550,
                "set_0", "ord_bivariate", FALSE, 50, "GAD_2", "ord_HSROC_random_C",     1750,
                "set_0", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  1900,
                "set_0", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2300,
                "set_0", "ord_bivariate", FALSE, 50, "HADS",  "ord_HSROC_random_C",     1400,
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  5650,
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 5150,
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "ord_HSROC_random_C",     4100,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  2050,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 2200,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "ord_HSROC_random_C",     3550,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  6850,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 6750,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "ord_HSROC_random_C",     3400,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  3550,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 2000,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "ord_HSROC_random_C",     2500,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  850,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 500,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "ord_HSROC_random_C",     500,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  3750,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 2800,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "ord_HSROC_random_C",     1350,
                "set_0", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_fixed_C",  6250,
                "set_0", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_random_C", 6250,
                "set_0", "ord_HSROC",           TRUE,  10, "BAI",   "ord_HSROC_random_C",     7350,
                "set_0", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_fixed_C",  6900,
                "set_0", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_random_C", 5800,
                "set_0", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_HSROC_random_C",     5150,
                "set_0", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_fixed_C",  6000,
                "set_0", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_random_C", 4950,
                "set_0", "ord_HSROC",           TRUE,  10, "HADS",  "ord_HSROC_random_C",     5950,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_fixed_C",  4550,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_random_C", 3400,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "ord_HSROC_random_C",     5000,
                "set_0", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_fixed_C",  1400,
                "set_0", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_random_C", 1550,
                "set_0", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_HSROC_random_C",     1050,
                "set_0", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_fixed_C",  2100,
                "set_0", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_random_C", 1900,
                "set_0", "ord_HSROC",           TRUE,  50, "HADS",  "ord_HSROC_random_C",     3050,
                "set_0", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  7800,
                "set_0", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 9150,
                "set_0", "Bivariate",     FALSE, 10, "BAI",   "ord_HSROC_random_C",     8500,
                "set_0", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  6050,
                "set_0", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 6850,
                "set_0", "Bivariate",     FALSE, 10, "GAD_2", "ord_HSROC_random_C",     4300,
                "set_0", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  7350,
                "set_0", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 6600,
                "set_0", "Bivariate",     FALSE, 10, "HADS",  "ord_HSROC_random_C",     7800,
                "set_0", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1850,
                "set_0", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1900,
                "set_0", "Bivariate",     FALSE, 50, "BAI",   "ord_HSROC_random_C",     1200,
                "set_0", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  700,
                "set_0", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 800,
                "set_0", "Bivariate",     FALSE, 50, "GAD_2", "ord_HSROC_random_C",     650,
                "set_0", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1950,
                "set_0", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1450,
                "set_0", "Bivariate",     FALSE, 50, "HADS",  "ord_HSROC_random_C",     1100,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  2500,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 3100,
                "set_1", "ord_bivariate", FALSE, 10, "BAI",   "ord_HSROC_random_C",     2750,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3200,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 2800,
                "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "ord_HSROC_random_C",     2750,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  2850,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 2350,
                "set_1", "ord_bivariate", FALSE, 10, "HADS",  "ord_HSROC_random_C",     2300,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2300,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2100,
                "set_1", "ord_bivariate", FALSE, 50, "BAI",   "ord_HSROC_random_C",     1850,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  1900,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2000,
                "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "ord_HSROC_random_C",     1600,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  1850,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2250,
                "set_1", "ord_bivariate", FALSE, 50, "HADS",  "ord_HSROC_random_C",     1150,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  3550,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 3200,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "ord_HSROC_random_C",     5050,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  2050,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 1850,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "ord_HSROC_random_C",     1900,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  4650,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 10250,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "ord_HSROC_random_C",     8900,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  3300,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 1750,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "ord_HSROC_random_C",     2000,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  650,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 500,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "ord_HSROC_random_C",     450,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  2850,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 2500,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "ord_HSROC_random_C",     1100,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_fixed_C",  4650,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_random_C", 4800,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "ord_HSROC_random_C",     5000,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_fixed_C",  4600,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_random_C", 3350,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_HSROC_random_C",     3750,
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_fixed_C",  3250,
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_random_C", 3100,
                "set_1", "ord_HSROC",           TRUE,  10, "HADS",  "ord_HSROC_random_C",     5200,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_fixed_C",  2600,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_random_C", 2600,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "ord_HSROC_random_C",     3100,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_fixed_C",  1050,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_random_C", 1000,
                "set_1", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_HSROC_random_C",     1200,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_fixed_C",  1500,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_random_C", 1800,
                "set_1", "ord_HSROC",           TRUE,  50, "HADS",  "ord_HSROC_random_C",     2300,
                "set_1", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  3850,
                "set_1", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 5300,
                "set_1", "Bivariate",     FALSE, 10, "BAI",   "ord_HSROC_random_C",     4100,
                "set_1", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  4050,
                "set_1", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 4400,
                "set_1", "Bivariate",     FALSE, 10, "GAD_2", "ord_HSROC_random_C",     4100,
                "set_1", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  5450,
                "set_1", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 3800,
                "set_1", "Bivariate",     FALSE, 10, "HADS",  "ord_HSROC_random_C",     4650,
                "set_1", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1600,
                "set_1", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1700,
                "set_1", "Bivariate",     FALSE, 50, "BAI",   "ord_HSROC_random_C",     1100,
                "set_1", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  700,
                "set_1", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 650,
                "set_1", "Bivariate",     FALSE, 50, "GAD_2", "ord_HSROC_random_C",     700,
                "set_1", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1750,
                "set_1", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1250,
                "set_1", "Bivariate",     FALSE, 50, "HADS",  "ord_HSROC_random_C",     1050,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  4800,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 7950,
                "set_2", "ord_bivariate", FALSE, 10, "BAI",   "ord_HSROC_random_C",     3100,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  4000,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 3250,
                "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "ord_HSROC_random_C",     3700,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  3600,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 4400,
                "set_2", "ord_bivariate", FALSE, 10, "HADS",  "ord_HSROC_random_C",     2950,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2200,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2300,
                "set_2", "ord_bivariate", FALSE, 50, "BAI",   "ord_HSROC_random_C",     1950,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  2000,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2100,
                "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "ord_HSROC_random_C",     1750,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  1800,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2100,
                "set_2", "ord_bivariate", FALSE, 50, "HADS",  "ord_HSROC_random_C",     1200,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  5850,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 8200,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "ord_HSROC_random_C",     4000,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  2700,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 2050,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "ord_HSROC_random_C",     2050,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  5850,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 5850,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "ord_HSROC_random_C",     4300,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  2950,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 2000,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "ord_HSROC_random_C",     2000,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  750,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 450,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "ord_HSROC_random_C",     400,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  3400,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 2400,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "ord_HSROC_random_C",     1600,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_fixed_C",  4000,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "ord_biv_random_C", 4550,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "ord_HSROC_random_C",     4700,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_fixed_C",  16050,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_biv_random_C", 4250,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "ord_HSROC_random_C",     4750,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_fixed_C",  4800,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "ord_biv_random_C", 5050,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "ord_HSROC_random_C",     4500,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_fixed_C",  2700,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "ord_biv_random_C", 3200,
                "set_2", "ord_HSROC",           TRUE,  50, "BAI",   "ord_HSROC_random_C",     3650,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_fixed_C",  1450,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_biv_random_C", 1150,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "ord_HSROC_random_C",     1350,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_fixed_C",  1450,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "ord_biv_random_C", 1600,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "ord_HSROC_random_C",     2350,
                "set_2", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  4900,
                "set_2", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 5150,
                "set_2", "Bivariate",     FALSE, 10, "BAI",   "ord_HSROC_random_C",     4800,
                "set_2", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3250,
                "set_2", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 4550,
                "set_2", "Bivariate",     FALSE, 10, "GAD_2", "ord_HSROC_random_C",     4950,
                "set_2", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  4850,
                "set_2", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 5100,
                "set_2", "Bivariate",     FALSE, 10, "HADS",  "ord_HSROC_random_C",     4550,
                "set_2", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1600,
                "set_2", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1800,
                "set_2", "Bivariate",     FALSE, 50, "BAI",   "ord_HSROC_random_C",     1150,
                "set_2", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  650,
                "set_2", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 700,
                "set_2", "Bivariate",     FALSE, 50, "GAD_2", "ord_HSROC_random_C",     800,
                "set_2", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1750,
                "set_2", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1300,
                "set_2", "Bivariate",     FALSE, 50, "HADS",  "ord_HSROC_random_C",     950,
                ## ---- ord-HSROC with random cutpoints on the DGMs Hoyer_2018, logit_GLMM and diagmeta: the cells still to finish,
                ##      from a pilot of 10 fits of 1000 iterations per chain of each cell with the partially centred files; the
                ##      rule: the fewest iterations (a multiple of 50) at which at least 90% of the pilot fits reach the min ESS
                ##      asked of the cell, so that between 85% and 95% of the fits are kept.
                "set_0", "ord_HSROC",           TRUE,  10, "BAI",   "diagmeta",         5550,
                "set_0", "ord_HSROC",           TRUE,  10, "BAI",   "logit_GLMM",       5350,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "Hoyer_2018",       2000,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "diagmeta",         2750,
                "set_0", "ord_HSROC",           TRUE,  50, "BAI",   "logit_GLMM",       2600,
                "set_0", "ord_HSROC",           TRUE,  50, "HADS",  "logit_GLMM",       1300,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "Hoyer_2018",       3800,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "diagmeta",         4500,
                "set_1", "ord_HSROC",           TRUE,  10, "BAI",   "logit_GLMM",       3700,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "Hoyer_2018",       3200,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "diagmeta",         2650,
                "set_1", "ord_HSROC",           TRUE,  10, "GAD_2", "logit_GLMM",       3400,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "Hoyer_2018",       1750,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "diagmeta",         2400,
                "set_1", "ord_HSROC",           TRUE,  50, "BAI",   "logit_GLMM",       2400,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "Hoyer_2018",       5600,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "diagmeta",         4100,
                "set_2", "ord_HSROC",           TRUE,  10, "BAI",   "logit_GLMM",       10100,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "Hoyer_2018",       3600,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "diagmeta",         3250,
                "set_2", "ord_HSROC",           TRUE,  10, "GAD_2", "logit_GLMM",       4950,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "Hoyer_2018",       5050,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "diagmeta",         7050,
                "set_2", "ord_HSROC",           TRUE,  10, "HADS",  "logit_GLMM",       3750,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "Hoyer_2018",       1000,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "diagmeta",         850,
                "set_2", "ord_HSROC",           TRUE,  50, "GAD_2", "logit_GLMM",       1050,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "Hoyer_2018",       650,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "diagmeta",         1150,
                "set_2", "ord_HSROC",           TRUE,  50, "HADS",  "logit_GLMM",       1900,
                ## ---- O-biv-RC (ord_bivariate with random cutpoints) in the partially centred form on the DGMs Hoyer_2018, logit_GLMM
                ##      and diagmeta, from a pilot of 10 fits of 1000 iterations per chain of each cell, by the same rule as the rows
                ##      of ord-HSROC above (at least 90% of the pilot fits reach the min ESS asked of the cell). The rows of O-biv-RC on
                ##      the three ordinal DGMs higher up in this table are from the same pilot and the same rule.
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "Hoyer_2018",       4350,
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "diagmeta",         5700,
                "set_0", "ord_bivariate", TRUE,  10, "BAI",   "logit_GLMM",       6250,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "Hoyer_2018",       2250,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "diagmeta",         2150,
                "set_0", "ord_bivariate", TRUE,  10, "GAD_2", "logit_GLMM",       2100,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "Hoyer_2018",       4500,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "diagmeta",         5200,
                "set_0", "ord_bivariate", TRUE,  10, "HADS",  "logit_GLMM",       6100,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "Hoyer_2018",       2250,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "diagmeta",         2350,
                "set_0", "ord_bivariate", TRUE,  50, "BAI",   "logit_GLMM",       3650,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "Hoyer_2018",       500,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "diagmeta",         500,
                "set_0", "ord_bivariate", TRUE,  50, "GAD_2", "logit_GLMM",       650,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "Hoyer_2018",       2350,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "diagmeta",         1600,
                "set_0", "ord_bivariate", TRUE,  50, "HADS",  "logit_GLMM",       3250,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "Hoyer_2018",       3750,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "diagmeta",         3700,
                "set_1", "ord_bivariate", TRUE,  10, "BAI",   "logit_GLMM",       4300,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "Hoyer_2018",       2000,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "diagmeta",         2000,
                "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "logit_GLMM",       2050,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "Hoyer_2018",       4150,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "diagmeta",         3150,
                "set_1", "ord_bivariate", TRUE,  10, "HADS",  "logit_GLMM",       3400,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "Hoyer_2018",       2150,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "diagmeta",         2350,
                "set_1", "ord_bivariate", TRUE,  50, "BAI",   "logit_GLMM",       2900,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "Hoyer_2018",       450,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "diagmeta",         500,
                "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "logit_GLMM",       700,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "Hoyer_2018",       2150,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "diagmeta",         1450,
                "set_1", "ord_bivariate", TRUE,  50, "HADS",  "logit_GLMM",       2400,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "Hoyer_2018",       5900,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "diagmeta",         5700,
                "set_2", "ord_bivariate", TRUE,  10, "BAI",   "logit_GLMM",       6200,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "Hoyer_2018",       1800,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "diagmeta",         1850,
                "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "logit_GLMM",       1850,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "Hoyer_2018",       6300,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "diagmeta",         5300,
                "set_2", "ord_bivariate", TRUE,  10, "HADS",  "logit_GLMM",       5200,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "Hoyer_2018",       2200,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "diagmeta",         2650,
                "set_2", "ord_bivariate", TRUE,  50, "BAI",   "logit_GLMM",       3250,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "Hoyer_2018",       450,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "diagmeta",         450,
                "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "logit_GLMM",       550,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "Hoyer_2018",       2000,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "diagmeta",         1400,
                "set_2", "ord_bivariate", TRUE,  50, "HADS",  "logit_GLMM",       2950),
                dimnames = list(NULL, c("prior_set", "model", "random_thresholds", "n_studies", "test", "DGM", "n_iter"))))
          row_of_this_cell <- which(n_iter_of_each_cell$prior_set == prior_set &
                                    n_iter_of_each_cell$model == model_parameterisation &
                                    n_iter_of_each_cell$random_thresholds == random_thresholds &
                                    n_iter_of_each_cell$n_studies == n_studies &
                                    n_iter_of_each_cell$test == test &
                                    n_iter_of_each_cell$DGM == DGM)
          if (length(row_of_this_cell) == 1 && N_per_study_mean == 500) {

                n_iter_for_4_chains <- n_iter_of_each_cell$n_iter[[row_of_this_cell]]

          } else if (is.null(n_iter_for_4_chains_of_this_run)) {

                stop("no sampling iterations are written for this cell: ", prior_set, ", ", model_parameterisation,
                     ", random_thresholds = ", random_thresholds, ", ", n_studies, " studies, ", test, ", DGM ", DGM)

          }
          ## (a pilot or a real-data fit gives its own number)
          if (!is.null(n_iter_for_4_chains_of_this_run)) {
                n_iter_for_4_chains <- as.integer(n_iter_for_4_chains_of_this_run)
          }
          n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
          ##
          ## (the 4 Set 0 O-biv-RC 10-study cells that use adapt_delta = 0.95 and n_burnin = 1000; every other cell uses 0.80 and 500)
          is_O_biv_RC_set_0_adapt_delta_095_cell <- (prior_set == "set_0" &&
                                                     model_parameterisation == "ord_bivariate" &&
                                                     isTRUE(random_thresholds) &&
                                                     n_studies == 10 &&
                                                     test %in% c("BAI", "HADS") &&
                                                     DGM %in% c("Hoyer_2018", "ord_HSROC_random_C"))
          ## (the same 4 kinds of cell at prior sets 1 and 2, 8 cells: adapt_delta = 0.95 and learning_rate = 0.025 with the usual
          ##  500 burn-in. With 0.80 and 0.05 they lost 14% to 36% of their fits to divergences; on the same seeds this setting
          ##  lost 1 of 64. Every other cell of sets 1 and 2 uses 0.80 and 0.05.)
          is_O_biv_RC_sets_1_and_2_adapt_delta_095_cell <- (prior_set %in% c("set_1", "set_2") &&
                                                            model_parameterisation == "ord_bivariate" &&
                                                            isTRUE(random_thresholds) &&
                                                            n_studies == 10 &&
                                                            test %in% c("BAI", "HADS") &&
                                                            DGM %in% c("Hoyer_2018", "ord_HSROC_random_C"))
          ## (4 more O-biv-RC cells that use adapt_delta = 0.95 and learning_rate = 0.025, and no other cell: BAI with 85% of its
          ##  thresholds missing, 10 studies, DGM Jones and DGM ord_biv_random_C, at prior sets 1 and 2. At set 2, with 0.80 and
          ##  0.05, they lost 22.2% (173 of 780) and 19.2% (136 of 710) of their fits to divergences; on 51 of those seeds fitted
          ##  again with this setting 0 were lost; their fits at 0.80 were taken out of the sim on 16 Nov and the two cells are
          ##  fitted again. The two cells of set 1 lost 2 of 10 and 4 of 10 fits in the pilot and get the setting before they run.)
          is_O_biv_RC_sets_1_and_2_BAI_85_pct_adapt_delta_095_cell <- (prior_set %in% c("set_1", "set_2") &&
                                                                model_parameterisation == "ord_bivariate" &&
                                                                isTRUE(random_thresholds) &&
                                                                n_studies == 10 &&
                                                                test == "BAI" &&
                                                                isTRUE(all.equal(target_proportion_missing_of_this_run, 0.85)) &&
                                                                DGM %in% c("Jones", "ord_biv_random_C"))
          ## (the 2 Set 0 Jones cells of the factorial design that use adapt_delta = 0.95: GAD-2 with 85% of its thresholds missing,
          ##  50 studies, DGM Hoyer_2018 and DGM ord_biv_fixed_C. With 0.80 they lost 11.8% and 10.3% of their fits to divergences;
          ##  with 0.95, 1 of 878 and 1 of 815 fits. Every other Jones cell uses 0.80.)
          is_Jones_set_0_adapt_delta_095_cell <- (prior_set == "set_0" &&
                                                  model_parameterisation == "Jones" &&
                                                  n_studies == 50 &&
                                                  test == "GAD_2" &&
                                                  isTRUE(all.equal(target_proportion_missing_of_this_run, 0.85)) &&
                                                  DGM %in% c("Hoyer_2018", "ord_biv_fixed_C"))
          ## (the cells sampled with NicoStan's auto-chunking, 2 chunks at the burn-in and at the sampling, 1 thread a chain: the
          ##  model's library built with the auto-chunking engine lies beside its Stan file in the cache as
          ##  <name>_nuisance_first_auto_chunking_on_model.so and gives the same log density and gradient as the library without it.
          ##  They are the cells where a gradient was faster with chunks than with the sim library's own build: the stratified
          ##  bivariate with 50 studies at BAI (1.43 times) and HADS (1.09 times), cut by its (threshold, study) pairs; O-biv-FC with
          ##  50 studies at BAI (1.56 times) and HADS (1.08 times), O-biv-RC at BAI with 50 studies (1.41 times, with 8 chunks),
          ##  Jones at BAI with 50 studies (1.18 times) and ord-HSROC at BAI with 50 studies (1.13 times, with 4 chunks), cut by
          ##  their studies. Every other cell is sampled without auto-chunking.)
          is_stratified_bivariate_auto_chunking_cell <- (model_parameterisation == "Bivariate" &&
                                                         !isTRUE(random_thresholds) &&
                                                         n_studies == 50 &&
                                                         test %in% c("BAI", "HADS"))
          ## (O-biv-FC is sampled with auto-chunking at prior set 0 only. At prior set 1, HADS, 50 studies, the workers of 776 of
          ##  844 lane runs died of memory corruption inside the chunked library on 16 Nov, 01:36 to 01:59 (munmap_chunk(): invalid
          ##  pointer; corrupted double-linked list; malloc(): smallbin double linked list corrupted); at set 0 no run ever did.
          ##  Sets 1 and 2 use the library without auto-chunking until the chunked one is proven on their priors.)
          is_O_biv_FC_auto_chunking_cell <- (model_parameterisation == "ord_bivariate" &&
                                             !isTRUE(random_thresholds) &&
                                             n_studies == 50 &&
                                             prior_set == "set_0" &&
                                             test %in% c("BAI", "HADS"))
          ## (O-biv-RC is sampled without auto-chunking in every cell: its fused file of 15 Nov without auto-chunking is 1.39
          ##  times faster an iteration at BAI with 50 studies than the fused file before it with 8 chunks)
          is_O_biv_RC_auto_chunking_cell <- FALSE
          is_Jones_auto_chunking_cell <- (model_parameterisation == "Jones" &&
                                          test == "BAI" && n_studies == 50)
          is_ord_HSROC_auto_chunking_cell <- (model_parameterisation == "ord_HSROC" &&
                                              test == "BAI" && n_studies == 50)
          is_auto_chunking_cell <- is_stratified_bivariate_auto_chunking_cell || is_O_biv_FC_auto_chunking_cell ||
                                   is_O_biv_RC_auto_chunking_cell || is_Jones_auto_chunking_cell ||
                                   is_ord_HSROC_auto_chunking_cell
          ## (the number of chunks of a cell sampled with auto-chunking, at the burn-in and at the sampling)
          N_chunks_of_this_cell <- if (is_O_biv_RC_auto_chunking_cell) 8 else if (is_ord_HSROC_auto_chunking_cell) 4 else 2
          ## (the burn-in iterations: given to this function by the driver, from the settings file of the run: one number for
          ##  every cell, or "the_rule_of_each_cell" = 1000 on the 4 Set 0 O-biv-RC 10-study cells that use adapt_delta = 0.95
          ##  and 500 on every other cell)
          stopifnot(length(n_burnin) == 1, is.numeric(n_burnin) || identical(n_burnin, "the_rule_of_each_cell"))
          if (identical(n_burnin, "the_rule_of_each_cell")) n_burnin <- if (is_O_biv_RC_set_0_adapt_delta_095_cell) 1000L else 500L
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
          ##
          ## ---- THE PRIORS OF EACH MODEL: the ONE place where they are set. One branch per model, every value written out,
          ##      nothing set twice, no environment variable. Not in this block: the Se / Sp prior of prior sets 1 and 2 of the
          ##      ordinal and Jones models, and the Normal on mu of the stratified bivariate (fn_induced_accuracy_priors, below).
          ##
          priors$use_probit_link     <- 1
          priors$K_fold_CV_indicator <- rep(1, n_studies_actual)
          ##
          if (model_parameterisation == "ord_bivariate" && random_thresholds == FALSE) {

                ## ---- O-biv-FC (ordinal bivariate, fixed cutpoints):
                priors$prior_dirichlet_alpha    <- list(rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))  ## Dirichlet(1, ..., 1), each class
                ##
                priors$prior_beta_mu_mean[1, 1] <- 0     ## intercept, class 1
                priors$prior_beta_mu_mean[2, 1] <- 0     ## intercept, class 2
                priors$prior_beta_mu_SD[1, 1]   <- 1.0
                priors$prior_beta_mu_SD[2, 1]   <- 1.0
                ##
                priors$prior_beta_SD_SD         <- c(0.5, 0.5)

          } else if (model_parameterisation == "ord_bivariate" && random_thresholds == TRUE) {

                ## ---- O-biv-RC (ordinal bivariate, random cutpoints): log(kappa) ~ Student-t(dof, mean, SD) in each class
                priors$prior_dirichlet_alpha    <- list(rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))  ## Dirichlet(1, ..., 1), each class
                ##
                priors$prior_beta_mu_mean[1, 1] <- 0     ## intercept, class 1
                priors$prior_beta_mu_mean[2, 1] <- 0     ## intercept, class 2
                priors$prior_beta_mu_SD[1, 1]   <- 1.0
                priors$prior_beta_mu_SD[2, 1]   <- 1.0
                ##
                priors$prior_beta_SD_SD         <- c(0.5, 0.5)
                ##
                priors$prior_kappa_mean         <- c(log(125), log(125))
                priors$prior_kappa_SD           <- c(1.0, 1.0)
                priors$prior_kappa_dof          <- c(1000, 1000)
                ##
                # priors$kappa_lb                 <- 0.001
                ## (initial values)
                for (kk in 1:n_chains) {
                      init_lists_per_chain[[kk]]$kappa <- c(50, 50)
                }

          } else if (model_parameterisation == "ord_HSROC" && random_thresholds == TRUE) {

                ## ---- O-HSROC-RC (ordinal HSROC, random cutpoints): log(kappa) ~ Student-t(dof, mean, SD)
                priors$prior_dirichlet_alpha    <- rep(1.0, n_thr + 1) ## Dirichlet(1, ..., 1)
                ##
                priors$prior_beta_mu_mean[1]    <- 0     ## intercept
                priors$prior_beta_mu_SD[1]      <- 1.0
                priors$prior_beta_SD_SD         <- 0.5
                ##
                ## (log-scale intercept, the same for every test: the ratio of the two classes' latent SDs is exp(2 * log-scale),
                ##  so N(0, 0.125) puts that ratio between 0.61 and 1.63 with 95% prior probability)
                priors$prior_raw_scale_mu_SD    <- 0.125
                priors$prior_raw_scale_SD_SD    <- 0.5
                ##
                priors$prior_kappa_mean         <- log(125)
                priors$prior_kappa_SD           <- 1.0
                priors$prior_kappa_dof          <- 1000
                # priors$kappa_lb                 <- 0.001
                ## (initial values)
                for (kk in 1:n_chains) {
                      init_lists_per_chain[[kk]]$kappa        <- 50
                      init_lists_per_chain[[kk]]$beta_SD      <- 0.25
                      init_lists_per_chain[[kk]]$raw_scale_SD <- 0.25
                }

          } else if (model_parameterisation == "Jones") {

                ## ---- Jones:
                priors$prior_beta_SD_SD         <- c(0.5, 0.5)
                ##
                priors$prior_raw_scale_SD_SD    <- c(0.5, 0.5)
                ##
                ## ---- TO FILL IN (Jones is done LAST, after every other model): the Box-Cox lambda prior, ONE value for every test.
                ##      NA until it is decided. (As the model stands, this lambda prior line is used only when the one-hill Se / Sp
                ##      prior is off: DTA_MA_Jones_Se_Sp_prior.stan, lines 442 to 444.)
                ## (Stan data cannot hold NA: the package's own N(0, 1) of R_fn_set_priors_MA.R is passed; under the one-hill
                ##  Se / Sp prior of prior sets 0, 1 and 2 the model never reaches its lambda prior line, so it is not used)
                priors$prior_boxcox_lambda_mean <- c(0, 0)
                priors$prior_boxcox_lambda_SD   <- c(1, 1)

          } else if (model_parameterisation == "Bivariate") {

                ## ---- Stratified bivariate: its ONE prior is the Normal on mu (probit Se, probit Sp) at each threshold,
                ##      prior_mu_mean and prior_mu_SD, set by fn_induced_accuracy_priors below; nothing else to set here.

          } else {

                stop("no priors are written for this model: ", model_parameterisation, ", random_thresholds = ", random_thresholds)

          }
          ##
          ## ---- The induced-accuracy priors of the prior set (set_0, set_1 or set_2) ------------------------------------------------------
          ##
          stopifnot(prior_set %in% c("set_0", "set_1", "set_2"))
          priors <- fn_induced_accuracy_priors(priors = priors,
                                               model_parameterisation = model_parameterisation,
                                               test = test,
                                               n_thr = n_thr,
                                               prior_set = as.numeric(sub("set_", "", prior_set)))
          ## ("Dirichlet": each ordinal model's own Dirichlet prior (prior_dirichlet_alpha, all ones), no Se / Sp prior term)
          if (prior_set == "set_0" && ordinal_prior_type == "Dirichlet" &&
              model_parameterisation %in% c("ord_bivariate", "ord_HSROC")) {
                priors$use_Se_Sp_prior <- 0
          }
          ## (set 0 never has a Se / Sp prior term in the ordinal models; the stratified bivariate has no such term at all)
          if (prior_set == "set_0" && model_parameterisation %in% c("ord_bivariate", "ord_HSROC")) stopifnot(priors$use_Se_Sp_prior == 0)
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



            ## (the model's name for its nuisance blocks: the Jones model has no custom file, its own file names it)
            name_of_the_model_of_this_fit <- if (is.null(custom_file_name)) "DTA_MA_Jones" else
                                             sub("[.]stan$", "", custom_file_name)
            ## (the nuisance blocks come from MetaOrdDTA itself, never written here: one source of truth)
            nuisance_params <- MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                  name_of_the_model      = name_of_the_model_of_this_fit,
                  model_parameterisation = model_parameterisation,
                  random_thresholds      = random_thresholds)
            ## (nuisance_params_of_this_run: only those blocks as nuisance; the rest go to the main block)
            if (!is.null(nuisance_params_of_this_run)) {
                  nuisance_params <- nuisance_params_of_this_run
            }
            ##
            ##
            interest_only_params <- if (model_parameterisation == "ord_bivariate") {
              
                                              if (random_thresholds) c("beta_mu", "dirichlet_phi")
                                              else c("beta_mu", "C_raw_vec")
              
                                    } else if (model_parameterisation == "ord_HSROC") {
                                      
                                              if (random_thresholds) c("beta_mu", "raw_scale_mu", "dirichlet_phi")
                                              else c("beta_mu", "raw_scale_mu", "C_raw_vec")
                                      
                                    } else {
                                          
                                              NULL
                                      
                                    }
            ##
            # if (any(candidate_route) && is_a_pop_median_model) {
            ## (set 0 only: the pop-median model's Se / Sp prior term off; sets 1 and 2 keep their Se / Sp prior fields)
            if (any(candidate_route) && is_a_pop_median_model && ordinal_prior_type == "Dirichlet" &&
                model_parameterisation == "ord_bivariate" && prior_set == "set_0") {
                priors[names(prior_fields)] <- prior_fields
            }
            ## (sets 1 and 2 always have their Se / Sp prior term in the ordinal and Jones models; the stratified bivariate's
            ##  one prior is the Normal on mu, with its mean and SD at each threshold)
            if (prior_set %in% c("set_1", "set_2") && model_parameterisation != "Bivariate") stopifnot(priors$use_Se_Sp_prior == 1)
            if (model_parameterisation == "Bivariate") stopifnot(!is.null(priors$prior_mu_mean), !is.null(priors$prior_mu_SD))
            time_after_the_model_prep <- Sys.time()
            ## (keep_the_data_file_in = a folder: the Stan data of this cmdstanr fit is saved there as it is written; NULL = not saved)
            if (!is.null(keep_the_data_file_in) && sampling_algorithm == "Stan") {
                  file_of_the_Stan_data <- file.path(keep_the_data_file_in,
                                                     paste0("seed_", seed, "_", file_name_string, ".rds"))
                  trace("write_stan_json", where = asNamespace("cmdstanr"), print = FALSE,
                        tracer = bquote(saveRDS(data, paste0(.(file_of_the_Stan_data),
                                                         if (grepl("init", basename(file))) "_init" else ""))))
            }
            ## (a test of sampler settings: the settings named in this file replace the values written below)
            ## (sampler_settings_to_test_file = the .rds file of a list of settings; NULL = no test, the values written below)
            sampler_settings_to_test <- if (!is.null(sampler_settings_to_test_file))
                  readRDS(sampler_settings_to_test_file) else list()
            ## (a cell sampled with auto-chunking: its settings, which a test file can still replace; and N, the number of rows the
            ##  chunks cut, given to NicoStan's initialise_model in this worker: Nico_model$new does not hand N on, and without N
            ##  NicoStan chunks nothing. N = the number of (threshold, study) pairs with a reported count, sum(n_obs_cutpoints), for
            ##  the stratified bivariate; the number of studies for O-biv-FC, O-biv-RC, Jones and ord-HSROC.)
            if (is_auto_chunking_cell) {
                  options(MetaOrdDTA.rows_of_auto_chunking = if (is_stratified_bivariate_auto_chunking_cell) "pairs" else "studies")
                  sampler_settings_to_test <- utils::modifyList(list(auto_chunking            = TRUE,
                                                                     run_in_fresh_R_process   = FALSE,
                                                                     n_threads_WCP_burnin     = 1,
                                                                     num_chunks_burnin        = N_chunks_of_this_cell,
                                                                     n_threads_WCP_sampling   = 1,
                                                                     num_chunks_sampling      = N_chunks_of_this_cell),
                                                                sampler_settings_to_test)
                  namespace_of_NicoStan <- asNamespace("NicoStan")
                  if (!isTRUE(attr(get("initialise_model", envir = namespace_of_NicoStan), "gives_N_of_auto_chunking"))) {
                        initialise_model_of_NicoStan <- get("initialise_model", envir = namespace_of_NicoStan)
                        initialise_model_giving_N <- function(...) {
                              arguments <- list(...)
                              if (is.null(arguments$N)) {
                                    arguments$N <- if (identical(getOption("MetaOrdDTA.rows_of_auto_chunking"), "pairs"))
                                          as.integer(sum(arguments$Stan_data_list$n_obs_cutpoints)) else
                                          as.integer(arguments$Stan_data_list$n_studies)
                              }
                              do.call(initialise_model_of_NicoStan, arguments)
                        }
                        attr(initialise_model_giving_N, "gives_N_of_auto_chunking") <- TRUE
                        unlockBinding("initialise_model", namespace_of_NicoStan)
                        assign("initialise_model", initialise_model_giving_N, envir = namespace_of_NicoStan)
                        lockBinding("initialise_model", namespace_of_NicoStan)
                  }
            }
            adapt_delta_of_this_cell <- if (is_O_biv_RC_set_0_adapt_delta_095_cell || is_O_biv_RC_sets_1_and_2_adapt_delta_095_cell ||
                                            is_Jones_set_0_adapt_delta_095_cell || is_O_biv_RC_sets_1_and_2_BAI_85_pct_adapt_delta_095_cell) 0.95 else 0.80
            learning_rate_of_this_cell <- if (is_O_biv_RC_sets_1_and_2_adapt_delta_095_cell || is_O_biv_RC_sets_1_and_2_BAI_85_pct_adapt_delta_095_cell) 0.025 else 0.05
            model_samples_obj <- model_prep_obj$sample( seed = MCMC_seed, ## MCMC seed
                                            n_burnin = n_burnin,
                                            n_iter = n_iter,
                                            adapt_delta = adapt_delta_of_this_cell,
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
                                                  ## (the settings of a sampler test, when a file of them is given, replace the values written here)
                                                  NicoStan_args = utils::modifyList(list(
                                            ## (every $sample setting that changes sampling, written out)
                                            adapt_delta                          = adapt_delta_of_this_cell,  ## 0.95 on the 4 Set 0 O-biv-RC 10s cells and the 2 Set 0 Jones cells of GAD-2 with 85% missing; 0.80 everywhere else
                                            ## (the paper's models, sampled with no auto-chunking library but in the cells of
                                            ##  is_auto_chunking_cell, whose settings replace this one; every model,
                                            ##  the Jones model included, runs its fused AVX file)
                                            auto_chunking                        = FALSE,
                                            debug                                = FALSE,
                                            Stan_file_of_the_fused_version       = Stan_file_of_the_fused_version,
                                            ###
                                            # metric_shape_main                    = "dense", ## --------
                                            metric_shape_main                    = "diag",    ## --------
                                            ##
                                            nuisance_params                      = nuisance_params,
                                            sample_nuisance                      = !is.null(nuisance_params),
                                            ## ("none": no nuisance values kept for the summary or the traces)
                                            nuisance_to_track                    = "none",
                                            diffusion_HMC                        = TRUE,
                                            diffusion_HMC_integrator             = "kick_flow_kick",  ## signature
                                            partitioned_HMC                      = FALSE,  ## R_fn_sample.R:1606
                                            ##
                                            burnin_algorithm                     = "ESJD", ## --------
                                            ##
                                            n_adapt                              = n_burnin - round(n_burnin / 10),
                                            clip_iter                            = 70L,
                                            clip_iter_tau                        = 200L,
                                            ##
                                            n_refresh                            = round(n_burnin / 20),
                                            parallel_method                      = "RcppParallel",
                                            ##
                                            n_superchains                        = 4,
                                            ##
                                            learning_rate                        = learning_rate_of_this_cell,  ## 0.025 on the 8 O-biv-RC 10s cells of sets 1 and 2 at 0.95; 0.05 everywhere else
                                            # learning_rate                        = 0.025,
                                            ##
                                            eps_initial                          = NULL,
                                            eps_initial_iter                     = NULL,
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
                                            share_tau_ii_across_chains_in_burnin = TRUE,
                                            force_multi_attempts_for_metric      = FALSE,
                                            ##
                                            tau_initial                          = pi,         ## pi ## --------
                                            tau_adaptation_block                 = "main",  ## --------
                                            interest_only                        = NULL,
                                            ##
                                            tau_sampling_scale                   = "none",  ## R_fn_sample.R:1446
                                            tau_jitter_burnin                    = "uniform",
                                            tau_jitter_sampling_type             = "uniform",  ## same as burn-in
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
                                            metric_update_rescales_eps_and_tau   = TRUE,  ## signature ## ----
                                            metric_update_rescales_tau           = TRUE,
                                            initial_eps_search_at_each_chain     = FALSE,
                                            ## (the sampling runs at the step size handed over by the burn-in)
                                            calibrate_eps_on_sampling_chains_before_sampling = FALSE,
                                            ##
                                            pooled_metric_warm_up                = 16, ## says how many draws it waits for before its first estimate replaces the starting identity metric.
                                            metric_start_iter                    = round(n_burnin/10) ## says WHEN the metric starts collecting draws.
                                            ), sampler_settings_to_test))
            ##SO
            ## ---- Get summaries:
            ##
            time_after_the_sample_call <- Sys.time()
            ns_sum <- if (identical(sampling_algorithm, "NicoStan")) {
                  model_samples_obj$internal_obj$outs_stan_sampling$NicoStan_summary_obj
            } else {
                  NULL
            }
            if (!is.null(ns_sum)) {
                  model_summary_obj <- NULL
                  summary_obj_for_getters <- ns_sum
                  tibble_main <- ns_sum$get_summary_main()
                  tibble_gq   <- ns_sum$get_summary_generated_quantities()
                  tibble_all  <- rbind(tibble_main, tibble_gq)
            } else {
                  model_summary_obj <- model_samples_obj$summary( compute_main_params = TRUE,
                                                                  compute_transformed_parameters = FALSE,
                                                                  compute_generated_quantities = TRUE,
                                                                  save_log_lik_trace = FALSE)
                  summary_obj_for_getters <- model_summary_obj
                  tibble_main <- model_summary_obj$get_summary_main()
                  tibble_gq   <- model_summary_obj$get_summary_generated_quantities()
                  tibble_all  <- rbind(tibble_main, tibble_gq)
            }
            time_after_the_summary <- Sys.time()
            ##
            # tibbles_list_testing[[MCMC_seed]] <- tibble_gq

            if (!is.null(model_summary_obj)) {
              try({
                 model_summary_obj$extract_params(params = "kappa")
                 model_summary_obj$extract_params(params = "C_raw_vec")
                 model_summary_obj$extract_params(params = "beta_mu")
                 model_summary_obj$extract_params(params = "beta_SD")
              })
            }

            summary_obj_for_getters$get_efficiency_metrics()$time_burnin
            summary_obj_for_getters$get_efficiency_metrics()$time_sampling


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
              ESS_vec_Se <- pmin(tibble_Se$n_eff, tibble_Se$n_eff_tail, tibble_Se$n_eff_sd) ## min ESS = min(bulk, tail, SD)
              ESS_vec_Sp <- pmin(tibble_Sp$n_eff, tibble_Sp$n_eff_tail, tibble_Sp$n_eff_sd)
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
              div_of_this_fit <- summary_obj_for_getters$get_divergences()
              n_divergent_per_chain <- div_of_this_fit$n_divs_per_chain
              if (!is.null(n_divergent_per_chain)) {
                n_divergent <- sum(n_divergent_per_chain)
              } else if (!is.null(div_of_this_fit$n_divs) && length(div_of_this_fit$n_divs) == 1L && !is.na(div_of_this_fit$n_divs)) {
                n_divergent <- div_of_this_fit$n_divs
              } else {
                stop("the number of divergent transitions of this fit was not found")
              }
              percent_divergent <- 100 * n_divergent / (n_chains * n_iter)
              ## (keep_the_divergence_trace_in = a folder: the divergence trace and the draws of this NicoStan fit are saved there; NULL = not saved)
              if (!is.null(keep_the_divergence_trace_in) && sampling_algorithm == "NicoStan") {
                    outs_of_this_fit <- model_samples_obj$internal_obj$outs_stan_sampling
                    saveRDS(list(div_trace = outs_of_this_fit$NicoStan_model_obj$result$sampling_object[[2]],
                                 ## (the draws of the nuisance blocks that nuisance_to_track keeps, as NicoStan holds them: one
                                 ##  matrix a chain, one row a nuisance value, on the scale it is sampled on; empty when none is kept)
                                 nuisance_trace = outs_of_this_fit$NicoStan_model_obj$result$sampling_object[[3]],
                                 draws_array = outs_of_this_fit$draws_array),
                            file.path(keep_the_divergence_trace_in,
                                      paste0("seed_", seed, "_", file_name_string, ".rds")))
              }
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
              efficiency_of_this_fit <- tryCatch(summary_obj_for_getters$get_efficiency_metrics(),
                                                 error = function(e) list())
              seconds_burnin         <- c(efficiency_of_this_fit$time_burnin, NA)[1]
              seconds_sampling       <- c(efficiency_of_this_fit$time_sampling, NA)[1]
              HMC_info_of_this_fit   <- tryCatch(summary_obj_for_getters$get_HMC_info(),
                                                 error = function(e) list())
              HMC_info_of_the_sampling <- HMC_info_of_this_fit$stan_HMC_info_list$out_list_sampling
              L_max_per_chain        <- HMC_info_of_the_sampling$L_max_per_chain_sampling
              eps_sampling_val       <- if (!is.null(ns_sum)) c(HMC_info_of_this_fit$eps_main, NA)[1] else c(HMC_info_of_this_fit$eps, NA)[1]
              tau_sampling_val       <- if (!is.null(ns_sum)) c(HMC_info_of_this_fit$tau_main, NA)[1] else c(HMC_info_of_this_fit$tau, NA)[1]
              L_mean_of_this_fit     <- if (!is.null(ns_sum)) tau_sampling_val / eps_sampling_val else c(HMC_info_of_this_fit$L_main_during_sampling, NA)[1]
              n_grad_evals_sampling  <- if (!is.null(ns_sum)) {
                L_mean_of_this_fit * n_iter * n_chains
              } else {
                c(efficiency_of_this_fit$Min_ESS_main /
                  efficiency_of_this_fit$Min_ess_per_grad_main_samp, NA)[1]
              }
              min_ESS_per_1000_grads <- 1000 * min_ESS / n_grad_evals_sampling
              seconds_of_each_step   <- tryCatch(diff(as.numeric(c(time_of_the_start_of_this_fit,
                                                     time_after_the_simulated_data,
                                                     time_after_the_model_prep,
                                                     time_after_the_sample_call,
                                                     time_after_the_summary,
                                                     Sys.time()))),
                                                 error = function(e) rep(NA_real_, 5))
              ##
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
              file_of_the_HMC_log <- file.path(output_dir,
                                               paste0("HMC_of_every_fit_", file_name_string, ".csv"))
              suppressWarnings(utils::write.table(data.frame(seed = seed,
                                                             sampling_algorithm = sampling_algorithm,
                                                             n_burnin = n_burnin,
                                                             eps_sampling = eps_sampling_val,
                                                             tau_sampling = tau_sampling_val,
                                                             L_mean_sampling = L_mean_of_this_fit,
                                                             L_max_sampling = if (length(L_max_per_chain) > 0)
                                                                                max(L_max_per_chain) else NA),
                                                  file = file_of_the_HMC_log,
                                                  sep = ",",
                                                  row.names = FALSE,
                                                  col.names = !file.exists(file_of_the_HMC_log),
                                                  append = TRUE))
              ##
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
            ## (the MCSE of the cell is worked out by one lane at a time: the check reads every results file of the cell and
            ##  works out the metrics of all of them, 16 s for 1,900 fits on one thread, and it ran after every fit in
            ##  every lane. A lane runs it only when no lane has run it for this cell in the last 30 seconds: the time of
            ##  the last check is the time of a file in the cell's folder. The lane whose check finds the target met
            ##  leaves a second file there; a lane that sees that file, written in the last 10 minutes, stops as if its
            ##  own check had found the target met. The check itself, and so the stopping rule, is unchanged.
            ##  The time now is read as the time of a file written now in the same folder, so that the three times are
            ##  those of one clock, the folder's, on every machine.)
            file_of_the_last_MCSE_check_of_this_cell <- file.path(output_dir, paste0("time_of_the_last_MCSE_check_", file_name_string, ".txt"))
            file_of_the_MCSE_target_met_for_every_lane <- file.path(output_dir, paste0("MCSE_target_met_", file_name_string, ".txt"))
            file_written_now <- file.path(output_dir, paste0("time_now_of_", Sys.info()[["nodename"]], "_process_", Sys.getpid(), ".txt"))
            writeLines("", file_written_now)
            time_now_in_seconds <- as.numeric(file.mtime(file_written_now))
            unlink(file_written_now)
            seconds_since_the_target_was_found_met <- if (file.exists(file_of_the_MCSE_target_met_for_every_lane))
                  time_now_in_seconds - as.numeric(file.mtime(file_of_the_MCSE_target_met_for_every_lane)) else Inf
            seconds_since_the_last_MCSE_check <- if (file.exists(file_of_the_last_MCSE_check_of_this_cell))
                  time_now_in_seconds - as.numeric(file.mtime(file_of_the_last_MCSE_check_of_this_cell)) else Inf
            if (stop_early_if_target_met && isTRUE(seconds_since_the_target_was_found_met < 600)) {
                  cat("TARGET MCSE ACHIEVED (the check of another lane found it)! Stopping early.\n")
                  return(list(stopped_early = TRUE,
                              final_n_sims = NA_integer_))
            }
            ## (never two lanes in the check of one cell at once: the lane that runs the check makes this folder, which only
            ##  one lane can make, and removes it when its check ends or when this function ends; a folder left by a lane
            ##  that died is removed once it is 5 minutes old)
            folder_of_the_lane_that_checks_the_MCSE <- file.path(output_dir, paste0("a_lane_is_working_out_the_MCSE_of_", file_name_string))
            if (dir.exists(folder_of_the_lane_that_checks_the_MCSE) &&
                isTRUE(time_now_in_seconds - as.numeric(file.mtime(folder_of_the_lane_that_checks_the_MCSE)) > 300)) {
                  unlink(folder_of_the_lane_that_checks_the_MCSE, recursive = TRUE)
            }
            if ((iii %% check_MCSE_interval == 0 || iii == length(sim_ids_chunk)) && !isTRUE(seconds_since_the_last_MCSE_check < 30) &&
                isTRUE(dir.create(folder_of_the_lane_that_checks_the_MCSE, showWarnings = FALSE))) {

                    ## (this lane removes the folder only while it is the lane that made it)
                    this_lane_made_the_folder_of_the_MCSE_check <- TRUE
                    on.exit(if (isTRUE(this_lane_made_the_folder_of_the_MCSE_check)) unlink(folder_of_the_lane_that_checks_the_MCSE, recursive = TRUE), add = TRUE)
                    writeLines("", file_of_the_last_MCSE_check_of_this_cell)
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
                                                              min_studies_per_threshold = 1)  ## (every threshold that at least one study reports counts; the value the results tables are made with)



                    ## (check_MCSE_status returns MCSE_val once it has min_runs_needed results, and
                    ##  mcse_value = NA before then; reading mcse_value alone stopped every worker at its first
                    ##  full check)
                    ## (the 30 seconds before the next check count from the end of this one)
                    writeLines("", file_of_the_last_MCSE_check_of_this_cell)
                    unlink(folder_of_the_lane_that_checks_the_MCSE, recursive = TRUE)
                    this_lane_made_the_folder_of_the_MCSE_check <- FALSE
                    MCSE_value_of_this_check <- c(current_mcse_status$MCSE_val, current_mcse_status$mcse_value)[1]
                    cat(sprintf("Completed sims: %d, Current MCSE: %.4f%%\n",
                                current_mcse_status$n_sims,
                                MCSE_value_of_this_check * 100))

                    # Check if MCSE is valid AND target is met
                    if (is.finite(MCSE_value_of_this_check) &&
                        current_mcse_status$mcse_met) {
                      cat("✓ TARGET MCSE ACHIEVED! Stopping early.\n")
                      if (stop_early_if_target_met) {
                        ## (the file holds the number of sims and the MCSE of this check: the driver prints them at the
                        ##  end of the cell, so that not every lane runs the full check again)
                        writeLines(paste(current_mcse_status$n_sims, format(MCSE_value_of_this_check, digits = 15)),
                                   file_of_the_MCSE_target_met_for_every_lane)
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





