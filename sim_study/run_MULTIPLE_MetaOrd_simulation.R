
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
                prior_set <- if (is.null(scenario_row$prior_set)) "set_0" else scenario_row$prior_set
                # box_cox <- scenario_row$box_cox
                softplus <- FALSE ## scenario_row$softplus
                random_thresholds <- scenario_row$random_thresholds
                Dirichlet_random_effects_type <- scenario_row$Dirichlet_random_effects_type
                ##
                ## ---- The sampler settings of this cell (SETTINGS_OF_EVERY_CELL = a csv file, one row per cell) ---------------------
                ##
                if (nzchar(Sys.getenv("SETTINGS_OF_EVERY_CELL"))) {

                      settings_of_every_cell <- read.csv(Sys.getenv("SETTINGS_OF_EVERY_CELL"),
                                                         stringsAsFactors = FALSE)
                      row_of_this_cell <- which(settings_of_every_cell$test == test &
                                                settings_of_every_cell$DGM == DGM &
                                                settings_of_every_cell$n_studies == n_studies &
                                                settings_of_every_cell$model == model_parameterisation &
                                                settings_of_every_cell$random_thresholds == random_thresholds &
                                                settings_of_every_cell$prior_set == prior_set)
                      if (length(row_of_this_cell) != 1) {
                            stop("the settings file needs one row for this cell")
                      }
                      ## (a setting of another cell never stays set)
                      Sys.unsetenv(unique(sub("=.*", "", unlist(strsplit(settings_of_every_cell$settings, " ")))))
                      settings_of_this_cell <- strsplit(settings_of_every_cell$settings[row_of_this_cell], " ")[[1]]
                      do.call(Sys.setenv, stats::setNames(as.list(sub("^[^=]*=", "", settings_of_this_cell)),
                                                          sub("=.*", "", settings_of_this_cell)))

                }

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
                      model_parameterisation %in% c("ord_bivariate", "R&G"),
                      isTRUE(random_thresholds), Dirichlet_random_effects_type == "kappa",
                      is.character(file_name_string), length(file_name_string) == 1,
                      !is.na(file_name_string),
                      grepl("popmed", file_name_string, fixed = TRUE) == is_a_pop_median_model)
            stopifnot(scenario_row$custom_file_name[[1]] %in%
                            c("DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan",
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

          if (nzchar(Sys.getenv("SIM_REAL_DATA_DIR"))) {
                real_x <- readRDS(file.path(Sys.getenv("SIM_REAL_DATA_DIR"),
                                           paste0("real_x_", c(GAD_2 = "GAD2", HADS = "HADS", BAI = "BAI")[[test]], ".RDS")))
                x <- list(real_x$nondiseased, real_x$diseased) ## the model takes the non-diseased class first
          }
          n_thr <- ncol(x[[1]]) - 1
          n_studies_actual <- nrow(x[[1]])

          use_custom_file <- FALSE
          custom_file_name <- NULL
          ##
          ## (one file per model: its Se_Sp_prior fused AVX file; use_Se_Sp_prior = 0 is the model's own Dirichlet prior)
          if (model_parameterisation == "Bivariate") {
              use_custom_file  <- TRUE
              custom_file_name <- "DTA_MA_strat_bivariate_Se_Sp_prior_NicoStan_fused_AVX.stan"
          }
          ## (the stratified bivariate is sampled from its fused file with no Se / Sp prior block)
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
          if (model_parameterisation == "R&G") {
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
          ## (STAN_MATH_ONLY_MODEL = 1: the model's own Stan file, Stan::math only, no user header)
          Stan_math_only_model <- Sys.getenv("STAN_MATH_ONLY_MODEL") == "1"
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

          ## (LOG_LIK_INDICATOR = 0: log_lik_study is not computed in generated quantities)
          if (nzchar(Sys.getenv("LOG_LIK_INDICATOR"))) {

                model_prep_obj$advanced_model_options$log_lik_indicator <- as.integer(Sys.getenv("LOG_LIK_INDICATOR"))
                model_prep_obj$internal_obj$outs_data$stan_data_list$log_lik_indicator <-
                      as.integer(Sys.getenv("LOG_LIK_INDICATOR"))

          }

          # Set priors
          priors <- model_prep_obj$priors

          if (prior_set == "set_0") {
                priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
          } else if (prior_set %in% c("set_1", "set_2")) {
                ## (the Jones and R&G priors of prior set 1 / 2 are matched to the ordinal-bivariate reference prior of that
                ##  set by a prior model, Stan sampling; the reference has alpha = (1, ..., 1))
                priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
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

          if (prior_set == "set_0") {

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
              prior_set == "set_0") {

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
              prior_set == "set_0") {

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
              prior_set == "set_0") {

                      n_iter_for_4_chains <- 3500

          }
          ##
          ## ---- The iterations per chain (4 chains) of each cell -------------------------------------------------------------------------
          ##
          ## (min ESS 2500 at 10 studies and 500 at 50 studies, from the 10th percentile of the min ESS
          ##  of the 10 pilot fits of 1000 iterations of each cell; the ordinal models: prior set 0 only)
          {

                      n_iter_of_each_cell <- as.data.frame(matrix(byrow = TRUE, ncol = 6, data = list(
                            ## model,          random_thresholds, n_studies, test, DGM, n_iter
                            "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  4350,
                            "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 4450,
                            "ord_bivariate", FALSE, 10, "BAI",   "R&G_random_C",     6650,
                            "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3800,
                            "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 3700,
                            "ord_bivariate", FALSE, 10, "GAD_2", "R&G_random_C",     3500,
                            "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  4450,
                            "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 4000,
                            "ord_bivariate", FALSE, 10, "HADS",  "R&G_random_C",     3550,
                            "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2350,
                            "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2650,
                            "ord_bivariate", FALSE, 50, "BAI",   "R&G_random_C",     2100,
                            "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  2100,
                            "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2200,
                            "ord_bivariate", FALSE, 50, "GAD_2", "R&G_random_C",     1700,
                            "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  2100,
                            "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2300,
                            "ord_bivariate", FALSE, 50, "HADS",  "R&G_random_C",     1600,
                            "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  5550,
                            "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 2350,
                            "ord_bivariate", TRUE,  10, "BAI",   "R&G_random_C",     2900,
                            "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  2700,
                            "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 2850,
                            "ord_bivariate", TRUE,  10, "GAD_2", "R&G_random_C",     3000,
                            "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  3300,
                            "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 3550,
                            "ord_bivariate", TRUE,  10, "HADS",  "R&G_random_C",     5250,
                            "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  850,
                            "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 800,
                            "ord_bivariate", TRUE,  50, "BAI",   "R&G_random_C",     650,
                            "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  650,
                            "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 2750,
                            "ord_bivariate", TRUE,  50, "GAD_2", "R&G_random_C",     650,
                            "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  850,
                            "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 600,
                            "ord_bivariate", TRUE,  50, "HADS",  "R&G_random_C",     1050,
                            "R&G",           TRUE,  10, "BAI",   "ord_biv_fixed_C",  5250,
                            "R&G",           TRUE,  10, "BAI",   "ord_biv_random_C", 6100,
                            "R&G",           TRUE,  10, "BAI",   "R&G_random_C",     6300,
                            "R&G",           TRUE,  10, "GAD_2", "ord_biv_fixed_C",  6700,
                            "R&G",           TRUE,  10, "GAD_2", "ord_biv_random_C", 4700,
                            "R&G",           TRUE,  10, "GAD_2", "R&G_random_C",     3950,
                            "R&G",           TRUE,  10, "HADS",  "ord_biv_fixed_C",  6200,
                            "R&G",           TRUE,  10, "HADS",  "ord_biv_random_C", 6200,
                            "R&G",           TRUE,  10, "HADS",  "R&G_random_C",     5850,
                            "R&G",           TRUE,  50, "BAI",   "ord_biv_fixed_C",  3250,
                            "R&G",           TRUE,  50, "BAI",   "ord_biv_random_C", 3350,
                            "R&G",           TRUE,  50, "BAI",   "R&G_random_C",     4900,
                            "R&G",           TRUE,  50, "GAD_2", "ord_biv_fixed_C",  1300,
                            "R&G",           TRUE,  50, "GAD_2", "ord_biv_random_C", 1250,
                            "R&G",           TRUE,  50, "GAD_2", "R&G_random_C",     1450,
                            "R&G",           TRUE,  50, "HADS",  "ord_biv_fixed_C",  2400,
                            "R&G",           TRUE,  50, "HADS",  "ord_biv_random_C", 1750,
                            "R&G",           TRUE,  50, "HADS",  "R&G_random_C",     2950,
                            "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  8800,
                            "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 9600,
                            "Bivariate",     FALSE, 10, "BAI",   "R&G_random_C",     7900,
                            "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  6450,
                            "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 3600,
                            "Bivariate",     FALSE, 10, "GAD_2", "R&G_random_C",     4000,
                            "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  7200,
                            "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 8800,
                            "Bivariate",     FALSE, 10, "HADS",  "R&G_random_C",     13600,
                            "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1800,
                            "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1850,
                            "Bivariate",     FALSE, 50, "BAI",   "R&G_random_C",     1400,
                            "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  700,
                            "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 1000,
                            "Bivariate",     FALSE, 50, "GAD_2", "R&G_random_C",     800,
                            "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1550,
                            "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1550,
                            "Bivariate",     FALSE, 50, "HADS",  "R&G_random_C",     1000),
                            dimnames = list(NULL, c("model", "random_thresholds", "n_studies", "test", "DGM",
                                                    "n_iter"))))
                      row_of_this_cell <- which(n_iter_of_each_cell$model == model_parameterisation &
                                                n_iter_of_each_cell$random_thresholds == random_thresholds &
                                                n_iter_of_each_cell$n_studies == n_studies &
                                                n_iter_of_each_cell$test == test &
                                                n_iter_of_each_cell$DGM == DGM)
                      if (length(row_of_this_cell) == 1 &&
                          N_per_study_mean == 500 &&
                          (prior_set == "set_0" || model_parameterisation == "Bivariate")) {

                            n_iter_for_4_chains <- n_iter_of_each_cell$n_iter[[row_of_this_cell]]

                      }

          }
          ##
          ## ---- The iterations per chain (4 chains) of each cell of prior sets 1 and 2 ---------------------------------------------------
          ##
          ## (from the 10 pilot fits of 1000 iterations of each cell that reach the min ESS within 5000 iterations: the
          ##  ord-biv models with the cutpoints relative to the intercept, and the stratified bivariate)
          {

                      n_iter_of_sets_1_and_2 <- as.data.frame(matrix(byrow = TRUE, ncol = 7, data = list(
                            ## prior_set, model, random_thresholds, n_studies, test, DGM, n_iter
                            "set_1", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  3200,
                            "set_1", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 3900,
                            "set_1", "ord_bivariate", FALSE, 10, "BAI",   "R&G_random_C",     4150,
                            "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3850,
                            "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 4800,
                            "set_1", "ord_bivariate", FALSE, 10, "GAD_2", "R&G_random_C",     4300,
                            "set_1", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  3400,
                            "set_1", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 3300,
                            "set_1", "ord_bivariate", FALSE, 10, "HADS",  "R&G_random_C",     3150,
                            "set_1", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2550,
                            "set_1", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2500,
                            "set_1", "ord_bivariate", FALSE, 50, "BAI",   "R&G_random_C",     1800,
                            "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  2350,
                            "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2800,
                            "set_1", "ord_bivariate", FALSE, 50, "GAD_2", "R&G_random_C",     2050,
                            "set_1", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  2050,
                            "set_1", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2350,
                            "set_1", "ord_bivariate", FALSE, 50, "HADS",  "R&G_random_C",     1550,
                            "set_1", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  4450,
                            "set_1", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 5000,
                            "set_1", "ord_bivariate", TRUE,  10, "BAI",   "R&G_random_C",     5000,
                            "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  4850,
                            "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 5000,
                            "set_1", "ord_bivariate", TRUE,  10, "GAD_2", "R&G_random_C",     3650,
                            "set_1", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  4500,
                            "set_1", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 4350,
                            "set_1", "ord_bivariate", TRUE,  10, "HADS",  "R&G_random_C",     4750,
                            "set_1", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  3600,
                            "set_1", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 4950,
                            "set_1", "ord_bivariate", TRUE,  50, "BAI",   "R&G_random_C",     4550,
                            "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  4550,
                            "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 3950,
                            "set_1", "ord_bivariate", TRUE,  50, "GAD_2", "R&G_random_C",     1150,
                            "set_1", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  1850,
                            "set_1", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 3000,
                            "set_1", "ord_bivariate", TRUE,  50, "HADS",  "R&G_random_C",     2650,
                            "set_1", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  3400,
                            "set_1", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 4550,
                            "set_1", "Bivariate",     FALSE, 10, "BAI",   "R&G_random_C",     4400,
                            "set_1", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3250,
                            "set_1", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 2550,
                            "set_1", "Bivariate",     FALSE, 10, "GAD_2", "R&G_random_C",     3400,
                            "set_1", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  3950,
                            "set_1", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 3650,
                            "set_1", "Bivariate",     FALSE, 10, "HADS",  "R&G_random_C",     4400,
                            "set_1", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1700,
                            "set_1", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1800,
                            "set_1", "Bivariate",     FALSE, 50, "BAI",   "R&G_random_C",     1150,
                            "set_1", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  900,
                            "set_1", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 1000,
                            "set_1", "Bivariate",     FALSE, 50, "GAD_2", "R&G_random_C",     700,
                            "set_1", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1500,
                            "set_1", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1300,
                            "set_1", "Bivariate",     FALSE, 50, "HADS",  "R&G_random_C",     950,
                            "set_2", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_fixed_C",  3550,
                            "set_2", "ord_bivariate", FALSE, 10, "BAI",   "ord_biv_random_C", 3600,
                            "set_2", "ord_bivariate", FALSE, 10, "BAI",   "R&G_random_C",     3450,
                            "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_fixed_C",  4550,
                            "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "ord_biv_random_C", 4550,
                            "set_2", "ord_bivariate", FALSE, 10, "GAD_2", "R&G_random_C",     4100,
                            "set_2", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_fixed_C",  3850,
                            "set_2", "ord_bivariate", FALSE, 10, "HADS",  "ord_biv_random_C", 3800,
                            "set_2", "ord_bivariate", FALSE, 10, "HADS",  "R&G_random_C",     3800,
                            "set_2", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_fixed_C",  2450,
                            "set_2", "ord_bivariate", FALSE, 50, "BAI",   "ord_biv_random_C", 2400,
                            "set_2", "ord_bivariate", FALSE, 50, "BAI",   "R&G_random_C",     1900,
                            "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_fixed_C",  2400,
                            "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "ord_biv_random_C", 2700,
                            "set_2", "ord_bivariate", FALSE, 50, "GAD_2", "R&G_random_C",     2100,
                            "set_2", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_fixed_C",  2200,
                            "set_2", "ord_bivariate", FALSE, 50, "HADS",  "ord_biv_random_C", 2650,
                            "set_2", "ord_bivariate", FALSE, 50, "HADS",  "R&G_random_C",     1850,
                            "set_2", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_fixed_C",  4700,
                            "set_2", "ord_bivariate", TRUE,  10, "BAI",   "ord_biv_random_C", 5000,
                            "set_2", "ord_bivariate", TRUE,  10, "BAI",   "R&G_random_C",     5000,
                            "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_fixed_C",  4800,
                            "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "ord_biv_random_C", 5000,
                            "set_2", "ord_bivariate", TRUE,  10, "GAD_2", "R&G_random_C",     4550,
                            "set_2", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_fixed_C",  4350,
                            "set_2", "ord_bivariate", TRUE,  10, "HADS",  "ord_biv_random_C", 5000,
                            "set_2", "ord_bivariate", TRUE,  10, "HADS",  "R&G_random_C",     4700,
                            "set_2", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_fixed_C",  3000,
                            "set_2", "ord_bivariate", TRUE,  50, "BAI",   "ord_biv_random_C", 4850,
                            "set_2", "ord_bivariate", TRUE,  50, "BAI",   "R&G_random_C",     3800,
                            "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_fixed_C",  4500,
                            "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "ord_biv_random_C", 4250,
                            "set_2", "ord_bivariate", TRUE,  50, "GAD_2", "R&G_random_C",     1100,
                            "set_2", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_fixed_C",  2100,
                            "set_2", "ord_bivariate", TRUE,  50, "HADS",  "ord_biv_random_C", 3850,
                            "set_2", "ord_bivariate", TRUE,  50, "HADS",  "R&G_random_C",     4600,
                            "set_2", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_fixed_C",  3950,
                            "set_2", "Bivariate",     FALSE, 10, "BAI",   "ord_biv_random_C", 4700,
                            "set_2", "Bivariate",     FALSE, 10, "BAI",   "R&G_random_C",     4050,
                            "set_2", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_fixed_C",  3650,
                            "set_2", "Bivariate",     FALSE, 10, "GAD_2", "ord_biv_random_C", 3250,
                            "set_2", "Bivariate",     FALSE, 10, "GAD_2", "R&G_random_C",     3600,
                            "set_2", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_fixed_C",  4200,
                            "set_2", "Bivariate",     FALSE, 10, "HADS",  "ord_biv_random_C", 4850,
                            "set_2", "Bivariate",     FALSE, 10, "HADS",  "R&G_random_C",     4750,
                            "set_2", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_fixed_C",  1850,
                            "set_2", "Bivariate",     FALSE, 50, "BAI",   "ord_biv_random_C", 1900,
                            "set_2", "Bivariate",     FALSE, 50, "BAI",   "R&G_random_C",     1200,
                            "set_2", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_fixed_C",  800,
                            "set_2", "Bivariate",     FALSE, 50, "GAD_2", "ord_biv_random_C", 950,
                            "set_2", "Bivariate",     FALSE, 50, "GAD_2", "R&G_random_C",     700,
                            "set_2", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_fixed_C",  1550,
                            "set_2", "Bivariate",     FALSE, 50, "HADS",  "ord_biv_random_C", 1400,
                            "set_2", "Bivariate",     FALSE, 50, "HADS",  "R&G_random_C",     950),
                            dimnames = list(NULL, c("prior_set", "model", "random_thresholds", "n_studies", "test", "DGM",
                                                    "n_iter"))))
                      row_of_this_cell <- which(n_iter_of_sets_1_and_2$prior_set == prior_set &
                                                n_iter_of_sets_1_and_2$model == model_parameterisation &
                                                n_iter_of_sets_1_and_2$random_thresholds == random_thresholds &
                                                n_iter_of_sets_1_and_2$n_studies == n_studies &
                                                n_iter_of_sets_1_and_2$test == test &
                                                n_iter_of_sets_1_and_2$DGM == DGM)
                      if (length(row_of_this_cell) == 1 &&
                          N_per_study_mean == 500) {

                            n_iter_for_4_chains <- n_iter_of_sets_1_and_2$n_iter[[row_of_this_cell]]

                      }

          }
          ##
          if (nzchar(Sys.getenv("N_ITER_FOR_4_CHAINS"))) {
                n_iter_for_4_chains <- as.integer(Sys.getenv("N_ITER_FOR_4_CHAINS"))
          }
          n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
          ##
          ## (NicoStan's burn-in is 125-250 iterations; Stan's is never fewer than 200)
          n_burnin <- if (nzchar(Sys.getenv("N_BURNIN"))) as.integer(Sys.getenv("N_BURNIN")) else
                      if (sampling_algorithm == "NicoStan") 500 else 500
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
                        if (Sys.getenv("KAPPA_PRIOR") == "current") {
                              priors$prior_kappa_mean <- c(log(300), log(300))
                              priors$prior_kappa_SD <- c(0.5, 0.5)
                              priors$prior_kappa_dof <- c(1000, 1000)
                        } else if (Sys.getenv("KAPPA_PRIOR") %in% c("", "wide")) {
                              priors$prior_kappa_mean <- c(log(150), log(150))
                              priors$prior_kappa_SD <- c(1.25, 1.25)
                              priors$prior_kappa_dof <- c(1000, 1000)
                        }
                        ##
                        for (kk in 1:n_chains) {
                          init_lists_per_chain[[kk]]$kappa <- c(50, 50)
                        }
                    } else {
                        priors$prior_beta_SD_SD <- c(0.5, 0.5)
                        priors$prior_kappa_mean <- log(100)
                        priors$prior_kappa_SD <- 1.0
                        priors$prior_kappa_dof <- 5
                        ## (KAPPA_PRIOR=current: the tighter kappa prior, for the sampling check of the pop-median model)
                        if (Sys.getenv("KAPPA_PRIOR") == "current") {
                              priors$prior_kappa_mean <- log(300)
                              priors$prior_kappa_SD <- 0.5
                              priors$prior_kappa_dof <- 1000
                        } else if (Sys.getenv("KAPPA_PRIOR") %in% c("", "wide")) {
                              priors$prior_kappa_mean <- log(150)
                              priors$prior_kappa_SD <- 1.25
                              priors$prior_kappa_dof <- 1000
                        }
                        ##
                        for (kk in 1:n_chains) {
                          init_lists_per_chain[[kk]]$kappa <- 50
                          init_lists_per_chain[[kk]]$beta_SD <- 0.25
                          init_lists_per_chain[[kk]]$raw_scale_SD <- 0.25
                        }
                    }
              
          }
          ##
          if (model_parameterisation == "ord_bivariate") {
            
                    priors$prior_beta_SD_SD <- c(0.5, 0.5)
                    ##
                    if (prior_set %in% c("set_1", "set_2")) {
                          priors$prior_beta_mu_mean[1, 1] <-  0
                          priors$prior_beta_mu_mean[2, 1] <-  0
                          priors$prior_beta_mu_SD[1, 1] <-  c(HADS = 0.20, GAD_2 = 0.25, BAI = 0.35)[[test]]
                          priors$prior_beta_mu_SD[2, 1] <-  c(HADS = 0.20, GAD_2 = 0.25, BAI = 0.35)[[test]]
                    }
                
          } else if (model_parameterisation == "R&G") {
                
                    priors$prior_beta_SD_SD      <- 0.5
                    priors$prior_raw_scale_SD_SD <- 0.5
                    ##
                    ## (alpha = (1, ..., 1): the SD of the intercept prior and of the log-scale prior matched to the
                    ##  ord_biv fixed-cutpoint prior by NicoStan prior-only fits of the pop-median model, 5 seeds, 0.05 rule;
                    ##  MATCHED_PRIORS=0 keeps the runner's own SDs, for the sampling checks)
                    if (Sys.getenv("MATCHED_PRIORS") != "0") {
                          SD_of_the_intercept_matched <-  c(HADS = 1.0, GAD_2 = 1.0, BAI = 1.0)[[test]]
                          SD_of_the_log_scale_matched <-  c(HADS = 0.25, GAD_2 = 0.25, BAI = 0.25)[[test]]
                          if (Sys.getenv("KAPPA_PRIOR") %in% c("", "wide")) SD_of_the_log_scale_matched <-  c(HADS = 0.125, GAD_2 = 0.1875, BAI = 0.1875)[[test]]
                          if (!is.null(SD_of_the_intercept_matched)) priors$prior_beta_mu_SD[1] <-  SD_of_the_intercept_matched
                          if (!is.null(SD_of_the_log_scale_matched)) priors$prior_raw_scale_mu_SD[1] <-  SD_of_the_log_scale_matched
                    }
                
          } else if (model_parameterisation == "Jones") {

                    priors$prior_beta_SD_SD      <- c(0.5, 0.5)
                    priors$prior_raw_scale_SD_SD <- c(0.5, 0.5)
                    ## (the Box-Cox lambda prior of each test, as in the paper)
                    priors$prior_boxcox_lambda_mean <- rep(c(HADS = 0.75, GAD_2 = 1.0, BAI = 0.5)[[test]], 2)
                    priors$prior_boxcox_lambda_SD   <- rep(c(HADS = 0.1, GAD_2 = 0.1, BAI = 0.05)[[test]], 2)

          } else if (model_parameterisation == "Bivariate") {
                
                    ##
                    ## (its Se / Sp priors come from fn_induced_accuracy_priors below)
            
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
              model_parameterisation %in% c("ord_bivariate", "R&G")) {
                priors$use_Se_Sp_prior <- 0
          }
          ## (set 0 never has a Se / Sp prior term, in any model)
          if (prior_set == "set_0" && model_parameterisation != "Jones") stopifnot(priors$use_Se_Sp_prior == 0)
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
            ## (NUISANCE_PARAMS = "beta_z,raw_scale_z" holds only those blocks as nuisance; the rest go to the main block)
            if (nzchar(Sys.getenv("NUISANCE_PARAMS"))) {
                  nuisance_params <- strsplit(Sys.getenv("NUISANCE_PARAMS"), ",")[[1]]
            }
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
            ## (set 0 only: the pop-median model's Se / Sp prior term off; sets 1 and 2 keep their Se / Sp prior fields)
            if (any(candidate_route) && is_a_pop_median_model && ordinal_prior_type == "Dirichlet" &&
                model_parameterisation == "ord_bivariate" && prior_set == "set_0") {
                priors[names(prior_fields)] <- prior_fields
            }
            ## (sets 1 and 2 always have their Se / Sp prior term, in every model)
            if (prior_set %in% c("set_1", "set_2")) stopifnot(priors$use_Se_Sp_prior == 1)
            time_after_the_model_prep <- Sys.time()
            ## (KEEP_THE_DATA_FILE_IN = a folder: the Stan data of this cmdstanr fit is saved there as it is written)
            if (nzchar(Sys.getenv("KEEP_THE_DATA_FILE_IN")) && sampling_algorithm == "Stan") {
                  file_of_the_Stan_data <- file.path(Sys.getenv("KEEP_THE_DATA_FILE_IN"),
                                                     paste0("seed_", seed, "_", file_name_string, ".rds"))
                  trace("write_stan_json", where = asNamespace("cmdstanr"), print = FALSE,
                        tracer = bquote(saveRDS(data, paste0(.(file_of_the_Stan_data),
                                                         if (grepl("init", basename(file))) "_init" else ""))))
            }
            ## (a test of sampler settings: the settings named in this file replace the values written below)
            sampler_settings_to_test <- if (nzchar(Sys.getenv("SAMPLER_SETTINGS_TO_TEST")))
                  readRDS(Sys.getenv("SAMPLER_SETTINGS_TO_TEST")) else list()
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
                                                  NicoStan_args = utils::modifyList(c(list(
                                            ## (every $sample setting that changes sampling, written out)
                                            adapt_delta                          = 0.80,  ## the PS7 value
                                            ## (the paper's models, sampled with no auto-chunking library; every model,
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
                                            learning_rate                        = 0.0625,
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
                                            share_tau_ii_across_chains_in_burnin = FALSE,
                                            force_multi_attempts_for_metric      = FALSE,
                                            ##
                                            tau_initial                          = "adaptive", ## pi ## --------
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
                                            metric_update_rescales_eps_and_tau   = FALSE, ## signature ## ----
                                            metric_update_rescales_tau           = TRUE,
                                            initial_eps_search_at_each_chain     = FALSE,
                                            ##
                                            pooled_metric_warm_up                = 8, ## says how many draws it waits for before its first estimate replaces the starting identity metric.
                                            metric_start_iter                    = round(n_burnin/10) ## says WHEN the metric starts collecting draws.
                                            ),
                                            ## (METRIC_POOLED_WINDOW_RESETS = none, or iterations split by commas)
                                            if (nzchar(Sys.getenv("METRIC_POOLED_WINDOW_RESETS")))
                                                  list(metric_pooled_window_resets =
                                                             if (Sys.getenv("METRIC_POOLED_WINDOW_RESETS") == "none") "none" else
                                                                   as.integer(strsplit(Sys.getenv("METRIC_POOLED_WINDOW_RESETS"), ",")[[1]])),
                                            ## (METRIC_ADAPTATION_END_ITER = the last metric update)
                                            if (nzchar(Sys.getenv("METRIC_ADAPTATION_END_ITER")))
                                                  list(metric_adaptation_end_iter = as.integer(Sys.getenv("METRIC_ADAPTATION_END_ITER"))),
                                            ## (TAU_SAMPLING_SCALE_WITH_A_JITTERED_BURNIN = 1: TAU_SAMPLING_SCALE also applies with a jittered burn-in tau)
                                            if (Sys.getenv("TAU_SAMPLING_SCALE_WITH_A_JITTERED_BURNIN") == "1")
                                                  list(tau_sampling_scale_with_a_jittered_burnin = TRUE),
                                            ## (RESTART_STUCK_CHAINS = 1: a stuck chain restarts from a chain that moves)
                                            if (Sys.getenv("RESTART_STUCK_CHAINS") == "1")
                                                  list(restart_stuck_chains_in_burnin = TRUE,
                                                       n_iter_of_a_stuck_chain = if (nzchar(Sys.getenv("N_ITER_OF_A_STUCK_CHAIN")))
                                                             as.integer(Sys.getenv("N_ITER_OF_A_STUCK_CHAIN")) else 30L),
                                            ## (TAU_WITHIN_THE_CAP_OF_STEPS = 1: tau is never longer than max_L * eps in burn-in)
                                            if (Sys.getenv("TAU_WITHIN_THE_CAP_OF_STEPS") == "1")
                                                  list(hold_tau_within_the_cap_of_steps_in_burnin = TRUE),
                                            ## (EPS_SEARCH_AT_HANDOVER_FROM_EACH_CHAIN = 1: the median of the searches from each chain)
                                            if (Sys.getenv("EPS_SEARCH_AT_HANDOVER_FROM_EACH_CHAIN") == "1")
                                                  list(eps_search_at_handover_from_each_chain = TRUE),
                                            ## (FREEZE_TAU_WHILE_THE_CAP_BINDS = 1: tau is held while it is longer than max_L * eps)
                                            if (Sys.getenv("FREEZE_TAU_WHILE_THE_CAP_BINDS") == "1")
                                                  list(freeze_tau_while_the_cap_of_steps_binds_in_burnin = TRUE),
                                            ## (EPS_SEARCH_WHILE_THE_ACCEPTANCE_IS_HIGH = 1: eps goes up to the search from each chain,
                                            ##  from FIRST_ITER_OF_THE_EPS_SEARCHES on, default 130)
                                            if (Sys.getenv("EPS_SEARCH_WHILE_THE_ACCEPTANCE_IS_HIGH") == "1")
                                                  list(eps_search_while_the_acceptance_is_high_in_burnin = TRUE,
                                                       first_iter_of_the_eps_searches = if (nzchar(Sys.getenv("FIRST_ITER_OF_THE_EPS_SEARCHES")))
                                                             as.integer(Sys.getenv("FIRST_ITER_OF_THE_EPS_SEARCHES")) else 130L,
                                                       n_iter_between_the_eps_searches = if (nzchar(Sys.getenv("N_ITER_BETWEEN_THE_EPS_SEARCHES")))
                                                             as.integer(Sys.getenv("N_ITER_BETWEEN_THE_EPS_SEARCHES")) else 10L),
                                            ## (EPS_DOUBLING_WHILE_THE_ACCEPTANCE_IS_HIGH = 1: eps is doubled while the acceptance is 0.95 or more,
                                            ##  from FIRST_ITER_OF_THE_EPS_DOUBLINGS on, default 60, every N_ITER_BETWEEN_THE_EPS_DOUBLINGS, default 5)
                                            if (Sys.getenv("EPS_DOUBLING_WHILE_THE_ACCEPTANCE_IS_HIGH") == "1")
                                                  list(eps_doubling_while_the_acceptance_is_high_in_burnin = TRUE,
                                                       first_iter_of_the_eps_doublings = if (nzchar(Sys.getenv("FIRST_ITER_OF_THE_EPS_DOUBLINGS")))
                                                             as.integer(Sys.getenv("FIRST_ITER_OF_THE_EPS_DOUBLINGS")) else 60L,
                                                       n_iter_between_the_eps_doublings = if (nzchar(Sys.getenv("N_ITER_BETWEEN_THE_EPS_DOUBLINGS")))
                                                             as.integer(Sys.getenv("N_ITER_BETWEEN_THE_EPS_DOUBLINGS")) else 5L)),
                                                                                    sampler_settings_to_test))
            ##SO
            ## ---- Get summaries:
            ##
            time_after_the_sample_call <- Sys.time()
            model_summary_obj <- model_samples_obj$summary( compute_main_params = TRUE,
                                                            compute_transformed_parameters = FALSE,
                                                            compute_generated_quantities = TRUE,
                                                            save_log_lik_trace = FALSE)
            time_after_the_summary <- Sys.time()
            ##
            tibble_main <- model_summary_obj$get_summary_main()
            tibble_gq   <- model_summary_obj$get_summary_generated_quantities()
            tibble_all  <- rbind(tibble_main, tibble_gq)
            if (nzchar(Sys.getenv("DUMP_SUMMARY"))) {
                  saveRDS(tibble_all, file.path(output_dir, paste0("summary_of_every_parameter_of_seed_", seed, "_",
                                                                   file_name_string, ".rds")))
            }
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
              n_divergent_per_chain <- model_summary_obj$get_divergences()$n_divs_per_chain
              if (is.null(n_divergent_per_chain)) {
                stop("the number of divergent transitions of this fit was not found")
              }
              n_divergent <- sum(n_divergent_per_chain)
              percent_divergent <- 100 * n_divergent / (n_chains * n_iter)
              ## (KEEP_THE_DIVERGENCE_TRACE_IN = a folder: the divergence trace and the draws of this NicoStan fit are saved there)
              if (nzchar(Sys.getenv("KEEP_THE_DIVERGENCE_TRACE_IN")) && sampling_algorithm == "NicoStan") {
                    outs_of_this_fit <- model_samples_obj$internal_obj$outs_stan_sampling
                    saveRDS(list(div_trace = outs_of_this_fit$NicoStan_model_obj$result$sampling_object[[2]],
                                 draws_array = outs_of_this_fit$draws_array),
                            file.path(Sys.getenv("KEEP_THE_DIVERGENCE_TRACE_IN"),
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
              if (Sys.getenv("DUMP_BURNIN_TRACE") == "1") {
                    NicoStan_result_of_this_fit <- tryCatch(
                          model_samples_obj$internal_obj$outs_stan_sampling$NicoStan_model_obj$result,
                          error = function(e) NULL)
                    if (!is.null(NicoStan_result_of_this_fit)) {
                          burnin_of_this_fit <- NicoStan_result_of_this_fit$burnin_object
                          saveRDS(list(seed = seed,
                                       n_burnin = n_burnin,
                                       burnin_tau = burnin_of_this_fit$burnin_tau_main_ii_all_chains,
                                       burnin_eps = burnin_of_this_fit$burnin_eps_main_all_chains,
                                       burnin_L = burnin_of_this_fit$burnin_L_main_ii_all_chains,
                                       burnin_div = burnin_of_this_fit$burnin_div_main_all_chains,
                                       burnin_p_jump = burnin_of_this_fit$burnin_p_jump_main_all_chains,
                                       burnin_tau_update_records = burnin_of_this_fit$burnin_tau_update_records,
                                       burnin_schedule = burnin_of_this_fit$burnin_schedule,
                                       burnin_start_main = burnin_of_this_fit$burnin_start_main_all_chains,
                                       restarts_of_stuck_chains = burnin_of_this_fit$restarts_of_stuck_chains,
                                       theta_main_end_of_burnin = burnin_of_this_fit$theta_main_vectors_all_chains_input_from_R,
                                       theta_us_end_of_burnin = burnin_of_this_fit$theta_us_vectors_all_chains_input_from_R,
                                       eps_calibration = NicoStan_result_of_this_fit$eps_calibration_on_sampling_chains,
                                       HMC_info = HMC_info_of_this_fit),
                                  file.path(output_dir, paste0("burnin_trace_of_seed_", seed, "_", file_name_string, ".rds")))
                    }
              }
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





