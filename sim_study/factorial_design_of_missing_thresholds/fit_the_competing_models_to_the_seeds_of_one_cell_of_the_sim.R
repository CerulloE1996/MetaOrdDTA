##
## ---- The competing models (Hoyer et al. 2018; the logit GLMM) fitted to seeds of one cell of the sim study ----------------------
##
## Each seed's data set is made exactly as the sim runner makes it (run_MULTIPLE_MetaOrd_simulation.R, the lines from
## set.seed(seed) to "Select data"): the same three package files are sourced with the same values, so the competing models
## are given, seed for seed, the data set the MetaOrdDTA models are given. The data set of a seed is made once and both
## models are fitted to it.
## Each fit is returned as a results list with the fields the sim's metrics read (thr_sim_study_helper_fns.R,
## compute_simulation_metrics): x, n_studies, tibble_all, the true Se / Sp and the Se / Sp point estimates and interval
## limits, all in percent.
## The competing models are competitors of the sim only: each is fitted by its own package and no prior of any model is
## taken from them.
##
## This file needs, sourced before it: fit_the_Hoyer_2018_time_to_event_model_to_one_data_set_of_the_sim.R and
## fit_the_logit_GLMM_to_one_data_set_of_the_sim.R; and the package MetaOrdDTA loaded (the runner's require(MetaOrdDTA)),
## whose apply_missingness_pattern() the data step calls.
##

##
## ---- The data set of one seed of a cell, as the sim runner makes it ---------------------------------------------------------------
##
simulated_data_set_of_one_seed_of_a_cell <- function(local_pkg_dir,
                                                     DGM,
                                                     test,
                                                     n_studies,
                                                     N_per_study_mean,
                                                     reduction_factor,
                                                     seed,
                                                     target_proportion_missing_of_this_run) {

      ## (target_proportion_missing_of_this_run: NULL = the test gets its real reporting pattern of missing thresholds, or the
      ##  proportion of the thresholds to make missing, e.g. 0.85; the data script of the package reads it, as in the runner)
      ## (the runner's values, its lines 78 to 86)
      missing_thr <-  TRUE
      N_per_study_SD <-  N_per_study_mean
      assume_perfect_GS <-  TRUE
      true_Mean_prev <-  0.20
      true_SD_probit_prev <-  0.25
      ## (the runner's lines 138 to 162)
      set.seed(seed)
      source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA_v2.R"), local = TRUE)
      source(file.path(local_pkg_dir, "inst", "examples", "MA_missing_thr_prep_data.R"), local = TRUE)
      source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"), local = TRUE)
      if (test == "GAD_2") {
            x <-  x_GAD2_w_missing_thr
      } else if (test == "HADS") {
            x <-  x_HADS_w_missing_thr
      } else if (test == "BAI") {
            x <-  x_BAI_w_missing_thr
      } else {
            stop(paste("Unknown test type:", test))
      }
      return(x)

}

##
## ---- The results list of one fit, in the sim's form, from the fit's table of Se and Sp -------------------------------------------
##
results_list_of_the_fit_of_a_competing_model <- function(name_of_the_model,
                                                         Se_and_Sp,
                                                         x,
                                                         seed,
                                                         n_studies,
                                                         true_DGM_Se,
                                                         true_DGM_Sp) {

      n_thr <-  nrow(Se_and_Sp)
      ## (these fits have no draws: the columns the metrics read for the min ESS are infinite, so no fit is dropped for a low
      ##  ESS)
      tibble_all <-  dplyr::tibble(parameter   = c(paste0("Se_baseline[", seq_len(n_thr), "]"), paste0("Sp_baseline[", seq_len(n_thr), "]")),
                                   mean        = c(Se_and_Sp$Sens, Se_and_Sp$Spec),
                                   sd          = NA_real_,
                                   `2.5%`      = c(Se_and_Sp$SensCI_Lower, Se_and_Sp$SpecCI_Lower),
                                   `50%`       = c(Se_and_Sp$Sens, Se_and_Sp$Spec),
                                   `97.5%`     = c(Se_and_Sp$SensCI_Upper, Se_and_Sp$SpecCI_Upper),
                                   n_eff       = Inf,
                                   Rhat        = NA_real_,
                                   n_Rhat      = NA_real_,
                                   n_eff_tail  = Inf,
                                   n_eff_sd    = Inf)
      return(list(seed                    = seed,
                  percent_divergent       = 0,
                  x                       = x,
                  n_studies               = n_studies,
                  model_parameterisation  = name_of_the_model,
                  softplus                = FALSE,
                  random_thresholds       = FALSE,
                  tibble_all              = tibble_all,
                  true_DGM_Se             = true_DGM_Se,
                  true_DGM_Sp             = true_DGM_Sp,
                  Se_median_vec           = 100 * Se_and_Sp$Sens,
                  Sp_median_vec           = 100 * Se_and_Sp$Spec,
                  Se_mean_vec             = 100 * Se_and_Sp$Sens,
                  Sp_mean_vec             = 100 * Se_and_Sp$Spec,
                  Se_lower_vec            = 100 * Se_and_Sp$SensCI_Lower,
                  Sp_lower_vec            = 100 * Se_and_Sp$SpecCI_Lower,
                  Se_upper_vec            = 100 * Se_and_Sp$SensCI_Upper,
                  Sp_upper_vec            = 100 * Se_and_Sp$SpecCI_Upper))

}

##
## ---- One seed: its data set, then both competing models fitted to it --------------------------------------------------------------
##
## (returns, for each model, its results list (NULL when the fit stopped with an error) and one row for the gate table of the
##  cell: the fit's code (0 = the fit ended without a failure code), the number of its warning messages, whether the fitted
##  random-effects covariance is singular (logit GLMM), its seconds and its error message)
fits_of_the_competing_models_to_one_seed <- function(local_pkg_dir,
                                                     DGM,
                                                     test,
                                                     n_studies,
                                                     N_per_study_mean,
                                                     reduction_factor,
                                                     seed,
                                                     true_DGM_Se,
                                                     true_DGM_Sp,
                                                     distribution_of_the_Hoyer_model,
                                                     smallest_of_the_Hoyer_model,
                                                     largest_of_the_Hoyer_model,
                                                     confidence_level,
                                                     target_proportion_missing_of_this_run) {

      x <-  simulated_data_set_of_one_seed_of_a_cell(local_pkg_dir       = local_pkg_dir,
                                                     DGM                 = DGM,
                                                     test                = test,
                                                     n_studies           = n_studies,
                                                     N_per_study_mean    = N_per_study_mean,
                                                     reduction_factor    = reduction_factor,
                                                     seed                = seed,
                                                     target_proportion_missing_of_this_run = target_proportion_missing_of_this_run)
      n_thr <-  ncol(x[[1]]) - 1
      counts <-  table_of_counts_of_each_study_and_threshold_from_the_real_data(real_x = list(nondiseased = x[[1]],
                                                                                              diseased    = x[[2]]))
      ##
      Hoyer <-  tryCatch(fit_the_Hoyer_2018_model_to_a_table_of_counts(counts            = counts,
                                                                       n_thr             = n_thr,
                                                                       distribution      = distribution_of_the_Hoyer_model,
                                                                       smallest          = smallest_of_the_Hoyer_model,
                                                                       largest           = largest_of_the_Hoyer_model,
                                                                       confidence_level  = confidence_level),
                         error = function(e) conditionMessage(e))
      GLMM <-  tryCatch(fit_the_logit_GLMM_to_a_table_of_counts(counts            = counts,
                                                                n_thr             = n_thr,
                                                                confidence_level  = confidence_level),
                        error = function(e) conditionMessage(e))
      ##
      fn_results <-  function(name_of_the_model, fit) {

            if (is.character(fit)) return(NULL)
            return(results_list_of_the_fit_of_a_competing_model(name_of_the_model  = name_of_the_model,
                                                                Se_and_Sp          = fit$Se_and_Sp,
                                                                x                  = x,
                                                                seed               = seed,
                                                                n_studies          = n_studies,
                                                                true_DGM_Se        = true_DGM_Se,
                                                                true_DGM_Sp        = true_DGM_Sp))

      }
      gate_Hoyer <-  if (is.character(Hoyer)) {
            data.frame(seed = seed, code = NA, n_messages = NA, singular = NA, seconds = NA, error = substr(Hoyer, 1, 120))
      } else {
            data.frame(seed = seed, code = Hoyer$convergence, n_messages = 0, singular = NA, seconds = Hoyer$seconds, error = "")
      }
      gate_GLMM <-  if (is.character(GLMM)) {
            data.frame(seed = seed, code = NA, n_messages = NA, singular = NA, seconds = NA, error = substr(GLMM, 1, 120))
      } else {
            data.frame(seed = seed, code = GLMM$code_of_the_fit, n_messages = length(GLMM$messages_of_the_fit),
                       singular = GLMM$singular, seconds = GLMM$seconds, error = "")
      }
      return(list(Hoyer_2018  = list(results = fn_results("Hoyer_2018", Hoyer), gate = gate_Hoyer),
                  logit_GLMM  = list(results = fn_results("logit_GLMM", GLMM), gate = gate_GLMM)))

}
