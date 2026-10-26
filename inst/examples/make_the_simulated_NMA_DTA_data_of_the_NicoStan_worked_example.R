##
## ==== make_the_simulated_NMA_DTA_data_of_the_NicoStan_worked_example.R =========================================
##
## ---- PRIVATE (never copied into NicoStan): the simulated network meta-analysis (NMA) data of NicoStan's worked
##      example 4 (MetaOrdDTA's NMA-DTA model with Xu's fixed cutpoints), made with MetaOrdDTA's own simulator,
##      R_fn_sim_NMA_data_with_covariates_varying() (R/R_fn_sim_data_ord_NMA.R), in exactly the flow and with
##      exactly the settings of inst/examples/NMA_example_main.R (Xu, fixed cutpoints, N_per_study_mean = 500,
##      N_per_study_SD = 500, missing thresholds, true_Mean_prev = 0.25, true_SD_probit_prev = 0.25, seed 1,
##      the covariate DGM MR_model_5), through inst/examples/NMA_missing_thr_prep_data_v2.R.
## ---- The simulator reads the real data of Klaus et al. (CONFIDENTIAL) for its design only: the patterns of
##      tests of the studies, the covariate rows it draws from, and the proportion of missing thresholds of each
##      test. Only SIMULATED quantities are saved: the simulated counts x (missing thresholds coded -1), which
##      index tests each simulated study has, their number and the tests' names. Nothing real (no covariate
##      value, no count, no study identifier) is saved, printed or copied.
##
## Writes:  NicoStan/inst/examples/data/simulated_network_meta_analysis_of_four_ordinal_screening_tests_with_
##          missing_thresholds.rds (never remade or overwritten: the script stops at once when it exists)
## Usage:   Rscript make_the_simulated_NMA_DTA_data_of_the_NicoStan_worked_example.R
##
##
## ------- Settings (those of NMA_example_main.R):  --------------------------------------------------------------
##
{
      local_pkg_dir <-  path.expand("~/Documents/Work/PhD_work/R_packages/MetaOrdDTA")
      output_file <-  file.path(path.expand("~/Documents/Work/PhD_work/R_packages/NicoStan/inst/examples/data"),
                                paste0("simulated_network_meta_analysis_of_four_ordinal_screening_tests_with_",
                                       "missing_thresholds.rds"))
      ##
      model_parameterisation <-  "Xu"
      box_cox <-  FALSE
      cts <-  FALSE
      random_thresholds <-  FALSE
      Dirichlet_random_effects_type <-  "fixed"
      ##
      N_per_study_mean <-  500
      N_per_study_SD <-  N_per_study_mean
      n_studies <-  25
      missing_thr <-  TRUE
      true_Mean_prev <-  0.25
      true_SD_probit_prev <-  0.25
      index_test_chosen_index <-  1 + 1
      network <-  TRUE
      softplus <-  FALSE
      ##
      seed <-  1
      intercept_only_for_sim_DGM <-  FALSE
      MR_model_for_sim_DGM <-  "MR_model_5"
      ##
      names_of_the_index_tests <-  c("GAD-2", "GAD-7", "HADS", "BAI")
}
##
## ---- A saved data file is never remade:
##
if (file.exists(output_file)) {
      message(paste0("\033[36m", "Already made (not remade): ", output_file, "\033[0m"))
      quit(save = "no", status = 0)
}
##
##
## ---- The simulation, in the flow of NMA_example_main.R (its sourced files use MetaOrdDTA's functions without
##      the namespace, hence MetaOrdDTA attached, as there): -----------------------------------------------------
##
start_time <-  Sys.time()
suppressPackageStartupMessages(library(MetaOrdDTA))
setwd(local_pkg_dir)
set.seed(seed)
source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_NMA.R"))
source(file.path(local_pkg_dir, "inst", "examples", "NMA_missing_thr_prep_data_v2.R"))
##
if (missing_thr == TRUE) {
      x <-  x_NMA_w_missing_thr
} else {
      x <-  x_NMA
}
n_studies <-  nrow(x[[1]][[1]])
##
## ---- x[[c]][[t]] (class, then test) reshaped to x[[t]][[c]] (test, then class), as in NMA_example_main.R:
##
x_old <-  x
x <-  list()
for (t in 1:n_index_tests) {
      x[[t]] <-  list()
}
for (t in 1:n_index_tests) {
      for (c in 1:2) {
            x[[t]][[c]] <-  x_old[[c]][[t]]
      }
}
indicator_index_test_in_study <-  MetaOrdDTA:::create_indicator_index_test_in_study(x)
##
##
## ------- Only the simulated quantities, with no names or attributes that could come from the real data:  -------
##
for (t in 1:n_index_tests) {
      for (c in 1:2) {
            x[[t]][[c]] <-  matrix(as.numeric(x[[t]][[c]]), nrow = nrow(x[[t]][[c]]), ncol = ncol(x[[t]][[c]]))
      }
      names(x[[t]]) <-  c("non_diseased", "diseased")
}
names(x) <-  names_of_the_index_tests
indicator_index_test_in_study <-  matrix(as.numeric(indicator_index_test_in_study),
                                         nrow = nrow(indicator_index_test_in_study),
                                         ncol = ncol(indicator_index_test_in_study),
                                         dimnames = list(NULL, names_of_the_index_tests))
simulated_NMA_data <-  list(x                              = x,
                            indicator_index_test_in_study  = indicator_index_test_in_study,
                            n_index_tests_per_study        = as.numeric(rowSums(indicator_index_test_in_study)),
                            names_of_the_index_tests       = names_of_the_index_tests,
                            n_studies                      = n_studies,
                            how_the_data_were_simulated    = paste0("MetaOrdDTA's R_fn_sim_NMA_data_with_",
                                                                    "covariates_varying(), seed ", seed,
                                                                    ", N_per_study_mean ", N_per_study_mean,
                                                                    ", N_per_study_SD ", N_per_study_SD,
                                                                    ", true_Mean_prev ", true_Mean_prev,
                                                                    ", true_SD_probit_prev ",
                                                                    true_SD_probit_prev,
                                                                    ", missing thresholds coded -1"))
##
## ---- The check that nothing real is saved: every count matrix plain numeric, with no names or dimnames, and
##      only the fields written here:
for (t in 1:n_index_tests) {
      for (c in 1:2) {
            stopifnot(is.numeric(x[[t]][[c]]),
                      is.null(dimnames(x[[t]][[c]])))
      }
}
stopifnot(identical(names(simulated_NMA_data), c("x", "indicator_index_test_in_study", "n_index_tests_per_study",
                                                 "names_of_the_index_tests", "n_studies",
                                                 "how_the_data_were_simulated")))
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(simulated_NMA_data, output_file)
##
message(paste0("\033[32m", "Saved: ", output_file, "\033[0m"))
message(paste0("\033[36m", "Fields and dims: ", "\033[0m"))
for (t in 1:n_index_tests) {
      message(paste0("\033[36m", "  x$", names(x)[t], ": 2 matrices of ", nrow(x[[t]][[1]]), " x ",
                     ncol(x[[t]][[1]]), " (", sum(x[[t]][[1]] == -1), " entries -1 in the non-diseased)",
                     "\033[0m"))
}
message(paste0("\033[36m", "  indicator_index_test_in_study: ", nrow(indicator_index_test_in_study), " x ",
               ncol(indicator_index_test_in_study), "; studies per test: ",
               paste0(colSums(indicator_index_test_in_study), collapse = ", "), "\033[0m"))
message(paste0("\033[36m", "Took ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2),
               " seconds", "\033[0m"))






















