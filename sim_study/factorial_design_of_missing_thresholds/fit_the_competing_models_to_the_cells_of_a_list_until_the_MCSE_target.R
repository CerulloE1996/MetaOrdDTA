options(warning.length = 8170, width = 250, error = function() { m <- geterrmessage(); cat("\nTAIL OF THE ERROR:", substr(m, max(1, nchar(m) - 2500), nchar(m)), "\n"); q(status = 1) })
##
## ---- The competing models (Hoyer et al. 2018, time-to-event; the logit GLMM) fitted to the cells of a list, each cell
##      until BOTH models are under the MCSE target of the sim (or the largest number of seeds is reached).
##      The same code as the job of 11 Nov that took the GAD-2, 10-study cells to the target; here the test, the number of
##      studies and the DGM of a cell are read from the list. Each seed's data set is the one the sim runner makes for that
##      cell and seed (fit_the_competing_models_to_the_seeds_of_one_cell_of_the_sim.R). The MCSE is the sim's own: the
##      thresholds of interest of the test, thresholds reported by at least 1 study.
##      The settings of the two models stand as on 11 Nov: Hoyer 2018 Weibull, ends 0.5 and (number of thresholds + 1);
##      the logit GLMM with lme4's own settings; a fit is kept when its code is 0; 95% intervals.
##
## ---- Usage: Rscript <this file> <list file> <number of cores> <folder of the results> <largest number of seeds of a cell>
##      The list file: one cell per line, fields separated by one space: <test> <number of studies> <DGM>
##                     <proportion of the thresholds missing: real_pattern | 0.15 | 0.55 | 0.85>
##
command_line_arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(command_line_arguments) == 4)
list_file <-  command_line_arguments[1]
n_cores <-  as.numeric(command_line_arguments[2])
folder_of_the_results <-  command_line_arguments[3]
largest_number_of_seeds <-  as.numeric(command_line_arguments[4])
##
folder_of_the_packages <-  "/home/enzocerullo/Documents/Work/PhD_work/R_packages"
folder_of_the_sim_outputs <-  file.path(folder_of_the_packages, "MetaOrdDTA_outputs", "sim_study")
library_of_the_competing_methods <-  file.path(folder_of_the_packages, "R_library_for_the_competing_methods_of_the_sim_study")
library_of_the_sim <-  file.path(folder_of_the_packages, "R_library_for_the_sim_on_the_local_HPC_with_the_final_O_biv_RC_of_11_Nov")
package_folder <-  file.path(folder_of_the_packages, "MetaOrdDTA")
stopifnot(file.exists(list_file), dir.exists(folder_of_the_results), dir.exists(library_of_the_competing_methods),
          dir.exists(library_of_the_sim), dir.exists(package_folder))
suppressMessages({
      library(dplyr)
      library(stringr)
      library(MetaOrdDTA, lib.loc = library_of_the_sim)
      library(dtametaTMB, lib.loc = library_of_the_competing_methods)
})
reduction_factor <-  1.0
for (f in c("fit_the_Hoyer_2018_time_to_event_model_to_one_data_set_of_the_sim.R",
            "fit_the_logit_GLMM_to_one_data_set_of_the_sim.R",
            "fit_the_competing_models_to_the_seeds_of_one_cell_of_the_sim.R")) source(file.path(folder_of_the_sim_outputs, f))
for (f in c("MetaOrdDTA_sim_helpers_1.R", "MetaOrdDTA_sim_helpers_2.R", "MetaOrdDTA_sim_helpers_3.R",
            "thr_sim_study_helper_fns.R")) source(file.path(package_folder, "sim_study", f))
target <-  0.125
seeds_per_batch <-  400
thresholds_of_interest_of_each_test <-  list(GAD_2 = 1:6, HADS = 3:17, BAI = 5:43)
cells <-  utils::read.table(list_file, header = FALSE, sep = " ", stringsAsFactors = FALSE, col.names = c("test", "n_studies", "DGM", "proportion_missing"),
                            colClasses = c("character", "numeric", "character", "character"),
                            comment.char = "")
stopifnot(all(cells$test %in% names(thresholds_of_interest_of_each_test)), all(cells$n_studies %in% c(10, 50)),
          all(cells$proportion_missing %in% c("real_pattern", "0.15", "0.55", "0.85")))
models <-  c("Hoyer_2018", "logit_GLMM")
for (row_of_the_cell in seq_len(nrow(cells))) {

      test <-  cells$test[row_of_the_cell]
      n_studies <-  cells$n_studies[row_of_the_cell]
      DGM <-  cells$DGM[row_of_the_cell]
      ## (NULL = the real reporting pattern of the test, as in the first sim; or the proportion of the thresholds to make missing)
      target_proportion_missing_of_this_run <-  if (cells$proportion_missing[row_of_the_cell] == "real_pattern") NULL else
                                                as.numeric(cells$proportion_missing[row_of_the_cell])
      thresholds_of_interest <-  thresholds_of_interest_of_each_test[[test]]
      min_ESS_of_the_cell <-  if (n_studies == 10) 2500 else 500
      folder_of_the_cell <-  file.path(folder_of_the_results, paste0("DGM_", DGM))
      dir.create(folder_of_the_cell, showWarnings = FALSE)
      truth <-  get_true_values(DGM = DGM, test = test)
      n_thr <-  length(truth$Se)
      stopifnot(n_thr == list(GAD_2 = 6, HADS = 21, BAI = 63)[[test]], length(truth$Sp) == n_thr)
      fn_name_of_the_cell <- function(model) {
            paste0("test_", test, "red_fact_1DGM_", DGM, "N_500n_studies_", n_studies, "param_", model, "soft_FALSErand_thr_FALSEset_0")
      }
      fn_file_of_the_results_of_a_seed <- function(seed, model) {
            file.path(folder_of_the_cell, paste0("seed_", seed, "_", fn_name_of_the_cell(model), "_results_list.RDS"))
      }
      fn_gate_file <- function(model) file.path(folder_of_the_cell, paste0("gate_of_every_fit_", fn_name_of_the_cell(model), ".csv"))
      fn_score <- function(list_of_results, label, seconds) {

            metrics <-  NULL
            invisible(capture.output(metrics <-  suppressWarnings(compute_simulation_metrics(test = test,
                                                                                            output_dir = folder_of_the_cell,
                                                                                            complete_results_list = list_of_results,
                                                                                            true_DGM_Se = truth$Se,
                                                                                            true_DGM_Sp = truth$Sp,
                                                                                            min_ESS_thr = min_ESS_of_the_cell,
                                                                                            vec_index_inner = thresholds_of_interest,
                                                                                            min_studies_per_threshold = 1))))
            by_threshold <-  suppressWarnings(create_metrics_summary_table(metrics = metrics, vec_index_inner = thresholds_of_interest))$by_threshold
            fn_number <- function(v) as.numeric(sub("%", "", v))
            ## (the sim's own rule, check_MCSE_status in MetaOrdDTA_sim_helpers_3.R: a threshold of interest with no value -
            ##  reported by no study in every data set - is left out)
            MCSE_of_the_cell <-  (max(max(fn_number(by_threshold$MCSE_RMSE_Se), na.rm = TRUE), max(fn_number(by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)) +
                                  max(mean(fn_number(by_threshold$MCSE_RMSE_Se), na.rm = TRUE), mean(fn_number(by_threshold$MCSE_RMSE_Sp), na.rm = TRUE))) / 2
            data.frame(test = test, n_studies = n_studies, DGM = DGM, model = label, data_sets = length(list_of_results),
                       MCSE = round(MCSE_of_the_cell, 3),
                       thresholds_of_interest = length(thresholds_of_interest),
                       thresholds_with_a_value = sum(!is.na(fn_number(by_threshold$MCSE_RMSE_Se))),
                       mean_abs_bias_Se = mean(abs(fn_number(by_threshold$Bias_Se)), na.rm = TRUE), mean_abs_bias_Sp = mean(abs(fn_number(by_threshold$Bias_Sp)), na.rm = TRUE),
                       mean_RMSE_Se = mean(fn_number(by_threshold$RMSE_Se), na.rm = TRUE), mean_RMSE_Sp = mean(fn_number(by_threshold$RMSE_Sp), na.rm = TRUE),
                       lowest_coverage_Se = min(fn_number(by_threshold$Coverage_Se), na.rm = TRUE), mean_coverage_Se = mean(fn_number(by_threshold$Coverage_Se), na.rm = TRUE),
                       lowest_coverage_Sp = min(fn_number(by_threshold$Coverage_Sp), na.rm = TRUE), mean_coverage_Sp = mean(fn_number(by_threshold$Coverage_Sp), na.rm = TRUE),
                       mean_width_Se = mean(fn_number(by_threshold$Width_Se), na.rm = TRUE), mean_width_Sp = mean(fn_number(by_threshold$Width_Sp), na.rm = TRUE),
                       median_seconds_per_fit = round(seconds, 2), cores = 1)

      }
      fn_one_seed <- function(seed) {

            fits <-  NULL
            invisible(capture.output(fits <-  tryCatch(suppressWarnings(suppressMessages(
                  fits_of_the_competing_models_to_one_seed(local_pkg_dir = package_folder,
                                                           DGM = DGM,
                                                           test = test,
                                                           n_studies = n_studies,
                                                           N_per_study_mean = 500,
                                                           reduction_factor = 1.0,
                                                           seed = seed,
                                                           true_DGM_Se = truth$Se,
                                                           true_DGM_Sp = truth$Sp,
                                                           distribution_of_the_Hoyer_model = "weibull",
                                                           smallest_of_the_Hoyer_model = 0.5,
                                                           largest_of_the_Hoyer_model = n_thr + 1,
                                                           confidence_level = 0.95,
                                                           target_proportion_missing_of_this_run = target_proportion_missing_of_this_run))),
                  error = function(e) conditionMessage(e))))
            if (is.character(fits)) {
                  gate <-  data.frame(seed = seed, code = NA, n_messages = NA, singular = NA, seconds = NA, error = substr(fits, 1, 120))
                  return(list(Hoyer_2018 = gate, logit_GLMM = gate))
            }
            for (model in models) {
                  if (!is.null(fits[[model]]$results) && isTRUE(fits[[model]]$gate$code == 0)) {
                        saveRDS(fits[[model]]$results, fn_file_of_the_results_of_a_seed(seed, model))
                  }
            }
            list(Hoyer_2018 = fits$Hoyer_2018$gate, logit_GLMM = fits$logit_GLMM$gate)

      }
      cat("CELL", test, n_studies, "studies, DGM", DGM, ": start", format(Sys.time(), "%d %b %H:%M:%S"), "\n")
      n_done <-  if (file.exists(fn_gate_file("Hoyer_2018"))) max(utils::read.csv(fn_gate_file("Hoyer_2018"))$seed) else 0
      state <-  list()
      repeat {
            if (n_done > 0) {
                  for (model in models) {
                        gate <-  utils::read.csv(fn_gate_file(model))
                        list_of_results <-  lapply(seq_len(n_done), function(seed) {
                              f <-  fn_file_of_the_results_of_a_seed(seed, model)
                              if (file.exists(f)) readRDS(f) else NULL
                        })
                        list_of_results <-  list_of_results[!vapply(list_of_results, is.null, logical(1))]
                        state[[model]] <-  cbind(fn_score(list_of_results, if (model == "Hoyer_2018") "Hoyer 2018" else "Logit GLMM",
                                                          median(gate$seconds, na.rm = TRUE)),
                                                 fits = nrow(gate), fits_not_kept = sum(!gate$kept))
                  }
                  cat("RESULT", test, n_studies, "studies, DGM", DGM, ": seeds 1 to", n_done, "| MCSE Hoyer 2018", state$Hoyer_2018$MCSE,
                      "(kept", state$Hoyer_2018$data_sets, ") | logit GLMM", state$logit_GLMM$MCSE, "(kept", state$logit_GLMM$data_sets, ") |",
                      format(Sys.time(), "%H:%M:%S"), "\n")
                  if (all(vapply(state, function(z) z$MCSE, numeric(1)) < target) || n_done >= largest_number_of_seeds) break
            }
            new_seeds <-  (n_done + 1):(n_done + seeds_per_batch)
            gates <-  parallel::mclapply(new_seeds, fn_one_seed, mc.cores = n_cores)
            for (model in models) {
                  rows_to_add <-  do.call(rbind, lapply(gates, function(z) z[[model]]))
                  rows_to_add$kept <-  !is.na(rows_to_add$code) & rows_to_add$code == 0
                  rows_before <-  if (file.exists(fn_gate_file(model))) utils::read.csv(fn_gate_file(model)) else NULL
                  write.csv(rbind(rows_before, rows_to_add), fn_gate_file(model), row.names = FALSE)
            }
            n_done <-  n_done + seeds_per_batch
      }
      summary_of_the_cell <-  do.call(rbind, state)
      rownames(summary_of_the_cell) <-  NULL
      saveRDS(summary_of_the_cell, file.path(folder_of_the_results, paste0("band_summary_of_the_two_competing_models_at_the_MCSE_target_test_", test,
                                                                           "_n_studies_", n_studies, "_DGM_", DGM, ".rds")))
      lines_of_the_summary <-  capture.output(print(format(summary_of_the_cell, digits = 3), row.names = FALSE))
      cat(paste("SUMMARY ", lines_of_the_summary), sep = "\n")
      cat("CELL", test, n_studies, "studies, DGM", DGM, ": end", format(Sys.time(), "%d %b %H:%M:%S"), "\n")

}
cat("EVERY CELL OF THE LIST IS DONE", format(Sys.time(), "%d %b %H:%M:%S"), "\n")
