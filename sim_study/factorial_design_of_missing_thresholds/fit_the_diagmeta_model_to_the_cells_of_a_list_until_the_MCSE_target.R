options(warning.length = 8170, width = 250, error = function() { m <- geterrmessage(); cat("\nTAIL OF THE ERROR:", substr(m, max(1, nchar(m) - 2500), nchar(m)), "\n"); q(status = 1) })
##
## ---- The model of Steinhauser, Schumacher and Rücker (2016) (R package diagmeta) fitted to the cells of a list, each cell until
##      the model is under the MCSE target of the sim (or the largest number of seeds is reached).
##      Each seed's data set is the one the sim runner makes for that cell and seed
##      (fit_the_competing_models_to_the_seeds_of_one_cell_of_the_sim.R). The MCSE is the sim's own: the thresholds of interest of
##      the test, thresholds reported by at least 1 study.
##      The model is frequentist: it has no prior, so each cell is fitted once and is filed with the cells of prior set 0.
##      A fit is kept when lmer's code is 0 and every Se, Sp and interval limit is a number; a fit whose covariance of the study
##      effects is singular is kept and counted. 95% intervals.
##
## ---- Usage: Rscript <this file> <list file> <number of cores> <folder of the results> <largest number of seeds of a cell>
##                                 <seeds per batch> <variant of diagmeta: e.g. DIDS> <distribution: logistic | normal>
##      The list file: one cell per line, fields separated by one space: <test> <number of studies> <DGM>
##                     <proportion of the thresholds missing: real_pattern | 0.15 | 0.55 | 0.85>
##
command_line_arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(command_line_arguments) == 7)
list_file <-  command_line_arguments[1]
n_cores <-  as.numeric(command_line_arguments[2])
folder_of_the_results <-  command_line_arguments[3]
largest_number_of_seeds <-  as.numeric(command_line_arguments[4])
seeds_per_batch <-  as.numeric(command_line_arguments[5])
model_of_diagmeta <-  command_line_arguments[6]
distr_of_diagmeta <-  command_line_arguments[7]
stopifnot(model_of_diagmeta %in% c("CI", "DI", "CS", "DS", "CICS", "DICS", "CIDS", "DIDS"), distr_of_diagmeta %in% c("logistic", "normal"),
          seeds_per_batch >= 1, largest_number_of_seeds >= seeds_per_batch)
##
folder_of_the_packages <-  "/home/enzocerullo/Documents/Work/PhD_work/R_packages"
folder_of_the_sim_outputs <-  file.path(folder_of_the_packages, "MetaOrdDTA_outputs", "sim_study")
library_of_the_competing_methods <-  file.path(folder_of_the_packages, "R_library_for_the_competing_methods_of_the_sim_study")
library_of_the_sim <-  file.path(folder_of_the_packages, "R_library_for_the_sim_on_the_local_HPC_with_the_final_O_biv_RC_of_11_Nov")
package_folder <-  file.path(folder_of_the_packages, "MetaOrdDTA")
stopifnot(file.exists(list_file), dir.exists(folder_of_the_results), dir.exists(library_of_the_competing_methods),
          dir.exists(library_of_the_sim), dir.exists(package_folder))
## (diagmeta is in the library of the competing methods; lme4, which it calls, is found there or in R's own libraries)
.libPaths(c(library_of_the_competing_methods, .libPaths()))
suppressMessages({
      library(dplyr)
      library(stringr)
      library(MetaOrdDTA, lib.loc = library_of_the_sim)
      library(lme4)
      library(diagmeta, lib.loc = library_of_the_competing_methods)
})
for (f in c("fit_the_Hoyer_2018_time_to_event_model_to_one_data_set_of_the_sim.R",
            "fit_the_diagmeta_model_to_one_data_set_of_the_sim.R",
            "fit_the_competing_models_to_the_seeds_of_one_cell_of_the_sim.R")) source(file.path(folder_of_the_sim_outputs, f))
for (f in c("MetaOrdDTA_sim_helpers_1.R", "MetaOrdDTA_sim_helpers_2.R", "MetaOrdDTA_sim_helpers_3.R",
            "thr_sim_study_helper_fns.R")) source(file.path(package_folder, "sim_study", f))
target <-  0.125
thresholds_of_interest_of_each_test <-  list(GAD_2 = 1:6, HADS = 3:17, BAI = 5:43)
cells <-  utils::read.table(list_file, header = FALSE, sep = " ", stringsAsFactors = FALSE, col.names = c("test", "n_studies", "DGM", "proportion_missing"),
                            colClasses = c("character", "numeric", "character", "character"),
                            comment.char = "")
stopifnot(all(cells$test %in% names(thresholds_of_interest_of_each_test)), all(cells$n_studies %in% c(10, 50)),
          all(cells$proportion_missing %in% c("real_pattern", "0.15", "0.55", "0.85")))
model <-  "diagmeta"
label_of_the_model <-  paste0("diagmeta (", model_of_diagmeta, ", ", distr_of_diagmeta, ")")
cat(label_of_the_model, "; diagmeta", as.character(utils::packageVersion("diagmeta")), "; cells", nrow(cells), "; cores", n_cores, ";", date(), "\n")
for (row_of_the_cell in seq_len(nrow(cells))) {

      test <-  cells$test[row_of_the_cell]
      n_studies <-  cells$n_studies[row_of_the_cell]
      DGM <-  cells$DGM[row_of_the_cell]
      ## (NULL = the real reporting pattern of the test, as in the first sim; or the proportion of the thresholds to make missing)
      target_proportion_missing_of_this_run <-  if (cells$proportion_missing[row_of_the_cell] == "real_pattern") NULL else
                                                as.numeric(cells$proportion_missing[row_of_the_cell])
      thresholds_of_interest <-  thresholds_of_interest_of_each_test[[test]]
      ## (these fits have no draws and their ESS is written as infinite, so this value drops no fit; the sim's values are kept)
      min_ESS_of_the_cell <-  if (n_studies == 10) 2500 else 500
      folder_of_the_cell <-  file.path(folder_of_the_results, paste0("DGM_", DGM))
      dir.create(folder_of_the_cell, showWarnings = FALSE)
      truth <-  get_true_values(DGM = DGM, test = test)
      n_thr <-  length(truth$Se)
      stopifnot(n_thr == list(GAD_2 = 6, HADS = 21, BAI = 63)[[test]], length(truth$Sp) == n_thr)
      name_of_the_cell <-  paste0("test_", test, "red_fact_1DGM_", DGM, "N_500n_studies_", n_studies, "param_", model, "soft_FALSErand_thr_FALSEset_0")
      fn_file_of_the_results_of_a_seed <- function(seed) {
            file.path(folder_of_the_cell, paste0("seed_", seed, "_", name_of_the_cell, "_results_list.RDS"))
      }
      gate_file <-  file.path(folder_of_the_cell, paste0("gate_of_every_fit_", name_of_the_cell, ".csv"))
      fn_score <- function(list_of_results, seconds) {

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
            data.frame(test = test, n_studies = n_studies, DGM = DGM, model = label_of_the_model, data_sets = length(list_of_results),
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
      ## (one seed: its data set, the fit, the results list saved when the fit is kept; returns the seed's row of the gate table)
      fn_one_seed <- function(seed) {

            fit <-  NULL
            invisible(capture.output(fit <-  tryCatch(suppressWarnings(suppressMessages({
                  x <-  simulated_data_set_of_one_seed_of_a_cell(local_pkg_dir     = package_folder,
                                                                 DGM               = DGM,
                                                                 test              = test,
                                                                 n_studies         = n_studies,
                                                                 N_per_study_mean  = 500,
                                                                 reduction_factor  = 1.0,
                                                                 seed              = seed,
                                                                 target_proportion_missing_of_this_run = target_proportion_missing_of_this_run)
                  counts <-  table_of_counts_of_each_study_and_threshold_from_the_real_data(real_x = list(nondiseased = x[[1]],
                                                                                                          diseased    = x[[2]]))
                  result <-  fit_the_diagmeta_model_to_a_table_of_counts(counts            = counts,
                                                                        n_thr             = n_thr,
                                                                        model_of_diagmeta = model_of_diagmeta,
                                                                        distr_of_diagmeta = distr_of_diagmeta,
                                                                        confidence_level  = 0.95)
                  result$x <-  x
                  result$fit <-  NULL
                  result
            })), error = function(e) conditionMessage(e))))
            if (is.character(fit)) {
                  return(data.frame(seed = seed, code = NA, n_messages = NA, singular = NA, seconds = NA, every_value_is_a_number = NA,
                                    kept = FALSE, error = substr(gsub("[\r\n,]", " ", fit), 1, 160)))
            }
            every_value_is_a_number <-  all(is.finite(as.matrix(fit$Se_and_Sp[, c("Sens", "SensCI_Lower", "SensCI_Upper", "Spec", "SpecCI_Lower", "SpecCI_Upper")])))
            kept <-  isTRUE(fit$code_of_the_fit == 0) && every_value_is_a_number
            if (kept) {
                  saveRDS(results_list_of_the_fit_of_a_competing_model(name_of_the_model  = model,
                                                                       Se_and_Sp          = fit$Se_and_Sp,
                                                                       x                  = fit$x,
                                                                       seed               = seed,
                                                                       n_studies          = n_studies,
                                                                       true_DGM_Se        = truth$Se,
                                                                       true_DGM_Sp        = truth$Sp),
                          fn_file_of_the_results_of_a_seed(seed))
            }
            return(data.frame(seed = seed, code = fit$code_of_the_fit, n_messages = length(fit$messages_of_the_fit), singular = fit$singular,
                              seconds = fit$seconds, every_value_is_a_number = every_value_is_a_number, kept = kept, error = ""))

      }
      cat("CELL", test, n_studies, "studies, DGM", DGM, ": start", date(), "\n")
      n_done <-  if (file.exists(gate_file)) max(utils::read.csv(gate_file)$seed) else 0
      state <-  NULL
      repeat {
            if (n_done > 0) {
                  gate <-  utils::read.csv(gate_file)
                  list_of_results <-  lapply(seq_len(n_done), function(seed) {
                        f <-  fn_file_of_the_results_of_a_seed(seed)
                        if (file.exists(f)) readRDS(f) else NULL
                  })
                  list_of_results <-  list_of_results[!vapply(list_of_results, is.null, logical(1))]
                  if (length(list_of_results) == 0) {
                        ## (no fit of the cell could be kept: the model cannot be fitted to these data sets, as when every study
                        ##  reports one threshold only; the cell is written down as not fitted and the next cell is taken)
                        cat("NOT FITTED", test, n_studies, "studies, DGM", DGM, ": no fit kept in seeds 1 to", n_done, "| error of the first seed:",
                            as.character(gate$error[1]), "|", date(), "\n")
                        state <-  NULL
                        break
                  }
                  state <-  cbind(fn_score(list_of_results, median(gate$seconds, na.rm = TRUE)),
                                  fits = nrow(gate),
                                  percent_of_fits_not_kept = round(100 * mean(!gate$kept), 2),
                                  percent_of_fits_with_an_error = round(100 * mean(is.na(gate$code)), 2),
                                  percent_of_kept_fits_singular = round(100 * mean(gate$singular[gate$kept]), 2))
                  cat("RESULT", test, n_studies, "studies, DGM", DGM, ": seeds 1 to", n_done, "| MCSE", state$MCSE, "| kept", state$data_sets,
                      "| not kept", state$percent_of_fits_not_kept, "% | singular among the kept", state$percent_of_kept_fits_singular, "% |", date(), "\n")
                  if (state$MCSE < target || n_done >= largest_number_of_seeds) break
            }
            new_seeds <-  (n_done + 1):(n_done + seeds_per_batch)
            gates <-  parallel::mclapply(new_seeds, fn_one_seed, mc.cores = n_cores)
            stopifnot(all(vapply(gates, is.data.frame, logical(1))))
            rows_before <-  if (file.exists(gate_file)) utils::read.csv(gate_file) else NULL
            write.csv(rbind(rows_before, do.call(rbind, gates)), gate_file, row.names = FALSE)
            n_done <-  n_done + seeds_per_batch
      }
      if (is.null(state)) {
            cat("CELL", test, n_studies, "studies, DGM", DGM, ": end, not fitted", date(), "\n")
            next
      }
      saveRDS(state, file.path(folder_of_the_results, paste0("band_summary_of_the_diagmeta_model_at_the_MCSE_target_test_", test,
                                                             "_n_studies_", n_studies, "_DGM_", DGM, ".rds")))
      lines_of_the_summary <-  capture.output(print(format(state, digits = 3), row.names = FALSE))
      cat(paste("SUMMARY ", lines_of_the_summary), sep = "\n")
      cat("CELL", test, n_studies, "studies, DGM", DGM, ": end", date(), "\n")

}
cat("EVERY CELL OF THE LIST IS DONE", date(), "\n")
