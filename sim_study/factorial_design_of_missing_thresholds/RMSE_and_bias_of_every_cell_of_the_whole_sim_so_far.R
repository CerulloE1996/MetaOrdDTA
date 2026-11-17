##
## ---- RMSE, bias, coverage and interval width of every cell of the whole sim so far, with the sim's own functions ------------------
##
## For each cell with 30 or more kept fits (the files seed_*_<cell>_results_list.RDS of the whole sim's folder, which
## holds the kept fits of the local HPC and of ALICE): compute_simulation_metrics and create_metrics_summary_table of
## thr_sim_study_helper_fns.R on the thresholds of interest, with the true values of get_true_values. One row per
## cell, with the columns of create_summary_row_MetaOrd (MetaOrdDTA_sim_helpers_2.R), which the plot functions of the
## paper read, plus the prior set, the number of kept fits, and MCSE_val = (largest MCSE of the RMSE over the
## thresholds + the larger of the two mean MCSEs) / 2, the number the sim compares with its MCSE target.
## Usage: Rscript <this file> <folder of the MetaOrdDTA package> <folder of the whole sim's results> <file to save, .rds>
##        <number of cells read at once> <fewest studies a threshold needs to be counted> <MCSE target of the sim>
##
arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 6)
local_pkg_dir <-  arguments[1]
folder <-  arguments[2]
file_to_save <-  arguments[3]
n_cells_at_once <-  as.integer(arguments[4])
min_studies_per_threshold <-  as.numeric(arguments[5])
target_MCSE_RMSE <-  as.numeric(arguments[6])
suppressPackageStartupMessages({ library(dplyr); library(parallel) })
source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
files <-  list.files(folder, pattern = "^seed_[0-9]+_.*_results_list[.]RDS$")
cell_of <-  sub("_results_list[.]RDS$", "", sub("^seed_[0-9]+_", "", files))
cells <-  names(which(table(cell_of) >= 30))
of_interest <-  list(GAD_2 = 1:6, HADS = 3:17, BAI = 5:43)
fn_one_cell <- function(cell) {

      test <-  sub("^test_([A-Z_0-9]+)red_fact.*", "\\1", cell)
      DGM <-  sub(".*1DGM_(.*)N_500n_studies.*", "\\1", cell)
      n_studies <-  as.numeric(sub(".*n_studies_([0-9]+)param.*", "\\1", cell))
      model_parameterisation <-  sub(".*param_(.*)soft_FALSE.*", "\\1", cell)
      random_thresholds <-  grepl("rand_thr_TRUE", cell)
      prior_set <-  sub(".*rand_thr_(TRUE|FALSE)(set_[0-9]+).*", "\\2", cell)
      results <-  lapply(file.path(folder, files[cell_of == cell]), function(f) tryCatch(readRDS(f), error = function(e) NULL))
      results <-  results[!vapply(results, is.null, logical(1))]
      true_values <-  get_true_values(test = test, DGM = DGM)
      metrics <-  tryCatch(compute_simulation_metrics(test = test,
                                                      complete_results_list = results,
                                                      output_dir = folder,
                                                      true_DGM_Se = true_values$Se,
                                                      true_DGM_Sp = true_values$Sp,
                                                      min_ESS_thr = if (n_studies == 10) 2500 else 500,
                                                      vec_index_inner = of_interest[[test]],
                                                      min_studies_per_threshold = min_studies_per_threshold),
                           error = function(e) NULL)
      if (is.null(metrics)) return(NULL)
      by_threshold <-  create_metrics_summary_table(metrics = metrics, vec_index_inner = of_interest[[test]])$by_threshold
      MCSE_Se <-  as.numeric(by_threshold$MCSE_RMSE_Se)
      MCSE_Sp <-  as.numeric(by_threshold$MCSE_RMSE_Sp)
      MCSE_val <-  (max(MCSE_Se, MCSE_Sp, na.rm = TRUE) + max(mean(MCSE_Se, na.rm = TRUE), mean(MCSE_Sp, na.rm = TRUE))) / 2
      o <-  metrics$overall
      data.frame(test = test, DGM = DGM, n_studies = n_studies, model_parameterisation = model_parameterisation,
                 prior_set = prior_set, random_thresholds = random_thresholds,
                 model = if (model_parameterisation == "ord_bivariate") { if (random_thresholds) "O-biv-RC" else "O-biv-FC" } else
                         if (model_parameterisation == "R&G") "ord-HSROC" else
                         if (model_parameterisation == "Bivariate") "strat. bivariate" else model_parameterisation,
                 kept_fits = length(results), n_sims = metrics$n_sims, n_excluded = o$n_excluded,
                 MCSE_val = MCSE_val, mcse_met = MCSE_val <= target_MCSE_RMSE,
                 avg_RMSE_Se = o$avg_RMSE_Se, avg_RMSE_Sp = o$avg_RMSE_Sp, avg_RMSE = o$avg_RMSE_Se + o$avg_RMSE_Sp,
                 max_RMSE_Se = o$max_RMSE_Se, max_RMSE_Sp = o$max_RMSE_Sp,
                 avg_MCSE_RMSE_Se = o$avg_MCSE_RMSE_Se, avg_MCSE_RMSE_Sp = o$avg_MCSE_RMSE_Sp,
                 max_MCSE_RMSE_Se = o$max_MCSE_RMSE_Se, max_MCSE_RMSE_Sp = o$max_MCSE_RMSE_Sp,
                 avg_abs_bias_Se = o$avg_abs_bias_Se, avg_abs_bias_Sp = o$avg_abs_bias_Sp,
                 avg_abs_bias = o$avg_abs_bias_Se + o$avg_abs_bias_Sp,
                 max_abs_bias_Se = o$max_abs_bias_Se, max_abs_bias_Sp = o$max_abs_bias_Sp,
                 avg_MCSE_bias_Se = o$avg_MCSE_bias_Se, avg_MCSE_bias_Sp = o$avg_MCSE_bias_Sp,
                 avg_coverage_Se = o$avg_coverage_prob_Se, avg_coverage_Sp = o$avg_coverage_prob_Sp,
                 avg_width_Se = o$avg_interval_width_Se, avg_width_Sp = o$avg_interval_width_Sp,
                 mean_min_ESS = metrics$mean_min_ESS,
                 stringsAsFactors = FALSE)

}
rows <-  mclapply(cells, function(cell) { sink(nullfile()); on.exit(sink()); fn_one_cell(cell) }, mc.cores = n_cells_at_once)
table_of_the_cells <-  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
saveRDS(table_of_the_cells, file_to_save)
cat("cells with 30 or more kept fits:", length(cells), "; cells with metrics:", nrow(table_of_the_cells), "; kept fits read:",
    sum(table_of_the_cells$kept_fits), "; cells with the MCSE target met:", sum(table_of_the_cells$mcse_met), "\n")
