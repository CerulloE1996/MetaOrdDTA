##
## ---- RMSE, bias, coverage, interval width and MCSE of every cell of every model of the whole sim, with the sim's own functions,
##      read from a LIST of folders (the kept fits of the sim sit in several folders, each model, prior set and DGM in one of them) -----
##
## The scoring of a cell is that of RMSE_and_bias_of_every_cell_of_the_whole_sim_so_far.R, unchanged: compute_simulation_metrics and
## create_metrics_summary_table of thr_sim_study_helper_fns.R on the thresholds of interest (GAD-2 1:6, HADS 3:17, BAI 5:43), with the
## true values of get_true_values; MCSE_val = (largest MCSE of the RMSE over the thresholds + the larger of the two mean MCSEs) / 2,
## the number the sim compares with its MCSE target.
## A cell = the files seed_*_<cell>_results_list.RDS of that name over the folders given; a file name is read once (the first folder
## of the list that holds it). Each folder is given with a pattern: only the cells whose name matches it are taken from that folder.
## Added to each row: the median seconds of a fit, from the gate table of the cell (gate_of_every_fit_<cell>.csv in the same
## folders) where it has the columns: burn-in + sampling, and the whole sample call.
##
## Usage: Rscript <this file> <folder of the MetaOrdDTA package> <text file: one line a folder, "<folder><TAB><pattern of the cells>">
##                <file to save, .rds> <number of cells read at once> <fewest studies a threshold needs to be counted>
##                <MCSE target of the sim> <fewest kept fits a cell needs to be scored>
##
arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 7)
local_pkg_dir <-  arguments[1]
list_of_folders <-  utils::read.delim(arguments[2], header = FALSE, stringsAsFactors = FALSE, col.names = c("folder", "pattern"), quote = "")
file_to_save <-  arguments[3]
n_cells_at_once <-  as.integer(arguments[4])
min_studies_per_threshold <-  as.numeric(arguments[5])
target_MCSE_RMSE <-  as.numeric(arguments[6])
fewest_kept_fits <-  as.integer(arguments[7])
stopifnot(all(dir.exists(list_of_folders$folder)))
suppressPackageStartupMessages({ library(dplyr); library(parallel) })
source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
##
## ---- The results files and the gate tables of every folder of the list:
##
files_of_every_folder <-  do.call(rbind, lapply(seq_len(nrow(list_of_folders)), function(i) {

      names_here <-  list.files(list_of_folders$folder[i], pattern = "^seed_[0-9]+_.*_results_list[.]RDS$")
      cell_here <-  sub("_results_list[.]RDS$", "", sub("^seed_[0-9]+_", "", names_here))
      taken <-  grepl(list_of_folders$pattern[i], cell_here)
      data.frame(folder = rep(list_of_folders$folder[i], sum(taken)), file = names_here[taken], cell = cell_here[taken], stringsAsFactors = FALSE)

}))
files_of_every_folder <-  files_of_every_folder[!duplicated(files_of_every_folder$file), ]
cells <-  names(which(table(files_of_every_folder$cell) >= fewest_kept_fits))
files_by_cell <-  split(file.path(files_of_every_folder$folder, files_of_every_folder$file), files_of_every_folder$cell)
of_interest <-  list(GAD_2 = 1:6, HADS = 3:17, BAI = 5:43)
name_of_the_model <-  function(model_parameterisation,
                               random_thresholds) {

      if (model_parameterisation == "ord_bivariate") return(if (random_thresholds) "O-biv-RC" else "O-biv-FC")
      switch(model_parameterisation, "R&G" = "ord-HSROC", "ord_HSROC" = "ord-HSROC", "Bivariate" = "strat. bivariate", "Jones" = "Jones", "Hoyer_2018" = "Hoyer 2018",
             "logit_GLMM" = "logit GLMM", "diagmeta" = "diagmeta", model_parameterisation)

}
fn_seconds_of_a_fit <- function(cell,
                                folders_of_the_cell) {

      gate_files <-  file.path(folders_of_the_cell, paste0("gate_of_every_fit_", cell, ".csv"))
      gate_files <-  gate_files[file.exists(gate_files)]
      nothing <-  c(seconds_burnin_and_sampling = NA_real_, seconds_sample_call = NA_real_)
      if (!length(gate_files)) return(nothing)
      gate <-  tryCatch(do.call(rbind, lapply(gate_files, function(f) {
            one <-  utils::read.csv(f, colClasses = "character")
            one[, intersect(c("seed", "seconds_burnin", "seconds_sampling", "seconds_sample_call"), names(one)), drop = FALSE]
      })), error = function(e) NULL)
      if (is.null(gate) || !("seed" %in% names(gate))) return(nothing)
      gate <-  gate[!is.na(suppressWarnings(as.numeric(gate$seed))), , drop = FALSE]
      fn_median <- function(x) suppressWarnings(stats::median(as.numeric(x), na.rm = TRUE))
      c(seconds_burnin_and_sampling = if (all(c("seconds_burnin", "seconds_sampling") %in% names(gate)))
                                            fn_median(as.numeric(gate$seconds_burnin) + as.numeric(gate$seconds_sampling)) else NA_real_,
        seconds_sample_call = if ("seconds_sample_call" %in% names(gate)) fn_median(gate$seconds_sample_call) else NA_real_)

}
fn_one_cell <- function(cell) {

      test <-  sub("^test_([A-Z_0-9]+)red_fact.*", "\\1", cell)
      DGM <-  sub(".*1DGM_(.*)N_500n_studies.*", "\\1", cell)
      n_studies <-  as.numeric(sub(".*n_studies_([0-9]+)param.*", "\\1", cell))
      model_parameterisation <-  sub(".*param_(.*)soft_FALSE.*", "\\1", cell)
      random_thresholds <-  grepl("rand_thr_TRUE", cell)
      prior_set <-  sub(".*rand_thr_(TRUE|FALSE)(set_[0-9]+).*", "\\2", cell)
      files_of_the_cell <-  files_by_cell[[cell]]
      results <-  lapply(files_of_the_cell, function(f) tryCatch(readRDS(f), error = function(e) NULL))
      n_unreadable <-  sum(vapply(results, is.null, logical(1)))
      results <-  results[!vapply(results, is.null, logical(1))]
      ## (the fits of the first sim carry the DGM name R&G_random_C in their file names; the package now calls that DGM ord_HSROC_random_C)
      true_values <-  get_true_values(test = test, DGM = if (DGM == "R&G_random_C") "ord_HSROC_random_C" else DGM)
      metrics <-  tryCatch(compute_simulation_metrics(test = test,
                                                      complete_results_list = results,
                                                      output_dir = dirname(files_of_the_cell[1]),
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
      seconds <-  fn_seconds_of_a_fit(cell, unique(dirname(files_of_the_cell)))
      o <-  metrics$overall
      data.frame(test = test, DGM = DGM, n_studies = n_studies, model_parameterisation = model_parameterisation,
                 prior_set = prior_set, random_thresholds = random_thresholds,
                 model = name_of_the_model(model_parameterisation, random_thresholds),
                 kept_fits = length(results), files_unreadable = n_unreadable, n_sims = metrics$n_sims, n_excluded = o$n_excluded,
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
                 seconds_burnin_and_sampling = unname(seconds["seconds_burnin_and_sampling"]),
                 seconds_sample_call = unname(seconds["seconds_sample_call"]),
                 folder = paste(basename(unique(dirname(files_of_the_cell))), collapse = " + "),
                 stringsAsFactors = FALSE)

}
rows <-  mclapply(cells, function(cell) { sink(nullfile()); on.exit(sink()); tryCatch(fn_one_cell(cell), error = function(e) NULL) }, mc.cores = n_cells_at_once)
table_of_the_cells <-  do.call(rbind, rows[!vapply(rows, function(r) is.null(r) || !is.data.frame(r), logical(1))])
saveRDS(table_of_the_cells, file_to_save)
cat("results files listed:", nrow(files_of_every_folder), "; cells with", fewest_kept_fits, "or more kept fits:", length(cells), "; cells with metrics:",
    nrow(table_of_the_cells), "; kept fits read:", sum(table_of_the_cells$kept_fits), "; files unreadable:", sum(table_of_the_cells$files_unreadable),
    "; cells with the MCSE target met:", sum(table_of_the_cells$mcse_met), "\n")
print(as.data.frame(table(model = table_of_the_cells$model, set = table_of_the_cells$prior_set)), row.names = FALSE)
