##
## ---- Set the sampling iterations of every cell from the 10-seeds-per-cell pilot ----------------------------------------------------
##
## (reads every gate_of_every_fit_*.csv of the pilot; for each cell with all 10 fits, the min ESS of Se and Sp per
##  transition of its fits that neither exceed 1% divergent nor had a stuck chain, r = min_ESS / n_transitions, and its
##  5th percentile q05 = quantile(r, 0.05) give the transitions that put about 95% of the fits over the cell's ESS
##  threshold (2500 at 10 studies, 500 at 50): N_ITER_FOR_4_CHAINS = threshold / q05 / 4, rounded up to 50, at most
##  5000; the settings table has one row per cell and the iterations only, as every cell runs the same sampler settings)
##
## ---- Controls ----------------------------------------------------------------------------------------------------------------------
##
arguments <- commandArgs(trailingOnly = TRUE)
root <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs/sim_study"
folders_of_the_pilot <- if (length(arguments) > 0) arguments else
                        file.path(root, "pilot_10_seeds_per_cell_confirmed_settings")
folder_of_the_tables <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/sim_study/sim_study_redo_on_ALICE"
settings_file_after_the_pilot <- file.path(folder_of_the_tables, "settings_of_every_cell_of_the_whole_sim_after_the_pilot.csv")
n_chains <- 4
n_seeds <- 10
max_n_iter <- 5000
##
## ---- One row per cell --------------------------------------------------------------------------------------------------------------
##
fn_iterations_of_one_cell <- function(files_of_the_gate_tables_of_the_cell) {

      gate <- do.call(rbind, lapply(files_of_the_gate_tables_of_the_cell, read.csv))
      gate <- gate[!is.na(gate$min_ESS), ]
      cell <- sub("gate_of_every_fit_test_", "", basename(files_of_the_gate_tables_of_the_cell[1]))
      ## (a fit over 1% divergent fails the divergence gate whatever its length, and a fit whose min ESS is under a tenth
      ##  of the cell's median had a stuck chain; neither is a question of length, so both are left out of the ESS rate)
      failed <- gate$percent_divergent > 1 | gate$min_ESS < 0.1 * median(gate$min_ESS)
      r <- (gate$min_ESS / gate$n_transitions)[!failed]
      q05 <- unname(quantile(r, 0.05))
      threshold <- gate$min_ESS_threshold[1]
      n_iter_new <- min(max_n_iter, 50 * ceiling(threshold / q05 / n_chains / 50))
      data.frame(test = sub("red_fact.*", "", cell),
                 DGM = sub(".*DGM_(.*)N_500n_studies.*", "\\1", cell),
                 n_studies = as.integer(sub(".*n_studies_([0-9]+)param.*", "\\1", cell)),
                 model = sub(".*param_(.*)soft_.*", "\\1", cell),
                 random_thresholds = as.logical(sub(".*rand_thr_(TRUE|FALSE).*", "\\1", cell)),
                 prior_set = if (grepl("alpha_pi_one", cell)) "set_0" else sub(".*(set_[0-2]).*", "\\1", cell),
                 fits = nrow(gate),
                 threshold = threshold,
                 min_ESS_median = round(median(gate$min_ESS)),
                 n_iter_new = n_iter_new,
                 failed_fits = sum(failed),
                 fits_over_1pct_divergent = sum(gate$percent_divergent > 1),
                 divergent_pct_max = round(max(gate$percent_divergent), 3))

}
##
files_of_the_gate_tables <- unlist(lapply(folders_of_the_pilot,
                                          list.files,
                                          pattern = "^gate_of_every_fit_.*\\.csv$",
                                          full.names = TRUE))
## (one cell = every gate table of that name)
iterations_of_every_cell <- do.call(rbind, lapply(split(files_of_the_gate_tables, basename(files_of_the_gate_tables)),
                                                  fn_iterations_of_one_cell))
## (only cells with all their fits)
iterations_of_every_cell <- iterations_of_every_cell[iterations_of_every_cell$fits >= n_seeds, ]
iterations_of_every_cell <- iterations_of_every_cell[order(iterations_of_every_cell$prior_set,
                                                           iterations_of_every_cell$model,
                                                           iterations_of_every_cell$random_thresholds,
                                                           iterations_of_every_cell$DGM,
                                                           iterations_of_every_cell$test,
                                                           iterations_of_every_cell$n_studies), ]
options(width = 250)
print(iterations_of_every_cell, row.names = FALSE)
##
## ---- The settings table after the pilot: the iterations of each cell only --------------------------------------------------------
##
settings <- data.frame(test = iterations_of_every_cell$test,
                       DGM = iterations_of_every_cell$DGM,
                       n_studies = iterations_of_every_cell$n_studies,
                       model = iterations_of_every_cell$model,
                       random_thresholds = iterations_of_every_cell$random_thresholds,
                       prior_set = iterations_of_every_cell$prior_set,
                       settings =paste0("N_ITER_FOR_4_CHAINS=", iterations_of_every_cell$n_iter_new),
                       stringsAsFactors = FALSE)
write.csv(settings, settings_file_after_the_pilot, row.names = FALSE)
cat("\ncells in the settings table:", nrow(settings), "\nwritten:", settings_file_after_the_pilot, "\n")
