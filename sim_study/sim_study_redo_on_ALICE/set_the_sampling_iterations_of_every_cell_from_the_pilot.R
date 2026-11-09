##
## ---- Set the sampling iterations of every cell from the 10-seeds-per-cell pilot ----------------------------------------------------
##
## (reads every gate_of_every_fit_*.csv of the pilot folders; for each cell the min ESS of Se and Sp per transition of its
##  fits that neither exceed 1% divergent nor had a stuck chain, r = min_ESS / n_transitions, is taken as log-normal and
##  its 5th percentile q05 = exp(mean(log r) - 1.645 sd(log r))
##  gives the transitions that put about 95% of the fits over the cell's ESS threshold: N_ITER_FOR_4_CHAINS = threshold /
##  q05 / 4, rounded up to 50; the per-cell settings table gets that value, and a row for every cell it lacks, with the
##  default settings of the launcher and the new iterations)
##
## ---- Controls ----------------------------------------------------------------------------------------------------------------------
##
arguments <- commandArgs(trailingOnly = TRUE)
root <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs/sim_study"
folders_of_the_pilot <- if (length(arguments) > 0) arguments else
                        file.path(root, c("pilot_10_seeds_per_cell_prior_mean_anchor_9_Nov",
                                          "pilot_10_seeds_per_cell_prior_mean_anchor_9_Nov_copy_from_ALICE"))
folder_of_the_tables <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/sim_study/sim_study_redo_on_ALICE"
settings_file <- file.path(folder_of_the_tables, "settings_of_every_cell_of_the_whole_sim.csv")
settings_file_after_the_pilot <- file.path(folder_of_the_tables, "settings_of_every_cell_of_the_whole_sim_after_the_pilot.csv")
n_chains <- 4
## (the launcher's defaults, for a cell the table lacks)
default_settings <- paste("TAU_INITIAL=pi LEARNING_RATE=0.05 CLIP_ITER=40 CLIP_ITER_TAU=120 DIFFUSION_HMC=1",
                          "TAU_JITTER_BURNIN=halton SHARE_TAU_II=1 ESJD_JUMP_POWER=2 METRIC_START_ITER=85",
                          "METRIC_UPDATE_RESCALES=0 EPS_REINIT_AT_HANDOVER=1 INITIAL_EPS_SEARCH_AT_EACH_CHAIN=1")
##
## ---- One row per cell --------------------------------------------------------------------------------------------------------------
##
fn_iterations_of_one_cell <- function(file_of_the_gate_table) {

      gate <- read.csv(file_of_the_gate_table)
      cell <- sub("gate_of_every_fit_test_", "", basename(file_of_the_gate_table))
      ## (a fit over 1% divergent fails the divergence gate whatever its length, and a fit whose min ESS is under a tenth
      ##  of the cell's median had a stuck chain; neither is a question of length, so both are left out of the ESS rate)
      failed <- gate$percent_divergent > 1 | gate$min_ESS < 0.1 * median(gate$min_ESS)
      r <- (gate$min_ESS / gate$n_transitions)[!failed]
      q05 <- exp(mean(log(r)) - 1.645 * sd(log(r)))
      threshold <- gate$min_ESS_threshold[1]
      n_iter_now <- round(gate$n_transitions[1] / n_chains)
      n_iter_new <- 50 * ceiling(threshold / q05 / n_chains / 50)
      data.frame(test = sub("red_fact.*", "", cell),
                 DGM = sub(".*DGM_(.*)N_500n_studies.*", "\\1", cell),
                 n_studies = as.integer(sub(".*n_studies_([0-9]+)param.*", "\\1", cell)),
                 model = sub(".*param_(.*)soft_.*", "\\1", cell),
                 random_thresholds = sub(".*rand_thr_(TRUE|FALSE).*", "\\1", cell),
                 fits = nrow(gate),
                 threshold = threshold,
                 min_ESS_lowest = round(min(gate$min_ESS)),
                 min_ESS_median = round(median(gate$min_ESS)),
                 passed_now = sum(gate$min_ESS >= threshold),
                 n_iter_now = n_iter_now,
                 n_iter_new = n_iter_new,
                 passed_at_new = sum(r * n_iter_new * n_chains >= threshold),
                 failed_fits = sum(failed),
                 divergent_pct_max = round(max(gate$percent_divergent), 3))

}
##
files_of_the_gate_tables <- unlist(lapply(folders_of_the_pilot,
                                          list.files,
                                          pattern = "^gate_of_every_fit_.*\\.csv$",
                                          full.names = TRUE))
iterations_of_every_cell <- do.call(rbind, lapply(files_of_the_gate_tables, fn_iterations_of_one_cell))
iterations_of_every_cell <- iterations_of_every_cell[order(iterations_of_every_cell$model,
                                                           iterations_of_every_cell$random_thresholds,
                                                           iterations_of_every_cell$DGM,
                                                           iterations_of_every_cell$test,
                                                           iterations_of_every_cell$n_studies), ]
options(width = 250)
print(iterations_of_every_cell, row.names = FALSE)
##
## ---- The settings table after the pilot --------------------------------------------------------------------------------------------
##
settings <- read.csv(settings_file, stringsAsFactors = FALSE)
settings$random_thresholds <- as.character(settings$random_thresholds)
for (i in seq_len(nrow(iterations_of_every_cell))) {

      cell <- iterations_of_every_cell[i, ]
      row <- which(settings$test == cell$test & settings$DGM == cell$DGM & settings$n_studies == cell$n_studies &
                   settings$model == cell$model & settings$random_thresholds == cell$random_thresholds)
      new_setting <- paste0("N_ITER_FOR_4_CHAINS=", cell$n_iter_new)
      if (length(row) == 1) {
            settings$settings[row] <- if (grepl("N_ITER_FOR_4_CHAINS=", settings$settings[row]))
                                            sub("N_ITER_FOR_4_CHAINS=[0-9]+", new_setting, settings$settings[row]) else
                                            paste(settings$settings[row], new_setting)
      } else {
            settings <- rbind(settings, data.frame(test = cell$test, DGM = cell$DGM, n_studies = cell$n_studies,
                                                   model = cell$model, random_thresholds = cell$random_thresholds,
                                                   settings = paste(default_settings, new_setting),
                                                   stringsAsFactors = FALSE))
      }

}
settings$random_thresholds <- as.logical(settings$random_thresholds)
write.csv(settings, settings_file_after_the_pilot, row.names = FALSE)
cat("\ncells in the pilot:", nrow(iterations_of_every_cell),
    " rows in the settings table after the pilot:", nrow(settings),
    "\nwritten:", settings_file_after_the_pilot, "\n")
