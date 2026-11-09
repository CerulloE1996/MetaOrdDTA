##
## ---- Summarise the gate tables of the 10-seeds-per-cell pilot ----------------------------------------------------------------------
##
## (one row per cell, from every gate_of_every_fit_*.csv of the folders given as arguments: fits, kept, min ESS of
##  Se and Sp, divergent transitions; the fused-model fits of the pilot of 9 Nov on ALICE and the local HPC)
##
## ---- Controls ----------------------------------------------------------------------------------------------------------------------
##
folders_of_the_pilot <- commandArgs(trailingOnly = TRUE)
if (length(folders_of_the_pilot) == 0) {
      root <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA_outputs/sim_study"
      folders_of_the_pilot <- file.path(root, c("pilot_10_seeds_per_cell_prior_mean_anchor_9_Nov",
                                                "pilot_10_seeds_per_cell_prior_mean_anchor_9_Nov_copy_from_ALICE"))
}
##
## ---- One row per cell --------------------------------------------------------------------------------------------------------------
##
fn_summarise_one_gate_table <- function(file_of_the_gate_table) {

      gate <- read.csv(file_of_the_gate_table)
      cell <- sub("gate_of_every_fit_test_", "", basename(file_of_the_gate_table))
      data.frame(test = sub("red_fact.*", "", cell),
                 DGM = sub(".*DGM_(.*)N_500n_studies.*", "\\1", cell),
                 n_studies = as.integer(sub(".*n_studies_([0-9]+)param.*", "\\1", cell)),
                 model = sub(".*param_(.*)soft_.*", "\\1", cell),
                 rand_thr = sub(".*rand_thr_(TRUE|FALSE).*", "\\1", cell),
                 fits = nrow(gate),
                 kept = sum(gate$kept),
                 ESS_threshold = gate$min_ESS_threshold[1],
                 min_ESS_lowest = round(min(gate$min_ESS)),
                 min_ESS_median = round(median(gate$min_ESS)),
                 under_threshold = sum(gate$min_ESS < gate$min_ESS_threshold),
                 divergent_total = sum(gate$n_divergent),
                 divergent_pct_max = round(max(gate$percent_divergent), 3),
                 over_1_pct_divergent = sum(gate$percent_divergent > 1),
                 seconds_per_fit = round(mean(gate$seconds_burnin + gate$seconds_sampling)))

}
##
## ---- The table ---------------------------------------------------------------------------------------------------------------------
##
files_of_the_gate_tables <- unlist(lapply(folders_of_the_pilot,
                                          list.files,
                                          pattern = "^gate_of_every_fit_.*\\.csv$",
                                          full.names = TRUE))
summary_of_every_cell <- do.call(rbind, lapply(files_of_the_gate_tables, fn_summarise_one_gate_table))
summary_of_every_cell <- summary_of_every_cell[order(summary_of_every_cell$model,
                                                     summary_of_every_cell$rand_thr,
                                                     summary_of_every_cell$DGM,
                                                     summary_of_every_cell$test,
                                                     summary_of_every_cell$n_studies), ]
options(width = 250)
print(summary_of_every_cell, row.names = FALSE)
cat("\ncells:", nrow(summary_of_every_cell),
    " fits:", sum(summary_of_every_cell$fits),
    " kept:", sum(summary_of_every_cell$kept),
    " under the ESS threshold:", sum(summary_of_every_cell$under_threshold),
    " with over 1% divergent:", sum(summary_of_every_cell$over_1_pct_divergent), "\n")
