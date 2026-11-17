##
## ---- The realised proportion of missing thresholds, and the mean number of thresholds reported per study, in the simulated data
## ---- of the nine sets, read from the data saved with each kept fit (x: the non-diseased and diseased counts, -1 = missing)
##
## One group a test, proportion of thresholds missing and number of studies: up to <files a group> results files of the O-biv-FC
## model under prior set 0 (the data of a seed do not depend on the model), all DGMs. GAD-2 at ~15% and BAI at ~55% are read
## from the folder of the earlier sim (the real reporting pattern of each test).
##
## Usage: Rscript <this file> <sim_study folder> <files a group> <csv to write>
##
arguments <-  commandArgs(trailingOnly = TRUE)
S <-  arguments[1]
files_a_group <-  as.integer(arguments[2])
out <-  arguments[3]
set.seed(1)
factorial <-  function(p) file.path(S, sprintf("TEST_sim_in_the_factorial_design_of_missing_thresholds_15_Nov_with_%d_pct_of_thresholds_missing", p))
earlier <-  file.path(S, "whole_sim_sets_0_1_2_MCSE_0.125_cutpoint_bounds_minus_10_and_2.5")
groups <-  list()
for (n in c(10, 50)) {
      groups[[length(groups) + 1]] <-  list(test = "GAD_2", pct = 15, n = n, folder = earlier)
      groups[[length(groups) + 1]] <-  list(test = "BAI", pct = 55, n = n, folder = earlier)
      for (p in c(55, 85)) groups[[length(groups) + 1]] <-  list(test = "GAD_2", pct = p, n = n, folder = factorial(p))
      for (p in c(15, 55, 85)) groups[[length(groups) + 1]] <-  list(test = "HADS", pct = p, n = n, folder = factorial(p))
      for (p in c(15, 85)) groups[[length(groups) + 1]] <-  list(test = "BAI", pct = p, n = n, folder = factorial(p))
}
rows <-  NULL
for (g in groups) {
      pattern <-  sprintf("^seed_[0-9]+_test_%sred.*N_500n_studies_%dparam_ord_bivariatesoft_FALSErand_thr_FALSEset_0_results_list[.]RDS$", g$test, g$n)
      files <-  list.files(g$folder, pattern = pattern, full.names = TRUE)
      if (!length(files)) { cat("no files:", g$test, g$pct, g$n, "\n"); next }
      files <-  files[sample.int(length(files), min(files_a_group, length(files)))]
      missing <-  reported_per_study <-  fewest <-  numeric(0)
      for (f in files) {
            x <-  tryCatch(readRDS(f)$x[[1]], error = function(e) NULL)
            if (is.null(x)) next
            thresholds <-  as.matrix(x[, -1, drop = FALSE])
            missing <-  c(missing, mean(thresholds == -1))
            reported <-  rowSums(thresholds != -1)
            reported_per_study <-  c(reported_per_study, mean(reported))
            fewest <-  c(fewest, min(reported))
      }
      rows <-  rbind(rows, data.frame(test = g$test, pct = g$pct, n_studies = g$n, n_thresholds = ncol(thresholds), data_sets = length(missing),
                                      pct_missing_mean = round(100 * mean(missing), 2), pct_missing_min = round(100 * min(missing), 1), pct_missing_max = round(100 * max(missing), 1),
                                      thresholds_reported_per_study_mean = round(mean(reported_per_study), 3), pct_of_data_sets_where_every_study_reports_one = round(100 * mean(reported_per_study == 1), 1),
                                      fewest_reported_by_a_study = min(fewest)))
}
utils::write.csv(rows, out, row.names = FALSE)
print(rows, row.names = FALSE)
