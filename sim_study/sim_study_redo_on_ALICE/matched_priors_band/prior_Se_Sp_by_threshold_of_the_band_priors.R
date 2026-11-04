##
## ---- The chosen setting of each model, test and kappa prior (closest_setting_on_the_band.csv): the pooled prior
##      median, 2.5% and 97.5% quantiles of the summary Se and Sp at EVERY threshold of the test, for the figure.
##
folder_of_the_fits <-  commandArgs(trailingOnly = TRUE)[1]
all_of_each_test <-  list(GAD_2 = 1:6, HADS = 1:21, BAI = 1:63)
##
closest <-  utils::read.csv(file.path(folder_of_the_fits, "closest_setting_on_the_band.csv"), colClasses = "character")
rows <-  list()
for (i in seq_len(nrow(closest))) {
      a <-  closest[i, ]
      test <-  a$test
      folder <-  file.path(folder_of_the_fits, paste0(a$model, "_", test, "_kappa_", a$kappa_prior))
      files <-  Sys.glob(file.path(folder, "*.rds"))
      start <-  paste0(a$model, "_", test, "_SDs_", a$SD_1, "_", a$SD_2, "_", a$SD_3)
      if (a$mean_1 != "runner") start <-  paste0(start, "_means_", a$mean_1, "_", a$mean_2)
      files <-  files[startsWith(basename(files), start)]
      settings <-  unique(sub("_seed_[0-9]+\\.rds$", "", basename(files)))
      if (a$model == "jones") settings <-  settings[grepl(paste0("_", a$mean_3, "$"), settings) |
                                                    !grepl("_means_[^_]+_[^_]+_[^_]+$", settings)]
      stopifnot(length(settings) >= 1)
      setting <-  settings[1]
      files_of_the_setting <-  files[startsWith(basename(files), paste0(setting, "_seed_"))]
      draws_of_the_files <-  lapply(files_of_the_setting, function(file_of_the_fit) readRDS(file_of_the_fit)$draws_of_Se_Sp_and_kappa)
      for (quantity in c("Se", "Sp")) for (threshold in all_of_each_test[[test]]) {
            name <-  paste0(quantity, "_baseline[", threshold, "]")
            draws <-  unlist(lapply(draws_of_the_files, function(d) as.vector(d[, , name])))
            q <-  stats::quantile(draws, c(0.5, 0.025, 0.975))
            rows[[length(rows) + 1]] <-  data.frame(model = a$model, test = test, kappa_prior = a$kappa_prior,
                                                    setting = setting, quantity = quantity, threshold = threshold,
                                                    score_at_or_above = threshold,
                                                    n_seeds = length(files_of_the_setting), n_draws = length(draws),
                                                    median = q[[1]], lower = q[[2]], upper = q[[3]])
      }
}
result <-  do.call(rbind, rows)
utils::write.csv(result, file.path(folder_of_the_fits, "prior_Se_Sp_by_threshold_of_the_band_priors.csv"), row.names = FALSE)
print(table(result$model, result$kappa_prior))
