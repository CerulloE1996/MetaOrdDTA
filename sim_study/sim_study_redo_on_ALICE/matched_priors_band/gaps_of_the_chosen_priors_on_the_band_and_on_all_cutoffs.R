##
## ---- The chosen setting of each model, test and kappa prior (closest_setting_on_the_band.csv): the largest and
##      mean absolute gap to the reference over the band AND over every threshold of the test.
##
folder_of_the_fits <-  commandArgs(trailingOnly = TRUE)[1]
band_of_each_test <-  list(GAD_2 = 1:6, HADS = 3:17, BAI = 5:43)
all_of_each_test <-  list(GAD_2 = 1:6, HADS = 1:21, BAI = 1:63)
##
fn_pooled_quantiles_of_the_band <- function(files,
                                            band) {

      draws_of_the_files <-  lapply(files, function(file_of_the_fit) readRDS(file_of_the_fit)$draws_of_Se_Sp_and_kappa)
      quantiles <-  list()
      for (quantity in c("Se", "Sp")) for (threshold in band) {
            name <-  paste0(quantity, "_baseline[", threshold, "]")
            draws <-  unlist(lapply(draws_of_the_files, function(a) as.vector(a[, , name])))
            quantiles[[paste(quantity, threshold)]] <-  stats::quantile(draws, c(0.5, 0.025, 0.975))
      }
      return(quantiles)

}
fn_gaps <- function(this_prior,
                    reference) {

      gaps <-  c()
      names_of_the_gaps <-  c()
      for (name in names(reference)) {
            gaps <-  c(gaps, abs(this_prior[[name]] - reference[[name]]))
            names_of_the_gaps <-  c(names_of_the_gaps, paste(name, c("median", "lower", "upper")))
      }
      return(list(largest = max(gaps), where = names_of_the_gaps[which.max(gaps)], mean = mean(gaps)))

}
closest <-  utils::read.csv(file.path(folder_of_the_fits, "closest_setting_on_the_band.csv"), colClasses = "character")
rows <-  list()
for (i in seq_len(nrow(closest))) {
      a <-  closest[i, ]
      test <-  a$test
      files_of_the_reference <-  Sys.glob(file.path(folder_of_the_fits, paste0("xu_fixed_", test, "_kappa_current"),
                                                    paste0("xu_fixed_", test, "_SDs_*_seed_*.rds")))
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
      on_the_band <-  fn_gaps(fn_pooled_quantiles_of_the_band(files_of_the_setting, band_of_each_test[[test]]),
                              fn_pooled_quantiles_of_the_band(files_of_the_reference, band_of_each_test[[test]]))
      on_all <-  fn_gaps(fn_pooled_quantiles_of_the_band(files_of_the_setting, all_of_each_test[[test]]),
                         fn_pooled_quantiles_of_the_band(files_of_the_reference, all_of_each_test[[test]]))
      rows[[length(rows) + 1]] <-  data.frame(model = a$model, test = test, kappa_prior = a$kappa_prior,
                                              setting = setting, n_fits = length(files_of_the_setting),
                                              band_largest = round(on_the_band$largest, 4),
                                              band_mean = round(on_the_band$mean, 4), band_where = on_the_band$where,
                                              all_largest = round(on_all$largest, 4),
                                              all_mean = round(on_all$mean, 4), all_where = on_all$where)
}
result <-  do.call(rbind, rows)
utils::write.csv(result, file.path(folder_of_the_fits, "gaps_of_the_chosen_priors_on_the_band_and_on_all_cutoffs.csv"),
                 row.names = FALSE)
print(result[, c("model", "test", "kappa_prior", "n_fits", "band_largest", "band_mean", "all_largest", "all_mean",
                 "all_where")], row.names = FALSE)
print(result$setting)
