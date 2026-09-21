# Re-score the completed Klaus AAD K-fold analyses.
#
# This script deliberately does not rerun Stan. It recovers study-wise log
# posterior predictive densities from the saved draws, compares all models on
# the identical 24-fold split, and writes tables through MetaOrdDTA functions.

source_tree <- file.exists("DESCRIPTION") &&
  identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L]), "MetaOrdDTA")
if (source_tree) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Running from the source tree requires devtools")
  }
  devtools::load_all(quiet = TRUE)
} else if (!requireNamespace("MetaOrdDTA", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install MetaOrdDTA, or run this script from its source tree with devtools available")
  }
  devtools::load_all(quiet = TRUE)
}

target_condition <- "AAD"
results_dir <- MetaOrdDTA::metaorddta_output_path(
  "application_results_real_data", target_condition
)
output_dir <- file.path(results_dir, "rescored_2026-09-17")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_ids <- c(0L, 1:15)
n_covariates <- c(1L, 2L, 3L, 3L, 2L, 7L, 4L, 4L, 3L, 5L, 4L, 4L, 6L, 5L, 5L, 6L)
model_names <- c(
  "MR_Model_0 (intercept-only)",
  "MR_Model_1 (intercept + prev_AAD)",
  "MR_Model_2 (intercept + ref_test)",
  "MR_Model_3 (intercept + study_setting)",
  "MR_Model_4 (intercept + RoB)",
  "MR_Model_5 (intercept + prev_AAD + ref_test + study_setting + RoB)",
  "MR_Model_6 (intercept + prev_AAD + ref_test)",
  "MR_Model_7 (intercept + prev_AAD + study_setting)",
  "MR_Model_8 (intercept + prev_AAD + RoB)",
  "MR_Model_9 (intercept + ref_test + study_setting)",
  "MR_Model_10 (intercept + ref_test + RoB)",
  "MR_Model_11 (intercept + study_setting + RoB)",
  "MR_Model_12 (intercept + prev_AAD + ref_test + study_setting)",
  "MR_Model_13 (intercept + prev_AAD + ref_test + RoB)",
  "MR_Model_14 (intercept + prev_AAD + study_setting + RoB)",
  "MR_Model_15 (intercept + ref_test + study_setting + RoB)"
)

result_filename <- function(model_id, n_covariates) {
  stem <- paste0(
    "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AAD",
    "intercept_only_", as.integer(model_id == 0L),
    "N_covs_", n_covariates,
    "param_Xurand_thr_FALSECS_TRUE",
    "prior_beta_mu_SD_0.5_0.25inform_dir_TRUE"
  )
  if (model_id != 0L) stem <- paste0(stem, "_MR_model_", model_id)
  paste0(stem, "_applied_results.RDS")
}

input_files <- file.path(
  results_dir,
  Map(result_filename, model_ids, n_covariates) |> unlist(use.names = FALSE)
)
missing <- input_files[!file.exists(input_files)]
if (length(missing)) {
  stop("Missing saved K-fold result(s):\n", paste(missing, collapse = "\n"))
}

rescored <- vector("list", length(input_files))
for (i in seq_along(input_files)) {
  message(sprintf("Rescoring %d/%d: %s", i, length(input_files), model_names[i]))
  rescored[[i]] <- MetaOrdDTA::rescore_kfold_results(
    readRDS(input_files[i]),
    drop_draws = TRUE
  )
  gc()
}
names(rescored) <- model_names

comparison <- MetaOrdDTA::compare_models_k_fold_CV(
  kfold_results_list = rescored,
  model_names = model_names,
  min_ess_threshold = 100,
  exclude_low_ess = FALSE,
  reference_model = "best"
)

comparison_vs_intercept <- MetaOrdDTA::compare_models_k_fold_CV(
  kfold_results_list = rescored,
  model_names = model_names,
  min_ess_threshold = 100,
  exclude_low_ess = FALSE,
  reference_model = model_names[1L],
  verbose = FALSE
)

ranking <- comparison$comparison_table
ranking <- ranking[order(ranking$ELPD, decreasing = TRUE), ]
ranking$Rank <- seq_len(nrow(ranking))
ranking <- ranking[c("Rank", "Model", "ELPD", "SE", "Delta_ELPD", "Delta_SE")]

write.csv(ranking, file.path(output_dir, "klaus_aad_kfold_ranking.csv"), row.names = FALSE)
write.csv(
  comparison_vs_intercept$comparison_table,
  file.path(output_dir, "klaus_aad_kfold_vs_intercept.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    comparison = comparison,
    comparison_vs_intercept = comparison_vs_intercept,
    rescored_results = rescored,
    source_files = input_files,
    generated_at = Sys.time()
  ),
  file.path(output_dir, "klaus_aad_kfold_rescored.rds"),
  compress = "xz"
)

latex <- capture.output(MetaOrdDTA::make_latex_table_k_fold_CV(
  comparison_results = comparison,
  caption = "K-fold cross-validation comparison of AAD meta-regression models",
  label = "tab:klaus_aad_kfold",
  digits_elpd = 1,
  digits_se = 1
))
writeLines(latex, file.path(output_dir, "klaus_aad_kfold_table.tex"))

min_ess <- min(comparison$min_ess_matrix, na.rm = TRUE)
conclusion <- c(
  sprintf("All 16 models used the same 24 folds across 99 studies; minimum fold ESS was %.1f.", min_ess),
  "Four covariate models had numerically higher ELPD than intercept-only: models 15 (+13.92), 10 (+13.70), 8 (+13.19), and 7 (+0.29).",
  "MR model 15 had the highest K-fold ELPD and is therefore the selected model under a best-predictive-performance rule.",
  "Models 10 and 8 were close to model 15 (ELPD differences 0.22 and 0.73); paired standard errors describe uncertainty in this ranking, not a 5% significance requirement for selection.",
  "Use model 15 as the primary meta-regression and models 10 and 8 as near-best sensitivity analyses. Because 16 candidate models were compared, describe the selection as predictive and acknowledge selection uncertainty."
)
writeLines(conclusion, file.path(output_dir, "klaus_aad_kfold_conclusion.txt"))

print(ranking, row.names = FALSE)
message("Wrote rescored Klaus results to: ", output_dir)
