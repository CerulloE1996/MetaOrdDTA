# Paper Code Boxes 1-2: model fitting and K-fold cross-validation
#
# This script is safe to source as a quick smoke test. Set
# METAORDDTA_RUN_SAMPLING=true to compile and sample the demonstration model,
# and METAORDDTA_RUN_KFOLD=true to additionally run the full K-fold analysis.

library(MetaOrdDTA)

example_file <- system.file(
  "examples", "data_example_1_NMA_list.RDS",
  package = "MetaOrdDTA"
)
if (!nzchar(example_file)) {
  example_file <- file.path("inst", "examples", "data_example_1_NMA_list.RDS")
}
example_data <- readRDS(example_file)

# The bundled legacy object stores disease class first; the public API stores
# test first, then non-diseased/diseased matrices.
x <- lapply(seq_along(example_data$x_NMA[[1L]]), function(test) {
  list(
    example_data$x_NMA[[1L]][[test]],
    example_data$x_NMA[[2L]][[test]]
  )
})
indicator_index_test_in_study <- example_data$indicator_index_test_in_study

run_sampling <- identical(
  tolower(Sys.getenv("METAORDDTA_RUN_SAMPLING", "false")),
  "true"
)
run_kfold <- identical(
  tolower(Sys.getenv("METAORDDTA_RUN_KFOLD", "false")),
  "true"
)

model <- MetaOrd_model$new(
  x = x,
  indicator_index_test_in_study = indicator_index_test_in_study,
  intercept_only = TRUE,
  cov_data = NULL,
  network = TRUE,
  compound_symmetry = TRUE,
  prior_only = FALSE,
  model_parameterisation = "ord_bivariate",
  random_thresholds = FALSE,
  n_chains = 2,
  compile = run_sampling,
  quiet = TRUE
)

fold_assignments <- model$make_folds_k_fold_CV(seed = 123, K = 5)

if (run_sampling) {
  model$sample(
    seed = 123,
    n_chains = 2,
    n_burnin = 500,
    n_iter = 500,
    adapt_delta = 0.8
  )
  model$summary(
    compute_main_params = TRUE,
    compute_transformed_parameters = TRUE,
    compute_generated_quantities = TRUE,
    use_BayesMVP_for_faster_summaries = FALSE
  )

  print(model$get_summary_main(), n = 50)
  print(model$extract_params(params = "Se_baseline"), n = 20)
  print(model$extract_params(params = "Sp_baseline"), n = 20)
}

if (run_kfold) {
  if (!run_sampling) {
    stop("Set METAORDDTA_RUN_SAMPLING=true when running K-fold CV")
  }
  cv <- model$run_k_fold_CV(
    seed = 123,
    fold_assignments = fold_assignments,
    n_burnin = 500,
    n_iter = 500,
    n_chains = 2,
    adapt_delta = 0.8,
    max_treedepth = 10,
    parallel = TRUE,
    n_workers = 5
  )
  saveRDS(
    cv,
    file.path(
      metaorddta_output_path("cv_results", create = TRUE),
      "paper_example_fixed_cutpoints_CS.rds"
    )
  )
}
