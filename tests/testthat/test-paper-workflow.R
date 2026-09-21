test_that("paper NMA example initializes without compiling", {
  example_file <- system.file(
    "examples", "data_example_1_NMA_list.RDS",
    package = "MetaOrdDTA"
  )
  expect_true(nzchar(example_file))

  example_data <- readRDS(example_file)
  x <- lapply(seq_along(example_data$x_NMA[[1L]]), function(test) {
    list(
      example_data$x_NMA[[1L]][[test]],
      example_data$x_NMA[[2L]][[test]]
    )
  })

  capture.output(
    model <- MetaOrd_model$new(
      x = x,
      indicator_index_test_in_study = example_data$indicator_index_test_in_study,
      intercept_only = TRUE,
      cov_data = NULL,
      network = TRUE,
      compound_symmetry = TRUE,
      prior_only = FALSE,
      model_parameterisation = "ord_bivariate",
      random_thresholds = FALSE,
      n_chains = 2,
      compile = FALSE,
      quiet = TRUE
    )
  )

  expect_null(model$internal_obj$outs_stan_compile$stan_model_obj)
  expect_equal(model$internal_obj$outs_data$stan_data_list$n_studies, 10)
  capture.output(folds <- model$make_folds_k_fold_CV(seed = 123, K = 5))
  expect_length(folds, 10)
  expect_identical(sort(unique(folds)), 1:5)
})
