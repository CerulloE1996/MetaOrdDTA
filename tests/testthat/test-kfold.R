make_log_lik_draws <- function(offset = 0) {
  probabilities <- c(
    0.1, 0.9, 0.2, 0.8,
    0.3, 0.7, 0.4, 0.6
  )
  array(
    log(probabilities) + offset,
    dim = c(2, 2, 2),
    dimnames = list(
      iteration = NULL,
      chain = NULL,
      variable = c("log_lik_study[1]", "log_lik_study[2]")
    )
  )
}


make_kfold_result <- function(offset = 0, fold_assignments = c(1L, 2L)) {
  draws <- make_log_lik_draws(offset)
  list(
    outs_kfold = list(
      K = 2L,
      fold_assignments = fold_assignments,
      min_ESS_vec = c(500, 500),
      fold_results = list(
        list(
          fold = 1L,
          status = "completed",
          holdout_studies = 1L,
          log_lik_study_array = draws
        ),
        list(
          fold = 2L,
          status = "completed",
          holdout_studies = 2L,
          log_lik_study_array = draws
        )
      )
    )
  )
}


test_that("fold ELPD uses log posterior predictive density", {
  metrics <- MetaOrdDTA:::kfold_fold_metrics(make_log_lik_draws(), 1:2)

  expected <- c(
    log(mean(c(0.1, 0.9, 0.2, 0.8))),
    log(mean(c(0.3, 0.7, 0.4, 0.6)))
  )
  expect_equal(unname(metrics$pointwise_elpd), expected)
  expect_equal(metrics$elpd, sum(expected))
})


test_that("model comparison uses paired studywise differences", {
  comparison <- compare_models_k_fold_CV(
    list(make_kfold_result(), make_kfold_result(0.1)),
    model_names = c("reference", "shifted"),
    reference_model = "reference",
    verbose = FALSE
  )

  expect_equal(comparison$comparison_table$Delta_ELPD, c(0, 0.2))
  expect_equal(comparison$comparison_table$Delta_SE, c(0, 0))
})


test_that("model comparison defaults to the model with highest ELPD", {
  comparison <- compare_models_k_fold_CV(
    list(make_kfold_result(), make_kfold_result(0.1)),
    model_names = c("reference", "shifted"),
    verbose = FALSE
  )

  expect_identical(comparison$reference_model, "shifted")
  expect_equal(comparison$comparison_table$Delta_ELPD, c(-0.2, 0))
})


test_that("comparison reporting does not turn paired SE into a significance gate", {
  comparison <- compare_models_k_fold_CV(
    list(make_kfold_result(), make_kfold_result(0.1)),
    model_names = c("reference", "shifted"),
    verbose = FALSE
  )

  summary_text <- capture.output(summarise_model_comparison_k_fold_CV(comparison))
  latex_text <- capture.output(make_latex_table_k_fold_CV(comparison))

  expect_true(any(grepl("not a selection threshold", summary_text, fixed = TRUE)))
  expect_false(any(grepl("Significantly different", latex_text, fixed = TRUE)))
  expect_true(any(grepl("not a hypothesis-test threshold", latex_text, fixed = TRUE)))
})


test_that("different fold assignments cannot be compared", {
  expect_error(
    compare_models_k_fold_CV(
      list(make_kfold_result(), make_kfold_result(fold_assignments = c(2L, 1L))),
      verbose = FALSE
    ),
    "different fold assignments"
  )
})


test_that("fold creation is reproducible and balanced", {
  indicators <- rbind(
    c(1, 0), c(1, 0), c(1, 0),
    c(0, 1), c(0, 1), c(0, 1),
    c(1, 1), c(1, 1), c(1, 1)
  )
  first <- create_folds(indicators, K = 3, seed = 42, network = TRUE, n_studies = 9)
  second <- create_folds(indicators, K = 3, seed = 42, network = TRUE, n_studies = 9)

  expect_identical(first, second)
  expect_identical(sort(unique(first)), 1:3)
  expect_lte(max(table(first)) - min(table(first)), 1)
})


test_that("sequential and R6 interfaces accept the same CV controls", {
  sequential_args <- names(formals(R_fn_run_kfold_cv))
  r6_args <- names(formals(MetaOrd_model$public_methods$run_k_fold_CV))
  expected <- c("seed", "fold_assignments", "cmdstanr_args", "output_dir")

  expect_true(all(expected %in% sequential_args))
  expect_true(all(expected %in% r6_args))
})


test_that("NULL initial values remain valid", {
  expect_null(resize_init_list(NULL, 2))
})
