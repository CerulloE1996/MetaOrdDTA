
#  
# resize_init_list <- function(init_lists_per_chain, 
#                              n_chains_new) {
#         
#         n_chains_old <- length(init_lists_per_chain)
#         
#         if (n_chains_new == n_chains_old) {
#           # No change needed
#           return(init_lists_per_chain)
#         } else if (n_chains_new < n_chains_old) {
#           # Fewer chains - just take first n_chains_new
#           return(init_lists_per_chain[1:n_chains_new])
#         } else {
#           # More chains - repeat cyclically
#           indices <- rep(1:n_chains_old, length.out = n_chains_new)
#           return(init_lists_per_chain[indices])
#         }
#         
# }





create_folds_new <- function( indicator_index_test_in_study, 
                          K, 
                          seed, 
                          network, 
                          n_studies
) {
  
        set.seed(seed, kind = "Mersenne-Twister")
        
        if (!network) {
          study_test_matrix <- matrix(rep(1, n_studies), ncol = 1)
        } else {
          study_test_matrix <- indicator_index_test_in_study
        }
        
        test_patterns <- apply(study_test_matrix, 1, paste, collapse = "")
        
        folds <- rep(NA, n_studies)
        fold_counts <- rep(0, K)
        
        # Sort patterns by frequency (rare patterns first for better balance)
        pattern_freq <- table(test_patterns)
        patterns_sorted <- names(sort(pattern_freq))
        
        for (pattern in patterns_sorted) {
          
                idx <- which(test_patterns == pattern)
                shuffled_idx <- sample(idx)
                
                for (study in shuffled_idx) {
                  # Assign to the fold with fewest studies (greedy balancing)
                  target_fold <- which.min(fold_counts)
                  folds[study] <- target_fold
                  fold_counts[target_fold] <- fold_counts[target_fold] + 1
                }
          
        }
        
        return(folds)
  
}







#' Stratified sampling by test availability
#' create_folds
#' @keywords internal
#' @export
create_folds <- function(indicator_index_test_in_study,
                         K, 
                         seed,
                         network,
                         n_studies
) {
        if (length(K) != 1L || is.na(K) || K < 2L || K > n_studies) {
          stop("K must be between 2 and the number of studies")
        }

        if (isTRUE(network) &&
            (is.null(indicator_index_test_in_study) ||
             nrow(indicator_index_test_in_study) != n_studies)) {
          stop("indicator_index_test_in_study must have one row per study")
        }

        create_folds_new(
          indicator_index_test_in_study = indicator_index_test_in_study,
          K = as.integer(K),
          seed = seed,
          network = isTRUE(network),
          n_studies = n_studies
        )
  
}


validate_kfold_assignments <- function(fold_assignments, K, n_studies) {
  if (length(fold_assignments) != n_studies) {
    stop("fold_assignments length must equal the number of studies")
  }
  if (anyNA(fold_assignments) || any(fold_assignments != as.integer(fold_assignments))) {
    stop("fold_assignments must contain integer fold labels without missing values")
  }

  labels <- sort(unique(as.integer(fold_assignments)))
  if (!identical(labels, seq_len(as.integer(K)))) {
    stop("fold_assignments must use every integer label from 1 through K")
  }

  as.integer(fold_assignments)
}


log_mean_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) {
    return(maximum)
  }
  maximum + log(mean(exp(x - maximum)))
}


kfold_fold_metrics <- function(log_lik_study_array, holdout_studies) {
  draws <- as.array(log_lik_study_array)
  if (length(dim(draws)) != 3L) {
    stop("log_lik_study_array must have iteration, chain, and study dimensions")
  }

  variable_names <- dimnames(draws)[[3L]]
  study_columns <- holdout_studies
  if (!is.null(variable_names)) {
    expected <- sprintf("log_lik_study[%d]", holdout_studies)
    matched <- match(expected, variable_names)
    if (!anyNA(matched)) {
      study_columns <- matched
    }
  }

  holdout_draws <- draws[, , study_columns, drop = FALSE]
  pointwise <- vapply(
    seq_along(holdout_studies),
    function(i) log_mean_exp(as.numeric(holdout_draws[, , i])),
    numeric(1L)
  )
  names(pointwise) <- as.character(holdout_studies)

  list(
    elpd = sum(pointwise),
    pointwise_elpd = pointwise
  )
}


kfold_pointwise_from_results <- function(fold_results, n_studies) {
  pointwise <- rep(NA_real_, n_studies)

  for (fold_result in fold_results) {
    if (is.null(fold_result) || !identical(fold_result$status, "completed")) {
      next
    }

    fold_pointwise <- fold_result$pointwise_elpd
    if (is.null(fold_pointwise) && !is.null(fold_result$log_lik_study_array)) {
      fold_pointwise <- kfold_fold_metrics(
        fold_result$log_lik_study_array,
        fold_result$holdout_studies
      )$pointwise_elpd
    }
    if (is.null(fold_pointwise)) {
      stop("Fold results do not contain pointwise ELPD or log-likelihood draws")
    }

    pointwise[fold_result$holdout_studies] <- unname(fold_pointwise)
  }

  pointwise
}


kfold_elpd_summary <- function(pointwise_elpd) {
  observed <- pointwise_elpd[is.finite(pointwise_elpd)]
  list(
    elpd_kfold = sum(observed),
    se_elpd = if (length(observed) > 1L) {
      sqrt(length(observed) * stats::var(observed))
    } else {
      NA_real_
    },
    n_scored = length(observed)
  )
}


#' Recompute and optionally compact saved K-fold results
#'
#' Re-scores legacy results using study-wise log posterior predictive density.
#' Large per-fold draw arrays can be dropped after their pointwise scores have
#' been recovered.
#'
#' @param kfold_results A K-fold result object, optionally wrapped in a list
#'   named `outs_kfold`.
#' @param drop_draws If `TRUE`, discard large log-likelihood arrays and fold
#'   summary tibbles after rescoring.
#' @return The result object with corrected pointwise ELPD, total ELPD and SE.
#' @export
rescore_kfold_results <- function(kfold_results, drop_draws = TRUE) {
  wrapped <- !is.null(kfold_results$outs_kfold)
  results <- if (wrapped) kfold_results$outs_kfold else kfold_results
  n_studies <- length(results$fold_assignments)
  validate_kfold_assignments(results$fold_assignments, results$K, n_studies)

  pointwise <- kfold_pointwise_from_results(results$fold_results, n_studies)
  if (any(!is.finite(pointwise))) {
    stop("K-fold results have incomplete held-out predictions")
  }
  summary <- kfold_elpd_summary(pointwise)
  results$pointwise_elpd <- pointwise
  results$elpd_kfold <- summary$elpd_kfold
  results$se_elpd <- summary$se_elpd

  for (fold in seq_along(results$fold_results)) {
    fold_result <- results$fold_results[[fold]]
    if (!is.null(fold_result$log_lik_study_array)) {
      metrics <- kfold_fold_metrics(
        fold_result$log_lik_study_array,
        fold_result$holdout_studies
      )
      fold_result$pointwise_elpd <- metrics$pointwise_elpd
      fold_result$elpd_fold <- metrics$elpd
    }
    if (isTRUE(drop_draws)) {
      fold_result$log_lik_study_array <- NULL
      fold_result$tibble_all <- NULL
    }
    results$fold_results[[fold]] <- fold_result
  }

  if (wrapped) {
    kfold_results$outs_kfold <- results
    return(kfold_results)
  }
  results
}



#' R_fn_run_kfold_cv
#' @keywords internal
#' @export
R_fn_run_kfold_cv <- function(   debugging = FALSE,
                                 ##
                                 K = 10,
                                 seed = 123,
                                 fold_assignments,
                                 ##
                                 model_prep_obj,
                                 ##
                                 stan_data_list,
                                 ##
                                 cmdstanr_args = NULL,
                                 ##
                                 priors,
                                 init_lists_per_chain = NULL,
                                 ##
                                 n_burnin,
                                 n_iter,
                                 n_chains,
                                 adapt_delta,
                                 max_treedepth,
                                 ##
                                 output_dir = metaorddta_output_path("cv_results"),
                                 ##
                                 use_BayesMVP_for_faster_summaries = FALSE
) {
          cat("Using", K, "-fold CV (sequential)\n")

          n_studies <- stan_data_list$n_studies
          fold_ids <- validate_kfold_assignments(fold_assignments, K, n_studies)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          cmdstanr_args <- if (is.null(cmdstanr_args)) list() else cmdstanr_args

          fold_results <- vector("list", K)

          for (k in seq_len(K)) {
            holdout_studies <- which(fold_ids == k)
            cat("\nFold", k, "- holding out", length(holdout_studies), "studies\n")

            fold_priors <- priors
            fold_priors$K_fold_CV_indicator <- as.integer(fold_ids != k)
            fold_cmdstanr_args <- utils::modifyList(
              cmdstanr_args,
              list(
                seed = seed + k,
                output_dir = file.path(output_dir, sprintf("fold_%02d", k)),
                output_basename = paste0("fold_", k, "_chain")
              )
            )

            start_time <- Sys.time()
            model_samples_obj <- model_prep_obj$sample(
              seed = seed + k,
              n_burnin = n_burnin,
              n_iter = n_iter,
              adapt_delta = adapt_delta,
              max_treedepth = max_treedepth,
              metric_shape = "diag_e",
              priors = fold_priors,
              n_chains = n_chains,
              init_lists_per_chain = init_lists_per_chain,
              cmdstanr_args = fold_cmdstanr_args
            )
            model_summary_and_trace_obj <- model_samples_obj$summary(
              compute_main_params = TRUE,
              compute_transformed_parameters = FALSE,
              compute_generated_quantities = TRUE,
              save_log_lik_trace = TRUE,
              use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries
            )

            trace_gq <- model_summary_and_trace_obj$get_all_traces_list()$traces_as_arrays$trace_gq
            log_lik_study_array <- filter_param_draws_string(
              named_draws_array = trace_gq,
              param_string = "log_lik_study",
              condition = "exact_match",
              exclude_NA = TRUE
            )
            metrics <- kfold_fold_metrics(log_lik_study_array, holdout_studies)

            fold_results[[k]] <- list(
              fold = k,
              status = "completed",
              elpd = metrics$elpd,
              pointwise_elpd = metrics$pointwise_elpd,
              n_holdout = length(holdout_studies),
              holdout_studies = holdout_studies,
              duration_mins = as.numeric(difftime(Sys.time(), start_time, units = "mins"))
            )

            rm(model_samples_obj, model_summary_and_trace_obj, trace_gq, log_lik_study_array)
            gc()
          }

          pointwise_elpd <- kfold_pointwise_from_results(fold_results, n_studies)
          summary <- kfold_elpd_summary(pointwise_elpd)

          list(
            elpd_kfold = summary$elpd_kfold,
            se_elpd = summary$se_elpd,
            pointwise_elpd = pointwise_elpd,
            K = K,
            fold_results = fold_results,
            fold_assignments = fold_ids,
            errors = list()
          )
}










##
## ---- NicoStan's summary of one fold of the K-fold CV (NicoStan or cmdstanr) ----------------------------------------------------------
##
#' fn_K_fold_summaries_of_a_draws_array
#' @keywords internal
fn_K_fold_summaries_of_a_draws_array <- function(draws_array) {

        names_of_the_variables <- dimnames(draws_array)[[3]]
        tibble_all <-  MetaOrdDTA:::generate_summary_tibble( n_threads = NULL,
                                                             input_array = draws_array,
                                                             iter_dim = 1,
                                                             chain_dim = 2,
                                                             param_names = names_of_the_variables,
                                                             compute_nested_rhat = FALSE,
                                                             n_superchains = NULL)
        fn_rows_of_one_variable <-  function(name_of_the_variable) {

              return(tibble_all[grepl(paste0("^", name_of_the_variable, "(\\[|$)"), tibble_all$parameter), ])

        }
        is_log_lik_study <-  grepl("^log_lik_study(\\[|$)", names_of_the_variables)
        log_lik_study_array <-  draws_array[, , is_log_lik_study, drop = FALSE]
        log_lik_study_tibble <-  tibble_all[is_log_lik_study, ]
        ##
        return(list( log_lik_study_array  = log_lik_study_array,
                     tibble_all           = tibble_all,
                     beta_mu              = fn_rows_of_one_variable("beta_mu"),
                     beta_sigma           = fn_rows_of_one_variable("beta_sigma"),
                     beta_tau             = fn_rows_of_one_variable("beta_tau"),
                     ## (min ESS = min of the bulk, tail and theta-squared ESS)
                     min_ESS              = min(log_lik_study_tibble$n_eff, log_lik_study_tibble$n_eff_tail,
                                                log_lik_study_tibble$n_eff_sd, na.rm = TRUE),
                     max_rhat             = max(log_lik_study_tibble$Rhat, na.rm = TRUE)))

}




##
## ---- One fold of the K-fold CV sampled with NicoStan: the same outputs as a cmdstanr fold -----------------------------------------
##
#' R_fn_fit_one_K_fold_CV_fold_with_NicoStan
#' @keywords internal
R_fn_fit_one_K_fold_CV_fold_with_NicoStan  <-  function(  stan_data_list,
                                                          stan_model_file_name,
                                                          n_chains,
                                                          init_lists_per_chain,
                                                          seed,
                                                          n_burnin,
                                                          n_iter,
                                                          NicoStan_args) {

        ## (the transformed parameters too, e.g. beta_tau, as in a cmdstanr fold)
        NicoStan_args$compute_transformed_parameters <-  TRUE
        outs <-  MetaOrdDTA:::R_fn_sample_model_with_NicoStan( stan_data_list        = stan_data_list,
                                                               stan_model_file_name  = stan_model_file_name,
                                                               n_chains              = n_chains,
                                                               init_lists_per_chain  = init_lists_per_chain,
                                                               seed                  = seed,
                                                               n_burnin              = n_burnin,
                                                               n_iter                = n_iter,
                                                               NicoStan_args         = NicoStan_args)
        ## (iterations x chains x variables, with the transformed parameters and generated quantities)
        summaries <-  fn_K_fold_summaries_of_a_draws_array(unclass(outs$draws_array))
        ##
        return(list( log_lik_study_array  = summaries$log_lik_study_array,
                     tibble_all           = summaries$tibble_all,
                     beta_mu              = summaries$beta_mu,
                     beta_sigma           = summaries$beta_sigma,
                     beta_tau             = summaries$beta_tau,
                     min_ESS              = summaries$min_ESS,
                     max_rhat             = summaries$max_rhat,
                     divergences          = outs$divergences,
                     seconds_of_the_fit   = outs$seconds_of_the_fit,
                     so_file              = outs$library_used))

}



#' R_fn_run_kfold_cv_parallel
#' Run k-fold cross-validation in parallel across multiple folds
#' @importFrom foreach %dopar%
#' @export
R_fn_run_kfold_cv_parallel <- function(  debugging,
                                         ##
                                         K = 10, 
                                         seed = 123,
                                         fold_assignments,
                                         ##
                                         model_prep_obj,
                                         stan_data_list, 
                                         ##
                                         cmdstanr_args = NULL, ## (Stan folds only)
                                         ##
                                         priors,
                                         init_lists_per_chain,
                                         ##
                                         n_burnin,
                                         n_iter,
                                         n_chains,
                                         adapt_delta,
                                         max_treedepth,
                                         ##
                                         n_workers = NULL,
                                         output_dir = metaorddta_output_path("cv_results"),
                                         ##
                                         use_BayesMVP_for_faster_summaries,
                                         ##
                                         ## ("Stan": cmdstanr; "NicoStan": each fold with NicoStan, NicoStan_args as in
                                         ##  MetaOrd_model$sample)
                                         algorithm = "Stan",
                                         NicoStan_args = list()
) {
  
        if (!requireNamespace("doParallel", quietly = TRUE) ||
            !requireNamespace("foreach", quietly = TRUE)) {
          stop("Parallel K-fold CV requires the doParallel and foreach packages")
        }
        ##
        if (use_BayesMVP_for_faster_summaries &&
            !requireNamespace("BayesMVP", quietly = TRUE)) {
          stop("BayesMVP is required when use_BayesMVP_for_faster_summaries = TRUE")
        }
        ##
        cat("Using", K, "-fold CV (LOO inappropriate for NMA models) - PARALLEL VERSION\n")
        ##
        # Calculate optimal number of workers
        if (is.null(n_workers)) {
          # Use all available cores divided by chains per model
          n_workers <- min(K, max(1L, floor(parallel::detectCores() / n_chains)))
        }
        cat(sprintf("Running %d folds in parallel using %d workers (%d chains each = %d cores)\n", 
                    K, n_workers, n_chains, n_workers * n_chains))
        ##
        ## Create output directory if it doesn't exist:
        ##
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        # Create a minimal environment
        minimal_env <- new.env()
        ##
        # Set up fold assignments first (before parallel)
        minimal_env$n_studies <- stan_data_list$n_studies
        minimal_env$n_index_tests <- stan_data_list$n_index_tests
        minimal_env$n_thr_max <- stan_data_list$n_thr_max
        minimal_env$n_covariates_max <- stan_data_list$n_covariates_max
        minimal_env$intercept_only <- ifelse(minimal_env$n_covariates_max == 1, 1, 0)
        ##
        minimal_env$model_parameterisation <- model_prep_obj$advanced_model_options$model_parameterisation
        minimal_env$random_thresholds <- model_prep_obj$advanced_model_options$random_thresholds
        minimal_env$Dirichlet_random_effects_type <- model_prep_obj$advanced_model_options$Dirichlet_random_effects_type
        minimal_env$box_cox <- model_prep_obj$advanced_model_options$box_cox
        minimal_env$softplus <- model_prep_obj$advanced_model_options$softplus
        ##
        minimal_env$network <- model_prep_obj$basic_model_options$network
        minimal_env$cts <- model_prep_obj$basic_model_options$cts
        minimal_env$prior_only <- model_prep_obj$basic_model_options$prior_only
        ##
        minimal_env$cov_data <- model_prep_obj$cov_data
        minimal_env$x <- model_prep_obj$x
        minimal_env$indicator_index_test_in_study <- model_prep_obj$indicator_index_test_in_study
        # ##
        if (is.null(fold_assignments)) {
          stop("fold_assignments must be provided! Create them externally with create_folds()")
        }
        fold_ids <- validate_kfold_assignments(
          fold_assignments,
          K,
          stan_data_list$n_studies
        )
        ##
        minimal_env$fold_ids <- fold_ids
        ##
        current_wd <- getwd()
        ##
        minimal_env$debugging <- debugging
        ##
        minimal_env$K <- K
        minimal_env$seed <- seed
        ##
        minimal_env$model_prep_obj <- model_prep_obj
        minimal_env$stan_data_list <- stan_data_list
        ##
        minimal_env$priors <- priors
        minimal_env$init_lists_per_chain <- init_lists_per_chain
        ##
        minimal_env$n_burnin <- n_burnin
        minimal_env$n_iter <- n_iter
        minimal_env$n_chains <- n_chains
        minimal_env$adapt_delta <- adapt_delta
        minimal_env$max_treedepth <- max_treedepth
        ##
        minimal_env$fold_ids <- fold_ids
        minimal_env$current_wd <- current_wd
        minimal_env$output_dir <- output_dir
        ##
        minimal_env$use_BayesMVP_for_faster_summaries <- use_BayesMVP_for_faster_summaries
        ## (here too, as foreach also exports this function's own cmdstanr_args to the workers; NULL for NicoStan folds)
        cmdstanr_args <- if (is.null(cmdstanr_args)) list() else cmdstanr_args
        minimal_env$cmdstanr_args <- cmdstanr_args
        ##
        stopifnot(algorithm %in% c("Stan", "NicoStan"))
        minimal_env$algorithm <- algorithm
        minimal_env$NicoStan_args <- NicoStan_args
        minimal_env$R_fn_fit_one_K_fold_CV_fold_with_NicoStan <- R_fn_fit_one_K_fold_CV_fold_with_NicoStan
        minimal_env$fn_K_fold_summaries_of_a_draws_array <- fn_K_fold_summaries_of_a_draws_array
        ##
        # Set up parallel cluster
        cl <- parallel::makeCluster(n_workers, outfile = "")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        doParallel::registerDoParallel(cl)
        ##
        # # Export everything to all workers
        # clusterExport(cl,
        #               varlist = ls(),
        #               envir = environment())
        # Export necessary objects to workers
        minimal_env$resize_init_list <- MetaOrdDTA::resize_init_list
        minimal_env$filter_param_draws_string <- MetaOrdDTA::filter_param_draws_string
        ##
        parallel::clusterExport(cl,
                                ls(minimal_env),
                                envir = minimal_env)
        parallel::clusterEvalQ(cl, {
          if (identical(algorithm, "Stan")) library(cmdstanr) ## (cmdstanr only for Stan folds)
          library(MetaOrdDTA)
          library(posterior)
          library(dplyr)
          library(stringr)
        })
        ##
        k <- 2
        ##
        # Run folds in parallel
        fold_results_list <- foreach::foreach(k = 1:K,
                                     .combine = 'c',
                                     # .export = ls(envir = asNamespace("MetaOrdDTA")),
                                     .packages = c("MetaOrdDTA",
                                                   if (identical(algorithm, "Stan")) "cmdstanr",
                                                   "posterior",
                                                   "dplyr",
                                                   "stringr"),
                                     .errorhandling = "pass") %dopar% {
                                                             
                                             fold_result <- list(
                                               fold = k,
                                               status = "starting", 
                                               start_time = Sys.time()
                                             )
                                             
                                              tryCatch({
                                                     ##
                                                     setwd(current_wd)
                                                     ##
                                                     ## Load required packages
                                                     require(MetaOrdDTA)
                                                     if (identical(algorithm, "Stan")) require(cmdstanr)
                                                     require(posterior)
                                                     require(dplyr)
                                                     ##
                                                     if (use_BayesMVP_for_faster_summaries) { 
                                                       require(RcppParallel)
                                                       require(BayesMVP)
                                                     }
                                                     ##
                                                     # Create unique directories for this fold
                                                     fold_id <- sprintf("fold_%02d_pid_%s", k, Sys.getpid())
                                                     fold_output_dir <- file.path(output_dir, fold_id)
                                                     dir.create(fold_output_dir, showWarnings = FALSE, recursive = TRUE)
                                                     # ##
                                                     # Add this RIGHT AFTER creating fold_output_dir (around line 136)
                                                     # Create unique temp directory
                                                     fold_temp_dir <- file.path(tempdir(), fold_id)
                                                     dir.create(fold_temp_dir, showWarnings = FALSE, recursive = TRUE)
                                                     ##
                                                     # Save current environment
                                                     old_tmpdir <- Sys.getenv("TMPDIR")
                                                     old_temp   <- Sys.getenv("TEMP")
                                                     old_tmp    <- Sys.getenv("TMP")
                                                     ##
                                                     # Set new temp directories
                                                     Sys.setenv(TMPDIR = fold_temp_dir)
                                                     Sys.setenv(TEMP = fold_temp_dir)
                                                     Sys.setenv(TMP = fold_temp_dir)
                                                     ## 
                                                     # Also set Stan-specific temp
                                                     if (identical(algorithm, "Stan")) Sys.setenv(CMDSTAN_TEMP_DIR = fold_temp_dir)
                                                     ##
                                                     # init_lists_per_chain <- model_prep_obj$init_lists_per_chain
                                                     # # priors <- model_prep_obj$priors
                                                     # stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
                                                     ##
                                                     cat("Checking what's available:\n")
                                                     cat("stan_data_list names:", names(stan_data_list), "\n")
                                                     if ("n_index_tests" %in% names(stan_data_list)) {
                                                       n_index_tests <- stan_data_list$n_index_tests
                                                       cat("Set n_index_tests =", n_index_tests, "\n")
                                                     }
                                                     if (any(sapply(stan_data_list, is.null))) {
                                                       print(names(stan_data_list)[sapply(stan_data_list, is.null)])
                                                       stop("NULL values in stan_data_list")
                                                     }
                                                     ##
                                                     cat(sprintf("\n===== Worker starting fold %d/%d =====\n", k, K))
                                                     ##
                                                     ##
                                                     ## Get fold-specific information
                                                     holdout_studies <- which(fold_ids == k)
                                                     remaining_studies <- which(fold_ids != k)
                                                     n_holdout <- length(holdout_studies)
                                                     ##
                                                     cat(sprintf("Fold %d - holding out %d studies\n", k, n_holdout))
                                                     ##
                                                     ## Create fold-specific priors
                                                     priors$K_fold_CV_indicator <- as.integer(fold_ids != k)
                                                     ##
                                                     init_lists_per_chain <- resize_init_list(
                                                       init_lists_per_chain = init_lists_per_chain,
                                                       n_chains_new = n_chains
                                                     )
                                                     ##
                                                     priors$softplus <- softplus
                                                     ##
                                                     final_stan_args <- list(
                                                       seed = seed,
                                                       # data = stan_data_list,
                                                       init = init_lists_per_chain,
                                                       chains = n_chains,
                                                       parallel_chains = n_chains,
                                                       iter_sampling = n_iter,
                                                       iter_warmup = n_burnin,
                                                       max_treedepth = max_treedepth,
                                                       adapt_delta = adapt_delta,
                                                       metric = "diag_e"
                                                     )
                                                     ##
                                                     stan_data_list <- remove_duplicates_by_name(c(stan_data_list, priors))
                                                     ##
                                                     print(paste("compound_symmetry = ", stan_data_list$compound_symmetry))
                                                     ##
                                                     final_stan_args <- utils::modifyList(
                                                       final_stan_args,
                                                       cmdstanr_args
                                                     )
                                                     final_stan_args <- utils::modifyList(
                                                       final_stan_args,
                                                       list(
                                                         data = stan_data_list,
                                                         seed = seed + k,
                                                         output_dir = fold_output_dir,
                                                         output_basename = paste0("fold_", k, "_chain"),
                                                         save_latent_dynamics = FALSE
                                                       )
                                                     )
                                                     # ##
                                                     # final_stan_args <- modifyList(final_stan_args, cmdstanr_args)
                                                     # ##
                                                     if (identical(algorithm, "NicoStan")) {
                                                       ##
                                                       ## ---- The fold sampled with NicoStan (the model's NicoStan version):
                                                       ##
                                                       outs_of_the_fold <- R_fn_fit_one_K_fold_CV_fold_with_NicoStan(
                                                         stan_data_list = stan_data_list,
                                                         stan_model_file_name =
                                                           model_prep_obj$internal_obj$outs_stan_model_name$stan_model_file_name,
                                                         n_chains = n_chains,
                                                         init_lists_per_chain = init_lists_per_chain,
                                                         seed = seed + k,
                                                         n_burnin = n_burnin,
                                                         n_iter = n_iter,
                                                         NicoStan_args = NicoStan_args)
                                                       log_lik_study_array <- outs_of_the_fold$log_lik_study_array
                                                       tibble_all <- outs_of_the_fold$tibble_all
                                                       beta_mu <- outs_of_the_fold$beta_mu
                                                       beta_sigma <- outs_of_the_fold$beta_sigma
                                                       beta_tau <- outs_of_the_fold$beta_tau
                                                       min_ESS <- outs_of_the_fold$min_ESS
                                                       max_rhat <- outs_of_the_fold$max_rhat
                                                       fold_result$divergences <- outs_of_the_fold$divergences
                                                       fold_result$seconds_of_the_fit <- outs_of_the_fold$seconds_of_the_fit
                                                       fold_result$so_file <- outs_of_the_fold$so_file
                                                     } else {
                                                     stan_mod_samples <- do.call(model_prep_obj$internal_obj$outs_stan_compile$stan_model_obj$sample, 
                                                                                 final_stan_args)
                                                     ##
                                                     summaries <- fn_K_fold_summaries_of_a_draws_array(
                                                       unclass(stan_mod_samples$draws(format = "array")))
                                                     log_lik_study_array <- summaries$log_lik_study_array
                                                     tibble_all <- summaries$tibble_all
                                                     beta_mu <- summaries$beta_mu
                                                     beta_sigma <- summaries$beta_sigma
                                                     beta_tau <- summaries$beta_tau
                                                     min_ESS  <- summaries$min_ESS
                                                     max_rhat <- summaries$max_rhat
                                                     }
                                                     ##
                                                     metrics <- MetaOrdDTA:::kfold_fold_metrics(
                                                       log_lik_study_array,
                                                       holdout_studies
                                                     )
                                                     elpd_holdout <- metrics$elpd
                                                     ##
                                                     ### Update fold result
                                                     fold_result$status <- "completed"
                                                     fold_result$elpd <- elpd_holdout
                                                     fold_result$pointwise_elpd <- metrics$pointwise_elpd
                                                     ##
                                                     fold_result$tibble_all <- tibble_all
                                                     fold_result$beta_mu <- beta_mu
                                                     fold_result$beta_sigma <- beta_sigma
                                                     fold_result$beta_tau <- beta_tau
                                                     fold_result$log_lik_study_array <- log_lik_study_array
                                                     ##
                                                     fold_result$min_ESS <- min_ESS
                                                     fold_result$max_rhat <- max_rhat
                                                     ##
                                                     fold_result$n_holdout <- n_holdout
                                                     fold_result$holdout_studies <- holdout_studies
                                                     fold_result$end_time <- Sys.time()
                                                     fold_result$duration_mins <- as.numeric(difftime(fold_result$end_time, 
                                                                                                      fold_result$start_time, 
                                                                                                      units = "mins"))
                                                     ##
                                                     cat(sprintf("\n===== Fold %d completed in %.1f minutes (elpd = %.2f) =====\n", 
                                                                 k, fold_result$duration_mins, elpd_holdout))
                                                     # ##
                                                     # ## Optionally save fold-specific results:
                                                     # ##
                                                     # saveRDS(fold_result,
                                                     #         file = file.path(output_dir, sprintf("fold_%02d_results.rds", k)))
                                                     ##
                                                     ## ADD THIS - delete the massive CSV files:
                                                     ##
                                                     csv_files <- list.files(fold_output_dir, pattern = "\\.csv$", full.names = TRUE)
                                                     file.remove(csv_files)
                                                     unlink(fold_output_dir, recursive = TRUE)
                                               
                                             }, error = function(e) {
                                               cat(sprintf("\nFold %d ERROR: %s\n", k, e$message))
                                               fold_result$status <- "error"
                                               fold_result$error_message <- e$message
                                               fold_result$error_call <- toString(e$call)

                                               # CLEANUP even on error
                                               if (exists("old_tmpdir")) {
                                                 Sys.setenv(TMPDIR = old_tmpdir)
                                                 Sys.setenv(TEMP = old_temp)
                                                 Sys.setenv(TMP = old_tmp)
                                               }
                                               if (exists("fold_temp_dir")) {
                                                 unlink(fold_temp_dir, recursive = TRUE, force = TRUE)
                                               }
                                             }, finally = {
                                               # ALWAYS cleanup
                                               if (exists("old_tmpdir")) {
                                                 Sys.setenv(TMPDIR = old_tmpdir)
                                                 Sys.setenv(TEMP = old_temp)
                                                 Sys.setenv(TMP = old_tmp)
                                               }
                                               if (exists("fold_temp_dir")) {
                                                 unlink(fold_temp_dir, recursive = TRUE, force = TRUE)
                                               }
                                             })
                                             
                                             # fold_result
                                             list(fold_result)
                                     }
  
  
        # Stop workers before aggregating the results.
        parallel::stopCluster(cl)
        on.exit(NULL, add = FALSE)
        
        # Clean up temporary files
        for (k in 1:K) {
          fold_file <- file.path(output_dir, paste0("cv_environment_fold_", k, ".RData"))
          if (file.exists(fold_file)) {
            file.remove(fold_file)
          }
        }
             
        # Process results, retaining fold positions even when a worker failed.
        fold_results <- vector("list", K)
        for (result in fold_results_list) {
          if (is.list(result) && !is.null(result$fold)) {
            fold_results[[result$fold]] <- result
          }
        }
        successful_folds <- Filter(
          function(x) !is.null(x) && identical(x$status, "completed"),
          fold_results
        )
        
        if (length(successful_folds) < K) {
          warning(sprintf("Only %d/%d folds completed successfully", 
                          length(successful_folds), K))
        }
        
        pointwise_elpd <- kfold_pointwise_from_results(fold_results, stan_data_list$n_studies)
        elpd_summary <- kfold_elpd_summary(pointwise_elpd)
        elpd_kfold <- elpd_summary$elpd_kfold
        se_elpd <- elpd_summary$se_elpd
        
        # Summary of timing
        total_time <- sum(sapply(successful_folds, function(x) x$duration_mins))
        parallel_time <- if (length(successful_folds)) {
          max(sapply(successful_folds, function(x) x$duration_mins))
        } else {
          NA_real_
        }
        ##
        min_ESS_vec <- max_rhat_vec <- c()
        min_ESS_vec <- max_rhat_vec <- rep(NA_real_, K)
        for (result in successful_folds) {
          min_ESS_vec[result$fold] <- result$min_ESS
          max_rhat_vec[result$fold] <- result$max_rhat
        }
        
        cat("\n===== K-FOLD CV SUMMARY =====\n")
        cat(sprintf("Successful folds: %d/%d\n", length(successful_folds), K))
        cat(sprintf("Total elpd_kfold: %.2f (SE: %.2f)\n", elpd_kfold, se_elpd))
        cat(sprintf("Total CPU time: %.1f minutes\n", total_time))
        cat(sprintf("Wall clock time: %.1f minutes\n", parallel_time))
        cat(sprintf("Speedup factor: %.1fx\n", total_time / parallel_time))
        ##
        print(paste("min_ESS:"))
        print(quantile(min_ESS_vec))
        ##
        print(paste("max_Rhat:"))
        print(quantile(max_rhat_vec))
        ##
        # Return results compatible with original function
        return(list(
          elpd_kfold = elpd_kfold,
          se_elpd = se_elpd,
          pointwise_elpd = pointwise_elpd,
          K = K,
          ##
          min_ESS_vec = min_ESS_vec,
          max_rhat_vec = max_rhat_vec,
          ##
          fold_results = fold_results,
          fold_assignments = fold_ids,
          n_workers = n_workers,
          timing = list(
            total_cpu_mins = total_time,
            wall_clock_mins = parallel_time,
            speedup = total_time / parallel_time
          ),
          errors = Filter(function(x) !is.null(x) && identical(x$status, "error"),
                          fold_results)
        ))
  
}















# R_fn_run_kfold_cv_parallel_BayesMVP <- function( debugging,
#                                                  ##
#                                                  K = 10, 
#                                                  seed = 123,
#                                                  ##
#                                                  model_prep_obj,
#                                                  Stan_data_list, 
#                                                  ##
#                                                  cmdstanr_args,
#                                                  ##
#                                                  priors,
#                                                  init_lists_per_chain,
#                                                  ##
#                                                  n_burnin,
#                                                  n_iter,
#                                                  ##
#                                                  n_chains_burnin,
#                                                  n_chains_sampling,
#                                                  ##
#                                                  adapt_delta,
#                                                  max_treedepth,
#                                                  ##
#                                                  n_workers = NULL,
#                                                  output_dir = metaorddta_output_path("cv_results"),
#                                                  ##
#                                                  use_BayesMVP_for_faster_summaries,
#                                                  stanc_args
# ) {
#   
#         require(doParallel)
#         require(foreach)
#         ##
#         if (use_BayesMVP_for_faster_summaries) { 
#           require(RcppParallel)
#           require(BayesMVP)
#         }
#         ##
#         cat("Using", K, "-fold CV (LOO inappropriate for NMA models) - PARALLEL VERSION\n")
#         
#         # Calculate optimal number of workers
#         if (is.null(n_workers)) {
#           # Use all available cores divided by chains per model
#           n_workers <- min(K, floor(parallel::detectCores() / n_chains_sampling))
#         }
#         cat(sprintf("Running %d folds in parallel using %d workers (%d chains each = %d cores)\n", 
#                     K, n_workers, n_chains_sampling, n_workers * n_chains_sampling))
#         
#         # Create output directory if it doesn't exist
#         if (!dir.exists(output_dir)) {
#           dir.create(output_dir, recursive = TRUE)
#         }
#         # Create a minimal environment
#         minimal_env <- new.env()
#         ##
#         # Set up fold assignments first (before parallel)
#         minimal_env$n_studies <- Stan_data_list$n_studies
#         minimal_env$n_index_tests <- Stan_data_list$n_index_tests
#         minimal_env$n_thr_max <- Stan_data_list$n_thr_max
#         ##
#         minimal_env$indicator_index_test_in_study <- Stan_data_list$indicator_index_test_in_study
#         ##
#         fold_ids <-  create_folds(K = K, 
#                                   study_test_matrix = Stan_data_list$indicator_index_test_in_study,
#                                   seed = seed)
#         minimal_env$fold_ids <- fold_ids
#         ##
#         current_wd <- getwd()
#         ##
#         minimal_env$debugging <- debugging
#         ##
#         minimal_env$K <- K
#         minimal_env$seed <- seed
#         ##
#         minimal_env$model_prep_obj <- model_prep_obj
#         ##
#         minimal_env$Stan_data_list <- Stan_data_list
#         ##
#         minimal_env$priors <- priors
#         minimal_env$init_lists_per_chain <- init_lists_per_chain
#         ##
#         minimal_env$n_burnin <- n_burnin
#         minimal_env$n_iter <- n_iter
#         minimal_env$n_chains <- n_chains
#         minimal_env$adapt_delta <- adapt_delta
#         minimal_env$max_treedepth <- max_treedepth
#         ##
#         minimal_env$fold_ids <- fold_ids
#         minimal_env$current_wd <- current_wd
#         minimal_env$output_dir <- output_dir
#         ##
#         minimal_env$use_BayesMVP_for_faster_summaries <- use_BayesMVP_for_faster_summaries
#         ##
#         minimal_env$n_chains_sampling <- n_chains_sampling
#         minimal_env$n_chains_burnin <- n_chains_burnin
#         ##
#         minimal_env$stanc_args <- stanc_args
#         ##
#         # Set up parallel cluster
#         cl <- makeCluster(n_workers, 
#                           outfile = "")
#         registerDoParallel(cl)
#         ##
#         minimal_env$resize_init_list <- MetaOrdDTA::resize_init_list
#         minimal_env$filter_param_draws_string <- MetaOrdDTA::filter_param_draws_string
#         ##
#         clusterExport(cl, 
#                       ls(minimal_env), 
#                       envir = minimal_env)
#         clusterEvalQ(cl, {
#           library(cmdstanr)
#           library(MetaOrdDTA)
#           library(posterior)
#           library(dplyr)
#           library(stringr)
#         })
#         ##
#         k <- 2
#         ##
#         # Run folds in parallel
#         fold_results_list <- foreach(k = 1:K,
#                                      .combine = 'c',
#                                      # .export = ls(envir = asNamespace("MetaOrdDTA")),
#                                      .packages = c("MetaOrdDTA",
#                                                    "cmdstanr", 
#                                                    "bridgestan", 
#                                                    "RcppParallel",
#                                                    "BayesMVPtest",
#                                                    "posterior",
#                                                    "dplyr",
#                                                    "stringr"),
#                                      .errorhandling = "pass") %dopar% {
#                                        
#                                        fold_result <- list(
#                                          fold = k,
#                                          status = "starting", 
#                                          start_time = Sys.time()
#                                        )
#                                        
#                                        tryCatch({
#                                          
#                                          ##
#                                          setwd(current_wd)
#                                          ##
#                                          ## Load required packages
#                                          require(MetaOrdDTA)
#                                          require(cmdstanr)
#                                          require(posterior)
#                                          require(dplyr)
#                                          ##
#                                          if (use_BayesMVP_for_faster_summaries) { 
#                                            require(RcppParallel)
#                                            require(BayesMVPtest)
#                                          }
#                                          ## set bs environment variable (otherwise it'll try downloading it even if already installed...)
#                                          bs_path <- bridgestan_path()
#                                          Sys.setenv(BRIDGESTAN = bs_path)
#                                          ##
#                                          # Create unique directories for this fold
#                                          fold_id <- sprintf("fold_%02d_pid_%s", k, Sys.getpid())
#                                          fold_output_dir <- file.path(output_dir, fold_id)
#                                          dir.create(fold_output_dir, showWarnings = FALSE, recursive = TRUE)
#                                          ##
#                                          # Add this RIGHT AFTER creating fold_output_dir (around line 136)
#                                          # Create unique temp directory
#                                          fold_temp_dir <- file.path(tempdir(), fold_id)
#                                          dir.create(fold_temp_dir, showWarnings = FALSE, recursive = TRUE)
#                                          ##
#                                          # Save current environment
#                                          old_tmpdir <- Sys.getenv("TMPDIR")
#                                          old_temp <- Sys.getenv("TEMP")
#                                          old_tmp <- Sys.getenv("TMP")
#                                          ##
#                                          # Set new temp directories
#                                          Sys.setenv(TMPDIR = fold_temp_dir)
#                                          Sys.setenv(TEMP = fold_temp_dir)
#                                          Sys.setenv(TMP = fold_temp_dir)
#                                          ##
#                                          # Also set Stan-specific temp
#                                          Sys.setenv(CMDSTAN_TEMP_DIR = fold_temp_dir)
#                                          ##
#                                          cat(sprintf("\n===== Worker starting fold %d/%d =====\n", k, K))
#                                          ##
#                                          init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain, 
#                                                                                    n_chains_new = n_chains_burnin)
#                                          ##
#                                          ## Create fold-specific priors:
#                                          ##
#                                          Stan_data_list$K_fold_CV_indicator <- as.integer(fold_ids != k)
#                                          ##
#                                          ## Get fold-specific information:
#                                          ##
#                                          holdout_studies <- which(fold_ids == k)
#                                          remaining_studies <- which(fold_ids != k)
#                                          n_holdout <- length(holdout_studies)
#                                          ##
#                                          cat(sprintf("Fold %d - holding out %d studies\n", k, n_holdout))
#                                          ##
#                                          partitioned_HMC <- TRUE ;    diffusion_HMC <- FALSE
#                                          metric_shape_main <- "diag" ; metric_type_main <- "Empirical"
#                                          ##
#                                          model_samples <-  model_prep_obj$sample(  Stan_data_list = Stan_data_list,
#                                                                                    init_lists_per_chain = init_lists_per_chain,
#                                                                                    n_chains_burnin = n_chains_burnin,
#                                                                                    ##
#                                                                                    partitioned_HMC = partitioned_HMC,
#                                                                                    diffusion_HMC = diffusion_HMC,
#                                                                                    seed = seed,
#                                                                                    n_burnin = n_burnin,
#                                                                                    n_iter = n_iter,
#                                                                                    ##
#                                                                                    n_chains_sampling = n_chains_sampling,
#                                                                                    n_superchains = n_superchains,
#                                                                                    ##
#                                                                                    adapt_delta = 0.80,
#                                                                                    learning_rate = 0.075,
#                                                                                    ##
#                                                                                    metric_shape_main = metric_shape_main,
#                                                                                    metric_type_main = metric_type_main,
#                                                                                    ##
#                                                                                    tau_mult = 1.6,
#                                                                                    ##
#                                                                                    clip_iter = 50,
#                                                                                    clip_iter_tau = 150,
#                                                                                    ##
#                                                                                    interval_width_main = 10,
#                                                                                    interval_width_nuisance = 10,
#                                                                                    ratio_M_us = 0.50,
#                                                                                    ratio_M_main = 0.50,
#                                                                                    ##
#                                                                                    beta1_adam = 0.0,
#                                                                                    beta2_adam = 0.95,
#                                                                                    eps_adam = 1e-8,
#                                                                                    ##
#                                                                                    parallel_method = "RcppParallel",
#                                                                                    stanc_args = stanc_args)
#                                          ##
#                                          model_fit <- model_samples$summary(save_log_lik_trace = FALSE, 
#                                                                             ##
#                                                                             compute_nested_rhat = FALSE,
#                                                                             ##
#                                                                             compute_main_params = TRUE,
#                                                                             compute_transformed_parameters = FALSE,
#                                                                             compute_generated_quantities = TRUE)
#       
#                                          trace_gq <- model_fit$traces$traces_as_arrays$trace_generated_quantities
#                                          ##
#                                          log_lik_study_array <- extract_params_from_array_batch(debugging = FALSE,
#                                                                                                 array_3d = trace_gq, 
#                                                                                                 combine = TRUE,
#                                                                                                 param_strings_vec = c("log_lik_study"),
#                                                                                                 condition = "containing") 
#                                          ##
#                                          
#                                          summary_tibbles <- model_fit$summaries$summary_tibbles
#                                          ##
#                                          log_lik_study_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
#                                                                                        tibble = tibble_gq, 
#                                                                                        param_strings_vec = c("log_lik_study"),
#                                                                                        condition = "containing") %>% print(n=100)
#                                          ##
#                                          log_lik_study_tibble <- stan_mod_samples$summary(variables = c("log_lik_study"))
#                                          min_ESS <- min(log_lik_study_tibble$n_eff, na.rm = TRUE)
#                                          max_rhat <- max(log_lik_study_tibble$Rhat, na.rm = TRUE)
#                                          ##
#                                          ## Calculate holdout log-likelihood
#                                          log_lik_holdout_samples <- log_lik_study_array[, , holdout_studies, drop = FALSE]
#                                          ##
#                                          ## Sum across holdout studies and average across iterations/chains
#                                          elpd_holdout <- mean(apply(log_lik_holdout_samples, c(1, 2), sum))
#                                          ##
#                                          ### Update fold result
#                                          fold_result$summary_tibbles <- summary_tibbles
#                                          ##
#                                          fold_result$status <- "completed"
#                                          fold_result$elpd <- elpd_holdout
#                                          ##
#                                          fold_result$min_ESS <- min_ESS
#                                          fold_result$max_rhat <- max_rhat
#                                          ##
#                                          fold_result$n_holdout <- n_holdout
#                                          fold_result$holdout_studies <- holdout_studies
#                                          fold_result$end_time <- Sys.time()
#                                          fold_result$duration_mins <- as.numeric(difftime(fold_result$end_time, 
#                                                                                           fold_result$start_time, 
#                                                                                           units = "mins"))
#                                          ##
#                                          cat(sprintf("\n===== Fold %d completed in %.1f minutes (elpd = %.2f) =====\n", 
#                                                      k, fold_result$duration_mins, elpd_holdout))
#                                          ##
#                                          ## Optionally save fold-specific results
#                                          saveRDS(fold_result, 
#                                                  file = file.path(output_dir, sprintf("fold_%02d_results.rds", k)))
#                                          
#                                        }, error = function(e) {
#                                          cat(sprintf("\nFold %d ERROR: %s\n", k, e$message))
#                                          fold_result$status <- "error"
#                                          fold_result$error_message <- e$message
#                                          fold_result$error_call <- toString(e$call)
#                                          
#                                          # CLEANUP even on error
#                                          if (exists("old_tmpdir")) {
#                                            Sys.setenv(TMPDIR = old_tmpdir)
#                                            Sys.setenv(TEMP = old_temp)
#                                            Sys.setenv(TMP = old_tmp)
#                                          }
#                                          if (exists("fold_temp_dir")) {
#                                            unlink(fold_temp_dir, recursive = TRUE, force = TRUE)
#                                          }
#                                        }, finally = {
#                                          # ALWAYS cleanup
#                                          if (exists("old_tmpdir")) {
#                                            Sys.setenv(TMPDIR = old_tmpdir)
#                                            Sys.setenv(TEMP = old_temp)
#                                            Sys.setenv(TMP = old_tmp)
#                                          }
#                                          if (exists("fold_temp_dir")) {
#                                            unlink(fold_temp_dir, recursive = TRUE, force = TRUE)
#                                          }
#                                        })
#                                        
#                                        # fold_result
#                                        list(fold_result)
#                                      }
#         
#         
#         # Stop cluster
#         stopCluster(cl)
#         
#         # Clean up temporary files
#         for (k in 1:K) {
#           fold_file <- file.path(output_dir, paste0("cv_environment_fold_", k, ".RData"))
#           if (file.exists(fold_file)) {
#             file.remove(fold_file)
#           }
#         }
#         
#         # Process results
#         successful_folds <- Filter(function(x) !is.null(x) && x$status == "completed", 
#                                    fold_results_list)
#         
#         if (length(successful_folds) < K) {
#           warning(sprintf("Only %d/%d folds completed successfully", 
#                           length(successful_folds), K))
#         }
#         
#         # Calculate total elpd_kfold
#         elpd_kfold <- sum(sapply(successful_folds, function(x) x$elpd))
#         
#         # Compute SE using fold-level results
#         fold_elpds <- sapply(successful_folds, function(x) x$elpd)
#         se_elpd <- sd(fold_elpds) * sqrt(K)
#         
#         # Summary of timing
#         total_time <- sum(sapply(successful_folds, function(x) x$duration_mins))
#         parallel_time <- max(sapply(successful_folds, function(x) x$duration_mins))
#         ##
#         min_ESS_vec <- max_rhat_vec <- c()
#         for (k in 1:length(successful_folds)) {
#           min_ESS_vec[k]  <- successful_folds[[k]]$min_ESS
#           max_rhat_vec[k] <- successful_folds[[k]]$max_rhat
#         }
#         
#         cat("\n===== K-FOLD CV SUMMARY =====\n")
#         cat(sprintf("Successful folds: %d/%d\n", length(successful_folds), K))
#         cat(sprintf("Total elpd_kfold: %.2f (SE: %.2f)\n", elpd_kfold, se_elpd))
#         cat(sprintf("Total CPU time: %.1f minutes\n", total_time))
#         cat(sprintf("Wall clock time: %.1f minutes\n", parallel_time))
#         cat(sprintf("Speedup factor: %.1fx\n", total_time / parallel_time))
#         ##
#         print(paste("min_ESS:"))
#         print(quantile(min_ESS_vec))
#         ##
#         print(paste("max_Rhat:"))
#         print(quantile(max_rhat_vec))
#         ##
#         # Return results compatible with original function
#         return(list(
#           elpd_kfold = elpd_kfold,
#           se_elpd = se_elpd,
#           K = K,
#           ##
#           min_ESS_vec = min_ESS_vec,
#           max_rhat_vec = max_rhat_vec,
#           ##
#           fold_results = successful_folds,
#           fold_assignments = fold_ids,
#           n_workers = n_workers,
#           timing = list(
#             total_cpu_mins = total_time,
#             wall_clock_mins = parallel_time,
#             speedup = total_time / parallel_time
#           ),
#           errors = Filter(function(x) !is.null(x) && x$status == "error", 
#                           fold_results_list)
#         ))
#         
# }








#' filter_kfold_by_ess
#' @export
filter_kfold_by_ess <- function(kfold_results,
                                min_ess_threshold = 400, 
                                verbose = TRUE) {
  
        # Extract ESS values and identify good folds
        ess_values <- kfold_results$min_ESS_vec
        good_folds <- which(ess_values >= min_ess_threshold)
        bad_folds <- which(ess_values < min_ess_threshold)
        
        if (verbose) {
          cat("Original K-fold results:\n")
          cat(sprintf("  Total folds: %d\n", kfold_results$K))
          cat(sprintf("  ELPD: %.1f (SE: %.1f)\n", kfold_results$elpd_kfold, kfold_results$se_elpd))
          
          if (length(bad_folds) > 0) {
            cat("\nExcluding folds with ESS <", min_ess_threshold, ":\n")
            for (i in bad_folds) {
              cat(sprintf("  Fold %d: ESS = %.1f\n", i, ess_values[i]))
            }
          }
        }
        
        if (!length(good_folds)) {
          stop("No folds meet the requested ESS threshold")
        }

        # Filter fold results
        filtered_fold_results <- kfold_results$fold_results[good_folds]
        pointwise <- kfold_results$pointwise_elpd
        if (is.null(pointwise)) {
          pointwise <- kfold_pointwise_from_results(
            kfold_results$fold_results,
            length(kfold_results$fold_assignments)
          )
        }
        studies_used <- which(kfold_results$fold_assignments %in% good_folds)
        pointwise <- pointwise[studies_used]
        summary <- kfold_elpd_summary(pointwise)
        
        # Create filtered results object
        filtered_results <- list(
          elpd_kfold = summary$elpd_kfold,
          se_elpd = summary$se_elpd,
          pointwise_elpd = pointwise,
          K_original = kfold_results$K,
          K_used = length(good_folds),
          min_ESS_vec = ess_values[good_folds],
          max_rhat_vec = if (!is.null(kfold_results$max_rhat_vec)) {
            kfold_results$max_rhat_vec[good_folds]
          } else {
            NULL
          },
          fold_results = filtered_fold_results,
          excluded_folds = bad_folds,
          min_ess_threshold = min_ess_threshold
        )
        
        if (verbose) {
          cat(sprintf("\nFiltered results (using %d/%d folds):\n", 
                      length(good_folds), kfold_results$K))
          cat(sprintf("  ELPD: %.1f (SE: %.1f)\n", summary$elpd_kfold, summary$se_elpd))
          cat(sprintf("  Change in ELPD: %.1f\n", summary$elpd_kfold - kfold_results$elpd_kfold))
        }
        
        return(filtered_results)
  
}



 


 
#' compare_kfold_thresholds
#' @export
compare_kfold_thresholds <- function(kfold_results, 
                                     thresholds = c(1, 50, 100, 200, 400, 600, 650, 700, 704, 706, 725, 750, 800, 1000)) {
  
        cat("ESS threshold comparison:\n")
        cat("Threshold | Folds used | ELPD (SE)\n")
        cat("----------|------------|-------------\n")
        
        for (thresh in thresholds) {
          filtered <- filter_kfold_by_ess(kfold_results, thresh, verbose = FALSE)
          cat(sprintf("%9d | %10d | %.0f (%.0f)\n", 
                      thresh, 
                      filtered$K_used,
                      filtered$elpd_kfold, 
                      filtered$se_elpd))
        }
  
}







#' compare_models_k_fold_CV
#' @export
compare_models_k_fold_CV <- function(  kfold_results_list, 
                                       model_names = NULL,
                                       min_ess_threshold = 100,
                                       exclude_low_ess = FALSE,
                                       allow_incomplete = FALSE,
                                       reference_model = "best",
                                       verbose = TRUE) {

        n_models <- length(kfold_results_list)
        if (n_models < 2L) {
          stop("At least two K-fold result objects are required")
        }

        results <- lapply(kfold_results_list, function(x) {
          if (!is.null(x$outs_kfold)) x$outs_kfold else x
        })

        if (is.null(model_names)) {
          model_names <- paste0("Model_", seq_len(n_models))
        }
        if (length(model_names) != n_models) {
          stop("model_names must have one entry per model")
        }

        if (verbose) cat("Checking fold consistency across models...\n")
        ref_fold_assignments <- results[[1L]]$fold_assignments
        ref_K <- results[[1L]]$K
        n_studies <- length(ref_fold_assignments)

        validate_kfold_assignments(ref_fold_assignments, ref_K, n_studies)
        for (i in seq_len(n_models)[-1L]) {
          if (!identical(results[[i]]$fold_assignments, ref_fold_assignments)) {
            stop(sprintf(
              "Model %d (%s) has different fold assignments from Model 1",
              i, model_names[i]
            ))
          }
          if (!identical(results[[i]]$K, ref_K)) {
            stop(sprintf("Model %d (%s) has a different K", i, model_names[i]))
          }
        }

        for (fold in seq_len(ref_K)) {
          ref_holdout <- results[[1L]]$fold_results[[fold]]$holdout_studies
          for (i in seq_len(n_models)[-1L]) {
            model_holdout <- results[[i]]$fold_results[[fold]]$holdout_studies
            if (!identical(sort(ref_holdout), sort(model_holdout))) {
              stop(sprintf("Fold %d has different holdout studies between models", fold))
            }
          }
        }
        if (verbose) cat("All models used identical fold assignments\n\n")

        all_min_ess <- matrix(NA_real_, nrow = ref_K, ncol = n_models)
        colnames(all_min_ess) <- model_names
        for (i in seq_len(n_models)) {
          if (!is.null(results[[i]]$min_ESS_vec)) {
            all_min_ess[, i] <- results[[i]]$min_ESS_vec
          }
        }
        min_ess_per_fold <- apply(all_min_ess, 1L, function(x) {
          if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
        })
        diagnostic_folds <- which(!is.na(min_ess_per_fold) &
                                    min_ess_per_fold < min_ess_threshold)

        folds_used <- if (isTRUE(exclude_low_ess)) {
          setdiff(seq_len(ref_K), diagnostic_folds)
        } else {
          seq_len(ref_K)
        }
        if (!length(folds_used)) {
          stop("No folds remain after applying the ESS threshold")
        }

        if (length(diagnostic_folds)) {
          warning(
            sprintf(
              "%d fold(s) have ESS below %s; rerun those folds before final inference%s",
              length(diagnostic_folds),
              min_ess_threshold,
              if (exclude_low_ess) " (they were excluded by explicit request)" else ""
            ),
            call. = FALSE
          )
        }

        studies_used <- which(ref_fold_assignments %in% folds_used)
        filtered_results <- vector("list", n_models)
        names(filtered_results) <- model_names
        pointwise_by_model <- vector("list", n_models)

        for (i in seq_len(n_models)) {
          pointwise <- results[[i]]$pointwise_elpd
          if (is.null(pointwise)) {
            pointwise <- kfold_pointwise_from_results(
              results[[i]]$fold_results,
              n_studies
            )
          }
          if (!allow_incomplete && any(!is.finite(pointwise[studies_used]))) {
            stop(sprintf("Model %s has incomplete held-out predictions", model_names[i]))
          }
          selected <- pointwise[studies_used]
          summary <- kfold_elpd_summary(selected)
          pointwise_by_model[[i]] <- selected
          fold_elpds <- vapply(folds_used, function(fold) {
            sum(pointwise[ref_fold_assignments == fold], na.rm = TRUE)
          }, numeric(1L))

          filtered_results[[i]] <- list(
            elpd_kfold = summary$elpd_kfold,
            se_elpd = summary$se_elpd,
            pointwise_elpd = selected,
            K_original = ref_K,
            K_used = length(folds_used),
            fold_elpds = fold_elpds,
            excluded_folds = setdiff(seq_len(ref_K), folds_used),
            diagnostic_folds = diagnostic_folds,
            min_ESS_vec = all_min_ess[folds_used, i]
          )
        }

        comparison_table <- data.frame(
          Model = model_names,
          ELPD = vapply(filtered_results, `[[`, numeric(1L), "elpd_kfold"),
          SE = vapply(filtered_results, `[[`, numeric(1L), "se_elpd"),
          Delta_ELPD = 0,
          Delta_SE = 0,
          stringsAsFactors = FALSE
        )

        reference_index <- if (identical(reference_model, "best")) {
          which.max(comparison_table$ELPD)
        } else if (is.character(reference_model) && length(reference_model) == 1L) {
          match(reference_model, model_names)
        } else if (is.numeric(reference_model) && length(reference_model) == 1L) {
          as.integer(reference_model)
        } else {
          NA_integer_
        }
        if (is.na(reference_index) || reference_index < 1L || reference_index > n_models) {
          stop("reference_model must be 'best', a model name, or a valid model index")
        }

        reference_pointwise <- pointwise_by_model[[reference_index]]
        for (i in seq_len(n_models)) {
          delta <- pointwise_by_model[[i]] - reference_pointwise
          delta <- delta[is.finite(delta)]
          comparison_table$Delta_ELPD[i] <- sum(delta)
          comparison_table$Delta_SE[i] <- if (length(delta) > 1L) {
            sqrt(length(delta) * stats::var(delta))
          } else {
            NA_real_
          }
        }

        if (verbose) {
          cat(sprintf("Final comparison using %d/%d folds:\n", length(folds_used), ref_K))
          print(comparison_table, digits = 1)
        }

        list(
          filtered_results = filtered_results,
          comparison_table = comparison_table,
          excluded_folds = setdiff(seq_len(ref_K), folds_used),
          diagnostic_folds = diagnostic_folds,
          K_used = length(folds_used),
          min_ess_matrix = all_min_ess,
          fold_assignments = ref_fold_assignments,
          reference_model = model_names[reference_index]
        )
  
}




#' summarise_model_comparison_k_fold_CV
#' @export
summarise_model_comparison_k_fold_CV <- function(comparison_results) {
  
        cat("\nModel ranking by ELPD:\n")
        table <- comparison_results$comparison_table
        table <- table[order(table$ELPD, decreasing = TRUE), ]
        
        for (i in 1:nrow(table)) {
          if (i == 1) {
            cat(sprintf("%d. %s: ELPD = %.0f (%.0f) [BEST]\n", 
                        i, table$Model[i], table$ELPD[i], table$SE[i]))
          } else {
            cat(sprintf("%d. %s: ELPD = %.0f (%.0f), Δ = %.0f (%.0f)\n",
                        i, table$Model[i], table$ELPD[i], table$SE[i],
                        table$Delta_ELPD[i], table$Delta_SE[i]))
          }
        }
        cat("\nModels are ranked by ELPD (higher is better). Delta SE describes\n")
        cat("uncertainty in the paired difference; it is not a selection threshold.\n")
  
}




#' Generate LaTeX table for k-fold CV model comparison
#' @export
kfold_CV_to_latex <- function(comparison_results,
                              target_condition = "GAD",  # "GAD" or "AAD"
                              K,
                              label_suffix = "Meta_reg_kfold_cv") {
  
  # Get comparison table and sort by ELPD (best first)
  table <- comparison_results$comparison_table
  table <- table[order(table$ELPD, decreasing = TRUE), ]
  
  # Function to format numbers with commas
  fmt_num <- function(x, digits = 0) {
    formatC(round(x, digits), format = "f", digits = digits, big.mark = ",")
  }
  
  # Function to clean model names and extract covariates
  parse_model <- function(model_name) {
    # Extract model number
    if (grepl("Model_0|intercept-only", model_name)) {
      return(list(num = 0, covs = "Intercept-only"))
    } else if (grepl("FULL|Model_5", model_name)) {
      return(list(num = 5, covs = paste0("+ ", target_condition, " prev. + ref. test + study setting + RoB")))
    }
    
    # Extract number
    num <- as.numeric(gsub(".*Model_?(\\d+).*", "\\1", model_name))
    
    # Parse covariates from model name
    covs <- ""
    if (grepl("prev", model_name, ignore.case = TRUE)) {
      covs <- paste0(covs, "+ ", target_condition, " prev. ")
    }
    if (grepl("ref_test|ref test", model_name, ignore.case = TRUE)) {
      covs <- paste0(covs, "+ ref. test ")
    }
    if (grepl("study_setting|study setting", model_name, ignore.case = TRUE)) {
      covs <- paste0(covs, "+ study setting ")
    }
    if (grepl("RoB", model_name)) {
      covs <- paste0(covs, "+ RoB ")
    }
    
    covs <- trimws(covs)
    if (covs == "") covs <- "Intercept-only"
    
    return(list(num = num, covs = covs))
  }
  
  # Start building LaTeX
  cat("%%%%\n%%%%\n")
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\caption{\n")
  cat(sprintf("K-fold cross-validation results comparing meta-regression models ($K = %d$);\n", K))
  cat("ranked in order from best to worst.\n")
  cat("}\n")
  cat("\\small\n")
  cat("\\begin{tabular}{llcccc}\n")
  cat("%%%%\n")
  cat("\\toprule\n")
  cat("%%%%\n")
  cat("\\textbf{Model$^*$} & \\textbf{Covariates} & \\textbf{ELPD} & \\textbf{SE} & \\textbf{$\\Delta$ELPD} & \\textbf{$\\Delta$SE} \\\\\n")
  cat("%%%%\n")
  cat("\\midrule\n")
  cat("%%%%\n")
  
  # Generate rows
  for (i in 1:nrow(table)) {
    parsed <- parse_model(table$Model[i])
    model_num <- parsed$num
    covs <- parsed$covs
    
    elpd <- table$ELPD[i]
    se <- table$SE[i]
    
    is_best <- (i == 1)
    is_intercept <- (model_num == 0)
    is_reference <- identical(table$Model[i], comparison_results$reference_model)
    
    # Format model name
    if (is_best) {
      model_str <- sprintf("\\textbf{MR model %d}", model_num)
      cov_str <- sprintf("\\textbf{%s}", covs)
      elpd_str <- sprintf("$\\mathbf{%s}$", fmt_num(elpd))
      se_str <- sprintf("$\\mathbf{%s}$", fmt_num(se))
      
      if (is_reference) {
        delta_elpd_str <- "\\textbf{--}"
        delta_se_str <- "\\textbf{--}"
      } else {
        delta_elpd_str <- sprintf("$\\mathbf{%s}$", fmt_num(table$Delta_ELPD[i]))
        delta_se_str <- sprintf("$\\mathbf{%s}$", fmt_num(table$Delta_SE[i]))
      }
    } else {
      model_str <- sprintf("MR model %d", model_num)
      cov_str <- covs
      elpd_str <- sprintf("$%s$", fmt_num(elpd))
      se_str <- sprintf("$%s$", fmt_num(se))
      
      if (is_reference) {
        delta_elpd_str <- "--"
        delta_se_str <- "--"
      } else {
        delta_elpd_str <- sprintf("$%s$", fmt_num(table$Delta_ELPD[i]))
        delta_se_str <- sprintf("$%s$", fmt_num(table$Delta_SE[i]))
      }
    }
    
    # Add %% before intercept-only row
    if (is_intercept) {
      cat("%%\n")
    }
    
    # Print row
    cat(sprintf("%s & %s & %s & %s & %s & %s \\\\\n",
                model_str, cov_str, elpd_str, se_str, delta_elpd_str, delta_se_str))
    
    # Add %% after intercept-only row
    if (is_intercept) {
      cat("%%\n")
    }
  }
  
  cat("%%%%\n")
  cat("\\bottomrule\n")
  cat("%%%%\n")
  cat("\\end{tabular}\n")
  cat("\\begin{tablenotes}\n")
  cat("\\footnotesize\n")
  cat("\\item $^*$ All models use fixed cutpoints with compound symmetry correlation structure (i.e. based on Model A from baseline analysis).\n")
  cat("\\item ELPD = expected log predictive density; SE = standard error\n")
  cat(sprintf(
    "\\item $\\Delta$ELPD and $\\Delta$SE use paired study-wise differences relative to %s.\n",
    gsub("_", "\\\\_", comparison_results$reference_model)
  ))
  cat("\\end{tablenotes}\n")
  cat(sprintf("\\label{table_application_%s_%s}\n", target_condition, label_suffix))
  cat("\\end{table}\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
}






#' Create LaTeX table for k-fold CV comparison
#' @export
make_latex_table_k_fold_CV <- function(comparison_results,
                                       caption = "K-fold cross-validation model comparison",
                                       label = "tab:kfold_comparison",
                                       digits_elpd = 0,
                                       digits_se = 0,
                                       show_rank = TRUE,
                                       escape_underscore = TRUE) {
          
          # Get the comparison table
          table_data <- comparison_results$comparison_table
          n_models <- nrow(table_data)
          
          # Sort by ELPD (best first)
          table_data <- table_data[order(table_data$ELPD, decreasing = TRUE), ]
          
          # Add rank column if requested
          if (show_rank) {
            table_data$Rank <- 1:n_models
          }
          
          # Format the values
          formatted_table <- data.frame(
            Rank = if(show_rank) table_data$Rank else NULL,
            Model = table_data$Model
          )
          
          # Escape underscores in model names if needed
          if (escape_underscore) {
            formatted_table$Model <- gsub("_", "\\\\_", formatted_table$Model)
          }
          
          # Format ELPD with SE
          formatted_table$ELPD <- sprintf("$%.*f$ ($%.*f$)", 
                                          digits_elpd, table_data$ELPD,
                                          digits_se, table_data$SE)
          
          # Format paired Delta ELPD with its SE. The SE quantifies uncertainty
          # in the estimated difference; it is not a hypothesis-test gate for
          # ranking models by predictive performance.
          formatted_table$`$\\Delta$ELPD` <- character(n_models)
          formatted_table$`$\\Delta$ELPD`[1] <- "---" # Best model
          
          for (i in 2:n_models) {
            formatted_table$`$\\Delta$ELPD`[i] <- sprintf("$%.*f$ ($%.*f$)",
                                                          digits_elpd, table_data$Delta_ELPD[i],
                                                          digits_se, table_data$Delta_SE[i])
          }
          
          # Start building LaTeX code
          latex_code <- character()
          
          # Table header
          latex_code <- c(latex_code, "\\begin{table}[H]")
          latex_code <- c(latex_code, "\\centering")
          latex_code <- c(latex_code, sprintf("\\caption{%s}", caption))
          latex_code <- c(latex_code, sprintf("\\label{%s}", label))
          
          # Determine column alignment
          if (show_rank) {
            col_align <- "clcc"
          } else {
            col_align <- "lcc"
          }
          
          latex_code <- c(latex_code, sprintf("\\begin{tabular}{%s}", col_align))
          latex_code <- c(latex_code, "\\hline")
          
          # Column headers
          if (show_rank) {
            headers <- "Rank & Model & ELPD (SE) & $\\Delta$ELPD (SE) \\\\"
          } else {
            headers <- "Model & ELPD (SE) & $\\Delta$ELPD (SE) \\\\"
          }
          latex_code <- c(latex_code, headers)
          latex_code <- c(latex_code, "\\hline")
          
          # Table rows
          for (i in 1:n_models) {
            if (show_rank) {
              row <- sprintf("%d & %s & %s & %s \\\\",
                             formatted_table$Rank[i],
                             formatted_table$Model[i],
                             formatted_table$ELPD[i],
                             formatted_table$`$\\Delta$ELPD`[i])
            } else {
              row <- sprintf("%s & %s & %s \\\\",
                             formatted_table$Model[i],
                             formatted_table$ELPD[i],
                             formatted_table$`$\\Delta$ELPD`[i])
            }
            
            # Add best model indicator
            if (i == 1) {
              row <- gsub("\\\\\\\\", " \\\\textbf{[Best]} \\\\\\\\", row)
            }
            
            latex_code <- c(latex_code, row)
          }
          
          latex_code <- c(latex_code, "\\hline")
          latex_code <- c(latex_code, "\\end{tabular}")
          
          # Add footnotes
          footnotes <- character()
          
          # Note about excluded folds
          if (length(comparison_results$excluded_folds) > 0) {
            excluded_str <- paste(comparison_results$excluded_folds, collapse = ", ")
            footnotes <- c(footnotes, 
                           sprintf("\\footnotesize Note: %d of %d folds used (excluded fold%s %s due to low ESS).",
                                   comparison_results$K_used,
                                   comparison_results$filtered_results[[1]]$K_original,
                                   ifelse(length(comparison_results$excluded_folds) > 1, "s", ""),
                                   excluded_str))
          }
          
          # Note about ELPD
          footnotes <- c(footnotes,
                         "\\footnotesize ELPD = expected log predictive density (higher is better). Models are ranked by ELPD; the paired $\\Delta$SE describes uncertainty and is not a hypothesis-test threshold.")
          
          # Add footnotes with line breaks
          if (length(footnotes) > 0) {
            latex_code <- c(latex_code, "\\vspace{0.2cm}")
            latex_code <- c(latex_code, "\\begin{flushleft}")
            for (fn in footnotes) {
              latex_code <- c(latex_code, fn, "\\\\")
            }
            latex_code <- c(latex_code, "\\end{flushleft}")
          }
          
          latex_code <- c(latex_code, "\\end{table}")
          
          # Combine all lines
          latex_output <- paste(latex_code, collapse = "\n")
          
          # Print to console
          cat(latex_output)
          cat("\n")
          
          # Return invisibly
          invisible(latex_output)
  
}
