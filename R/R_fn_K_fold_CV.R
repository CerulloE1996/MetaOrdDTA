
#' resize_init_list
#' @keywords internal
#' @export
resize_init_list <- function(init_lists_per_chain, 
                             n_chains_new) {
        
        n_chains_old <- length(init_lists_per_chain)
        
        if (n_chains_new == n_chains_old) {
          # No change needed
          return(init_lists_per_chain)
        } else if (n_chains_new < n_chains_old) {
          # Fewer chains - just take first n_chains_new
          return(init_lists_per_chain[1:n_chains_new])
        } else {
          # More chains - repeat cyclically
          indices <- rep(1:n_chains_old, length.out = n_chains_new)
          return(init_lists_per_chain[indices])
        }
        
}



#' Stratified sampling by test availability
#' create_folds
#' @keywords internal
#' @export
create_folds <- function(study_test_matrix, 
                         K, 
                         seed) {
        
        set.seed(seed)
  
        # Group studies by test pattern
        test_patterns <- apply(study_test_matrix, 1, paste, collapse = "")
        
        # Create folds within each pattern
        folds <- rep(NA, nrow(study_test_matrix))
        for (pattern in unique(test_patterns)) {
          idx <- which(test_patterns == pattern)
          folds[idx] <- sample(rep(1:K, length.out = length(idx)))
        }
        
        return(folds)
  
}



#' R_fn_run_kfold_cv
#' @keywords internal
#' @export
R_fn_run_kfold_cv <- function(   debugging,
                                 ##
                                 K = 10, 
                                 ##
                                 model_prep_obj,
                                 ##
                                 stan_data_list, 
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
                                 use_BayesMVP_for_faster_summaries
) {
  
          cat("Using", K, "-fold CV (LOO inappropriate for NMA models)\n")
          
          # Try to keep network connected in each fold
          n_studies <- stan_data_list$n_studies
          ##
          fold_size <- ceiling(n_studies / K)
          fold_ids <- rep(1:K, each = fold_size)[1:n_studies]
          fold_ids <- sample(fold_ids)  # Randomize
          
          elpd_kfold <- 0
          fold_results <- list()
          
          k <- 1
          
          for (k in 1:K) {
                  
                  cat("\nFold", k, "- holding out", sum(fold_ids == k), "studies\n")
                  
                  # Check network connectivity
                  holdout_studies   <- which(fold_ids == k)
                  remaining_studies <- which(fold_ids != k)
                  ##
                  priors$K_fold_CV_indicator <- as.integer(fold_ids != k)
                  ##
                  model_samples_obj <-  model_prep_obj$sample(   
                    n_burnin = n_burnin,
                    n_iter   = n_iter,
                    adapt_delta = adapt_delta, 
                    max_treedepth = max_treedepth,
                    metric_shape = "diag_e",
                    ##
                    priors = priors,
                    ##
                    n_chains = n_chains,
                    ##
                    init_lists_per_chain = init_lists_per_chain
                  )
                  ##
                  model_summary_and_trace_obj <- model_samples_obj$summary(
                    compute_main_params = TRUE,
                    compute_transformed_parameters = FALSE, 
                    compute_generated_quantities = TRUE,
                    ##
                    save_log_lik_trace = TRUE,
                    ##
                    use_BayesMVP_for_faster_summaries = TRUE)
                  ##
                  rm(model_samples_obj)
                  gc() ; gc()
                  ##
                  # posterior_draws <- model_summary_and_trace_obj$get_all_traces_list()
                  trace_gq <- model_summary_and_trace_obj$get_all_traces_list()$traces_as_arrays$trace_gq
                  ##
                  rm(model_summary_and_trace_obj)
                  gc() ; gc()
                  ##
                  ## Extract log_lik_study for ALL studies
                  log_lik_study_array <- filter_param_draws_string(
                    named_draws_array = trace_gq,
                    param_string = "log_lik_study",
                    condition = "exact_match",
                    exclude_NA = TRUE
                  )
                  ##
                  rm(trace_gq)
                  gc() ; gc()
                  ##
                  # Calculate holdout log-likelihood (only for holdout studies)
                  holdout_indices <- which(fold_ids == k)
                  log_lik_holdout_samples <- log_lik_study_array[, , holdout_indices, drop = FALSE]
                  sum(log_lik_study_array)
                  sum(log_lik_holdout_samples)
                  ##
                  # Sum across holdout studies and average across iterations/chains
                  elpd_holdout <- mean(apply(log_lik_holdout_samples, c(1, 2), sum))
                  ##
                  fold_results[[k]] <- list(
                    elpd = elpd_holdout,
                    n_holdout = length(holdout_indices),
                    holdout_studies = holdout_indices
                  )
                  ##
                  elpd_kfold <- elpd_kfold + elpd_holdout
                  ##
                  fold_results[[k]]
                  elpd_kfold
            
          }
          
          # Compute SE using fold-level results
          fold_elpds <- sapply(fold_results, function(x) x$elpd)
          se_elpd <- sd(fold_elpds) * sqrt(K)
          
          return(list(
            elpd_kfold = elpd_kfold,
            se_elpd = se_elpd,
            K = K,
            fold_results = fold_results,
            fold_assignments = fold_ids
          ))
  
}










#' R_fn_run_kfold_cv_parallel
#' Run k-fold cross-validation in parallel across multiple folds
#' @export
R_fn_run_kfold_cv_parallel <- function(debugging,
                                       ##
                                       K = 10, 
                                       seed = 123,
                                       fold_assignments,
                                       ##
                                       model_prep_obj,
                                       stan_data_list, 
                                       ##
                                       cmdstanr_args,
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
                                       output_dir = "cv_results",
                                       ##
                                       use_BayesMVP_for_faster_summaries
) {
  
        require(doParallel)
        require(foreach)
        ##
        if (use_BayesMVP_for_faster_summaries) { 
          require(RcppParallel)
          require(BayesMVP)
        }
        ##
        cat("Using", K, "-fold CV (LOO inappropriate for NMA models) - PARALLEL VERSION\n")
        
        # Calculate optimal number of workers
        if (is.null(n_workers)) {
          # Use all available cores divided by chains per model
          n_workers <- min(K, floor(parallel::detectCores() / n_chains))
        }
        cat(sprintf("Running %d folds in parallel using %d workers (%d chains each = %d cores)\n", 
                    K, n_workers, n_chains, n_workers * n_chains))
        
        # Create output directory if it doesn't exist
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
        minimal_env$minimal_env$prior_only <- model_prep_obj$basic_model_options$prior_only
        ##
        minimal_env$cov_data <- model_prep_obj$cov_data
        minimal_env$x <- model_prep_obj$x
        minimal_env$indicator_index_test_in_study <- model_prep_obj$indicator_index_test_in_study
        # ##
        # With this:
        if (is.null(fold_assignments)) {
          stop("fold_assignments must be provided! Create them externally with create_folds()")
        }
        fold_ids <- fold_assignments
        ##
        # Validate they match
        if (length(fold_ids) != model_prep_obj$indicator_index_test_in_study %>% nrow()) {
          stop("fold_assignments length doesn't match number of studies!")
        }
        if (max(fold_ids) != K) {
          stop("fold_assignments has different K than specified!")
        }
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
      
        ##
       
        # minimal_env$cmdstanr_args <- cmdstanr_args
        ##
        # Set up parallel cluster
        cl <- makeCluster(n_workers, 
                          outfile = "")
        registerDoParallel(cl)
        ##
        # # Export everything to all workers
        # clusterExport(cl,
        #               varlist = ls(),
        #               envir = environment())
        # Export necessary objects to workers
        minimal_env$resize_init_list <- MetaOrdDTA::resize_init_list
        minimal_env$filter_param_draws_string <- MetaOrdDTA::filter_param_draws_string
        ##
        clusterExport(cl, 
                      ls(minimal_env), 
                      envir = minimal_env)
        clusterEvalQ(cl, {
          library(cmdstanr)
          library(MetaOrdDTA)
          library(posterior)
          library(dplyr)
          library(stringr)
        })
        ##
        k <- 2
        ##
        # Run folds in parallel
        fold_results_list <- foreach(k = 1:K,
                                     .combine = 'c',
                                     # .export = ls(envir = asNamespace("MetaOrdDTA")),
                                     .packages = c("MetaOrdDTA",
                                                   "cmdstanr", 
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
                                               require(cmdstanr)
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
                                               cmdstanr_args <- list(    seed = seed + k,
                                                                         output_dir = fold_output_dir,  # CSV files go here
                                                                         output_basename = paste0("fold_", k, "_chain"),  # Unique prefix
                                                                         save_latent_dynamics = FALSE  # Don't save extra diagnostic files
                                               )
                                               ##
                                               # Add this RIGHT AFTER creating fold_output_dir (around line 136)
                                               # Create unique temp directory
                                               fold_temp_dir <- file.path(tempdir(), fold_id)
                                               dir.create(fold_temp_dir, showWarnings = FALSE, recursive = TRUE)
                                               
                                               # Save current environment
                                               old_tmpdir <- Sys.getenv("TMPDIR")
                                               old_temp <- Sys.getenv("TEMP")
                                               old_tmp <- Sys.getenv("TMP")
                                               
                                               # Set new temp directories
                                               Sys.setenv(TMPDIR = fold_temp_dir)
                                               Sys.setenv(TEMP = fold_temp_dir)
                                               Sys.setenv(TMP = fold_temp_dir)
                                               
                                               # Also set Stan-specific temp
                                               Sys.setenv(CMDSTAN_TEMP_DIR = fold_temp_dir)
                                               ##
                                               init_lists_per_chain <- model_prep_obj$init_lists_per_chain
                                               # priors <- model_prep_obj$priors
                                               stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
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
                                               init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain, 
                                                                                         n_chains_new = n_chains)
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
                                               cmdstanr_args <- list(data = c(stan_data_list, priors),
                                                                     seed = seed + k,
                                                                     output_dir = fold_output_dir,  # CSV files go here
                                                                     output_basename = paste0("fold_", k, "_chain"),  # Unique prefix
                                                                     save_latent_dynamics = FALSE  # Don't save extra diagnostic files
                                                                     )
                                               final_stan_args <- modifyList(final_stan_args, cmdstanr_args)
                                               # ##
                                               # final_stan_args <- modifyList(final_stan_args, cmdstanr_args)
                                               # ##
                                               stan_mod_samples <- do.call(model_prep_obj$internal_obj$outs_stan_compile$stan_model_obj$sample, 
                                                                           final_stan_args)
                                               ##
                                               log_lik_study_array <- stan_mod_samples$draws(variables = c("log_lik_study"))
                                               str(log_lik_study_array)
                                               ##
                                               log_lik_study_tibble <- stan_mod_samples$summary(variables = c("log_lik_study"))
                                               min_ESS <- min(log_lik_study_tibble$ess_bulk, na.rm = TRUE)
                                               max_rhat <- max(log_lik_study_tibble$rhat, na.rm = TRUE)
                                               ##
                                               ## Calculate holdout log-likelihood
                                               log_lik_holdout_samples <- log_lik_study_array[, , holdout_studies, drop = FALSE]
                                               ##
                                               ## Sum across holdout studies and average across iterations/chains
                                               elpd_holdout <- mean(apply(log_lik_holdout_samples, c(1, 2), sum))
                                               ##
                                               ### Update fold result
                                               fold_result$status <- "completed"
                                               fold_result$elpd <- elpd_holdout
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
                                               ##
                                               ## Optionally save fold-specific results
                                               saveRDS(fold_result, 
                                                       file = file.path(output_dir, sprintf("fold_%02d_results.rds", k)))
                                               
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
  
  
        # Stop cluster
        stopCluster(cl)
        
        # Clean up temporary files
        for (k in 1:K) {
          fold_file <- file.path(output_dir, paste0("cv_environment_fold_", k, ".RData"))
          if (file.exists(fold_file)) {
            file.remove(fold_file)
          }
        }
             
        # Process results
        successful_folds <- Filter(function(x) !is.null(x) && x$status == "completed", 
                                   fold_results_list)
        
        if (length(successful_folds) < K) {
          warning(sprintf("Only %d/%d folds completed successfully", 
                          length(successful_folds), K))
        }
        
        # Calculate total elpd_kfold
        elpd_kfold <- sum(sapply(successful_folds, function(x) x$elpd))
        
        # Compute SE using fold-level results
        fold_elpds <- sapply(successful_folds, function(x) x$elpd)
        se_elpd <- sd(fold_elpds) * sqrt(K)
        
        # Summary of timing
        total_time <- sum(sapply(successful_folds, function(x) x$duration_mins))
        parallel_time <- max(sapply(successful_folds, function(x) x$duration_mins))
        ##
        min_ESS_vec <- max_rhat_vec <- c()
        for (k in 1:length(successful_folds)) {
          min_ESS_vec[k]  <- successful_folds[[k]]$min_ESS
          max_rhat_vec[k] <- successful_folds[[k]]$max_rhat
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
          K = K,
          ##
          min_ESS_vec = min_ESS_vec,
          max_rhat_vec = max_rhat_vec,
          ##
          fold_results = successful_folds,
          fold_assignments = fold_ids,
          n_workers = n_workers,
          timing = list(
            total_cpu_mins = total_time,
            wall_clock_mins = parallel_time,
            speedup = total_time / parallel_time
          ),
          errors = Filter(function(x) !is.null(x) && x$status == "error", 
                          fold_results_list)
        ))
  
}













#' R_fn_run_kfold_cv_parallel_BayesMVP
#' Run k-fold cross-validation in parallel across multiple folds
#' @export
R_fn_run_kfold_cv_parallel_BayesMVP <- function( debugging,
                                                 ##
                                                 K = 10, 
                                                 seed = 123,
                                                 ##
                                                 model_prep_obj,
                                                 Stan_data_list, 
                                                 ##
                                                 cmdstanr_args,
                                                 ##
                                                 priors,
                                                 init_lists_per_chain,
                                                 ##
                                                 n_burnin,
                                                 n_iter,
                                                 ##
                                                 n_chains_burnin,
                                                 n_chains_sampling,
                                                 ##
                                                 adapt_delta,
                                                 max_treedepth,
                                                 ##
                                                 n_workers = NULL,
                                                 output_dir = "cv_results",
                                                 ##
                                                 use_BayesMVP_for_faster_summaries,
                                                 stanc_args
) {
  
  require(doParallel)
  require(foreach)
  ##
  if (use_BayesMVP_for_faster_summaries) { 
    require(RcppParallel)
    require(BayesMVP)
  }
  ##
  cat("Using", K, "-fold CV (LOO inappropriate for NMA models) - PARALLEL VERSION\n")
  
  # Calculate optimal number of workers
  if (is.null(n_workers)) {
    # Use all available cores divided by chains per model
    n_workers <- min(K, floor(parallel::detectCores() / n_chains_sampling))
  }
  cat(sprintf("Running %d folds in parallel using %d workers (%d chains each = %d cores)\n", 
              K, n_workers, n_chains_sampling, n_workers * n_chains_sampling))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Create a minimal environment
  minimal_env <- new.env()
  ##
  # Set up fold assignments first (before parallel)
  minimal_env$n_studies <- Stan_data_list$n_studies
  minimal_env$n_index_tests <- Stan_data_list$n_index_tests
  minimal_env$n_thr_max <- Stan_data_list$n_thr_max
  ##
  minimal_env$indicator_index_test_in_study <- Stan_data_list$indicator_index_test_in_study
  ##
  fold_ids <-  create_folds(K = K, 
                            study_test_matrix = Stan_data_list$indicator_index_test_in_study,
                            seed = seed)
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
  ##
  minimal_env$Stan_data_list <- Stan_data_list
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
  ##
  minimal_env$n_chains_sampling <- n_chains_sampling
  minimal_env$n_chains_burnin <- n_chains_burnin
  ##
  minimal_env$stanc_args <- stanc_args
  ##
  # Set up parallel cluster
  cl <- makeCluster(n_workers, 
                    outfile = "")
  registerDoParallel(cl)
  ##
  minimal_env$resize_init_list <- MetaOrdDTA::resize_init_list
  minimal_env$filter_param_draws_string <- MetaOrdDTA::filter_param_draws_string
  ##
  clusterExport(cl, 
                ls(minimal_env), 
                envir = minimal_env)
  clusterEvalQ(cl, {
    library(cmdstanr)
    library(MetaOrdDTA)
    library(posterior)
    library(dplyr)
    library(stringr)
  })
  ##
  k <- 2
  ##
  # Run folds in parallel
  fold_results_list <- foreach(k = 1:K,
                               .combine = 'c',
                               # .export = ls(envir = asNamespace("MetaOrdDTA")),
                               .packages = c("MetaOrdDTA",
                                             "cmdstanr", 
                                             "bridgestan", 
                                             "RcppParallel",
                                             "BayesMVPtest",
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
                                   require(cmdstanr)
                                   require(posterior)
                                   require(dplyr)
                                   ##
                                   if (use_BayesMVP_for_faster_summaries) { 
                                     require(RcppParallel)
                                     require(BayesMVPtest)
                                   }
                                   ## set bs environment variable (otherwise it'll try downloading it even if already installed...)
                                   bs_path <- bridgestan_path()
                                   Sys.setenv(BRIDGESTAN = bs_path)
                                   ##
                                   # Create unique directories for this fold
                                   fold_id <- sprintf("fold_%02d_pid_%s", k, Sys.getpid())
                                   fold_output_dir <- file.path(output_dir, fold_id)
                                   dir.create(fold_output_dir, showWarnings = FALSE, recursive = TRUE)
                                   ##
                                   # Add this RIGHT AFTER creating fold_output_dir (around line 136)
                                   # Create unique temp directory
                                   fold_temp_dir <- file.path(tempdir(), fold_id)
                                   dir.create(fold_temp_dir, showWarnings = FALSE, recursive = TRUE)
                                   ##
                                   # Save current environment
                                   old_tmpdir <- Sys.getenv("TMPDIR")
                                   old_temp <- Sys.getenv("TEMP")
                                   old_tmp <- Sys.getenv("TMP")
                                   ##
                                   # Set new temp directories
                                   Sys.setenv(TMPDIR = fold_temp_dir)
                                   Sys.setenv(TEMP = fold_temp_dir)
                                   Sys.setenv(TMP = fold_temp_dir)
                                   ##
                                   # Also set Stan-specific temp
                                   Sys.setenv(CMDSTAN_TEMP_DIR = fold_temp_dir)
                                   ##
                                   cat(sprintf("\n===== Worker starting fold %d/%d =====\n", k, K))
                                   ##
                                   init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain, 
                                                                             n_chains_new = n_chains_burnin)
                                   ##
                                   ## Create fold-specific priors:
                                   ##
                                   Stan_data_list$K_fold_CV_indicator <- as.integer(fold_ids != k)
                                   ##
                                   ## Get fold-specific information:
                                   ##
                                   holdout_studies <- which(fold_ids == k)
                                   remaining_studies <- which(fold_ids != k)
                                   n_holdout <- length(holdout_studies)
                                   ##
                                   cat(sprintf("Fold %d - holding out %d studies\n", k, n_holdout))
                                   ##
                                   partitioned_HMC <- TRUE ;    diffusion_HMC <- FALSE
                                   metric_shape_main <- "diag" ; metric_type_main <- "Empirical"
                                   ##
                                   model_samples <-  model_prep_obj$sample(  Stan_data_list = Stan_data_list,
                                                                             init_lists_per_chain = init_lists_per_chain,
                                                                             n_chains_burnin = n_chains_burnin,
                                                                             ##
                                                                             partitioned_HMC = partitioned_HMC,
                                                                             diffusion_HMC = diffusion_HMC,
                                                                             seed = seed,
                                                                             n_burnin = n_burnin,
                                                                             n_iter = n_iter,
                                                                             ##
                                                                             n_chains_sampling = n_chains_sampling,
                                                                             n_superchains = n_superchains,
                                                                             ##
                                                                             adapt_delta = 0.80,
                                                                             learning_rate = 0.075,
                                                                             ##
                                                                             metric_shape_main = metric_shape_main,
                                                                             metric_type_main = metric_type_main,
                                                                             ##
                                                                             tau_mult = 1.6,
                                                                             ##
                                                                             clip_iter = 50,
                                                                             clip_iter_tau = 150,
                                                                             ##
                                                                             interval_width_main = 10,
                                                                             interval_width_nuisance = 10,
                                                                             ratio_M_us = 0.50,
                                                                             ratio_M_main = 0.50,
                                                                             ##
                                                                             beta1_adam = 0.0,
                                                                             beta2_adam = 0.95,
                                                                             eps_adam = 1e-8,
                                                                             ##
                                                                             parallel_method = "RcppParallel",
                                                                             stanc_args = stanc_args)
                                   ##
                                   model_fit <- model_samples$summary(save_log_lik_trace = FALSE, 
                                                                      ##
                                                                      compute_nested_rhat = FALSE,
                                                                      ##
                                                                      compute_main_params = TRUE,
                                                                      compute_transformed_parameters = FALSE,
                                                                      compute_generated_quantities = TRUE)

                                   trace_gq <- model_fit_object$traces$traces_as_arrays$trace_generated_quantities
                                   ##
                                   log_lik_study_array <- extract_params_from_array_batch(debugging = FALSE,
                                                                                          array_3d = trace_gq, 
                                                                                          combine = TRUE,
                                                                                          param_strings_vec = c("log_lik_study"),
                                                                                          condition = "containing") 
                                   ##
                                   
                                   tibble_gq <- model_fit_object$summaries$summary_tibbles$summary_tibble_generated_qua
                                   ##
                                   log_lik_study_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
                                                                                 tibble = tibble_gq, 
                                                                                 param_strings_vec = c("log_lik_study"),
                                                                                 condition = "containing") %>% print(n=100)
                                   ##
                                   log_lik_study_tibble <- stan_mod_samples$summary(variables = c("log_lik_study"))
                                   min_ESS <- min(log_lik_study_tibble$n_eff, na.rm = TRUE)
                                   max_rhat <- max(log_lik_study_tibble$Rhat, na.rm = TRUE)
                                   ##
                                   ## Calculate holdout log-likelihood
                                   log_lik_holdout_samples <- log_lik_study_array[, , holdout_studies, drop = FALSE]
                                   ##
                                   ## Sum across holdout studies and average across iterations/chains
                                   elpd_holdout <- mean(apply(log_lik_holdout_samples, c(1, 2), sum))
                                   ##
                                   ### Update fold result
                                   fold_result$status <- "completed"
                                   fold_result$elpd <- elpd_holdout
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
                                   ##
                                   ## Optionally save fold-specific results
                                   saveRDS(fold_result, 
                                           file = file.path(output_dir, sprintf("fold_%02d_results.rds", k)))
                                   
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
  
  
  # Stop cluster
  stopCluster(cl)
  
  # Clean up temporary files
  for (k in 1:K) {
    fold_file <- file.path(output_dir, paste0("cv_environment_fold_", k, ".RData"))
    if (file.exists(fold_file)) {
      file.remove(fold_file)
    }
  }
  
  # Process results
  successful_folds <- Filter(function(x) !is.null(x) && x$status == "completed", 
                             fold_results_list)
  
  if (length(successful_folds) < K) {
    warning(sprintf("Only %d/%d folds completed successfully", 
                    length(successful_folds), K))
  }
  
  # Calculate total elpd_kfold
  elpd_kfold <- sum(sapply(successful_folds, function(x) x$elpd))
  
  # Compute SE using fold-level results
  fold_elpds <- sapply(successful_folds, function(x) x$elpd)
  se_elpd <- sd(fold_elpds) * sqrt(K)
  
  # Summary of timing
  total_time <- sum(sapply(successful_folds, function(x) x$duration_mins))
  parallel_time <- max(sapply(successful_folds, function(x) x$duration_mins))
  ##
  min_ESS_vec <- max_rhat_vec <- c()
  for (k in 1:length(successful_folds)) {
    min_ESS_vec[k]  <- successful_folds[[k]]$min_ESS
    max_rhat_vec[k] <- successful_folds[[k]]$max_rhat
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
    K = K,
    ##
    min_ESS_vec = min_ESS_vec,
    max_rhat_vec = max_rhat_vec,
    ##
    fold_results = successful_folds,
    fold_assignments = fold_ids,
    n_workers = n_workers,
    timing = list(
      total_cpu_mins = total_time,
      wall_clock_mins = parallel_time,
      speedup = total_time / parallel_time
    ),
    errors = Filter(function(x) !is.null(x) && x$status == "error", 
                    fold_results_list)
  ))
  
}















#' R_fn_run_kfold_cv_auto
#' Wrapper function to choose between sequential and parallel k-fold CV
#' @export
R_fn_run_kfold_cv_auto <- function(  debugging,
                                     ##
                                     K = 10, 
                                     ##
                                     model_prep_obj,
                                     stan_data_list, 
                                     ##
                                     priors,
                                     init_lists_per_chain,
                                     n_burnin,
                                     n_iter,
                                     n_chains,
                                     adapt_delta,
                                     max_treedepth,
                                     ##
                                     parallel,
                                     n_workers,
                                     output_dir,
                                     ##
                                     use_BayesMVP_for_faster_summaries
) {
  
        # Decide whether to use parallel version
        if (parallel == "auto") {
          # Use parallel if we have enough cores for at least 2 workers
          available_cores <- parallel::detectCores()
          potential_workers <- floor(available_cores / n_chains)
          use_parallel <- potential_workers >= 2 && K >= 4
        } else {
          use_parallel <- parallel
        }
        
        if (use_parallel) {
              cat("Using PARALLEL k-fold CV\n")
              return(R_fn_run_kfold_cv_parallel(
                debugging = debugging,
                K = K,
                model_prep_obj = model_prep_obj,
                stan_data_list = stan_data_list,
                priors = priors,
                init_lists_per_chain = init_lists_per_chain,
                n_burnin = n_burnin,
                n_iter = n_iter,
                n_chains = n_chains,
                adapt_delta = adapt_delta,
                max_treedepth = max_treedepth,
                n_workers = n_workers,
                output_dir = output_dir,
                use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries
              ))
        } else {
              cat("Using SEQUENTIAL k-fold CV\n")
              return(R_fn_run_kfold_cv(
                debugging = debugging,
                K = K,
                model_prep_obj = model_prep_obj,
                stan_data_list = stan_data_list,
                priors = priors,
                init_lists_per_chain = init_lists_per_chain,
                n_burnin = n_burnin,
                n_iter = n_iter,
                n_chains = n_chains,
                adapt_delta = adapt_delta,
                max_treedepth = max_treedepth,
                use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries
              ))
        }
  
}















##
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
  
  # Filter fold results
  filtered_fold_results <- kfold_results$fold_results[good_folds]
  
  # Recompute ELPD (sum of fold ELPDs)
  fold_elpds <- sapply(filtered_fold_results, function(x) x$elpd)
  new_elpd_kfold <- sum(fold_elpds)
  
  # Recompute SE 
  # Note: This is tricky - we're using fewer folds now
  # Standard approach: SE = sd(fold_elpds) * sqrt(n_folds)
  new_se_elpd <- sd(fold_elpds) * sqrt(length(good_folds))
  
  # Create filtered results object
  filtered_results <- list(
    elpd_kfold = new_elpd_kfold,
    se_elpd = new_se_elpd,
    K_original = kfold_results$K,
    K_used = length(good_folds),
    min_ESS_vec = ess_values[good_folds],
    max_rhat_vec = kfold_results$max_rhat_vec[good_folds],
    fold_results = filtered_fold_results,
    excluded_folds = bad_folds,
    min_ess_threshold = min_ess_threshold
  )
  
  if (verbose) {
    cat(sprintf("\nFiltered results (using %d/%d folds):\n", 
                length(good_folds), kfold_results$K))
    cat(sprintf("  ELPD: %.1f (SE: %.1f)\n", new_elpd_kfold, new_se_elpd))
    cat(sprintf("  Change in ELPD: %.1f\n", new_elpd_kfold - kfold_results$elpd_kfold))
  }
  
  return(filtered_results)
}












# Quick comparison function:
compare_kfold_thresholds <- function(kfold_results, thresholds = c(1, 50, 100, 200, 400, 600, 650, 700, 704, 706, 725, 750, 800, 1000)) {
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








compare_kfold_models <- function(kfold_results_list, 
                                 model_names = NULL,
                                 min_ess_threshold = 100,
                                 verbose = TRUE) {
  
  n_models <- length(kfold_results_list)
  
  # Set model names if not provided
  if (is.null(model_names)) {
    model_names <- paste0("Model_", 1:n_models)
  }
  
  # 1. Check that all models used the same fold assignments
  if (verbose) cat("Checking fold consistency across models...\n")
  
  # Get fold assignments from first model
  ref_fold_assignments <- kfold_results_list[[1]]$fold_assignments
  ref_K <- kfold_results_list[[1]]$outs_kfold$K
  
  for (i in 2:n_models) {
    if (!identical(kfold_results_list[[i]]$fold_assignments, ref_fold_assignments)) {
      stop(sprintf("Model %d (%s) has different fold assignments than Model 1!", 
                   i, model_names[i]))
    }
    print(kfold_results_list[[i]]$outs_kfold$K)
    ##
    if (kfold_results_list[[i]]$outs_kfold$K != ref_K) {
      stop(sprintf("Model %d (%s) has different K (%d) than Model 1 (%d)!", 
                   i, model_names[i], kfold_results_list[[i]]$K, ref_K))
      ##
      
    }
  }
  
  if (verbose) cat("✓ All models used identical fold assignments\n\n")
  
  # Check held-out studies match (extra safety check)
  for (fold in 1:ref_K) {
    ref_holdout <- kfold_results_list[[1]]$outs_kfold$fold_results[[fold]]$holdout_studies
    for (i in 2:n_models) {
      model_holdout <- kfold_results_list[[i]]$outs_kfold$fold_results[[fold]]$holdout_studies
      if (!identical(sort(ref_holdout), sort(model_holdout))) {
        stop(sprintf("Fold %d has different holdout studies between models!", fold))
      }
    }
  }
  
  # 2. Find ALL folds with low ESS across ANY model
  all_min_ess <- matrix(NA, nrow = ref_K, ncol = n_models)
  colnames(all_min_ess) <- model_names
  
  for (i in 1:n_models) {
    all_min_ess[, i] <- kfold_results_list[[i]]$outs_kfold$min_ESS_vec
  }
  
  # Find minimum ESS for each fold across all models
  min_ess_per_fold <- apply(all_min_ess, 1, min, na.rm = TRUE)
  
  # Identify folds to exclude (any model has ESS below threshold)
  bad_folds <- which(min_ess_per_fold < min_ess_threshold)
  good_folds <- which(min_ess_per_fold >= min_ess_threshold)
  
  if (verbose) {
    cat("ESS summary across models:\n")
    print(round(all_min_ess, 1))
    cat("\nMinimum ESS per fold:", round(min_ess_per_fold, 1), "\n")
    
    if (length(bad_folds) > 0) {
      cat(sprintf("\nExcluding %d fold(s) with ESS < %d in ANY model:\n", 
                  length(bad_folds), min_ess_threshold))
      for (fold in bad_folds) {
        cat(sprintf("  Fold %d: min ESS = %.1f\n", fold, min_ess_per_fold[fold]))
      }
    }
  }
  
  # 3. Filter all models using the SAME good folds
  filtered_results <- list()
  
  for (i in 1:n_models) {
    # Get fold results for good folds only
    filtered_fold_results <- kfold_results_list[[i]]$outs_kfold$fold_results[good_folds]
    
    # Recompute ELPD and SE
    fold_elpds <- sapply(filtered_fold_results, function(x) x$elpd)
    new_elpd <- sum(fold_elpds)
    new_se <- sd(fold_elpds) * sqrt(length(good_folds))
    
    filtered_results[[model_names[i]]] <- list(
      elpd_kfold = new_elpd,
      se_elpd = new_se,
      K_original = ref_K,
      K_used = length(good_folds),
      fold_elpds = fold_elpds,
      excluded_folds = bad_folds,
      min_ESS_vec = kfold_results_list[[i]]$min_ESS_vec[good_folds]
    )
  }
  
  # 4. Create comparison table
  comparison_table <- data.frame(
    Model = model_names,
    ELPD = sapply(filtered_results, function(x) x$elpd_kfold),
    SE = sapply(filtered_results, function(x) x$se_elpd),
    Delta_ELPD = NA,
    Delta_SE = NA,
    stringsAsFactors = FALSE
  )
  
  # Calculate differences from first model
  for (i in 2:n_models) {
    comparison_table$Delta_ELPD[i] <- comparison_table$ELPD[i] - comparison_table$ELPD[1]
    # SE of difference (assuming independence - conservative)
    comparison_table$Delta_SE[i] <- sqrt(comparison_table$SE[i]^2 + comparison_table$SE[1]^2)
  }
  
  if (verbose) {
    cat(sprintf("\n\nFinal comparison using %d/%d folds:\n", length(good_folds), ref_K))
    print(comparison_table, digits = 1)
  }
  
  # Return everything
  return(list(
    filtered_results = filtered_results,
    comparison_table = comparison_table,
    excluded_folds = bad_folds,
    K_used = length(good_folds),
    min_ess_matrix = all_min_ess
  ))
}




# Quick summary function
summarize_model_comparison <- function(comparison_results) {
  cat("\nModel ranking by ELPD:\n")
  table <- comparison_results$comparison_table
  table <- table[order(table$ELPD, decreasing = TRUE), ]
  
  for (i in 1:nrow(table)) {
    if (i == 1) {
      cat(sprintf("%d. %s: ELPD = %.0f (%.0f) [BEST]\n", 
                  i, table$Model[i], table$ELPD[i], table$SE[i]))
    } else {
      sig <- ifelse(abs(table$Delta_ELPD[i]) > 2*table$Delta_SE[i], "*", "")
      cat(sprintf("%d. %s: ELPD = %.0f (%.0f), Δ = %.0f (%.0f) %s\n", 
                  i, table$Model[i], table$ELPD[i], table$SE[i],
                  table$Delta_ELPD[i], table$Delta_SE[i], sig))
    }
  }
  cat("\n* = Significantly different from best model (|Δ| > 2*SE)\n")
}

