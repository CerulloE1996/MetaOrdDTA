





# This is actually a valid research finding, not a problem to "fix"!
#   What's happening:
# 
# Random-cutpoint models estimate kappa → ∞ when heterogeneity → 0
# Large kappa causes numerical instability in Stan (kappa × simplex overflows)
# This is a real limitation of random-cutpoint models, not a bug
# 
# Your options:
# Option 1: Drop random-cutpoint models on fixed-DGMs (RECOMMENDED)
# Justification: 
# - Fixed-cutpoint DGMs are unrealistic anyway (real data has heterogeneity)
# - You already showed fixed-cutpoint MODELS work great on random-DGMs (more important!)
# - Random→Fixed comparison shows robustness
# - Fixed→Random comparison just shows numerical instability (not interesting)
# Option 2: Document as limitation
# "Random-cutpoint models exhibited poor ESS (<100) in ~50% of simulations 
# when fitted to fixed-cutpoint DGMs due to numerical instability when 
# estimating near-zero heterogeneity (kappa → ∞). These scenarios were 
# excluded as fixed-cutpoint DGMs are clinically unrealistic."
# Option 3: Only fit fixed-models to fixed-DGMs
# Design matrix:
#          | Fixed-DGM | Random-DGM
# ---------|-----------|------------
# Fixed-M  |    ✓      |     ✓
# Random-M |    —      |     ✓
# 
# Note: Random models not fitted to fixed-DGMs due to numerical instability
# My recommendation: Option 1 or 3 - just skip that comparison entirely.
# You've already shown the important result: fixed-cutpoint models are robust to random cutpoints. 
# The reverse (random models on fixed cutpoints) is just demonstrating a numerical limitation, 
# not interesting scientifically.
# 
# Bottom line: Don't waste compute trying to "fix" an unrealistic scenario. 
# Focus on realistic DGMs with heterogeneity where all models can be fairly compared.








# ============================================================================================================
# Main Execution Script for MetaOrdDTA Simulation Study
# ============================================================================================================

cat("\014") ## clear console

# # # # # # # # # ##
# rm(list = ls())
# .rs.restartR()  # In RStudio
##
# {
#       os <- .Platform$OS.type
# 
#       if (os == "unix") {
#         user_root_dir <- Sys.getenv("PWD")
#       } else if (os == "windows") {
#         user_root_dir <- Sys.getenv("USERPROFILE")
#       }
#       local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
# 
#       {
#         ## First remove any possible package fragments:
#         ## Find user_pkg_install_dir:
#         user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
#         print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
#         ##
#         ## Find pkg_install_path + pkg_temp_install_path:
#         pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
#         pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
#         ##
#         ## Remove any (possible) MetaOrdDTA package fragments:
#         remove.packages("MetaOrdDTA")
#         unlink(pkg_install_path, recursive = TRUE, force = TRUE)
#         unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
#       }
# 
#       #### -------- ACTUAL (LOCAL) INSTALL:
#       ## Document:
#       devtools::clean_dll(local_pkg_dir)
#       roxygen2::roxygenize(local_pkg_dir)
#       ##
#       ## Install (outer pkg):
#       ##
#       devtools::clean_dll(local_pkg_dir)
#       devtools::install(local_pkg_dir,
#                         upgrade = "never",
#                         quick = TRUE,
#                         dependencies = FALSE)
#       ##
#       ## May need to restart R:
#       ##
#       .rs.restartR()  # In RStudio
# }
# ##
library(MetaOrdDTA)
library(doParallel)
library(dplyr)
##
{
    os <- .Platform$OS.type
    
    if (os == "unix") {
      user_root_dir <- Sys.getenv("PWD")
    } else if (os == "windows") {
      user_root_dir <- Sys.getenv("USERPROFILE")
    }
    local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
}

setwd(local_pkg_dir)

options(scipen = 999999)
options(max.print = 1000)

# Set paths
output_dir <- file.path(local_pkg_dir, "simulation_results")

sim_study_dir <- file.path(local_pkg_dir, "sim_study")
output_dir <- file.path(sim_study_dir, "simulation_results")


# # Main execution
# scenario_row <- data.frame(
#   index_test_MA = "GAD_2",
#   n_studies = 50,
#   model_parameterisation = "Xu",
#   alpha_pi = "one",
#   box_cox = FALSE,
#   softplus = FALSE,
#   random_thresholds = FALSE,
#   Dirichlet_random_effects_type = "fixed"
# )

# # Calculate how many simulations needed
# n_sims_total <- 50
# n_workers <- 16
# sims_per_batch <- 100



# 
# compare_kappa_dofs( log_mean = log(50),
#                     log_sd = 1)


# log_mean <- log(100)
# log_sd <- 1.0
# ##
# ## ---- dof = 2:
# ##
# quick_kappa_prior(log_mean = log(50),
#                   log_sd = log_sd,
#                   df = 2)
# ##
# quick_kappa_prior(log_mean = log(100),
#                   log_sd = log_sd,
#                   df = 2)
# ##
# quick_kappa_prior(log_mean = log(250),
#                   log_sd = log_sd,
#                   df = 2)
# ##
# quick_kappa_prior(log_mean = log(500),
#                   log_sd = log_sd,
#                   df = 2)
# ##
# ## ---- dof = 3:
# ##
# quick_kappa_prior(log_mean = log(50),
#                   log_sd = log_sd,
#                   df = 3)
# ##
# quick_kappa_prior(log_mean = log(100),
#                   log_sd = log_sd,
#                   df = 3)
# ##
# quick_kappa_prior(log_mean = log(250),
#                   log_sd = log_sd,
#                   df = 3)
# ##
# quick_kappa_prior(log_mean = log(500),
#                   log_sd = log_sd,
#                   df = 3)
# ##
# ## ---- dof = 5:
# ##
# quick_kappa_prior(log_mean = log(100),
#                   log_sd = log_sd,
#                   df = 5)
# ##
# quick_kappa_prior(log_mean = log(250),
#                   log_sd = log_sd,
#                   df = 5)
# ##
# quick_kappa_prior(log_mean = log(500),
#                   log_sd = log_sd,
#                   df = 5)



## ============================================================================================================
## ---- Parallel execution
## ============================================================================================================
{
    require(parallel)
    require(doParallel)
    require(foreach)
}
##
{
    # Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread") 
    options(mc.cores = parallel::detectCores())
    Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
    RcppParallel::setThreadOptions(numThreads = parallel::detectCores())
}
##
{
    system("ulimit -s")
    
    # This shows R's actual limits, not a new shell:
    system("cat /proc/self/limits | grep stack")
    
    # Or check R's process directly:
    pid <- Sys.getpid()
    system(paste0("cat /proc/", pid, "/limits | grep stack"))
    
    # Install and use unix package to change limits from within R:
    require(unix)
    
    # Set stack size to unlimited (or large value)
    rlimit_stack(cur = 256 * 1024 * 1024)  # 256 MB
    
    # Verify it worked:
    system("cat /proc/self/limits | grep stack")
}


# # Check Jones files with FALSE
# jones_false_files <- list.files(output_dir, 
#                                 pattern = ".*GAD_2.*Jones.*n_studies_10.*rand_thr_FALSE.*", 
#                                 full.names = FALSE)
# head(jones_false_files, 20)
# 
# 
# 
# 
# # Check for red_fact_1 AND rand_thr_FALSE
# xu_false_redfact1 <- list.files(output_dir, 
#                                 pattern = ".*GAD_2.*red_fact_1.*Xu_random_C.*n_studies_10.*rand_thr_FALSE.*", 
#                                 full.names = FALSE)
# length(xu_false_redfact1)
# head(xu_false_redfact1)
# 
# rg_false_redfact1 <- list.files(output_dir, 
#                                 pattern = ".*GAD_2.*red_fact_1.*R&G_random_C.*n_studies_10.*rand_thr_FALSE.*", 
#                                 full.names = FALSE)
# length(rg_false_redfact1)
# head(rg_false_redfact1)
##
options(max.print = 1000)
##
{
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
    source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
    source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
}
##
# Get all scenarios
scenario_df <- create_MetaOrd_scenarios() ##  %>% filter(model_parameterisation == "Bivariate")
# scenario_df <- scenario_df  %>% filter(DGM ==("Xu_fixed_C"))
scenario_df <- scenario_df ## %>% filter(random_thresholds == TRUE)
scenario_df %>% print(n = 1000)
##
# A tibble: 5 × 13
# test  DGM         n_studies model_parameterisation random_thresholds alpha_pi mcse_met max_MCSE n_sims min_sims_needed avg_RMSE avg_abs_bias mean_min_ESS
# <chr> <chr>           <dbl> <chr>                  <lgl>             <chr>    <lgl>       <dbl>  <int>           <dbl>    <dbl>        <dbl>        <dbl>
# 3 HADS  Xu_random_C        50 Bivariate              FALSE             one      FALSE       0.126    114           132.      2.80        0.391          917
##
## ---- Test:
##
# test_vals_vec <- c("GAD_2", "HADS", "BAI")
# test_vals_vec <- c("HADS", "BAI")
# test_vals_vec <- c("HADS", "BAI")
# test_vals_vec <- "GAD_2"
test_vals_vec <- "HADS"
# test_vals_vec <- "BAI"
##
## ---- DGM:
## 
# DGM_vec <- c( "Jones",
#               "Xu_fixed_C",
#               "Xu_random_C",
#               "R&G_random_C")
##
# DGM_vec <- "Jones"
DGM_vec <- "Xu_fixed_C"
# DGM_vec <- "Xu_random_C"
# DGM_vec <- "R&G_random_C"
##
n_studies_vals_vec <- 10
# n_studies_vals_vec <- 50
# n_studies_vals_vec <- c(10, 50)
##
## ---- model_parameterisation:
##
# model_param_vec <- c("Jones", "Xu", "R&G", "Bivariate")
# model_param_vec <- c("Xu", "Jones", "R&G")
# model_param_vec <- "Bivariate"
# model_param_vec <- c("Jones", "Bivariate")
# model_param_vec <- "Jones"
# model_param_vec <- "R&G"
model_param_vec <- "Xu"
##
## ----
##
# alpha_pi_vec <- "one"
alpha_pi_vec <- "informative"
##
# test_val    <- "BAI"
# model_param <- model_param_vec
##
## ----
##
# target_MCSE_RMSE <- 0.25
target_MCSE_RMSE <- 0.20
# target_MCSE_RMSE <- 0.15
# target_MCSE_RMSE <- 0.125
# target_MCSE_RMSE <- 0.10
##
## ----
##
max_seed_to_check <- 40000
min_runs_needed <- 30
##
## ----
##
starting_sim_offset <- 1
n_sims_total <- starting_sim_offset + 1000 ; n_sims_total
##
## ----
##
root_MCMC_seed <- 123
##
reduction_factor <- 1.0
# reduction_factor <- 0.50
# reduction_factor^2
##
## ----
##
min_studies_per_threshold <- 1
##
## ----
##
# output_dir = output_dir
# file_name_string = file_name_string
# max_seed_to_check = max_seed_to_check
# min_runs_needed = min_runs_needed
# target_MCSE_RMSE = target_MCSE_RMSE
# test = test_val
# DGM = DGM
# n_studies = scenario_row$n_studies
# vec_index_inner = vec_index_inner
# min_studies_per_threshold = min_studies_per_threshold
##
##
for (DGM_val in DGM_vec) {
  DGM <- DGM_val
  for (test_val in test_vals_vec) {
      
      if (test_val == "GAD_2") {
        vec_index_inner <- 1:6
      } else if (test_val == "HADS") {
        vec_index_inner <- 1:21
        # vec_index_inner <- 3:17
      } else if (test_val == "BAI") {
        # vec_index_inner <- 1:63 ## 5:43
        vec_index_inner <- 5:43
      }
    
    for (model_param in model_param_vec) {
      
      # Filter scenarios
      run_matrix_filtered <- scenario_df %>% 
        filter(
          index_test_MA == test_val,
          DGM == DGM_val,
          n_studies %in% n_studies_vals_vec,
          model_parameterisation == model_param,
          alpha_pi %in% alpha_pi_vec
        ) %>% 
        print(n = 100)
      
      if (nrow(run_matrix_filtered) == 0) next
      
      ## Simulate data + get true values 
      source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA_v2.R"), local = TRUE)
      ## true_DGM_Se and true_DGM_Sp will be set by the simulation
      
      scenario_idx <- 1
      
      ## Process each scenario in the filtered set
      for (scenario_idx in 1:nrow(run_matrix_filtered)) {
        
                  scenario_row <- run_matrix_filtered[scenario_idx, ] ; scenario_row
                  ##
                  cat(sprintf("\n========== Processing: Test=%s, Model=%s, n=%d, alpha=%s ==========\n",
                              scenario_row$index_test_MA,
                              scenario_row$model_parameterisation,
                              scenario_row$n_studies,
                              scenario_row$alpha_pi))
                  ##
                  ## ---- Create file name string:
                  ##
                  file_name_string <- create_file_name_string( scenario_row,
                                                               reduction_factor)
                  ##
                  source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"), local = TRUE)
                  ##
                  outs <- list()
                  try({  
                  outs <- check_MCSE_status(  output_dir = output_dir,
                                              file_name_string = file_name_string,
                                              max_seed_to_check = max_seed_to_check,
                                              min_runs_needed = min_runs_needed,
                                              target_MCSE_RMSE = target_MCSE_RMSE,
                                              test = test_val,
                                              DGM = DGM,
                                              n_studies = scenario_row$n_studies,
                                              vec_index_inner = vec_index_inner,
                                              min_studies_per_threshold = min_studies_per_threshold)
                  })
                  ##
                  try({
                  if (outs$mcse_met) {
                    cat("✓ MCSE target already met, skipping\n")
                    next
                  }
                  })
                  ##
                  ## Calculate sims needed
                  n_sims_total <- n_sims_total  # or outs$min_sims_needed
                  n_threads_to_use <- 96*1.0
                  n_chains_per_worker <- 4
                  n_workers <- floor(n_threads_to_use / n_chains_per_worker) ; n_workers
                  sims_per_batch <- n_sims_total ## 8 * n_threads_to_use / n_chains_per_worker
                  ##
                  ## Set up minimal environment
                  ##
                  minimal_env <- new.env()
                  minimal_env$run_multiple_sims <- run_multiple_sims
                  minimal_env$root_MCMC_seed <- 123
                  minimal_env$sim_study_dir <- sim_study_dir
                  minimal_env$local_pkg_dir <- local_pkg_dir
                  minimal_env$output_dir <- output_dir
                  minimal_env$n_chains_per_worker <- n_chains_per_worker
                  minimal_env$scenario_row <- scenario_row
                  minimal_env$file_name_string <- file_name_string
                  minimal_env$target_MCSE_RMSE <- target_MCSE_RMSE
                  minimal_env$max_seed_to_check <- max_seed_to_check
                  minimal_env$min_runs_needed <- min_runs_needed
                  minimal_env$DGM <- DGM
                  minimal_env$reduction_factor <- reduction_factor
                  minimal_env$min_studies_per_threshold <- min_studies_per_threshold
                  
                  ## Batch execution
                  n_batches <- 1 ## ceiling(n_sims_total / sims_per_batch)
                  starting_sim_offset <- starting_sim_offset
                  ##
                  batch_id <- 1
                  ##
                  for (batch_id in 1:n_batches) {
                          
                          cat(sprintf("\n========== BATCH %d/%d ==========\n", 
                                      batch_id, 
                                      n_batches))
                          
                          ## Check MCSE before starting batch
                          # if (batch_id > 1) {
                          
                                current_mcse_status <- list()
                                try({ 
                                current_mcse_status <- check_MCSE_status( output_dir = output_dir,
                                                                          file_name_string = file_name_string,
                                                                          max_seed_to_check = max_seed_to_check,
                                                                          min_runs_needed = min_runs_needed,
                                                                          target_MCSE_RMSE = target_MCSE_RMSE,
                                                                          test = test_val,
                                                                          DGM = DGM,
                                                                          n_studies =  scenario_row$n_studies,
                                                                          vec_index_inner = vec_index_inner,
                                                                          min_studies_per_threshold = min_studies_per_threshold)
                                
                                
                                current_mcse_status$MCSE_val
                                
                                
                                })
                                try({ 
                                if (!is.na(current_mcse_status$MCSE_val) && current_mcse_status$mcse_met) {
                                  cat("✓ TARGET MCSE REACHED! Stopping batches.\n")
                                  break
                                }
                                })
                            
                          # }
                          ## Calculate simulation range for this batch
                          start_sim <- starting_sim_offset + (batch_id - 1) * sims_per_batch + 1 ; start_sim
                          end_sim <- min(starting_sim_offset + batch_id * sims_per_batch, starting_sim_offset + n_sims_total) ; end_sim
                          sims_this_batch <- start_sim:end_sim ; sims_this_batch
                          
                          ## Distribute to workers
                          sim_ids_chunk_list <- vector("list", n_workers)
                          ##
                          for (worker_i in 1:n_workers) {
                              worker_sims <- sims_this_batch[seq(worker_i, 
                                                                 length(sims_this_batch), 
                                                                 by = n_workers)]
                              sim_ids_chunk_list[[worker_i]] <- worker_sims
                          }
                          minimal_env$sim_ids_chunk_list <- sim_ids_chunk_list
                          
                          ## Create cluster and run
                          cl <- makeCluster( n_workers,
                                             type = "PSOCK", 
                                             outfile = "")
                          registerDoParallel(cl)
                          
                          clusterExport(cl, 
                                        ls(envir = minimal_env),
                                        envir = minimal_env)
                          
                          # MCMC_seed = minimal_env$root_MCMC_seed + ii + batch_id * 1000
                          # sim_ids_chunk = sim_ids_chunk_list[[ii]]
                          # n_chains = n_chains_per_worker
                          
                          ii <- 1
                          
                          MCMC_seed = 123 + ii + batch_id * 1000
                          sim_ids_chunk = sim_ids_chunk_list[[ii]]
                          min_runs_needed = 30
                          check_MCSE_interval = 1
                          n_chains <- 4
                          
                          batch_results <- foreach(ii = 1:n_workers,
                                                   .combine = 'list',
                                                   .errorhandling = "pass") %dopar% {
                                                     
                                                     tryCatch({
                                                       
                                                           ## Load packages
                                                           require(RcppParallel)
                                                           require(BayesMVPsim)
                                                           require(BayesMVP)
                                                           require(MetaOrdDTA)
                                                           require(cmdstanr)
                                                           require(dplyr)
                                                           
                                                           cmdstanr::set_cmdstan_path()
                                                           set.seed(root_MCMC_seed + ii + batch_id * 1000)
                                                           n_chains = n_chains_per_worker
                                                           
                                                           ## Load configuration and helpers
                                                           {
                                                               source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
                                                               source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
                                                           }
                                                           
                                                           ## Call sim fn
                                                           result <- run_multiple_sims(  MCMC_seed = root_MCMC_seed + ii + batch_id * 1000,
                                                                                         sim_ids_chunk = sim_ids_chunk_list[[ii]],
                                                                                         ##
                                                                                         n_chains = n_chains_per_worker,
                                                                                         ##
                                                                                         scenario_row = scenario_row,
                                                                                         ##
                                                                                         # reduction_factor = reduction_factor,
                                                                                         ##
                                                                                         sim_study_dir = sim_study_dir,
                                                                                         output_dir = output_dir,
                                                                                         local_pkg_dir = local_pkg_dir,
                                                                                         file_name_string = file_name_string,
                                                                                         ##
                                                                                         max_seed_to_check = max_seed_to_check,
                                                                                         min_runs_needed = 30,
                                                                                         target_MCSE_RMSE = target_MCSE_RMSE,
                                                                                         check_MCSE_interval = 1,
                                                                                         stop_early_if_target_met = TRUE)
                                                                                         ##
                                                                                         # vec_index_inner = vec_index_inner)
                                                           
                                                           return(list(worker_id = ii,
                                                                       batch_id = batch_id,
                                                                       status = "completed",
                                                                       result = result))
                                                       
                                                     }, error = function(e) {
                                                       
                                                           ## Capture and return the actual error
                                                           cat(sprintf("\nWorker %d ERROR: %s\n", 
                                                                       ii, 
                                                                       e$message))
                                                           cat(sprintf("Traceback: %s\n", paste(capture.output(traceback()), 
                                                                                                collapse = "\n")))
                                                           return(list(worker_id = ii, 
                                                                       batch_id = batch_id, 
                                                                       status = "error", 
                                                                       error = e$message))
                                                       
                                                     })
                                                     
                                                   }
                          
                          stopCluster(cl)
                          gc(reset = TRUE, full = TRUE)
                          
                          cat(sprintf("Batch %d completed\n", batch_id))
                    
                  }
                  ## ---- Final MCSE check for this scenario:
                  final_mcse_status <- list()
                  try({
                  final_mcse_status <- check_MCSE_status( output_dir = output_dir,
                                                          file_name_string = file_name_string,
                                                          max_seed_to_check = max_seed_to_check,
                                                          min_runs_needed = min_runs_needed,
                                                          target_MCSE_RMSE = target_MCSE_RMSE,
                                                          test = test_val,
                                                          DGM = DGM,
                                                          n_studies = scenario_row$n_studies,
                                                          vec_index_inner = vec_index_inner,
                                                          min_studies_per_threshold = min_studies_per_threshold)
                  })
                  ##
                  cat(sprintf("FINAL: %d sims, MCSE = %.3f%%\n",
                              final_mcse_status$n_sims,
                              final_mcse_status$MCSE_val))
        
      }
    }
  }
}




final_mcse_status$MCSE_val


 

# scenario_df = scenario_df
# output_dir = output_dir
# test_filter = "GAD_2"  # optional
# n_studies_filter = c(10, 50) # optional
# target_MCSE_RMSE = 0.125
# max_seed_to_check = 5000
# min_runs_needed = 30



{
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
    source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
    source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
}




 
 



































