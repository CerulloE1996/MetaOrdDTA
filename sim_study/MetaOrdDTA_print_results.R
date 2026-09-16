





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



# # # # # # # # ##
# rm(list = ls())
# .rs.restartR()  # In RStudio
# 
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



library(MetaOrdDTA)
library(doParallel)
library(dplyr)

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
sim_study_dir <- file.path(local_pkg_dir, "sim_study")
output_dir <- MetaOrdDTA::metaorddta_output_path(
  "sim_study", "simulation_results"
)


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
require(parallel)
require(doParallel)
require(foreach)
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
scenario_df %>% print(n = 100)
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
test_vals_vec <- "BAI"
##
## ---- DGM:
## 
DGM_vec <- c("Jones", "Xu_fixed_C", "Xu_random_C", "R&G_random_C")
# DGM_vec <- "Jones"
# DGM_vec <- "Xu_fixed_C"
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
model_param_vec <- "Bivariate"
# model_param_vec <- c("Jones", "Bivariate")
# model_param_vec <- "Jones"
# model_param_vec <- "R&G"
##
## ----
##
alpha_pi_vec <- "one"
##
# test_val    <- "BAI"
# model_param <- model_param_vec
##
## ----
##
target_MCSE_RMSE <- 0.25
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
starting_sim_offset <- 10000
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
## ---- Create scenarios:
##
scenario_df <- create_MetaOrd_scenarios() #### %>% filter(model_parameterisation == "Bivariate")
scenario_df
##
scenario_df <- scenario_df %>% filter(model_parameterisation == "Xu") ## --------------------------
##
reduction_factor   <- 1.0
# reduction_factor <- 0.5
##
output_dir <- MetaOrdDTA::metaorddta_output_path(
  "sim_study", "simulation_results"
)
##
## ---- GAD-2; n_studies = 10: -----------------------------------------------------
##
{
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
    ##
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
    ##
    source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
    source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
    source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
}
##
# min_studies_per_threshold <- 0
min_studies_per_threshold <- 3
##
results_GAD_2_n_studies_10 <- process_all_scenarios_MetaOrd(
                              reduction_factor = reduction_factor,
                              scenario_df = scenario_df,
                              output_dir  = output_dir,
                              test_filter = "GAD_2",
                              n_studies_filter = c(10),
                              target_MCSE_RMSE = target_MCSE_RMSE,
                              max_seed_to_check = max_seed_to_check,
                              min_runs_needed = 5,
                              min_ESS_thr = 3000,
                              vec_index_inner = 1:6,
                              min_studies_per_threshold = min_studies_per_threshold)
##
results_GAD_2_n_studies_10_summary <- filter( results_GAD_2_n_studies_10$summary_table,
                                              n_studies == 10)
##
results_GAD_2_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
##
## ---- GAD-2; n_studies = 50:  -----------------------------------------------------
##
##
results_GAD_2_n_studies_50 <- process_all_scenarios_MetaOrd(
                              reduction_factor = reduction_factor,
                              scenario_df = scenario_df,
                              output_dir  = output_dir,
                              test_filter = "GAD_2",
                              n_studies_filter  = c(50),
                              target_MCSE_RMSE  = target_MCSE_RMSE,
                              max_seed_to_check = max_seed_to_check,
                              min_runs_needed = 5,
                              min_ESS_thr = 400,
                              vec_index_inner = 1:6,
                              min_studies_per_threshold = min_studies_per_threshold)
##
results_GAD_2_n_studies_50_summary <- filter( results_GAD_2_n_studies_50$summary_table,
                                              n_studies == 50)
##
results_GAD_2_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_GAD_2_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE, avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
##
##  --------------------------------------------------------------------
##
# metric <- "Bias"
metric <- "RMSE"
##
create_MetaOrd_summary_3( ## ---- Jones > (Bivariate, Xu)  ## ---- RMSE diff (Jones vs Xu): ~ 0.50
  data =  filter(results_GAD_2_n_studies_10_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- (Xu, Jones) > Bivariate
  data =  filter(results_GAD_2_n_studies_10_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- (Xu, Xu, Jones, R&G) > Bivariate
  data =  filter(results_GAD_2_n_studies_10_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- (Xu, Jones, R&G, Xu) > Bivariate
  data =  filter(results_GAD_2_n_studies_10_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
## ----
##
create_MetaOrd_summary_3( ## ---- (Jones, Xu) > Bivariate
  data =  filter(results_GAD_2_n_studies_50_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- (Xu, Bivariate) > Jones ## ---- RMSE diff (Jones vs Xu): ~ 0.56
  data =  filter(results_GAD_2_n_studies_50_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- (Xu, Xu, R&G) > (Jones, Bivariate) ## ---- RMSE diff (Jones vs Xu): ~ 0.60 (0.44 for Xu-Rand vs. Jones)
  data =  filter(results_GAD_2_n_studies_50_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3( ## ---- All 5 models tied
  data =  filter(results_GAD_2_n_studies_50_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)






##
## ---- Create scenarios:
##
scenario_df <- create_MetaOrd_scenarios() ### %>% filter(random_thresholds == TRUE)
# ##
# reduction_factor <- 1.0
# # reduction_factor <- 0.5
##
## ----HADS; n_studies = 10: -----------------------------------------------------------------------------------------------------------
##
# Process all or filtered scenarios
results_HADS_n_studies_10 <-  process_all_scenarios_MetaOrd(
                              reduction_factor = reduction_factor,
                              scenario_df = scenario_df,
                              output_dir  = output_dir,
                              test_filter = "HADS",
                              n_studies_filter = c(10),
                              target_MCSE_RMSE = 0.20,
                              max_seed_to_check = max_seed_to_check,
                              min_runs_needed = 5,
                              min_ESS_thr = 3000,
                              vec_index_inner = 1:21,
                              min_studies_per_threshold = min_studies_per_threshold)
##
results_HADS_n_studies_10_summary <- filter( results_HADS_n_studies_10$summary_table,
                                             n_studies == 10,
                                             alpha_pi != "flat",
                                             alpha_pi != "custom")
##
results_HADS_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
# ##
# # A tibble: 2 × 13
#   test  DGM        n_studies model_parameterisation random_thresholds alpha_pi    mcse_met MCSE_val n_sims min_sims_needed avg_RMSE avg_abs_bias mean_min_ESS
#   <chr> <chr>          <dbl> <chr>                  <lgl>             <chr>       <lgl>       <dbl>  <int>           <dbl>    <dbl>        <dbl>        <dbl>
# 1 HADS  Xu_fixed_C        10 Xu                     FALSE             one         TRUE        0.123    576             440     4.62        0.621         6190
# 2 HADS  Xu_fixed_C        10 Xu                     FALSE             informative FALSE       0.172    314             470     4.84        0.582         5924
##
## informative:
##
#  A tibble: 2 × 13
#   test  DGM        n_studies model_parameterisation random_thresholds alpha_pi    mcse_met MCSE_val n_sims min_sims_needed avg_RMSE avg_abs_bias mean_min_ESS
#   <chr> <chr>          <dbl> <chr>                  <lgl>             <chr>       <lgl>       <dbl>  <int>           <dbl>    <dbl>        <dbl>        <dbl>
# 1 HADS  Xu_fixed_C        10 Xu                     FALSE             one         TRUE        0.123    576             248     4.62        0.621         6190
# 2 HADS  Xu_fixed_C        10 Xu                     FALSE             informative TRUE        0.184    287             279     4.68        1.11          6810



##
## ---- HADS; n_studies = 50:  --------------------------------------------------------------------------------------------------------
##
# scenario_df <- filter(scenario_df, model_parameterisation == "Bivariate")
# scenario_df
##
results_HADS_n_studies_50 <- process_all_scenarios_MetaOrd(
  reduction_factor = reduction_factor,
  scenario_df = scenario_df,
  output_dir  = output_dir,
  test_filter = "HADS",
  n_studies_filter = c(50),
  target_MCSE_RMSE = target_MCSE_RMSE,
  max_seed_to_check = max_seed_to_check,
  min_runs_needed = 5,
  min_ESS_thr = 500,
  vec_index_inner = 1:21,
  min_studies_per_threshold = min_studies_per_threshold)
##
results_HADS_n_studies_50_summary <- filter( results_HADS_n_studies_50$summary_table,
                                             n_studies == 50,
                                             alpha_pi == "one")
##
results_HADS_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_HADS_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
##
##  --------------------------------------------------------------------
##
# metric <- "Bias"
metric <- "RMSE"
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_10_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_10_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_10_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_10_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
## ----
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_50_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_50_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_50_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_HADS_n_studies_50_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
## ---- BAI; n_studies = 10: --------------------------------------------------------------------------------------------------------
##
reduction_factor <- 1.0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <- 0.5
##
scenario_df <- create_MetaOrd_scenarios() %>% print(n = 100)
##
## ---- Process all or filtered scenarios:
##
results_BAI_n_studies_10 <- process_all_scenarios_MetaOrd(
                            scenario_df = scenario_df,
                            output_dir = output_dir,
                            test_filter = "BAI",
                            n_studies_filter = 10,
                            target_MCSE_RMSE = target_MCSE_RMSE,
                            max_seed_to_check = max_seed_to_check,
                            min_runs_needed = 3,
                            min_ESS_thr = 1000,
                            vec_index_inner = 5:43,
                            ######## vec_index_inner = 1:63,
                            min_studies_per_threshold = min_studies_per_threshold)
##
results_BAI_n_studies_10_summary <- filter( results_BAI_n_studies_10$summary_table,
                                            n_studies == 10)
##
results_BAI_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi, 
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_10_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_10_summary %>% 
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi, 
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_10_summary %>% 
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi, 
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
## 
## ---- BAI; n_studies = 50: --------------------------------------------------------------------------------------------------------
## 
results_BAI_n_studies_50 <- process_all_scenarios_MetaOrd(
                            scenario_df = scenario_df,
                            output_dir = output_dir,
                            test_filter = "BAI",
                            n_studies_filter = 50,
                            target_MCSE_RMSE = target_MCSE_RMSE,
                            max_seed_to_check = max_seed_to_check,
                            min_runs_needed = 3,
                            min_ESS_thr = 250,
                            vec_index_inner = 5:43,
                            ######## vec_index_inner = 1:63,
                            min_studies_per_threshold = min_studies_per_threshold)
##
results_BAI_n_studies_50_summary <- filter( results_BAI_n_studies_50$summary_table,
                                            n_studies == 50)
##
results_BAI_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Jones") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_fixed_C") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "Xu_random_C") %>%
  print(n = 100, width = 9000)
##
results_BAI_n_studies_50_summary %>%
  select(test, DGM, n_studies, model_parameterisation, random_thresholds, alpha_pi,
         mcse_met, MCSE_val, n_sims, min_sims_needed,
         avg_RMSE,
         avg_abs_bias,
         mean_min_ESS) %>%
  filter(DGM == "R&G_random_C") %>%
  print(n = 100, width = 9000)
##
##  --------------------------------------------------------------------
##
metric <- "Bias"
# metric <- "RMSE"
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_10_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_10_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_10_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_10_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
## ----
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_50_summary, DGM == "Jones"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_50_summary, DGM == "Xu_fixed_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_50_summary, DGM == "Xu_random_C"),
  metric = metric,
  rel_diff_threshold = 10)
##
create_MetaOrd_summary_3(
  data =  filter(results_BAI_n_studies_50_summary, DGM == "R&G_random_C"),
  metric = metric,
  rel_diff_threshold = 10)








# You said to split into two analyses:
#   Analysis 1 (Main - Ordinal models only):
#   
#   Models: Xu, Jones, R&G (NO bivariate)
# Thresholds: ALL thresholds (including those with 0, 1, 2 studies)
# Purpose: Shows ordinal models' full capabilities (borrowing + interpolation)
# 
# Analysis 2 (Secondary - Fair comparison):
# 
# Models: Bivariate, Xu, Jones, R&G (ALL 4 models)
# Thresholds: Only ≥3 studies per threshold
# Purpose: Fair comparison when bivariate has adequate data
# 
# Rationale:
# 
# Bivariate CANNOT estimate thresholds with <3 studies (needs data at each threshold independently)
# Ordinal models CAN because they borrow information across thresholds
# Including bivariate on sparse thresholds would be unfair
# The main scientific point is showing ordinal models work on ALL thresholds, while bivariate needs restriction
# 


mat <- readRDS("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/real_x_GAD2.RDS")
mat <- readRDS("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/real_x_HADS.RDS")
mat <- readRDS("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/real_x_BAI.RDS")
##
nrow(mat$diseased) 
##
## GAD_2 = 29 studies, HADS = 36 studies, BAI = 7 studies.



# rename_files_add_red_fact(directory = MetaOrdDTA::metaorddta_output_path("sim_study", "simulation_results"),
#                           reduction_factor = 1.00,
#                           dry_run = F)


{
      ##
      source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
      source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
      source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
      source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
      source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
      source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
      ##
      # data_0 <- sim_data_GAD_2_summary_0
      ##
      sim_data_GAD_2_summary_0 <- rbind( results_GAD_2_n_studies_10_summary,
                                         results_GAD_2_n_studies_50_summary)
      # sim_data_GAD_2_summary <- sim_data_GAD_2_summary_0 %>% mutate( model_parameterisation = ifelse( model_parameterisation %in% c("Xu"),
      #                                                                                                 "Ord-Biv",
      #                                                                                                 model_parameterisation ))
      # 
      # data <- data_0 %>% 
      #   dplyr::mutate( model_parameterisation = case_when( model_parameterisation == "R&G" ~ "Ord-HSROC",
      #                                                      model_parameterisation == "Xu"  ~ "Ord-Biv",
      #                                                      TRUE ~ model_parameterisation ))
      # data
      # 
      # 
      # ##
      sim_data_HADS_summary_0 <- rbind( results_HADS_n_studies_10_summary,
                                        results_HADS_n_studies_50_summary)
      # sim_data_HADS_summary <- sim_data_HADS_summary_0 %>% mutate( model_parameterisation = ifelse( model_parameterisation %in% c("Xu"),
      #                                                                                               "Ord-Biv",
      #                                                                                               model_parameterisation ))
      ##
      sim_data_BAI_summary_0 <- rbind( results_BAI_n_studies_10_summary,
                                       results_BAI_n_studies_50_summary)
      # sim_data_BAI_summary <- sim_data_BAI_summary_0 %>% mutate( model_parameterisation = ifelse( model_parameterisation %in% c("Xu"),
      #                                                                                             "Ord-Biv",
      #                                                                                             model_parameterisation ))
      ##
      sim_data_GAD_2_summary   <- R_fn_edit_sim_data_tibble(sim_data_GAD_2_summary_0)
      sim_data_HADS_summary    <- R_fn_edit_sim_data_tibble(sim_data_HADS_summary_0)
      sim_data_BAI_summary     <- R_fn_edit_sim_data_tibble(sim_data_BAI_summary_0)
      ##
      sim_data_GAD_2_summary
}
##
sim_data_GAD_2_summary_0$model_parameterisation
sim_data_GAD_2_summary$model_parameterisation
##
sim_data_GAD_2_summary_0
sim_data_GAD_2_summary
##
sim_data_GAD_2_summary$model_parameterisation
sim_data_GAD_2_summary$Model
## 
## ---------------------------------------------------------------------------------------------------
##
require(patchwork)
##
## ---- Create plots for each DGM:
##
# Get unique DGMs
# DGM_names <- unique(overall_sim_data$DGM)
##
DGM_names <- c("Jones",
               "Xu_fixed_C",
               "Xu_random_C",
               "R&G_random_C")
##
test_vec <- c("GAD-2", 
              "HADS",
              "BAI")
##
# Create lists to store plots
rmse_plots_Se <- rmse_plots_Sp <- list()
bias_plots_Se <- bias_plots_Sp <- list()
coverage_plots_Se <- coverage_plots_Sp <- list()
width_plots_Se <- width_plots_Sp <- list()
##
base_size <- 14
##
# dgm <- "Jones"
# dgm <- "Xu_fixed_C"
# dgm <- "Xu_random_C"
dgm <- "R&G_random_C"
##
test <- "GAD-2"
# test <- "HADS"
# test <- "BAI"
##
for (test in test_vec) {
  
            if (test == "GAD-2") { 
              overall_sim_data <- sim_data_GAD_2_summary
            } else if (test == "HADS") { 
              overall_sim_data <- sim_data_HADS_summary
            } else if (test == "BAI") { 
              overall_sim_data <- sim_data_BAI_summary
            }
            ##
            ## ---- RMSE plots:
            ##
            rmse_plots_Se[[test]] <- create_RMSE_plot( data = overall_sim_data,
                                                       metric_type = "Se",
                                                       base_size = base_size)
            rmse_plots_Se[[test]]
            ##
            rmse_plots_Sp[[test]] <- create_RMSE_plot( data = overall_sim_data,
                                                       metric_type = "Sp",
                                                       base_size = base_size)
            rmse_plots_Sp[[test]]
            ##
            ## ---- Bias plots:
            ##
            bias_plots_Se[[test]] <- create_bias_plot( data = overall_sim_data,
                                                       metric_type = "Se",
                                                       base_size = base_size)
            bias_plots_Se[[test]]
            ##
            bias_plots_Sp[[test]] <- create_bias_plot( data = overall_sim_data,
                                                       metric_type = "Sp",
                                                       base_size = base_size)
            bias_plots_Sp[[test]]
            ##
            ## ---- Coverage plots:
            ##
            coverage_plots_Se[[test]] <- create_coverage_plot( data = overall_sim_data,
                                                               metric_type = "Se",
                                                               base_size = base_size)
            coverage_plots_Se[[test]]
            ##
            coverage_plots_Sp[[test]] <- create_coverage_plot( data = overall_sim_data,
                                                               metric_type = "Sp",
                                                               base_size = base_size)
            coverage_plots_Sp[[test]]
            ##
            ## ---- Width plots:
            ##
            width_plots_Se[[test]] <- create_width_plot(   data = overall_sim_data,
                                                           metric_type = "Se",
                                                           base_size = base_size)
            width_plots_Se[[test]] 
            ##
            width_plots_Sp[[test]] <- create_width_plot(  data = overall_sim_data,
                                                          metric_type = "Sp",
                                                          base_size = base_size)
            width_plots_Sp[[test]]
            
  
}
##
## ----
##
overall_sim_data %>% print(n = 1000, width = 9000)
##
sim_data_GAD_2_summary  %>% 
  # select(-avg_RMSE_Se, -avg_RMSE_Sp, -avg_RMSE, -max_RMSE_Se, -max_RMSE_Sp, 
  #        -avg_MCSE_RMSE_Se, -avg_MCSE_RMSE_Sp, -max_MCSE_RMSE_Se, -max_MCSE_RMSE_Sp, 
  #        -avg_abs_bias_Se, -avg_abs_bias_Sp) %>%
  print(n = 1000, width = 9000)
##
sim_data_HADS_summary  %>% 
  # select(-avg_RMSE_Se, -avg_RMSE_Sp, -avg_RMSE, -max_RMSE_Se, -max_RMSE_Sp, 
  #        -avg_MCSE_RMSE_Se, -avg_MCSE_RMSE_Sp, -max_MCSE_RMSE_Se, -max_MCSE_RMSE_Sp, 
  #        -avg_abs_bias_Se, -avg_abs_bias_Sp) %>%
  print(n = 1000, width = 9000)
##
sim_data_BAI_summary   %>% 
  # select(-avg_RMSE_Se, -avg_RMSE_Sp, -avg_RMSE, -max_RMSE_Se, -max_RMSE_Sp, 
  #        -avg_MCSE_RMSE_Se, -avg_MCSE_RMSE_Sp, -max_MCSE_RMSE_Se, -max_MCSE_RMSE_Sp, 
  #        -avg_abs_bias_Se, -avg_abs_bias_Sp) %>%
  print(n = 1000, width = 9000)
##
sim_data_GAD_2_summary %>% filter(n_studies == 10) %>% pull(avg_RMSE) %>% mean() ## 8.916097
sim_data_GAD_2_summary %>% filter(n_studies == 50) %>% pull(avg_RMSE) %>% mean() ## 4.206083
##
sim_data_HADS_summary %>% filter(n_studies == 10) %>% pull(avg_RMSE) %>% mean() ## 5.034427
sim_data_HADS_summary %>% filter(n_studies == 50) %>% pull(avg_RMSE) %>% mean() ## 2.31455
##
sim_data_BAI_summary %>% filter(n_studies == 10) %>% pull(avg_RMSE) %>% mean() ## 6.504302
sim_data_BAI_summary %>% filter(n_studies == 50) %>% pull(avg_RMSE) %>% mean() ## 3.232993
##
## ----
##
plot_RMSE_GAD_2 <-  rmse_plots_Se[["GAD-2"]] +
                    rmse_plots_Sp[["GAD-2"]] +
                    plot_layout(ncol = 1)
##
plot_RMSE_HADS <-   rmse_plots_Se[["HADS"]] +
                    rmse_plots_Sp[["HADS"]] +
                    plot_layout(ncol = 1)
##
plot_RMSE_BAI <-    rmse_plots_Se[["BAI"]] +
                    rmse_plots_Sp[["BAI"]] +
                    plot_layout(ncol = 1)
##
## ----
##
plot_Bias_GAD_2 <-  bias_plots_Se[["GAD-2"]] +
                    bias_plots_Sp[["GAD-2"]] +
                    plot_layout(ncol = 1)
##
plot_Bias_HADS <-   bias_plots_Se[["HADS"]] +
                    bias_plots_Sp[["HADS"]] +
                    plot_layout(ncol = 1)
##
plot_Bias_BAI <-    bias_plots_Se[["BAI"]] +
                    bias_plots_Sp[["BAI"]] +
                    plot_layout(ncol = 1)
##
## ----
##
plot_Cvg_GAD_2 <-   coverage_plots_Se[["GAD-2"]] +
                    coverage_plots_Sp[["GAD-2"]] +
                    plot_layout(ncol = 1)
##
plot_Cvg_HADS <-    coverage_plots_Se[["HADS"]] +
                    coverage_plots_Sp[["HADS"]] +
                    plot_layout(ncol = 1)
##
plot_Cvg_BAI <-     coverage_plots_Se[["BAI"]] +
                    coverage_plots_Sp[["BAI"]] +
                    plot_layout(ncol = 1)
##
dpi <- 200
mult <- 1.1
##
plot_RMSE_GAD_2
plot_Bias_GAD_2
plot_Cvg_GAD_2
##
ggsave("Sim_study_RMSE_GAD_2.png", plot_RMSE_GAD_2, width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_RMSE_HADS.png",  plot_RMSE_HADS,  width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_RMSE_BAI.png",   plot_RMSE_BAI,   width = 11*mult, height = 16*mult, dpi = dpi)
##
ggsave("Sim_study_Bias_GAD_2.png", plot_Bias_GAD_2, width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_Bias_HADS.png",  plot_Bias_HADS,  width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_Bias_BAI.png",   plot_Bias_BAI,   width = 11*mult, height = 16*mult, dpi = dpi)
##
ggsave("Sim_study_Coverage_GAD_2.png", plot_Cvg_GAD_2, width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_Coverage_HADS.png",  plot_Cvg_HADS,  width = 11*mult, height = 16*mult, dpi = dpi)
ggsave("Sim_study_Coverage_BAI.png",   plot_Cvg_BAI,   width = 11*mult, height = 16*mult, dpi = dpi)
##
## ----
##
# ggsave("Sim_study_Bias_GAD_2_DGM_1_Jones.png", bias_plots[["GAD-2"]][["Jones"]], width = 9, height = 16, dpi = dpi)
# ggsave("Sim_study_Bias_GAD_2_DGM_2_Xu_fixed.png", bias_plots[["GAD-2"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_GAD_2_DGM_3_Xu_random.png", bias_plots[["GAD-2"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_GAD_2_DGM_4_R&G_random.png", bias_plots[["GAD-2"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# ##
# ggsave("Sim_study_Bias_HADS_DGM_1_Jones.png", bias_plots[["HADS"]][["Jones"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_HADS_DGM_2_Xu_fixed.png", bias_plots[["HADS"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_HADS_DGM_3_Xu_random.png", bias_plots[["HADS"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_HADS_DGM_4_R&G_random.png", bias_plots[["HADS"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# ##
# ggsave("Sim_study_Bias_BAI_DGM_1_Jones.png", bias_plots[["BAI"]][["Jones"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_BAI_DGM_2_Xu_fixed.png", bias_plots[["BAI"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_BAI_DGM_3_Xu_random.png", bias_plots[["BAI"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Bias_BAI_DGM_4_R&G_random.png", bias_plots[["BAI"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# ##
# ## ----
# ##
# ggsave("Sim_study_Coverage_GAD_2_DGM_1_Jones.png", coverage_plots[["GAD-2"]][["Jones"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_GAD_2_DGM_2_Xu_fixed.png", coverage_plots[["GAD-2"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_GAD_2_DGM_3_Xu_random.png", coverage_plots[["GAD-2"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_GAD_2_DGM_4_R&G_random.png", coverage_plots[["GAD-2"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# ##
# ggsave("Sim_study_Coverage_HADS_DGM_1_Jones.png", coverage_plots[["HADS"]][["Jones"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_HADS_DGM_2_Xu_fixed.png", coverage_plots[["HADS"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_HADS_DGM_3_Xu_random.png", coverage_plots[["HADS"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_HADS_DGM_4_R&G_random.png", coverage_plots[["HADS"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# ##
# ggsave("Sim_study_Coverage_BAI_DGM_1_Jones.png", coverage_plots[["BAI"]][["Jones"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_BAI_DGM_2_Xu_fixed.png", coverage_plots[["BAI"]][["Xu_fixed_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_BAI_DGM_3_Xu_random.png", coverage_plots[["BAI"]][["Xu_random_C"]], width = 21, height = 9, dpi = dpi)
# ggsave("Sim_study_Coverage_BAI_DGM_4_R&G_random.png", coverage_plots[["BAI"]][["R&G_random_C"]], width = 21, height = 9, dpi = dpi)
# 

 
 



ggsave("Sim_study_Bias_GAD_2.png", bias_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_GAD_2.png", coverage_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)

##
## ---- Plots for HADS:
##
# sROC_plots[["HADS"]]
rmse_plots[["HADS"]]
bias_plots[["HADS"]]
coverage_plots[["HADS"]]
##
ggsave("Sim_study_RMSE_HADS.png", rmse_plots[["HADS"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Bias_HADS.png", bias_plots[["HADS"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_HADS.png", coverage_plots[["HADS"]], width = 21, height = 9, dpi = 400)

##
## ---- Plots for BAI:j
##
# sROC_plots[["BAI"]]
rmse_plots[["BAI"]]
bias_plots[["BAI"]]
coverage_plots[["BAI"]]
##
# ggsave("Sim_study_sROC_BAI.png", sROC_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_RMSE_BAI.png", rmse_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Bias_BAI.png", bias_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_BAI.png", coverage_plots[["BAI"]], width = 21, height = 9, dpi = 400)




graphics.off()   # Close all devices
.rs.restartR()   # Restart R session in RStudio

dev.off()  # Close current device

X11()  # For Linux


library(ggplot2)
library(patchwork)
# ## --- Plot 4: Coverage probability -------------------------------------------------------------------
# ##
# plot_coverage_Se <- overall_sim_data %>%
#   pivot_longer(cols = c(avg_coverage_Se, avg_coverage_Sp),
#                names_to = "metric",
#                values_to = "coverage") %>%
#   mutate(
#     metric = factor(metric,
#                     levels = c("avg_coverage_Se"),
#                     labels = c("Se"))
#   ) %>%
#   ggplot(aes(x = Model, y = coverage, color = Model)) +
#   geom_point(size = 4) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      limits = c(0.7, 1)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(y = "Coverage Probability",
#        title = "95% CI Coverage Across Models",
#        subtitle = "Red line = Nominal 95% coverage") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   facet_wrap(~ DGM + n_studies, scales = "free_y", nrow = 3)
# ##
# plot_coverage_Se
# ##
# plot_coverage_Sp <- overall_sim_data %>%
#   pivot_longer(cols = c(avg_coverage_Se, avg_coverage_Sp),
#                names_to = "metric",
#                values_to = "coverage") %>%
#   mutate(
#     metric = factor(metric,
#                     levels = c("avg_coverage_Sp"),
#                     labels = c("Sp"))
#   ) %>%
#   ggplot(aes(x = Model, y = coverage, color = Model)) +
#   geom_point(size = 4) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      limits = c(0.7, 1)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(y = "Coverage Probability",
#        title = "95% CI Coverage Across Models",
#        subtitle = "Red line = Nominal 95% coverage") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   facet_wrap(~ DGM + n_studies, scales = "free_y", nrow = 3)
# ##
# plot_coverage_Sp
# ##
# ##
# require(patchwork)
# plot_coverage_Se + plot_coverage_Sp
# 
# 
# 
# true_mean_GAD_2_params_list$Ord_bivariate$beta_SD <- c(0.341, 0.559)
# true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.172
# ##
# true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.235,  0.392)
# true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.302
# ##
# true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.257,  0.368)
# true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.005574
# 

# ##
# ## --- Combine plots
# ##
# # Using patchwork to combine
# combined_plot <- (plot_rmse | plot_bias) / 
#   plot_coverage /
#   plot_accuracy +
#   plot_layout(heights = c(1, 1, 2))
# 
# # Save
# ggsave("simulation_results_combined.png", combined_plot, 
#        width = 12, height = 14, dpi = 300)
# 
































