





# This is actually a valid research finding, not a problem to "fix"!
#   What's happening:
#
# Random-cutpoint models estimate kappa → ∞ when heterogeneity → 0
# Large kappa causes numerical instability in Stan (kappa × simplex overflows)
# This is a real limitation of random-cutpoint models, not a bug
#
# Options:
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

# 
# 
# system2("/home/enzocerullo/R/R-4.3.3/bin/R",
#         args = c("CMD",
#                  "INSTALL",
#                  "--library=/home/enzocerullo/R/R-4.3.3/lib/R/library",
#                  "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA"))
##
# rm(list = ls())
# .rs.restartR()  # In RStudio
##
# ============================================================================================================
# Main Execution Script for MetaOrdDTA Simulation Study
# ============================================================================================================

cat("\014") ## clear console

# # # # # # # # # ##
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
#         ## First remove any possible packagDe fragments:
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
##
## ---- The frozen NicoStan and NicoSIMD of this sim: installs by other work never reach its workers ----------------------------------
##
## (copies of the installs of 24 Oct, 21:37 and 21:42; R_LIBS puts them first in every worker and every
##  fresh-process fit)
frozen_library_of_this_sim <- "/home/enzocerullo/R/frozen_library_of_the_MetaOrdDTA_sim"
##
## ---- Unload any MetaOrdDTA / NicoStan / NicoSIMD loaded from the main library (e.g. by the restored .RData) ------------------------
##
for (package_to_reload in c("MetaOrdDTA", "NicoStan", "NicoSIMD")) {

      if (package_to_reload %in% loadedNamespaces() &&
          !startsWith(getNamespaceInfo(package_to_reload, "path"), frozen_library_of_this_sim)) {
            unloadNamespace(package_to_reload)
      }

}
.libPaths(c(frozen_library_of_this_sim, .libPaths()))
Sys.setenv(R_LIBS = frozen_library_of_this_sim)
library(MetaOrdDTA)
library(doParallel)
stopifnot(startsWith(find.package("NicoStan"), frozen_library_of_this_sim),
          startsWith(find.package("NicoSIMD"), frozen_library_of_this_sim))
# library(dplyr)
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

# # Main execution
# scenario_row <- data.frame(
#   index_test_MA = "GAD_2",
#   n_studies = 50,
#   model_parameterisation = "ord_bivariate",
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
##
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
##
n_threads_to_use_HPC <- 160
# n_threads_to_use_HPC <- 100
##
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
#                                 pattern = ".*GAD_2.*red_fact_1.*ord_biv_random_C.*n_studies_10.*rand_thr_FALSE.*",
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
# scenario_df <- scenario_df  %>% filter(DGM ==("ord_biv_fixed_C"))
## (the ord_bivariate fixed-cutpoint AND the ord_bivariate random-kappa models; the random one runs on the ord_bivariate random-cutpoint data
##  only, by create_MetaOrd_scenarios(); the ordinal-HSROC random-cutpoint model is added once its informative
##  priors from the literature exist)
scenario_df %>% print(n = 1000)
##
# A tibble: 5 × 13
# test  DGM         n_studies model_parameterisation random_thresholds alpha_pi mcse_met max_MCSE n_sims min_sims_needed avg_RMSE avg_abs_bias mean_min_ESS
# <chr> <chr>           <dbl> <chr>                  <lgl>             <chr>    <lgl>       <dbl>  <int>           <dbl>    <dbl>        <dbl>        <dbl>
# 3 HADS  ord_biv_random_C        50 Bivariate              FALSE             one      FALSE       0.126    114           132.      2.80        0.391          917
##
## ---- Test:
##
test_vals_vec <- c("GAD_2", "HADS", "BAI") ## ----
##
# test_vals_vec <- "GAD_2"
# test_vals_vec <- "HADS"
# test_vals_vec <- "BAI"
# test_vals_vec <- c("HADS", "BAI")
##
## ---- Select algorithm:
##
sampling_algorithm <- "NicoStan"
# sampling_algorithm <- "Stan"
## 
## ---- Sim study type ----------------------------------------------------------------------------------------------------------------
##
## ("paper" = the sim study of the paper; "RC_ord_biv_median_vs_pop_median" = the tests and numbers of studies
##  set below, two random-cutpoint ord_biv models: the median model of the paper,
##  DTA_MA_ord_biv_RANDthr_kappa_NicoStan_fused_AVX, and the pop-median model,
##  DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX, its Se / Sp prior off; each with its
##  off-centred cutpoints (C_raw) and with its centred cutpoints (C_v):
##  DTA_MA_ord_biv_RANDthr_kappa_centred_C_NicoStan_fused_AVX and
##  DTA_MA_ord_biv_RANDthr_kappa_popmed_centred_C_Se_Sp_prior_NicoStan_fused_AVX;
##  "RC_ord_HSROC_median_vs_pop_median" = the same for the two random-cutpoint ord_HSROC models: the median model,
##  DTA_MA_ord_HSROC_RANDthr_kappa_NicoStan_fused_AVX, and the pop-median model,
##  DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX, its Se / Sp prior off; off-centred
##  cutpoints only)
##
# sim_study_type <- "paper" ## ----
sim_study_type <- "RC_ord_biv_median_vs_pop_median"
# sim_study_type <- "RC_ord_HSROC_median_vs_pop_median"
##
## ---- The prior of the ordinal models (ord_bivariate, R&G): their own Dirichlet prior on the cutpoints or the Se / Sp prior --------------------
##
## ("Dirichlet" = priors$use_Se_Sp_prior <- 0, each model's own Dirichlet prior on the cutpoints; "Se_Sp" =
##  priors$use_Se_Sp_prior <- 1, the direct Se / Sp prior, type 4: N(0, 1) on probit Sp and probit Se at every
##  cutoff (run_MULTIPLE_MetaOrd_simulation.R); the Se / Sp fits are saved as "..._Se_Sp_prior")
ordinal_prior_type <- "Dirichlet"
# ordinal_prior_type <- "Se_Sp"
stopifnot(ordinal_prior_type %in% c("Dirichlet", "Se_Sp"))
##
# Set paths
sim_study_dir <- file.path(local_pkg_dir, "sim_study")
##
{
  
    if (sim_study_type == "paper") { 
      
            ## NicoStan's folder of the sim study: "simulation_results_NicoStan_fused_AVX"
            ## Stan's folder of the divergence test: "simulation_results_Stan_divergence_test"
            if (sampling_algorithm == "NicoStan") { 
              # inner_path <-  "fixed_raw_C_J_simu_results_NicoStan_diag_burnin_250"
              # inner_path <-   "fixed_raw_C_J_sim_results_NicoStan_test_MCSE_0.5_test_2"
              # inner_path <-   "fixed_raw_C_J_sim_results_NicoStan_MCSE_0.25_Se_Sp_prior_flat_1"
                inner_path <-   "fixed_raw_C_J_sim_results_NicoStan_MCSE_0.25_alpha_prior_flat_1"
            } else { 
              inner_path <-  "fixed_raw_C_J_sim_results_Stan_reference_grid"
            }
      
    } else if (sim_study_type == "RC_ord_biv_median_vs_pop_median") {
      
            inner_path <- paste0("RC_ord_biv_med_vs_popmed_", inner_path)
          
    } else if (sim_study_type == "RC_ord_HSROC_median_vs_pop_median") {
      
            inner_path <- paste0("RC_ord_HSROC_med_vs_popmed_", inner_path)
            
    }
  
}
##
## ---- The folder of the sim study type ----------------------------------------------------------------------------------------------
##
## ("paper" keeps the folder above; every other sim study type gets its own folder: the same name + its type)
if (sim_study_type != "paper") {
      inner_path <- paste0(inner_path, "_", sim_study_type)
}
##
output_dir <- MetaOrdDTA::metaorddta_output_path(
  "sim_study",
  inner_path,
  create = TRUE
)
##
# test_vals_vec <- "HADS"
## (all 6 cells: HADS, GAD-2 and BAI at 10 and 50 studies)o
# test_vals_vec <- c("HADS", "GAD_2", "BAI")
# test_vals_vec <- c("HADS", "GAD_2")
# test_vals_vec <- "BAI"
##
## ---- DGM:
##
# DGM_vec <- c("ord_biv_fixed_C", "Jones", "ord_biv_random_C", "R&G_random_C") ## ----
DGM_vec <- c("ord_biv_random_C", "R&G_random_C")
# DGM_vec <- "Jones"
# DGM_vec <- "ord_biv_random_C"
# DGM_vec <- "R&G_random_C"
##
n_studies_vals_vec <- c(10, 50) ## ----
##
## ---- model_parameterisation:
##
# model_param_vec <- c("Jones", "ord_bivariate", "R&G", "Bivariate")
# model_param_vec <- c("ord_bivariate", "Jones", "R&G")
# model_param_vec <- "Bivariate"
# model_param_vec <- c("Jones", "Bivariate")
# model_param_vec <- "Jones"
# model_param_vec <- "R&G"
# model_param_vec <- c("ord_bivariate", "R&G")
# model_param_vec <- "ord_bivariate"
# model_param_vec <- c("ord_bivariate", "R&G", "Bivariate")
# model_param_vec <- c("ord_bivariate", "Bivariate")
model_param_vec <- "ord_bivariate"
##
## ----
##
alpha_pi_vec <- "one"
# alpha_pi_vec <- "informative_from_the_literature"
##
## ---- Settings of the RC ord_biv median vs pop-median sim study ---------------------------------------------------------------------
##
random_thresholds_vec <- c(FALSE, TRUE)
##
if (sim_study_type == "RC_ord_biv_median_vs_pop_median") {
      model_param_vec       <- "ord_bivariate"
      random_thresholds_vec <- TRUE
      RC_ord_biv_model_vec  <- c("median", "pop_median")
      # RC_cutpoints_vec      <- c("off_centred", "centred")
      ## (off-centred only: the centred C_v fits keep none, every one 8-12% divergent with min ESS 4)
      RC_cutpoints_vec      <- "off_centred"
}
##
## ---- Settings of the RC ord_HSROC median vs pop-median sim study -------------------------------------------------------------------
##
if (sim_study_type == "RC_ord_HSROC_median_vs_pop_median") {
      model_param_vec          <- "R&G"
      random_thresholds_vec    <- TRUE
      RC_ord_HSROC_model_vec   <- c("median", "pop_median")
}
##
## ---- sampler and the min ESS of Se / Sp a fit needs to be kept:
##
min_ESS_of_10_studies <- 1000
min_ESS_of_50_studies <- 200
##
# test_val    <- "BAI"
# model_param <- model_param_vec
##
## ----
##
# target_MCSE_RMSE <- 1.0
# target_MCSE_RMSE <- 0.50
target_MCSE_RMSE <- 0.25
# target_MCSE_RMSE <- 0.125 ## ----
##
## ----
##
max_seed_to_check <- 40000
min_runs_needed <- 30
##
## ----
##
starting_sim_offset <- 1
## (the laptop's share: 180 chains in all on the local-HPC (4 per fit) and 4 fits on the laptop (all its 16
##  threads), the two never overlapping: the laptop's own seeds from 20001 (below max_seed_to_check, which
##  the MCSE checks read up to), its 16 threads, and the
##  scenarios in the reverse order of the local-HPC's (BAI first), so that the two machines meet in the middle; the
##  results of both are copied to each other, so that each machine's MCSE checks count them all)
machine_of_this_run <- if (any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))) "local_HPC" else
                         "laptop"
if (machine_of_this_run == "laptop") {
  starting_sim_offset <- 20000
  test_vals_vec       <- rev(test_vals_vec)
  n_studies_vals_vec  <- rev(n_studies_vals_vec)
  DGM_vec             <- rev(DGM_vec)
}
# n_sims_total <- starting_sim_offset + 30 ; n_sims_total
n_sims_total <- starting_sim_offset + 1000 ; n_sims_total
##
## ----
##
root_MCMC_seed <- 123
##
reduction_factor <- 1.0
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
## ---- the fixed-cutpoint cells first, then the random-cutpoint cells:
for (random_thresholds_of_this_pass in random_thresholds_vec) {
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
            alpha_pi %in% alpha_pi_vec,
            random_thresholds == random_thresholds_of_this_pass
          ) %>%
          print(n = 100)

        if (nrow(run_matrix_filtered) == 0) next
        ##
        ## ---- Each scenario once per model (median, pop-median) and cutpoints (off-centred, centred) --------------------------------
        ##
        if (sim_study_type == "RC_ord_biv_median_vs_pop_median") {
              run_matrix_filtered <- tidyr::expand_grid(run_matrix_filtered,
                                                        RC_ord_biv_model = RC_ord_biv_model_vec,
                                                        RC_cutpoints     = RC_cutpoints_vec)
        }
        ##
        ## ---- Each scenario once per ord_HSROC model (median, pop-median) -----------------------------------------------------------
        ##
        if (sim_study_type == "RC_ord_HSROC_median_vs_pop_median") {
              run_matrix_filtered <- tidyr::expand_grid(run_matrix_filtered,
                                                        RC_ord_HSROC_model = RC_ord_HSROC_model_vec)
        }

        ## Simulate data + get true values
        source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA_v2.R"), local = TRUE)
        ## true_DGM_Se and true_DGM_Sp will be set by the simulation

        scenario_idx <- 1

        ## Process each scenario in the filtered set
        for (scenario_idx in 1:nrow(run_matrix_filtered)) {

                    scenario_row <- run_matrix_filtered[scenario_idx, ] ; scenario_row
                    min_ESS_of_this_cell <-  if (scenario_row$n_studies == 10) min_ESS_of_10_studies else
                                             min_ESS_of_50_studies
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
                    ## ---- The pop-median model, its user header and its Se / Sp prior fields, the prior off -------------------------
                    ##
                    ## (the fields of R_fn_set_Se_Sp_prior_MA.R:44-65 with n_chosen_max = 1, nothing chosen and
                    ##  use_Se_Sp_prior = 0; the worker takes them at run_MULTIPLE_MetaOrd_simulation.R:62-81,
                    ##  142-145, 157-159 and 540-542)
                    if (sim_study_type == "RC_ord_biv_median_vs_pop_median" &&
                        scenario_row$RC_ord_biv_model == "pop_median") {
                          scenario_row$custom_file_name <- list(
                                "DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan")
                          # scenario_row$custom_cpp_user_header_file_path <- list(
                          #       file.path(local_pkg_dir, "inst", "stan_functions", "RC_popmed_user_header.hpp"))
                          ## (the same two includes, from the frozen NicoSIMD)
                          scenario_row$custom_cpp_user_header_file_path <- list(
                                file.path(frozen_library_of_this_sim, "RC_popmed_user_header.hpp"))
                          scenario_row$prior_fields <- list(list(
                                use_Se_Sp_prior          = 0,
                                Se_Sp_prior_type         = 4,
                                n_chosen_max             = 1,
                                n_chosen                 = c(0L, 0L),
                                chosen_cutoffs           = matrix(1L, 2, 1),
                                prior_Se_Sp_dirichlet_A  = list(rep(0, 2), rep(0, 2)),
                                prior_Se_Sp_beta         = list(matrix(0, 1, 2), matrix(0, 1, 2)),
                                prior_Se_Sp_link_mean_SD = list(matrix(0, 1, 2), matrix(0, 1, 2))))
                          file_name_string <- paste0(file_name_string, "_popmed")
                    }
                    ##
                    ## ---- The ord_HSROC pop-median model, its user header and its Se / Sp prior fields, the prior off ---------------
                    ##
                    ## (the same fields as the ord_biv pop-median model above)
                    if (sim_study_type == "RC_ord_HSROC_median_vs_pop_median" &&
                        scenario_row$RC_ord_HSROC_model == "pop_median") {
                          scenario_row$custom_file_name <- list(
                                "DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan")
                          scenario_row$custom_cpp_user_header_file_path <- list(
                                file.path(frozen_library_of_this_sim, "RC_popmed_user_header.hpp"))
                          scenario_row$prior_fields <- list(list(
                                use_Se_Sp_prior          = 0,
                                Se_Sp_prior_type         = 4,
                                n_chosen_max             = 1,
                                n_chosen                 = c(0L, 0L),
                                chosen_cutoffs           = matrix(1L, 2, 1),
                                prior_Se_Sp_dirichlet_A  = list(rep(0, 2), rep(0, 2)),
                                prior_Se_Sp_beta         = list(matrix(0, 1, 2), matrix(0, 1, 2)),
                                prior_Se_Sp_link_mean_SD = list(matrix(0, 1, 2), matrix(0, 1, 2))))
                          file_name_string <- paste0(file_name_string, "_popmed")
                    }
                    ##
                    ## ---- The centred cutpoints: the centred file of each model (the median model on the NicoSIMD header) -----------
                    ##
                    if (sim_study_type == "RC_ord_biv_median_vs_pop_median" &&
                        scenario_row$RC_cutpoints == "centred") {
                          if (scenario_row$RC_ord_biv_model == "pop_median") {
                                scenario_row$custom_file_name <- list(paste0(
                                      "DTA_MA_ord_biv_RANDthr_kappa_popmed_centred_C_Se_Sp_prior_",
                                      "NicoStan_fused_AVX.stan"))
                          } else {
                                scenario_row$custom_file_name <- list(
                                      "DTA_MA_ord_biv_RANDthr_kappa_centred_C_NicoStan_fused_AVX.stan")
                                # scenario_row$custom_cpp_user_header_file_path <- list(
                                #     file.path(dirname(dirname(local_pkg_dir)), "R_packages", "NicoStan", "inst",
                                #               "NicoSIMD", "inst", "include", "NicoSIMD",
                                #               "stan_external_functions_of_NicoSIMD.hpp"))
                                scenario_row$custom_cpp_user_header_file_path <- list(
                                      file.path(find.package("NicoSIMD"), "include", "NicoSIMD",
                                                "stan_external_functions_of_NicoSIMD.hpp"))
                          }
                          file_name_string <- paste0(file_name_string, "_centred_C")
                    }
                    ##
                    ## ---- The Se / Sp prior: its fits saved under names of their own ------------------------------------------------
                    ##
                    ## (the same test as the worker's: run_MULTIPLE_MetaOrd_simulation.R, "The direct Se/Sp prior")
                    uses_the_Se_Sp_prior <- ordinal_prior_type == "Se_Sp" &&
                                            scenario_row$model_parameterisation %in% c("ord_bivariate", "R&G") &&
                                            scenario_row$alpha_pi == "one"
                    if (uses_the_Se_Sp_prior &&
                        identical(scenario_row$custom_file_name[[1]],
                                  "DTA_MA_ord_biv_RANDthr_kappa_centred_C_NicoStan_fused_AVX.stan")) {
                          stop("The centred median model has no Se / Sp prior file: set ordinal_prior_type ",
                               "to \"Dirichlet\" or leave the centred median model out.")
                    }
                    if (uses_the_Se_Sp_prior) {
                          file_name_string <- paste0(file_name_string, "_Se_Sp_prior")
                    }
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
                                                # min_ESS_thr = if (scenario_row$n_studies == 10) 1000 else 200,
                                                  min_ESS_thr = min_ESS_of_this_cell,
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
                    # n_threads_to_use <- 180
                    n_threads_to_use <- if (machine_of_this_run == "laptop") 16 else n_threads_to_use_HPC
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
                    minimal_env$sampling_algorithm <- sampling_algorithm
                    minimal_env$min_ESS_of_10_studies <- min_ESS_of_10_studies
                    minimal_env$min_ESS_of_50_studies <- min_ESS_of_50_studies
                    minimal_env$ordinal_prior_type <- ordinal_prior_type

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
                                                                            # min_ESS_thr = if (scenario_row$n_studies == 10) 1000 else 200,
                                                                            min_ESS_thr = min_ESS_of_this_cell,
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

                            ## ---- Hand the seeds to the workers one at a time:
                            ## (foreach gives each seed to whichever worker is free, through doParallel's
                            ##  clusterApplyLB, so no worker waits for another at the end of a cell)
                            minimal_env$sims_this_batch <- sims_this_batch
                            ##
                            ## (the file a worker writes once this cell's MCSE target is met: every seed not yet
                            ##  started then returns at once instead of being fitted)
                            file_of_the_MCSE_target_met <- file.path(tempdir(),
                                                                     paste0("MCSE_target_met_",
                                                                            file_name_string,
                                                                            "_batch_",
                                                                            batch_id,
                                                                            "_",
                                                                            format(Sys.time(), "%Y%m%d_%H%M%S"),
                                                                            ".flag"))
                            minimal_env$file_of_the_MCSE_target_met <- file_of_the_MCSE_target_met

                            ##
                            ## ---- Build this scenario's model ONCE, before the workers start: ---------------------------------------
                            ## (a run that started 45 workers, each of which
                            ##  built the same cmdstanr executable at the same moment: one worker read the
                            ##  executable whilst another was writing it, which raised "Permission denied",
                            ##  "No such file or directory" and a missing file time. The model is built here once,
                            ##  with the same arguments MetaOrd_model$new gives it
                            ##  (run_MULTIPLE_MetaOrd_simulation.R:131-160), so that no worker compiles it; a
                            ##  build that fails stops the run here, before any worker starts.)
                            ## (the same test as the fits' own, run_MULTIPLE_MetaOrd_simulation.R:110-114, never the
                            ##  machine's name: the name and the CPU can disagree, and then this build and the
                            ##  workers' would use different flags and every worker would rebuild)
                            CPU_has_AVX512_of_this_machine <- any(grepl("\\bavx512f\\b",
                                                                         readLines("/proc/cpuinfo", warn = FALSE)))
                            AVX_FLAGS_of_this_machine <- if (CPU_has_AVX512_of_this_machine) {
                                  "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
                            } else {
                                  "-march=native -mtune=native -mavx2"
                            }
                            # NicoSIMD_user_header <- file.path(dirname(dirname(local_pkg_dir)), "R_packages",
                            #                                   "NicoStan", "inst", "NicoSIMD", "inst", "include",
                            #                                   "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp")
                            NicoSIMD_user_header <- file.path(find.package("NicoSIMD"), "include", "NicoSIMD",
                                                              "stan_external_functions_of_NicoSIMD.hpp")
                            use_custom_file_of_this_scenario <- scenario_row$model_parameterisation %in%
                                  c("Bivariate", "ord_bivariate", "R&G")
                            custom_file_name_of_this_scenario <- if (scenario_row$model_parameterisation ==
                                                                   "Bivariate") {
                                  "DTA_MA_strat_bivariate_NicoStan_fused_AVX.stan"
                            } else if (scenario_row$model_parameterisation %in% c("ord_bivariate", "R&G")) {
                                  paste0(if (scenario_row$model_parameterisation == "ord_bivariate") "DTA_MA_ord_biv_" else
                                         "DTA_MA_ord_HSROC_",
                                         if (scenario_row$random_thresholds) "RANDthr_kappa" else "FIXEDthr",
                                         "_NicoStan_fused_AVX.stan")
                            } else {
                                  NULL
                            }
                            custom_cpp_user_header_of_this_scenario <-
                                  if (use_custom_file_of_this_scenario) {
                                        NicoSIMD_user_header
                                  } else {
                                        NULL
                                  }
                            ##
                            ## ---- The pop-median model, the centred cutpoints and their user headers in the build -------------------
                            ##
                            if (sim_study_type == "RC_ord_biv_median_vs_pop_median" &&
                                (scenario_row$RC_ord_biv_model == "pop_median" ||
                                 scenario_row$RC_cutpoints == "centred")) {
                                  custom_file_name_of_this_scenario <- scenario_row$custom_file_name[[1]]
                                  custom_cpp_user_header_of_this_scenario <-
                                        scenario_row$custom_cpp_user_header_file_path[[1]]
                            }
                            if (sim_study_type == "RC_ord_HSROC_median_vs_pop_median" &&
                                scenario_row$RC_ord_HSROC_model == "pop_median") {
                                  custom_file_name_of_this_scenario <- scenario_row$custom_file_name[[1]]
                                  custom_cpp_user_header_of_this_scenario <-
                                        scenario_row$custom_cpp_user_header_file_path[[1]]
                            }
                            Stan_model_name_of_this_scenario <- MetaOrdDTA:::R_fn_get_stan_model_file_name(
                                  custom_file_name              = custom_file_name_of_this_scenario,
                                  # use_Se_Sp_prior               = scenario_row$model_parameterisation %in%
                                  #                                 c("ord_bivariate", "R&G"),
                                  ## (as the worker: priors$use_Se_Sp_prior, run_MULTIPLE_MetaOrd_simulation)
                                  # use_Se_Sp_prior               = FALSE,
                                  use_Se_Sp_prior               = uses_the_Se_Sp_prior,
                                  cts                           = scenario_row$model_parameterisation == "Jones",
                                  network                       = FALSE,
                                  prior_only                    = FALSE,
                                  model_parameterisation        = scenario_row$model_parameterisation,
                                  random_thresholds             = scenario_row$random_thresholds,
                                  Dirichlet_random_effects_type = scenario_row$Dirichlet_random_effects_type)
                            model_file_name_of_this_scenario <-
                                  Stan_model_name_of_this_scenario$stan_model_file_name
                            user_header_of_this_scenario <- custom_cpp_user_header_of_this_scenario
                            start_of_the_build <- Sys.time()
                            if (use_custom_file_of_this_scenario) {
                                  build_of_this_scenario_advanced <-
                                        MetaOrdDTA:::R_fn_compile_stan_advanced_given_file_name(
                                              stan_model_file_name             = model_file_name_of_this_scenario,
                                              cts                              = FALSE,
                                              network                          = FALSE,
                                              prior_only                       = FALSE,
                                              debugging                        = FALSE,
                                              force_recompile                  = FALSE,
                                              quiet                            = TRUE,
                                              set_custom_CXX_CPP_flags         = TRUE,
                                              CCACHE_PATH                      = " ",
                                              custom_cpp_user_header_file_path = user_header_of_this_scenario,
                                              CXX_COMPILER_PATH                = NULL,
                                              CPP_COMPILER_PATH                = NULL,
                                              MATH_FLAGS                       = NULL,
                                              FMA_FLAGS                        = "-mfma",
                                              AVX_FLAGS                        = AVX_FLAGS_of_this_machine,
                                              THREAD_FLAGS                     = NULL)
                                  exe_file_of_this_scenario <-
                                        build_of_this_scenario_advanced$stan_model_obj$exe_file()
                            } else {
                                  build_of_this_scenario_basic <-
                                        MetaOrdDTA:::R_fn_compile_stan_model_basic_given_file_name(
                                              stan_model_file_name = model_file_name_of_this_scenario,
                                              cts                  = FALSE,
                                              network              = FALSE,
                                              prior_only           = FALSE,
                                              debugging            = FALSE,
                                              force_recompile      = FALSE,
                                              quiet                = TRUE,
                                              compile              = TRUE)
                                  exe_file_of_this_scenario <-
                                        build_of_this_scenario_basic$stan_model_obj$exe_file()
                            }
                            if (!file.exists(exe_file_of_this_scenario)) {
                                  stop(paste0("The model of this scenario was not built: ",
                                              exe_file_of_this_scenario, ": no worker is started."))
                            }
                            seconds_of_the_build <- as.numeric(difftime(Sys.time(), start_of_the_build,
                                                                        units = "secs"))
                            message(paste0("\033[36m", "The model of this scenario, built once for the workers: ",
                                           exe_file_of_this_scenario, " (",
                                           formatC(seconds_of_the_build, format = "f", digits = 1), " s)",
                                           "\033[0m"))
                            ## Create cluster and run

                                                      ## (a free port for the workers, never one another program holds, e.g. Ollama's 11434)
                            port_of_the_workers <-  NA
                            for (port in 11500:11999) {
                                  socket <-  tryCatch(serverSocket(port), error = function(e) NULL)
                                  if (!is.null(socket)) {
                                        close(socket)
                                        port_of_the_workers <-  port
                                        break
                                  }
                            }
                            ## Create cluster and run
                            cl <- makeCluster( n_workers,
                                               type = "PSOCK",
                                               # outfile = "",
                                               outfile = file.path(output_dir, "log_of_the_workers.txt"),
                                               port = port_of_the_workers)
                            # cl <- makeCluster( n_workers,
                            #                    type = "PSOCK",
                            #                    outfile = "")
                            registerDoParallel(cl)

                            clusterExport(cl,
                                          ls(envir = minimal_env),
                                          envir = minimal_env)
                            ## (every worker on the frozen NicoStan and NicoSIMD)
                            stopifnot(all(unlist(clusterEvalQ(cl, .libPaths()[1])) == frozen_library_of_this_sim))

                            # MCMC_seed = minimal_env$root_MCMC_seed + ii + batch_id * 1000
                            # sim_ids_chunk = sim_ids_chunk_list[[ii]]
                            # n_chains = n_chains_per_worker

                            ii <- 1

                            MCMC_seed = 123 + sims_this_batch[ii] + batch_id * 1000
                            sim_ids_chunk = sims_this_batch[ii]
                            min_runs_needed = 30
                            check_MCSE_interval = 1
                            n_chains <- 4

                            batch_results <- foreach(ii = seq_along(sims_this_batch),
                                                     .errorhandling = "pass") %dopar% {

                                                       tryCatch({

                                                             ## ---- A seed not yet started returns at once when this cell's MCSE target is met:
                                                             if (file.exists(file_of_the_MCSE_target_met)) {
                                                                   return(list(seed = sims_this_batch[ii],
                                                                               batch_id = batch_id,
                                                                               status = "not fitted: MCSE target met"))
                                                             }

                                                             ## Load packages
                                                             require(RcppParallel)
                                                             require(BayesMVPsim)
                                                             require(BayesMVP)
                                                             require(MetaOrdDTA)
                                                             require(cmdstanr)
                                                             require(dplyr)

                                                             set.seed(root_MCMC_seed + sims_this_batch[ii] + batch_id * 1000)
                                                             n_chains = n_chains_per_worker

                                                             ## Load configuration and helpers (once per worker: each cell starts new workers)
                                                             if (!exists("helpers_of_the_sim_study_are_sourced", envir = globalenv())) {
                                                                 cmdstanr::set_cmdstan_path()
                                                                 source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_config.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_1.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_2.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_helpers_3.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "run_MULTIPLE_MetaOrd_simulation.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "MetaOrdDTA_sim_run_batch_MetaOrd.R"))
                                                                 source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
                                                                 assign("helpers_of_the_sim_study_are_sourced", TRUE, envir = globalenv())
                                                             }

                                                             ## Call sim fn
                                                             result <- run_multiple_sims(  MCMC_seed = root_MCMC_seed + sims_this_batch[ii] + batch_id * 1000,
                                                                                           sim_ids_chunk = sims_this_batch[ii],
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
                                                                                           stop_early_if_target_met = TRUE,
                                                                                           sampling_algorithm = sampling_algorithm,
                                                                                           min_ESS_of_10_studies = min_ESS_of_10_studies,
                                                                                           min_ESS_of_50_studies = min_ESS_of_50_studies,
                                                                                           ordinal_prior_type = ordinal_prior_type)
                                                                                           ##
                                                                                           # vec_index_inner = vec_index_inner)

                                                             ## ---- Once this cell's MCSE target is met, no worker starts another of its seeds:
                                                             if (isTRUE(result$stopped_early)) {
                                                                   file.create(file_of_the_MCSE_target_met)
                                                             }

                                                             return(list(seed = sims_this_batch[ii],
                                                                         batch_id = batch_id,
                                                                         status = "completed",
                                                                         result = result))

                                                       }, error = function(e) {

                                                             ## Capture and return the actual error
                                                             cat(sprintf("\nSeed %d ERROR: %s\n",
                                                                         sims_this_batch[ii],
                                                                         e$message))
                                                             cat(sprintf("Traceback: %s\n", paste(capture.output(traceback()),
                                                                                                  collapse = "\n")))
                                                             return(list(seed = sims_this_batch[ii],
                                                                         batch_id = batch_id,
                                                                         status = "error",
                                                                         error = e$message))

                                                       })

                                                     }

                            stopCluster(cl)
                            unlink(file_of_the_MCSE_target_met)
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
                                                            # min_ESS_thr = if (scenario_row$n_studies == 10) 1000 else 200,
                                                            min_ESS_thr = min_ESS_of_this_cell,
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








































