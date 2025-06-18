

 
 
# # # # # # # ##
rm(list = ls())


{

os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
#


 {
  ## First remove any possible package fragments:
  ## Find user_pkg_install_dir:
  user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
  print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
  ##
  ## Find pkg_install_path + pkg_temp_install_path:
  pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
  pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
  ##
  ## Remove any (possible) MetaOrdDTA package fragments:
  remove.packages("MetaOrdDTA")
  unlink(pkg_install_path, recursive = TRUE, force = TRUE)
  unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
}

#


#
#### -------- ACTUAL (LOCAL) INSTALL:
## Document:
devtools::clean_dll(local_pkg_dir)
roxygen2::roxygenize(local_pkg_dir)
##
## Install (outer pkg):
##
devtools::clean_dll(local_pkg_dir)
devtools::install(local_pkg_dir,
                  upgrade = "never",
                  quick = TRUE)
##
## May need to restart R:
##
.rs.restartR()  # In RStudio

# ?devtools::install


}
#
# # 
# 



  

start_sim_index <- 1



  # 
#### ----
{
  model_parameterisation = "Bivariate"
  ##
  softplus <- FALSE
  ##
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
}


{
  model_parameterisation = "Jones"
  ##
  softplus <- FALSE
  ##
  box_cox <- TRUE
  cts <- TRUE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
}


{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
}


{
  model_parameterisation = "R&G"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
}


{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  TRUE
  ##"
  Dirichlet_random_effects_type <- "kappa"
}


{
  model_parameterisation = "R&G"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  TRUE
  ##"
  Dirichlet_random_effects_type <- "kappa"
}











{

os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")


setwd(local_pkg_dir)

{
  require(MetaOrdDTA)
  require(RcppParallel)
  require(BayesMVP)
}

source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))



{
    
    n_chains <-  8 ##  parallel::detectCores() / 2
    ##
    if (random_thresholds)    n_iter   <- ceiling(1000/(n_chains/4)) ## rand-thr actually have better ESS/sec 
    if (!(random_thresholds)) n_iter   <- ceiling(1000/(n_chains/4))
    n_iter
    ##
    n_burnin <- 500
    ##
    ##  max_treedepth <- 9
     max_treedepth <- 10
    # 2^max_treedepth
      ##  adapt_delta <- 0.80
      adapt_delta <- 0.65

}



##  N_per_study_mean <- 500  ; n_studies <- 10
 N_per_study_mean <- 500  ; n_studies <- 50
##
N_per_study_SD <- N_per_study_mean
## N_per_study_SD <- 1 ## N_per_study_mean




##
    index_test_MA <- "GAD_2" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "GAD_2" ;  missing_thr <- FALSE   ## FOR TESTING
##
   ##    index_test_MA <- "HADS" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "HADS" ;  missing_thr <- FALSE   ## FOR TESTING
##
##  index_test_MA <- "GAD_7" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "GAD_7" ;  missing_thr <- FALSE   ## FOR TESTING


alpha_pi <- "flat"
## alpha_pi <- "weak_info"


{
  true_Mean_prev <- 0.25 ;  true_SD_probit_prev <- 0.25
  ##
   ############# scale_ALL_bs_SDs_by <- 0.5
   ###   scale_ALL_bs_SDs_by <- 1.0
  ##
 ###   bs_het_C_nd_prob_scale <- 0.03932276 ; bs_het_C_d_prob_scale <- 0.07712564 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
}


if (index_test_MA == "GAD_2") {
        ##
        ## ---- based on "kappa" Xu model fitted to real GAD-2 data:
        ##
        # bivariate_locations_nd_bs_het_SD <- 0.343
        # bivariate_locations_d_bs_het_SD  <- 0.553 
        ##
        index_test_chosen_index <- 1 + 1
        n_thr <- 6
        ##
        vec_index_inner <- 1:n_thr ## NULL
        compute_sim_study_metrics <- NULL
 
          
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
 
          
} else if (index_test_MA == "HADS") {
          ##
          ## ---- based on "kappa" Xu model fitted to real HADS data
          ##
          # bivariate_locations_nd_bs_het_SD <- 0.235
          # bivariate_locations_d_bs_het_SD  <- 0.396
          index_test_chosen_index <- 3 + 1
          n_thr <- 21
          ##
          vec_index_inner <- 3:17
          n_inner_thr <- length(vec_index_inner)
          vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
          n_outer_thr <- length(vec_index_outer_thr)
          # compute_sim_study_metrics <- 1
 
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
  
} else if (index_test_MA == "GAD_7") {
          ##
          ## ---- based on "kappa" Xu model fitted to real GAD-7 data:
          ##
          # bivariate_locations_nd_bs_het_SD <- 0.488
          # bivariate_locations_d_bs_het_SD  <- 0.669 
          index_test_chosen_index <- 4 + 1
          n_thr <- 21
          ##
          vec_index_inner <- 3:17
          n_inner_thr <- length(vec_index_inner)
          vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
          n_outer_thr <- length(vec_index_outer_thr)
          # compute_sim_study_metrics <- 1
          ##
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
  
}




 

## 
if (n_studies < 50) { 
    n_sims <- 250
} else { 
    n_sims <- 50
}

##



{
    print(paste("--------------------------------------------------------------"))
    print(paste("n_sims = ", n_sims))
    print(paste("n_studies = ", n_studies))
    print(paste("N_per_study_mean = ", N_per_study_mean))
    ##
    # print(paste("--------------------------------------------------------------"))
    # print(paste("bivariate_locations_nd_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD))
    # print(paste("bivariate_locations_d_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD))
    # ##
    # print(paste("--------------------------------------------------------------"))
    # if (bs_het_C_nd_prob_scale > 0.01) { 
    #     cutpoint_SD <- "significant het"
    # } else { 
    #     cutpoint_SD <- "No het"
    # }
    # print(paste("cutpoint_SD = ", cutpoint_SD))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("model_parameterisation = ", model_parameterisation))
    print(paste("--------------------------------------------------------------"))
    print(paste("random_thresholds = ", random_thresholds))
    print(paste("--------------------------------------------------------------"))
    print(paste("index_test_MA = ", index_test_MA))
    print(paste("--------------------------------------------------------------"))
    print(paste("Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))
    print(paste("alpha_pi = ", alpha_pi))
    ##
    print(paste("--------------------------------------------------------------"))
}






{
  
  output_dir = "simulation_results"
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # # Initialize storage for results
  # # Set up summary storage
  # summary_df <- create_summary_df(n_sims = n_sims,
  #                                 min_k = vec_index_inner[1],
  #                                 max_k = tail(vec_index_inner, 1),
  #                                 start_index = start_sim_index)
  
  network <- FALSE
  ##
  softplus <- FALSE
  # softplus <- TRUE
  
  
}



# file_name <- "DTA_MA_Xu_FIXEDthr_daig_J.stan"
##
# use_custom_file <- TRUE
use_custom_file <- FALSE




}




start_sim_index
n_sims <- 500


# 
source(file.path("inst", "examples", "MetaOrdDTA_MA_eg_1_debug_v4.R"))
# 
min_ESS
# # priors
# alpha_pi
# 
# 
# 
# 
# 
# 
# 
# 



dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)


dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline"))) %>% print(n = 1000)



quantile(exp(rnorm(n = 10000, mean = log(2*100), sd = 1)), c(0.025, 0.50, 0.975))
quantile(exp(rnorm(n = 10000, mean = log(100), sd = 1)), c(0.025, 0.50, 0.975))

# Generate samples from the prior
n_samples <- 10000

# Sample log(kappa) from student_t(df=3, location=log(50), scale=1)
log_kappa_samples <- rt(n_samples, df = 5) * 1 + log(200)
# Note: rt() generates standard t-distribution, so we scale and shift

# Transform to get kappa
kappa_samples <- exp(log_kappa_samples)

# Quick summary
summary(kappa_samples)
quantile(kappa_samples, probs = c(0.025, 0.5, 0.975))

# readRDS("simulation_results/test_GAD_7N_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flatscale_het_1het_nd_0.488het_d_0.669cut_SD_1_complete_results_list.RDS")
# 
# readRDS("simulation_results/test_GAD_7N_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flatscale_het_1het_nd_0.488het_d_0.669cut_SD_1.RDS")
# 
# 



index_test_MA
##
## ---- For GAD-2:
##
test_GAD_2N_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_flat.RDS
test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_flat.RDS
##
## ---- Bivariate:
##
file_name <- "test_GAD_2N_500n_studies_50param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_metrics_overall.csv"
##
## ---- Jones:
##
file_name <- "test_GAD_2N_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_metrics_overall.csv"
##
## ---- Xu, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_flat_metrics_overall.csv"
##
## ---- R&G, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flat_metrics_overall.csv"
##
## ---- Xu, Random C:
##
file_name <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_flat_metrics_overall.csv"
##
## ---- R&G, Random C:
##
file_name <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_flat_metrics.RDS"
csv_name  <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_flat_metrics_overall.csv"

 
## need to run for GAD-2, N_s 50:
##  Bivariate
##  Jones
##  Xu, fixed
##  R&G, random


# test_HADSN_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
# test_HADSN_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
# test_HADSN_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flat.RDS
# test_HADSN_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_flat.RDS
# test_HADSN_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_flat.RDS


index_test_MA



complete_results_list <- readRDS(file.path("simulation_results", file_name))
complete_results_list

csv_file <- read.csv(file.path("simulation_results", csv_name))
csv_file



# Usage example:
metrics <- compute_simulation_metrics(complete_results_list = complete_results_list,
                                      true_DGM_Se = true_DGM_Se, 
                                      true_DGM_Sp = true_DGM_Sp,
                                      min_k = vec_index_inner[1],
                                      max_k = tail(vec_index_inner, 1))

# Create nice tables
tables <- create_metrics_summary_table(metrics = metrics,
                                       min_k = vec_index_inner[1], 
                                       max_k = tail(vec_index_inner, 1))


# Now you can use calc_min_sim_size



{
    # Print results
    cat("\n========== SIMULATION METRICS SUMMARY ==========\n")
    cat(sprintf("Number of simulations: %d\n", metrics$n_sims))
    cat(sprintf("Thresholds analyzed: %d to %d\n\n", vec_index_inner[1], tail(vec_index_inner, 1)))
    
    cat("--- Overall Performance (averaged across thresholds) ---\n")
    print(tables$overall, row.names = FALSE)
    
    cat("\n--- Performance by Threshold ---\n")
    print(tables$by_threshold, row.names = FALSE)
}


summary_df$model_parameterisation[1]
summary_df$random_thresholds[1]
(summary_df$seed)

##
## ---- 
 

{
  
  # cat(paste("\n--------------------------------------------------\n"))
  # message("if min_MCSE_pct = 0.5%:")
  # print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec,
  #                         min_MCSE_pct = 0.50,
  #                         min_k = vec_index_inner[1],
  #                         max_k = tail(vec_index_inner, 1)))
  # cat(paste("\n--------------------------------------------------\n"))
  cat(paste("\n--------------------------------------------------\n"))
  message("if min_MCSE_pct = 0.25%:")
  print(calc_min_sim_size(
    SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
    min_MCSE_pct = 0.25, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if min_MCSE_pct = 0.125%:")
  print(calc_min_sim_size(
    SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
    min_MCSE_pct = 0.125, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if min_MCSE_pct = 0.1%:")
  print(calc_min_sim_size(
    SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
    min_MCSE_pct = 0.10, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if min_MCSE_pct = 0.05%:")
  print(calc_min_sim_size(
    SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
    min_MCSE_pct = 0.05, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  
  cat(paste("\n--------------------------------------------------\n"))
}


# ##
{
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.25%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.25, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.125%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.125, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.1%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.10, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.05%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.05, 
    min_k = vec_index_inner[1],
    max_k = tail(vec_index_inner, 1)
  ))
  
  cat(paste("\n--------------------------------------------------\n"))
}




