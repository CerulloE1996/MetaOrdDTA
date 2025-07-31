


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

 
require(dplyr)
require(cmdstanr)


n_tests <- 5
n_index_tests <- n_tests - 1
##
##
assume_perfect_GS <- 1
 
set.seed(seed)
##
missing_indicator <- -1
##
n_covariates <- 1 ## intercept only !!
##
options(max.print = 100000)
##

##
## test_names <- c("GAD_2", "GAD_7", "HADS", "BAI")
## real: 29 39 36  7 (w/ 76 studies in total)
n_studies <- 80
n_studies_per_test <-   c( n_studies,
                           round(n_studies/2), round(n_studies/2), round(n_studies/2), round(n_studies/5))
# n_studies_per_test <-   c( 25, 
#                            15, 20, 15, 5)

 
# SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = 1.0, SD_prob_scale = 0.04102194)
# SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = 2.0, SD_prob_scale = 0.07832947)
# 
# seq_probit_scale <- seq(from = 0.01, to = 1.00, by = 0.005)
# for (i in 1:100) { 
#   SD_probit_scale <- seq_probit_scale[i]
#   print(paste("SD on prob scale for SD_probit_scale = ",  SD_probit_scale, "is equal to: ", SD_approx_probit_to_prob(SD_probit_scale)))
# }
# 

 
if (network == FALSE) {
        ##
        # Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
        sim_results <- R_fn_sim_data_ord_MA(          n_studies = n_studies,
                                                      N_per_study_mean = N_per_study_mean,
                                                      N_per_study_SD = N_per_study_SD,
                                                      assume_perfect_GS = assume_perfect_GS,
                                                      ##
                                                      seed = seed,
                                                      ##
                                                      true_Mean_prev = true_Mean_prev,
                                                      true_SD_probit_prev = true_SD_probit_prev)

} else { ## This is: NMA + covariates.
        
        if (intercept_only_for_sim_DGM == FALSE) {
          
                  outs_covs_subset <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_cov_data_4_tests"))
                  real_data <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_4_tests"))
                  ##
                  X <- list(outs_covs_subset$X_nd, 
                            outs_covs_subset$X_d)
                  ##
                  real_data$cov_data$X <- X
                  ##
                  # real_data$indicator_index_test_in_study
                  ##
                  if (MR_model_for_sim_DGM == "MR_model_1") { 
                    outs_covs_subset <- outs_covs_subset$MR_model_1
                  } else if (MR_model_for_sim_DGM == "MR_model_2") {
                    outs_covs_subset <- outs_covs_subset$MR_model_2
                  } else if (MR_model_for_sim_DGM == "MR_model_3") {
                    outs_covs_subset <- outs_covs_subset$MR_model_3
                  } else if (MR_model_for_sim_DGM == "MR_model_4") {
                    outs_covs_subset <- outs_covs_subset$MR_model_4
                  } else if (MR_model_for_sim_DGM == "MR_model_5") {
                    outs_covs_subset <- outs_covs_subset$MR_model_5
                  }
                  ##
                  n_covariates <- outs_covs_subset$n_covariates_max
                  ##
                  covariates_same_between_classes <- rep(TRUE, n_covariates)
                  ##
                  covariate_names <- outs_covs_subset$covariate_names
                  ##
                  continuous_covs <- outs_covs_subset$continuous_covs
                  binary_covs <- outs_covs_subset$binary_covs
                  categorical_covs <- outs_covs_subset$categorical_covs
                  ##
                  covariate_settings = list(
                      continuous_vars = continuous_covs,
                      binary_vars = binary_covs,
                      categorical_vars = categorical_covs,
                      covariate_names = covariate_names,
                      n_covariates = n_covariates
                  )
                  ##
                  X_nd <- outs_covs_subset$X_nd ## this is the real covariate data
                  X_d  <- outs_covs_subset$X_d  ## this is the real covariate data
                  X <- list(X_nd, X_d)
                  ##
                  X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies( X,
                                                                           indicator_index_test_in_study = real_data$indicator_index_test_in_study)
                  real_data$cov_data$X <- X
                  original_cov_data <- list(X = X)
                  ##
                  # n_studies = n_studies
                  N_per_study_mean = N_per_study_mean
                  N_per_study_SD = N_per_study_SD
                  assume_perfect_GS = assume_perfect_GS
                  ##
                  seed = seed
                  ##
                  true_Mean_prev = true_Mean_prev
                  true_SD_probit_prev = true_SD_probit_prev
                  ##
                  covariate_settings = covariate_settings
                  ##
                  simulate_covariates <- TRUE
              
        } else {
              
                  covariate_settings = list(
                    covariate_names = c("intercept"),
                    n_covariates = 1
                  )
                  ##
                  X_nd_t <- array(1, dim = c(n_studies, 1))
                  X_nd <- rep(list(X_nd_t), n_index_tests)
                  X_d  <- X_nd
                  X <- list(X_nd, X_d)
                  ##
                  original_cov_data <- list(X = X)
                  ##
                  # n_studies = n_studies
                  N_per_study_mean = N_per_study_mean
                  N_per_study_SD = N_per_study_SD
                  assume_perfect_GS = assume_perfect_GS
                  ##
                  seed = seed
                  ##
                  true_Mean_prev = true_Mean_prev
                  true_SD_probit_prev = true_SD_probit_prev
                  ##
                  covariate_settings = covariate_settings
                  ##
                  simulate_covariates <- FALSE
              
        }

          ##
        require_index_test = TRUE 
        missing_value_covs = -999
        ##
        sim_results <- R_fn_sim_NMA_data_with_covariates_varying( 
                                           # n_studies = n_studies,
                                           n_studies_per_test = n_studies_per_test,
                                           N_per_study_mean = N_per_study_mean,
                                           N_per_study_SD = N_per_study_SD,
                                           assume_perfect_GS = assume_perfect_GS,
                                           ##
                                           seed = seed,
                                           ##
                                           true_Mean_prev = true_Mean_prev,
                                           true_SD_probit_prev = true_SD_probit_prev,
                                           ##
                                           covariate_settings = covariate_settings,
                                           original_cov_data = original_cov_data,
                                           simulate_covariates = simulate_covariates,
                                           missing_value_covs = missing_value_covs)
        
        n_studies <- sim_results$n_studies
                                         
        if (intercept_only_for_sim_DGM == FALSE) {
                ##
                ## ---- Extract covariate data (SAME for each test):
                ##
                cov_data_list <- list()
                ##
                X <- sim_results$X_mat
                # X_list <- rep(list(X), n_index_tests)
                cov_data_list$X_nd <- X[[1]]
                cov_data_list$X_d  <- X[[2]]
                ##
                cov_data_list$n_covariates_max <- ncol(X)
                cov_data_list$n_covariates_nd <- rep(cov_data_list$n_covariates_max, n_index_tests)
                cov_data_list$n_covariates_d <- cov_data_list$n_covariates_nd 
                ##
                # X[1:10,]
                baseline_case <- c(1,              ## intercept
                                   ##
                                   0.0,   ## prev_GAD
                                   # mean(X[, 3]),   ## prev_AAD
                                   ##
                                   0, ## low_RoB_QUADAS
                                   # 0, ## low_RoB_liberal
                                   ##
                                   0, ## Ref_test_SCID
                                   0, ## Ref_test_Structured
                                   ##
                                   0, ## study_setting_2
                                   1) ## study_setting_3 (most common)
                ##
                cov_data_list$baseline_case_nd <- baseline_case
                cov_data_list$baseline_case_d  <- baseline_case
              
        } else { 
                ##
                ## ---- Extract covariate data (SAME for each test):
                ##
                cov_data_list <- list()
                ##
                X <- sim_results$X_mat
                X_list <- rep(list(X), n_index_tests)
                cov_data_list$X_nd <- cov_data_list$X_d <- X_list
                ##
                cov_data_list$n_covariates_max <- ncol(X)
                cov_data_list$n_covariates_nd <- rep(cov_data_list$n_covariates_max, n_index_tests)
                cov_data_list$n_covariates_d <- cov_data_list$n_covariates_nd 
                ##
                # X[1:10,]
                baseline_case <- c(1)
                ##
                cov_data_list$baseline_case_nd <- baseline_case
                cov_data_list$baseline_case_d  <- baseline_case
                
        }
        

}


# sim_results$thr_for_all_tests_for_all_studies_nd[1, 2, 1:6]
# sim_results$thr_for_all_tests_for_all_studies_nd[2, 2, 1:6]
# sim_results$thr_for_all_tests_for_all_studies_nd[3, 2, 1:6]
# sim_results$thr_for_all_tests_for_all_studies_nd[4, 2, 1:6]
# sim_results$thr_for_all_tests_for_all_studies_nd[5, 2, 1:6]

             
sim_results$true_DGM_Se                                            
 # index_test_chosen_index <- 1 + 1
 # #index_test_chosen_index <- 4 + 1

{
  
    y_list <- sim_results$y_list
    str(y_list)
    
    sim_results$n_thr_per_test
    sim_results$N_per_study_vec

    # n_thr <- 10
    n_thr_vec <- sim_results$n_thr_per_test_excl_ref ; n_thr_vec
    # ## subset to remove ref test:
    # n_thr <- n_thr[-c(1)]
    n_cat_vec <- n_thr_vec + 1 ; n_cat_vec
    n_thr_max <- max(n_thr_vec) ; n_thr_max
    n_cat_max <- max(n_cat_vec) ; n_cat_max
    
}




 
 sim_results$true_DGM_Fp
 sim_results$true_DGM_Sp[[1]]
 sim_results$true_DGM_Se[[1]]

{
  true_DGM_Fp <- 100 * sim_results$true_DGM_Fp[[index_test_chosen_index - 1]]
  true_DGM_Sp <- 100 * sim_results$true_DGM_Sp[[index_test_chosen_index - 1]]
  true_DGM_Se <- 100 * sim_results$true_DGM_Se[[index_test_chosen_index - 1]]
}
  
 
 
 


x_nd <- x_d <- list()
for (t in 2:n_tests) {
    n_thr_chosen_index <- n_thr_vec[t - 1]
    x_nd[[t - 1]] <- sim_results$n_FP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)] ## apply(sim_results$n_FP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)], c(1, 2), mean)
    x_d[[t - 1]]  <- sim_results$n_TP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)] ## apply(sim_results$n_TP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)], c(1, 2), mean)
}


x_nd

str(sim_results$n_FP_at_each_threshold)

sim_results$n_FP_at_each_threshold
 
###
### Prep / save final NMA simulated data outputs:
x_NMA <- list(x_nd, x_d)
x <- x_NMA

x[[2]][[2]]
 
###
### ---- Extract specific index tests:
####
x_GAD2 <- x_HADS <- x_BAI <- x_GAD7 <- list()
##
GAD2_index <- 1
GAD7_index <- 2
HADS_index <- 3
BAI_index  <- 4

##
for (c in 1:2) {
  x_GAD2[[c]] <- x[[c]][[GAD2_index]]
  x_GAD7[[c]] <- x[[c]][[GAD7_index]]
  x_HADS[[c]] <- x[[c]][[HADS_index]]
  x_BAI[[c]]  <- x[[c]][[BAI_index]]
}
x_GAD7
x_GAD2
x_HADS
BAI_index
 


#### ------------------- Apply ** MISSING THRESHOLDS ** to GAD-2 data/test: -----------------------------------------------------------------------------



source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
# source(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "make_Klaus_et_al_data.R")) ## to get "x_gad2" dataset
##
x_GAD2_REAL <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_x_GAD2.RDS"))
##
x_GAD2_w_missing_thr <- apply_missingness_pattern_NMA(   x_complete = x_GAD2,
                                                         x_pattern = x_GAD2_REAL,
                                                         enforce_consecutive_missingness = FALSE,
                                                         x_complete_missing_indicator = -1,
                                                         x_pattern_missing_indicator = -1)
##
x_GAD2_w_missing_thr
x_GAD2
x_GAD2_REAL


## Compute the % of missing thr data:
100 * sum(x_GAD2_w_missing_thr[[1]] == -1) / (n_studies * (ncol(x_GAD2_w_missing_thr[[1]]) - 1))
100 * sum(x_GAD2_REAL[[1]] == -1) / (n_studies * (-1 + ncol(x_GAD2_REAL[[1]])))

 

#### ------------------- Apply ** MISSING THRESHOLDS ** to HADS data/test: -----------------------------------------------------------------------------



source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
# source(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "make_Klaus_et_al_data.R")) ## to get "x_gad2" dataset
##
x_HADS_REAL <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_x_HADS.RDS"))
##
x_HADS_w_missing_thr <- apply_missingness_pattern_NMA( x_complete = x_HADS,
                                                   x_pattern = x_HADS_REAL,
                                                   enforce_consecutive_missingness = FALSE,
                                                   x_complete_missing_indicator = -1,
                                                   x_pattern_missing_indicator = -1)
##
x_HADS_w_missing_thr
x_HADS
x_HADS_REAL


## Compute the % of missing thr data:
100 * sum(x_HADS_w_missing_thr[[1]] == -1) / (n_studies * (ncol(x_HADS_w_missing_thr[[1]]) - 1))
100 * sum(x_HADS_REAL[[1]] == -1) / (n_studies * (ncol(x_HADS_REAL[[1]]) - 1))

#### ------------------- Apply ** MISSING THRESHOLDS ** to BAI data/test: -----------------------------------------------------------------------------



source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
# source(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "make_Klaus_et_al_data.R")) ## to get "x_gad2" dataset
##
x_BAI_REAL <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_x_BAI.RDS"))
##
x_BAI_w_missing_thr <- apply_missingness_pattern_NMA(  x_complete = x_BAI,
                                                   x_pattern = x_BAI_REAL,
                                                   enforce_consecutive_missingness = FALSE,
                                                   x_complete_missing_indicator = -1,
                                                   x_pattern_missing_indicator = -1)
##
x_BAI_w_missing_thr
x_BAI
x_BAI_REAL


## Compute the % of missing thr data:
100 * sum(x_BAI_w_missing_thr[[1]] == -1) / (n_studies * (ncol(x_BAI_w_missing_thr[[1]]) - 1))
100 * sum(x_BAI_REAL[[1]] == -1) / (n_studies * (ncol(x_BAI_REAL[[1]]) - 1))


#### ------------------- Apply ** MISSING THRESHOLDS ** to GAD-7 data/test: -----------------------------------------------------------------------------



source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
# source(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "make_Klaus_et_al_data.R")) ## to get "x_gad2" dataset
##
x_GAD7_REAL <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_x_GAD7.RDS"))
##
x_GAD7_w_missing_thr <- apply_missingness_pattern_NMA( x_complete = x_GAD7,
                                                   x_pattern  = x_GAD7_REAL,
                                                   enforce_consecutive_missingness = FALSE,
                                                   x_complete_missing_indicator = -1,
                                                   x_pattern_missing_indicator = -1)
##
x_GAD7_w_missing_thr
x_GAD7
x_GAD7_REAL


## Compute the % of missing thr data:
100 * sum(x_GAD7_w_missing_thr[[1]] == -1) / (n_studies * (ncol(x_GAD7_w_missing_thr[[1]]) - 1))
100 * sum(x_GAD7_REAL[[1]] == -1) / (n_studies * ncol(x_GAD7_REAL[[1]]))

#### ------------------- Make NMA dataset w/ missing thresholds for some tests --------------------------------------------------------------------------



###
### Prep / save final NMA simulated data outputs:
x_NMA <- list(x_nd, x_d)
x_NMA_w_missing_thr <- x_NMA
##
## Fpr GAD-2 (index test #1)
##
str(x_NMA)
for (c in 1:2) { 
  x_NMA_w_missing_thr[[c]][[1]] <- x_GAD2_w_missing_thr[[c]]
}
##
## Fpr BAI (index test #2)
##
for (c in 1:2) { 
  x_NMA_w_missing_thr[[c]][[2]] <- x_GAD7_w_missing_thr[[c]]
}
##
## Fpr HADS (index test #3)
##
for (c in 1:2) { 
  x_NMA_w_missing_thr[[c]][[3]] <- x_HADS_w_missing_thr[[c]]
}
##
## For GAD-7 (index test #4)
##
for (c in 1:2) { 
  x_NMA_w_missing_thr[[c]][[4]] <- x_BAI_w_missing_thr[[c]]
}

 











