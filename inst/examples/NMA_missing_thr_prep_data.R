


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

 
require(dplyr)
require(cmdstanr)


n_tests <- 5
n_index_tests <- n_tests - 1
##
##
assume_perfect_GS <- 1
# ##
# true_Mean_prev <- 0.25
# ##
# true_SD_probit_prev <- 0.25
# # true_SD_probit_prev <- 0.10
# # true_SD_probit_prev <- 0.05
# # true_SD_probit_prev <- 0.025
# ## true_SD_probit_prev <- 0.01
##
set.seed(seed)
##
missing_indicator <- -1
##
n_covariates <- 1 ## intercept only !!
##
options(max.print = 100000)
##
# 

# pnorm(qnorm(true_Mean_prev) + 2.0*true_SD_probit_prev)
# pnorm(qnorm(true_Mean_prev) - 2.0*true_SD_probit_prev)
# 
# mean_test <- -0.1
# pnorm(mean_test + 2*bivariate_locations_d_bs_het_SD)
# pnorm(mean_test)
# pnorm(mean_test - 2*bivariate_locations_d_bs_het_SD)


# SD_approx_probit_to_prob(0.50)
# SD_approx_probit_to_prob(0.25)
# SD_approx_probit_to_prob(0.16)
# SD_approx_probit_to_prob(0.125)
# SD_approx_probit_to_prob(0.10)
##
 
SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = 1.0, SD_prob_scale = 0.04102194)
SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = 2.0, SD_prob_scale = 0.07832947)

seq_probit_scale <- seq(from = 0.01, to = 1.00, by = 0.005)
for (i in 1:100) { 
  SD_probit_scale <- seq_probit_scale[i]
  print(paste("SD on prob scale for SD_probit_scale = ",  SD_probit_scale, "is equal to: ", SD_approx_probit_to_prob(SD_probit_scale)))
}



# n_studies = n_studies
# N_per_study_mean = N_per_study_mean
# N_per_study_SD = N_per_study_SD
# assume_perfect_GS = assume_perfect_GS
# ##
# seed = seed
# ##
# true_Mean_prev = true_Mean_prev
# true_SD_probit_prev = true_SD_probit_prev
 

## mean(dirichlet_cat_SDs_sigma_nd) =  0.04387889 -> SD on probit scale ~ 0.17 (O-biv-RC model fitted to GAD-7 data)
## mean(dirichlet_cat_SDs_sigma_nd) =  0.06703344 -> SD on probit scale ~ 0.225 (O-biv-RC model fitted to GAD-7 data)


if (covariates == FALSE) {
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

        # n_studies = n_studies
        # N_per_study_mean = N_per_study_mean
        # N_per_study_SD = N_per_study_SD
        # assume_perfect_GS = assume_perfect_GS
        # ##
        # seed = seed
        # ##
        # true_Mean_prev = true_Mean_prev
        # true_SD_probit_prev = true_SD_probit_prev
        ##
        ##
        covariate_settings = list(
            continuous_vars = c("logit_prev_GAD"),
            ##
            binary_vars = c("low_RoB_QUADAS_clean"),  # Match actual column name
            ##
            categorical_vars = c("Ref_test_clean", "study_setting"),  # Match actual pattern
            ##
            covariate_names = c("intercept", 
                                ##
                                "logit_prev_GAD", 
                                ##
                                "low_RoB_QUADAS_clean", 
                                ##
                                "study_setting_2",  
                                "study_setting_3", 
                                ##
                                "Ref_test_clean_SCID",
                                "Ref_test_clean_Structured"),
            n_covariates = 7
        )
        ##
        X_nd <- outs_covs_subset$X_nd
        X_d  <- outs_covs_subset$X_d
        X <- list(X_nd, X_d)
        original_cov_data <- list(X = X)
        ##
        # n_studies = n_studies
        # N_per_study_mean = N_per_study_mean
        # N_per_study_SD = N_per_study_SD
        # assume_perfect_GS = assume_perfect_GS
        # ##
        # seed = seed
        # ##
        # true_Mean_prev = true_Mean_prev
        # true_SD_probit_prev = true_SD_probit_prev
        # ##
        # covariate_settings = covariate_settings
        # ##
        # simulate_covariates <- TRUE
        ##
        ## test_names <- c("GAD_2", "BAI", "HADS", "GAD_7")
        n_studies_per_test <- c(50, 
                                40, 10, 25, 25)
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
                                           original_cov_data = original_cov_data)
        
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
        X[1:10,]
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
                           0) ## study_setting_3
        ##
        cov_data_list$baseline_case_nd <- baseline_case
        cov_data_list$baseline_case_d  <- baseline_case
        ##

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
 
 # round(true_DGM_Sp, 3)
 # round(true_DGM_Se, 3)
 # 
 # sim_results$n_TP_at_current_threshold_OVERALL
 # sim_results$n_FP_at_current_threshold_OVERALL
 # 
 # sim_results$n_TP_at_each_threshold[, 1 + 1, 1:7]
 # sim_results$n_FP_at_each_threshold[, 1 + 1, 1:7]
 # 
 # 
 # sim_results$Se_all_tests_all_thresholds[, 1 + 1, 1:7]
 # sim_results$Sp_all_tests_all_thresholds[, 1 + 1, 1:7]
 # 
 # 
 # 
 # 
 # apply(sim_results$n_TP_at_each_threshold[, 1 + 1, 1:7], c(1, 2), mean)
 # apply(sim_results$n_FP_at_each_threshold[, 1 + 1, 1:7], c(1, 2), mean)
 # 
 # 
 # apply(sim_results$Se_all_tests_all_thresholds[, 1 + 1, 1:7], c(1, 2), mean)
 # 
 # apply(sim_results$Sp_all_tests_all_thresholds[, 1 + 1, 1:7], c(1, 2), mean)
 # 
 # 
 # 
 
 
 
 
# 
# sim_results$Se_OVERALL_all_tests_all_thresholds[2, ]
# 
# sim_results$Se_OVERALL_all_tests_all_thresholds
# # sim_results$Se_OVERALL_all_tests_all_thresholds[1, 1:n_thr_vec[1]]
# sim_results$Se_OVERALL_all_tests_all_thresholds[2, 1:n_thr_vec[2]]
 
 # true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr_vec[index_test_chosen_index - 1]]
 # true_Se_OVERALL
 # ##
 # true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr_vec[index_test_chosen_index - 1]]
 # true_Sp_OVERALL

  # plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

# sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
##


# min(y_list[[1]][, 5])
# max(y_list[[1]][, 5])
# 
# ##
# y_list_example_1 <- list()
# for (s in 1:n_studies) { 
#     N <- nrow(y_list[[s]])
#     y <- array(NA, dim = c(N, n_index_tests + 1))
#     for (t in 1:(n_index_tests + 1)) { 
#       y[, t] <- y_list[[s]][, t]
#     }
#     y_list_example_1[[s]] <- y
# }
# 
# min(y[,index_test_chosen_index])
# max(y[,index_test_chosen_index])
# 
# 
# str(y_list_example_1)
# 
# str(y_list)
# # y_list_example_1
# 

# ##
# plot(density(y[,index_test_chosen_index]))
# 
# y_nd <- y[,index_test_chosen_index][y[, 1] == 0]
# y_d <- y[,index_test_chosen_index][y[, 1] == 1]
# ##
# par(mfrow  = c(2, 1))
# plot(density(y_nd))
# plot(density(y_d))

## plot(density(log(y[,index_test_chosen_index])))
# 
# ## Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
# agg_data_cumulative <- convert_to_aggregate_counts( y_list = y_list_example_1, 
#                                                     n_studies = n_studies,
#                                                     n_tests_inc_ref = n_index_tests + 1,
#                                                     n_thr = n_thr_vec)
# 
# 
# 
# agg_data_cumulative$Se_per_study_list[[1]]
# agg_data_cumulative$x_nd_list[[1]]
# 
# 
# ## Get OBSERVED "true" values BEFORE you remove any thresholds and BEFORE you remove any tests for a given study s!!!
# Se_per_study_list <- agg_data_cumulative$Se_per_study_list
# Sp_per_study_list <- agg_data_cumulative$Sp_per_study_list
# 
# true_Se_OVERALL <- true_Sp_OVERALL <- list() 
# true_Se_OVERALL_weighted <- true_Sp_OVERALL_weighted <- list() 
# 
# agg_data_cumulative$x_d_list
# 
# 
# sim_results$n_TP_at_each_threshold[, 1 + 1, 1:7]
# sim_results$n_FP_at_each_threshold[, 1 + 1, 1:7]
# 
# agg_data_cumulative$x_nd_list[[1]]


# 0.9897257 0.9566693 0.8985080 0.7359064 0.3932815 
# 0.2462610 0.5247035 0.7138989 0.9045126 0.9878030  # Sp

# for (t in 1:n_index_tests) {
#   
#       ##
#       ## this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
#       ##
#       # true_Se_OVERALL[[t]] <- 100 * colMeans(Se_per_study_list[[t]]) ; true_Se_OVERALL[[t]]
#       # true_Sp_OVERALL[[t]] <- 100 * colMeans(Sp_per_study_list[[t]]) ; true_Sp_OVERALL[[t]]
#       ##
#       ## this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
#       ##
#       true_Se_OVERALL_weighted[[t]] <- sim_results$Se_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Se_OVERALL_weighted[[t]]
#       true_Sp_OVERALL_weighted[[t]] <- sim_results$Sp_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Sp_OVERALL_weighted[[t]]
#       ##
#       ##
#       # print(true_Se_OVERALL[[t]] - true_Se_OVERALL_weighted[[t]])
#       # print(true_Sp_OVERALL[[t]] - true_Sp_OVERALL_weighted[[t]])
# 
# }
# 
# 
# true_Se_OVERALL[[3]] - true_Se_OVERALL_weighted[[3]]
# true_Sp_OVERALL[[3]] - true_Sp_OVERALL_weighted[[3]]
# 
# 
# # true_Sp_OVERALL_weighted
# # true_Sp_OVERALL
# 
# 
# true_Se_OVERALL_weighted[[1]]
# true_Se_OVERALL[[1]]
 
 
 
 
 
 
 
 

 
 
 
 

# #### ------------------- Apply ** MISSING TESTS ** (hence "NMA" - optional): ---------------------------------------------------------

outs_NMA <- MetaOrdDTA:::create_test_in_study_for_NMA( seed = seed,
                                          n_studies = n_studies,
                                          n_tests_inc_ref = n_index_tests + 1,
                                          prob_present = 0.40,
                                          min_index_tests_per_study = 1)

indicator_index_test_in_study <- outs_NMA$indicator_index_test_in_study ; indicator_index_test_in_study
n_index_tests_per_study       <- outs_NMA$n_index_tests_per_study ; n_index_tests_per_study

# 
# ## Then update the "true values" list:
# for (t in 1:n_index_tests) {
#   indicator_mat <- array(1.0, dim = c(n_studies, n_thr_vec[t]))
#   for (k in 1:n_thr_vec[t]) {
#      indicator_mat[,k] <- ifelse(indicator_index_test_in_study[,t] == 0, NA, 1)
#   }
#   Se_per_study_list[[t]] <- Se_per_study_list[[t]] * indicator_mat
# }
# # 
# 
# ## Then re-calculate:
# for (t in 1:n_index_tests) {
#   ##
#   ## this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
#   ##
#   true_Se_OVERALL[[t]] <- 100 * colMeans(Se_per_study_list[[t]], na.rm = TRUE) ; true_Se_OVERALL[[t]]
#   true_Sp_OVERALL[[t]] <- 100 * colMeans(Sp_per_study_list[[t]], na.rm = TRUE) ; true_Sp_OVERALL[[t]]
#   ##
#   ## this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
#   ##
#   true_Se_OVERALL_weighted[[t]] <- sim_results$Se_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Se_OVERALL_weighted[[t]]
#   true_Sp_OVERALL_weighted[[t]] <- sim_results$Sp_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Sp_OVERALL_weighted[[t]]
#   ##
#   ##
#   print(true_Se_OVERALL[[t]] - true_Se_OVERALL_weighted[[t]])
#   print(true_Sp_OVERALL[[t]] - true_Sp_OVERALL_weighted[[t]])
#   
# }


x_nd <- x_d <- list()
for (t in 2:n_tests) {
    n_thr_chosen_index <- n_thr_vec[t - 1]
    x_nd[[t - 1]] <- apply(sim_results$n_FP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)], c(1, 2), mean)
    x_d[[t - 1]]  <- apply(sim_results$n_TP_at_each_threshold[, t, 1:(n_thr_chosen_index + 1)], c(1, 2), mean)
}


 

####
#### Non-diseased group:
####
print(paste("x_non_diseased_cumulative = ")) ; print(x_nd)
####
#### Diseased group:
####
print(paste("x_diseased_cumulative = ")) ; print(x_d)

###
### Prep / save final NMA simulated data outputs:
x_NMA <- list(x_nd, x_d)
x <- x_NMA

data_example_1_NMA_list <- list()
data_example_1_NMA_list$x_NMA <- x_NMA
data_example_1_NMA_list$indicator_index_test_in_study <- indicator_index_test_in_study
data_example_1_NMA_list$n_index_tests_per_study <- n_index_tests_per_study

saveRDS(object = data_example_1_NMA_list, file = file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))


x_nd[[4]]

x_d[[4]]

###
### ---- Extract specific index tests:
####

try({  
 
  
  
  
  
x_GAD2 <- x_HADS <- x_BAI <- x_GAD7 <- list()
##
GAD2_index <- 1
BAI_index  <- 2
HADS_index <- 3
GAD7_index <- 4
##
for (c in 1:2) {
  x_GAD2[[c]] <- x[[c]][[GAD2_index]]
  x_HADS[[c]] <- x[[c]][[HADS_index]]
  x_BAI[[c]]  <- x[[c]][[BAI_index]]
  x_GAD7[[c]] <- x[[c]][[GAD7_index]]
}
x_GAD7
x_GAD2
x_HADS
BAI_index

x_GAD7[[1]]


# #### ------------------- Apply ** MISSING THRESHOLDS ** to PHQ-9 data/test: -----------------------------------------------------------------------------
# # 
# # agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative
# # 
# # agg_data_cumulative_with_missing_thresholds$x_d_list
# 
# 
# 
# ## Apply PHQ-9 "missing thresholds":
# ## (1) 1st quintile (e.g., studies 1-5   if 25 total studies)  0% of studies report @ ALL threshodls, so do nothing for these.
# ## (2) 2nd quintile (e.g., studies 6-10  if 25 total studies) report 
# ## (3) 3rd quintile (e.g., studies 11-15 if 25 total studies) report ...
# ## (3) 4th quintile (e.g., studies 16-20 if 25 total studies) report ...
# ## (3) 5th quintile (e.g., studies 21-25 if 25 total studies) report ...
# 
# # n_groups <- 5
# # n_studies_div_frac <- ceiling(n_studies/n_groups)
# # n_groups * n_studies_div_frac
# # n_studies
# # last_group <- n_studies - frac*n_studies_div_frac
# # last_group
# 
# 
# n_subgroups <- 3
# ##
# study_thr_groups <- MetaOrdDTA:::divide_studies(n_studies = n_studies, n_groups = n_subgroups)
# study_thr_groups
# ##
# first_3rd  <- study_thr_groups[[1]] ## ~ 33% of studies provide ONLY ONE threshold (cases 1 and 3)
# first_3rd_split <- MetaOrdDTA:::divide_studies(n_studies = length(first_3rd), n_groups = 2)
# subset_case_1 <- first_3rd_split[[1]]
# subset_case_3 <- first_3rd_split[[2]]
# ##
# second_3rd <- study_thr_groups[[2]] ## ~ 33% of studies provide ALL thr (case 5)
# subset_case_5 <- second_3rd
# ##
# start_index <- tail(subset_case_5, 1) + 1
# third_3rd  <- study_thr_groups[[3]] ## ~ 33% of studies provide SOME thr (between 2 and < n_thr, i.e. cases 2 and 4)
# third_3rd_split <- MetaOrdDTA:::divide_studies(n_studies = length(third_3rd), n_groups = 2)
# subset_case_2 <- third_3rd_split[[1]] + start_index - 1
# subset_case_4 <- third_3rd_split[[2]] + start_index - 1
# ##
# set.seed(seed)
# ##
# {
#     ## When a study ONLY reports at the single cut-off of 10
#     ## (i.e. the "standard" PHQ-9 screening cut-off):
#     case_1 <- list()
#     ##
#     case_1$studies <- subset_case_1
#     ##
#     case_1$thr_combo_vec_1 <- c(10)
#     ##
#     print(case_1)
# }
# {
#     ## This is based off of the thresholds reported in Levis et al, 2019 
#     ## (BMJ PHQ-9 paper):
#     case_2 <- list()
#     ##
#     case_2$studies <- subset_case_2
#     ##
#     case_2$thr_combo_vec_1 <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
#     ##
#     print(case_2)
# }
# {
#     ## When a study reports at a SINGLE cut-off, cut NOT the standard cut-off of 10
#     ## (can be anything between 8 and 15):
#     case_3 <- list()
#     ##
#     case_3$studies <- subset_case_3
#     ##
#     random_single_thr <- sample(8:15, 1)
#     case_3$thr_combo_vec_1 <- c(random_single_thr)
#     ##
#     print(case_3)
# }
# {
#     case_4 <- list()
#     ##
#     case_4$studies <- subset_case_4
#     ##
#     case_4$thr_combo_vec_1  <- c(8, 9, 10, 11, 12)
#     case_4$thr_combo_vec_2  <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
#     case_4$thr_combo_vec_3  <- c(9, 10, 11, 12, 13)
#     case_4$thr_combo_vec_4  <- c(10, 11, 12, 13, 14)
#     case_4$thr_combo_vec_5  <- c(9, 10, 11, 12, 13, 14, 15)
#     case_4$thr_combo_vec_6  <- c(6, 7, 8)
#     case_4$thr_combo_vec_7  <- c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
#     case_4$thr_combo_vec_8  <- c(11, 12, 13)
#     case_4$thr_combo_vec_9  <- c(5, 6, 7, 8, 9, 10, 11, 12)
#     case_4$thr_combo_vec_10 <- c(4, 5, 6, 7, 8, 9, 10)
#     ## e.g. a study did 10 for screening and 15 for DIAGNOSIS (not common w/ PHQ-9 
#     ## but some studies do this or investigate it):
#     case_4$thr_combo_vec_11 <- c(10, 15) 
#     ## this is also based on the study that looked at PHQ-9 for both screening (i.e., using the
#     ## standard/lower thresholds) AND DIAGNOSIS (higher thresholds):
#     case_4$thr_combo_vec_12 <- c(9, 10, 11, 14, 15, 16) 
#     ##
#     print(case_4)
# }
# {
#     case_5 <- list()
#     ##
#     case_5$studies <- subset_case_5
#     ##
#     case_5$thr_combo_vec_1 <- c(1:27) ## i.e. these studies report at ALL thresholds
#     ##
#     print(case_5)
# }
# 
# # Combine all cases into a single list
# all_cases <- list(
#   case_1 = case_1,
#   case_2 = case_2,
#   case_3 = case_3,
#   case_4 = case_4,
#   case_5 = case_5
# )
# 
# 
# x_PHQ9
# 
# x_PHQ9_w_missing_thr <-  R_fn_apply_missingness_pattern_PHQ_9( x_PHQ = x_PHQ9,
#                                                               case_list = all_cases)
# 
# x_PHQ9_w_missing_thr
#  
# # x_PHQ9_w_missing_thr_old <- x_PHQ9_w_missing_thr
# # x_PHQ9_w_missing_thr[[1]] - x_PHQ9_w_missing_thr_old[[1]]
# 
# ## Compute the % of missing thr data:
# 100 * sum(x_PHQ9_w_missing_thr[[1]] == -1) / (n_studies * ncol(x_PHQ9_w_missing_thr[[1]]))
# 
# save_list_obj <- list()
# x_NMA <- x
# save_list_obj$indicator_index_test_in_study <- indicator_index_test_in_study
# save_list_obj$n_index_tests_per_study <- n_index_tests_per_study
# save_list_obj$x_NMA <- x_NMA
# save_list_obj$x_PHQ9 <- x_PHQ9
# save_list_obj$x_PHQ9_w_missing_thr <- x_PHQ9_w_missing_thr
# saveRDS(object = save_list_obj,
#         file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))
# 



})




#### ------------------- Apply ** MISSING THRESHOLDS ** to GAD-2 data/test: -----------------------------------------------------------------------------



source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
# source(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "make_Klaus_et_al_data.R")) ## to get "x_gad2" dataset
##
x_GAD2_REAL <- readRDS(file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_x_GAD2.RDS"))
##
x_GAD2_w_missing_thr <- apply_missingness_pattern( x_complete = x_GAD2,
                                                   x_pattern = x_GAD2_REAL,
                                                   enforce_consecutive_missingness = FALSE)
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
x_HADS_w_missing_thr <- apply_missingness_pattern( x_complete = x_HADS,
                                                   x_pattern = x_HADS_REAL,
                                                   enforce_consecutive_missingness = FALSE)
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
x_BAI_w_missing_thr <- apply_missingness_pattern(  x_complete = x_BAI,
                                                   x_pattern = x_BAI_REAL,
                                                   enforce_consecutive_missingness = FALSE)
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
x_GAD7_w_missing_thr <- apply_missingness_pattern( x_complete = x_GAD7,
                                                   x_pattern  = x_GAD7_REAL,
                                                   enforce_consecutive_missingness = FALSE)
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
  x_NMA_w_missing_thr[[c]][[2]] <- x_BAI_w_missing_thr[[c]]
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
  x_NMA_w_missing_thr[[c]][[4]] <- x_GAD7_w_missing_thr[[c]]
}

data_example_1_NMA_list <- list()
data_example_1_NMA_list$x_NMA <- x_NMA
##
data_example_1_NMA_list$x_NMA_w_missing_thr <- x_NMA_w_missing_thr
##
data_example_1_NMA_list$indicator_index_test_in_study <- indicator_index_test_in_study
data_example_1_NMA_list$n_index_tests_per_study <- n_index_tests_per_study

saveRDS(object = data_example_1_NMA_list, file = file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))


















