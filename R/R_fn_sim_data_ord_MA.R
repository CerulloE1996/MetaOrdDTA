
 
 



generate_ordered_thresholds <- function(mean_thresholds, 
                                        sd_thresholds) {
  
          n_thr <- length(mean_thresholds)
          
          if (n_thr <= 0) {
            return(numeric(0))
          }
          
          # Create thresholds that maintain the same pattern as the means
          ordered_thresholds <- numeric(n_thr)
          
          # Generate correlated noise
          # Start with initial random offset
          offset <- rnorm(1, 0, sd_thresholds[1] * 0.5)
          
          # Apply offset to first threshold
          ordered_thresholds[1] <- mean_thresholds[1] + offset
          
          # For each subsequent threshold
          if (n_thr > 1) {
            
            for (i in 2:n_thr) {
              # Calculate the original spacing between consecutive thresholds
              orig_diff <- mean_thresholds[i] - mean_thresholds[i-1]
              
              # Add some noise to the spacing, but ensure it stays positive
              noise_sd <- min(sd_thresholds[i], sd_thresholds[i-1]) * 0.5
              diff_noise <- max(0.01, rnorm(1, orig_diff, noise_sd))
              
              # Set the threshold based on previous threshold plus noisy difference
              ordered_thresholds[i] <- ordered_thresholds[i-1] + diff_noise
            }
            
          }
          
          ordered_thresholds
          
          return(ordered_thresholds)
  
}




 
 

## -| ------------------ R function to simulate a meta-analysis dataset for (binary + ordinal) LC_MVP  --------------------------------------------
 
 
# seed = 123
# n_studies = 50
# N_per_study_mean = 2500
# N_per_study_SD = 500
# assume_perfect_GS = 1
# true_Mean_prev = 0.2
# true_SD_probit_prev = 0.25
# HSROC_locations_bs_het_SD = 0.5
# HSROC_raw_scales_bs_het_SD = 0.25
# bivariate_locations_nd_bs_het_SD = 0.25
# bivariate_locations_d_bs_het_SD = 0.5
# bs_het_C_nd_prob_scale
# bs_het_C_d_prob_scale
# scale_ALL_bs_SDs_by = 0.75


#' R_fn_sim_data_ord_MA
#' @keywords internal
#' @export
R_fn_sim_data_ord_MA <- function(  seed = 123, 
                                   n_studies = 50, 
                                   N_per_study_mean = 2500, 
                                   N_per_study_SD = 500, 
                                   assume_perfect_GS = 1, 
                                   true_Mean_prev = 0.2, 
                                   true_SD_probit_prev = 0.25, 
                                   # HSROC_locations_bs_het_SD = 0.5, 
                                   # HSROC_raw_scales_bs_het_SD = 0.25, 
                                   bivariate_locations_nd_bs_het_SD = 0.25, 
                                   bivariate_locations_d_bs_het_SD = 0.5, 
                                   bs_het_C_nd_prob_scale, 
                                   bs_het_C_d_prob_scale,
                                   scale_ALL_bs_SDs_by = 1
) {
  
  
          bivariate_locations_nd_bs_het_SD <- scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD
          bivariate_locations_d_bs_het_SD  <- scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD
          
          # Setup parameters
          n_binary_tests <- 1
          n_ordinal_tests <- 4
          n_tests <- n_binary_tests + n_ordinal_tests
          set.seed(seed, kind = "L'Ecuyer-CMRG")
          
          # Initialize lists and arrays
          y_ord_and_bin_list <- list()
          Sigma_nd_true_observed_list <- Sigma_d_true_observed_list <- list()
          prev_true_observed_list <- Se_true_observed_list <- Sp_true_observed_list <- list()
          Phi_Se_observed_list <- Phi_Fp_observed_list <- list()
          true_correlations_observed_vec_list <- observed_table_probs_list <- true_estimates_observed_list <- observed_cell_counts_list <- list()
          ii_dataset <- 0
          
          # Generate study sizes
          N_per_study_vec <- round(TruncatedNormal::rtnorm(n = n_studies, mu = N_per_study_mean, sd = N_per_study_SD, lb = 100), 0)
          
          # Initialize location parameters
          bivariate_locations_nd <- rep(-1, n_tests)
          bivariate_locations_d <- rep(+1, n_tests)
          if (assume_perfect_GS == 1) {
            bivariate_locations_nd[1] <- -100
            bivariate_locations_d[1] <- +100
          }
          
          # Generate prevalence parameters
          true_Mean_probit_prev <- qnorm(true_Mean_prev)
          true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev, 
                                              sd = true_SD_probit_prev)
          true_prev_per_study <- pnorm(true_probit_prev_per_study)
          
          # Set up test parameters
          {
            ##
            ## ---- PHQ-9 parameters (from single-study REAL PHQ-9 dataset from Baron et al).
            ##
            true_mean_PHQ_9_Cerullo_Gat_params_list <- list()
            ##
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$beta_mu <- c(-1.17, -0.118)
            ##
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd <- c(-2.2994594075, -1.85076751625, -1.4823266075, -1.191850410875, 
                                                                            -0.9185224348125, -0.6824948028125, -0.465187581006875, -0.28524418401225, 
                                                                            -0.123824854394419, 0.0714019907404125, 0.244086194109869, 0.408806356951388, 
                                                                            0.5236473534, 0.61907566584375, 0.880868595, 0.93506378875, 1.13891996525, 
                                                                            1.362710241125, 1.5555927000625, 1.71732335, 1.73495572625, 1.753785303125, 
                                                                            1.7730980725, 2.15678496, 2.221392155, 2.293717775, 2.38682278375)
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd_SD_prob_scale <- bs_het_C_nd_prob_scale
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd_SDs <- SD_approx_ID_ord_prob_to_C_probit(
                                                                                         mean_C_cutpoint_scale = true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd,
                                                                                         SD_prob_scale =  true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd_SD_prob_scale)
            ##
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d <- c(-2.917814866875, -2.309792485, -2.073766918125, -1.915478330625, 
                                                                           -1.70158212875, -1.364530803125, -1.089750625625, -0.8431635128125, 
                                                                           -0.63443191095, -0.471749624497031, -0.195740285686575, -0.047765369834075, 
                                                                           0.07715920840186, 0.104641682855775, 0.285555863388069, 0.396088297391706, 
                                                                           0.572971156206875, 0.6713609324625, 1.022715328875, 1.126140787875, 
                                                                           1.2449444535, 1.5585936949375, 1.6775888, 1.83199706625, 2.053982441875, 
                                                                           2.078204196875, 2.10569573375)
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d_SD_prob_scale <- bs_het_C_d_prob_scale
            true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d_SDs <- SD_approx_ID_ord_prob_to_C_probit(
                                                                                          mean_C_cutpoint_scale = true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d,
                                                                                          SD_prob_scale =  true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d_SD_prob_scale)
            ##
            ## ---- GAD-2 parameters:
            ## (from fitting Xu-Rand-Kappa model to REAL GAD-2 MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_GAD_2_params_list <- list()
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$beta_mu <- c(-0.550, 0.184)
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-0.8394009, -0.0902575, 0.5646964, 0.9078822, 1.3215716, 1.6462413)
            true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale <- bs_het_C_nd_prob_scale
            true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SDs <- SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale =  true_mean_GAD_2_params_list$Ord_bivariate$C_nd,
                                                                                                   SD_prob_scale = true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale)
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$C_d <- c(-1.3322775, -0.7746882, -0.0929791, 0.3271562, 0.8693018, 1.3237361)
            true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale <- bs_het_C_d_prob_scale
            true_mean_GAD_2_params_list$Ord_bivariate$C_d_SDs <- SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale =  true_mean_GAD_2_params_list$Ord_bivariate$C_d,
                                                                                                   SD_prob_scale = true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale)
          }
          
          # Set model type and parameters
          bs_model <- "ord_bivariate"
          if (bs_model == "ord_HSROC") {
            # HSROC model code would go here
          } else {
            # Get parameters for bivariate model
            true_beta_PHQ_9 <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$beta_mu
            print(paste("true_beta_PHQ_9 = "))
            print(true_beta_PHQ_9)
            true_C_nd_PHQ_9 <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd
            print(paste("true_C_nd_PHQ_9 = "))
            print(true_C_nd_PHQ_9)
            true_C_d_PHQ_9 <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d
            print(paste("true_C_d_PHQ_9 = "))
            print(true_C_d_PHQ_9)
            true_beta_GAD_2 <- true_mean_GAD_2_params_list$Ord_bivariate$beta_mu
            print(paste("true_beta_GAD_2 = "))
            print(true_beta_GAD_2)
            true_C_nd_GAD_2 <- true_mean_GAD_2_params_list$Ord_bivariate$C_nd
            print(paste("true_C_nd_GAD_2 = "))
            print(true_C_nd_GAD_2)
            true_C_d_GAD_2 <- true_mean_GAD_2_params_list$Ord_bivariate$C_d
            print(paste("true_C_d_GAD_2 = "))
            print(true_C_d_GAD_2)
            
            # Set bivariate locations
            bivariate_locations_nd <- c(-1, true_beta_GAD_2[1], -1, -1, true_beta_PHQ_9[1])
            bivariate_locations_d <- c(+1, true_beta_GAD_2[2], +1, +1, true_beta_PHQ_9[2])
          }
          
          # Override PHQ-9 location parameters
          bivariate_locations_nd[5] <- true_beta_PHQ_9[1]
          bivariate_locations_d[5] <- true_beta_PHQ_9[2]
          
          # Set up threshold parameters
          n_thr_per_test <- c(1, 6, 10, 25, 27)
          max_threshold_across_all_tests <- max(n_thr_per_test)
          threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_nd <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_d <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
          SD_of_thr_for_all_tests_array_nd <- array( 0.001, dim = c(n_tests, max_threshold_across_all_tests))
          SD_of_thr_for_all_tests_array_d <- array( 0.001, dim = c(n_tests, max_threshold_across_all_tests))
          ##
          ## ---- Set thresholds for test 1 (binary reference test):
          ##
          test <- 1
          Mean_thresholds_for_test_1 <- rep(NA, n_thr_per_test[1])
          Mean_thresholds_for_test_1 <- c(0)
          Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]] <- Mean_thresholds_for_test_1
          Mean_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]] <- Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]
          ##
          ## ---- Set thresholds for test 2 (GAD-2):
          ##
          test <- 2
          Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]] <- true_mean_GAD_2_params_list$Ord_bivariate$C_nd
          SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]   <- true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SDs
          ##
          Mean_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]  <- true_mean_GAD_2_params_list$Ord_bivariate$C_d
          SD_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]    <- true_mean_GAD_2_params_list$Ord_bivariate$C_d_SDs
          ##
          ## ---- Set thresholds for test 3:
          ##
          test <- 3
          Mean_thresholds_for_test_3 <- rep(NA, n_thr_per_test[3])
          Mean_thresholds_for_test_3 <- c(-1.75, -1.5, -1, -0.8, -0.5, 0.25, 1, 1.5, 1.8, 2.2)
          length(Mean_thresholds_for_test_3)
          Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]] <- Mean_thresholds_for_test_3
          SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]   <- SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = Mean_thresholds_for_test_3,
                                                                                                                SD_prob_scale = bs_het_C_d_prob_scale)
          ##
          Mean_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]  <- Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]
          SD_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]   <- SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]
          ##
          ## ---- Set thresholds for test 4:
          ##
          test <- 4
          Mean_thresholds_for_test_4 <- rep(NA, n_thr_per_test[4])
          Mean_thresholds_for_test_4 <- c(
            -2.5, -2.2, -2, -1.9, -1.75, -1.5, -1.35, -1, -0.8, -0.5, -0.25, -0.1, 
            +0, +0.1, +0.25, +0.4, +0.8, +1, +1.5, +1.8, +2.2, +2.4, +2.5, +2.6, +2.7
          )
          length(Mean_thresholds_for_test_4)
          Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]] <- Mean_thresholds_for_test_4
          SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]   <- SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale = Mean_thresholds_for_test_4,
                                                                                                                SD_prob_scale = bs_het_C_d_prob_scale)
          ##
          Mean_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]  <- Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]
          SD_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]   <- SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]
          ##
          ## ---- Set thresholds for test 5 (PHQ-9):
          ##
          test <- 5
          Mean_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]] <- true_C_nd_PHQ_9
          SD_of_thr_for_all_tests_array_nd[test, 1:n_thr_per_test[test]]   <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd_SDs
          ##
          Mean_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]  <- true_C_d_PHQ_9
          SD_of_thr_for_all_tests_array_d[test, 1:n_thr_per_test[test]]    <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d_SDs
          ##
          ## -----------------------------------------------
          ##
          if (bs_het_C_nd_prob_scale < 0.02) {
            SD_of_thr_for_all_tests_array_nd[,] <- 0.001
          }
          if (bs_het_C_d_prob_scale < 0.02) {
            SD_of_thr_for_all_tests_array_d[,] <- 0.001
          }
          ##
          SD_of_thr_for_all_tests_array_nd[,] <- scale_ALL_bs_SDs_by * SD_of_thr_for_all_tests_array_nd[,] 
          SD_of_thr_for_all_tests_array_d[,]  <- scale_ALL_bs_SDs_by * SD_of_thr_for_all_tests_array_d[,] 
          ##
          # Initialize more lists and arrays
          Se_for_current_study_at_threshold_0_list <- Sp_for_current_study_at_threshold_0_list <- Fp_for_current_study_at_threshold_0_list <- list()
          Se_per_study_all_tests_all_thresholds_list <- Sp_per_study_all_tests_all_thresholds_list <- list()
          prev_true_observed_list <- list()
          Se_per_study_ref <- Sp_per_study_ref <- list()
          thresholds_for_all_tests_for_current_study_array_list <- list()
          y_list <- list()
          y_df_list <- list()
          
          # Initialize arrays for threshold calculations
          n_TP_at_each_threshold_OVERALL <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          n_FP_at_each_threshold_OVERALL <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          n_TP_at_each_threshold <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          n_FP_at_each_threshold <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Se_all_tests_all_thresholds <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Sp_all_tests_all_thresholds <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Fp_all_tests_all_thresholds <- array(0, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          
          # Initialize overall counters
          n_pos_OVERALL <- 0
          n_neg_OVERALL <- 0
          n_total_OVERALL <- 0
          Se_OVERALL_all_tests_all_thresholds <- array(0, dim = c(n_tests, max_threshold_across_all_tests + 1))
          Sp_OVERALL_all_tests_all_thresholds <- array(0, dim = c(n_tests, max_threshold_across_all_tests + 1))
          
          # Initialize study counter
          s <- 1
          
          # Calculate true DGM values
          true_DGM_Se <- true_DGM_Fp <- true_DGM_Sp <- list()
          for (t in 2:n_tests) {
            true_DGM_Se[[t - 1]] <- 1 - pnorm(Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] - bivariate_locations_d[t])
            true_DGM_Fp[[t - 1]] <- 1 - pnorm(Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] - bivariate_locations_nd[t])
            true_DGM_Sp[[t - 1]] <- 1 - true_DGM_Fp[[t - 1]]
          }
          
          ##
          ## ---- Between-study correlation (and var-cov) matrix (using same one for each test)
          ##
          if (bs_model == "ord_bivariate") {
              ##
              Omega_bs <- matrix(1.0, nrow = 2, ncol = 2)
              Sigma_bs <- matrix(NA, nrow = 2, ncol = 2)
              ##
              bs_corr <- 0.122 ## this bs-corr value was obtained from GAD-2 (Klaus et al) data fitting the Xu-Rand model to it.
              Omega_bs[1, 2] <- bs_corr
              Omega_bs[2, 1] <- bs_corr
              ##
              Sigma_bs[1, 1] <- bivariate_locations_nd_bs_het_SD^2
              Sigma_bs[2, 2] <- bivariate_locations_d_bs_het_SD^2
              Sigma_bs[1, 2] <- bs_corr * bivariate_locations_nd_bs_het_SD * bivariate_locations_d_bs_het_SD
              Sigma_bs[2, 1] <- Sigma_bs[1, 2]
              ##
          }
          # Loop through studies
          for (s in 1:n_studies) {
            
                    # Set up study parameters
                    N <- N_per_study_vec[s]
                    ##
                    true_prev <- true_prev_per_study[s]
                    d_ind <- sort(rbinom(n = N, size = 1, prob = true_prev))
                    ##
                    n_pos <- sum(d_ind)
                    n_pos_OVERALL <- n_pos_OVERALL + n_pos
                    ##
                    n_neg <- N - sum(d_ind)
                    n_neg_OVERALL <- n_neg_OVERALL + n_neg
                    ##
                    n_total <- n_pos + n_neg
                    n_total_OVERALL <- n_total_OVERALL + n_total
                    
                    # Initialize threshold arrays for current study
                    thr_for_all_tests_for_current_study_nd <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                    thr_for_all_tests_for_current_study_d <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                    thr_for_all_tests_for_current_study_per_n <- array(NA, dim = c(N, n_tests, max_threshold_across_all_tests))
                    
                    # Generate thresholds for all tests
                    for (t in 2:n_tests) {
                      n_thr <- n_thr_per_test[t]
                      if (n_thr > 0) {
                        thr_for_all_tests_for_current_study_nd[t, 1:n_thr] <- generate_ordered_thresholds(
                          mean_thresholds = Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr], 
                          sd_thresholds = SD_of_thr_for_all_tests_array_nd[t, 1:n_thr]
                        )
                        thr_for_all_tests_for_current_study_d[t, 1:n_thr] <- generate_ordered_thresholds(
                          mean_thresholds = Mean_of_thr_for_all_tests_array_d[t, 1:n_thr], 
                          sd_thresholds = SD_of_thr_for_all_tests_array_d[t, 1:n_thr]
                        )
                      }
                      
                      # Assign thresholds to diseased and non-diseased subjects
                      diseased_indices <- which(d_ind == 1)
                      nondiseased_indices <- which(d_ind == 0)
                      
                      for (i in diseased_indices) {
                        thr_for_all_tests_for_current_study_per_n[i, t, 1:n_thr] <- thr_for_all_tests_for_current_study_d[t, 1:n_thr]
                      }
                      for (i in nondiseased_indices) {
                        thr_for_all_tests_for_current_study_per_n[i, t, 1:n_thr] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr]
                      }
                    }
                    
                    # Generate study-specific parameters
                    if (bs_model == "ord_bivariate") {
                      ##
                      # location_nd_study_s <- rnorm(n = n_tests, mean = bivariate_locations_nd, sd = bivariate_locations_nd_bs_het_SD)
                      # location_d_study_s  <- rnorm(n = n_tests, mean = bivariate_locations_d, sd = bivariate_locations_d_bs_het_SD)
                      ##
                      ##
                      location_nd_study_s <- c()
                      location_d_study_s <- c()
                      ##
                      for (t in 1:n_tests) {
                        locations <- LaplacesDemon::rmvn(n = 1, mu = c(bivariate_locations_nd[t], bivariate_locations_d[t]), Sigma = Sigma_bs)
                        location_nd_study_s[t] <- locations[1]
                        location_d_study_s[t]  <- locations[2]
                      }
                      ##
                      scale_nd_study_s <- rep(1, n_tests)
                      scale_d_study_s <- rep(1, n_tests)
                    } else if (bs_model == "ord_HSROC") {
                      HSROC_location_study_s <- rnorm(n = n_tests, mean = HSROC_locations, sd = HSROC_locations_bs_het_SD)
                      location_nd_study_s <- -1 * HSROC_location_study_s
                      location_d_study_s <- +1 * HSROC_location_study_s
                      HSROC_raw_scale_study_s <- rnorm(n = n_tests, mean = HSROC_raw_scales, sd = HSROC_raw_scales_bs_het_SD)
                      scale_nd_study_s <- exp(-1 * HSROC_raw_scale_study_s)
                      scale_d_study_s <- exp(+1 * HSROC_raw_scale_study_s)
                    }
                    
                    # Set correlation parameters
                    if (assume_perfect_GS == TRUE) {
                      rho1 <- 0.0
                    }
                    else {
                      rho1 <- 0.2
                    }
                    ##
                    ## Set up correlation matrices:
                    Omega_highly_varied <- matrix(c(1,     rho1,  rho1,     rho1,     rho1,
                                                    rho1,  1.00,  0.50,     0.20,     0.10,
                                                    rho1,  0.50,  1.00,     0.40,     0.40,
                                                    rho1,  0.20,  0.40,     1.00,     0.70,
                                                    rho1,  0.10,  0.40,     0.70,     1.00),
                                                  nrow = n_tests,
                                                  ncol = n_tests)
                    
                    {
                      Omega_d <- Omega_highly_varied
                      diag(Omega_d) <- rep(1, n_tests)
                      Omega_nd <- 0.5 * Omega_highly_varied
                      diag(Omega_nd) <- rep(1, n_tests)
                    }
                    
                    # Ensure correlation matrices are proper
                    Omega_nd <- as.matrix(Matrix::nearPD(Omega_nd, keepDiag = TRUE)$mat)
                    Omega_d <- as.matrix(Matrix::nearPD(Omega_d, keepDiag = TRUE)$mat)
                    L_Omega_nd <- t(chol(Omega_nd))
                    L_Omega_d <- t(chol(Omega_d))
                    
                    # Calculate covariance matrices
                    Sigma_nd <- diag(scale_nd_study_s) %*% Omega_nd %*% diag(scale_nd_study_s)
                    Sigma_d <- diag(scale_d_study_s) %*% Omega_d %*% diag(scale_d_study_s)
                    Sigma_nd <- as.matrix(Matrix::nearPD(Sigma_nd)$mat)
                    Sigma_d <- as.matrix(Matrix::nearPD(Sigma_d)$mat)
                    Sigma_nd <- as.matrix(Matrix::forceSymmetric(Sigma_nd))
                    Sigma_d <- as.matrix(Matrix::forceSymmetric(Sigma_d))
                    L_Sigma_nd <- t(chol(Sigma_nd))
                    L_Sigma_d <- t(chol(Sigma_d))
                    
                    # Generate latent continuous results
                    latent_cts_results_neg <- LaplacesDemon::rmvn(n = n_neg, mu = location_nd_study_s, Sigma = Sigma_nd)
                    latent_cts_results_pos <- LaplacesDemon::rmvn(n = n_pos, mu = location_d_study_s, Sigma = Sigma_d)
                    latent_results <- rbind(latent_cts_results_neg, latent_cts_results_pos)
                    
                    # Initialize result arrays
                    results_pos <- array(NA, dim = c(n_pos, n_tests))
                    results_neg <- array(NA, dim = c(n_neg, n_tests))
                    y <- array(NA, dim = c(n_neg + n_pos, n_tests))
                    
                    # Set binary test results (test 1)
                    y[, 1] <- ifelse(latent_results[, 1] > 0, 1, 0)
                    disease_status <- d_ind
                    print(sum(disease_status == 1))
                    print(sum(d_ind == 1))
                    
                    # Generate ordinal test results (tests 2-5)
                    for (t in 2:n_tests) {
                            n_thr <- n_thr_per_test[t]
                            n_cat <- n_thr + 1
                            
                            # First category (below first threshold)
                            threshold_lower <- rep(-9999, N)
                            threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, 1]
                            first_category <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
                            y[first_category, t] <- 0
                            
                            # Middle categories
                            for (k in 2:(n_thr)) {
                              threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, k - 1]
                              threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, k]
                              category_k <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
                              y[category_k, t] <- k - 1
                            }
                            
                            # Last category (above last threshold)
                            threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, n_thr]
                            threshold_upper <- rep(9999, N)
                            last_category <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
                            y[last_category, t] <- n_cat - 1
                    }
                    
                    # Create data frames
                    df <- dplyr::tibble(y, latent_results, d_ind)
                    df_pos <- dplyr::filter(df, d_ind == 1)
                    df_neg <- dplyr::filter(df, d_ind == 0)
                    
                    # Calculate prevalence
                    prev_true_observed <- print(round(sum(d_ind)/N, 3))
                    prev_true_observed_list[[s]] <- prev_true_observed
                    
                    # Calculate test performance for each threshold
                    for (t in 2:n_tests) {
                      n_thr <- n_thr_per_test[t]
                      n_cat <- n_thr + 1
                      
                      # Initialize counts for first threshold
                      n_TP_at_each_threshold[s, t, 1] <- n_pos
                      n_FP_at_each_threshold[s, t, 1] <- n_neg
                      
                      # Calculate counts for each threshold
                      for (k in 2:n_cat) {
                        df_pos_y <- y[, t][disease_status == 1]
                        df_neg_y <- y[, t][disease_status == 0]
                        
                        if (t == 2 && s == 1) {
                          print(paste("k =", k, "| Unique values in df_pos_y:", paste(sort(unique(df_pos_y)), collapse = ", ")))
                          print(paste("Count of df_pos_y >= (k-1):", sum(df_pos_y >= (k - 1)), "| Total positive cases:", n_pos))
                        }
                        
                        positive_count <- sum(df_pos_y >= (k - 1))
                        negative_count <- sum(df_neg_y >= (k - 1))
                        n_TP_at_each_threshold[s, t, k] <- sum(df_pos_y >= (k - 1))
                        n_FP_at_each_threshold[s, t, k] <- sum(df_neg_y >= (k - 1))
                      }
                      
                      # Calculate sensitivity and specificity
                      Se_all_tests_all_thresholds[s, t, 1:n_cat] <- n_TP_at_each_threshold[s, t, 1:n_cat]/n_pos
                      Fp_all_tests_all_thresholds[s, t, 1:n_cat] <- n_FP_at_each_threshold[s, t, 1:n_cat]/n_neg
                      Sp_all_tests_all_thresholds[s, t, 1:n_cat] <- 1 - Fp_all_tests_all_thresholds[s, t, 1:n_cat]
                    }
                    
                    # Calculate reference test performance (test 1)
                    t <- 1
                    Phi_Se_ref <- qnorm(sum(df_pos$results[, t])/nrow(df_pos))
                    Se_ref <- pnorm(Phi_Se_ref)
                    Se_per_study_ref[[s]] <- Se_ref
                    Phi_Fp_ref <- qnorm(1 - ((nrow(df_neg) - sum(df_neg$results[, t]))/nrow(df_neg)))
                    Fp_ref <- pnorm(Phi_Fp_ref)
                    Sp_ref <- 1 - Fp_ref
                    Sp_per_study_ref[[s]] <- Sp_ref
                    
                    # Store results
                    y_list[[s]] <- y
                    y_df_list[[s]] <- data.frame(y)
            
          }
          
          # Combine results into a tibble
          y_tibble <- NULL
          y_tibble <- tibble(data.table::rbindlist(y_df_list, idcol = "Study"))
          
          # Return results
          return(list(
            N_per_study_vec = N_per_study_vec, 
            y_list = y_list, 
            y_df_list = y_df_list, 
            y_tibble = y_tibble, 
            true_Mean_prev = true_Mean_prev, 
            n_thr_per_test = n_thr_per_test, 
            n_thr_per_test_excl_ref = n_thr_per_test[-c(1)], 
            n_TP_at_each_threshold = n_TP_at_each_threshold, 
            n_FP_at_each_threshold = n_FP_at_each_threshold, 
            Se_all_tests_all_thresholds = Se_all_tests_all_thresholds, 
            Sp_all_tests_all_thresholds = Sp_all_tests_all_thresholds, 
            Fp_all_tests_all_thresholds = Fp_all_tests_all_thresholds, 
            prev_true_observed_list = prev_true_observed_list, 
            Se_per_study_ref = Se_per_study_ref, 
            Sp_per_study_ref = Sp_per_study_ref, 
            true_DGM_Se = true_DGM_Se, 
            true_DGM_Fp = true_DGM_Fp, 
            true_DGM_Sp = true_DGM_Sp
          ))
  
}



 








  
  
  
  
