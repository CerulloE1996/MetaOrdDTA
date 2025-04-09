
 
 
 
 
 # 
 # n_studies <- 10
 # N_per_study_mean <- 500
 # N_per_study_SD <- 1000
 # assume_perfect_GS <- 1
 # seed <- 123
 # true_Mean_prev <- 0.20



## -| ------------------ R function to simulate a meta-analysis dataset for (binary + ordinal) LC_MVP  --------------------------------------------
 
 
# n_studies <- 25
# N_per_study_mean <- 2500
# N_per_study_SD <- 250
# assume_perfect_GS <- 1
# seed <- 123


#' R_fn_sim_data_ordinal_MA
#' @keywords internal
#' @export
R_fn_sim_data_ordinal_MA <- function(   seed = 123,
                                        n_studies = 30,
                                        N_per_study_mean = 2000,
                                        N_per_study_SD   = 500,
                                        assume_perfect_GS = 1,
                                        true_Mean_prev = 0.20,
                                        true_SD_probit_prev = 0.25,
                                        location_nd_between_study_SD = 0.25,
                                        location_d_between_study_SD  = 0.50,
                                        raw_scale_nd_between_study_SD =  0.125,
                                        raw_scale_d_between_study_SD  =  0.250,
                                        threshsolds_latent_between_study_SD = 0.10
) {
   

     
     n_binary_tests <- 1 ## Single binary "reference" test - an imperfect gold standard (GS).
     n_ordinal_tests <- 4 # fixed to 4
     n_tests <- n_binary_tests + n_ordinal_tests
     
     ## set the seed (keep this OUTSIDE the N loop):
     set.seed(seed, kind = "L'Ecuyer-CMRG") 
     
     ## Make storage lists:
     y_ord_and_bin_list <- list()
     Sigma_nd_true_observed_list <- Sigma_d_true_observed_list <- list()
     prev_true_observed_list <- Se_true_observed_list <- Sp_true_observed_list <- list()
     Phi_Se_observed_list <- Phi_Fp_observed_list <- list()
     true_correlations_observed_vec_list <- observed_table_probs_list <- true_estimates_observed_list <- observed_cell_counts_list <- list()
     
     ii_dataset <- 0 ## counter
     
     # N_per_study_SD <- 250 ##  N_per_study_mean / 2
     N_per_study_vec <- round(TruncatedNormal::rtnorm(n = n_studies, mu = N_per_study_mean, sd = N_per_study_SD, lb = 100), 0)
  
     # ## True values for locations in both groups:
     location_nd <- rep(-1.00, n_tests)
     location_d  <- rep(+1.00, n_tests)
     
     ## For ref test:
     if (assume_perfect_GS == 1) { 
         location_nd[1] <- -100
         location_d[1]  <- +100
     }
     
     ## Set true values for disease prevelance:
     true_Mean_probit_prev <- qnorm(true_Mean_prev)
     true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev, sd = true_SD_probit_prev)
     true_prev_per_study <- pnorm(true_probit_prev_per_study)
     
     # ## Set true calues for the mean scales:
     # scale_nd <- c(1.00, 1.10, 1.20, 1.30, 1.40) ##  rep(0.25, n_tests)
     # scale_d  <- c(1.00, 0.90, 0.80, 0.70, 0.60) ## rep(0.50, n_tests)
     
     ##
     ## For PHQ-9-based data (test #5):
     ##
     ## ----------------------
     #### file_path <- file.path(MetaOrdinal_admin_path_examples, "true_mean_PHQ_9_Cerullo_Gat_params_list.RDS")
     true_mean_PHQ_9_Cerullo_Gat_params_list <- list() #### readRDS(file_path)
     true_mean_PHQ_9_Cerullo_Gat_params_list$beta_mu <- +1.67
     true_mean_PHQ_9_Cerullo_Gat_params_list$raw_scale_mu <- +0.061
     
     true_mean_PHQ_9_Cerullo_Gat_params_list$C <- c(-2.2000, -2.0000, -1.5300, -1.1600, -0.8590,
                                                    -0.5780 ,-0.3300, -0.1020,  0.0908,  0.2640 ,
                                                    0.4590,   0.6600,  0.8240,  0.9440,  1.0200 ,
                                                    1.2500,   1.3300,  1.5200,  1.6700,  1.9800 ,
                                                    2.1100,   2.2000,  2.4500,  2.5400,  2.8000  ,
                                                    3.0500,   3.0800)
     
     ## ----------------------
     true_beta_mu_PHQ_9_Cerullo_Gat      <- true_mean_PHQ_9_Cerullo_Gat_params_list$beta_mu
     print(paste("true_beta_mu_PHQ_9_Cerullo_Gat = "))
     print(true_beta_mu_PHQ_9_Cerullo_Gat)
     ##
     true_raw_scale_mu_PHQ_9_Cerullo_Gat <- true_mean_PHQ_9_Cerullo_Gat_params_list$raw_scale_mu
     print(paste("true_raw_scale_mu_PHQ_9_Cerullo_Gat = "))
     print(true_raw_scale_mu_PHQ_9_Cerullo_Gat)
     ##
     true_C_mu_PHQ_9_Cerullo_Gat         <- true_mean_PHQ_9_Cerullo_Gat_params_list$C
     print(paste("true_C_mu_PHQ_9_Cerullo_Gat = "))
     print(true_C_mu_PHQ_9_Cerullo_Gat)
     ##
     location_nd[5] <- -0.5*true_beta_mu_PHQ_9_Cerullo_Gat ; location_nd[5]
     location_d[5]  <- +0.5*true_beta_mu_PHQ_9_Cerullo_Gat ;  location_d[5]
     # ##
     # scale_nd[5] <- exp(-0.5*true_raw_scale_mu_PHQ_9_Cerullo_Gat);
     # scale_d[5]  <- exp(+0.5*true_raw_scale_mu_PHQ_9_Cerullo_Gat);
 
     
     ##
     ## get the "raw" scale params (using R&G HSROC-based param.):
     ##
     # raw_scale_nd <- log(scale_nd)*(-2.0)
     # raw_scale_d  <- log(scale_d)*(+2.0)
     # raw_scale_d
     ##
     raw_RG_scale <- c(0.00, -0.05, -0.10, -0.15, true_raw_scale_mu_PHQ_9_Cerullo_Gat)
     raw_scale_nd_MEDIAN <- -0.5*raw_RG_scale ## since median of log-normal is equal to exp(mu), and mean is exp(mu + 0.5*sd)
     raw_scale_d_MEDIAN  <- +0.5*raw_RG_scale ## since median of log-normal is equal to exp(mu), and mean is exp(mu + 0.5*sd)
     # scale_nd_MEDIAN <- exp(raw_scale_nd_MEDIAN)
     # scale_d_MEDIAN  <- exp(raw_scale_d_MEDIAN)
     # ## check:
     # scale_nd_MEDIAN
     # scale_d_MEDIAN
     
     
     ##
     ## Threshold info (test 1 is binary so only 1 threshold @0 on the "latent" scale):
     ##
     n_thr_per_test <- c(1, 5, 10, 25, 27)
     max_threshold_across_all_tests <- max(n_thr_per_test)
     threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
     Mean_of_thr_for_all_tests_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
     SD_of_thr_for_all_tests_array <- array(threshsolds_latent_between_study_SD, dim = c(n_tests, max_threshold_across_all_tests)) ## between-study heterogenity for thresholds 
     
     ## Set the "true values" of the mean threshold between studies for each test (on the * latent * scale):
     ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! 
     ## (i.e. the TRUE threshold parameters are "homogenous between classes")
     ## Binary reference test:
     Mean_thresholds_for_test_1 <- rep(NA, n_thr_per_test[1])
     Mean_thresholds_for_test_1 <- c(0.0)
     Mean_of_thr_for_all_tests_array[1, 1:n_thr_per_test[1]] <- Mean_thresholds_for_test_1
     ## Test 2:
     Mean_thresholds_for_test_2 <- rep(NA, n_thr_per_test[2])
     Mean_thresholds_for_test_2 <- c(-1.75, -1.00, -0.50, 0.25, 1.25)
     length(Mean_thresholds_for_test_2)
     Mean_of_thr_for_all_tests_array[2, 1:n_thr_per_test[2]] <- Mean_thresholds_for_test_2
     ## Test 3:
     Mean_thresholds_for_test_3 <- rep(NA, n_thr_per_test[3])
     Mean_thresholds_for_test_3 <- c(-1.75, -1.50, -1.00, -0.80, -0.50, 0.25, 1.00, 1.50, 1.80, 2.20)
     length(Mean_thresholds_for_test_3)
     Mean_of_thr_for_all_tests_array[3, 1:n_thr_per_test[3]] <- Mean_thresholds_for_test_3
     ## Test 4:
     Mean_thresholds_for_test_4 <- rep(NA, n_thr_per_test[4])
     Mean_thresholds_for_test_4 <- c(-2.50, -2.20, -2.00, -1.90, -1.75, 
                                     -1.50, -1.35, -1.00, -0.80, -0.50, 
                                     -0.25, -0.10, +0.00, +0.10, +0.25, 
                                     +0.40, +0.80, +1.00, +1.50, +1.80, 
                                     +2.20, +2.40, +2.50, +2.60, +2.70)
     length(Mean_thresholds_for_test_4)
     Mean_of_thr_for_all_tests_array[4, 1:n_thr_per_test[4]] <- Mean_thresholds_for_test_4
     ##
     ## Test 5 (BASED ON REAL_LIFE PHQ-9 DATA):
     ##
     Mean_thresholds_for_test_5 <- true_C_mu_PHQ_9_Cerullo_Gat ### c(-2.2, true_C_mu_PHQ_9_Cerullo_Gat)
     Mean_of_thr_for_all_tests_array[5, 1:n_thr_per_test[5]] <- Mean_thresholds_for_test_5
     

     
     Se_for_current_study_at_threshold_0_list <- Sp_for_current_study_at_threshold_0_list <- Fp_for_current_study_at_threshold_0_list <- list()
     Se_per_study_all_tests_all_thresholds_list <- Sp_per_study_all_tests_all_thresholds_list <- list()
     prev_true_observed_list <- list()
     Se_per_study_ref <- Sp_per_study_ref <- list()
     thresholds_for_all_tests_for_current_study_array_list <- list()
     y_list <- list()
     y_df_list <- list()
     ##
     n_TP_at_current_threshold_OVERALL <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     n_FP_at_current_threshold_OVERALL <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     n_pos_OVERALL <- 0 
     n_neg_OVERALL <- 0 
     n_total_OVERALL <- 0
     Se_OVERALL_all_tests_all_thresholds <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     Sp_OVERALL_all_tests_all_thresholds <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     ##
     # true_DGP_one_m_Se <- pnorm(location_d - Mean_of_thr_for_all_tests_array)
     # true_DGP_Sp       <- pnorm(location_nd - Mean_of_thr_for_all_tests_array)
     ##
     location_nd_mu <- location_nd
     location_d_mu  <- location_d
    
    s <- 1
    
    ## Compute true values of Se and Sp for the DGP (these will NOT vary between identical simulations with different seeds!):
    true_DGM_Se <- true_DGM_Fp <- true_DGM_Sp <- list()
    true_DGM_Se[[5]] <- 1 - pnorm((true_C_mu_PHQ_9_Cerullo_Gat - location_d[5])/exp(raw_scale_d_MEDIAN[5]))
    true_DGM_Fp[[5]] <- 1 - pnorm((true_C_mu_PHQ_9_Cerullo_Gat - location_nd[5])/exp(raw_scale_nd_MEDIAN[5]))
    true_DGM_Sp[[5]] <- 1 - true_DGM_Fp[[5]]
    

    
    
   for (s in 1:n_studies) {
             
             N <- N_per_study_vec[s]
             true_prev <- true_prev_per_study[s]
             
             ## Generate "true" thresholds for study s:
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             thr_for_all_tests_for_current_study_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               thr_for_all_tests_for_current_study_array[t, 1:n_thr_per_test[t]] <- rnorm( n = n_thr_per_test[t], 
                                                                                           mean = Mean_of_thr_for_all_tests_array[t, 1:n_thr_per_test[t]], 
                                                                                           sd = SD_of_thr_for_all_tests_array[t, 1:n_thr_per_test[t]])
             }
             
             thresholds_for_all_tests_for_current_study_array_list[[s]] <- thr_for_all_tests_for_current_study_array
             
             ##
             ## Draw study-specific locations:
             ##
             location_nd_study_s <- rnorm(n = n_tests, mean = location_nd_mu,  sd = location_nd_between_study_SD)
             location_d_study_s  <- rnorm(n = n_tests, mean = location_d_mu,   sd = location_d_between_study_SD)
             ##
             ## Draw study-specific RAW scale params::
             ##
             raw_scale_nd_study_s <- rnorm(n = n_tests, mean = raw_scale_nd_MEDIAN,  sd = raw_scale_nd_between_study_SD)
             raw_scale_d_study_s  <- rnorm(n = n_tests, mean = raw_scale_d_MEDIAN,   sd = raw_scale_d_between_study_SD)
             ##
             scale_nd_study_s <- exp(-0.5*raw_scale_nd_study_s)
             scale_d_study_s  <- exp(+0.5*raw_scale_d_study_s)
             
             
             if (assume_perfect_GS == TRUE) {
                 rho1 <- 1.00
             } else { 
                 rho1 <- 0.20
             }
             ##
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
                   Omega_nd <-  0.5 * Omega_highly_varied ## Corr(D-) is HALF that of the Corr(D+)
                   diag(Omega_nd) <- rep(1, n_tests)
             }
             ##
             Omega_nd <- as.matrix(Matrix::nearPD(Omega_nd, keepDiag = TRUE)$mat)
             Omega_d  <- as.matrix(Matrix::nearPD(Omega_d, keepDiag = TRUE)$mat)
             ##
             L_Omega_nd   <- t(chol(Omega_nd)) # PD check (fails if not PD)
             L_Omega_d    <- t(chol(Omega_d))   ## BayesMVP:::Rcpp_Chol(Sigma_nd) # PD check (fails if not PD)
             ##
             Sigma_nd <- diag(scale_nd_study_s) %*% Omega_nd %*% diag(scale_nd_study_s)
             Sigma_d  <- diag(scale_d_study_s)  %*% Omega_d  %*% diag(scale_d_study_s)
             ##
             Sigma_nd <- as.matrix(Matrix::nearPD(Sigma_nd)$mat)
             Sigma_d  <- as.matrix(Matrix::nearPD(Sigma_d)$mat)
             ##
             Sigma_nd <- as.matrix(Matrix::forceSymmetric(Sigma_nd))
             Sigma_d  <- as.matrix(Matrix::forceSymmetric(Sigma_d))
             ##
             L_Sigma_nd   <- t(chol(Sigma_nd)) # PD check (fails if not PD)
             L_Sigma_d    <- t(chol(Sigma_d))   ## BayesMVP:::Rcpp_Chol(Sigma_nd) # PD check (fails if not PD)
             ##
             d_ind <- sort(rbinom(n = N,
                                  size = 1, 
                                  prob = true_prev))
             ##
             n_pos <- sum(d_ind)
             n_pos_OVERALL <- n_pos_OVERALL + n_pos
             ##
             n_neg <- N - sum(d_ind)
             n_neg_OVERALL <- n_neg_OVERALL + n_neg
             ##
             n_total <- n_pos + n_neg
             n_total_OVERALL <- n_total_OVERALL + n_total
             ##
             ## Simulate the underlying "latent" continuous test responses which are distributed multivariate normal:
             latent_cts_results_neg <- LaplacesDemon::rmvn(n = n_neg, mu = location_nd_study_s, Sigma = Sigma_nd)
             latent_cts_results_pos <- LaplacesDemon::rmvn(n = n_pos, mu = location_d_study_s,  Sigma = Sigma_d)
             latent_results <- rbind(latent_cts_results_neg, latent_cts_results_pos)
             
             ## Now simulate the BINARY data for the reference test (i.e. the "imperfect gold standard" - test # 1):
             results_pos <- array(NA, dim = c(n_pos, n_tests))
             results_neg <- array(NA, dim = c(n_neg, n_tests))
             y <-     array(NA, dim = c(n_neg + n_pos, n_tests))
             ## BINARY reference test results:
             y[, 1] <- ifelse(latent_results[, 1] > 0.0, 1, 0)
             
             ## Now simulate the ORDINAL data for tests 2-5:
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             ## Hence we can just use the "latent_results" array directly (no need to split into +'ve and -'ve)
             ## Negatives / Sp:
             for (t in 2:n_tests) {
               
                       n_thr <- n_thr_per_test[t] 
                       
                       threshold_lower <- -9999
                       threshold_upper <- head(thr_for_all_tests_for_current_study_array[t, 1:n_thr], 1)
                       y[, t] <- ifelse(  (latent_results[, t] > threshold_lower) & (latent_results[, t] <= threshold_upper),
                                                1, 
                                                y[, t])
                    
                       for (k in 2:(n_thr)) {
                           threshold_lower <- thr_for_all_tests_for_current_study_array[t, k - 1] ; threshold_lower
                           threshold_upper <- thr_for_all_tests_for_current_study_array[t, k]     ; threshold_upper
                           y[, t] <- ifelse(  (latent_results[, t] > threshold_lower) & (latent_results[, t] <= threshold_upper), 
                                                     k, 
                                                     y[, t])
                       }
                       
                       threshold_lower <- tail(thr_for_all_tests_for_current_study_array[t, 1:n_thr], 1)  
                       threshold_upper <- +9999
                       y[, t] <- ifelse(     (latent_results[, t] > threshold_lower) & (latent_results[, t] < threshold_upper),
                                                   n_thr + 1, 
                                                   y[, t])
                 
             }
             y
             
             df <- dplyr::tibble(y, latent_results, d_ind)
             df_pos <- dplyr::filter(df, d_ind == 1)
             df_neg <- dplyr::filter(df, d_ind == 0)

             # prev:
             prev_true_observed <-  print(round(sum(d_ind)/N, 3))
             prev_true_observed_list[[s]] <- prev_true_observed

             ## Compute OVERALL observed TRUE Se:
             Se_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               n_thr <- n_thr_per_test[t] 
               for (k in 1:n_thr) {
                 n_TP_at_current_threshold <- sum(df_pos$y[, t] >= k + 1)
                 Se_current_study_all_tests_all_thresholds[t, k] <- n_TP_at_current_threshold/n_pos
                 ## Increment comulative counters:
                 n_TP_at_current_threshold_OVERALL[t, k] <- n_TP_at_current_threshold_OVERALL[t, k] + n_TP_at_current_threshold
               }
             }
             Se_per_study_all_tests_all_thresholds_list[[s]] <- Se_current_study_all_tests_all_thresholds
             
             ## Compute OVERALL observed TRUE Fp / Sp:
             Sp_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               n_thr <- n_thr_per_test[t] 
               for (k in 1:n_thr) {
                 n_FP_at_current_threshold <- sum(df_neg$y[, t] >= k + 1)
                 Fp_at_threshold_k <- n_FP_at_current_threshold/n_neg
                 Sp_at_threshold_k <- 1.0 - Fp_at_threshold_k
                 Sp_current_study_all_tests_all_thresholds[t, k] <- Sp_at_threshold_k
                 ## Increment comulative counters:
                 n_FP_at_current_threshold_OVERALL[t, k] <- n_FP_at_current_threshold_OVERALL[t, k] + n_FP_at_current_threshold
               }
             }
             Sp_per_study_all_tests_all_thresholds_list[[s]] <- Sp_current_study_all_tests_all_thresholds
             
             ## For reference test:
             t <- 1
             ## Se:
             Phi_Se_ref <- qnorm(sum(df_pos$results[ ,t])/nrow(df_pos))
             Se_ref <- pnorm(Phi_Se_ref)
             Se_per_study_ref[[s]] <- Se_ref
             
             ## Sp:
             Phi_Fp_ref <-  qnorm( 1.0 - ((nrow(df_neg) - sum(df_neg$results[ ,t]))/nrow(df_neg))  )
             Fp_ref <- pnorm(Phi_Fp_ref)
             Sp_ref <- 1.0 - Fp_ref
             Sp_per_study_ref[[s]] <- Sp_ref

             y_list[[s]] <- y
             y_df_list[[s]] <- data.frame(y)
             
     
     
   }
   

       ## Compute OVERALL observed TRUE Se:
       for (t in 2:n_tests) {
         n_thr <- n_thr_per_test[t]
         for (k in 1:n_thr) {
           Se_OVERALL_all_tests_all_thresholds[t, k] <- n_TP_at_current_threshold_OVERALL[t, k]/n_pos_OVERALL
         }
       }
       # Se_OVERALL_all_tests_all_thresholds <- Reduce(`+`, Se_per_study_all_tests_all_thresholds_list) / length(Se_per_study_all_tests_all_thresholds_list)
      
       
       ## Compute OVERALL observed TRUE Fp / Sp:
       for (t in 2:n_tests) {
         n_thr <- n_thr_per_test[t]
         for (k in 1:n_thr) {
           Sp_OVERALL_all_tests_all_thresholds[t, k] <- 1.0 - (n_FP_at_current_threshold_OVERALL[t, k]/n_neg_OVERALL)
         }
       }
       # Sp_OVERALL_all_tests_all_thresholds <- Reduce(`+`, Sp_per_study_all_tests_all_thresholds_list) / length(Sp_per_study_all_tests_all_thresholds_list)
       
       y_tibble <- NULL
       y_tibble <- tibble(data.table::rbindlist(y_df_list, idcol = "Study"))
   
   return(list(
     N_per_study_vec = N_per_study_vec,
     y_list = y_list,
     y_df_list = y_df_list,
     y_tibble = y_tibble,
     ## Between-study true params:
     # Mean_Se_at_threshold_0 = Mean_Se_at_threshold_0,
     # Mean_Sp_at_threshold_0 = Mean_Sp_at_threshold_0,
     true_Mean_prev = true_Mean_prev,
     ## Between-study true params:
     # SD_of_Phi_Se_at_threshold_0 = SD_of_Phi_Se_at_threshold_0,
     # SD_of_Phi_Fp_at_threshold_0 = SD_of_Phi_Fp_at_threshold_0,
     ## Between-study true params:
     n_thr_per_test = n_thr_per_test,
     n_thr_per_test_excl_ref = n_thr_per_test[-c(1)],
     ## Within-study true params:
     Se_for_current_study_at_threshold_0_list = Se_for_current_study_at_threshold_0_list,
     Sp_for_current_study_at_threshold_0_list = Sp_for_current_study_at_threshold_0_list,
     thresholds_for_all_tests_for_current_study_array_list = thresholds_for_all_tests_for_current_study_array_list,
     Se_per_study_all_tests_all_thresholds_list = Se_per_study_all_tests_all_thresholds_list,
     Sp_per_study_all_tests_all_thresholds_list = Sp_per_study_all_tests_all_thresholds_list,
     prev_true_observed_list = prev_true_observed_list,
     Se_per_study_ref = Se_per_study_ref,
     Sp_per_study_ref = Sp_per_study_ref,
     Se_OVERALL_all_tests_all_thresholds = Se_OVERALL_all_tests_all_thresholds,
     Sp_OVERALL_all_tests_all_thresholds = Sp_OVERALL_all_tests_all_thresholds,
     ##
     true_DGM_Se = true_DGM_Se,
     true_DGM_Fp = true_DGM_Fp,
     true_DGM_Sp = true_DGM_Sp
     # true_DGP_one_m_Se = true_DGP_one_m_Se,
     # true_DGP_Sp = true_DGP_Sp
     # ## 
     # Sigma_nd_true_observed_list = Sigma_nd_true_observed_list,
     # Sigma_d_true_observed_list = Sigma_d_true_observed_list,
     # ##
     # Phi_Se_observed_list = Phi_Se_observed_list,
     # Phi_Fp_observed_list = Phi_Fp_observed_list,
     # ##
     # prev_true_observed_list = prev_true_observed_list,
     # Se_true_observed_list = Se_true_observed_list,
     # Sp_true_observed_list = Sp_true_observed_list,
     # ##
     # true_correlations_observed_vec_list = true_correlations_observed_vec_list,
     # observed_table_probs_list = observed_table_probs_list,
     # true_estimates_observed_list = true_estimates_observed_list,
     # observed_cell_counts_list = observed_cell_counts_list
   ))
   
   
 }
 
 
 

 
 
 
 
 
 





 








  
  
  
  
