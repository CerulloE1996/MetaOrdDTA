
 
 
#  
# n_studies = n_studies
# N_per_study_mean = N_per_study_mean
# N_per_study_SD = N_per_study_SD
# assume_perfect_GS = assume_perfect_GS
# # ##
# # seed = 123
# # ##
# # true_Mean_prev = 0.20
# # true_SD_probit_prev = 0.25
# # ##
# # bivariate_locations_d_bs_het_SD = 0.50
# # bivariate_locations_nd_bs_het_SD = 0.25
# # ##
# # threshsolds_latent_between_study_SD = bs_het_C
# 
# 
# 
# t = 5
# mean_thresholds = Mean_of_thr_for_all_tests_array_d[t, 1:n_thr]
# sd_thresholds = SD_of_thr_for_all_tests_array_d[t, 1:n_thr]


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
 
 
# n_studies <- 25
# N_per_study_mean <- 2500
# N_per_study_SD <- 250
# assume_perfect_GS <- 1
# seed <- 123


#' @keywords internal
#' @export
R_fn_sim_data_ord_MA_PHQ9 <- function(   seed = 123,
                                        ##
                                        n_studies = 50,
                                        ##
                                        N_per_study_mean = 2500,
                                        N_per_study_SD   = 500,
                                        ##
                                        assume_perfect_GS = 1,
                                        ##
                                        true_Mean_prev = 0.20,
                                        true_SD_probit_prev = 0.25,
                                        ##
                                        HSROC_locations_bs_het_SD = 0.50,
                                        HSROC_raw_scales_bs_het_SD = 0.25,
                                        ##
                                        bivariate_locations_nd_bs_het_SD = 0.25,
                                        bivariate_locations_d_bs_het_SD  = 0.50,
                                        ##
                                        bs_het_C_nd,
                                        bs_het_C_d
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
     # true_mean_PHQ_9_Cerullo_Gat_params_list <- list() #### readRDS(file_path)
     # true_mean_PHQ_9_Cerullo_Gat_params_list$beta_mu <- +1.67
     # true_mean_PHQ_9_Cerullo_Gat_params_list$raw_scale_mu <- +0.061
     # 
     # true_mean_PHQ_9_Cerullo_Gat_params_list$C <- c(-2.2000, -2.0000, -1.5300, -1.1600, -0.8590,
     #                                                -0.5780 ,-0.3300, -0.1020,  0.0908,  0.2640 ,
     #                                                0.4590,   0.6600,  0.8240,  0.9440,  1.0200 ,
     #                                                1.2500,   1.3300,  1.5200,  1.6700,  1.9800 ,
     #                                                2.1100,   2.2000,  2.4500,  2.5400,  2.8000  ,
     #                                                3.0500,   3.0800)
     
     ##
     ## Store + save "true" parameter values to use for simulation:
     ##
     {
       true_mean_PHQ_9_Cerullo_Gat_params_list <- list()
       ##
       # true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC <- list()
       # true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate <- list()
       # ##
       # ## Params for ord-HSROC DGP:
       # ##
       # true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$beta_mu <-  0.837
       # true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$raw_scale_mu  <- -0.0397 
       # true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$C <-  c(-5.057555, -2.012125, -1.54322, -1.161175, -0.8611485, 
       #                                                           -0.578003, -0.327348, -0.09665075, 0.09790845, 0.272485, 0.4688355, 
       #                                                           0.6720595, 0.8375195, 0.958285, 1.031275, 1.265085, 1.347275, 
       #                                                           1.54239, 1.6897, 2.010015, 2.132755, 2.23139, 2.48997, 2.57877, 
       #                                                           2.74645, 3.00994, 3.01703)
       ##
       ## Params for ord-bivariate DGP:
       ##
       true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$beta_mu <-  c(-1.72,  -0.159)
       ##
       true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd <- c(-2.84831821, -2.3991021125, -2.03067032625, -1.7397899325, 
                                                                       -1.466111744125, -1.229449821875, -1.011687069, -0.831343162875, 
                                                                       -0.669146864435, -0.473289126897625, -0.299543294298375, -0.1337316803615, 
                                                                       -0.01839557547548, 0.078408769757035, 0.344478924181375, 0.399385793905875, 
                                                                       0.609455893075, 0.843181055, 1.052238411, 1.237164372875, 1.25896475775, 
                                                                       1.28145171775, 1.30471162275, 1.83039288375, 1.96358992625, 2.139608805, 
                                                                       2.43692291875)
       ##
       true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d  <- c(-3.32756737125, -2.42484669625, -2.16246129, -1.99281247125, 
                                                                       -1.76690405, -1.419330548, -1.139900163625, -0.889595967125, 
                                                                       -0.6783053727825, -0.514200418774488, -0.235775588161387, -0.0867989636942887, 
                                                                       0.0385369718472375, 0.065950386820325, 0.248325709728404, 0.35965754351225, 
                                                                       0.538172728131, 0.63804592431375, 0.994727939625, 1.100288967375, 
                                                                       1.22190518575, 1.54635508675, 1.673320262875, 1.839718557375, 
                                                                       2.09151451375, 2.12236033875, 2.15612568625)
     }
     ##
     bs_model <- "ord_bivariate"
     ##
     if (bs_model == "ord_HSROC") {
             
               # ## ----------------------
               # true_beta_mu_PHQ_9_Cerullo_Gat      <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$beta_mu
               # print(paste("true_beta_mu_PHQ_9_Cerullo_Gat = "))
               # print(true_beta_mu_PHQ_9_Cerullo_Gat)
               # ##
               # true_raw_scale_mu_PHQ_9_Cerullo_Gat <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$raw_scale_mu
               # print(paste("true_raw_scale_mu_PHQ_9_Cerullo_Gat = "))
               # print(true_raw_scale_mu_PHQ_9_Cerullo_Gat)
               # ##
               # true_C_mu_PHQ_9_Cerullo_Gat         <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$C
               # print(paste("true_C_mu_PHQ_9_Cerullo_Gat = "))
               # print(true_C_mu_PHQ_9_Cerullo_Gat)
               ## HSROC_locations <- c(rep(1.0, 4), true_beta_mu_PHQ_9_Cerullo_Gat)
               ##
               # location_nd[5] <- -1.0*true_beta_mu_PHQ_9_Cerullo_Gat ; location_nd[5]
               # location_d[5]  <- +1.0*true_beta_mu_PHQ_9_Cerullo_Gat ;  location_d[5]
               # ##
               # scale_nd[5] <- exp(-0.5*true_raw_scale_mu_PHQ_9_Cerullo_Gat);
               # scale_d[5]  <- exp(+0.5*true_raw_scale_mu_PHQ_9_Cerullo_Gat);
               ##
               ## get the "raw" scale params (if using R&G HSROC-based param.):
               ##
               # raw_scale_nd <- log(scale_nd)*(-2.0)
               # raw_scale_d  <- log(scale_d)*(+2.0)
               # raw_scale_d
               ##
               HSROC_raw_scales <- c(0.00, -0.025, -0.05, -0.075, true_raw_scale_mu_PHQ_9_Cerullo_Gat)
               raw_scale_nd_MEDIAN <- -1.0*HSROC_raw_scales ## since median of log-normal is equal to exp(mu), and mean is exp(mu + 0.5*sd)
               raw_scale_d_MEDIAN  <- +1.0*HSROC_raw_scales ## since median of log-normal is equal to exp(mu), and mean is exp(mu + 0.5*sd)
               # scale_nd_MEDIAN <- exp(raw_scale_nd_MEDIAN)
               # scale_d_MEDIAN  <- exp(raw_scale_d_MEDIAN)
               # ## check:
               # scale_nd_MEDIAN
               # scale_d_MEDIAN
               
     } else { 
       
               ## ----------------------
               true_beta_PHQ_9 <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$beta_mu
               print(paste("true_beta_PHQ_9 = "))
               print(true_beta_PHQ_9)
               ##
               true_C_nd_PHQ_9  <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd
               print(paste("true_C_nd_PHQ_9 = "))
               print(true_C_nd_PHQ_9)
               ##
               true_C_d_PHQ_9  <- true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d
               print(paste("true_C_d_PHQ_9 = "))
               print(true_C_d_PHQ_9)
               ##
               bivariate_locations_nd <- c(rep(-1.0, 4), true_beta_PHQ_9[1])
               bivariate_locations_d  <- c(rep(+1.0, 4), true_beta_PHQ_9[2])
       
     }
     location_nd[5] <- true_beta_PHQ_9[1]
     location_d[5]  <- true_beta_PHQ_9[2]
     ##
     ## ---- Threshold info (test 1 is binary so only 1 threshold @0 on the "latent" scale):
     ##
     n_thr_per_test <- c(1, 5, 10, 25, 27)
     max_threshold_across_all_tests <- max(n_thr_per_test)
     threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
     ##
     ## Mean_of_thr_for_all_tests_array
     Mean_of_thr_for_all_tests_array_nd <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
     Mean_of_thr_for_all_tests_array_d  <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
     ##
     ## ---- SD_of_thr_for_all_tests_array <- array(threshsolds_latent_between_study_SD, dim = c(n_tests, max_threshold_across_all_tests)) ## between-study heterogenity for thresholds 
     SD_of_thr_for_all_tests_array_nd <- array(bs_het_C_nd, dim = c(n_tests, max_threshold_across_all_tests))
     SD_of_thr_for_all_tests_array_d  <- array(bs_het_C_d,  dim = c(n_tests, max_threshold_across_all_tests))
     ##
     ## Set the "true values" of the mean threshold between studies for each test (on the * latent * scale):
     ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! 
     ## (i.e. the TRUE threshold parameters are "homogenous between classes")
     ## Binary reference test:
     Mean_thresholds_for_test_1 <- rep(NA, n_thr_per_test[1])
     Mean_thresholds_for_test_1 <- c(0.0)
     Mean_of_thr_for_all_tests_array_nd[1, 1:n_thr_per_test[1]] <- Mean_thresholds_for_test_1
     Mean_of_thr_for_all_tests_array_d[1, 1:n_thr_per_test[1]] <- Mean_of_thr_for_all_tests_array_nd[1, 1:n_thr_per_test[1]]
     ## Test 2:
     Mean_thresholds_for_test_2 <- rep(NA, n_thr_per_test[2])
     Mean_thresholds_for_test_2 <- c(-1.75, -1.00, -0.50, 0.25, 1.25)
     length(Mean_thresholds_for_test_2)
     Mean_of_thr_for_all_tests_array_nd[2, 1:n_thr_per_test[2]] <- Mean_thresholds_for_test_2
     Mean_of_thr_for_all_tests_array_d[2, 1:n_thr_per_test[2]] <- Mean_of_thr_for_all_tests_array_nd[2, 1:n_thr_per_test[2]]
     ## Test 3:
     Mean_thresholds_for_test_3 <- rep(NA, n_thr_per_test[3])
     Mean_thresholds_for_test_3 <- c(-1.75, -1.50, -1.00, -0.80, -0.50, 0.25, 1.00, 1.50, 1.80, 2.20)
     length(Mean_thresholds_for_test_3)
     Mean_of_thr_for_all_tests_array_nd[3, 1:n_thr_per_test[3]] <- Mean_thresholds_for_test_3
     Mean_of_thr_for_all_tests_array_d[3, 1:n_thr_per_test[3]] <- Mean_of_thr_for_all_tests_array_nd[3, 1:n_thr_per_test[3]]
     ## Test 4:
     Mean_thresholds_for_test_4 <- rep(NA, n_thr_per_test[4])
     Mean_thresholds_for_test_4 <- c(-2.50, -2.20, -2.00, -1.90, -1.75, 
                                     -1.50, -1.35, -1.00, -0.80, -0.50, 
                                     -0.25, -0.10, +0.00, +0.10, +0.25, 
                                     +0.40, +0.80, +1.00, +1.50, +1.80, 
                                     +2.20, +2.40, +2.50, +2.60, +2.70)
     length(Mean_thresholds_for_test_4)
     Mean_of_thr_for_all_tests_array_nd[4, 1:n_thr_per_test[4]] <- Mean_thresholds_for_test_4
     Mean_of_thr_for_all_tests_array_d[4, 1:n_thr_per_test[4]] <- Mean_of_thr_for_all_tests_array_nd[4, 1:n_thr_per_test[4]]
     ##
     ## Test 5 (BASED ON REAL_LIFE PHQ-9 DATA):
     ##
     # Mean_thresholds_for_test_5 <- true_C_mu_PHQ_9_Cerullo_Gat ### c(-2.2, true_C_mu_PHQ_9_Cerullo_Gat)
     # Mean_of_thr_for_all_tests_array[5, 1:n_thr_per_test[5]] <- Mean_thresholds_for_test_5
     ##
     Mean_thresholds_for_test_5_nd <- true_C_nd_PHQ_9
     Mean_thresholds_for_test_5_d  <- true_C_d_PHQ_9
     ##
     Mean_of_thr_for_all_tests_array_nd[5, 1:n_thr_per_test[5]] <- Mean_thresholds_for_test_5_nd
     Mean_of_thr_for_all_tests_array_d[5, 1:n_thr_per_test[5]]  <- Mean_thresholds_for_test_5_d
     
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
     # ##
     # location_nd_mu <- location_nd
     # location_d_mu  <- location_d
    
     s <- 1
    
     ## Compute true values of Se and Sp for the DGP (these will NOT vary between identical simulations with different seeds!):
     true_DGM_Se <- true_DGM_Fp <- true_DGM_Sp <- list()
     ##
     for (t in 2:n_tests) {
           true_DGM_Se[[t - 1]] <- 1.0 - pnorm(Mean_of_thr_for_all_tests_array_d[t,  1:n_thr_per_test[t]] - location_d[t])
           true_DGM_Fp[[t - 1]] <- 1.0 - pnorm(Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] - location_nd[t])
           true_DGM_Sp[[t - 1]] <- 1.0 - true_DGM_Fp[[t - 1]]
     }
     # ##
     # if (bs_model == "ord_HSROC") {
     # 
     #       true_DGM_Se[[5]] <- 1 - pnorm((true_C_mu_PHQ_9_Cerullo_Gat - location_d[5])/exp(raw_scale_d_MEDIAN[5]))
     #       true_DGM_Fp[[5]] <- 1 - pnorm((true_C_mu_PHQ_9_Cerullo_Gat - location_nd[5])/exp(raw_scale_nd_MEDIAN[5]))
     #       true_DGM_Sp[[5]] <- 1 - true_DGM_Fp[[5]]
     #         
     # } else { 
     #   
     #       true_DGM_Se[[5]] <- 1 - pnorm(true_C_d_PHQ_9  - true_beta_PHQ_9[2])
     #       true_DGM_Fp[[5]] <- 1 - pnorm(true_C_nd_PHQ_9 - true_beta_PHQ_9[1])
     #       true_DGM_Sp[[5]] <- 1 - true_DGM_Fp[[5]]
     #   
     # }
     
     
   for (s in 1:n_studies) {
             
             N <- N_per_study_vec[s]
             true_prev <- true_prev_per_study[s]
             ##
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
             
             ## Generate "true" thresholds for study s:
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             ##
             ## thr_for_all_tests_for_current_study_array ----
             thr_for_all_tests_for_current_study_nd <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             thr_for_all_tests_for_current_study_d  <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             thr_for_all_tests_for_current_study_per_n <- array(NA, dim = c(N, n_tests, max_threshold_across_all_tests))
             ##
             # ##
             for (t in 2:n_tests) {
               n_thr <- n_thr_per_test[t]
               
               if (n_thr > 0) {
                 # Generate thresholds for non-diseased group
                 thr_for_all_tests_for_current_study_nd[t, 1:n_thr] <- generate_ordered_thresholds(
                   mean_thresholds = Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr],
                   sd_thresholds = SD_of_thr_for_all_tests_array_nd[t, 1:n_thr]
                 )
                 
                 # Generate thresholds for diseased group
                 thr_for_all_tests_for_current_study_d[t, 1:n_thr] <- generate_ordered_thresholds(
                   mean_thresholds = Mean_of_thr_for_all_tests_array_d[t, 1:n_thr],
                   sd_thresholds = SD_of_thr_for_all_tests_array_d[t, 1:n_thr]
                 )
               }
               
               # Assign to individual subjects based on disease status
               # for (n in 1:N) {
               #   if (d_ind[n] == 1) {
               #     thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr] <- thr_for_all_tests_for_current_study_d[t, 1:n_thr]
               #   } else {
               #     thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr]
               #   }
               # }
               # 
               # Get indices of diseased and non-diseased subjects
               diseased_indices <- which(d_ind == 1)
               nondiseased_indices <- which(d_ind == 0)
               
               # Assign thresholds to diseased subjects (all at once)
               for (i in diseased_indices) {
                 thr_for_all_tests_for_current_study_per_n[i, t, 1:n_thr] <- thr_for_all_tests_for_current_study_d[t, 1:n_thr]
               }
               
               # Assign thresholds to non-diseased subjects (all at once)
               for (i in nondiseased_indices) {
                 thr_for_all_tests_for_current_study_per_n[i, t, 1:n_thr] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr]
               }
               
             }
             ##
             if (bs_model == "ord_bivariate") { 
               
                         ##
                         ## ---- Draw study-specific locations:
                         ##
                         location_nd_study_s <- rnorm(n = n_tests, mean = bivariate_locations_nd,  sd = bivariate_locations_nd_bs_het_SD)
                         location_d_study_s  <- rnorm(n = n_tests, mean = bivariate_locations_d,   sd = bivariate_locations_d_bs_het_SD)
                         ####
                         ##
                         ## ---- Draw study-specific RAW scale params:
                         ##
                         scale_nd_study_s <- rep(1.0, n_tests)
                         scale_d_study_s  <- rep(1.0, n_tests)
                         ##
               
             } else if (bs_model == "ord_HSROC") { 
               
                         ##
                         ## ---- Draw study-specific locations:
                         ##
                         HSROC_location_study_s <- rnorm(n = n_tests, mean = HSROC_locations, sd = HSROC_locations_bs_het_SD)
                         location_nd_study_s <- -1.0*HSROC_location_study_s
                         location_d_study_s  <- +1.0*HSROC_location_study_s
                         ##
                         ## ---- Draw study-specific RAW scale params:
                         ##
                         HSROC_raw_scale_study_s <- rnorm(n = n_tests, mean = HSROC_raw_scales, sd = HSROC_raw_scales_bs_het_SD)
                         scale_nd_study_s <- exp(-1.0*HSROC_raw_scale_study_s)
                         scale_d_study_s  <- exp(+1.0*HSROC_raw_scale_study_s)
               
             }
             ##
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
             ## Simulate the underlying "latent" continuous test responses which are distributed multivariate normal:
             latent_cts_results_neg <- LaplacesDemon::rmvn(n = n_neg, mu = location_nd_study_s, Sigma = Sigma_nd)
             latent_cts_results_pos <- LaplacesDemon::rmvn(n = n_pos, mu = location_d_study_s,  Sigma = Sigma_d)
             
             ## Logistic:
             for (t in 1:n_tests) {
                location_nd_study_s[t] <- rlogis(n = n_neg, location = location_nd_study_s[t], scale = scale_nd_study_s[t])
                location_d_study_s[t]  <- rlogis(n = n_pos, location = location_d_study_s[t],  scale = scale_d_study_s[t])
             }
             latent_results <- rbind(latent_cts_results_neg, latent_cts_results_pos)
             # latent_results <- latent_results*1.702 ## approx logistic adjustment
             
             ## Now simulate the BINARY data for the reference test (i.e. the "imperfect gold standard" - test # 1):
             results_pos <- array(NA, dim = c(n_pos, n_tests))
             results_neg <- array(NA, dim = c(n_neg, n_tests))
             y <-     array(NA, dim = c(n_neg + n_pos, n_tests))
             ## BINARY reference test results:
             y[, 1] <- ifelse(latent_results[, 1] > 0.0, 1, 0)
             ##
             ## Now simulate the ORDINAL data for tests 2-5:
             ##
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             ## Hence we can just use the "latent_results" array directly (no need to split into +'ve and -'ve)
             ## Negatives / Sp:
             for (t in 2:n_tests) {
               n_thr <- n_thr_per_test[t]
               
               # First category (values <= first threshold)
               threshold_lower <- rep(-9999, N)
               threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, 1]
               
               # Vectorized comparison for first category
               first_category <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
               y[first_category, t] <- 1
               
               # Middle categories (between consecutive thresholds)
               for (k in 2:n_thr) {
                 threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, k - 1]
                 threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, k]
                 
                 # Vectorized comparison for category k
                 category_k <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
                 y[category_k, t] <- k
               }
               
               # Last category (values > last threshold)
               threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, n_thr]
               threshold_upper <- rep(9999, N)
               
               # Vectorized comparison for last category
               last_category <- (latent_results[1:N, t] > threshold_lower) & (latent_results[1:N, t] <= threshold_upper)
               y[last_category, t] <- n_thr + 1
             }
             ##
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
                 df_pos_y_m1 <- df_pos$y - 1
                 n_TP_at_current_threshold <- sum(df_pos_y_m1[, t] >= k)
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
                 df_neg_y_m1 <- df_neg$y - 1
                 n_FP_at_current_threshold <- sum(df_neg_y_m1[, t] >= k)
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
     ##
     # thresholds_for_all_tests_for_current_study_array_nd_list = thresholds_for_all_tests_for_current_study_array_nd_list,
     # thresholds_for_all_tests_for_current_study_array_d_list = thresholds_for_all_tests_for_current_study_array_d_list,
     # thresholds_for_all_tests_for_current_study_array_list = thresholds_for_all_tests_for_current_study_array_list,
     ##
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
 
 
 

 
 
 
 
 
 





 








  
  
  
  
