
 
#' SD_approx_probit_to_prob
#' @export
SD_approx_probit_to_prob <- function(SD_probit_scale) { 
  
          SD_prob_scale <- SD_probit_scale * dnorm(qnorm(SD_probit_scale))
          ##
          return(signif(SD_prob_scale, 3))
  
}


#' SD_approx_ID_ord_prob_to_C_probit
#' @export
SD_approx_ID_ord_prob_to_C_probit <- function(mean_C_cutpoint_scale, 
                                              SD_prob_scale) {
  
          SD_probit_scale <- (1.0 / dnorm(mean_C_cutpoint_scale)) * SD_prob_scale
          return(SD_probit_scale)
          
}


# Use this in the simulation
generate_ordered_thresholds_delta_method <- function(mean_cutpoints, 
                                                     SD_prob_scale_vec,
                                                     use_probit_link = TRUE) {
  
          n_thr <- length(mean_cutpoints)
          
          if (n_thr <= 0) return(numeric(0))
          
          # Transform SDs from probability to cutpoint scale
          sigma_C <- SD_approx_ID_ord_prob_to_C_probit(mean_cutpoints, SD_prob_scale_vec)
          
          # Generate cutpoints
          cutpoints <- numeric(n_thr)
          for (k in 1:n_thr) {
            cutpoints[k] <- rnorm(1, mean = mean_cutpoints[k], sd = sigma_C[k])
          }
          
          # Ensure ordering (this is the key addition for simulation)
          cutpoints <- sort(cutpoints)
          
          return(cutpoints)
  
}


#' generate_ordered_thresholds_delta_method_correlated
#' @export
generate_ordered_thresholds_delta_method_correlated <- function( mean_cutpoints, 
                                                                 SD_prob_scale_vec,
                                                                 use_probit_link = TRUE,
                                                                 base_correlation = 0.5) {
  
          n_thr <- length(mean_cutpoints)
          if (n_thr <= 0) return(numeric(0))
          ##
          ## Transform SDs from probability to cutpoint scale
          ##
          sigma_C <- SD_approx_ID_ord_prob_to_C_probit(mean_cutpoints, SD_prob_scale_vec)
          
          if (n_thr == 1) {
            return(rnorm(1, mean = mean_cutpoints[1], sd = sigma_C[1]))
          }
          ##
          # Create correlation matrix - higher correlation for adjacent cutpoints
          # This mimics the induced-Dirichlet structure
          cor_mat <- matrix(0, n_thr, n_thr)
          for (i in 1:n_thr) {
            for (j in 1:n_thr) {
              # Correlation decays with distance
              cor_mat[i,j] <- base_correlation^abs(i-j)
            }
          }
          ##
          # Covariance matrix
          cov_mat <- diag(sigma_C) %*% cor_mat %*% diag(sigma_C)
          ##
          # Generate correlated cutpoints
          cutpoints <- MASS::mvrnorm(1, mu = mean_cutpoints, Sigma = cov_mat)
          ##
          # Ensure ordering
          cutpoints <- sort(cutpoints)
          ##
          return(cutpoints)
  
}


#' generate_ordered_thresholds_induced_dirichlet
#' @export
generate_ordered_thresholds_induced_dirichlet <- function(alpha_vec, 
                                                          use_probit_link = TRUE) {
  
          n_cat <- length(alpha_vec)
          n_thr <- n_cat - 1
          
          if (n_thr <= 0) {
            return(numeric(0))
          }
          
          # Step 1: Generate ordinal probabilities from Dirichlet
          ord_probs <- as.vector(MCMCpack::rdirichlet(1, alpha_vec))
          
          # Step 2: Convert to cumulative probabilities
          cumul_probs <- cumsum(ord_probs)[1:n_thr]
          
          # Step 3: Convert to cutpoints
          if (use_probit_link) {
            cutpoints <- qnorm(cumul_probs)
          } else {
            cutpoints <- qlogis(cumul_probs)  # logit link
          }
          
          return(cutpoints)
          
}






#' generate_ordered_thresholds_induced_dirichlet_sample
#' @export
generate_ordered_thresholds_induced_dirichlet_sample <- function(n_samples = 100000,
                                                                 alpha_vec, 
                                                                 use_probit_link = TRUE) { 
  
        
          n_cat <- length(alpha_vec)
          n_thr <- n_cat - 1
          
          samples <- array(NA, dim = c(n_samples, ncol = n_thr))
          
          for (i in 1:n_samples) { 
              samples[i, ] <- generate_ordered_thresholds_induced_dirichlet(alpha_vec = alpha_vec,
                                                            use_probit_link = use_probit_link)
          }
          
          return(apply(samples, FUN = median, c(2)))
    
  
}




# Function to compute Dirichlet SDs from alpha
compute_dirichlet_SDs <- function(alpha_vec) {
  
  alpha_0 <- sum(alpha_vec)
  n_cat <- length(alpha_vec)
  
  # SD for each category probability
  sds <- sqrt((alpha_vec * (alpha_0 - alpha_vec)) / (alpha_0^2 * (alpha_0 + 1)))
  
  return(sds)
  
}


 
 

## -| ------------------ R function to simulate a meta-analysis dataset for (binary + ordinal) LC_MVP  --------------------------------------------
#  
#  
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
# bivariate_locations_d_bs_het_SD  = bivariate_locations_d_bs_het_SD
# bivariate_locations_nd_bs_het_SD = bivariate_locations_nd_bs_het_SD
# ##
# bs_het_C_nd_prob_scale = bs_het_C_nd_prob_scale
# bs_het_C_d_prob_scale = bs_het_C_d_prob_scale
# ##
# scale_ALL_bs_SDs_by = scale_ALL_bs_SDs_by





#' R_fn_sim_data_ord_MA_v2
#' @keywords internal
#' @export
R_fn_sim_data_ord_MA_v2 <- function(   seed = 123, 
                                       DGM,
                                       n_studies, 
                                       N_per_study_mean, 
                                       N_per_study_SD, 
                                       assume_perfect_GS = 1, 
                                       true_Mean_prev = 0.2, 
                                       true_SD_probit_prev = 0.25, 
                                       # # HSROC_locations_bs_het_SD = 0.5, 
                                       # # HSROC_raw_scales_bs_het_SD = 0.25, 
                                       # bivariate_locations_nd_bs_het_SD = 0.25, 
                                       # bivariate_locations_d_bs_het_SD = 0.5, 
                                       bs_het_C_nd_prob_scale, 
                                       bs_het_C_d_prob_scale,
                                       reduction_factor = 1.00
                                       # scale_ALL_bs_SDs_by = 1
) {
  
  
          # bivariate_locations_nd_bs_het_SD <- scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD
          # bivariate_locations_d_bs_het_SD  <- scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD
          
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
          
          # # Initialize location parameters
          # bivariate_locations_nd <- rep(-1, n_tests)
          # bivariate_locations_d <- rep(+1, n_tests)
          # if (assume_perfect_GS == 1) {
          #   bivariate_locations_nd[1] <- -100
          #   bivariate_locations_d[1]  <- +100
          # }
          
          # Generate prevalence parameters
          true_Mean_probit_prev <- qnorm(true_Mean_prev)
          true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev,
                                              sd = true_SD_probit_prev)
          true_prev_per_study <- pnorm(true_probit_prev_per_study)
          ##
          # reduction_factor <- 1.00
          # reduction_factor <- 0.70
          # reduction_factor <- 0.50
          # reduction_factor^2
          
          
          
          # ##
          # ## ---- GAD-2: --------------------------------------------------------------------------------------------------------------------
          # ##
          # ## ---- ord_biv_random_C:
          # ##
          # round(mean(c(true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale)), 4) ## 0.036
          # round(mean(c(true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale)),  4) ## 0.068
          # ##
          # ## ---- R&G_random_C:
          # ##
          # round(mean(c(true_mean_GAD_2_params_list$Ord_HSROC$C_SD_prob_scale)), 4) ## 0.0685
          # ##
          # ## ---- HADS: ---------------------------------------------------------------------------------------------------------------------
          # ##
          # ##
          # ## ---- ord_biv_random_C:
          # ##
          # round(mean(c(true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale)), 4) ## 0.0104
          # round(mean(c(true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale)),  4) ## 0.0108
          # ##
          # ## ---- R&G_random_C:
          # ##
          # round(mean(c(true_mean_HADS_params_list$Ord_HSROC$C_SD_prob_scale)), 4) ## 0.0095
          # ## 
          # ## ---- BAI: ----------------------------------------------------------------------------------------------------------------------
          # ##
          # ## ---- ord_biv_random_C:
          # ##
          # round(mean(c(true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale)), 4)  ## 0.0045
          # round(mean(c(true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale)),  4)  ## 0.0026
          # ##
          # ## ---- R&G_random_C:
          # ##
          # round(mean(c(true_mean_BAI_params_list$Ord_HSROC$C_SD_prob_scale)), 4) ## 0.0031
          ##
          ##
          ## Set up test parameters
          {
            ##
            ## ---- GAD-2 parameters: ----------------------------------------------------------------------------------------------------------------------
##
            ## ---- (from fitting ord_biv-Rand-Kappa model to REAL GAD-2 MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_GAD_2_params_list <- list()
            ##
            if (DGM == "ord_biv_random_C") {
              
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_mu   <- c(-0.762, 0.197)
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_mu   <- c(-0.685885, 0.936699)
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_SD   <- c(0.339, 0.561)*reduction_factor
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_SD   <- c(0.290853, 0.500444)*reduction_factor
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.176
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.3644
                  ##
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-1.002235, -0.3171695,  0.3326495, 0.6686015, 1.08083,   1.427205)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-0.929478, -0.238893, 0.412653, 0.753993, 1.18966, 1.56409)
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_d  <- c(-1.47935,  -0.8745465, -0.133177,  0.2706515, 0.7900285, 1.21108)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_d  <- c(-0.993491, -0.238206, 0.542912, 0.980165, 1.55188, 2.06426)
                  ##
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale <- c(0.0381236, 0.0424788, 0.0446779, 0.033099, 0.0326942, 0.0250071)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale <- c(0.0377096, 0.0412537, 0.0424379, 0.0311462, 0.0303931, 0.0228643)
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale  <- c(0.0555101, 0.0628668, 0.0843183, 0.0712172, 0.0745728, 0.0592809)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale  <- c(0.0408704, 0.04722, 0.0500518, 0.0367689, 0.0334796, 0.0216572)
                  ##
                  # true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd          <- c(14.20995, 19.4291, 22.71955, 
                  #                                                                  10.478, 10.03115, 5.48333, 6.93372) / (reduction_factor^2)
                  # true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd          <- c(18.5311, 23.5948, 25.6252, 11.6937, 11.0579, 5.9268, 6.4349) / (reduction_factor^2)
                  true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd          <- c(18.6967, 23.7175, 26.2555, 11.8253, 11.1197, 5.96799,
                                                                                   6.40415) / (reduction_factor^2)
                  # true_mean_GAD_2_params_list$Ord_bivariate$alpha_d           <- c(2.03021, 2.755295, 5.821855, 
                  #                                                                  3.763645, 4.173785, 2.400225, 2.800985) / (reduction_factor^2)
                  # true_mean_GAD_2_params_list$Ord_bivariate$alpha_d           <- c(13.5591, 19.9794, 24.0159, 10.5039, 8.4667, 3.31011, 1.9566) / (reduction_factor^2)
                  true_mean_GAD_2_params_list$Ord_bivariate$alpha_d           <- c(13.4193, 20.1435, 24.5148, 10.549, 8.43175, 3.29664,
                                                                                   1.93082) / (reduction_factor^2)
                
            } else if (DGM == "ord_biv_fixed_C") {
                
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_mu <- c(-0.987, 0.157)
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_mu <- c(-0.987855, 0.245392)
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_SD <- c(0.389, 0.668)*reduction_factor
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_SD <- c(0.334615, 0.5002)*reduction_factor
                  # true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.0819
                  true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.3225
                  ##
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-1.27417, -0.5334115,  0.155574,  0.493848,  0.916881,  1.275435)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-1.25265, -0.533424, 0.144763, 0.481504, 0.904145, 1.26096)
                  # true_mean_GAD_2_params_list$Ord_bivariate$C_d  <- c(-1.40413, -0.8865765, -0.1816915, 0.2231125, 0.7619075, 1.275885)
                  true_mean_GAD_2_params_list$Ord_bivariate$C_d  <- c(-1.66601, -0.93275, -0.151223, 0.278304, 0.83734, 1.36536)
                
            } else if (DGM == "Jones") {
              
                  true_mean_GAD_2_params_list$Jones$beta_mu <- c(0.279, 2.02)
                  true_mean_GAD_2_params_list$Jones$beta_SD <- c(0.415, 1.06)*reduction_factor
                  true_mean_GAD_2_params_list$Jones$beta_corr <- 0.00403
                  ##
                  # true_mean_GAD_2_params_list$Jones$cts_thr_values <- seq(from = 1, to = 21, by = 1)
                  # true_mean_GAD_2_params_list$Jones$lambda <- c(0.673, 0.936)
                  ##
                  true_mean_GAD_2_params_list$Jones$C_nd <- c(0, 0.800205,  1.38365, 1.85971, 2.26912, 2.632215)
                  true_mean_GAD_2_params_list$Jones$C_d  <- c(0, 0.8786875, 1.61281, 2.26642, 2.86603, 3.42584)
                  ##
                  true_mean_GAD_2_params_list$Jones$raw_scale_mu <- c(0.0609, 0.183)
                  true_mean_GAD_2_params_list$Jones$raw_scale_SD <- c(0.0853, 0.283)*reduction_factor
                  true_mean_GAD_2_params_list$Jones$raw_scale_corr <- -0.247
              
            } else if (DGM == "R&G_random_C") {
              
                  # true_mean_GAD_2_params_list$Ord_HSROC$beta_mu <- -0.710
                  true_mean_GAD_2_params_list$Ord_HSROC$beta_mu <- -0.772
                  # true_mean_GAD_2_params_list$Ord_HSROC$beta_SD <-  0.392*reduction_factor
                  true_mean_GAD_2_params_list$Ord_HSROC$beta_SD <-  0.2078*reduction_factor
                  ##
                  # true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_mu <-  -0.0438
                  true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_mu <-  -0.07801
                  # true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_SD <-   0.171*reduction_factor
                  true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_SD <-   0.09337*reduction_factor
                  ##
                  # true_mean_GAD_2_params_list$Ord_HSROC$C <- c(-0.989595, -0.2895845, 0.393812, 0.770447, 1.207325, 1.62195)
                  true_mean_GAD_2_params_list$Ord_HSROC$C <- c(-1.04768, -0.301202, 0.386064, 0.760449, 1.23366, 1.70676)
                  true_mean_GAD_2_params_list$Ord_HSROC$C_nd <- true_mean_GAD_2_params_list$Ord_HSROC$C
                  true_mean_GAD_2_params_list$Ord_HSROC$C_d  <- true_mean_GAD_2_params_list$Ord_HSROC$C
                  ##
                  # true_mean_GAD_2_params_list$Ord_HSROC$C_SD_prob_scale <-  c(0.0739518, 0.0813187, 0.0850478, 0.0632803, 0.0603268, 0.0471866)
                  true_mean_GAD_2_params_list$Ord_HSROC$C_SD_prob_scale <-  c(0.061047, 0.07108, 0.0744359, 0.0556163, 0.0534102, 0.0411834)
                  ##
                  # true_mean_GAD_2_params_list$Ord_HSROC$alpha <- c(4.35868, 5.585555, 6.34091, 2.985595, 2.68411, 1.56714, 1.502915) / (reduction_factor^2)
                  # true_mean_GAD_2_params_list$Ord_HSROC$alpha <- c(5.25336, 7.81922, 8.96809, 4.20731, 3.83192, 2.1587, 
                  #                                                  1.79068) / (reduction_factor^2)
                  true_mean_GAD_2_params_list$Ord_HSROC$alpha <- c(5.2323, 7.78426, 8.94391, 4.20263, 3.8345, 2.16141,
                                                                   1.79032) / (reduction_factor^2)
              
            }
            ##
            ## ---- GAD-7 parameters. -----------------------------------------------------------------------------------------------------------------------
            ##
            true_mean_GAD_7_params_list <- list()
            # ##
            # if (DGM == "ord_biv_random_C") {
                true_mean_GAD_7_params_list$Ord_bivariate$beta_mu <- c(-1.41, -0.283)
                true_mean_GAD_7_params_list$Ord_bivariate$beta_SD <- c(0.488, 0.674)*reduction_factor
                true_mean_GAD_7_params_list$Ord_bivariate$beta_corr <- 0.307
                ##
                true_mean_GAD_7_params_list$Ord_bivariate$C_nd <- c(-2.231685, -1.79927, -1.482965, -1.220955, -0.98261, -0.780394, 
                                                                    -0.5859865, -0.360594, -0.217805, -0.070975, 0.0544958, 0.176212, 
                                                                    0.283997, 0.407415, 0.5439385, 0.6842815, 0.8175955, 0.913633, 
                                                                    1.067695, 1.22144, 1.35017)
                true_mean_GAD_7_params_list$Ord_bivariate$C_nd_SD_prob_scale <-   c(0.0084992, 0.0102436, 0.012722, 0.014581, 0.0155888, 0.0165293, 
                                                                                    0.0171685, 0.0193348, 0.0161019, 0.0165393, 0.0154, 0.0152028, 
                                                                                    0.0142087, 0.01495, 0.0151354, 0.0148238, 0.0136074, 0.0112217, 
                                                                                    0.0134593, 0.0115781, 0.0109631)
                true_mean_GAD_7_params_list$Ord_bivariate$alpha_nd          <-  c(2.79575, 4.099455, 6.4108, 8.50915, 9.799435, 11.0924, 12.0361, 
                                                                                  15.5634, 10.525, 11.1519, 9.61439, 9.385075, 8.130495, 9.049675, 
                                                                                  9.303845, 8.918155, 7.438585, 5.00854, 7.273505, 5.341545, 4.76307, 
                                                                                  17.65125) / (reduction_factor^2)
                ##
                true_mean_GAD_7_params_list$Ord_bivariate$C_d <-  c(-2.209915, -2.02756, -1.88334, -1.705805, -1.51828, -1.286595, 
                                                                    -1.09908, -0.8463045, -0.6829905, -0.516839, -0.3368345, -0.150944, 
                                                                    0.0289825, 0.1866205, 0.3744105, 0.5237585, 0.667981, 0.8760615, 
                                                                    1.09799, 1.36988, 1.631645)
                true_mean_GAD_7_params_list$Ord_bivariate$C_d_SD_prob_scale <-     c(0.013356, 0.0087849, 0.0094477, 0.0116721, 0.013411, 0.0183166, 
                                                                                     0.0185152, 0.0237383, 0.0212544, 0.0226872, 0.0244091, 0.0257238, 
                                                                                     0.025128, 0.0240619, 0.0253914, 0.0219318, 0.0213946, 0.0227724, 
                                                                                     0.0224695, 0.021033, 0.018146)
                true_mean_GAD_7_params_list$Ord_bivariate$alpha_d           <-      c(1.799325, 0.767976, 0.887956, 1.371725, 1.82478, 3.448515, 
                                                                                      3.540835, 5.93037, 4.7235, 5.419735, 6.327325, 7.101355, 6.7747, 
                                                                                      6.191155, 6.943915, 5.119395, 4.84974, 5.52035, 5.36885, 4.6842, 
                                                                                      3.43053, 5.521015) / (reduction_factor^2)
            ##
            ## ---- HADS parameters: ----------------------------------------------------------------------------------------------------------------------
##
            ## (from fitting ord_biv-Rand-Kappa model to REAL HADS MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_HADS_params_list <- list()
            ##
            if (DGM == "ord_biv_random_C") {
              
                  # true_mean_HADS_params_list$Ord_bivariate$beta_mu <- c(-0.644, 0.0866)
                  true_mean_HADS_params_list$Ord_bivariate$beta_mu   <- c(-0.717607, 0.0205349)
                  # true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.233,  0.392)*reduction_factor
                  true_mean_HADS_params_list$Ord_bivariate$beta_SD   <- c(0.236112, 0.41698)*reduction_factor
                  # true_mean_HADS_params_list$Ord_bivariate$beta_corr <-  0.305
                  true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.2935
                  ##
                  # true_mean_HADS_params_list$Ord_bivariate$C_nd <-   c(-2.10529, -1.667765, -1.32778, -1.034325, -0.7666715, -0.5273575, 
                  #                                                      -0.2762745, -0.0323531, 0.2050175, 0.4217955, 0.599618, 0.7960495, 
                  #                                                      0.973136, 1.15867, 1.36864, 1.540585, 1.73854, 1.914335, 2.056615, 
                  #                                                      2.215475, 2.34737)
                  true_mean_HADS_params_list$Ord_bivariate$C_nd <- c(-2.20996, -1.75548, -1.41123, -1.11558, -0.845998, -0.605184,
                                                                     -0.353195, -0.107899, 0.12932, 0.346172, 0.523601, 0.722416, 
                                                                      0.90491, 1.09915, 1.32351, 1.50889, 1.73446, 1.94807,
                                                                      2.15298, 2.44522, 2.99679)
                  # true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale <-    c(0.0072208, 0.0090472, 0.0108992, 0.0124052, 0.0135303, 0.0141023, 
                  #                                                                     0.0152486, 0.0155346, 0.0152954, 0.0144312, 0.0127017, 0.0125942, 
                  #                                                                     0.0112547, 0.0106024, 0.010087, 0.0080056, 0.0075404, 0.0059428, 
                  #                                                                     0.0046284, 0.0041917, 0.0032816)  
                  true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale <- c(0.0059234, 0.00790024, 0.00966879, 0.0111565, 0.0123586, 0.0129809,
                                                                                   0.0141602, 0.0145371, 0.014444, 0.0137594, 0.0121893, 0.012229,
                                                                                   0.0110463, 0.0105015, 0.0100548, 0.00804148, 0.00764826, 0.00613696,
                                                                                   0.00491777, 0.00449456, 0.00369894)
                  # true_mean_HADS_params_list$Ord_bivariate$alpha_nd         <-       c(6.6585, 10.5856, 15.55235, 20.61605, 24.84645, 27.17665, 32.4214, 
                  #                                                                      33.84225, 32.6674, 28.78185, 21.8639, 21.4429, 16.7675, 14.7908, 
                  #                                                                      13.35695, 8.33584, 7.33371, 4.54828, 2.714125, 2.219215, 1.346705, 
                  #                                                                      3.62417) / (reduction_factor^2)
                  # true_mean_HADS_params_list$Ord_bivariate$alpha_nd          <- c(5.85941, 10.5468, 16.0195, 21.6411, 26.9269, 29.9459,
                  #                                                                 36.2425, 38.4251, 37.8778, 34.0143, 26.1402, 26.3234,
                  #                                                                 21.1908, 19.0457, 17.385,  10.9381, 9.86775, 6.29642,
                  #                                                                 4.02025, 3.3525, 2.2645, 0.864862) / (reduction_factor^2)
                  true_mean_HADS_params_list$Ord_bivariate$alpha_nd          <- c(5.90532, 10.6291, 16.0131, 21.5625, 26.6274, 29.6131,
                                                                                  35.6167, 38.0392, 37.0889, 33.1407, 25.5686, 25.3074,
                                                                                  20.6517, 18.3002, 16.6949, 10.4899, 9.33653, 6.03551,
                                                                                  3.90606, 3.25259, 2.24677, 0.812404) / (reduction_factor^2)
                  ##
                  # true_mean_HADS_params_list$Ord_bivariate$C_d <-      c(-2.014315, -1.936835, -1.7171, -1.584575, -1.4278, -1.23411, 
                  #                                                        -0.9858285, -0.7802815, -0.542741, -0.26083, -0.020179, 0.1654705, 
                  #                                                        0.3911145, 0.656948, 0.895045, 1.12504, 1.284525, 1.41728, 1.552285, 
                  #                                                        1.73874, 1.829675)
                  true_mean_HADS_params_list$Ord_bivariate$C_d  <- c(-2.60921, -2.44061, -2.01294, -1.79114, -1.58059, -1.3589, 
                                                                     -1.08735, -0.86826, -0.621757, -0.329962, -0.0798263, 0.111696, 
                                                                     0.350494, 0.629013, 0.890788, 1.15068, 1.34336, 1.5157,
                                                                     1.7088, 2.06551, 2.29048)
                  # true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale <-    c(0.0085058, 0.0037001, 0.0069388, 0.0063387, 0.0077475, 0.0097957, 
                  #                                                                    0.0125642, 0.0127646, 0.014842, 0.0171012, 0.0165019, 0.0146033, 
                  #                                                                    0.0155868, 0.0161474, 0.0141258, 0.012533, 0.0095831, 0.0078529, 
                  #                                                                    0.0071434, 0.0074656, 0.0045503) 
                  true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale  <- c(0.00316948, 0.00218316, 0.00519178, 0.00515816, 0.00611069, 0.00744644, 
                                                                                   0.00959548, 0.00987227, 0.0114361,  0.0132667,  0.0128782,  0.0115446,
                                                                                   0.0125575,  0.0129566,  0.0116694,  0.0104686,  0.00798533, 0.00666878, 
                                                                                   0.00616692, 0.00665284, 0.0038261)
                  # true_mean_HADS_params_list$Ord_bivariate$alpha_d           <-      c(7.12456, 1.25392, 4.960385, 4.089665, 6.108865, 9.94842, 16.71395, 
                  #                                                                      16.8367, 23.2625, 31.2592, 28.8158, 22.2174, 25.895, 27.872, 
                  #                                                                      21.18535, 16.4624, 9.543875, 6.33563, 5.19942, 5.683845, 2.03261, 
                  #                                                                      10.6996) / (reduction_factor^2)
                  # true_mean_HADS_params_list$Ord_bivariate$alpha_d           <- c(2.82742, 1.33769, 7.65691, 7.55662, 10.6689, 16.008,
                  #                                                                 27.173,  28.8604, 39.5765, 54.9949, 51.4344, 40.3991, 
                  #                                                                 48.6174, 52.1399, 41.3584, 32.7035, 18.4987, 12.7581, 
                  #                                                                 10.8704, 12.6957, 4.13052, 6.33936) / (reduction_factor^2)
                  true_mean_HADS_params_list$Ord_bivariate$alpha_d           <- c(3.16713, 1.74906, 9.1371, 9.03822, 12.3724, 18.419,
                                                                                  31.1213, 33.1508, 45.16, 62.7256, 58.451, 46.2714,
                                                                                  55.1144, 59.401, 46.4115, 36.4505, 21.1506, 14.3548,
                                                                                  12.2054, 14.4073, 4.84256, 6.79661) / (reduction_factor^2)
            
            } else if (DGM == "ord_biv_fixed_C") {
                 
                  # true_mean_HADS_params_list$Ord_bivariate$beta_mu <- c(-1.13, 0.136)
                  true_mean_HADS_params_list$Ord_bivariate$beta_mu <- c(-1.13167, 0.140572)
                  # true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.248, 0.423)*reduction_factor
                  true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.246743, 0.417525)*reduction_factor
                  # true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.295
                  true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.3015
                  ##
                  # true_mean_HADS_params_list$Ord_bivariate$C_nd <- c(-2.63468, -2.195915, -1.85045, -1.54265, -1.2776, -1.029515,
                  #                                                    -0.752784, -0.4990785, -0.2554335, -0.0422625, 0.1345545, 0.330531,
                  #                                                    0.5208185, 0.731341, 0.966377, 1.151985, 1.37807, 1.58074, 1.782565,
                  #                                                    2.0821, 2.582115)
                  true_mean_HADS_params_list$Ord_bivariate$C_nd <- c(-2.64194, -2.20272, -1.85626, -1.55009, -1.28585, -1.03681, 
                                                                     -0.760486, -0.506951, -0.26366, -0.0494368, 0.127302, 0.322418, 
                                                                     0.513391, 0.721886, 0.95937, 1.14338, 1.36828, 1.57127,
                                                                     1.77154, 2.06873, 2.57198)
                  # true_mean_HADS_params_list$Ord_bivariate$C_d  <- c(-2.4333, -2.294975, -1.901415, -1.65488, -1.44765, -1.234255,
                  #                                                    -0.9666575, -0.7418255, -0.49973, -0.2126005, 0.0370325, 0.2299425,
                  #                                                    0.4643065, 0.742412, 1.008775, 1.28073, 1.464895, 1.62725, 1.802565,
                  #                                                    2.17232, 2.39557)
                  true_mean_HADS_params_list$Ord_bivariate$C_d  <- c(-2.44652, -2.31948, -1.92772, -1.68312, -1.47943, -1.26532, 
                                                                     -0.997867, -0.775762, -0.532074, -0.246346, 0.00157962, 0.195771, 
                                                                     0.427589, 0.703408, 0.96664, 1.23754, 1.42198, 1.58223, 
                                                                     1.75609, 2.12563, 2.3441)
                  # ##
                  # true_mean_HADS_params_list$Ord_bivariate$beta_mu <- c(-1.55, 0.567)
                  # true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.252, 0.424)*reduction_factor
                  # true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.286
                  # ##
                  # true_mean_HADS_params_list$Ord_bivariate$C_nd <- c(-3.0599237, -2.6212467, -2.275351, -1.9692706, -1.7051523, 
                  #                                                    -1.4556182, -1.1794026, -0.9259213, -0.6828896, -0.4686991, -0.2913007, 
                  #                                                    -0.0959831, 0.095133, 0.3035555, 0.5400546, 0.7236019, 0.9469161, 
                  #                                                    1.1510247, 1.3530062, 1.6493928, 2.107553)
                  # true_mean_HADS_params_list$Ord_bivariate$C_d  <- c(-2.0010324, -1.846439, -1.4744186, -1.2333598, -1.028526, -0.8166647, 
                  #                                                    -0.5535559, -0.3306749, -0.0889998, 0.1935292, 0.4385729, 0.6326791, 
                  #                                                    0.8622773, 1.1349685, 1.397572, 1.6652311, 1.8549768, 2.0227153, 
                  #                                                    2.2059966, 2.5734272, 2.8180254)
              
            } else if (DGM == "Jones") {
              
                  true_mean_HADS_params_list$Jones$beta_mu <- c(3.20, 9.38)
                  true_mean_HADS_params_list$Jones$beta_SD <- c(0.516, 1.38)*reduction_factor
                  true_mean_HADS_params_list$Jones$beta_corr <- 0.211
                  ##
                  # true_mean_HADS_params_list$Jones$cts_thr_values <- seq(from = 1, to = 21, by = 1)
                  # true_mean_HADS_params_list$Jones$lambda <- c(0.673, 0.936)
                  ##
                  true_mean_HADS_params_list$Jones$C_nd <-  c(0, 0.8830485, 1.626035, 2.290555, 2.90222, 3.474825, 4.01685, 
                                                              4.53401, 5.03037, 5.508955, 5.9721, 6.421635, 6.859055, 7.285585, 
                                                              7.702245, 8.1099, 8.50929, 8.90106, 9.285765, 9.663885, 10.03585)
                  true_mean_HADS_params_list$Jones$C_d  <- c(0, 0.975819, 1.919665, 2.843325, 3.75226, 4.64964, 5.53751, 
                                                             6.417325, 7.29016, 8.156835, 9.018015, 9.87422, 10.72585, 11.5734, 
                                                             12.4171, 13.2572, 14.09395, 14.9276, 15.7583, 16.58625, 17.4115)
                  ##
                  true_mean_HADS_params_list$Jones$raw_scale_mu <- c(0.776, 1.23)
                  true_mean_HADS_params_list$Jones$raw_scale_SD <- c(0.116, 0.116)*reduction_factor
                  true_mean_HADS_params_list$Jones$raw_scale_corr <- 0.269
              
            } else if (DGM == "R&G_random_C") {
                  
                  # true_mean_HADS_params_list$Ord_HSROC$beta_mu <- -0.745
                  true_mean_HADS_params_list$Ord_HSROC$beta_mu <- -0.7677
                  # true_mean_HADS_params_list$Ord_HSROC$beta_SD <-  0.204*reduction_factor
                  true_mean_HADS_params_list$Ord_HSROC$beta_SD <-  0.2103*reduction_factor
                  ##
                  # true_mean_HADS_params_list$Ord_HSROC$raw_scale_mu <-  -0.000988
                  true_mean_HADS_params_list$Ord_HSROC$raw_scale_mu <-  -0.0222
                  # true_mean_HADS_params_list$Ord_HSROC$raw_scale_SD <-   0.117*reduction_factor
                  true_mean_HADS_params_list$Ord_HSROC$raw_scale_SD <-   0.1287*reduction_factor
                  ##
                  # true_mean_HADS_params_list$Ord_HSROC$C <- c(-2.208575, -1.775145, -1.43052, -1.139615, -0.8782825, -0.6452505, 
                  #                                            -0.397182, -0.1591875, 0.0759022, 0.3002705, 0.488208, 0.6809105, 
                  #                                            0.878373, 1.088115, 1.31888, 1.511295, 1.71618, 1.89627, 2.062955, 
                  #                                            2.260425, 2.385345)
                  true_mean_HADS_params_list$Ord_HSROC$C <- c(-2.29649, -1.83486, -1.47787, -1.17736, -0.907769, -0.665501,
                                                              -0.409621, -0.163745, 0.0770738, 0.307219, 0.499123, 0.696593,
                                                              0.90068, 1.1191, 1.35903, 1.56404, 1.78452, 1.98979, 
                                                              2.19635, 2.50281, 2.80724)
                  true_mean_HADS_params_list$Ord_HSROC$C_nd <- true_mean_HADS_params_list$Ord_HSROC$C
                  true_mean_HADS_params_list$Ord_HSROC$C_d  <- true_mean_HADS_params_list$Ord_HSROC$C
                  ##
                  # true_mean_HADS_params_list$Ord_HSROC$C_SD_prob_scale <-  c(0.0057212, 0.007391, 0.0091481, 0.0104899, 0.0115506, 0.0122124, 
                  #                                                           0.0134513, 0.013816, 0.0139312, 0.0135597, 0.0121398, 0.0117441, 
                  #                                                           0.011176, 0.0105672, 0.0098244, 0.0079068, 0.0070315, 0.0055991, 
                  #                                                           0.0045222, 0.004101, 0.0027446)
                  true_mean_HADS_params_list$Ord_HSROC$C_SD_prob_scale <-  c(0.00473229, 0.00658042, 0.00832652, 0.00965628, 0.0107397, 0.0113984, 
                                                                             0.0126041, 0.0129554, 0.0130747, 0.0126968, 0.0113753, 0.0110126,
                                                                             0.0104644, 0.00989028, 0.00914358, 0.00733038, 0.00645658, 0.00517164,
                                                                             0.00422725, 0.00389456, 0.00262681)
                  ##
                  # true_mean_HADS_params_list$Ord_HSROC$alpha <-  c(6.316175, 10.6233, 16.44925, 21.88265, 26.8957, 30.27135, 37.25105, 
                  #                                                 39.6049, 40.49045, 37.9255, 29.80595, 27.82585, 25.12695, 22.33755, 
                  #                                                 19.17145, 12.2595, 9.60213, 6.064095, 3.979375, 3.21631, 1.397125, 
                  #                                                 4.045135) / (reduction_factor^2)
                  # true_mean_HADS_params_list$Ord_HSROC$alpha <- c(5.8078, 11.3559, 18.4465, 25.1551, 31.5348, 35.8482, 
                  #                                                 44.672, 47.4861, 48.469, 45.4034, 35.6906, 33.2806, 
                  #                                                 29.8319, 26.461, 22.4276, 14.1726, 10.9229, 6.95217,
                  #                                                 4.62334, 3.91874, 1.77516, 1.58676) / (reduction_factor^2)
                  true_mean_HADS_params_list$Ord_HSROC$alpha <- c(5.81448, 11.3416, 18.5038, 25.1729, 31.5061, 35.8092,
                                                                  44.5978, 47.3811, 48.4394, 45.3039, 35.5748, 33.3057,
                                                                  29.8969, 26.5368, 22.4179, 14.3058, 10.9962, 6.96191,
                                                                  4.73241, 3.99418, 1.84492, 1.59377) / (reduction_factor^2)
              
            }
            ##
            ## ---- BAI parameters: ---------------------------------------------------------------------------------------------------
            ##
            ## (from fitting ord_biv-Rand-Kappa model to BAI HADS MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_BAI_params_list <- list()
            ##
            if (DGM == "ord_biv_random_C") {
              
                  # true_mean_BAI_params_list$Ord_bivariate$beta_mu <- c(-1.51, -0.347)
                  true_mean_BAI_params_list$Ord_bivariate$beta_mu   <- c(-1.54901, -0.400199)
                  # true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.373,  0.250)*reduction_factor
                  true_mean_BAI_params_list$Ord_bivariate$beta_SD   <- c(0.334592, 0.274101)*reduction_factor
                  # true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.0262
                  true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.1504
                  ##
                  # true_mean_BAI_params_list$Ord_bivariate$C_nd <-   c(-2.7115, -2.36566, -2.12313, -1.966295, -1.80715, -1.677715, 
                  #                                                     -1.52838, -1.422835, -1.315275, -1.2169, -1.119355, -1.035275, 
                  #                                                     -0.926281, -0.858484, -0.7834725, -0.6889895, -0.615242, -0.5419855, 
                  #                                                     -0.483087, -0.436301, -0.379436, -0.334049, -0.270132, -0.2054205, 
                  #                                                     -0.1521625, -0.1003645, -0.0553764, -0.0227455, 0.0165322, 0.0494109, 
                  #                                                     0.1047035, 0.1430875, 0.166283, 0.2047725, 0.237076, 0.2956025, 
                  #                                                     0.3482815, 0.4203975, 0.449386, 0.4788115, 0.5025445, 0.5480385, 
                  #                                                     0.575396, 0.628353, 0.7026125, 0.7394385, 0.8209555, 0.847542, 
                  #                                                     0.874193, 0.9086765, 0.9443915, 0.981077, 1.05217, 1.090975, 
                  #                                                     1.135785, 1.182495, 1.251555, 1.307605, 1.36542, 1.43101, 1.544135, 
                  #                                                     1.675555, 1.77604)
                  true_mean_BAI_params_list$Ord_bivariate$C_nd <- c(-2.75987, -2.38058, -2.13854, -1.96505, -1.79258, -1.66302, 
                                                                    -1.51412, -1.40725, -1.29856, -1.2135, -1.11691, -1.03217, 
                                                                    -0.923345, -0.857743, -0.783115, -0.687798, -0.613648, -0.538607,
                                                                    -0.480235, -0.430874, -0.374557, -0.328894, -0.262853, -0.196633, 
                                                                    -0.14101, -0.0876693, -0.0432382, -0.00814693, 0.0321933, 0.0675187,
                                                                    0.124321, 0.163258, 0.185753, 0.226458, 0.259769, 0.317945,
                                                                    0.37409, 0.448396, 0.476878, 0.506001, 0.532138, 0.580097, 
                                                                    0.609855, 0.666107, 0.744115, 0.784592, 0.871337, 0.896762, 
                                                                    0.923035, 0.958802, 0.998007, 1.03457, 1.11736, 1.1601,
                                                                    1.21116, 1.26416, 1.34189, 1.40711, 1.47668, 1.5532, 
                                                                    1.70537, 1.91261, 2.08012)
                  # true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale <-     c(0.002358, 0.0029225, 0.003367, 0.0035171, 0.0041432, 0.0041321, 
                  #                                                                     0.0049617, 0.0045639, 0.0049522, 0.0051236, 0.0054301, 0.0052936, 
                  #                                                                     0.0063219, 0.0051523, 0.0055877, 0.006542, 0.0058929, 0.0060159, 
                  #                                                                     0.0054111, 0.0049745, 0.0054706, 0.0049258, 0.0059544, 0.006023, 
                  #                                                                     0.0055156, 0.0053552, 0.0049767, 0.0043495, 0.004699, 0.0043609, 
                  #                                                                     0.0054868, 0.0046102, 0.0034478, 0.0046133, 0.0041473, 0.0055905, 
                  #                                                                     0.0053579, 0.0061172, 0.0036363, 0.0036762, 0.0033588, 0.0046603, 
                  #                                                                     0.0034891, 0.0049027, 0.0056647, 0.003851, 0.0056406, 0.0030387, 
                  #                                                                     0.0029721, 0.0032668, 0.0034083, 0.0032773, 0.0047283, 0.0032744, 
                  #                                                                     0.0034341, 0.0034647, 0.0039899, 0.0034031, 0.0033762, 0.0034169, 
                  #                                                                     0.0041718, 0.0042559, 0.0032206)
                  true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale <- c(0.00263281, 0.00336518, 0.0038723, 0.00404291, 0.00480719, 0.00478042,
                                                                                  0.0057275, 0.00532531, 0.00583639, 0.005454, 0.00619545, 0.00603648,
                                                                                  0.00718662, 0.00586598, 0.00639007, 0.00747298, 0.00676336, 0.00694946,
                                                                                  0.00625344, 0.00575804, 0.00629602, 0.00571784, 0.00691553, 0.00698501,
                                                                                  0.00637935, 0.00621736, 0.00581655, 0.00510076, 0.0055022, 0.00509098,
                                                                                  0.0064092, 0.00529466, 0.00392785, 0.00537956, 0.00479168, 0.00651336, 
                                                                                  0.00625356, 0.00709735, 0.00421312, 0.0042349, 0.00378844, 0.00547343, 
                                                                                  0.00410325, 0.00578289, 0.00664699, 0.00446988, 0.00665797, 0.00318332,
                                                                                  0.00320855, 0.00369558, 0.00378554, 0.00364327, 0.00563777, 0.00363642, 
                                                                                  0.00390008, 0.00393149, 0.00472193, 0.00375125, 0.003775, 0.00376099,
                                                                                  0.00502073, 0.00514945, 0.003564)
                  # true_mean_BAI_params_list$Ord_bivariate$alpha_nd         <-    c(2.30978, 3.574155, 4.67112, 5.126565, 7.194915, 7.083335, 10.26365, 
                  #                                                                  8.71072, 10.35995, 10.97255, 12.383, 11.82935, 16.83, 11.1916, 
                  #                                                                  13.16695, 18.13275, 14.68875, 15.34365, 12.48285, 10.49845, 12.6014, 
                  #                                                                  10.4, 15.1731, 15.54295, 13.01245, 12.2505, 10.4292, 8.026365, 
                  #                                                                  9.29906, 8.06027, 13.0912, 9.00764, 4.96975, 8.9957, 7.190805, 
                  #                                                                  13.23915, 12.08715, 15.98565, 5.567645, 5.60625, 4.704465, 9.25797, 
                  #                                                                  5.18286, 10.29925, 13.7561, 6.267155, 13.52375, 3.701475, 3.57462, 
                  #                                                                  4.373335, 4.747075, 4.370125, 9.416955, 4.34419, 4.775305, 4.843085, 
                  #                                                                  6.740945, 4.787045, 4.656905, 4.841155, 7.229655, 7.43725, 4.21635, 
                  #                                                                  24.05) / (reduction_factor^2)
                  # true_mean_BAI_params_list$Ord_bivariate$alpha_nd          <- c(2.01752, 3.30398, 4.38366, 4.78201, 6.78643, 6.7101,
                  #                                                                9.68648, 8.35278, 10.0655, 8.76822, 11.3702, 10.7822,
                  #                                                                15.4178, 10.1698, 12.1128, 16.7123, 13.6081, 14.3887,
                  #                                                                11.5888, 9.79201, 11.7508, 9.65324, 14.2446, 14.5405,
                  #                                                                12.0713, 11.4525, 9.99589, 7.65311, 8.92657, 7.62336,
                  #                                                                12.1872, 8.25539, 4.51142, 8.52663, 6.74214, 12.5963, 
                  #                                                                11.5893, 15.0259, 5.19719, 5.25159, 4.19436, 8.83187,
                  #                                                                4.92717, 9.87832, 13.1319, 5.8572, 13.1764, 2.9546, 
                  #                                                                3.00188, 3.98973, 4.18789, 3.87676, 9.37991, 3.86208,
                  #                                                                4.44732, 4.51985, 6.54485, 4.11178, 4.16442, 4.13333,
                  #                                                                7.41145, 7.8021, 3.70874, 11.2956) / (reduction_factor^2)
                  true_mean_BAI_params_list$Ord_bivariate$alpha_nd          <- c(2.00899, 3.39931, 4.48468, 4.91602, 6.928, 6.91689,
                                                                                 9.93505, 8.55651, 10.409, 9.09332, 11.6806, 11.1328,
                                                                                 15.8604, 10.6821, 12.5374, 17.2988, 14.2462, 14.9918,
                                                                                 12.1563, 10.3636, 12.5113, 10.4074, 14.989, 15.3816,
                                                                                 12.8416, 12.3251, 10.6864, 8.41658, 9.64992, 8.51466,
                                                                                 13.2948, 9.17625, 5.25954, 9.42787, 7.60462, 13.5917,
                                                                                 12.4927, 16.3391, 6.25261, 6.26359, 5.13549, 9.97028,
                                                                                 5.93225, 11.124, 14.8791, 7.01249, 14.6879, 4.20015,
                                                                                 4.27554, 5.57956, 5.6577, 5.43396, 10.954, 5.36349,
                                                                                 6.02921, 5.85821, 8.28843, 6.00386, 5.98133, 5.91644,
                                                                                 9.87254, 9.75062, 5.61821, 11.6388) / (reduction_factor^2)
                  ##
                  # true_mean_BAI_params_list$Ord_bivariate$C_d <-     c(-2.438265, -2.32965, -2.1953, -2.089345, -1.89032, -1.79675, 
                  #                                                      -1.71302, -1.60725, -1.53188, -1.441385, -1.396705, -1.27583, 
                  #                                                      -1.206545, -1.13503, -1.02006, -0.8973925, -0.8280695, -0.694393, 
                  #                                                      -0.625531, -0.564262, -0.5104715, -0.4560765, -0.3664795, -0.314171, 
                  #                                                      -0.245043, -0.217004, -0.150391, -0.0858726, -0.0442238, -0.0225385, 
                  #                                                      0.0681136, 0.1401925, 0.1774055, 0.2155285, 0.2551105, 0.3271175, 
                  #                                                      0.3494405, 0.4037965, 0.480177, 0.519243, 0.5864515, 0.6856125, 
                  #                                                      0.7236915, 0.7617215, 0.8003375, 0.860833, 0.925063, 0.9863845, 
                  #                                                      1.051545, 1.0928, 1.13591, 1.174145, 1.21501, 1.2753, 1.33788, 
                  #                                                      1.417865, 1.46565, 1.52343, 1.582555, 1.65265, 1.727495, 1.815865, 
                  #                                                      2.044695)
                  true_mean_BAI_params_list$Ord_bivariate$C_d  <- c(-2.84163, -2.63838, -2.42567, -2.28053, -2.05173, -1.9445, 
                                                                    -1.85112, -1.74589, -1.66714, -1.55206, -1.50248, -1.3751, 
                                                                    -1.29873, -1.22632, -1.10884, -0.981541, -0.894072, -0.73193,
                                                                    -0.645891, -0.580412, -0.519586, -0.461155, -0.362704, -0.304221,
                                                                    -0.229496, -0.194299, -0.12002, -0.0478964, -0.000312958, 0.0260897,
                                                                    0.115449, 0.187284, 0.223866, 0.261018, 0.300172, 0.373481,
                                                                    0.394434, 0.449233, 0.525329, 0.5672, 0.635365, 0.734499,
                                                                    0.771781, 0.812666, 0.856606, 0.919871, 0.989615, 1.05387,
                                                                    1.12522, 1.17135, 1.21735, 1.26008, 1.30284, 1.37285,
                                                                    1.44418, 1.53938, 1.60374, 1.67454, 1.75343, 1.84607,
                                                                    1.95442, 2.09384, 2.62229)
                  # true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale <-   c(0.0020563, 0.0010478, 0.0014252, 0.0014155, 0.0023965, 0.0018293, 
                  #                                                                  0.0018722, 0.002288, 0.0020586, 0.0024877, 0.001788, 0.0032144, 
                  #                                                                  0.0025958, 0.0026254, 0.0036675, 0.004052, 0.0031747, 0.0046054, 
                  #                                                                  0.0034303, 0.0031963, 0.0031116, 0.0031368, 0.0041784, 0.0031939, 
                  #                                                                  0.0036703, 0.0023667, 0.0036828, 0.003605, 0.0028875, 0.0020361, 
                  #                                                                  0.0043061, 0.0038247, 0.0026328, 0.0026575, 0.0027337, 0.0037605, 
                  #                                                                  0.0018628, 0.0031602, 0.0037008, 0.0025156, 0.0033907, 0.0039843, 
                  #                                                                  0.002241, 0.0022437, 0.0023091, 0.0028717, 0.0028866, 0.0026356, 
                  #                                                                  0.0025216, 0.0019958, 0.0020339, 0.001818, 0.0017884, 0.0021723, 
                  #                                                                  0.0020995, 0.0023874, 0.001653, 0.0017341, 0.0017044, 0.0017621, 
                  #                                                                  0.0017089, 0.0017028, 0.0024921)
                  true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale  <- c(0.00204623, 0.00147964, 0.00212767, 0.00215935, 0.00351402, 0.00283989,
                                                                                  0.00287449, 0.00338529, 0.00316408, 0.00419069, 0.00283558, 0.00512108, 
                                                                                  0.00417727, 0.00425657, 0.00586845, 0.00645617, 0.00569667, 0.00818151, 
                                                                                  0.00620487, 0.00555997, 0.00541135, 0.00543983, 0.00716381, 0.00562081,
                                                                                  0.00641342, 0.00427394, 0.0064894,  0.00635618, 0.00520147, 0.0036376,
                                                                                  0.00710496, 0.00631586, 0.00443879, 0.00439093, 0.00448195, 0.00619094,
                                                                                  0.0030139,  0.00529599, 0.00609203, 0.00427549, 0.00561541, 0.00661821,
                                                                                  0.00375507, 0.00381271, 0.003917,   0.0048498,  0.00484544, 0.00441006, 
                                                                                  0.00438091, 0.00336106, 0.00327413, 0.002922,   0.00292073, 0.00370182,
                                                                                  0.0035439,  0.00400926, 0.00285343, 0.00284943, 0.00280725, 0.00280539,
                                                                                  0.00279911, 0.00277995, 0.00423732)
                  # true_mean_BAI_params_list$Ord_bivariate$alpha_d           <-    c(13.0355, 3.031365, 5.812175, 5.896445, 17.5918, 10.33855, 10.70605, 
                  #                                                                   16.1598, 12.732, 19.037, 9.291165, 31.2761, 20.65685, 21.02305, 
                  #                                                                   41.2895, 50.29585, 30.8706, 64.7911, 36.12735, 31.221, 29.44055, 
                  #                                                                   30.11085, 53.0706, 31.44215, 42.4023, 16.68195, 41.44105, 39.4167, 
                  #                                                                   25.69055, 12.2199, 57.55395, 44.49145, 21.39975, 22.1469, 22.3883, 
                  #                                                                   44.14655, 10.11145, 31.1491, 41.368, 18.7385, 35.0389, 49.74355, 
                  #                                                                   14.86755, 14.886, 16.3358, 25.62505, 25.6765, 20.551, 20.1573, 
                  #                                                                   11.47565, 11.4116, 9.53656, 8.99987, 13.751, 13.53585, 17.1601, 
                  #                                                                   8.0716, 8.64498, 8.12516, 8.99806, 8.50395, 8.33532, 19.12505, 
                  #                                                                   33.42485) / (reduction_factor^2)
                  # true_mean_BAI_params_list$Ord_bivariate$alpha_d           <- c(2.16084, 1.12824, 2.33685, 2.4072,  6.41083, 4.17389, 
                  #                                                                4.27684, 5.94583, 5.18865, 9.15284, 4.16117, 13.7576,
                  #                                                                9.09355, 9.44678, 18.1806, 22.1295, 17.1054, 36.2764,
                  #                                                                20.3893, 16.2751, 15.3973, 15.5635, 27.457,  16.6419,
                  #                                                                21.828,  9.52513, 22.3656, 21.4277, 14.2019, 6.87415, 
                  #                                                                26.9895, 21.1482, 10.2851, 10.0613, 10.4891, 20.2951,
                  #                                                                4.70457, 14.7339, 19.6331, 9.53213, 16.6091, 23.2934,
                  #                                                                7.33,    7.55918, 7.98317, 12.3133, 12.2908, 10.1505, 
                  #                                                                10.0148, 5.8603,  5.55874, 4.42028, 4.41643, 7.12147, 
                  #                                                                6.52133, 8.36819, 4.21405, 4.20216, 4.07797, 4.07252, 
                  #                                                                4.05421, 3.99857, 9.36042, 3.8589) / (reduction_factor^2)
                  true_mean_BAI_params_list$Ord_bivariate$alpha_d           <- c(2.59618, 1.85877, 3.32728, 3.46307, 8.44884, 5.73802,
                                                                                 5.90851, 8.07855, 7.00664, 12.0885, 5.97491, 17.8712,
                                                                                 12.2939, 12.4797, 23.2808, 28.6875, 21.8652, 45.9727,
                                                                                 26.2045, 21.4794, 19.7809, 20.4896, 35.1546, 21.5541,
                                                                                 28.3403, 13.2437, 28.6595, 27.9058, 18.9864, 9.73877,
                                                                                 35.1387, 27.5108, 14.3272, 14.2004, 14.7024, 26.93,
                                                                                 7.47345, 20.012, 26.4276, 13.8518, 22.8029, 30.7472,
                                                                                 11.4051, 11.4962, 11.5193, 16.7147, 17.2678, 15.0079,
                                                                                 14.635, 9.45618, 9.161, 7.56052, 7.39364, 10.7696,
                                                                                 10.8349, 12.1793, 7.12504, 7.10511, 7.2552, 7.14676,
                                                                                 7.18202, 7.03166, 13.7592, 4.9269) / (reduction_factor^2)
            
            } else if (DGM == "ord_biv_fixed_C") {
                  
                  # true_mean_BAI_params_list$Ord_bivariate$beta_mu <- c(-1.49, -0.358)
                  true_mean_BAI_params_list$Ord_bivariate$beta_mu <- c(-1.53432, -0.406964)
                  # true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.385, 0.268)*reduction_factor
                  true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.347129, 0.273508)*reduction_factor
                  # true_mean_BAI_params_list$Ord_bivariate$beta_corr <- -0.0666
                  true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.06378
                  ##
                  # true_mean_BAI_params_list$Ord_bivariate$C_nd <- c(-2.697165, -2.366735, -2.118675, -1.95712, -1.78154, -1.658235, 
                  #                                                    -1.509775, -1.404975, -1.30039, -1.199635, -1.099655, -1.010015, 
                  #                                                    -0.902952, -0.840962, -0.766334, -0.678289, -0.609417, -0.539141, 
                  #                                                    -0.4829525, -0.434512, -0.382427, -0.3349525, -0.270738, -0.203197, 
                  #                                                    -0.14692, -0.0978044, -0.0560808, -0.0241374, 0.0147843, 0.0497706, 
                  #                                                    0.1088645, 0.1515355, 0.1719505, 0.210362, 0.242857, 0.301478, 
                  #                                                    0.361064, 0.4373065, 0.4652915, 0.4935335, 0.5156475, 0.5640035, 
                  #                                                    0.592413, 0.6513545, 0.7369215, 0.775247, 0.8616, 0.886257, 0.9118475, 
                  #                                                    0.9457825, 0.980719, 1.016785, 1.096365, 1.137085, 1.184045, 
                  #                                                    1.235355, 1.313385, 1.372655, 1.43508, 1.50447, 1.650595, 1.837665, 
                  #                                                    1.98404)
                  true_mean_BAI_params_list$Ord_bivariate$C_nd <- c(-2.69437,  -2.36391,  -2.11668,   -1.95548,   -1.78004,  -1.65627, 
                                                                    -1.50794,  -1.40213,  -1.29505,   -1.21028,   -1.11065,  -1.02151,
                                                                    -0.91535,  -0.854205, -0.780985,  -0.693652,  -0.625437, -0.55653,
                                                                    -0.501414, -0.453847, -0.402047,  -0.355303,  -0.292546, -0.225689,
                                                                    -0.170168, -0.121495, -0.0811179, -0.0498905, -0.0117294, 0.0222546,
                                                                    0.079519,   0.119483,  0.139067,   0.174589,   0.203929,  0.260448,
                                                                    0.316813,   0.39027,   0.414157,   0.438953,   0.458875,  0.503368,
                                                                    0.527106,   0.582149,  0.661236,   0.694387,   0.77607,   0.793218, 
                                                                    0.81075,    0.836378,  0.862687,   0.890012,   0.962354,  0.990933,
                                                                    1.02849,    1.06799,   1.13473,    1.17951,    1.22883,   1.28392,
                                                                    1.41102, 1.57869, 1.69594)
                  # true_mean_BAI_params_list$Ord_bivariate$C_d  <- c(-2.690125, -2.54255, -2.35055, -2.21676, -1.98042, -1.87441, 
                  #                                                    -1.78187, -1.664335, -1.5829, -1.486535, -1.44155, -1.31506, 
                  #                                                    -1.24019, -1.171065, -1.04944, -0.923497, -0.851394, -0.7143255, 
                  #                                                    -0.6413, -0.5821935, -0.5266825, -0.4708045, -0.378942, -0.3236045, 
                  #                                                    -0.251782, -0.221604, -0.15677, -0.0905013, -0.0465355, -0.025449, 
                  #                                                    0.0696832, 0.1426975, 0.1820825, 0.221577, 0.2602225, 0.3367695, 
                  #                                                    0.357437, 0.4147505, 0.4931115, 0.5361075, 0.6090605, 0.71257, 
                  #                                                    0.7507885, 0.790541, 0.830872, 0.8966195, 0.971349, 1.0344, 1.10528, 
                  #                                                    1.1516, 1.197425, 1.23752, 1.27867, 1.34473, 1.41797, 1.51321, 
                  #                                                    1.57391, 1.639725, 1.710355, 1.79398, 1.892395, 2.016135, 2.457965)
                  true_mean_BAI_params_list$Ord_bivariate$C_d  <- c(-2.8142,   -2.69419,   -2.5164,     -2.39809,  -2.17338,  -2.07801,
                                                                    -1.99709,  -1.88897,   -1.82124,    -1.7113,   -1.67002,  -1.54451,
                                                                    -1.47217,  -1.40547,   -1.28565,    -1.16094,  -1.07677,  -0.914916, 
                                                                    -0.831223, -0.768277,  -0.710391,   -0.653142, -0.556331, -0.499331,
                                                                    -0.424771, -0.39468,   -0.323484,   -0.255233, -0.210543, -0.189597, 
                                                                    -0.101793, -0.0333945, -0.000248207, 0.0334523, 0.0677788, 0.138055,
                                                                    0.152175,   0.204094,   0.276681,    0.311134,  0.378554,  0.475018,
                                                                    0.506227,   0.539071,   0.574464,    0.633392,  0.697522,  0.752151,
                                                                    0.810939,   0.845225,   0.881949,    0.910993,  0.940501,  0.993736, 
                                                                    1.0516,     1.13243,    1.17552,     1.22403,   1.27868,   1.34159,
                                                                    1.41693, 1.51298, 1.96908)
              
            } else if (DGM == "Jones") {
              
                  true_mean_BAI_params_list$Jones$beta_mu <- c(2.88, 8.42)
                  true_mean_BAI_params_list$Jones$beta_SD <- c(0.904, 0.898)*reduction_factor
                  true_mean_BAI_params_list$Jones$beta_corr <- -0.120
                  ##
                  true_mean_BAI_params_list$Jones$C_nd <-  c(0, 0.7964125, 1.373, 1.84136, 2.2428, 2.59789, 2.918545, 3.212405, 
                                                             3.484685, 3.739125, 3.978525, 4.205025, 4.420305, 4.625745, 4.82244, 
                                                             5.011315, 5.193135, 5.36856, 5.538155, 5.702405, 5.861735, 6.01651, 
                                                             6.16707, 6.31369, 6.45665, 6.59617, 6.73246, 6.86571, 6.9961, 
                                                             7.12379, 7.24891, 7.3716, 7.49198, 7.610165, 7.72625, 7.840335, 
                                                             7.95251, 8.062855, 8.171455, 8.278365, 8.383665, 8.487415, 8.589675, 
                                                             8.690495, 8.78993, 8.888035, 8.984845, 9.080415, 9.17478, 9.267985, 
                                                             9.360055, 9.45104, 9.540965, 9.629865, 9.71777, 9.804705, 9.8907, 
                                                             9.97579, 10.06, 10.1433, 10.2258, 10.3075, 10.38835)
                  true_mean_BAI_params_list$Jones$C_d  <- c(0, 0.8439165, 1.509195, 2.079825, 2.58904, 3.054115, 3.48543, 
                                                             3.889805, 4.27201, 4.63554, 4.983065, 5.316645, 5.63794, 5.9483, 
                                                             6.24884, 6.5405, 6.824075, 7.10024, 7.369585, 7.63263, 7.88983, 
                                                             8.141575, 8.38822, 8.6301, 8.867485, 9.10063, 9.329785, 9.555155, 
                                                             9.77693, 9.995295, 10.2104, 10.4224, 10.6315, 10.8377, 11.04125, 
                                                             11.24215, 11.4406, 11.6366, 11.83035, 12.02185, 12.2112, 12.3985, 
                                                             12.58385, 12.76715, 12.94865, 13.1284, 13.30635, 13.48265, 13.65725, 
                                                             13.83025, 14.00175, 14.1717, 14.3402, 14.5073, 14.673, 14.83735, 
                                                             15.0004, 15.16215, 15.3227, 15.4819, 15.64, 15.797, 15.9528)
                  ##
                  true_mean_BAI_params_list$Jones$raw_scale_mu <- c(0.881, 1.14)
                  true_mean_BAI_params_list$Jones$raw_scale_SD <- c(0.166, 0.358)*reduction_factor
                  true_mean_BAI_params_list$Jones$raw_scale_corr <- 0.299
                  
            } else if (DGM == "R&G_random_C") {
                  
                  # true_mean_BAI_params_list$Ord_HSROC$beta_mu <- -0.764
                  true_mean_BAI_params_list$Ord_HSROC$beta_mu <- -0.7887
                  # true_mean_BAI_params_list$Ord_HSROC$beta_SD <-  0.333*reduction_factor
                  true_mean_BAI_params_list$Ord_HSROC$beta_SD <-  0.2983*reduction_factor
                  ##
                  # true_mean_BAI_params_list$Ord_HSROC$raw_scale_mu <-  -0.00606
                  true_mean_BAI_params_list$Ord_HSROC$raw_scale_mu <-  -0.03466
                  # true_mean_BAI_params_list$Ord_HSROC$raw_scale_SD <-   0.149*reduction_factor
                  true_mean_BAI_params_list$Ord_HSROC$raw_scale_SD <-   0.1332*reduction_factor
                  ##
                  # true_mean_BAI_params_list$Ord_HSROC$C <- c(-2.00031, -2.0209, -2.0261, -2.09921, -1.992025, -2.011525, 
                  #                                               -2.00169, -1.6515, -1.678785, -1.698, -1.74414, -1.65547, -1.676645, 
                  #                                               -1.635135, -1.400465, -1.430955, -1.4554, -1.491285, -1.387285, 
                  #                                               -1.41443, -1.406005, -1.230205, -1.266045, -1.294405, -1.296525, 
                  #                                               -1.235635, -1.257335, -1.23759, -1.04746, -1.08462, -1.105655, 
                  #                                               -1.10106, -1.07395, -1.069125, -1.060525, -0.916482, -0.957749, 
                  #                                               -0.977008, -0.9675875, -0.9465195, -0.9471125, -0.9343785, -0.764437, 
                  #                                               -0.809607, -0.818484, -0.821121, -0.796313, -0.7966565, -0.79358, 
                  #                                               -0.6526225, -0.698821, -0.7049935, -0.7104555, -0.700862, -0.6848455, 
                  #                                               -0.6813505, -0.5443485, -0.5909885, -0.5917325, -0.5981345, -0.603597, 
                  #                                               -0.584624, -0.5752315)
                  true_mean_BAI_params_list$Ord_HSROC$C <- c(-2.06206, -1.69185, -1.43144, -1.25641, -1.06578, -0.93126, -0.773938, -0.659471, -0.544386, -0.455821, -0.360645, -0.267539, -0.161271, -0.0952302, -0.0142139, 0.0874505, 0.165217, 0.260516, 0.32489, 0.378592, 0.435582, 0.486651, 0.564139, 0.630071, 0.692746, 0.738611, 0.794555, 0.842266, 0.885529, 0.917694, 0.984457, 1.03412, 1.05972, 1.09575, 1.12922, 1.19055, 1.22922, 1.28988, 1.33294, 1.36292, 1.4021, 1.46718, 1.49483, 1.54282, 1.60194, 1.6484, 1.7197, 1.75629, 1.79685, 1.83223, 1.86781, 1.90014, 1.96098, 2.0079, 2.06278, 2.13009, 2.19448, 2.24932, 2.31215, 2.38432, 2.50081, 2.66129, 2.97197)
                  true_mean_BAI_params_list$Ord_HSROC$C_nd <- true_mean_BAI_params_list$Ord_HSROC$C
                  true_mean_BAI_params_list$Ord_HSROC$C_d  <- true_mean_BAI_params_list$Ord_HSROC$C
                  ##
                  # true_mean_BAI_params_list$Ord_HSROC$C_SD_prob_scale <- c(0.0041522, 0.0044421, 0.0048395, 0.0045842, 0.0053558, 0.004881, 
                  #                                                          0.0056063, 0.0050442, 0.0051693, 0.0051662, 0.0051382, 0.0051939, 
                  #                                                          0.0055983, 0.0044233, 0.0049543, 0.0055254, 0.0046796, 0.0051713, 
                  #                                                          0.0042187, 0.0038642, 0.0039497, 0.0037365, 0.0045396, 0.0040922, 
                  #                                                          0.0039323, 0.0032558, 0.0035338, 0.0031628, 0.003003, 0.0024837, 
                  #                                                          0.0038071, 0.0031435, 0.0021583, 0.0025582, 0.0023562, 0.0032352, 
                  #                                                          0.0024672, 0.0030622, 0.0024177, 0.0018974, 0.002189, 0.0027602, 
                  #                                                          0.0016566, 0.0021732, 0.0023364, 0.0019803, 0.0023578, 0.0015359, 
                  #                                                          0.0015725, 0.0014093, 0.0013785, 0.0013047, 0.001633, 0.0014093, 
                  #                                                          0.0013918, 0.0014907, 0.0013435, 0.0011966, 0.0011561, 0.0011473, 
                  #                                                          0.0013133, 0.0013248, 0.0013153)
                  true_mean_BAI_params_list$Ord_HSROC$C_SD_prob_scale <-  c(0.00529953, 0.00595171, 0.00649938, 0.00625365, 0.00726905, 0.00669293, 0.00767, 0.00696593, 0.00720422, 0.00655392, 0.00689347, 0.00692299, 0.00750302, 0.00597515, 0.00665384, 0.00739675, 0.00650085, 0.00711127, 0.00583602, 0.00524506, 0.00533923, 0.00502816, 0.00603297, 0.00549033, 0.0052491, 0.00438612, 0.00480233, 0.00429751, 0.00401657, 0.00337551, 0.0048028, 0.00402643, 0.00282754, 0.00329816, 0.00308528, 0.00412297, 0.00315414, 0.00384245, 0.00309634, 0.00249663, 0.00282314, 0.00350731, 0.00213343, 0.00282605, 0.00300289, 0.00256438, 0.00299784, 0.00198066, 0.0020053, 0.00174521, 0.00172, 0.00154732, 0.00217199, 0.00173998, 0.00179503, 0.00192303, 0.00171473, 0.00144007, 0.00143016, 0.00140908, 0.00164287, 0.00162519, 0.00172041)
                  ##
                  # true_mean_BAI_params_list$Ord_HSROC$alpha <- c(26.8374, 30.3177, 36.38685, 32.395, 44.41205, 36.57935, 48.847, 
                  #                                                39.63305, 41.5899, 41.4163, 41.16915, 41.49585, 49.04665, 29.9765, 
                  #                                                38.2134, 47.50965, 34.04035, 41.7983, 27.6553, 23.5528, 24.39075, 
                  #                                                21.86925, 32.267, 26.2921, 24.17695, 16.9432, 19.5393, 15.76975, 
                  #                                                14.26195, 9.94385, 22.96825, 15.84405, 7.69775, 10.4792, 9.00526, 
                  #                                                16.77205, 10.0927, 14.92775, 9.793295, 5.806525, 7.903205, 12.31325, 
                  #                                                4.496865, 7.799915, 8.995255, 6.373515, 9.105355, 3.86026, 3.99084, 
                  #                                                3.19067, 2.9886, 2.740055, 4.27295, 3.112475, 2.99751, 3.59191, 
                  #                                                2.9406, 2.192515, 2.102795, 2.071835, 2.699475, 2.761735, 2.810295, 
                  #                                                5.680975) / (reduction_factor^2)
                  # true_mean_BAI_params_list$Ord_HSROC$alpha <- c(14.7534, 18.7136, 22.4357, 20.7201, 28.3033, 23.8401, 31.6664, 25.9017, 27.7796, 22.8267, 25.3451, 25.571, 30.2396, 18.8654, 23.5526, 29.3512, 22.4462, 27.0382, 17.9741, 14.4453, 14.9801, 13.2528, 19.2427, 15.8599, 14.4681, 10.0386, 12.0688, 9.63149, 8.39876, 5.91099, 12.0712, 8.44053, 4.1373, 5.64104, 4.93143, 8.8553, 5.15563, 7.67854, 4.9671, 3.22144, 4.12435, 6.38585, 2.34947, 4.13292, 4.66984, 3.39949, 4.65405, 2.02409, 2.07492, 1.57048, 1.52534, 1.23394, 2.43545, 1.56107, 1.66163, 1.90771, 1.51598, 1.06856, 1.05389, 1.023, 1.39134, 1.36151, 1.52606, 1.41464) / (reduction_factor^2)
                  true_mean_BAI_params_list$Ord_HSROC$alpha <- c(14.4083, 18.4071, 21.9376, 20.2943, 27.6195, 23.3319,
                                                                 30.9705, 25.1424, 27.3267, 22.3746, 24.7933, 25.2362,
                                                                 29.6786, 18.586, 22.9406, 28.8046, 22.1909, 26.634,
                                                                 17.6295, 14.3911, 14.8432, 13.3731, 18.8943, 15.7842,
                                                                 14.2463, 10.1119, 11.8737, 9.62399, 8.48367, 6.17505,
                                                                 12.1466, 8.4787, 4.29961, 5.73965, 5.15759, 8.92505,
                                                                 5.3244, 7.7928, 5.18834, 3.44257, 4.37109, 6.39004,
                                                                 2.62691, 4.38675, 4.81803, 3.59881, 4.92485, 2.33708,
                                                                 2.39442, 1.90837, 1.89408, 1.60579, 2.67114, 1.92655,
                                                                 1.89542, 2.14739, 1.70413, 1.34222, 1.30965, 1.30336,
                                                                 1.73933, 1.68131, 1.68598, 1.35361) / (reduction_factor^2)
 
            }
            
          }
          
          # # Set model type and parameters
          # bs_model <- "ord_bivariate"
          # if (bs_model == "ord_HSROC") {
          #   # HSROC model code would go here
          # } else {
          # 
          # }
          if (DGM %in% c("ord_biv_fixed_C",
                         "ord_biv_random_C")) {
                  
                  # Set bivariate locations
                  beta_mu_nd <- c(-100, -1, -1, -1, -1)
                  beta_mu_d  <- c(+100, 1, 1, 1, 1)
                  beta_SD_nd <- rep(0.25, 5)
                  beta_SD_d  <- rep(0.5, 5)
                  beta_corr  <- rep(0.0, 5)
                  ##
                  raw_scale_mu_nd <- raw_scale_mu_d  <- rep(0.125, 5)
                  raw_scale_SD_nd <- raw_scale_SD_d  <- rep(0.25, 5)
                  raw_scale_corr <- rep(0.0, 5)
                  ##
                  ## GAD-2 location parameters
                  ##
                  beta_mu_nd[2] <- true_mean_GAD_2_params_list$Ord_bivariate$beta_mu[1]
                  beta_SD_nd[2] <- true_mean_GAD_2_params_list$Ord_bivariate$beta_SD[1]
                  beta_mu_d[2]  <- true_mean_GAD_2_params_list$Ord_bivariate$beta_mu[2]
                  beta_SD_d[2]  <- true_mean_GAD_2_params_list$Ord_bivariate$beta_SD[2]
                  beta_corr[2]  <- true_mean_GAD_2_params_list$Ord_bivariate$beta_corr
                  ##
                  ## HADS location parameters
                  ##
                  beta_mu_nd[4] <- true_mean_HADS_params_list$Ord_bivariate$beta_mu[1]
                  beta_SD_nd[4] <- true_mean_HADS_params_list$Ord_bivariate$beta_SD[1]
                  beta_mu_d[4]  <- true_mean_HADS_params_list$Ord_bivariate$beta_mu[2]
                  beta_SD_d[4]  <- true_mean_HADS_params_list$Ord_bivariate$beta_SD[2]
                  beta_corr[4]  <- true_mean_HADS_params_list$Ord_bivariate$beta_corr
                  ##
                  ## BAI location parameters
                  ##
                  beta_mu_nd[3] <- true_mean_BAI_params_list$Ord_bivariate$beta_mu[1]
                  beta_SD_nd[3] <- true_mean_BAI_params_list$Ord_bivariate$beta_SD[1]
                  beta_mu_d[3]  <- true_mean_BAI_params_list$Ord_bivariate$beta_mu[2]
                  beta_SD_d[3]  <- true_mean_BAI_params_list$Ord_bivariate$beta_SD[2]
                  beta_corr[3]  <- true_mean_BAI_params_list$Ord_bivariate$beta_corr
                  ##
                  ## GAD-7 location parameters
                  ##
                  beta_mu_nd[5] <- true_mean_GAD_7_params_list$Ord_bivariate$beta_mu[1]
                  beta_SD_nd[5] <- true_mean_GAD_7_params_list$Ord_bivariate$beta_SD[1]
                  beta_mu_d[5]  <- true_mean_GAD_7_params_list$Ord_bivariate$beta_mu[2]
                  beta_SD_d[5]  <- true_mean_GAD_7_params_list$Ord_bivariate$beta_SD[2]
                  beta_corr[5]  <- true_mean_GAD_7_params_list$Ord_bivariate$beta_corr
              
          } else if (DGM == "Jones") { 
            
                  # Set bivariate locations
                  beta_mu_nd <- c(-100, -1, -1, -1, -1)
                  beta_mu_d  <- c(+100, 1, 1, 1, 1)
                  beta_SD_nd <- rep(0.25, 5)
                  beta_SD_d  <- rep(0.5, 5)
                  beta_corr  <- rep(0.0, 5)
                  ##
                  raw_scale_mu_nd <- raw_scale_mu_d  <- rep(0.125, 5)
                  raw_scale_SD_nd <- raw_scale_SD_d  <- rep(0.25, 5)
                  raw_scale_corr <- rep(0.0, 5)
                  ##
                  ## GAD-2 location parameters
                  ##
                  beta_mu_nd[2] <- true_mean_GAD_2_params_list$Jones$beta_mu[1]
                  beta_SD_nd[2] <- true_mean_GAD_2_params_list$Jones$beta_SD[1]
                  beta_mu_d[2]  <- true_mean_GAD_2_params_list$Jones$beta_mu[2]
                  beta_SD_d[2]  <- true_mean_GAD_2_params_list$Jones$beta_SD[2]
                  beta_corr[2]  <- true_mean_GAD_2_params_list$Jones$beta_corr
                  ##
                  raw_scale_mu_nd[2]        <- true_mean_GAD_2_params_list$Jones$raw_scale_mu[1]
                  raw_scale_SD_nd[2]        <- true_mean_GAD_2_params_list$Jones$raw_scale_SD[1]
                  raw_scale_mu_d[2]         <- true_mean_GAD_2_params_list$Jones$raw_scale_mu[2]
                  raw_scale_SD_d[2]         <- true_mean_GAD_2_params_list$Jones$raw_scale_SD[2]
                  raw_scale_corr[2]         <- true_mean_GAD_2_params_list$Jones$raw_scale_corr
                  ##
                  ## HADS location parameters
                  ##
                  beta_mu_nd[4] <- true_mean_HADS_params_list$Jones$beta_mu[1]
                  beta_SD_nd[4] <- true_mean_HADS_params_list$Jones$beta_SD[1]
                  beta_mu_d[4]  <- true_mean_HADS_params_list$Jones$beta_mu[2]
                  beta_SD_d[4]  <- true_mean_HADS_params_list$Jones$beta_SD[2]
                  beta_corr[4]  <- true_mean_HADS_params_list$Jones$beta_corr
                  ##
                  raw_scale_mu_nd[4]        <- true_mean_HADS_params_list$Jones$raw_scale_mu[1]
                  raw_scale_SD_nd[4]        <- true_mean_HADS_params_list$Jones$raw_scale_SD[1]
                  raw_scale_mu_d[4]         <- true_mean_HADS_params_list$Jones$raw_scale_mu[2]
                  raw_scale_SD_d[4]         <- true_mean_HADS_params_list$Jones$raw_scale_SD[2]
                  raw_scale_corr[4]         <- true_mean_HADS_params_list$Jones$raw_scale_corr
                  ##
                  ## BAI location parameters
                  ##
                  beta_mu_nd[3] <- true_mean_BAI_params_list$Jones$beta_mu[1]
                  beta_SD_nd[3] <- true_mean_BAI_params_list$Jones$beta_SD[1]
                  beta_mu_d[3]  <- true_mean_BAI_params_list$Jones$beta_mu[2]
                  beta_SD_d[3]  <- true_mean_BAI_params_list$Jones$beta_SD[2]
                  beta_corr[3]  <- true_mean_BAI_params_list$Jones$beta_corr
                  ##
                  raw_scale_mu_nd[3]        <- true_mean_BAI_params_list$Jones$raw_scale_mu[1]
                  raw_scale_SD_nd[3]        <- true_mean_BAI_params_list$Jones$raw_scale_SD[1]
                  raw_scale_mu_d[3]         <- true_mean_BAI_params_list$Jones$raw_scale_mu[2]
                  raw_scale_SD_d[3]         <- true_mean_BAI_params_list$Jones$raw_scale_SD[2]
                  raw_scale_corr[3]         <- true_mean_BAI_params_list$Jones$raw_scale_corr
            
          } else if (DGM == "R&G_random_C") { 
            
                  # Set bivariate locations
                  location_mu  <- c(0, 0, 0, 0, 0)
                  location_SD <- rep(0.25, 5)
                  location_corr  <- rep(0.0, 5)
                  ##
                  raw_scale_mu  <- rep(0.0, 5)
                  raw_scale_SD  <- rep(0.125, 5)
                  raw_scale_corr <- rep(0.0, 5)
                  ##
                  mult_nd <- -1.0
                  mult_d  <- +1.0
                  ##
                  ## GAD-2 parameters
                  ##
                  location_mu[2]   <- true_mean_GAD_2_params_list$Ord_HSROC$beta_mu
                  location_SD[2]   <- true_mean_GAD_2_params_list$Ord_HSROC$beta_SD
                  location_corr[2] <- 0.0
                  ##
                  raw_scale_mu[2]   <- true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_mu
                  raw_scale_SD[2]   <- true_mean_GAD_2_params_list$Ord_HSROC$raw_scale_SD
                  raw_scale_corr[2] <- 0.0
                  ##
                  ## HADS parameters
                  ##
                  location_mu[4]   <- true_mean_HADS_params_list$Ord_HSROC$beta_mu
                  location_SD[4]   <- true_mean_HADS_params_list$Ord_HSROC$beta_SD
                  location_corr[4] <- 0.0
                  ##
                  raw_scale_mu[4]   <- true_mean_HADS_params_list$Ord_HSROC$raw_scale_mu
                  raw_scale_SD[4]   <- true_mean_HADS_params_list$Ord_HSROC$raw_scale_SD
                  raw_scale_corr[4] <- 0.0
                  ##
                  ## BAI parameters
                  ##
                  location_mu[3]   <- true_mean_BAI_params_list$Ord_HSROC$beta_mu
                  location_SD[3]   <- true_mean_BAI_params_list$Ord_HSROC$beta_SD
                  location_corr[3] <- 0.0
                  ##
                  raw_scale_mu[3]   <- true_mean_BAI_params_list$Ord_HSROC$raw_scale_mu
                  raw_scale_SD[3]   <- true_mean_BAI_params_list$Ord_HSROC$raw_scale_SD
                  raw_scale_corr[3] <- 0.0
            
          }
          ##
          ## Set up threshold parameters
          ##
          n_thr_per_test <- c(1, 
                              6,  ## GAD-2
                              63, ## BAI
                              21, ## HADS
                              21) ## GAD-7
          ##
          max_threshold_across_all_tests <- max(n_thr_per_test)
          Mean_of_thr_for_all_tests_array_per_study_nd <- array(NA, dim = c(n_studies, n_tests, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_per_study_d  <- array(NA, dim = c(n_studies, n_tests, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_nd <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_d  <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
          # SD_of_thr_for_all_tests_array_nd <- array( 0.001, dim = c(n_tests, max_threshold_across_all_tests))
          # SD_of_thr_for_all_tests_array_d <- array( 0.001, dim = c(n_tests, max_threshold_across_all_tests))
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
 
          # ## sd_thr_jitter_nd <- 0.001
          # sd_thr_jitter_nd <- 0.01
          # sd_thr_jitter_d  <- 2.0*sd_thr_jitter_nd
          
          s <- 1
          ##
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
                    ##
                    ##
                    ## ---- GAD-2 (test #2):
                    ##
                    t <- 2
                    ##
                    if (DGM == "ord_biv_random_C") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd,
                                                                                            use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_d,
                                                                                            use_probit_link =  TRUE)
                      
                    } else if (DGM == "ord_biv_fixed_C") { 
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_GAD_2_params_list$Ord_bivariate$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_GAD_2_params_list$Ord_bivariate$C_d
                      
                    } else if (DGM == "Jones") { 
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_GAD_2_params_list$Jones$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_GAD_2_params_list$Jones$C_d
                      
                    } else if (DGM == "R&G_random_C") { 
                      
                          thresholds <- generate_ordered_thresholds_induced_dirichlet( alpha_vec = true_mean_GAD_2_params_list$Ord_HSROC$alpha,
                                                                                       use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- thresholds
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- thresholds
                          
                    }
                    ##
                    ##
                    ## ---- BAI (test #3):
                    ##
                    t <- 3
                    ##
                    if (DGM == "ord_biv_random_C") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_nd,
                                                                                            use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_d,
                                                                                            use_probit_link =  TRUE)
                          
                    } else if (DGM == "ord_biv_fixed_C") { 
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_BAI_params_list$Ord_bivariate$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_BAI_params_list$Ord_bivariate$C_d
                      
                    } else if (DGM == "Jones") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_BAI_params_list$Jones$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_BAI_params_list$Jones$C_d
                      
                    } else if (DGM == "R&G_random_C") { 
                      
                          thresholds <- generate_ordered_thresholds_induced_dirichlet( alpha_vec = true_mean_BAI_params_list$Ord_HSROC$alpha,
                                                                                       use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- thresholds
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- thresholds
                          
                    }
                    ##
                    ## ---- HADS (test #4):
                    ##
                    t <- 4
                    ##
                    if (DGM == "ord_biv_random_C") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_nd,
                                                                                            use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_d,
                                                                                            use_probit_link =  TRUE)
                      
                    } else if (DGM == "ord_biv_fixed_C") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_HADS_params_list$Ord_bivariate$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_HADS_params_list$Ord_bivariate$C_d
                      
                    } else if (DGM == "Jones") {
                      
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_HADS_params_list$Jones$C_nd
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_HADS_params_list$Jones$C_d
                      
                    } else if (DGM == "R&G_random_C") { 
                      
                          thresholds <- generate_ordered_thresholds_induced_dirichlet( alpha_vec = true_mean_HADS_params_list$Ord_HSROC$alpha,
                                                                                       use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- thresholds
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- thresholds
                          
                    }
                    ##
                    ## ---- GAD-7 (test #5):
                    ##
                    t <- 5
                    ##
                    # if (DGM == "ord_biv_random_C") {
                          ##
                          thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_nd,
                                                                                            use_probit_link =  TRUE)
                          thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet(
                                                                                            alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_d,
                                                                                            use_probit_link =  TRUE)
                      
                    # } else if (DGM == "ord_biv_fixed_C") {
                    #       
                    #       thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_GAD_7_params_list$Ord_bivariate$C_nd
                    #       thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_GAD_7_params_list$Ord_bivariate$C_d
                    #   
                    # } else if (DGM == "Jones") {
                    #   
                    #       thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_mean_GAD_7_params_list$Jones$C_nd
                    #       thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_mean_GAD_7_params_list$Jones$C_d
                    #   
                    # }
                    # ##
                    for (t in 2:n_tests) {
                      
                            Mean_of_thr_for_all_tests_array_per_study_nd[s, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]]
                            Mean_of_thr_for_all_tests_array_per_study_d[s, t, 1:n_thr_per_test[t]]  <- thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]
                            ##
                            # Assign thresholds to diseased and non-diseased subjects
                            diseased_indices    <- which(d_ind == 1)
                            nondiseased_indices <- which(d_ind == 0)
                            ##
                            for (n in diseased_indices) {
                              thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]
                            }
                            for (n in nondiseased_indices) {
                              thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]]
                            }
                    }
                    ##
                    # Generate study-specific parameters
                    if (DGM %in% c("ord_biv_random_C", "ord_biv_fixed_C")) {
                      
                            Sigma_bs_list <- list()
                            ##
                            for (t in 1:n_tests) {
                                  Omega_bs <- matrix(1.0, nrow = 2, ncol = 2)
                                  Sigma_bs <- matrix(NA, nrow = 2, ncol = 2)
                                  ##
                                  Omega_bs[1, 2] <- beta_corr[t]
                                  Omega_bs[2, 1] <- beta_corr[t]
                                  ##
                                  Sigma_bs[1, 1] <- beta_SD_nd[t]^2
                                  Sigma_bs[2, 2] <- beta_SD_d[t]^2
                                  Sigma_bs[1, 2] <- beta_corr[t] * beta_SD_nd[t] * beta_SD_d[t]
                                  Sigma_bs[2, 1] <- Sigma_bs[1, 2]
                                  Sigma_bs_list[[t]] <- Sigma_bs
                            }
                            ##
                            location_nd_study_s <- c()
                            location_d_study_s <- c()
                            ##
                            for (t in 1:n_tests) {
                              locations <- LaplacesDemon::rmvn(n = 1, 
                                                               mu = c(beta_mu_nd[t], beta_mu_d[t]), 
                                                               Sigma = Sigma_bs_list[[t]])
                              location_nd_study_s[t] <- locations[1]
                              location_d_study_s[t]  <- locations[2]
                            }
                            ##
                            scale_nd_study_s <- rep(1, n_tests)
                            scale_d_study_s  <- rep(1, n_tests)
                          
                    } else if (DGM == "R&G_random_C") {
                      
                            location_nd_study_s <- location_d_study_s <- c()
                            scale_nd_study_s <- scale_d_study_s <- c()
                            HSROC_location_study_s <- c()
                            HSROC_raw_scale_study_s <- c()
                            
                            for (t in 1:n_tests) {
                                  HSROC_location_study_s[t] <- rnorm(n = 1, 
                                                                     mean = location_mu[t], 
                                                                     sd = location_SD[t])
                                  ## ---- To match Stan model, need to FLIP signs but ONLY for locations - NOT for scales!
                                  location_nd_study_s[t] <- -1.0 * mult_nd * HSROC_location_study_s[t] 
                                  location_d_study_s[t]  <- -1.0 * mult_d  * HSROC_location_study_s[t]
                                  ##
                                  HSROC_raw_scale_study_s[t] <- rnorm(n = 1, 
                                                                   mean = raw_scale_mu[t], 
                                                                   sd = raw_scale_SD[t])
                                  scale_nd_study_s[t] <- exp(mult_nd * HSROC_raw_scale_study_s[t])
                                  scale_d_study_s[t]  <- exp(mult_d  * HSROC_raw_scale_study_s[t])
                            }
                      
                    } else if (DGM == "Jones") {
                      
                            Sigma_bs_list <- list()
                            ##
                            for (t in 1:n_tests) {
                              Omega_bs <- matrix(1.0, nrow = 2, ncol = 2)
                              Sigma_bs <- matrix(NA, nrow = 2, ncol = 2)
                              ##
                              Omega_bs[1, 2] <- beta_corr[t]
                              Omega_bs[2, 1] <- beta_corr[t]
                              ##
                              Sigma_bs[1, 1] <- beta_SD_nd[t]^2
                              Sigma_bs[2, 2] <- beta_SD_d[t]^2
                              Sigma_bs[1, 2] <- beta_corr[t] * beta_SD_nd[t] * beta_SD_d[t]
                              Sigma_bs[2, 1] <- Sigma_bs[1, 2]
                              Sigma_bs_list[[t]] <- Sigma_bs
                            }
                            ##
                            location_nd_study_s <- c()
                            location_d_study_s <- c()
                            ##
                            for (t in 1:n_tests) {
                              locations <- LaplacesDemon::rmvn(n = 1, 
                                                               mu = c(beta_mu_nd[t], beta_mu_d[t]), 
                                                               Sigma = Sigma_bs_list[[t]])
                              location_nd_study_s[t] <- locations[1]
                              location_d_study_s[t]  <- locations[2]
                            }
                            ##
                            raw_scale_Sigma_bs_list <- list()
                            ##
                            for (t in 1:n_tests) {
                              raw_scale_Omega_bs <- matrix(1.0, nrow = 2, ncol = 2)
                              raw_scale_Sigma_bs <- matrix(NA, nrow = 2, ncol = 2)
                              ##
                              raw_scale_Omega_bs[1, 2] <- raw_scale_corr[t]
                              raw_scale_Omega_bs[2, 1] <- raw_scale_corr[t]
                              ##
                              raw_scale_Sigma_bs[1, 1] <- raw_scale_SD_nd[t]^2
                              raw_scale_Sigma_bs[2, 2] <- raw_scale_SD_d[t]^2
                              raw_scale_Sigma_bs[1, 2] <- raw_scale_corr[t] * raw_scale_SD_nd[t] * raw_scale_SD_d[t]
                              raw_scale_Sigma_bs[2, 1] <- raw_scale_Sigma_bs[1, 2]
                              raw_scale_Sigma_bs_list[[t]] <- raw_scale_Sigma_bs
                            }
                            ## 
                            scale_nd_study_s <- scale_d_study_s <- c()
                            ##
                            for (t in 1:n_tests) {
                              raw_scales <- LaplacesDemon::rmvn(n = 1, 
                                                               mu = c(raw_scale_mu_nd[t], raw_scale_mu_d[t]), 
                                                               Sigma = raw_scale_Sigma_bs_list[[t]])
                              scale_nd_study_s[t] <- exp(raw_scales[1])
                              scale_d_study_s[t]  <- exp(raw_scales[2])
                            }
                      
                    }
                    
                    
                    # Set correlation parameters
                    if (assume_perfect_GS == TRUE) {
                      
                          location_nd_study_s[1] <- -100
                          location_d_study_s[1]  <- +100
                          ##
                          scale_nd_study_s[1] <- 1
                          scale_d_study_s[1]  <- 1
                          ##
                          rho1 <- 0.0
                        
                    } else {
                      
                          rho1 <- 0.2
                        
                    }
                    ##
                    ## Set up correlation matrices:
                    Omega_highly_varied <- matrix(c(1,     rho1,  rho1,     rho1,     rho1,
                                                    rho1,  1.00,  0,     0,     0,
                                                    rho1,  0,     1.00,  0,     0,
                                                    rho1,  0,     0,     1.00,  0,
                                                    rho1,  0,     0,     0,     1.00),
                                                  nrow = n_tests,
                                                  ncol = n_tests)
                    ##
                    {
                      Omega_d <- Omega_highly_varied
                      diag(Omega_d) <- rep(1, n_tests)
                      Omega_nd <- 0.5 * Omega_highly_varied
                      diag(Omega_nd) <- rep(1, n_tests)
                    }
                    ##
                    # Ensure correlation matrices are proper
                    Omega_nd <- as.matrix(Matrix::nearPD(Omega_nd, keepDiag = TRUE)$mat)
                    Omega_d  <- as.matrix(Matrix::nearPD(Omega_d, keepDiag = TRUE)$mat)
                    L_Omega_nd <- t(chol(Omega_nd))
                    L_Omega_d  <- t(chol(Omega_d))
                    
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
                            
                            # if (t == 2 && s == 1) {
                            #   print(paste("k =", k, "| Unique values in df_pos_y:", paste(sort(unique(df_pos_y)), collapse = ", ")))
                            #   print(paste("Count of df_pos_y >= (k-1):", sum(df_pos_y >= (k - 1)), "| Total positive cases:", n_pos))
                            # }
                            
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
          # ##
          # ## ----------------- Calculate true DGM values:
          # ##
          # true_DGM_Se <- true_DGM_Fp <- true_DGM_Sp <- list()
          # ##
          # ## ---- GAD-2 (test #2):
          # ##
          # t <- 2
          # if (DGM == "ord_biv_random_C") { 
          #   
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd)
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_d)
          #   
          # } else if (DGM == "ord_biv_fixed_C") {
          #   
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_GAD_2_params_list$Ord_bivariate$C_nd
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_GAD_2_params_list$Ord_bivariate$C_d
          #   
          # } else if (DGM == "Jones") {
          #   
          # }
          # ##
          # ## ---- BAI (test #3):
          # ##
          # t <- 3
          # if (DGM == "ord_biv_random_C") { 
          #   
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_nd)
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_d)
          #   
          # } else if (DGM == "ord_biv_fixed_C") {
          #   
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_BAI_params_list$Ord_bivariate$C_nd
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_BAI_params_list$Ord_bivariate$C_d
          #   
          # } else if (DGM == "Jones") {
          #   
          #   
          # }
          # 
          # ##
          # ## ---- HADS (test #4):
          # ##
          # t <- 4
          # if (DGM == "ord_biv_random_C") { 
          #   
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_nd)
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet_sample(
          #                                                                     n_samples = 1000,
          #                                                                     alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_d)
          #   
          # } else if (DGM == "ord_biv_fixed_C") {
          #       
          #       Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_HADS_params_list$Ord_bivariate$C_nd
          #       Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_HADS_params_list$Ord_bivariate$C_d
          #       
          # } else if (DGM == "Jones") {
          #   
          #   
          # }
          # ##
          # ## ---- GAD-7 (test #5):
          # ##
          # t <- 5
          # ##
          # Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
          #   alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_nd)
          # Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
          #   alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_d)
          # 
          
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
            Sp_per_study_ref = Sp_per_study_ref
            # true_DGM_Se = true_DGM_Se, 
            # true_DGM_Fp = true_DGM_Fp, 
            # true_DGM_Sp = true_DGM_Sp
          ))
  
}














 
