
 

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




# Use this in your simulation
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
          
          # Transform SDs from probability to cutpoint scale
          sigma_C <- SD_approx_ID_ord_prob_to_C_probit(mean_cutpoints, SD_prob_scale_vec)
          
          if (n_thr == 1) {
            return(rnorm(1, mean = mean_cutpoints[1], sd = sigma_C[1]))
          }
          
          # Create correlation matrix - higher correlation for adjacent cutpoints
          # This mimics the induced-Dirichlet structure
          cor_mat <- matrix(0, n_thr, n_thr)
          for (i in 1:n_thr) {
            for (j in 1:n_thr) {
              # Correlation decays with distance
              cor_mat[i,j] <- base_correlation^abs(i-j)
            }
          }
          
          # Covariance matrix
          cov_mat <- diag(sigma_C) %*% cor_mat %*% diag(sigma_C)
          
          # Generate correlated cutpoints
          cutpoints <- MASS::mvrnorm(1, mu = mean_cutpoints, Sigma = cov_mat)
          
          # Ensure ordering
          cutpoints <- sort(cutpoints)
          
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
generate_ordered_thresholds_induced_dirichlet_sample <- function(n_samples = 10000000,
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
                                   # # HSROC_locations_bs_het_SD = 0.5, 
                                   # # HSROC_raw_scales_bs_het_SD = 0.25, 
                                   # bivariate_locations_nd_bs_het_SD = 0.25, 
                                   # bivariate_locations_d_bs_het_SD = 0.5, 
                                   bs_het_C_nd_prob_scale, 
                                   bs_het_C_d_prob_scale
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
          
          # Initialize location parameters
          bivariate_locations_nd <- rep(-1, n_tests)
          bivariate_locations_d <- rep(+1, n_tests)
          if (assume_perfect_GS == 1) {
            bivariate_locations_nd[1] <- -100
            bivariate_locations_d[1]  <- +100
          }
          
          # Generate prevalence parameters
          true_Mean_probit_prev <- qnorm(true_Mean_prev)
          true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev, 
                                              sd = true_SD_probit_prev)
          true_prev_per_study <- pnorm(true_probit_prev_per_study)
          
          # Set up test parameters
          {

            ##
            ## ---- GAD-2 parameters: ----------------------------------------------------------------------------------------------------------------------
            ## (from fitting Xu-Rand-Kappa model to REAL GAD-2 MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            ##
            true_mean_GAD_2_params_list <- list()
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$beta_mu <- c(-0.735, 0.197)
            true_mean_GAD_2_params_list$Ord_bivariate$beta_SD <- c(0.341, 0.559)
            true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.172
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$C_nd <- c(-0.9755975, -0.292786, 0.358953, 0.689446, 1.106, 1.46094)
            true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale <-  c(0.0378236, 0.041786, 0.043734, 0.0322398, 0.0317818, 0.0241089)
            true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd           <-  c(15.3304, 20.5346, 23.54825, 10.74365, 10.2649, 5.50998, 6.860215)
            ##
            true_mean_GAD_2_params_list$Ord_bivariate$C_d <-   c(-1.48302, -0.88127, -0.134779, 0.273662, 0.785527, 1.23541)
            true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale <-    c(0.0561697, 0.063388, 0.0852722, 0.0722046, 0.0756352, 0.0602498)
            true_mean_GAD_2_params_list$Ord_bivariate$alpha_d           <-    c(1.97502, 2.659315, 5.66787, 3.665955, 4.08459, 2.35475, 2.787845)
            ##
            true_mean_GAD_2_params_list$Ord_bivariate
            ##
            ## ---- GAD-7 parameters. -----------------------------------------------------------------------------------------------------------------------
            ##
            true_mean_GAD_7_params_list <- list()
            ##
            true_mean_GAD_7_params_list$Ord_bivariate$beta_mu <- c(-1.41, -0.283)
            true_mean_GAD_7_params_list$Ord_bivariate$beta_SD <- c(0.488, 0.674)
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
                                                                              17.65125)
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
                                                                                  3.43053, 5.521015)
            ##
            true_mean_GAD_7_params_list$Ord_bivariate
            ##
            ## ---- HADS parameters: ----------------------------------------------------------------------------------------------------------------------
            ## (from fitting Xu-Rand-Kappa model to REAL HADS MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_HADS_params_list <- list()
            ##
            true_mean_HADS_params_list$Ord_bivariate$beta_mu <- c(-0.640, 0.128)
            true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.235,  0.392)
            true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.302
            ##
            true_mean_HADS_params_list$Ord_bivariate$C_nd <-    c(-2.104755, -1.669755, -1.32896, -1.03661, -0.768047, -0.525203, 
                                                                  -0.276447, -0.0324133, 0.2037715, 0.4196795, 0.599935, 0.794462, 
                                                                  0.974998, 1.159995, 1.372315, 1.545365, 1.75456, 1.94278, 2.09185, 
                                                                  2.315215, 2.537795)
            true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale <-      c(0.0071724, 0.0089856, 0.0107984, 0.0122948, 0.0134223, 0.0139609, 
                                                                                  0.0150858, 0.0153606, 0.015128, 0.0142778, 0.0125449, 0.0124485, 
                                                                                  0.0111248, 0.0104644, 0.0099706, 0.0078793, 0.0074393, 0.0058618, 
                                                                                  0.0045857, 0.0041258, 0.0032296)
            true_mean_HADS_params_list$Ord_bivariate$alpha_nd         <-       c(6.80893, 10.8174, 15.9235, 21.05695, 25.51105, 27.81725, 33.0955, 
                                                                                 34.55905, 33.4195, 29.4049, 22.24045, 21.9021, 17.2082, 15.12155, 
                                                                                 13.60915, 8.403455, 7.452145, 4.61551, 2.79908, 2.250785, 1.35282, 
                                                                                 3.59519)
            ##
            true_mean_HADS_params_list$Ord_bivariate$C_d <-     c(-1.969345, -1.892385, -1.674995, -1.539365, -1.38404, -1.189135, 
                                                                  -0.94256, -0.736189, -0.501034, -0.2186955, 0.021306, 0.2022715, 
                                                                  0.430632, 0.6941875, 0.9345995, 1.16245, 1.328875, 1.470575, 
                                                                  1.610625, 1.826915, 1.944015)
            true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale <-     c(0.0093823, 0.0041046, 0.0076414, 0.0069644, 0.0083598, 0.0105113, 
                                                                                0.013535, 0.0137341, 0.0158818, 0.0182312, 0.0174599, 0.0153124, 
                                                                                0.0166366, 0.0168667, 0.0148842, 0.0130484, 0.0099909, 0.008245, 
                                                                                0.0074475, 0.0075942, 0.0047928)
            true_mean_HADS_params_list$Ord_bivariate$alpha_d           <-       c(7.05851, 1.270245, 4.69393, 3.93293, 5.774175, 8.96514, 15.0838, 
                                                                                  15.38695, 20.98205, 28.36065, 25.7017, 19.4406, 23.16915, 24.05625, 
                                                                                  18.4379, 14.0733, 8.107345, 5.52107, 4.50281, 4.663425, 1.740885, 
                                                                                  8.888425)
            ##
            true_mean_HADS_params_list$Ord_bivariate
            ##
            ## ---- BAI parameters: ----------------------------------------------------------------------------------------------------------------------
            ## (from fitting Xu-Rand-Kappa model to BAI HADS MA Klaus data (w/o imputing missing threhsolds - "-1"'s).)
            ##
            true_mean_BAI_params_list <- list()
            ##
            true_mean_BAI_params_list$Ord_bivariate$beta_mu <- c(-1.514, -0.3478)
            true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.2567,  0.3682)
            true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.005574
            ##
            true_mean_BAI_params_list$Ord_bivariate$C_nd <-    c(-2.732245, -2.387375, -2.131155, -1.9783, -1.82696, -1.695705, 
                                                                  -1.54511, -1.43848, -1.30658, -1.22135, -1.13454, -1.047695, 
                                                                  -0.9358325, -0.8690215, -0.7929205, -0.6970965, -0.6215295, -0.548447, 
                                                                  -0.4977135, -0.4503225, -0.3935105, -0.3476205, -0.2838145, -0.2169855, 
                                                                  -0.164063, -0.1078915, -0.0683982, -0.0317889, 0.0057544, 0.0355495, 
                                                                  0.1026155, 0.1354665, 0.156525, 0.195646, 0.2264185, 0.289278, 
                                                                  0.3289735, 0.4163955, 0.437654, 0.4579155, 0.483268, 0.543683, 
                                                                  0.54685, 0.6149605, 0.702533, 0.7322465, 0.8139215, 0.8669755, 
                                                                  0.864729, NA, 0.945261, 0.960884, 1.05953, 1.06125, 1.147085, 
                                                                  1.170595, 1.22346, NA, NA, 1.39561, NA, 1.65157, 1.743205)
            true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale <-      c(0.0022258, 0.0027771, 0.0031775, 0.003317, 0.0039419, 0.0039193, 
                                                                                  0.004729, 0.0043634, 0.0047205, 0.0048735, 0.0051855, 0.005064, 
                                                                                  0.0060246, 0.0049014, 0.0053449, 0.0062395, 0.0056115, 0.0057743, 
                                                                                  0.0051958, 0.0047255, 0.0052106, 0.0047283, 0.0056988, 0.0058049, 
                                                                                  0.0052882, 0.0051256, 0.0047954, 0.0042163, 0.0045332, 0.0042103, 
                                                                                  0.0053407, 0.0043758, 0.0032516, 0.0043956, 0.00397, 0.0053348, 
                                                                                  0.00517, 0.005898, 0.0034738, 0.0035204, 0.003252, 0.0045065, 
                                                                                  0.0033919, 0.0047635, 0.005488, 0.0037337, 0.0054154, 0.0028868, 
                                                                                  0.0028825, 0.0031611, 0.0031779, 0.0031876, 0.0045017, 0.0031124, 
                                                                                  0.0032843, 0.0033042, 0.00389, 0.0032669, 0.0031819, 0.0032341, 
                                                                                  0.0041153, 0.0039976, 0.0030255)
            true_mean_BAI_params_list$Ord_bivariate$alpha_nd         <-       c(2.45424, 3.80029, 5.00511, 5.489855, 7.745875, 7.67605, 11.1866, 
                                                                                 9.504145, 11.21095, 11.8736, 13.32125, 12.88525, 18.2642, 12.0509, 
                                                                                 14.4168, 19.6752, 15.92475, 16.74555, 13.59255, 11.30405, 13.6974, 
                                                                                 11.3949, 16.347, 17.1107, 14.12495, 13.37935, 11.4802, 8.81512, 
                                                                                 10.425, 8.9104, 14.3322, 9.8222, 5.27282, 9.798375, 7.925145, 
                                                                                 14.51595, 13.5244, 17.55315, 5.958355, 6.22533, 5.194005, 10.3437, 
                                                                                 5.693255, 11.47655, 15.37945, 6.96457, 14.9418, 4.02763, 3.97505, 
                                                                                 4.998775, 4.96461, 4.91987, 10.3152, 4.765575, 5.29702, 5.409755, 
                                                                                 7.550835, 5.157935, 4.94576, 5.101375, 8.279975, 8.05832, 4.51388, 
                                                                                 26.69655)
            ##
            true_mean_BAI_params_list$Ord_bivariate$C_d <-    c(-2.444395, -2.334045, -2.19856, -2.09331, -1.89484, -1.79838, 
                                                                 -1.721745, -1.61288, -1.53503, -1.445505, -1.40317, -1.28093, 
                                                                 -1.21059, -1.140995, -1.028165, -0.9017905, -0.834665, -0.70046, 
                                                                 -0.6317085, -0.571049, -0.516496, -0.4639385, -0.3725345, -0.317404, 
                                                                 -0.24743, -0.2175055, -0.150963, -0.0889879, -0.0431702, -0.0292832, 
                                                                 0.0655418, 0.139562, 0.1779475, 0.2172455, 0.2560985, 0.331448, 
                                                                 0.344908, 0.4035175, 0.475689, 0.513673, 0.584926, 0.6877265, 
                                                                 0.717862, 0.7596745, 0.804047, 0.8609135, 0.926506, 0.9966605, 
                                                                 1.052645, NA, 1.14518, 1.17595, 1.2203, 1.272135, 1.349195, 1.420015, 
                                                                 1.46297, NA, NA, 1.64526, NA, 1.804715, 2.044935)
            true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale <-    c(0.0021626, 0.0010989, 0.0014715, 0.0014889, 0.0024731, 0.0019332, 
                                                                               0.0019524, 0.0024103, 0.0021557, 0.0025875, 0.0018551, 0.0033631, 
                                                                               0.0027134, 0.0027322, 0.0038526, 0.0042342, 0.0033131, 0.0048345, 
                                                                               0.0035799, 0.0033566, 0.0032567, 0.0032841, 0.0043462, 0.0033693, 
                                                                               0.0038492, 0.0025173, 0.003852, 0.0037676, 0.0030316, 0.002122, 
                                                                               0.004543, 0.0040062, 0.0027845, 0.0027904, 0.002795, 0.0039004, 
                                                                               0.0019638, 0.0033492, 0.0038479, 0.0026635, 0.0035557, 0.004181, 
                                                                               0.0023594, 0.0023861, 0.0024149, 0.0029876, 0.003015, 0.0027166, 
                                                                               0.0027265, 0.0021096, 0.002106, 0.0019114, 0.0018586, 0.0022434, 
                                                                               0.0022594, 0.0024407, 0.0017831, 0.0018096, 0.0018111, 0.0017819, 
                                                                               0.0017882, 0.0017825, 0.002642)
            true_mean_BAI_params_list$Ord_bivariate$alpha_d           <-       c(11.8051, 2.807205, 5.439405, 5.486265, 15.6018, 9.42304, 9.71292, 
                                                                                  14.93015, 11.7485, 17.2174, 8.438585, 28.7092, 18.18885, 19.0014, 
                                                                                  37.7393, 45.7156, 27.78875, 59.584, 32.69305, 28.70255, 27.04855, 
                                                                                  27.35105, 48.0756, 29.16545, 37.90675, 15.93665, 37.776, 35.9761, 
                                                                                  23.57735, 10.8956, 53.23125, 41.3796, 19.69095, 19.66845, 20.0654, 
                                                                                  39.447, 9.469035, 28.359, 37.93865, 18.40605, 32.9332, 44.97685, 
                                                                                  13.879, 14.4692, 14.581, 23.5568, 23.43985, 18.25095, 18.3415, 
                                                                                  10.80625, 10.9846, 8.70559, 8.45787, 12.9152, 12.86445, 15.1903, 
                                                                                  7.95394, 7.97458, 8.06929, 7.90479, 7.955495, 7.93187, 18.1678, 
                                                                                  30.70595)
            ##
            true_mean_BAI_params_list$Ord_bivariate
          }
          
          # Set model type and parameters
          bs_model <- "ord_bivariate"
          if (bs_model == "ord_HSROC") {
            # HSROC model code would go here
          } else {
            # Set bivariate locations
            bivariate_locations_nd <- c(-1, -1, -1, -1, -1)
            bivariate_locations_d  <- c(1, 1, 1, 1, 1)
            ##
            bivariate_SD_nd <- rep(0.25, 5)
            bivariate_SD_d  <- rep(0.5, 5)
            ##
            bivariate_corr <- rep(0.0, 5)
          }
          ##
          ## GAD-2 location parameters
          ##
          bivariate_locations_nd[2] <- true_mean_GAD_2_params_list$Ord_bivariate$beta_mu[1]
          bivariate_SD_nd[2]        <- true_mean_GAD_2_params_list$Ord_bivariate$beta_SD[1]
          bivariate_locations_d[2]  <- true_mean_GAD_2_params_list$Ord_bivariate$beta_mu[2]
          bivariate_SD_d[2]         <- true_mean_GAD_2_params_list$Ord_bivariate$beta_SD[2]
          bivariate_corr[2]         <- true_mean_GAD_2_params_list$Ord_bivariate$beta_corr
          ##
          ## HADS location parameters
          ##
          bivariate_locations_nd[4] <- true_mean_HADS_params_list$Ord_bivariate$beta_mu[1]
          bivariate_SD_nd[4]        <- true_mean_HADS_params_list$Ord_bivariate$beta_SD[1]
          bivariate_locations_d[4]  <- true_mean_HADS_params_list$Ord_bivariate$beta_mu[2]
          bivariate_SD_d[4]         <- true_mean_HADS_params_list$Ord_bivariate$beta_SD[2]
          bivariate_corr[4]         <- true_mean_HADS_params_list$Ord_bivariate$beta_corr
          ##
          ## BAI location parameters
          ##
          bivariate_locations_nd[3] <- true_mean_BAI_params_list$Ord_bivariate$beta_mu[1]
          bivariate_SD_nd[3]        <- true_mean_BAI_params_list$Ord_bivariate$beta_SD[1]
          bivariate_locations_d[3]  <- true_mean_BAI_params_list$Ord_bivariate$beta_mu[2]
          bivariate_SD_d[3]         <- true_mean_BAI_params_list$Ord_bivariate$beta_SD[2]
          bivariate_corr[3]         <- true_mean_BAI_params_list$Ord_bivariate$beta_corr
          ##
          ## GAD-7 location parameters
          ##
          bivariate_locations_nd[5] <- true_mean_GAD_7_params_list$Ord_bivariate$beta_mu[1]
          bivariate_SD_nd[5]        <- true_mean_GAD_7_params_list$Ord_bivariate$beta_SD[1]
          bivariate_locations_d[5]  <- true_mean_GAD_7_params_list$Ord_bivariate$beta_mu[2]
          bivariate_SD_d[5]         <- true_mean_GAD_7_params_list$Ord_bivariate$beta_SD[2]
          bivariate_corr[5]         <- true_mean_GAD_7_params_list$Ord_bivariate$beta_corr
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
 
          s <- 1
          ##
          ## ---- Between-study correlation (and var-cov) matrix (using same one for each test)
          ##
          Sigma_bs_list <- list()
          for (t in 1:n_tests) {
            if (bs_model == "ord_bivariate") {
                ##
                Omega_bs <- matrix(1.0, nrow = 2, ncol = 2)
                Sigma_bs <- matrix(NA, nrow = 2, ncol = 2)
                ##
                Omega_bs[1, 2] <- bivariate_corr[t]
                Omega_bs[2, 1] <- bivariate_corr[t]
                ##
                Sigma_bs[1, 1] <- bivariate_SD_nd[t]^2
                Sigma_bs[2, 2] <- bivariate_SD_d[t]^2
                Sigma_bs[1, 2] <- bivariate_corr[t] * bivariate_SD_nd[t] * bivariate_SD_d[t]
                Sigma_bs[2, 1] <- Sigma_bs[1, 2]
                Sigma_bs_list[[t]] <- Sigma_bs
                ##
            }
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
                    ##
                    # Generate thresholds for all tests
                    # for (t in 2:n_tests) {
                    #         n_thr <- n_thr_per_test[t]
                    #         if (n_thr > 0) {
                    #           thr_for_all_tests_for_current_study_nd[t, 1:n_thr] <- generate_ordered_thresholds(
                    #             mean_thresholds = Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr],
                    #             sd_thresholds = SD_of_thr_for_all_tests_array_nd[t, 1:n_thr]
                    #           )
                    #           thr_for_all_tests_for_current_study_d[t, 1:n_thr] <- generate_ordered_thresholds(
                    #             mean_thresholds = Mean_of_thr_for_all_tests_array_d[t, 1:n_thr], 
                    #             sd_thresholds = SD_of_thr_for_all_tests_array_d[t, 1:n_thr]
                    #           )
                    #         }
                    # }
                    ##
                    ## ---- GAD-2 (test #2):
                    ##
                    t <- 2
                    # thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                         mean_cutpoints =  true_mean_GAD_2_params_list$Ord_bivariate$C_nd,
                    #                                                         SD_prob_scale_vec =  true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale,
                    #                                                         use_probit_link =  TRUE,
                    #                                                         base_correlation =  0.50)
                    # thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                       mean_cutpoints =  true_mean_GAD_2_params_list$Ord_bivariate$C_d,
                    #                                                       SD_prob_scale_vec =  true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale,
                    #                                                       use_probit_link =  TRUE,
                    #                                                       base_correlation =  0.50)
                    ##
                    thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd,
                      use_probit_link =  TRUE)
                    thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_d,
                      use_probit_link =  TRUE)
                    ##
                    ##
                    ## ---- BAI (test #3):
                    ##
                    t <- 3
                    thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                                                                      alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_nd,
                                                                                      use_probit_link =  TRUE)
                    thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet(
                                                                                      alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_d,
                                                                                      use_probit_link =  TRUE)
                    ##
                    ## ---- HADS (test #4):
                    ##
                    t <- 4
                    # thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                                     mean_cutpoints =  true_mean_HADS_params_list$Ord_bivariate$C_nd,
                    #                                                                     SD_prob_scale_vec =  true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale,
                    #                                                                     use_probit_link =  TRUE,
                    #                                                                     base_correlation =  0.50)
                    # thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                                   mean_cutpoints =  true_mean_HADS_params_list$Ord_bivariate$C_d,
                    #                                                                   SD_prob_scale_vec =  true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale,
                    #                                                                   use_probit_link =  TRUE,
                    #                                                                   base_correlation =  0.50)
                    ##
                    thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_nd,
                      use_probit_link =  TRUE)
                    thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_d,
                      use_probit_link =  TRUE)
                    ##
                    ## ---- GAD-7 (test #5):
                    ##
                    t <- 5
                    # thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                                     mean_cutpoints =  true_mean_GAD_7_params_list$Ord_bivariate$C_nd,
                    #                                                                     SD_prob_scale_vec =  true_mean_GAD_7_params_list$Ord_bivariate$C_nd_SD_prob_scale,
                    #                                                                     use_probit_link =  TRUE,
                    #                                                                     base_correlation =  0.50)
                    # thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_delta_method_correlated(
                    #                                                                   mean_cutpoints =  true_mean_GAD_7_params_list$Ord_bivariate$C_d,
                    #                                                                   SD_prob_scale_vec =  true_mean_GAD_7_params_list$Ord_bivariate$C_d_SD_prob_scale,
                    #                                                                   use_probit_link =  TRUE,
                    #                                                                   base_correlation =  0.50)
                    ##
                    thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_nd,
                      use_probit_link =  TRUE)
                    thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                      alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_d,
                      use_probit_link =  TRUE)
                    ##
                    for (t in 2:n_tests) {
                            ##
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
                        locations <- LaplacesDemon::rmvn(n = 1, mu = c(bivariate_locations_nd[t], bivariate_locations_d[t]), Sigma = Sigma_bs_list[[t]])
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
                    # ## Set up correlation matrices:
                    # Omega_highly_varied <- matrix(c(1,     rho1,  rho1,     rho1,     rho1,
                    #                                 rho1,  1.00,  0.50,     0.20,     0.10,
                    #                                 rho1,  0.50,  1.00,     0.40,     0.40,
                    #                                 rho1,  0.20,  0.40,     1.00,     0.70,
                    #                                 rho1,  0.10,  0.40,     0.70,     1.00),
                    #                               nrow = n_tests,
                    #                               ncol = n_tests)
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
          
          
          ##
          ## ----------------- Calculate true DGM values:
          ##
          true_DGM_Se <- true_DGM_Fp <- true_DGM_Sp <- list()
          ##
          ## ---- GAD-2 (test #2):
          ##
          t <- 2
          Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_GAD_2_params_list$Ord_bivariate$C_nd
          Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_GAD_2_params_list$Ord_bivariate$C_d
          # ##
          # Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample( n_samples = 10000000,
          #   alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_nd)
          # Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(
          #   alpha_vec =  true_mean_GAD_2_params_list$Ord_bivariate$alpha_d)
          # ##
          ##
          ## ---- BAI (test #3):
          ##
          t <- 3
          Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
                                                                        alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_nd)
          Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <-  generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
                                                                        alpha_vec =  true_mean_BAI_params_list$Ord_bivariate$alpha_d)
          ##
          ## ---- HADS (test #4):
          ##
          t <- 4
          Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
            alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_nd)
          Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
            alpha_vec =  true_mean_HADS_params_list$Ord_bivariate$alpha_d)
          # Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_HADS_params_list$Ord_bivariate$C_nd
          # Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_HADS_params_list$Ord_bivariate$C_d
          ##
          ## ---- GAD-7 (test #5):
          ##
          t <- 5
          # Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <-  true_mean_GAD_7_params_list$Ord_bivariate$C_nd
          # Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]]  <-  true_mean_GAD_7_params_list$Ord_bivariate$C_d
          Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
            alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_nd)
          Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet_sample(n_samples = 100,
            alpha_vec =  true_mean_GAD_7_params_list$Ord_bivariate$alpha_d)
          ##
          ## ---- Compute summary Se and Sp:
          ##
          # Mean_of_thr_for_all_tests_array_nd <- apply(Mean_of_thr_for_all_tests_array_per_study_nd, FUN = median, c(2, 3))
          # Mean_of_thr_for_all_tests_array_d <- apply(Mean_of_thr_for_all_tests_array_per_study_d, FUN = median, c(2, 3))
          ##
          for (t in 2:n_tests) {
            true_DGM_Se[[t - 1]] <- 1 - pnorm(Mean_of_thr_for_all_tests_array_d[t, 1:n_thr_per_test[t]] - bivariate_locations_d[t])
            true_DGM_Fp[[t - 1]] <- 1 - pnorm(Mean_of_thr_for_all_tests_array_nd[t, 1:n_thr_per_test[t]] - bivariate_locations_nd[t])
            true_DGM_Sp[[t - 1]] <- 1 - true_DGM_Fp[[t - 1]]
          }
          ##
          ## ---- GAD-2
          ##
          t <- 1 + 1
          true_DGM_Se[[t - 1]] <- c(95.563, 85.862, 63.12, 47.106, 27.594, 15.446)/100
          true_DGM_Fp[[t - 1]] <- c(59.611, 32.915, 13.738, 7.671, 3.31, 1.423)/100
          true_DGM_Sp[[t - 1]] <- c(40.389, 67.085, 86.262, 92.329, 96.69, 98.577)/100
          ##
          ## ---- HADS
          ##
          t <- 1 + 3
          true_DGM_Se[[t - 1]] <-         c(98.2, 97.834, 96.429, 95.211, 93.447, 90.584, 85.738, 80.568, 
                                            73.458, 63.496, 54.196, 46.979, 38.097, 28.518, 20.962, 14.992, 
                                            11.453, 8.952, 6.883, 4.428, 3.444)/100
          ##
          true_DGM_Fp[[t - 1]] <-   c(92.855, 84.843, 75.437, 65.381, 55.052, 45.399, 35.804, 27.157, 
                                      19.931, 14.467, 10.752, 7.575, 5.316, 3.598, 2.212, 1.44, 0.83, 
                                      0.491, 0.313, 0.155, 0.073)/100
          ##
          true_DGM_Sp[[t - 1]] <- c(7.145, 15.157, 24.563, 34.619, 44.948, 54.601, 64.196, 72.843, 
                                    80.069, 85.533, 89.248, 92.425, 94.684, 96.402, 97.788, 98.56, 
                                    99.17, 99.509, 99.687, 99.845, 99.927)/100
          ##
          ## ---- BAI:
          ##
          t <- 1 + 2
          true_DGM_Se[[t - 1]] <-    c(98.121, 97.614, 96.772, 95.936, 93.841, 92.638, 91.416, 89.622, 
                                       88.139, 86.268, 85.294, 82.32, 80.428, 78.427, 74.974, 70.799, 
                                       68.462, 63.491, 60.883, 58.539, 56.383, 54.35, 50.759, 48.56, 
                                       45.896, 44.695, 42.095, 39.667, 37.971, 37.413, 33.802, 31.229, 
                                       29.892, 28.623, 27.3, 24.923, 24.46, 22.646, 20.613, 19.532, 
                                       17.637, 15.099, 14.384, 13.499, 12.561, 11.41, 10.189, 8.997, 
                                       8.138, NA, 6.857, 6.437, 5.878, 5.314, 4.552, 3.904, 3.548, 
                                       NA, NA, 2.332, NA, 1.587, 0.86)/100
          ##
          true_DGM_Fp[[t - 1]] <-    c(88.795, 80.837, 73.135, 67.786, 62.114, 57.167, 51.226, 46.925, 
                                       41.544, 38.31, 35.111, 31.948, 28.056, 25.841, 23.451, 20.557, 
                                       18.526, 16.594, 15.377, 14.258, 13.021, 12.074, 10.828, 9.646, 
                                       8.761, 7.911, 7.347, 6.858, 6.377, 6.006, 5.26, 4.914, 4.696, 
                                       4.335, 4.058, 3.531, 3.247, 2.673, 2.538, 2.406, 2.264, 1.972, 
                                       1.943, 1.65, 1.329, 1.237, 0.998, 0.864, 0.868, NA, 0.696, 
                                       0.659, 0.502, 0.494, 0.387, 0.361, 0.308, NA, NA, 0.178, 
                                       NA, 0.075, 0.055)/100
          ##
          true_DGM_Sp[[t - 1]] <-    c(11.205, 19.163, 26.865, 32.213, 37.886, 42.833, 48.774, 53.075, 
                                       58.456, 61.69, 64.889, 68.052, 71.944, 74.159, 76.549, 79.443, 
                                       81.474, 83.406, 84.623, 85.742, 86.979, 87.926, 89.172, 90.354, 
                                       91.239, 92.089, 92.653, 93.142, 93.623, 93.994, 94.74, 95.086, 
                                       95.304, 95.665, 95.942, 96.469, 96.753, 97.327, 97.462, 97.594, 
                                       97.737, 98.028, 98.057, 98.35, 98.671, 98.763, 99.002, 99.136, 
                                       99.132, NA, 99.304, 99.341, 99.498, 99.506, 99.614, 99.639, 
                                       99.692, NA, NA, 99.822, NA, 99.925, 99.945)/100
          ##
          ## ---- GAD-7
          ##
          t <- 1 + 4
          true_DGM_Se[[t - 1]] <-            c(97.301, 95.945, 94.527, 92.268, 89.168, 84.232, 79.27, 71.388, 
                                               65.608, 59.321, 52.195, 44.815, 37.824, 31.994, 25.598, 21.024, 
                                               17.119, 12.366, 8.378, 4.91, 2.771)/100
          ##
          true_DGM_Fp[[t - 1]] <-  c(79.396, 65.052, 52.811, 42.426, 33.381, 26.387, 20.438, 14.645, 
                                     11.63, 9.01, 7.124, 5.612, 4.498, 3.445, 2.525, 1.805, 1.29, 
                                     1.003, 0.659, 0.422, 0.287)/100
          ##
          true_DGM_Sp[[t - 1]] <-    c(20.604, 34.948, 47.189, 57.574, 66.619, 73.613, 79.562, 85.355, 
                                       88.37, 90.99, 92.876, 94.388, 95.502, 96.555, 97.475, 98.195, 
                                       98.71, 98.997, 99.341, 99.578, 99.713)/100
          
          
          
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











 
