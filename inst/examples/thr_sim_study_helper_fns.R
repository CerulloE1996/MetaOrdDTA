



create_summary_df <- function(n_sims, 
                              min_k, 
                              max_k) {
  
        # Create vector of column names
        base_cols <- c(
          "seed", 
          "n_studies", 
          "model_parameterisation",
          "box_cox", 
          "softplus",
          "random_thresholds"
        )
        
        # Create empty vectors for threshold-specific column names
        se_cols <- character(max_k - min_k + 1)
        sp_cols <- character(max_k - min_k + 1)
        
        # Fill the threshold column names
        for (i in 1:length(se_cols)) {
          k <- min_k + i - 1
          se_cols[i] <- paste0("Se_mean_", k)
          sp_cols[i] <- paste0("Sp_mean_", k)
        }
        
        # Combine all column names
        all_cols <- c(base_cols, se_cols, sp_cols)
        
        # Create empty data frame with the right structure
        df <- data.frame(matrix(NA, nrow = n_sims, ncol = length(all_cols)))
        names(df) <- all_cols
        
        return(df)
  
}





# # Print summary:
# print_se_sp_stats <- function(k, 
#                               summary_stats, 
#                               true_DGM_Se,
#                               true_DGM_Sp, 
#                               output_dir) {
#   
#         cat(paste("--------------------------------------------------"))
#         cat(sprintf("\nSe (@thr = %d):\n", k))
#         cat(sprintf("  Mean: %s\n", formatC(summary_stats[[paste0("mean_Se_", k)]], digits=3, format="g")))
#         cat(sprintf("  SD: %s\n", formatC(summary_stats[[paste0("SD_Se_", k)]], digits=3, format="g")))
#         cat(sprintf("  Bias: %s\n", formatC(summary_stats[[paste0("mean_Se_", k)]] - true_DGM_Se[k], digits=3, format="g")))
#         ##
#         cat(sprintf("\nSp (@thr = %d):\n", k))
#         cat(sprintf("  Mean: %s\n", formatC(summary_stats[[paste0("mean_Sp_", k)]], digits=3, format="g")))
#         cat(sprintf("  SD: %s\n", formatC(summary_stats[[paste0("SD_Sp_", k)]], digits=3, format="g")))
#         cat(sprintf("  Bias: %s\n", formatC(summary_stats[[paste0("mean_Sp_", k)]] - true_DGM_Sp[k], digits=3, format="g")))
#         ##
#         cat(sprintf("\nResults saved to: %s\n", output_dir))
#   
# }



compute_average_bias <- function(summary_df, 
                                 true_DGM_Se, 
                                 true_DGM_Sp, 
                                 min_k,
                                 max_k) {
  
        # Initialize sums
        total_bias_Se <- 0
        total_bias_Sp <- 0
        total_SD_Se <- 0
        total_SD_Sp <- 0
        count <- 0
        ##
        # length_vec <- max_k - min_k + 1
        ##
        bias_Se_vec <- c()
        bias_Sp_vec <- c()
        SD_Se_vec <- c()
        SD_Sp_vec <- c()
        MCSE_bias_Se <- c()
        MCSE_bias_Sp <- c()
        
        N_sims <- nrow(summary_df)
        
        # summary_df[, paste0("Se_mean_", k)]
        # Compute bias for each k and add to total
        for (k in min_k:max_k) {
          
              ##
              ## ---- Se:
              ##
              bias_Se_vec[k] <- abs( mean(summary_df[, paste0("Se_mean_", k)], na.rm = TRUE) - true_DGM_Se[k])
              SD_Se_vec[k] <- sd(summary_df[, paste0("Se_mean_", k)], na.rm = TRUE)
              MCSE_bias_Se[k] <- sqrt(SD_Se_vec[k]^2/N_sims)
              ##
              ## ---- Sp
              ##
              bias_Sp_vec[k] <- abs( mean(summary_df[, paste0("Sp_mean_", k)], na.rm = TRUE) - true_DGM_Sp[k])
              SD_Sp_vec[k] <- sd(summary_df[, paste0("Sp_mean_", k)], na.rm = TRUE)
              MCSE_bias_Sp[k] <- sqrt(SD_Sp_vec[k]^2/N_sims)
              ##
              
        }
        
        # Calculate average biases
        mean_bias_Se <- mean(bias_Se_vec, na.rm = TRUE)
        max_bias_Se <- max(bias_Se_vec, na.rm = TRUE)
        mean_SD_Se   <- mean(SD_Se_vec, na.rm = TRUE)
        max_SD_Se   <- max(SD_Se_vec, na.rm = TRUE)
        mean_MCSE_bias_Se   <- mean(MCSE_bias_Se, na.rm = TRUE)
        max_MCSE_bias_Se   <- max(MCSE_bias_Se, na.rm = TRUE)
        ##
        mean_bias_Sp <- mean(bias_Sp_vec, na.rm = TRUE)
        max_bias_Sp <- max(bias_Sp_vec, na.rm = TRUE)
        mean_SD_Sp   <- mean(SD_Sp_vec, na.rm = TRUE)
        max_SD_Sp   <- max(SD_Sp_vec, na.rm = TRUE)
        mean_MCSE_bias_Sp   <- mean(MCSE_bias_Sp, na.rm = TRUE)
        max_MCSE_bias_Sp   <- max(MCSE_bias_Sp, na.rm = TRUE)
        
        # Return results
        return(list(
          bias_Se_vec = bias_Se_vec,
          SD_Se_vec = SD_Se_vec,
          ##
          bias_Sp_vec = bias_Sp_vec,
          SD_Sp_vec = SD_Sp_vec,
          ##
          mean_bias_Se = mean_bias_Se,
          max_bias_Se = max_bias_Se,
          mean_SD_Se = mean_SD_Se,
          max_SD_Se = max_SD_Se,
          mean_MCSE_bias_Se = mean_MCSE_bias_Se,
          max_MCSE_bias_Se = max_MCSE_bias_Se,
          ##
          mean_bias_Sp = mean_bias_Sp,
          max_bias_Sp = max_bias_Sp,
          mean_SD_Sp = mean_SD_Sp,
          max_SD_Sp = max_SD_Sp,
          mean_MCSE_bias_Sp = mean_MCSE_bias_Sp,
          max_MCSE_bias_Sp = max_MCSE_bias_Sp,
          ##
          summary_df = summary_df  # Return updated summary_df with Bias values
        ))
  
}




calc_min_sim_size <- function(SD_Se_or_Sp_vec, 
                              min_MCSE_pct,
                              min_k, 
                              max_k) {
  
        min_N_sim <- c()
        index <- 1
        
        for (k in min_k:max_k) {
          # Calculate minimum N required for this threshold
            min_N_sim[index] <- ceiling((SD_Se_or_Sp_vec[k]/min_MCSE_pct)^2)
            index <- index + 1
        }
        
        # Return both the vector and the maximum value
        return(list(
          ## min_N_sim = min_N_sim,
          min_required = min(min_N_sim, na.rm = TRUE),
          avg_required = mean(min_N_sim, na.rm = TRUE),
          max_required = max(min_N_sim, na.rm = TRUE)
        ))
  
}


















