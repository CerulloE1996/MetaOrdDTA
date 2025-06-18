



create_summary_df <- function(n_sims, 
                              min_k, 
                              max_k,
                              start_index = 1) {  # New parameter with default
  
          # Create vector of column names
          base_cols <- c(
            "seed", 
            # "sim_index",  
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
          
          # # Pre-populate the sim_index column
          # df$sim_index <- start_index:(start_index + n_sims - 1)
          
          return(df)
  
}





compute_simulation_metrics <- function(complete_results_list, 
                                       true_DGM_Se, 
                                       true_DGM_Sp,
                                       min_k,
                                       max_k) {
  
        # Filter out NULL results AND incomplete results (low ESS)
        valid_results <- complete_results_list[!sapply(complete_results_list, is.null)]
        
        # Only keep results that have credible intervals (i.e., successful runs)
        complete_results <- valid_results[sapply(valid_results, function(x) {
          !is.null(x$Se_lower_vec) && !is.null(x$Se_upper_vec) && 
            length(x$Se_lower_vec) > 0 && length(x$Se_upper_vec) > 0
        })]
        
        n_sims <- length(complete_results)
        n_excluded <- length(valid_results) - n_sims
        
        if(n_sims == 0) stop("No valid simulation results with adequate ESS found!")
        
        cat(sprintf("Using %d complete results (excluded %d with low ESS)\n", 
                    n_sims, n_excluded))
        
        # Initialize storage matrices
        coverage_Se <- coverage_Sp <- matrix(NA, n_sims, max_k)
        interval_width_Se <- interval_width_Sp <- matrix(NA, n_sims, max_k)
        bias_Se <- bias_Sp <- matrix(NA, n_sims, max_k)
        squared_error_Se <- squared_error_Sp <- matrix(NA, n_sims, max_k)
        
        # Store point estimates for SD calculation
        point_est_Se <- point_est_Sp <- matrix(NA, n_sims, max_k)
        
        # Loop through simulations (only complete ones now)
        for (i in 1:n_sims) {
          
          result <- complete_results[[i]]
          
          for (k in min_k:max_k) {
            
            # Store point estimates
            point_est_Se[i, k] <- result$Se_median_vec[k]
            point_est_Sp[i, k] <- result$Sp_median_vec[k]
            
            # Check coverage (is true value in CI?)
            coverage_Se[i, k] <- (true_DGM_Se[k] >= result$Se_lower_vec[k]) & 
              (true_DGM_Se[k] <= result$Se_upper_vec[k])
            coverage_Sp[i, k] <- (true_DGM_Sp[k] >= result$Sp_lower_vec[k]) & 
              (true_DGM_Sp[k] <= result$Sp_upper_vec[k])
            
            # Calculate interval width
            interval_width_Se[i, k] <- result$Se_upper_vec[k] - result$Se_lower_vec[k]
            interval_width_Sp[i, k] <- result$Sp_upper_vec[k] - result$Sp_lower_vec[k]
            
            # Calculate bias (using median as point estimate)
            bias_Se[i, k] <- result$Se_median_vec[k] - true_DGM_Se[k]
            bias_Sp[i, k] <- result$Sp_median_vec[k] - true_DGM_Sp[k]
            
            # Calculate squared error for RMSE
            squared_error_Se[i, k] <- (result$Se_median_vec[k] - true_DGM_Se[k])^2
            squared_error_Sp[i, k] <- (result$Sp_median_vec[k] - true_DGM_Sp[k])^2
            
          }
          
        }
        
        # Compute summary metrics for each threshold
        results_by_threshold <- list()
        ##
        bias_Se_vec <- c()
        bias_Sp_vec <- c()
        abs_bias_Se_vec <- c()
        abs_bias_Sp_vec <- c()
        MCSE_bias_Se_vec <- c()
        MCSE_bias_Sp_vec <- c()
        ##
        SD_Se_vec <- c()
        SD_Sp_vec <- c()
        ##
        RMSE_Se_vec <- c()
        RMSE_Sp_vec <- c()
        MCSE_RMSE_Se_vec <- c()
        MCSE_RMSE_Sp_vec <- c()
        
        
        for (k in min_k:max_k) {
          
          # Calculate bias metrics
          mean_bias_Se_k <- mean(bias_Se[, k], na.rm = TRUE)
          mean_bias_Sp_k <- mean(bias_Sp[, k], na.rm = TRUE)
          
          # Calculate SD of point estimates across simulations
          SD_Se_vec[k] <- sd(point_est_Se[, k], na.rm = TRUE)
          SD_Sp_vec[k] <- sd(point_est_Sp[, k], na.rm = TRUE)
          
          # Store for overall calculations
          bias_Se_vec[k] <- mean_bias_Se_k
          bias_Sp_vec[k] <- mean_bias_Sp_k
          abs_bias_Se_vec[k] <- abs(mean_bias_Se_k)
          abs_bias_Sp_vec[k] <- abs(mean_bias_Sp_k)
          
          # MCSE of bias
          MCSE_bias_Se_vec[k] <- SD_Se_vec[k] / sqrt(n_sims)
          MCSE_bias_Sp_vec[k] <- SD_Sp_vec[k] / sqrt(n_sims)
          
          # RMSE calculation
          RMSE_Se_k <- sqrt(mean(squared_error_Se[, k], na.rm = TRUE))
          RMSE_Sp_k <- sqrt(mean(squared_error_Sp[, k], na.rm = TRUE))
          
          RMSE_Se_vec[k] <- RMSE_Se_k
          RMSE_Sp_vec[k] <- RMSE_Sp_k
          
          # MCSE of RMSE using Delta method
          var_squared_error_Se <- var(squared_error_Se[, k], na.rm = TRUE)
          var_squared_error_Sp <- var(squared_error_Sp[, k], na.rm = TRUE)
          
          MCSE_RMSE_Se_k <- ifelse(RMSE_Se_k > 0, 
                                   sqrt(var_squared_error_Se) / (2 * sqrt(n_sims) * RMSE_Se_k),
                                   NA)
          MCSE_RMSE_Sp_k <- ifelse(RMSE_Sp_k > 0,
                                   sqrt(var_squared_error_Sp) / (2 * sqrt(n_sims) * RMSE_Sp_k),
                                   NA)
          
          MCSE_RMSE_Se_vec[k] <- MCSE_RMSE_Se_k
          MCSE_RMSE_Sp_vec[k] <- MCSE_RMSE_Sp_k
          
          results_by_threshold[[k]] <- list(
            # Coverage probability
            coverage_prob_Se = mean(coverage_Se[, k], na.rm = TRUE),
            coverage_prob_Sp = mean(coverage_Sp[, k], na.rm = TRUE),
            
            # Mean interval width
            mean_interval_width_Se = mean(interval_width_Se[, k], na.rm = TRUE),
            mean_interval_width_Sp = mean(interval_width_Sp[, k], na.rm = TRUE),
            
            # RMSE
            RMSE_Se = RMSE_Se_k,
            RMSE_Sp = RMSE_Sp_k,
            
            # MCSE of RMSE
            MCSE_RMSE_Se = MCSE_RMSE_Se_k,
            MCSE_RMSE_Sp = MCSE_RMSE_Sp_k,
            
            # Mean bias
            mean_bias_Se = mean_bias_Se_k,
            mean_bias_Sp = mean_bias_Sp_k,
            
            # Mean absolute bias
            mean_abs_bias_Se = abs(mean_bias_Se_k),
            mean_abs_bias_Sp = abs(mean_bias_Sp_k),
            
            # MCSE of bias
            MCSE_bias_Se = MCSE_bias_Se_vec[k],
            MCSE_bias_Sp = MCSE_bias_Sp_vec[k],
            
            # SD
            SD_Se = SD_Se_vec[k],
            SD_Sp = SD_Sp_vec[k]
          )
          
        }
        
        
        # Overall summaries (averaged across thresholds min_k:max_k)
        overall_results <- list(
          # Number of excluded results
          n_excluded = n_excluded,
          
          # Average RMSE across thresholds
          avg_RMSE_Se = mean(RMSE_Se_vec[min_k:max_k], na.rm = TRUE),
          avg_RMSE_Sp = mean(RMSE_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Maximum RMSE
          max_RMSE_Se = max(RMSE_Se_vec[min_k:max_k], na.rm = TRUE),
          max_RMSE_Sp = max(RMSE_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Average MCSE of RMSE
          avg_MCSE_RMSE_Se = mean(MCSE_RMSE_Se_vec[min_k:max_k], na.rm = TRUE),
          avg_MCSE_RMSE_Sp = mean(MCSE_RMSE_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # MCSE of max RMSE (using the MCSE at the threshold where max occurs)
          max_RMSE_Se_index = which.max(RMSE_Se_vec[min_k:max_k]) + min_k - 1,
          max_RMSE_Sp_index = which.max(RMSE_Sp_vec[min_k:max_k]) + min_k - 1,
          max_MCSE_RMSE_Se = MCSE_RMSE_Se_vec[which.max(RMSE_Se_vec[min_k:max_k]) + min_k - 1],
          max_MCSE_RMSE_Sp = MCSE_RMSE_Sp_vec[which.max(RMSE_Sp_vec[min_k:max_k]) + min_k - 1],
          
          # Average absolute bias across thresholds
          avg_abs_bias_Se = mean(abs_bias_Se_vec[min_k:max_k], na.rm = TRUE),
          avg_abs_bias_Sp = mean(abs_bias_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Maximum absolute bias
          max_abs_bias_Se = max(abs_bias_Se_vec[min_k:max_k], na.rm = TRUE),
          max_abs_bias_Sp = max(abs_bias_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Average MCSE of bias
          avg_MCSE_bias_Se = mean(MCSE_bias_Se_vec[min_k:max_k], na.rm = TRUE),
          avg_MCSE_bias_Sp = mean(MCSE_bias_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Maximum MCSE of bias
          max_MCSE_bias_Se = max(MCSE_bias_Se_vec[min_k:max_k], na.rm = TRUE),
          max_MCSE_bias_Sp = max(MCSE_bias_Sp_vec[min_k:max_k], na.rm = TRUE),
          
          # Average coverage across thresholds
          avg_coverage_prob_Se = mean(sapply(min_k:max_k, function(k) 
            results_by_threshold[[k]]$coverage_prob_Se), na.rm = TRUE),
          avg_coverage_prob_Sp = mean(sapply(min_k:max_k, function(k) 
            results_by_threshold[[k]]$coverage_prob_Sp), na.rm = TRUE),
          
          # Average interval width across thresholds
          avg_interval_width_Se = mean(sapply(min_k:max_k, function(k) 
            results_by_threshold[[k]]$mean_interval_width_Se), na.rm = TRUE),
          avg_interval_width_Sp = mean(sapply(min_k:max_k, function(k) 
            results_by_threshold[[k]]$mean_interval_width_Sp), na.rm = TRUE),
          
          # Average SD
          avg_SD_Se = mean(SD_Se_vec[min_k:max_k], na.rm = TRUE),
          avg_SD_Sp = mean(SD_Sp_vec[min_k:max_k], na.rm = TRUE),
          max_SD_Se = max(SD_Se_vec[min_k:max_k], na.rm = TRUE),
          max_SD_Sp = max(SD_Sp_vec[min_k:max_k], na.rm = TRUE)
        )
        
        return(list(
          by_threshold = results_by_threshold,
          overall = overall_results,
          n_sims = n_sims,
          bias_vectors = list(
            bias_Se_vec = bias_Se_vec,
            bias_Sp_vec = bias_Sp_vec,
            abs_bias_Se_vec = abs_bias_Se_vec,
            abs_bias_Sp_vec = abs_bias_Sp_vec,
            MCSE_bias_Se_vec = MCSE_bias_Se_vec,
            MCSE_bias_Sp_vec = MCSE_bias_Sp_vec,
            SD_Se_vec = SD_Se_vec,
            SD_Sp_vec = SD_Sp_vec,
            RMSE_Se_vec = RMSE_Se_vec,
            RMSE_Sp_vec = RMSE_Sp_vec,
            MCSE_RMSE_Se_vec = MCSE_RMSE_Se_vec,
            MCSE_RMSE_Sp_vec = MCSE_RMSE_Sp_vec
          )
        ))
  
}




# Updated summary table function
create_metrics_summary_table <- function(metrics, 
                                         min_k, 
                                         max_k) {
          
          # Create overall summary WITH RMSE first
          overall_df <- data.frame(
            Metric = c("Mean RMSE Se", "Mean RMSE Sp",
                       "Max RMSE Se", "Max RMSE Sp",
                       "Mean MCSE(RMSE) Se", "Mean MCSE(RMSE) Sp",
                       "Max MCSE(RMSE) Se", "Max MCSE(RMSE) Sp",
                       "Mean |Bias| Se", "Mean |Bias| Sp",
                       "Max |Bias| Se", "Max |Bias| Sp",
                       "Mean MCSE(Bias) Se", "Mean MCSE(Bias) Sp",
                       "Max MCSE(Bias) Se", "Max MCSE(Bias) Sp",
                       "Coverage Se", "Coverage Sp", 
                       "Mean Width Se", "Mean Width Sp"),
            Value = c(
              sprintf("%.3f", metrics$overall$avg_RMSE_Se),
              sprintf("%.3f", metrics$overall$avg_RMSE_Sp),
              sprintf("%.3f", metrics$overall$max_RMSE_Se),
              sprintf("%.3f", metrics$overall$max_RMSE_Sp),
              sprintf("%.4f", metrics$overall$avg_MCSE_RMSE_Se),
              sprintf("%.4f", metrics$overall$avg_MCSE_RMSE_Sp),
              sprintf("%.4f", metrics$overall$max_MCSE_RMSE_Se),
              sprintf("%.4f", metrics$overall$max_MCSE_RMSE_Sp),
              sprintf("%.3f", metrics$overall$avg_abs_bias_Se),
              sprintf("%.3f", metrics$overall$avg_abs_bias_Sp),
              sprintf("%.3f", metrics$overall$max_abs_bias_Se),
              sprintf("%.3f", metrics$overall$max_abs_bias_Sp),
              sprintf("%.4f", metrics$overall$avg_MCSE_bias_Se),
              sprintf("%.4f", metrics$overall$avg_MCSE_bias_Sp),
              sprintf("%.4f", metrics$overall$max_MCSE_bias_Se),
              sprintf("%.4f", metrics$overall$max_MCSE_bias_Sp),
              sprintf("%.1f%%", metrics$overall$avg_coverage_prob_Se * 100),
              sprintf("%.1f%%", metrics$overall$avg_coverage_prob_Sp * 100),
              sprintf("%.1f", metrics$overall$avg_interval_width_Se),
              sprintf("%.1f", metrics$overall$avg_interval_width_Sp)
            )
          )
          
          # Threshold table stays similar but with RMSE
          threshold_df <- do.call(rbind, lapply(min_k:max_k, function(k) {
            data.frame(
              Threshold = k,
              RMSE_Se = sprintf("%.3f", metrics$by_threshold[[k]]$RMSE_Se),
              RMSE_Sp = sprintf("%.3f", metrics$by_threshold[[k]]$RMSE_Sp),
              MCSE_RMSE_Se = sprintf("%.4f", metrics$by_threshold[[k]]$MCSE_RMSE_Se),
              MCSE_RMSE_Sp = sprintf("%.4f", metrics$by_threshold[[k]]$MCSE_RMSE_Sp),
              Bias_Se = sprintf("%.3f", metrics$by_threshold[[k]]$mean_bias_Se),
              Bias_Sp = sprintf("%.3f", metrics$by_threshold[[k]]$mean_bias_Sp),
              abs_Bias_Se = sprintf("%.3f", metrics$by_threshold[[k]]$mean_abs_bias_Se),
              abs_Bias_Sp = sprintf("%.3f", metrics$by_threshold[[k]]$mean_abs_bias_Sp),
              MCSE_Bias_Se = sprintf("%.4f", metrics$by_threshold[[k]]$MCSE_bias_Se),
              MCSE_Bias_Sp = sprintf("%.4f", metrics$by_threshold[[k]]$MCSE_bias_Sp),
              Coverage_Se = sprintf("%.1f%%", metrics$by_threshold[[k]]$coverage_prob_Se * 100),
              Coverage_Sp = sprintf("%.1f%%", metrics$by_threshold[[k]]$coverage_prob_Sp * 100),
              Width_Se = sprintf("%.1f", metrics$by_threshold[[k]]$mean_interval_width_Se),
              Width_Sp = sprintf("%.1f", metrics$by_threshold[[k]]$mean_interval_width_Sp)
            )
          }))
          
          return(list(
            by_threshold = threshold_df,
            overall = overall_df
          ))
  
}


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


calc_min_sim_size_RMSE <- function(metrics, 
                                   target_MCSE_RMSE,
                                   min_k,
                                   max_k) {
  
  min_N_sim <- c()
  index <- 1
  
  for (k in min_k:max_k) {
    # Get RMSE and its MCSE for this threshold
    RMSE_k <- metrics$by_threshold[[k]]$RMSE_Se  # or RMSE_Sp
    MCSE_RMSE_k <- metrics$by_threshold[[k]]$MCSE_RMSE_Se
    
    # From MCSE_RMSE_k at current n_sims, calculate SD of squared errors
    current_n_sims <- metrics$n_sims
    
    # MCSE(RMSE) ≈ SD(squared_errors) / (2 * sqrt(n) * RMSE)
    # So: SD(squared_errors) ≈ MCSE(RMSE) * 2 * sqrt(n) * RMSE
    SD_squared_errors <- MCSE_RMSE_k * 2 * sqrt(current_n_sims) * RMSE_k
    
    # Calculate n needed for target MCSE
    # target_MCSE = SD_squared_errors / (2 * sqrt(n_target) * RMSE)
    # Solving for n_target:
    n_target <- ceiling((SD_squared_errors / (2 * target_MCSE_RMSE * RMSE_k))^2)
    
    min_N_sim[index] <- n_target
    index <- index + 1
  }
  
  return(list(
    min_required = min(min_N_sim, na.rm = TRUE),
    avg_required = mean(min_N_sim, na.rm = TRUE),
    max_required = max(min_N_sim, na.rm = TRUE)
  ))
}




















