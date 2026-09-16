


#' @export
summarise_beta_mu <- function(beta_mu, 
                              test_names = c("GAD-2", "GAD-7", "HADS", "BAI", 
                                             "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                              covariate_names,  # e.g., c("intercept", "SCID", "Structured")
                              magnitude_threshold = 0.25,
                              skip_intercept = TRUE) {
  
          # Parse beta_mu indices
          parsed <- beta_mu %>%
            mutate(
              # Extract indices from variable name like "beta_mu[1,2,3]"
              indices = str_extract(variable, "\\[.*\\]"),
              test_idx = as.integer(str_extract(str_extract(indices, "\\[\\d+"), "\\d+")),
              ds_idx = as.integer(str_extract(str_extract(indices, ",\\d+,"), "\\d+")),
              cov_idx = as.integer(str_extract(str_extract(indices, ",\\d+\\]"), "\\d+"))
            ) %>%
            mutate(
              test = test_names[test_idx],
              disease_status = ifelse(ds_idx == 1, "ND (Sp)", "D (Se)"),
              covariate = covariate_names[cov_idx],
              # Check significance (95% CI excludes 0)
              significant = (`2.5%` > 0 | `97.5%` < 0),
              # Check large magnitude
              large_magnitude = abs(mean) >= magnitude_threshold,
              # Direction
              direction = case_when(
                mean > 0 ~ "+",
                mean < 0 ~ "-",
                TRUE ~ "0"
              )
            )
          
          # Skip intercept if requested
          if (skip_intercept) {
            parsed <- parsed %>% filter(cov_idx > 1)
          }
          
          # Filter to significant OR large magnitude
          notable <- parsed %>%
            filter(significant | large_magnitude) %>%
            arrange(covariate, disease_status, test) %>%
            select(test, disease_status, covariate, mean, `2.5%`, `97.5%`, 
                   significant, large_magnitude, direction)
          
          # Print summary
          cat("\n========== STATISTICALLY SIGNIFICANT EFFECTS (95% CI excludes 0) ==========\n\n")
          sig_only <- parsed %>% filter(significant)
          if (nrow(sig_only) > 0) {
            sig_only %>%
              arrange(covariate, disease_status, test) %>%
              select(test, disease_status, covariate, mean, `2.5%`, `97.5%`, direction) %>%
              print(n = 100)
          } else {
            cat("None\n")
          }
          
          cat("\n\n========== LARGE MAGNITUDE ('PRACTICALLY SIGNIFICANT'), BUT NOT STATISTICALLY SIGNIFICANT (|β| ≥",
              magnitude_threshold, ") ==========\n\n")
          large_only <- parsed %>% filter(large_magnitude & !significant)
          if (nrow(large_only) > 0) {
            large_only %>%
              arrange(covariate, disease_status, test) %>%
              select(test, disease_status, covariate, mean, `2.5%`, `97.5%`, direction) %>%
              print(n = 100)
          } else {
            cat("None\n")
          }
          
          # Return full parsed tibble invisibly
          invisible(parsed)
  
}



#' quick_kappa_prior
#' @export
quick_kappa_prior <- function(log_mean = 0, 
                              log_sd = 1, 
                              df = 5) {
  
        # Calculate key statistics on log scale
        log_quants <- c(0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995)
        log_q_vals <- log_mean + log_sd * qt(log_quants, df = df)
        
        # Transform to kappa scale
        kappa_q_vals <- exp(log_q_vals)
        
        # Simulate for mean and SD on kappa scale
        set.seed(123)
        log_samples <- log_mean + log_sd * rt(10000, df = df)
        kappa_samples <- exp(log_samples)
        kappa_samples <- kappa_samples[kappa_samples < quantile(kappa_samples, 0.9999)]
        
        # Print summary
        cat(sprintf("\nPrior: log(κ) ~ Student-t(df=%g, location=%.3f, scale=%.3f)\n", 
                    df, log_mean, log_sd))
        # cat("=" * 55, "\n")
        cat("\nKappa scale statistics:\n")
        cat(sprintf("Mean:     %.3f\n", mean(kappa_samples)))
        cat(sprintf("SD:       %.3f\n", sd(kappa_samples)))
        cat(sprintf("Median:   %.3f\n", kappa_q_vals[5]))
        ##
        # cat(sprintf("99%% CI:   [%.3f, %.3f]\n", kappa_q_vals[1], kappa_q_vals[9]))
        cat(sprintf("95%% CI:   [%.3f, %.3f]\n", kappa_q_vals[2], kappa_q_vals[8]))
        cat(sprintf("90%% CI:   [%.3f, %.3f]\n", kappa_q_vals[3], kappa_q_vals[7]))
        cat(sprintf("IQR:      [%.3f, %.3f]\n",  kappa_q_vals[4], kappa_q_vals[6]))
        
        # cat("\nProbabilities:\n")
        # cat(sprintf("P(κ < 0.5):  %.1f%%\n", mean(kappa_samples < 0.5) * 100))
        # cat(sprintf("P(κ > 1):    %.1f%%\n", mean(kappa_samples > 1) * 100))
        # cat(sprintf("P(κ > 2):    %.1f%%\n", mean(kappa_samples > 2) * 100))
        # cat(sprintf("P(κ > 3):    %.1f%%\n", mean(kappa_samples > 3) * 100))
        
        # # Quick plot on kappa scale only
        # log_x <- seq(log(0.01), log(10), length.out = 500)
        # kappa <- exp(log_x)
        # log_dens <- dt((log_x - log_mean)/log_sd, df = df)/log_sd
        # kappa_dens <- log_dens / kappa
        # 
        # plot(kappa, kappa_dens, type = "l", lwd = 3, col = "darkblue",
        #      main = sprintf("Kappa prior (df=%g)", df),
        #      xlab = "κ (kappa)", ylab = "Density",
        #      xlim = c(0, 5), ylim = c(0, max(kappa_dens[kappa < 5])))
        # 
        # # Add reference lines
        # abline(v = 1, lty = 2, col = "red", lwd = 2)
        # abline(v = kappa_q_vals[c(2, 5, 8)], lty = c(3, 1, 3), 
        #        col = c("gray40", "black", "gray40"), lwd = c(1, 2, 1))
        # 
        # # Add labels
        # text(1, max(kappa_dens[kappa < 5]) * 0.9, "κ=1", col = "red", pos = 4)
        # text(kappa_q_vals[4], max(kappa_dens[kappa < 5]) * 0.8, 
        #      sprintf("Median=%.2f", kappa_q_vals[5]), pos = 4)
        # 
        # # Add shaded 95% CI
        # kappa_ci <- kappa[kappa >= kappa_q_vals[2] & kappa <= kappa_q_vals[8]]
        # dens_ci <- kappa_dens[kappa >= kappa_q_vals[2] & kappa <= kappa_q_vals[8]]
        # polygon(c(kappa_ci, rev(kappa_ci)), c(dens_ci, rep(0, length(dens_ci))),
        #         col = rgb(0, 0, 1, 0.2), border = NA)
        
        invisible(kappa_q_vals)
}

      






 
#' compare_kappa_dofs
#' @export
compare_kappa_dofs <- function(log_mean = 0,
                               log_sd = 1, 
                               dofs = c(2, 3, 5, 10, 30)) {
        
            cat("\nComparing kappa priors with different DOF\n")
            cat(sprintf("Prior: log(κ) ~ t(df, location=%.2f, scale=%.2f)\n", log_mean, log_sd))
            # cat("=" * 65, "\n")
            
            results <- data.frame()
            
            for (df in dofs) {
              # Quantiles on log scale
              log_q <- log_mean + log_sd * qt(c(0.025, 0.5, 0.975), df = df)
              kappa_q <- exp(log_q)
              
              # Simulate for mean and SD
              set.seed(123)
              log_sim <- log_mean + log_sd * rt(10000, df = df)
              kappa_sim <- exp(log_sim)
              kappa_sim <- kappa_sim[kappa_sim < quantile(kappa_sim, 0.999)]
              
              row <- data.frame(
                DOF = df,
                Kappa_Mean = round(mean(kappa_sim), 3),
                Kappa_SD = round(sd(kappa_sim), 3),
                Kappa_Median = round(kappa_q[2], 3),
                Kappa_CI_Lower = round(kappa_q[1], 3),
                Kappa_CI_Upper = round(kappa_q[3], 3),
                CI_Width = round(kappa_q[3] - kappa_q[1], 3),
                P_gt_1 = round(mean(kappa_sim > 1) * 100, 1),
                P_gt_2 = round(mean(kappa_sim > 2) * 100, 1)
              )
              results <- rbind(results, row)
            }
            
            print(results)
            
            # Visual comparison
            par(mfrow = c(1, 2))
            
            # Density plot
            plot(NULL, xlim = c(0, 5), ylim = c(0, 2),
                 main = "Kappa prior densities", 
                 xlab = "κ (kappa)", ylab = "Density")
            
            colors <- c("red", "blue", "darkgreen", "purple", "orange", "black")
            
            for (i in 1:length(dofs)) {
              df <- dofs[i]
              log_x <- seq(log(0.01), log(10), length.out = 500)
              kappa <- exp(log_x)
              log_dens <- dt((log_x - log_mean)/log_sd, df = df)/log_sd
              kappa_dens <- log_dens / kappa
              
              lines(kappa[kappa <= 5], kappa_dens[kappa <= 5], 
                    col = colors[i], lwd = 2)
            }
            
            abline(v = 1, lty = 2, col = "gray40")
            legend("topright", legend = paste("df =", dofs), 
                   col = colors[1:length(dofs)], lwd = 2, cex = 0.8)
            
            # CDF plot  
            plot(NULL, xlim = c(0, 5), ylim = c(0, 1),
                 main = "Kappa prior CDFs", 
                 xlab = "κ (kappa)", ylab = "P(κ ≤ x)")
            
            for (i in 1:length(dofs)) {
              df <- dofs[i]
              kappa_seq <- seq(0.01, 5, length.out = 200)
              log_kappa_seq <- log(kappa_seq)
              cdf_vals <- pt((log_kappa_seq - log_mean)/log_sd, df = df)
              
              lines(kappa_seq, cdf_vals, col = colors[i], lwd = 2)
            }
            
            abline(v = 1, lty = 2, col = "gray40")
            abline(h = 0.5, lty = 3, col = "gray60")
            legend("bottomright", legend = paste("df =", dofs), 
                   col = colors[1:length(dofs)], lwd = 2, cex = 0.8)
            
            par(mfrow = c(1, 1))
            
            cat("\nInterpretation:\n")
            cat("- Lower DOF → heavier tails → more mass on extreme κ values\n")
            cat("- P(κ>1): Probability of better than random discrimination\n")
            cat("- P(κ>2): Probability of near-perfect discrimination\n")
            
            invisible(results)
  
}

# Example:
# compare_dofs(mean = 0, sd = 1)

# 
# mean_nd <- log(200)
# mean_d <- log(50)
# ##
# quick_kappa_prior(log_mean = mean_nd, 
#                   log_sd = 1, 
#                   df = 5)
# 
# 
# quick_kappa_prior(log_mean = mean_d, 
#                   log_sd = 1, 
#                   df = 5)
# 
# ##
# mean_nd <- log(100)
# mean_d <- log(25)
# 
# quick_kappa_prior(log_mean = mean_nd, 
#                       log_sd = 1, 
#                       df = 5)
# 












## method is either "sigma" or kappa" or "alpha":
#' induced_Dirichlet_ppc_plot
#' @export
induced_Dirichlet_ppc_plot <- function(  
                                         N = 5000,
                                         n_cat,
                                         other_args_list = NULL,
                                         log_mean,
                                         log_sd,
                                         df
                                         
) {
  
    require(ggplot2)
  
    p_if_uniform <- 1/n_cat
    
    message(paste("p_if_uniform = ", p_if_uniform))
    
    kappa <- NULL
    res   <- array(NA, dim = c(n_cat, N))
    alpha <- array(NA, dim = c(n_cat, N))
    log_alpha <- array(NA, dim = c(n_cat, N))
  
    # Simulate for mean and SD on kappa scale
    set.seed(123)
    ##
    log_kappa <- log_mean + log_sd * rt(10000, df = df)
    kappa <- exp(log_kappa)
    kappa <- kappa[kappa < quantile(kappa, 0.9999)]
    
    # use_log_kappa <- other_args_list$use_log_kappa
    # prior_dirichlet_cat_means_alpha <- other_args_list$prior_dirichlet_cat_means_alpha
    
  
    phi <- gtools::rdirichlet(N, rep(1, n_cat)) 
    
    log_phi <- log(phi)
    
    for (i in 1:N) {
          
          log_alpha[,i] <- log_kappa[i] + log_phi[i, 1:n_cat]
          ## Compute alpha:
          alpha[,i] <- exp(log_alpha[,i])
          ##
          res[,i] <- gtools::rdirichlet(n = 1, alpha[,i])
      
    }
 
  
  df <- data.frame(p1 = res[1,], 
                   p2 = res[2,], dist = 1) %>% 
    dplyr::filter(!is.na(p1), !is.na(p2))     
  
  g1 <- ggplot(df, aes(p1, p2)) + 
    geom_hex() + 
    scale_fill_continuous(trans = "log10") +
    theme_bw() + 
    # xlab(TeX("$P_{i}$")) + 
    # ylab(TeX("$P_{j}$")) + 
    xlim(0, 3.0*p_if_uniform) + 
    ylim(0, 3.0*p_if_uniform) + 
    geom_hline(yintercept = p_if_uniform) + 
    geom_vline(xintercept = p_if_uniform)
  
  print(g1)
  
  return(list(res = res, 
              alpha = alpha,
              kappa = kappa,
              log_alpha = log_alpha))
  
}


# N = 5000
# n_cat = 63
# log_mean = log(200)
# log_sd = 1
# df = 5
# 
# 
# quick_kappa_prior(log_mean = log(200),
#                   log_sd = 1,
#                   df = 5)
# quick_kappa_prior(log_mean = log(200),
#                   log_sd = 1,
#                   df = 15)
# 
# 
# outs <- induced_Dirichlet_ppc_plot(N = 5000,
#                            n_cat = 63,
#                            log_mean = log(200),
#                            log_sd = 1,
#                            df = 15)
# # 
# 



assess_interval_data_recovery <- function(x_matrix, 
                                          study_names = NULL) {
  
        # x_matrix: matrix with studies as rows, thresholds as columns
        # values are cumulative counts, -1 indicates missing
        
        n_studies <- nrow(x_matrix)
        n_thresholds <- ncol(x_matrix) - 1  # First col is study ID
        
        if(is.null(study_names)) {
          study_names <- x_matrix[, 1]
        }
        
        # Remove study ID column if present
        if(ncol(x_matrix) > n_thresholds) {
          counts <- x_matrix[, -1]
        } else {
          counts <- x_matrix
        }
        
        # Initialize results
        results <- list()
        
        # For each study
        for(s in 1:n_studies) {
          study_counts <- counts[s, ]
          n_intervals <- length(study_counts) - 1
          
          # Track which intervals are recoverable
          intervals_recoverable <- rep(FALSE, n_intervals)
          interval_values <- rep(NA, n_intervals)
          
          for(i in 1:n_intervals) {
            # Can only compute interval if both consecutive cumulative counts exist
            if(study_counts[i] != -1 && study_counts[i+1] != -1) {
              intervals_recoverable[i] <- TRUE
              interval_values[i] <- study_counts[i] - study_counts[i+1]
            }
          }
          
          results[[s]] <- list(
            study_id = study_names[s],
            cumulative_missing = sum(study_counts == -1),
            cumulative_available = sum(study_counts != -1),
            intervals_recoverable = sum(intervals_recoverable),
            intervals_missing = sum(!intervals_recoverable),
            interval_recovery_rate = mean(intervals_recoverable),
            interval_values = interval_values,
            recoverable_positions = which(intervals_recoverable)
          )
        }
        
        # Overall summary
        total_cumulative_data_points <- nrow(counts) * ncol(counts) ##  n_studies * length(counts[1,])
        total_cumulative_missing <- sum(counts == -1)
        total_interval_data_points <- n_studies * n_intervals
        total_intervals_recoverable <- sum(sapply(results, function(x) x$intervals_recoverable))
        ##
        interval_missing_rate = 1 - (total_intervals_recoverable / total_interval_data_points)
        cumulative_missing_rate = total_cumulative_missing / total_cumulative_data_points
        
        summary_stats <- list(
          cumulative_missing_rate = cumulative_missing_rate,
          interval_missing_rate = interval_missing_rate,
          cumulative_missing_pct = 100 * total_cumulative_missing / total_cumulative_data_points,
          interval_missing_pct = 100 * (1 - total_intervals_recoverable / total_interval_data_points),
          data_loss_amplification = (1 - total_intervals_recoverable / total_interval_data_points) / 
            (total_cumulative_missing / total_cumulative_data_points),
          data_loss_difference = cumulative_missing_rate - interval_missing_rate
        )
        
        return(list(
          by_study = results,
          summary = summary_stats
        ))
  
}


# More detailed comparison function
compare_data_approaches <- function(diseased_data, 
                                    nondiseased_data) {
  
  dis_recovery <- assess_interval_data_recovery(diseased_data)
  nondis_recovery <- assess_interval_data_recovery(nondiseased_data)
  
  cat("=== DATA RECOVERY COMPARISON ===\n\n")
  cat("CUMULATIVE DATA (Binomial Factorization Approach):\n")
  cat(sprintf("  Diseased missing rate: %.1f%%\n", dis_recovery$summary$cumulative_missing_pct))
  cat(sprintf("  Non-diseased missing rate: %.1f%%\n", nondis_recovery$summary$cumulative_missing_pct))
  cat(sprintf("  Overall missing rate: %.1f%%\n", 
              mean(c(dis_recovery$summary$cumulative_missing_pct, 
                     nondis_recovery$summary$cumulative_missing_pct))))
  
  cat("\nINTERVAL DATA (LC-MVP Approach):\n")
  cat(sprintf("  Diseased missing rate: %.1f%%\n", dis_recovery$summary$interval_missing_pct))
  cat(sprintf("  Non-diseased missing rate: %.1f%%\n", nondis_recovery$summary$interval_missing_pct))
  cat(sprintf("  Overall missing rate: %.1f%%\n",
              mean(c(dis_recovery$summary$interval_missing_pct,
                     nondis_recovery$summary$interval_missing_pct))))
  
  cat("\nDATA LOSS AMPLIFICATION:\n")
  cat(sprintf("  Diseased: %.2fx more missing data\n", dis_recovery$summary$data_loss_amplification))
  cat(sprintf("  Non-diseased: %.2fx more missing data\n", nondis_recovery$summary$data_loss_amplification))
  
  # Study-level breakdown
  cat("\n=== STUDY-LEVEL BREAKDOWN ===\n")
  worst_studies <- which(sapply(dis_recovery$by_study, 
                                function(x) x$interval_recovery_rate) < 0.3)
  if(length(worst_studies) > 0) {
    cat("\nStudies with <30% interval recovery:\n")
    for(s in worst_studies) {
      cat(sprintf("  Study %d: Only %d/%d intervals recoverable (%.1f%%)\n",
                  dis_recovery$by_study[[s]]$study_id,
                  dis_recovery$by_study[[s]]$intervals_recoverable,
                  dis_recovery$by_study[[s]]$intervals_recoverable + 
                    dis_recovery$by_study[[s]]$intervals_missing,
                  100 * dis_recovery$by_study[[s]]$interval_recovery_rate))
    }
  }
  
  return(list(diseased = dis_recovery, 
              nondiseased = nondis_recovery))
  
}













