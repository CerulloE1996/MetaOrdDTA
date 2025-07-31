



# 
# ##
# ## ---- sROC Plot: True vs Model-Predicted --------------------------------------------
# ##
# 
# create_sROC_plot <- function(data, 
#                              DGM_name, 
#                              base_size = 16) {
#   
#       require(dplyr)
#       require(ggplot2)
#       require(purrr)
#       require(tidyr)
#       require(patchwork)
#   
#       # Filter for specific DGM
#       plot_data <- data %>%
#         filter(DGM == DGM_name)
#       
#       plot_sroc <- plot_data %>%
#         # Calculate FPR (1 - Specificity) for both true and model
#         dplyr::mutate(
#           true_FPR = (100 - true_Sp)/100,  # Convert to proportion
#           model_FPR = (100 - Sp_median)/100,
#           true_Se_prop = true_Se/100,
#           model_Se_prop = Se_median/100
#         ) %>%
#         # Arrange by threshold to ensure proper line connections
#         dplyr::arrange(DGM, n_studies, Model, threshold) %>%
#         ggplot() +
#         ##
#         ## True sROC curve (black) - plot first
#         ##
#         geom_path(data = . %>% select(DGM, n_studies, threshold, true_FPR, true_Se_prop) %>% distinct(),
#                   aes(x = true_FPR, y = true_Se_prop),
#                   color = "black", linewidth = 1.5) +
#         ##
#         geom_point(data = . %>% select(DGM, n_studies, threshold, true_FPR, true_Se_prop) %>% distinct(),
#                    aes(x = true_FPR, y = true_Se_prop),
#                    color = "black", size = 3) +
#         ##
#         ## Model-predicted sROC curves (colored by model)
#         ##
#         geom_path(aes(x = model_FPR, y = model_Se_prop, color = Model, group = Model),
#                   linewidth = 1.0) +
#         ##
#         geom_point(aes(x = model_FPR, y = model_Se_prop, color = Model),
#                    size = 2, alpha = 0.5) +
#         ##
#         ## Diagonal reference line
#         geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
#         ##
#         scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                            limits = c(0, 1)) +
#         ##
#         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                            limits = c(0, 1)) +
#         ##
#         scale_color_brewer(palette = "Set1") +
#         ##
#         labs(x = "False Positive Rate (1 - Specificity)",
#              y = "Sensitivity",
#              title = "Summary ROC Curves: True vs Model Predictions",
#              subtitle = "Black = True sROC | Colored = Model estimates") +
#         ##
#         coord_equal() +
#         theme_bw() +
#         # theme(legend.position = "bottom",
#         #       strip.text = element_text(size = 8)) + 
#         theme(legend.position = "bottom",
#               # axis.text.x = element_text(angle = 45, hjust = 1),
#               ##
#               # Axis text
#               axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.9),
#               axis.text.y = element_text(size = base_size * 0.9),
#               # Axis titles
#               axis.title.x = element_text(size = base_size * 1.1),
#               axis.title.y = element_text(size = base_size * 1.1),
#               # Plot title and subtitle
#               plot.title = element_text(size = base_size * 1.3, face = "bold"),
#               plot.subtitle = element_text(size = base_size * 1.0),
#               # Facet labels
#               strip.text = element_text(size = base_size * 1.0, face = "bold"),
#               # If you have axis tick labels
#               axis.ticks = element_line(size = 0.5),
#               # LEGEND ELEMENTS
#               legend.title = element_text(size = base_size * 1.1, face = "bold"),  # Legend title
#               legend.text = element_text(size = base_size * 1.0),                  # Legend labels
#               legend.key.size = unit(1.2, "lines"),                                # Size of legend keys
#               legend.key.width = unit(1.5, "lines"),                               # Width of legend keys
#               legend.key.height = unit(1.2, "lines")                               # Height of legend keys
#         ) + 
#         # Facet
#         facet_wrap(~ DGM + n_studies, 
#                    nrow = 1)
#       
#       return(plot_sroc)
# 
# }


##
## ---- Function to create RMSE plots with MCSE error bars ---------------------------------------------------
##
create_RMSE_plot <- function(data, 
                             DGM_name, 
                             metric_type = c("Se", "Sp"),
                             base_size = 16) {
  
        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        
        metric_type <- match.arg(metric_type)
        
        # Filter for specific DGM
        plot_data <- data %>%
          filter(DGM == DGM_name)
        
        # Set up column names based on metric type
        if (metric_type == "Se") {
          rmse_col <- "avg_RMSE_Se"
          mcse_col <- "avg_MCSE_RMSE_Se"
          title_suffix <- "Sensitivity (Se)"
        } else {
          rmse_col <- "avg_RMSE_Sp"
          mcse_col <- "avg_MCSE_RMSE_Sp"
          title_suffix <- "Specificity (Sp)"
        }
        
        # Calculate y-axis limits
        plot_data_long <- plot_data %>%
          pivot_longer(cols = all_of(rmse_col),
                       names_to = "metric",
                       values_to = "RMSE") %>%
          mutate(
            MCSE = .data[[mcse_col]],
            metric = metric_type,
            lower_bound = RMSE - 1.96*MCSE,
            upper_bound = RMSE + 1.96*MCSE
          )
        
        # Calculate min and max for y-axis
        y_min <- min(plot_data_long$lower_bound, na.rm = TRUE) * 0.9  # 10% buffer below lowest point
        y_max <- max(plot_data_long$upper_bound, na.rm = TRUE) * 1.1  # 10% buffer above highest point
        
        # Create plot
        plot_data_long %>%
          ggplot(aes(x = Model, y = RMSE, colour = Model)) +
          geom_point(size = 5) +  # Use points instead of bar
          geom_col(alpha = 0.7) +
          geom_errorbar(aes(ymin = lower_bound, 
                            ymax = upper_bound),
                        width = 0.4) +
          scale_fill_brewer(palette = "Set1") +
          scale_y_continuous(limits = c(y_min, y_max),
                             expand = c(0, 0)) +  # Remove default expansion
          ##
          labs(y = "Root Mean Square Error",
               title = paste0("RMSE Comparison: ", DGM_name, " - ", title_suffix),
               subtitle = "Error bars = ± 1.96×MCSE") +
          theme_bw() +
          theme(legend.position = "none",
                # axis.text.x = element_text(angle = 45, hjust = 1),
                ##
                # Axis text
                axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.9),
                axis.text.y = element_text(size = base_size * 0.9),
                # Axis titles
                axis.title.x = element_text(size = base_size * 1.1),
                axis.title.y = element_text(size = base_size * 1.1),
                # Plot title and subtitle
                plot.title = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                # Facet labels
                strip.text = element_text(size = base_size * 1.0, face = "bold"),
                # If you have axis tick labels
                axis.ticks = element_line(size = 0.5)
          ) + 
          ##
          facet_wrap(~ n_studies, 
                     scales = "free_y", 
                     nrow = 1,
                     labeller = labeller(n_studies = function(x) paste(x, "studies")))
  
}
##
## ---- Function to create Bias plots for a specific DGM ---------------------------------------------------
##
create_bias_plot <- function(data, 
                             DGM_name, 
                             metric_type = c("Se", "Sp"),
                             base_size = 16) {
  
        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        
        metric_type <- match.arg(metric_type)
        
        # Filter for specific DGM
        plot_data <- data %>%
          filter(DGM == DGM_name)
        
        # Set up column names based on metric type
        if (metric_type == "Se") {
          bias_col <- "avg_abs_bias_Se"
          mcse_col <- "avg_MCSE_bias_Se"
          title_suffix <- "Sensitivity (Se)"
        } else {
          bias_col <- "avg_abs_bias_Sp"
          mcse_col <- "avg_MCSE_bias_Sp"
          title_suffix <- "Specificity (Sp)"
        }
        
        # Calculate y-axis limits
        plot_data_long <- plot_data %>%
          pivot_longer(cols = all_of(bias_col),
                       names_to = "metric",
                       values_to = "bias") %>%
          mutate(
            MCSE = .data[[mcse_col]],
            metric = metric_type,
            lower_bound = bias - 1.96*MCSE,
            upper_bound = bias + 1.96*MCSE
          )
        
        # Calculate min and max for y-axis
        # For bias, we might want to include 0 as the minimum if all values are positive
        y_min <- min(0, min(plot_data_long$lower_bound, na.rm = TRUE) * 0.9)
        y_max <- max(plot_data_long$upper_bound, na.rm = TRUE) * 1.1
        
        # Create plot
        plot_data_long %>%
          ggplot(aes(x = Model, y = bias, colour = Model)) +
          geom_point(size = 5) +  # Use points instead of bars
          geom_errorbar(aes(ymin = lower_bound, 
                            ymax = upper_bound),
                        width = 0.4) +
          scale_colour_brewer(palette = "Set1") +  # Changed to colour for points
          scale_y_continuous(limits = c(y_min, y_max),
                             expand = c(0, 0)) +  # Remove default expansion
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +  # Add reference line at 0
          ##
          labs(y = "Absolute Bias (%)",
               title = paste0("Bias Comparison: ", DGM_name, " - ", title_suffix),
               subtitle = "Error bars = ± 1.96×MCSE") +
          theme_bw(base_size = base_size) +
          theme(legend.position = "none",
                # Axis text
                axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.9),
                axis.text.y = element_text(size = base_size * 0.9),
                # Axis titles
                axis.title.x = element_text(size = base_size * 1.1),
                axis.title.y = element_text(size = base_size * 1.1),
                # Plot title and subtitle
                plot.title = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                # Facet labels
                strip.text = element_text(size = base_size * 1.0, face = "bold"),
                # If you have axis tick labels
                axis.ticks = element_line(size = 0.5)
          ) + 
          ##
          facet_wrap(~ n_studies, 
                     scales = "free_y", 
                     nrow = 1,
                     labeller = labeller(n_studies = function(x) paste(x, "studies")))
  
}
##
## ---- Function to create Coverage plots for a specific DGM ---------------------------------------------------
##
create_coverage_plot <- function(data, 
                                 DGM_name, 
                                 metric_type = c("Se", "Sp"),
                                 base_size = 16) {
  
      require(dplyr)
      require(ggplot2)
      require(purrr)
      require(tidyr)
      require(patchwork)
  
      metric_type <- match.arg(metric_type)
      
      # Filter for specific DGM
      plot_data <- data %>%
        filter(DGM == DGM_name)
      
      # Set up column names based on metric type
      if (metric_type == "Se") {
        coverage_col <- "avg_coverage_Se"
        title_suffix <- "Sensitivity (Se)"
      } else {
        coverage_col <- "avg_coverage_Sp"
        title_suffix <- "Specificity (Sp)"
      }
      
      # Create plot
      plot_data %>%
        pivot_longer(cols = all_of(coverage_col),
                     names_to = "metric",
                     values_to = "coverage") %>%
        mutate(metric = metric_type) %>%
        ggplot(aes(x = Model, y = coverage, color = Model)) +
        geom_point(size = 4) +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0.2, 1)) +  # Adjusted to capture O-HSROC's poor coverage
        scale_color_brewer(palette = "Set1") +
        labs(y = "Coverage Probability",
             title = paste0("95% CI Coverage: ", DGM_name, " - ", title_suffix),
             subtitle = "Red line = Nominal 95% coverage") +
        theme_bw() +
        theme(legend.position = "none",
              # axis.text.x = element_text(angle = 45, hjust = 1),
              ##
              # Axis text
              axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.9),
              axis.text.y = element_text(size = base_size * 0.9),
              # Axis titles
              axis.title.x = element_text(size = base_size * 1.1),
              axis.title.y = element_text(size = base_size * 1.1),
              # Plot title and subtitle
              plot.title = element_text(size = base_size * 1.3, face = "bold"),
              plot.subtitle = element_text(size = base_size * 1.0),
              # Facet labels
              strip.text = element_text(size = base_size * 1.0, face = "bold"),
              # If you have axis tick labels
              axis.ticks = element_line(size = 0.5)
        ) +  
        facet_wrap(~ n_studies, nrow = 1,
                   labeller = labeller(n_studies = function(x) paste(x, "studies")))
      
}




count_non_missing_by_column <- function(data, 
                                        missing_value = -1, 
                                        skip_first_col = TRUE,
                                        return_df = FALSE) {
        
        # Handle column selection
        if (skip_first_col) {
          if (ncol(data) < 2) {
            stop("Data must have at least 2 columns if skipping first column")
          }
          data_subset <- data[, -1, drop = FALSE]
          col_offset <- 1
        } else {
          data_subset <- data
          col_offset <- 0
        }
        
        # Count non-missing values
        counts <- colSums(data_subset != missing_value, na.rm = TRUE)
        
        # Return based on preference
        if (return_df) {
          result <- data.frame(
            original_col = (1:length(counts)) + col_offset,
            threshold = 1:length(counts),
            n_non_missing = as.numeric(counts),
            n_missing = nrow(data) - as.numeric(counts),
            pct_non_missing = round(100 * as.numeric(counts) / nrow(data), 1)
          )
          return(result)
        } else {
          return(as.numeric(counts))
        }
  
}
  
  
  
  
# Convert metrics list to tibble format
convert_metrics_to_tibble <- function(complete_results_list, 
                                      output_dir,
                                      test,
                                      true_DGM_Se,
                                      true_DGM_Sp,
                                      model_name,
                                      n_studies, 
                                      N,
                                      min_k,
                                      max_k) {

        metrics_list <- compute_simulation_metrics(complete_results_list = complete_results_list,
                                                   output_dir = output_dir,
                                                   test = test,
                                                   true_DGM_Se = true_DGM_Se, 
                                                   true_DGM_Sp = true_DGM_Sp,
                                                   min_k = min_k,
                                                   max_k = max_k)
        
        {
          MCSE_outs_1 <- calc_min_sim_size_RMSE(
            metrics = metrics_list,
            target_MCSE_RMSE = 0.125, 
            min_k = min_k,
            max_k = max_k
          )
          n_sims_0.125 <- MCSE_outs_1$avg_required
        }
        
        {
          MCSE_outs_2 <-  calc_min_sim_size_RMSE(
            metrics = metrics_list,
            target_MCSE_RMSE = 0.10, 
            min_k = min_k,
            max_k = max_k
          )
          n_sims_0.10 <- MCSE_outs_2$avg_required
        }
        
        tibble(
          DGM = test,
          Model = model_name,
          n_studies = n_studies,
          N = N,
          n_sims = metrics_list$n_sims,
          n_sims_0.125 = n_sims_0.125,
          n_sims_0.10 = n_sims_0.10,
          # Overall metrics:
          avg_RMSE_Se = metrics_list$overall$avg_RMSE_Se,
          avg_RMSE_Sp = metrics_list$overall$avg_RMSE_Sp,
          max_RMSE_Se = metrics_list$overall$max_RMSE_Se,
          max_RMSE_Sp = metrics_list$overall$max_RMSE_Sp,
          ##
          avg_MCSE_RMSE_Se = metrics_list$overall$avg_MCSE_RMSE_Se,
          avg_MCSE_RMSE_Sp = metrics_list$overall$avg_MCSE_RMSE_Sp,
          max_MCSE_RMSE_Se = metrics_list$overall$max_MCSE_RMSE_Se,
          max_MCSE_RMSE_Sp = metrics_list$overall$max_MCSE_RMSE_Sp,
          ##
          avg_abs_bias_Se = metrics_list$overall$avg_abs_bias_Se,
          avg_abs_bias_Sp = metrics_list$overall$avg_abs_bias_Sp,
          avg_MCSE_bias_Se = metrics_list$overall$avg_MCSE_bias_Se,
          avg_MCSE_bias_Sp = metrics_list$overall$avg_MCSE_bias_Sp,
          ##
          avg_coverage_Se = metrics_list$overall$avg_coverage_prob_Se,
          avg_coverage_Sp = metrics_list$overall$avg_coverage_prob_Sp,
          ##
          avg_width_Se = metrics_list$overall$avg_interval_width_Se,
          avg_width_Sp = metrics_list$overall$avg_interval_width_Sp
        )
          
}


get_threshold_data <- function(complete_results_list, 
                               output_dir,
                               test,
                               model_name, 
                               n_studies,
                               N,
                               min_k,
                               max_k,
                               true_DGM_Se, 
                               true_DGM_Sp) {
          
          metrics_list <- compute_simulation_metrics(complete_results_list = complete_results_list,
                                                     output_dir = output_dir,
                                                     test = test,
                                                     true_DGM_Se = true_DGM_Se,
                                                     true_DGM_Sp = true_DGM_Sp,
                                                     min_k = min_k,
                                                     max_k = max_k)
          
          n_thr <- length(metrics_list$by_threshold)
          
          # metrics_list$by_threshold
          # 
          # complete_results_list[[1]][[1]]
          
          # Extract metrics for each threshold
          threshold_data <- purrr::map_df(1:n_thr, function(k) {
                tibble(
                  DGM = test,
                  Model = model_name,
                  n_studies = n_studies,
                  N = N,
                  threshold = k,
                  # True values
                  true_Se = true_DGM_Se[k],
                  true_Sp = true_DGM_Sp[k],
                  # RMSE
                  RMSE_Se = metrics_list$by_threshold[[k]]$RMSE_Se,
                  RMSE_Sp = metrics_list$by_threshold[[k]]$RMSE_Sp,
                  # MCSE for RMSE
                  MCSE_RMSE_Se = metrics_list$by_threshold[[k]]$MCSE_RMSE_Se,
                  MCSE_RMSE_Sp = metrics_list$by_threshold[[k]]$MCSE_RMSE_Sp,
                  # Bias
                  Se_bias = metrics_list$by_threshold[[k]]$mean_bias_Se,
                  Sp_bias = metrics_list$by_threshold[[k]]$mean_bias_Sp,
                  # Coverage
                  coverage_Se = metrics_list$by_threshold[[k]]$coverage_prob_Se,
                  coverage_Sp = metrics_list$by_threshold[[k]]$coverage_prob_Sp,
                  # CI width
                  width_Se = metrics_list$by_threshold[[k]]$mean_interval_width_Se,
                  width_Sp = metrics_list$by_threshold[[k]]$mean_interval_width_Sp,
                  # Standard deviations
                  SD_Se = metrics_list$by_threshold[[k]]$SD_Se,
                  SD_Sp = metrics_list$by_threshold[[k]]$SD_Sp,
                  ##
                  Se_median = metrics_list$by_threshold[[k]]$Se_median,
                  Sp_median = metrics_list$by_threshold[[k]]$Sp_median,
                  ##
                  Se_mean = metrics_list$by_threshold[[k]]$Se_mean,
                  Sp_mean = metrics_list$by_threshold[[k]]$Sp_mean,
                  ##
                  Se_lower = metrics_list$by_threshold[[k]]$Se_lower,
                  Sp_lower = metrics_list$by_threshold[[k]]$Sp_lower,
                  ##
                  Se_upper = metrics_list$by_threshold[[k]]$Se_upper,
                  Sp_upper = metrics_list$by_threshold[[k]]$Sp_upper,
                  ##
                  Se_pred_lower = metrics_list$by_threshold[[k]]$Se_pred_lower,
                  Sp_pred_lower = metrics_list$by_threshold[[k]]$Sp_pred_lower,
                  ##
                  Se_pred_upper = metrics_list$by_threshold[[k]]$Se_pred_upper,
                  Sp_pred_upper = metrics_list$by_threshold[[k]]$Sp_pred_upper
                )
          })
          
          return(threshold_data)
  
}



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





compute_simulation_metrics <- function(test,
                                       complete_results_list, 
                                       output_dir,
                                       true_DGM_Se, 
                                       true_DGM_Sp,
                                       min_k,
                                       max_k) {
  
        n_studies <- complete_results_list[[1]]$n_studies
        
        # 
        file_name_string <- paste0( "sim_data_list_",
                                    "N_", 500, "_",
                                    "n_studies_", n_studies)
        ##
        data_list_file <- file.path(output_dir,
                                       paste0(file_name_string, ".RDS"))
        ##
        datasets_sim_list <- readRDS(data_list_file)
        
  
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
        
        
        non_missing_index_list <- missing_index_list <- list()
        
        for (i in 1:n_sims) {
          
            if (test == "GAD-2") { 
              
                 x <- datasets_sim_list[[i]]$x_GAD2_w_missing_thr
              
            } else if (test == "HADS") { 
              
                 x <- datasets_sim_list[[i]]$x_HADS_w_missing_thr
                 
            } else if (test == "BAI") { 
              
                 x <- datasets_sim_list[[i]]$x_BAI_w_missing_thr
              
            } else if (test == "GAD-7") { 
              
                 x <- datasets_sim_list[[i]]$x_GAD7_w_missing_thr
                 
            }
          
            non_missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) >= 2)
            missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) < 2)
          
        }
          
          
        
        # Initialize storage matrices
        coverage_Se <- coverage_Sp <- matrix(NA, n_sims, max_k)
        interval_width_Se <- interval_width_Sp <- matrix(NA, n_sims, max_k)
        bias_Se <- bias_Sp <- matrix(NA, n_sims, max_k)
        squared_error_Se <- squared_error_Sp <- matrix(NA, n_sims, max_k)
        ##
        Se_median_vec <- Sp_median_vec <- c()
        Se_mean_vec   <- Sp_mean_vec <- c()
        Se_lower_vec  <- Sp_lower_vec <- c()
        Se_upper_vec  <- Sp_upper_vec <- c()
        Se_pred_lower_vec <- Sp_pred_lower_vec <- c()
        Se_pred_upper_vec <- Sp_pred_upper_vec <- c()
        
        # Store point estimates for SD calculation
        point_est_Se <- point_est_Sp <- matrix(NA, n_sims, max_k)
        
        # Loop through simulations (only complete ones now)
        for (i in 1:n_sims) {
          
                result <- complete_results[[i]]
                
                ## replace values w/ less than 2 studies with NA's (otherwise can't compare to stratified bivariate fairly):
                result$Se_median_vec[missing_index_list[[i]]] <- NA
                result$Sp_median_vec[missing_index_list[[i]]] <- NA
                ##
                result$Se_mean_vec[missing_index_list[[i]]] <- NA
                result$Sp_mean_vec[missing_index_list[[i]]] <- NA
                ##
                result$Se_lower_vec[missing_index_list[[i]]] <- NA
                result$Sp_lower_vec[missing_index_list[[i]]] <- NA
                ##
                result$Se_upper_vec[missing_index_list[[i]]] <- NA
                result$Sp_upper_vec[missing_index_list[[i]]] <- NA
                ##
                result$Se_pred_lower_vec[missing_index_list[[i]]] <- NA
                result$Sp_pred_lower_vec[missing_index_list[[i]]] <- NA
                ##
                result$Se_pred_upper_vec[missing_index_list[[i]]] <- NA
                result$Sp_pred_upper_vec[missing_index_list[[i]]] <- NA
                
                for (k in min_k:max_k) {
                  
                        Se_median_vec[k] <- result$Se_median_vec[k]
                        Sp_median_vec[k] <- result$Sp_median_vec[k]
                        ##
                        Se_mean_vec[k] <- result$Se_mean_vec[k]
                        Sp_mean_vec[k] <- result$Sp_mean_vec[k]
                        ##
                        Se_lower_vec[k] <- result$Se_lower_vec[k]
                        Sp_lower_vec[k] <- result$Sp_lower_vec[k]
                        ##
                        Se_upper_vec[k] <- result$Se_upper_vec[k]
                        Sp_upper_vec[k] <- result$Sp_upper_vec[k]
                        ##
                        Se_pred_lower_vec[k] <- result$Se_pred_lower_vec[k]
                        Sp_pred_lower_vec[k] <- result$Sp_pred_lower_vec[k]
                        ##
                        Se_pred_upper_vec[k] <- result$Se_pred_upper_vec[k]
                        Sp_pred_upper_vec[k] <- result$Sp_pred_upper_vec[k]
                        
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
            ##
            # Mean interval width
            mean_interval_width_Se = mean(interval_width_Se[, k], na.rm = TRUE),
            mean_interval_width_Sp = mean(interval_width_Sp[, k], na.rm = TRUE),
            ##
            # RMSE
            RMSE_Se = RMSE_Se_k,
            RMSE_Sp = RMSE_Sp_k,
            ##
            # MCSE of RMSE
            MCSE_RMSE_Se = MCSE_RMSE_Se_k,
            MCSE_RMSE_Sp = MCSE_RMSE_Sp_k,
            ##
            # Mean bias
            mean_bias_Se = mean_bias_Se_k,
            mean_bias_Sp = mean_bias_Sp_k,
            ##
            # Mean absolute bias
            mean_abs_bias_Se = abs(mean_bias_Se_k),
            mean_abs_bias_Sp = abs(mean_bias_Sp_k),
            ##
            # MCSE of bias
            MCSE_bias_Se = MCSE_bias_Se_vec[k],
            MCSE_bias_Sp = MCSE_bias_Sp_vec[k],
            ##
            # SD
            SD_Se = SD_Se_vec[k],
            SD_Sp = SD_Sp_vec[k],
            ##
            ##
            ##
            Se_median = Se_median_vec[k],
            Sp_median = Sp_median_vec[k],
            ##
            Se_mean = Se_mean_vec[k],
            Sp_mean = Sp_mean_vec[k],
            ##
            Se_lower = Se_lower_vec[k],
            Sp_lower = Sp_lower_vec[k],
            ##
            Se_upper = Se_upper_vec[k],
            Sp_upper = Sp_upper_vec[k],
            ##
            Se_pred_lower = Se_pred_lower_vec[k],
            Sp_pred_lower = Sp_pred_lower_vec[k],
            ##
            Se_pred_upper = Se_pred_upper_vec[k],
            Sp_pred_upper = Sp_pred_upper_vec[k]
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
                       # "Max |Bias| Se", "Max |Bias| Sp",
                       "Mean MCSE(Bias) Se", "Mean MCSE(Bias) Sp",
                       # "Max MCSE(Bias) Se", "Max MCSE(Bias) Sp",
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
              # sprintf("%.3f", metrics$overall$max_abs_bias_Se),
              # sprintf("%.3f", metrics$overall$max_abs_bias_Sp),
              sprintf("%.4f", metrics$overall$avg_MCSE_bias_Se),
              sprintf("%.4f", metrics$overall$avg_MCSE_bias_Sp),
              # sprintf("%.4f", metrics$overall$max_MCSE_bias_Se),
              # sprintf("%.4f", metrics$overall$max_MCSE_bias_Sp),
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




















