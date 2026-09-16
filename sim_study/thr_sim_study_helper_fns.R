


 


##
## ---- Function to create RMSE plots with MCSE error bars ---------------------------------------------------
##
create_RMSE_plot <- function(  data, 
                               DGM_name = NULL,  
                               metric_type = c("Se", "Sp"),
                               base_size = 16
) {

        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        ##
        metric_type <- match.arg(metric_type)
        ##
        ## ---- Filter for specific DGM only if provided:
        ##
        if (!is.null(DGM_name)) {
                plot_data <- data %>% filter(DGM == DGM_name)
        } else {
                plot_data <- data
        }
        ##
        ## ---- Set up column names based on metric type:
        ##
        if (metric_type == "Se") {
          
                rmse_col <- "avg_RMSE_Se"
                mcse_col <- "avg_MCSE_RMSE_Se"
                title_suffix <- "Sensitivity (Se)"
            
        } else {
          
                rmse_col <- "avg_RMSE_Sp"
                mcse_col <- "avg_MCSE_RMSE_Sp"
                title_suffix <- "Specificity (Sp)"
            
        }
        ##
        ## ---- Calculate y-axis limits:
        ##
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
        ##
        best_models <-  plot_data_long %>%
                        group_by(DGM, n_studies) %>%
                        summarise(
                          min_rmse = min(RMSE),  # First find minimum RMSE
                          best_upper = max(upper_bound[RMSE == min_rmse]),  # Then get its upper bound
                          .groups = "drop"
                        )
        ##
        ## ---- Calculate min and max for y-axis:
        ##
        y_min <- min(plot_data_long$lower_bound, na.rm = TRUE) * 0.90
        y_max <- max(plot_data_long$upper_bound, na.rm = TRUE) * 1.10
        ##
        ## ---- Create plot:
        ##
        p <- plot_data_long %>%
          left_join(best_models, by = c("DGM", "n_studies")) %>%  # Add best_upper to data
          ggplot(aes(x = Model, y = RMSE, colour = Model)) +
          geom_hline(aes(yintercept = best_upper), 
                     linetype = "dashed", 
                     color = "darkgreen", 
                     alpha = 0.6) +  # Reference line for "best" group
          geom_point(size = 5) +
          geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.4) +
          ##
          labs(y = "Root Mean Square Error",
               title = paste0("RMSE: ", title_suffix, " - Test: ", unique(plot_data$test)[1]),
               subtitle = "Error bars = ± 1.96×MCSE") +
          theme_bw() +
          theme(legend.position = "none",
                axis.text.x   = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
                axis.text.y   = element_text(size = base_size * 0.9),
                axis.title.x  = element_text(size = base_size * 1.1),
                axis.title.y  = element_text(size = base_size * 1.1),
                plot.title    = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                ##
                strip.text = element_text(size = base_size * 0.65, face = "bold"),  # <- Make this smaller
                strip.background = element_rect(fill = "grey90"),
                ##
                axis.ticks = element_line(size = 0.5)
          ) +
          ##
          scale_x_discrete(limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          scale_color_brewer(palette = "Set1",
                             limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          ##
          facet_grid(  DGM ~ n_studies,
                       scales = "free_y",
                       labeller = labeller(
                         n_studies = function(x) paste(x, "studies"),
                         DGM = function(x) {
                           case_when(
                             x == "Jones"        ~ "DGM: Jones",
                             x == "Xu_fixed_C"   ~ "DGM: Ord-Biv-FC",
                             x == "Xu_random_C"  ~ "DGM: Ord-Biv-RC", 
                             x == "R&G_random_C" ~ "DGM: Ord-HSROC-RC",
                             TRUE ~ x
                           )
                         }
                       ))
          ##
          return(p)
  
}
##
## ---- Function to create Bias plots (modified for multiple DGMs) ---------------------------------------------------
##
create_bias_plot <- function(  data, 
                               DGM_name = NULL,  # Make optional - if NULL, plot all DGMs
                               metric_type = c("Se", "Sp"),
                               base_size = 16
) {
  
        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        ##
        metric_type <- match.arg(metric_type)
        ##
        # Filter for specific DGM only if provided
        if (!is.null(DGM_name)) {
          plot_data <- data %>% filter(DGM == DGM_name)
        } else {
          plot_data <- data
        }
        
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
        ##
        ## ---- Calculate y-axis limits:
        ##
        plot_data_long <- plot_data %>%
          pivot_longer(cols = all_of(bias_col),
                       names_to = "metric",
                       values_to = "bias") %>%
          mutate(
            MCSE = .data[[mcse_col]],
            metric = metric_type,
            lower_bound = bias - 1.96*MCSE,
            upper_bound = bias + 1.96*MCSE)
        ##
        best_models <- plot_data_long %>%
          group_by(DGM, n_studies) %>%
          summarise(
            best_upper = min(upper_bound[bias == min(bias)]),  # Upper bound of best model
            .groups = "drop"
          )
        ##
        # Calculate min and max for y-axis
        # For bias, we might want to include 0 as the minimum if all values are positive
        ##
        y_min <- min(0, min(plot_data_long$lower_bound, na.rm = TRUE) * 0.9)
        y_max <- max(plot_data_long$upper_bound, na.rm = TRUE) * 1.1
        ##
        ## ---- Create plot:
        ##
        p <- plot_data_long %>%
          left_join(best_models, by = c("DGM", "n_studies")) %>%  # Add best_upper to data
          ggplot(aes(x = Model, y = bias, colour = Model)) +
          geom_hline(aes(yintercept = best_upper),
                     linetype = "dashed",
                     color = "darkgreen",
                     alpha = 0.6) +  # Reference line for "best" group
          geom_point(size = 5) +
          geom_errorbar(aes( ymin = lower_bound, 
                             ymax = upper_bound), 
                        width = 0.4) +
        # p <- plot_data_long %>%
        #   ggplot(aes(x = Model, y = bias, colour = Model)) +
        #   geom_point(size = 5) +  
        #   geom_errorbar(aes(ymin = lower_bound, 
        #                     ymax = upper_bound),
        #                 width = 0.4) +
          scale_colour_brewer(palette = "Set1") +
          geom_hline( yintercept = 0, 
                      linetype = "dashed", 
                      color = "gray50",
                      alpha = 0.5) +
          ##
          labs(y = "Absolute Bias (%)",
               title = paste0("Bias: ", title_suffix, " - Test: ", unique(plot_data$test)[1]),
               subtitle = "Error bars = ± 1.96×MCSE") +
          theme_bw(base_size = base_size) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
                axis.text.y = element_text(size = base_size * 0.9),
                axis.title.x = element_text(size = base_size * 1.1),
                axis.title.y = element_text(size = base_size * 1.1),
                plot.title = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                ##
                strip.text = element_text(size = base_size * 0.65, face = "bold"),
                strip.background = element_rect(fill = "grey90"),
                ##
                axis.ticks = element_line(size = 0.5)
          )  +
          scale_x_discrete(limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          scale_color_brewer(palette = "Set1",
                             limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          ##
          facet_grid(DGM ~ n_studies, 
                     scales = "free_y",
                     labeller = labeller(
                       n_studies = function(x) paste(x, "studies"),
                       DGM = function(x) {
                         case_when(
                           x == "Jones"        ~ "DGM: Jones",
                           x == "Xu_fixed_C"   ~ "DGM: Ord-Biv-FC",
                           x == "Xu_random_C"  ~ "DGM: Ord-Biv-RC", 
                           x == "R&G_random_C" ~ "DGM: Ord-HSROC-RC",
                           TRUE ~ x
                         )
                       }
                     ))
        
        return(p)
        
}
##
## ---- Function to create Coverage plots for a specific DGM ---------------------------------------------------
##
create_coverage_plot <- function(data, 
                                 DGM_name = NULL,  # Make optional - if NULL, plot all DGMs
                                 metric_type = c("Se", "Sp"),
                                 base_size = 16
) {
  
        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        ##
        metric_type <- match.arg(metric_type)
        
        ## Filter for specific DGM only if provided
        if (!is.null(DGM_name)) {
          plot_data <- data %>% filter(DGM == DGM_name)
        } else {
          plot_data <- data
        }
        ##
        ## Set up column names based on metric type
        ##
        if (metric_type == "Se") {
          coverage_col <- "avg_coverage_Se"
          title_suffix <- "Sensitivity (Se)"
        } else {
          coverage_col <- "avg_coverage_Sp"
          title_suffix <- "Specificity (Sp)"
        }
        ##
        ## Create plot:
        ##
        p <- plot_data %>%
          pivot_longer(cols = all_of(coverage_col),
                       names_to = "metric",
                       values_to = "coverage") %>%
          mutate(metric = metric_type) %>%
          ggplot(aes(x = Model,
                     y = coverage, 
                     color = Model)) +
          geom_point(size = 4) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
          scale_color_brewer(palette = "Set1") +
          labs(y = "Coverage Probability",
               title = paste0("Coverage: ", title_suffix, " - Test: ", unique(plot_data$test)[1]),
               subtitle = "Red line = Nominal 95% coverage") +
          theme_bw(base_size = base_size) +
          theme(legend.position = "none",
                axis.text.x   = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
                axis.text.y   = element_text(size = base_size * 0.9),
                axis.title.x  = element_text(size = base_size * 1.1),
                axis.title.y  = element_text(size = base_size * 1.1),
                plot.title    = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                ##
                strip.text = element_text(size = base_size * 0.65, face = "bold"),  # <- Make this smaller
                strip.background = element_rect(fill = "grey90"),
                ##
                axis.ticks = element_line(size = 0.5)
          ) +
          scale_x_discrete(limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          scale_color_brewer(palette = "Set1",
                             limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          facet_grid(DGM ~ n_studies, 
                     scales = "free_y",
                     labeller = labeller(
                       n_studies = function(x) paste(x, "studies"),
                       DGM = function(x) {
                         case_when(
                           x == "Jones"        ~ "DGM: Jones",
                           x == "Xu_fixed_C"   ~ "DGM: Ord-Biv-FC",
                           x == "Xu_random_C"  ~ "DGM: Ord-Biv-RC",
                           x == "R&G_random_C" ~ "DGM: Ord-HSROC-RC",
                           TRUE ~ x
                         )
                       }
                     ))
        
        return(p)
  
}

##
## ---- Function to create Width plots for a specific DGM ---------------------------------------------------
##
create_width_plot <- function(data, 
                                 DGM_name = NULL,  # Make optional - if NULL, plot all DGMs
                                 metric_type = c("Se", "Sp"),
                                 base_size = 16
) {
  
        require(dplyr)
        require(ggplot2)
        require(purrr)
        require(tidyr)
        require(patchwork)
        ##
        metric_type <- match.arg(metric_type)
        
        ## Filter for specific DGM only if provided
        if (!is.null(DGM_name)) {
          plot_data <- data %>% filter(DGM == DGM_name)
        } else {
          plot_data <- data
        }
        ##
        ## Set up column names based on metric type
        ##
        if (metric_type == "Se") {
          width_col <- "avg_width_Se"
          title_suffix <- "Sensitivity (Se)"
        } else {
          width_col <- "avg_width_Sp"
          title_suffix <- "Specificity (Sp)"
        }
        ##
        ## Create plot:
        ##
        p <- plot_data %>%
          pivot_longer(cols = all_of(width_col),
                       names_to = "metric",
                       values_to = "width") %>%
          mutate(metric = metric_type) %>%
          ggplot(aes(x = Model,
                     y = width, 
                     color = Model)) +
          geom_point(size = 4) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
          scale_color_brewer(palette = "Set1") +
          labs(y = "Interval width",
               title = paste0("Width: ", title_suffix, " - Test: ", unique(plot_data$test)[1]),
               subtitle = "Red line = Nominal 95% width") +
          theme_bw(base_size = base_size) +
          theme(legend.position = "none",
                axis.text.x   = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
                axis.text.y   = element_text(size = base_size * 0.9),
                axis.title.x  = element_text(size = base_size * 1.1),
                axis.title.y  = element_text(size = base_size * 1.1),
                plot.title    = element_text(size = base_size * 1.3, face = "bold"),
                plot.subtitle = element_text(size = base_size * 1.0),
                ##
                strip.text = element_text(size = base_size * 0.65, face = "bold"),  # <- Make this smaller
                strip.background = element_rect(fill = "grey90"),
                ##
                axis.ticks = element_line(size = 0.5)
          ) +
          scale_x_discrete(limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          scale_color_brewer(palette = "Set1",
                             limits = c("Bivariate", "Jones-FC", "Ord-Biv-FC", "Ord-Biv-RC", "Ord-HSROC-RC")) +
          facet_grid(DGM ~ n_studies, 
                     scales = "free_y",
                     labeller = labeller(
                       n_studies = function(x) paste(x, "studies"),
                       DGM = function(x) {
                         case_when(
                           x == "Jones"        ~ "DGM: Jones",
                           x == "Xu_fixed_C"   ~ "DGM: Ord-Biv-FC",
                           x == "Xu_random_C"  ~ "DGM: Ord-Biv-RC",
                           x == "R&G_random_C" ~ "DGM: Ord-HSROC-RC",
                           TRUE ~ x
                         )
                       }
                     ))
        
        return(p)
  
}



R_fn_edit_sim_data_tibble <- function(data_0) { 
  
      data <- data_0 %>% 
        dplyr::mutate( model_parameterisation = case_when( model_parameterisation == "R&G" ~ "Ord-HSROC",
                                                           model_parameterisation == "Xu"  ~ "Ord-Biv",
                                                           TRUE ~ model_parameterisation ))
      data
      ##
      data <-  data %>% 
               dplyr::mutate( random_thr = case_when( random_thresholds == TRUE  ~ "-RC", 
                                                      random_thresholds == FALSE ~ "-FC",
                                                      TRUE ~ model_parameterisation),
                              Model = paste0(model_parameterisation, random_thr))
      ##
      data <- data %>% mutate(Model = ifelse(model_parameterisation == "Bivariate", "Bivariate", Model))
      ##
      data$Model <- factor(data$Model, levels = c(  "Bivariate", 
                                                    "Jones-FC",
                                                    "Ord-Biv-FC",
                                                    "Ord-Biv-RC",
                                                    "Ord-HSROC-RC"))
      ##
      data$DGM <- factor(data$DGM, levels = c(  "Jones", 
                                                "Xu_fixed_C", 
                                                "Xu_random_C", 
                                                "R&G_random_C"))
      
      return(data)
      
      
}








count_non_missing_by_column <- function(data, 
                                        missing_value = -1, 
                                        skip_first_col = TRUE,
                                        return_df = FALSE
) {
        
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
  
  
  

create_summary_df <- function(n_sims, 
                              min_k, 
                              max_k,
                              start_index = 1
) {

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
        
        # Fill the threshold column namesthat
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




# 
# compute_simulation_metrics <- function(test,
#                                        complete_results_list,
#                                        output_dir,
#                                        true_DGM_Se,
#                                        true_DGM_Sp,
#                                        min_ESS_thr,
#                                        vec_index_inner
# ) {
#   
#         min_k <- head(vec_index_inner, 1)
#         max_k <- tail(vec_index_inner, 1)
#         ##
#         datasets_sim_list <- complete_results_list ## readRDS(data_list_file)
#   
#         # Filter out NULL results AND incomplete results (low ESS)
#         valid_results <- complete_results_list[!sapply(complete_results_list, is.null)]
#         
#         # Only keep results that have credible intervals (i.e., successful runs)
#         complete_results <- valid_results[sapply(valid_results, function(x) {
#           !is.null(x$Se_lower_vec) && !is.null(x$Se_upper_vec) && 
#             length(x$Se_lower_vec) > 0 && length(x$Se_upper_vec) > 0
#         })]
#         
#         # *** NEW: Filter by ESS BEFORE main processing ***
#         complete_results_high_ESS <- Filter(function(result) {
#           
#             tryCatch({
#               Se_Sp_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch(
#                 tibble = result$tibble_all, 
#                 param_strings_vec = c("Se_baseline", "Sp_baseline"))
#               min_ESS <- min(Se_Sp_tibble$n_eff, na.rm = TRUE)
#               return(min_ESS >= min_ESS_thr)
#             }, error = function(e) {
#               return(FALSE)
#             })
#             
#         }, complete_results)
#         
#         n_studies <- valid_results[[1]]$n_studies ; n_studies
#         
#         n_sims <- length(complete_results_high_ESS)
#         n_excluded <- length(valid_results) - n_sims
#         
#         if(n_sims == 0) stop("No valid simulation results with adequate ESS found!")
#         
#         cat(sprintf("Using %d complete results (excluded %d with low ESS)\n", 
#                     n_sims, n_excluded))
#         
#         # Now use the filtered list
#         complete_results_all_ESS <- complete_results_high_ESS
#         n_studies <- complete_results_all_ESS[[1]]$n_studies
#         
#         non_missing_index_list <- missing_index_list <- list()
#         
#         for (i in 1:n_sims) {
#           
#             x <- valid_results[[i]]$x
#             
#             non_missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) >= 2)
#             missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) < 2)
#           
#         }
#         
#         # Initialize storage matrices
#         coverage_Se <- coverage_Sp <- matrix(NA, n_sims, max_k)
#         interval_width_Se <- interval_width_Sp <- matrix(NA, n_sims, max_k)
#         bias_Se <- bias_Sp <- matrix(NA, n_sims, max_k)
#         squared_error_Se <- squared_error_Sp <- matrix(NA, n_sims, max_k)
#         ##
#         Se_median_vec <- Sp_median_vec <- c()
#         Se_mean_vec   <- Sp_mean_vec <- c()
#         Se_lower_vec  <- Sp_lower_vec <- c()
#         Se_upper_vec  <- Sp_upper_vec <- c()
#         Se_pred_lower_vec <- Sp_pred_lower_vec <- c()
#         Se_pred_upper_vec <- Sp_pred_upper_vec <- c()
#         
#         # Store point estimates for SD calculation
#         point_est_Se <- point_est_Sp <- matrix(NA, n_sims, max_k)
#         ##
#         min_ESS_vec <- c()
#         
#         for (i in 1:n_sims) {
#           
#           try({ 
#           
#               result <- complete_results_all_ESS[[i]]
#               ##
#               Se_Sp_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch( tibble = result$tibble_all, 
#                                                                              param_strings_vec = c("Se_baseline", "Sp_baseline"))
#               min_ESS <- min(Se_Sp_tibble$n_eff, na.rm = TRUE)
#               min_ESS_vec[i] <- min_ESS
#               
#           })
#               
#         }
#         ##
#         mean_min_ESS <- round(mean(min_ESS_vec))
#         
#         # Loop through simulations (only complete ones now)
#         for (i in 1:n_sims) {
#           
#                 result <- complete_results_high_ESS[[i]]
#                 
#                 
#                 # Se_Sp_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch( 
#                 #                              tibble = result$tibble_all, 
#                 #                              param_strings_vec = c("Se_baseline", "Sp_baseline"))
#                 # min_ESS <- min(Se_Sp_tibble$n_eff, na.rm = TRUE)
#                 # min_ESS_vec[i] <- min_ESS
#                 
#                 ## replace values w/ less than 2 studies with NA's (otherwise can't compare to stratified bivariate fairly):
#                 result$Se_median_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_median_vec[missing_index_list[[i]]] <- NA
#                 ##
#                 result$Se_mean_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_mean_vec[missing_index_list[[i]]] <- NA
#                 ##
#                 result$Se_lower_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_lower_vec[missing_index_list[[i]]] <- NA
#                 ##
#                 result$Se_upper_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_upper_vec[missing_index_list[[i]]] <- NA
#                 ##
#                 result$Se_pred_lower_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_pred_lower_vec[missing_index_list[[i]]] <- NA
#                 ##
#                 result$Se_pred_upper_vec[missing_index_list[[i]]] <- NA
#                 result$Sp_pred_upper_vec[missing_index_list[[i]]] <- NA
#                 
#                 for (k in min_k:max_k) {
#                   
#                         Se_median_vec[k] <- result$Se_median_vec[k]
#                         Sp_median_vec[k] <- result$Sp_median_vec[k]
#                         ##
#                         Se_mean_vec[k] <- result$Se_mean_vec[k]
#                         Sp_mean_vec[k] <- result$Sp_mean_vec[k]
#                         ##
#                         Se_lower_vec[k] <- result$Se_lower_vec[k]
#                         Sp_lower_vec[k] <- result$Sp_lower_vec[k]
#                         ##
#                         Se_upper_vec[k] <- result$Se_upper_vec[k]
#                         Sp_upper_vec[k] <- result$Sp_upper_vec[k]
#                         ##
#                         Se_pred_lower_vec[k] <- result$Se_pred_lower_vec[k]
#                         Sp_pred_lower_vec[k] <- result$Sp_pred_lower_vec[k]
#                         ##
#                         Se_pred_upper_vec[k] <- result$Se_pred_upper_vec[k]
#                         Sp_pred_upper_vec[k] <- result$Sp_pred_upper_vec[k]
#                         
#                         # Store point estimates
#                         point_est_Se[i, k] <- result$Se_median_vec[k]
#                         point_est_Sp[i, k] <- result$Sp_median_vec[k]
#                         
#                         # Check coverage (is true value in CI?)
#                         coverage_Se[i, k] <- (true_DGM_Se[k] >= result$Se_lower_vec[k]) & 
#                                              (true_DGM_Se[k] <= result$Se_upper_vec[k])
#                         coverage_Sp[i, k] <- (true_DGM_Sp[k] >= result$Sp_lower_vec[k]) & 
#                                              (true_DGM_Sp[k] <= result$Sp_upper_vec[k])
#                         
#                         # Calculate interval width
#                         interval_width_Se[i, k] <- result$Se_upper_vec[k] - result$Se_lower_vec[k]
#                         interval_width_Sp[i, k] <- result$Sp_upper_vec[k] - result$Sp_lower_vec[k]
#                         
#                         # Calculate bias (using median as point estimate)
#                         bias_Se[i, k] <- result$Se_median_vec[k] - true_DGM_Se[k]
#                         bias_Sp[i, k] <- result$Sp_median_vec[k] - true_DGM_Sp[k]
#                         
#                         # Calculate squared error for RMSE
#                         squared_error_Se[i, k] <- (result$Se_median_vec[k] - true_DGM_Se[k])^2
#                         squared_error_Sp[i, k] <- (result$Sp_median_vec[k] - true_DGM_Sp[k])^2
#                   
#                 }
#           
#         }
#         
#        
#         
#         # Compute summary metrics for each threshold
#         results_by_threshold <- list()
#         ##
#         bias_Se_vec <- c()
#         bias_Sp_vec <- c()
#         abs_bias_Se_vec <- c()
#         abs_bias_Sp_vec <- c()
#         MCSE_bias_Se_vec <- c()
#         MCSE_bias_Sp_vec <- c()
#         ##
#         SD_Se_vec <- c()
#         SD_Sp_vec <- c()
#         ##
#         RMSE_Se_vec <- c()
#         RMSE_Sp_vec <- c()
#         MCSE_RMSE_Se_vec <- c()
#         MCSE_RMSE_Sp_vec <- c()
#         
#         
#         for (k in min_k:max_k) {
#           
#           # Calculate bias metrics
#           mean_bias_Se_k <- mean(bias_Se[, k], na.rm = TRUE)
#           mean_bias_Sp_k <- mean(bias_Sp[, k], na.rm = TRUE)
#           
#           # Calculate SD of point estimates across simulations
#           SD_Se_vec[k] <- sd(point_est_Se[, k], na.rm = TRUE)
#           SD_Sp_vec[k] <- sd(point_est_Sp[, k], na.rm = TRUE)
#           
#           # Store for overall calculations
#           bias_Se_vec[k] <- mean_bias_Se_k
#           bias_Sp_vec[k] <- mean_bias_Sp_k
#           ##
#           abs_bias_Se_vec[k] <- abs(mean_bias_Se_k)
#           abs_bias_Sp_vec[k] <- abs(mean_bias_Sp_k)
#           
#           # MCSE of bias
#           MCSE_bias_Se_vec[k] <- SD_Se_vec[k] / sqrt(n_sims)
#           MCSE_bias_Sp_vec[k] <- SD_Sp_vec[k] / sqrt(n_sims)
#           
#           # RMSE calculation
#           RMSE_Se_k <- sqrt(mean(squared_error_Se[, k], na.rm = TRUE))
#           RMSE_Sp_k <- sqrt(mean(squared_error_Sp[, k], na.rm = TRUE))
#           
#           RMSE_Se_vec[k] <- RMSE_Se_k
#           RMSE_Sp_vec[k] <- RMSE_Sp_k
#           
#           # MCSE of RMSE using Delta method
#           var_squared_error_Se <- var(squared_error_Se[, k], na.rm = TRUE)
#           var_squared_error_Sp <- var(squared_error_Sp[, k], na.rm = TRUE)
#           
#           MCSE_RMSE_Se_k <- ifelse(RMSE_Se_k > 0, 
#                                    sqrt(var_squared_error_Se) / (2 * sqrt(n_sims) * RMSE_Se_k),
#                                    NA)
#           MCSE_RMSE_Sp_k <- ifelse(RMSE_Sp_k > 0,
#                                    sqrt(var_squared_error_Sp) / (2 * sqrt(n_sims) * RMSE_Sp_k),
#                                    NA)
#           
#           MCSE_RMSE_Se_vec[k] <- MCSE_RMSE_Se_k
#           MCSE_RMSE_Sp_vec[k] <- MCSE_RMSE_Sp_k
#           
#           results_by_threshold[[k]] <- list(
#             # Coverage probability
#             coverage_prob_Se = mean(coverage_Se[, k], na.rm = TRUE),
#             coverage_prob_Sp = mean(coverage_Sp[, k], na.rm = TRUE),
#             ##
#             # Mean interval width
#             mean_interval_width_Se = mean(interval_width_Se[, k], na.rm = TRUE),
#             mean_interval_width_Sp = mean(interval_width_Sp[, k], na.rm = TRUE),
#             ##
#             # RMSE
#             RMSE_Se = RMSE_Se_k,
#             RMSE_Sp = RMSE_Sp_k,
#             ##
#             # MCSE of RMSE
#             MCSE_RMSE_Se = MCSE_RMSE_Se_k,
#             MCSE_RMSE_Sp = MCSE_RMSE_Sp_k,
#             ##
#             # Mean bias
#             mean_bias_Se = mean_bias_Se_k,
#             mean_bias_Sp = mean_bias_Sp_k,
#             ##
#             # Mean absolute bias
#             mean_abs_bias_Se = abs(mean_bias_Se_k),
#             mean_abs_bias_Sp = abs(mean_bias_Sp_k),
#             ##
#             # MCSE of bias
#             MCSE_bias_Se = MCSE_bias_Se_vec[k],
#             MCSE_bias_Sp = MCSE_bias_Sp_vec[k],
#             ##
#             # SD
#             SD_Se = SD_Se_vec[k],
#             SD_Sp = SD_Sp_vec[k],
#             ##
#             ##
#             ##
#             Se_median = Se_median_vec[k],
#             Sp_median = Sp_median_vec[k],
#             ##
#             Se_mean = Se_mean_vec[k],
#             Sp_mean = Sp_mean_vec[k],
#             ##
#             Se_lower = Se_lower_vec[k],
#             Sp_lower = Sp_lower_vec[k],
#             ##
#             Se_upper = Se_upper_vec[k],
#             Sp_upper = Sp_upper_vec[k],
#             ##
#             Se_pred_lower = Se_pred_lower_vec[k],
#             Sp_pred_lower = Sp_pred_lower_vec[k],
#             ##
#             Se_pred_upper = Se_pred_upper_vec[k],
#             Sp_pred_upper = Sp_pred_upper_vec[k]
#           )
#           
#         }
#         
#         # vec_index <- min_k:max_k
#         
#         # Overall summaries (averaged across thresholds min_k:max_k)
#         overall_results <- list(
#               # Number of excluded results
#               n_excluded = n_excluded,
#               
#               # Average RMSE across thresholds
#               avg_RMSE_Se = mean(RMSE_Se_vec, na.rm = TRUE),
#               avg_RMSE_Sp = mean(RMSE_Sp_vec, na.rm = TRUE),
#               
#               # Maximum RMSE
#               max_RMSE_Se = max(RMSE_Se_vec, na.rm = TRUE),
#               max_RMSE_Sp = max(RMSE_Sp_vec, na.rm = TRUE),
#               
#               # Average MCSE of RMSE
#               avg_MCSE_RMSE_Se = mean(MCSE_RMSE_Se_vec, na.rm = TRUE),
#               avg_MCSE_RMSE_Sp = mean(MCSE_RMSE_Sp_vec, na.rm = TRUE),
#               
#               # MCSE of max RMSE (using the MCSE at the threshold where max occurs)
#               max_RMSE_Se_index = which.max(RMSE_Se_vec),
#               max_RMSE_Sp_index = which.max(RMSE_Sp_vec),
#               max_MCSE_RMSE_Se = MCSE_RMSE_Se_vec[which.max(RMSE_Se_vec)],
#               max_MCSE_RMSE_Sp = MCSE_RMSE_Sp_vec[which.max(RMSE_Sp_vec)],
#               
#               # Average absolute bias across thresholds
#               avg_abs_bias_Se = mean(abs_bias_Se_vec, na.rm = TRUE),
#               avg_abs_bias_Sp = mean(abs_bias_Sp_vec, na.rm = TRUE),
#               
#               # Maximum absolute bias
#               max_abs_bias_Se = max(abs_bias_Se_vec, na.rm = TRUE),
#               max_abs_bias_Sp = max(abs_bias_Sp_vec, na.rm = TRUE),
#               
#               # Average MCSE of bias
#               avg_MCSE_bias_Se = mean(MCSE_bias_Se_vec, na.rm = TRUE),
#               avg_MCSE_bias_Sp = mean(MCSE_bias_Sp_vec, na.rm = TRUE),
#               
#               # Maximum MCSE of bias
#               max_MCSE_bias_Se = max(MCSE_bias_Se_vec, na.rm = TRUE),
#               max_MCSE_bias_Sp = max(MCSE_bias_Sp_vec, na.rm = TRUE),
#               
#               # Average coverage across thresholds
#               avg_coverage_prob_Se = mean(sapply(results_by_threshold, function(x) x$coverage_prob_Se), na.rm = TRUE),
#               avg_coverage_prob_Sp = mean(sapply(results_by_threshold, function(x) x$coverage_prob_Sp), na.rm = TRUE),
#               
#               # Average interval width across thresholds
#               avg_interval_width_Se = mean(sapply(results_by_threshold, function(x) x$mean_interval_width_Se), na.rm = TRUE),
#               avg_interval_width_Sp = mean(sapply(results_by_threshold, function(x) x$mean_interval_width_Sp), na.rm = TRUE),
#               
#               # Average SD
#               avg_SD_Se = mean(SD_Se_vec, na.rm = TRUE),
#               avg_SD_Sp = mean(SD_Sp_vec, na.rm = TRUE),
#               max_SD_Se = max(SD_Se_vec, na.rm = TRUE),
#               max_SD_Sp = max(SD_Sp_vec, na.rm = TRUE)
#             )
#         
#         return(list(
#           mean_min_ESS = mean_min_ESS,
#           by_threshold = results_by_threshold,
#           overall = overall_results,
#           n_sims = n_sims,
#           bias_vectors = list(
#             bias_Se_vec = bias_Se_vec,
#             bias_Sp_vec = bias_Sp_vec,
#             abs_bias_Se_vec = abs_bias_Se_vec,
#             abs_bias_Sp_vec = abs_bias_Sp_vec,
#             MCSE_bias_Se_vec = MCSE_bias_Se_vec,
#             MCSE_bias_Sp_vec = MCSE_bias_Sp_vec,
#             SD_Se_vec = SD_Se_vec,
#             SD_Sp_vec = SD_Sp_vec,
#             RMSE_Se_vec = RMSE_Se_vec,
#             RMSE_Sp_vec = RMSE_Sp_vec,
#             MCSE_RMSE_Se_vec = MCSE_RMSE_Se_vec,
#             MCSE_RMSE_Sp_vec = MCSE_RMSE_Sp_vec
#           )
#         ))
#   
# }




compute_simulation_metrics <- function(test,
                                       complete_results_list,
                                       output_dir,
                                       true_DGM_Se,
                                       true_DGM_Sp,
                                       min_ESS_thr,
                                       vec_index_inner,
                                       min_studies_per_threshold = NULL
) {
        
        min_k <- head(vec_index_inner, 1)
        max_k <- tail(vec_index_inner, 1)
        ##
        datasets_sim_list <- complete_results_list ## readRDS(data_list_file)
        
        # Filter out NULL results AND incomplete results (low ESS)
        valid_results <- complete_results_list[!sapply(complete_results_list, is.null)]
        
        # Only keep results that have credible intervals (i.e., successful runs)
        complete_results <- valid_results[sapply(valid_results, function(x) {
          !is.null(x$Se_lower_vec) && !is.null(x$Se_upper_vec) && 
            length(x$Se_lower_vec) > 0 && length(x$Se_upper_vec) > 0
        })]
        
        # *** NEW: Filter by ESS BEFORE main processing ***
        complete_results_high_ESS <- Filter(function(result) {
          
          tryCatch({
            Se_Sp_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch(
              tibble = result$tibble_all, 
              param_strings_vec = c("Se_baseline", "Sp_baseline"))
            min_ESS <- min(Se_Sp_tibble$n_eff, na.rm = TRUE)
            return(min_ESS >= min_ESS_thr)
          }, error = function(e) {
            return(FALSE)
          })
          
        }, complete_results)
        
        n_studies <- valid_results[[1]]$n_studies ; n_studies
        
        n_sims <- length(complete_results_high_ESS)
        n_excluded <- length(valid_results) - n_sims
        
        # CHANGED: Return empty results instead of stopping
        if(n_sims == 0) {
          cat(sprintf("WARNING: No valid simulation results with adequate ESS found! (excluded %d)\n", 
                      n_excluded))
          
          # Return empty structure that matches expected output
          return(list(
            mean_min_ESS = NA,
            by_threshold = list(),
            overall = list(
              n_excluded = n_excluded,
              avg_RMSE_Se = NA,
              avg_RMSE_Sp = NA,
              max_RMSE_Se = NA,
              max_RMSE_Sp = NA,
              avg_MCSE_RMSE_Se = NA,
              avg_MCSE_RMSE_Sp = NA,
              max_MCSE_RMSE_Se = NA,
              max_MCSE_RMSE_Sp = NA,
              avg_abs_bias_Se = NA,
              avg_abs_bias_Sp = NA,
              max_abs_bias_Se = NA,
              max_abs_bias_Sp = NA,
              avg_MCSE_bias_Se = NA,
              avg_MCSE_bias_Sp = NA,
              max_MCSE_bias_Se = NA,
              max_MCSE_bias_Sp = NA,
              avg_coverage_prob_Se = NA,
              avg_coverage_prob_Sp = NA,
              avg_interval_width_Se = NA,
              avg_interval_width_Sp = NA,
              avg_SD_Se = NA,
              avg_SD_Sp = NA,
              max_SD_Se = NA,
              max_SD_Sp = NA,
              included_thresholds = numeric(0),
              n_valid_thresholds = 0
            ),
            n_sims = 0,
            bias_vectors = list(
              bias_Se_vec = numeric(0),
              bias_Sp_vec = numeric(0),
              abs_bias_Se_vec  = numeric(0),
              abs_bias_Sp_vec  = numeric(0),
              MCSE_bias_Se_vec = numeric(0),
              MCSE_bias_Sp_vec = numeric(0),
              SD_Se_vec = numeric(0),
              SD_Sp_vec = numeric(0),
              RMSE_Se_vec = numeric(0),
              RMSE_Sp_vec = numeric(0),
              MCSE_RMSE_Se_vec = numeric(0),
              MCSE_RMSE_Sp_vec = numeric(0)
            )
          ))
        }
        
        cat(sprintf("Using %d complete results (excluded %d with low ESS)\n", 
                    n_sims, n_excluded))
        
        if(n_sims == 0) stop("No valid simulation results with adequate ESS found!")
        
        # cat(sprintf("Using %d complete results (excluded %d with low ESS)\n", 
        #             n_sims, n_excluded))
        
        # Now use the filtered list
        complete_results_all_ESS <- complete_results_high_ESS
        n_studies <- complete_results_all_ESS[[1]]$n_studies
        
        non_missing_index_list <- missing_index_list <- list()
        
        # NEW: Determine minimum threshold for missing data
        min_non_missing_threshold <- ifelse(!is.null(min_studies_per_threshold), 
                                            min_studies_per_threshold, 
                                            2)  # Default to 2 if not specified
        
        for (i in 1:n_sims) {
          
          x <- valid_results[[i]]$x
          
          # UPDATED: Use the configurable threshold
          non_missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) >= min_non_missing_threshold)
          missing_index_list[[i]] <- which(count_non_missing_by_column(x[[1]]) < min_non_missing_threshold)
          
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
        ##
        min_ESS_vec <- c()
        
        for (i in 1:n_sims) {
          
          try({ 
            
            result <- complete_results_all_ESS[[i]]
            ##
            Se_Sp_tibble <- MetaOrdDTA:::extract_params_from_tibble_batch( tibble = result$tibble_all, 
                                                                           param_strings_vec = c("Se_baseline", "Sp_baseline"))
            min_ESS <- min(Se_Sp_tibble$n_eff, na.rm = TRUE)
            min_ESS_vec[i] <- min_ESS
            
          })
          
        }
        ##
        mean_min_ESS <- round(mean(min_ESS_vec))
        
        # Loop through simulations (only complete ones now)
        for (i in 1:n_sims) {
          
          result <- complete_results_high_ESS[[i]]
          
          ## UPDATED COMMENT: replace values with less than min_studies_per_threshold with NA's
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
        
        # NEW: Create a vector of valid thresholds to include in overall summaries
        # Only include thresholds where we have enough data
        valid_thresholds <- c()
        for (k in min_k:max_k) {
          n_valid_k <- sum(!is.na(point_est_Se[, k]))
          if (n_valid_k > 0) {  # At least some valid data for this threshold
            valid_thresholds <- c(valid_thresholds, k)
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
          ##
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
        
        # UPDATED: Overall summaries only for valid thresholds
        overall_results <- list(
          # Number of excluded results
          n_excluded = n_excluded,
          
          # Average RMSE across VALID thresholds only
          avg_RMSE_Se = mean(RMSE_Se_vec[valid_thresholds], na.rm = TRUE),
          avg_RMSE_Sp = mean(RMSE_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Maximum RMSE (from valid thresholds)
          max_RMSE_Se = max(RMSE_Se_vec[valid_thresholds], na.rm = TRUE),
          max_RMSE_Sp = max(RMSE_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Average MCSE of RMSE
          avg_MCSE_RMSE_Se = mean(MCSE_RMSE_Se_vec[valid_thresholds], na.rm = TRUE),
          avg_MCSE_RMSE_Sp = mean(MCSE_RMSE_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # MCSE of max RMSE (using the MCSE at the threshold where max occurs)
          max_RMSE_Se_index = valid_thresholds[which.max(RMSE_Se_vec[valid_thresholds])],
          max_RMSE_Sp_index = valid_thresholds[which.max(RMSE_Sp_vec[valid_thresholds])],
          max_MCSE_RMSE_Se = MCSE_RMSE_Se_vec[valid_thresholds[which.max(RMSE_Se_vec[valid_thresholds])]],
          max_MCSE_RMSE_Sp = MCSE_RMSE_Sp_vec[valid_thresholds[which.max(RMSE_Sp_vec[valid_thresholds])]],
          
          # Average absolute bias across VALID thresholds
          avg_abs_bias_Se = mean(abs_bias_Se_vec[valid_thresholds], na.rm = TRUE),
          avg_abs_bias_Sp = mean(abs_bias_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Maximum absolute bias
          max_abs_bias_Se = max(abs_bias_Se_vec[valid_thresholds], na.rm = TRUE),
          max_abs_bias_Sp = max(abs_bias_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Average MCSE of bias
          avg_MCSE_bias_Se = mean(MCSE_bias_Se_vec[valid_thresholds], na.rm = TRUE),
          avg_MCSE_bias_Sp = mean(MCSE_bias_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Maximum MCSE of bias
          max_MCSE_bias_Se = max(MCSE_bias_Se_vec[valid_thresholds], na.rm = TRUE),
          max_MCSE_bias_Sp = max(MCSE_bias_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # Average coverage across VALID thresholds
          avg_coverage_prob_Se = mean(sapply(results_by_threshold[valid_thresholds], function(x) x$coverage_prob_Se), na.rm = TRUE),
          avg_coverage_prob_Sp = mean(sapply(results_by_threshold[valid_thresholds], function(x) x$coverage_prob_Sp), na.rm = TRUE),
          
          # Average interval width across VALID thresholds
          avg_interval_width_Se = mean(sapply(results_by_threshold[valid_thresholds], function(x) x$mean_interval_width_Se), na.rm = TRUE),
          avg_interval_width_Sp = mean(sapply(results_by_threshold[valid_thresholds], function(x) x$mean_interval_width_Sp), na.rm = TRUE),
          
          # Average SD
          avg_SD_Se = mean(SD_Se_vec[valid_thresholds], na.rm = TRUE),
          avg_SD_Sp = mean(SD_Sp_vec[valid_thresholds], na.rm = TRUE),
          max_SD_Se = max(SD_Se_vec[valid_thresholds], na.rm = TRUE),
          max_SD_Sp = max(SD_Sp_vec[valid_thresholds], na.rm = TRUE),
          
          # NEW: Track which thresholds were included
          included_thresholds = valid_thresholds,
          n_valid_thresholds = length(valid_thresholds)
        )
        
        return(list(
          mean_min_ESS = mean_min_ESS,
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

create_metrics_summary_table <- function( metrics, 
                                          vec_index_inner
) {
  
          min_k <- head(vec_index_inner, 1)
          max_k <- tail(vec_index_inner, 1)
          ##
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






calc_min_sim_size_RMSE <- function( metrics, 
                                    target_MCSE_RMSE,
                                    vec_index_inner
) {
  
        min_k <- head(vec_index_inner, 1)
        max_k <- tail(vec_index_inner, 1)
        ##
        min_N_sim <- c()
        index <- 1
        ##
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
        ##
        min_required = min(min_N_sim, na.rm = TRUE)
        avg_required = mean(min_N_sim, na.rm = TRUE)
        max_required = max(min_N_sim, na.rm = TRUE)
        ##
        n_sims_required <- round((max_required + avg_required)/2.0)
        ##
        return(list(
          min_required = min_required,
          avg_required = avg_required,
          max_required = max_required,
          n_sims_required = n_sims_required
        ))
  
}







 
 




