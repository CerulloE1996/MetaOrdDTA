







#' R_fn_sROC_plot_NMA
#' @keywords internal
#' @export
R_fn_sROC_plot_NMA <- function( stan_model_file_name,
                                ##
                                tibbles_provided = FALSE,
                                tibbles = NULL,
                                stan_mod_samples = NULL,
                                ##
                                # covariates,
                                # new_cov_data,
                                ##
                                df_true = NULL,
                                ##
                                conf_region_colour = "blue", 
                                pred_region_colour = "green",
                                ##
                                conf_region_alpha = 0.15,
                                pred_region_alpha = 0.30,
                                ##
                                base_size = 16,
                                point_size = 3,
                                linewidth = 0.5,
                                ##
                                n_index_tests,
                                n_thr,
                                ##
                                test_names = NULL,
                                relevant_thresholds = NULL
) {

                require(ggplot2)
                require(dplyr)
                require(patchwork)
                ##
                model_name <- stan_model_file_name
                n_thr_max <- max(n_thr)
                ##
                # Set default test names if not provided
                if (is.null(test_names)) {
                  test_names <- paste0("Test ", 1:n_index_tests)
                }
                ##
                Se_median_array <- Sp_median_array <- Fp_median_array <- array(NA, dim = c(n_index_tests, n_thr_max))
                Se_mean_array   <- Sp_mean_array   <- Fp_mean_array   <- array(NA, dim = c(n_index_tests, n_thr_max))
                ##
                Se_lower_array <- Sp_lower_array <- Fp_lower_array <- array(NA, dim = c(n_index_tests, n_thr_max))
                Se_upper_array <- Sp_upper_array <- Fp_upper_array <- array(NA, dim = c(n_index_tests, n_thr_max))
                ##
                Se_pred_lower_array <- Sp_pred_lower_array <- Fp_pred_lower_array <- array(NA, dim = c(n_index_tests, n_thr_max))
                Se_pred_upper_array <- Sp_pred_upper_array <- Fp_pred_upper_array <- array(NA, dim = c(n_index_tests, n_thr_max))
                ##
                ## Extract summaries as tibbles:
                ##
                if (tibbles_provided == TRUE) {
                  
                    Se <- tibbles$Se_summaries
                    Sp <- tibbles$Sp_summaries
                    Fp <- tibbles$Fp_summaries
                    ##
                    Se_pred <- tibbles$Se_pred_summaries
                    Sp_pred <- tibbles$Sp_pred_summaries
                    Fp_pred <- tibbles$Fp_pred_summaries
                    
                } else {
                  
                    Se <- stan_mod_samples$summary(c("Se_baseline"), mean, 
                                                   quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    Sp <- stan_mod_samples$summary(c("Sp_baseline"), mean, 
                                                   quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    Fp <- stan_mod_samples$summary(c("Fp_baseline"), mean, 
                                                   quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    ##
                    Se_pred <- stan_mod_samples$summary(c("Se_baseline_pred"), mean, 
                                                        quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    Sp_pred <- stan_mod_samples$summary(c("Sp_baseline_pred"), mean, 
                                                        quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    Fp_pred <- stan_mod_samples$summary(c("Fp_baseline_pred"), mean, 
                                                        quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
                    
                }
                ##
                counter <- 1
                for (k in 1:n_thr_max) {
                  for (t in 1:n_index_tests) {
                    # test_vec[counter] <- t
                    ## Posterior medians (of pooled estimates):
                    Se_median_array[t, k] <- Se$`50%`[counter]
                    Sp_median_array[t, k] <- Sp$`50%`[counter]
                    Fp_median_array[t, k] <- Fp$`50%`[counter]
                    ## Posterior means (of pooled estimates):
                    Se_mean_array[t, k] <- Se$mean[counter]
                    Sp_mean_array[t, k] <- Sp$mean[counter]
                    Fp_mean_array[t, k] <- Fp$mean[counter]
                    ## Posterior lower 95% (of pooled estimates):
                    Se_lower_array[t, k] <- Se$`2.5%`[counter]
                    Sp_lower_array[t, k] <- Sp$`2.5%`[counter]
                    Fp_lower_array[t, k] <- Fp$`2.5%`[counter]
                    ## Posterior upper 95% (of pooled estimates):
                    Se_upper_array[t, k] <- Se$`97.5%`[counter]
                    Sp_upper_array[t, k] <- Sp$`97.5%`[counter]
                    Fp_upper_array[t, k] <- Fp$`97.5%`[counter]
                    ## Posterior lower prediction 95% (of pooled estimates):
                    Se_pred_lower_array[t, k] <- Se_pred$`2.5%`[counter]
                    Sp_pred_lower_array[t, k] <- Sp_pred$`2.5%`[counter]
                    Fp_pred_lower_array[t, k] <- Fp_pred$`2.5%`[counter]
                    ## Posterior upper prediction 95% (of pooled estimates):
                    Se_pred_upper_array[t, k] <- Se_pred$`97.5%`[counter]
                    Sp_pred_upper_array[t, k] <- Sp_pred$`97.5%`[counter]
                    Fp_pred_upper_array[t, k] <- Fp_pred$`97.5%`[counter]
                    ##
                    counter <- counter + 1
                  }
                }
                ##
                test_vec <- thr_vec <- test_name_vec <- c()
                counter <- 1 
                for (t in 1:n_index_tests) {
                  for (k in 1:n_thr_max) {
                    test_vec[counter] <- t
                    test_name_vec[counter] <- test_names[t]  # Use actual test names
                    thr_vec[counter] <- k
                    
                    counter <- counter + 1
                  }
                }
                ##
                n_rows_total <- length(test_vec)
                ##
                ## Make tibble with test names
                tibble_NMA <- tibble(  Model = rep(model_name, n_rows_total),
                                       ##
                                       test = factor(test_vec),
                                       test_name = factor(test_name_vec, levels = test_names),  # Use actual test names
                                       threshold = thr_vec,
                                       ##
                                       Se_median = c(t(Se_median_array)),
                                       Sp_median = c(t(Sp_median_array)),
                                       Fp_median = c(t(Fp_median_array)),
                                       ##
                                       Se_mean = c(t(Se_mean_array)),
                                       Sp_mean = c(t(Sp_mean_array)),
                                       Fp_mean = c(t(Fp_mean_array)),
                                       ##
                                       Se_lower = c(t(Se_lower_array)),
                                       Sp_lower = c(t(Sp_lower_array)),
                                       Fp_lower = c(t(Fp_lower_array)),
                                       ##
                                       Se_upper = c(t(Se_upper_array)),
                                       Sp_upper = c(t(Sp_upper_array)),
                                       Fp_upper = c(t(Fp_upper_array)),
                                       ##
                                       Se_pred_lower = c(t(Se_pred_lower_array)),
                                       Sp_pred_lower = c(t(Sp_pred_lower_array)),
                                       Fp_pred_lower = c(t(Fp_pred_lower_array)),
                                       ##
                                       Se_pred_upper = c(t(Se_pred_upper_array)),
                                       Sp_pred_upper = c(t(Sp_pred_upper_array)),
                                       Fp_pred_upper = c(t(Fp_pred_upper_array)))
                ##
                tibble_NMA <- dplyr::filter(tibble_NMA, Se_median != -1)
                ##
                if (!(is.null(relevant_thresholds))) { 
                  tibble_NMA <-  filter_tibble_by_thresholds(  tibble_data = tibble_NMA,
                                                               relevant_thresholds = relevant_thresholds)
                }
                ##
                ## Update polygon creation to use test names
                polygon_Conf_list <- polygon_Pred_list <- list()
                
                for (t in 1:n_index_tests) {
                  
                    tibble_test_t <- dplyr::filter(tibble_NMA, test == t)
                    
                    # Confidence region
                    tibble_Conf_polygon_test_t <- create_confidence_polygon(tibble_test_t, 
                                                                            model_name = model_name)
                    n_rows = nrow(tibble_Conf_polygon_test_t)
                    tibble_Conf_polygon_test_t <- tibble_Conf_polygon_test_t %>% 
                      dplyr::mutate(test = rep(t, n_rows),
                                    test_name = factor(rep(test_names[t], n_rows), levels = test_names))
                    polygon_Conf_list[[t]] <- tibble_Conf_polygon_test_t
                    
                    # Prediction region
                    tibble_Pred_polygon_test_t <- create_prediction_region_smooth(tibble_test_t,
                                                                            model_name = model_name)
                    n_rows = nrow(tibble_Pred_polygon_test_t)
                    tibble_Pred_polygon_test_t <- tibble_Pred_polygon_test_t %>% 
                      dplyr::mutate(test = rep(t, n_rows), 
                                    test_name = factor(rep(test_names[t], n_rows), levels = test_names))
                    polygon_Pred_list[[t]] <- tibble_Pred_polygon_test_t
                    
                }
                
                polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% dplyr::filter(!(is.na(x)))
                polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% dplyr::filter(!(is.na(x)))
                
                  plot_list <- list()
                  
                  # Plot 1: All tests on one plot (use test names for color legend)
                  plot_1 <-  ggplot(tibble_NMA, 
                                    mapping = aes(x = Fp_median, 
                                                  y = Se_median,
                                                  colour = test_name)) + 
                    geom_line(linewidth = linewidth) + 
                    geom_point(size = point_size) + 
                    theme_bw(base_size = base_size) + 
                    xlab("False positive rate (Fp)") + 
                    ylab("Sensitivity (Se)") + 
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) + 
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    labs(colour = "Test")
                  print(plot_1)
                  plot_list[[1]] <- plot_1
                  
                  # Plot 2: Each test on separate panel (use test names for facets)
                  plot_2 <-  ggplot(tibble_NMA, 
                                    mapping = aes(x = Fp_median,
                                                  y = Se_median,
                                                  colour = Model)) + 
                    geom_line(linewidth = linewidth) + 
                    geom_point(size = point_size) + 
                    theme_bw(base_size = base_size) + 
                    facet_wrap(~ test_name) +  # Use test_name instead of test_char
                    xlab("False positive rate (Fp)") + 
                    ylab("Sensitivity (Se)") + 
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) + 
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) + 
                    theme(legend.position = "none")
                  print(plot_2)
                  plot_list[[2]] <- plot_2
                  
                  ##
                  ## ---- Plot 3 (plot w/ 95% confidence region):
                  ##
                  plot_3 <-   plot_2  +
                    ##
                    geom_polygon(data = polygon_Conf_tibble,
                                 aes(x = x, y = y, colour = Model),
                                 fill = conf_region_colour,
                                 alpha = conf_region_alpha)  +
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1))
                  # xlab("False positive rate (Fp)") +
                  # ylab("Sensitivity (Se)")
                  # geom_polygon(data = polygon_Pred_tibble, aes(x = x, y = y), fill = pred_region_colour, alpha = 0.25) +
                  ##
                  print(plot_3)
                  plot_list[[3]] <- plot_3
                  ##
                  ## ---- Plot 4 (plot w/ 95% prediction region):
                  ##
                  plot_4 <-  plot_2 +  
                    geom_polygon(
                                data = polygon_Pred_tibble,
                                aes(x = x, y = y, colour = Model),
                                fill = pred_region_colour,
                                alpha = conf_region_alpha)  +
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1))
                  # xlab("False positive rate (Fp)") +
                  # ylab("Sensitivity (Se)")
                  print(plot_4)
                  plot_list[[4]] <- plot_4
                  ##
                  ## ---- Plot 5 (plot w/ BOTH the 95% confidence + prediction regions):
                  ##
                  plot_5 <-  plot_2 +
                    geom_polygon(
                                data = polygon_Conf_tibble,
                                aes(x = x, y = y, colour = Model),
                                fill = conf_region_colour,
                                alpha = conf_region_alpha) +
                    ##
                    geom_polygon(
                      data = polygon_Pred_tibble,
                      aes(x = x, y = y, colour = Model),
                      fill = pred_region_colour,
                      alpha = pred_region_alpha) +
                    # geom_ribbon(aes(x = Fp_pred_lower, 
                    #                 ymin = Se_pred_lower, 
                    #                 ymax = Se_pred_upper),
                    #             alpha = 0.2, 
                    #             fill = conf_region_colour) +
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1))
                  # xlab("False positive rate (Fp)") +
                  # ylab("Sensitivity (Se)")
                  print(plot_5)
                  plot_list[[5]] <- plot_5
                  ##
                  ## ---- Plot 6 (plot w/ BOTH the 95% confidence + prediction regions):
                  ##
                  plot_6 <-  plot_1 +
                    geom_polygon(
                      data = polygon_Conf_tibble,
                      aes(x = x, 
                          y = y,
                          # colour = test_name,
                          fill = test_name),
                      alpha = pred_region_alpha) +
                    ##
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1))
                  # xlab("False positive rate (Fp)") +
                  # ylab("Sensitivity (Se)")
                  print(plot_6)
                  plot_list[[6]] <- plot_6
                  ##
                  ##
                  ## ---- Plot 7 (plots of Se and Sp vs. threshold for each test - base plot)
                  ##
                  # Normalize thresholds for each test to 0-1 range
                  tibble_NMA_normalized <- tibble_NMA %>%
                    group_by(test_name) %>%
                    mutate(
                      # Find actual threshold range for each test
                      min_thr = min(threshold[!is.na(Se_median)]),
                      max_thr = max(threshold[!is.na(Se_median)]),
                      # Normalize to 0-1
                      threshold_norm = (threshold - min_thr) / (max_thr - min_thr),
                      # Create labels for key points
                      threshold_label = ifelse(
                        threshold %in% c(min_thr, max_thr, round((min_thr + max_thr)/2)),
                        as.character(threshold),
                        ""
                      )
                    ) %>%
                    ungroup()
                  ##
                  plot_Se_vs_threshold_base <-   ggplot( tibble_NMA_normalized, 
                                                    mapping = aes(x = threshold_norm, 
                                                                  y = Se_median, 
                                                                  colour = test_name)) + 
                    geom_line(linewidth = linewidth) + 
                    geom_point(size = point_size) + 
                    theme_bw(base_size = base_size) + 
                    xlab("Test threshold") + 
                    ylab("Sensitivity (Se)") + 
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    theme(
                      axis.text.x = element_blank(),  # Remove tick labels
                      axis.ticks.x = element_blank()  # Remove tick marks
                    ) + 
                    labs(colour = "Test") + 
                    theme(legend.position = "bottom")
                  ##
                  plot_Se_vs_threshold_base
                  ##
                  plot_Se_vs_threshold_base_w_CrI <- plot_Se_vs_threshold_base +   
                                                 geom_ribbon(aes( x = threshold_norm, 
                                                                  ymin = Se_lower, 
                                                                  ymax = Se_upper,
                                                                  fill = test_name),
                                                             alpha = conf_region_alpha) 
                  ##
                  plot_Se_vs_threshold_panel <- plot_Se_vs_threshold_base_w_CrI + 
                    facet_wrap(~ test_name, 
                               scales = "free", 
                               ncol = 1) + 
                    geom_ribbon(aes( x = threshold_norm, 
                                     ymin = Se_pred_lower, 
                                     ymax = Se_pred_upper,
                                     fill = test_name),
                                alpha = pred_region_alpha) + 
                    theme(legend.position = "none")
                  plot_Se_vs_threshold_panel
                  ##
                  ##
                  plot_Sp_vs_threshold_base <-   ggplot( tibble_NMA_normalized, 
                                                         mapping = aes(x = threshold_norm, 
                                                                       y = Sp_median, 
                                                                       colour = test_name)) + 
                    geom_line(linewidth = linewidth) + 
                    geom_point(size = point_size) + 
                    theme_bw(base_size = base_size) + 
                    xlab("Test threshold") + 
                    ylab("Specificity (Sp)") + 
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    theme(
                      axis.text.x = element_blank(),  # Remove tick labels
                      axis.ticks.x = element_blank()  # Remove tick marks
                    ) + 
                    labs(colour = "Test") + 
                    theme(legend.position = "bottom")
                  ##
                  plot_Sp_vs_threshold_base
                  ##
                  plot_Sp_vs_threshold_base_w_CrI <- plot_Sp_vs_threshold_base +   
                                                  geom_ribbon(aes( x = threshold_norm, 
                                                                   ymin = Sp_lower, 
                                                                   ymax = Sp_upper,
                                                                   fill = test_name),
                                                                   alpha = conf_region_alpha)
                  plot_Sp_vs_threshold_base_w_CrI
                  ##
                  plot_Sp_vs_threshold_panel <- plot_Sp_vs_threshold_base_w_CrI + 
                    facet_wrap(~ test_name, 
                               scales = "free", 
                               ncol = 1) + 
                    geom_ribbon(aes( x = threshold_norm, 
                                     ymin = Sp_pred_lower, 
                                     ymax = Sp_pred_upper,
                                     fill = test_name),
                                alpha = pred_region_alpha) + 
                    theme(legend.position = "none")
                  plot_Sp_vs_threshold_panel
                  ##
                  plot_7 <- plot_Se_vs_threshold_base + 
                    theme(legend.position = "none") + 
                    xlab(" ") +  
                    plot_Sp_vs_threshold_base + 
                    patchwork::plot_layout(ncol = 1)
                  ##
                  print(plot_7)
                  plot_list[[7]] <- plot_7
                  ##
                  plot_8 <- plot_Se_vs_threshold_base_w_CrI + 
                    theme(legend.position = "none") + 
                    xlab(" ") + 
                    plot_Sp_vs_threshold_base_w_CrI + 
                    patchwork::plot_layout(ncol = 1)
                  ##
                  print(plot_8)
                  plot_list[[8]] <- plot_8
                  ##
                  plot_9 <- plot_Se_vs_threshold_panel + plot_Sp_vs_threshold_panel
                  print(plot_9)
                  plot_list[[9]] <- plot_9
                  ##
                
                  return(list(  plot_list = plot_list,
                                ##
                                tibble_NMA = tibble_NMA,
                                ##
                                polygon_Conf_list   = polygon_Conf_list,
                                polygon_Conf_tibble = polygon_Conf_tibble,
                                polygon_Pred_list   = polygon_Pred_list,
                                polygon_Pred_tibble = polygon_Pred_tibble))
  
}






 



#' plot_se_differences
#' @export
plot_se_differences <- function(filtered_tibble, 
                                show_top_n = 30,
                                order_by = "magnitude",
                                plot_title = "Sensitivity Differences Between Tests/Thresholds"
) {
  
        # Order by absolute difference magnitude or significance
        if (order_by == "magnitude") {
          plot_data <- filtered_tibble %>%
            mutate(abs_diff = abs(diff_Se_median)) %>%
            arrange(desc(abs_diff)) %>%
            slice_head(n = show_top_n)
        } else if (order_by == "significance") {
          # Order by significance (distance from 0 relative to SE)
          plot_data <- filtered_tibble %>%
            mutate(z_score = abs(diff_Se_median) / diff_Se_sd) %>%
            arrange(desc(z_score)) %>%
            slice_head(n = show_top_n)
        } else { 
          plot_data <- filtered_tibble %>%
            # arrange(test1, threshold1, test2, threshold2) %>%
            mutate(
              # Create factor with levels in this order
              comparison_label_ordered = factor(comparison_label, 
                                                levels = unique(comparison_label))
            )
        }
        
        ggplot(plot_data, aes(y = comparison_label, 
                              x = diff_Se_median)) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          geom_errorbarh(aes(xmin = diff_Se_lower, 
                             xmax = diff_Se_upper), 
                         height = 0.2, alpha = 0.5) +
          geom_point(size = 3) +
          facet_wrap(~ comparison_type, scales = "free_y", ncol = 1) +
          labs(x = "Difference in Sensitivity", 
               y = "",
               title = plot_title) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 16))
  
}





#' plot_se_differences
#' @export
plot_sp_differences <- function(filtered_tibble, 
                                show_top_n = 30,
                                order_by = "magnitude",
                                plot_title = "Specificity Differences Between Tests/Thresholds"
) {
        
        # Order by absolute difference magnitude or significance
        if (order_by == "magnitude") {
          plot_data <- filtered_tibble %>%
            mutate(abs_diff = abs(diff_Sp_median)) %>%
            arrange(desc(abs_diff)) %>%
            slice_head(n = show_top_n)
        } else if (order_by == "significance") {
          # Order by significance (distance from 0 relative to SE)
          plot_data <- filtered_tibble %>%
            mutate(z_score = abs(diff_Sp_median) / diff_Sp_sd) %>%
            arrange(desc(z_score)) %>%
            slice_head(n = show_top_n)
        } else { 
          plot_data <- filtered_tibble %>%
            # arrange(test1, threshold1, test2, threshold2) %>%
            mutate(
              # Create factor with levels in this order
              comparison_label_ordered = factor(comparison_label, 
                                                levels = unique(comparison_label))
            )
        }
        
        ggplot(plot_data, aes(y = comparison_label, 
                              x = diff_Sp_median)) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          geom_errorbarh(aes(xmin = diff_Sp_lower, 
                             xmax = diff_Sp_upper), 
                         height = 0.2, alpha = 0.5) +
          geom_point(size = 3) +
          facet_wrap(~ comparison_type, scales = "free_y", ncol = 1) +
          labs(x = "Difference in Specificity", 
               y = "",
               title = plot_title) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 16))
  
}





#' Compute analysis for multiple baseline scenarios
#' @export
R_fn_compute_MR_multiple_scenarios <- function(
                                                model_prep_obj,
                                                model_summary_and_trace_obj,
                                                baseline_case_nd_list,
                                                baseline_case_d_list,
                                                scenario_names = NULL,
                                                use_probit_link = TRUE,
                                                test_names = NULL,
                                                compute_comparisons = TRUE,
                                                save_results = TRUE,
                                                output_dir = "MR_results"
) {
        
        require(dplyr)
        require(purrr)
        
        n_scenarios <- length(baseline_case_nd_list)
        
        # Use provided names or create default ones
        if (is.null(scenario_names)) {
          scenario_names <- names(baseline_case_nd_list)
          if (is.null(scenario_names)) {
            scenario_names <- paste0("Scenario_", 1:n_scenarios)
          }
        }
        
        # Create output directory if saving
        if (save_results && !dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        # Process each scenario
        scenario_results <- list()
        
        cat("Processing", n_scenarios, "meta-regression scenarios...\n")
        
        for (i in 1:n_scenarios) {
          
              cat("\nScenario", i, ":", scenario_names[i], "\n")
              
              # Compute results for this scenario
              results <- R_fn_compute_MR_complete_new_baseline_analysis_v1(
                debugging = FALSE,
                model_prep_obj = model_prep_obj,
                model_summary_and_trace_obj = model_summary_and_trace_obj,
                baseline_case_nd = baseline_case_nd_list[[i]],
                baseline_case_d = baseline_case_d_list[[i]],
                use_probit_link = use_probit_link,
                test_names = test_names,
                compute_comparisons = compute_comparisons
              )
              
              # Add scenario info
              results$scenario_name <- scenario_names[i]
              results$scenario_index <- i
              
              scenario_results[[scenario_names[i]]] <- results
              
              # Save individual results if requested
              if (save_results) {
                saveRDS(results, file.path(output_dir, paste0("scenario_", i, "_", 
                                                              gsub("[^A-Za-z0-9]", "_", scenario_names[i]), 
                                                              ".rds")))
              }
          
        }
        
        return(scenario_results)
  
}



  




#' Create sROC plots for meta-regression scenarios
#' @export
R_fn_plot_MR_scenarios_sROC <- function(
                                        scenario_results,
                                        test_names = NULL,
                                        plot_type = c("tests_by_scenario", "scenarios_by_test", "individual"),
                                        show_CI = FALSE,
                                        show_PI = FALSE,
                                        base_size = 14,
                                        point_size = 2,
                                        line_size = 0.8,
                                        save_plots = TRUE,
                                        output_dir = "MR_plots",
                                        plot_width = 12,
                                        plot_height = 10
) {
  
        require(ggplot2)
        require(dplyr)
        require(purrr)
        require(patchwork)
        
        plot_type <- match.arg(plot_type, several.ok = TRUE)
        
        # Get n_thr from model_prep_obj
        n_thr <- model_prep_obj$internal_obj$outs_data$stan_data_list$n_thr
        n_index_tests <- model_prep_obj$internal_obj$outs_data$stan_data_list$n_index_tests
        
        # Extract test names if not provided
        if (is.null(test_names)) {
          if (!is.null(model_prep_obj$stan_data_list$test_names)) {
            test_names <- model_prep_obj$stan_data_list$test_names
          } else {
            test_names <- paste0("Test", 1:n_index_tests)
          }
        }
        
        n_scenarios <- length(scenario_results)
        n_tests <- length(test_names)
        scenario_names <- names(scenario_results)
  
        
        # Create output directory
        if (save_plots && !dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        # Prepare combined data for plotting  
        combined_data <- map_df( 1:length(scenario_results), 
                                 function(i) {
                
                res <- scenario_results[[i]]
                scenario_name <- names(scenario_results)[i]
                
                # Create tibble data structure
                tibble_data <- list(
                  Se_summaries = res$Se_summaries,
                  Sp_summaries = res$Sp_summaries,
                  Fp_summaries = res$Fp_summaries,
                  Se_pred_summaries = res$Se_pred_summaries,
                  Sp_pred_summaries = res$Sp_pred_summaries,
                  Fp_pred_summaries = res$Fp_pred_summaries
                )
                
                # Use existing plotting function to get formatted data
                plot_results <- R_fn_sROC_plot_NMA(
                  stan_model_file_name = scenario_name,
                  tibbles_provided = TRUE,
                  tibbles = tibble_data,
                  n_index_tests = n_index_tests,
                  n_thr = n_thr,  # Use n_thr from model_prep_obj
                  test_names = test_names
                )
                
                plot_results$tibble_NMA %>%
                  mutate(
                    scenario = scenario_name,
                    scenario_index = i
                  )
          
        })
        
        
        plot_list <- list()
        
        # 1. Tests by Scenario: 2x2 panels (one per test), 9 scenarios per panel
        if ("tests_by_scenario" %in% plot_type) {
          
          # Create color palette for scenarios
          scenario_colors <- scales::hue_pal()(n_scenarios)
          
          # Create plot for each test
          test_plots <- map(1:n_tests, function(t) {
            
            test_data <- combined_data %>%
              filter(test == t)
            
            p <- ggplot(test_data, 
                        aes(x = Fp_median, y = Se_median, 
                            color = scenario,
                            group = scenario)) +
              geom_line(size = line_size) +
              geom_point(size = point_size) +
              scale_color_manual(values = scenario_colors,
                                 name = "Scenario") +
              scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1)) +
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1)) +
              labs(x = "False positive rate",
                   y = "Sensitivity",
                   title = test_names[t]) +
              theme_bw(base_size = base_size) +
              theme(legend.position = "bottom",
                    legend.text = element_text(size = 8))
            
            p
          })
          
          # Combine into 2x2 layout
          combined_plot <- wrap_plots(test_plots, ncol = 2) + 
            plot_layout(guides = "collect") & 
            theme(legend.position = "bottom")
          
          plot_list[["tests_by_scenario"]] <- combined_plot
          
          if (save_plots) {
            ggsave(file.path(output_dir, "sROC_tests_by_scenario.pdf"),
                   combined_plot, 
                   width = plot_width, 
                   height = plot_height)
          }
        }
        
        # 2. Scenarios by Test: 3x3 panels (one per scenario), 4 tests per panel
        if ("scenarios_by_test" %in% plot_type) {
          
          # Create color palette for tests
          test_colors <- scales::hue_pal()(n_tests)
          
          # Create plot for each scenario
          scenario_plots <- map(1:n_scenarios, function(s) {
            
            scenario_data <- combined_data %>%
              filter(scenario_index == s)
            
            p <- ggplot(scenario_data, 
                        aes(x = Fp_median, y = Se_median, 
                            color = test_name,
                            group = test_name)) +
              geom_line(size = line_size) +
              geom_point(size = point_size) +
              scale_color_manual(values = test_colors,
                                 name = "Test") +
              scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1)) +
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1)) +
              labs(x = "False positive rate",
                   y = "Sensitivity",
                   title = scenario_names[s]) +
              theme_bw(base_size = base_size) +
              theme(legend.position = "none",
                    plot.title = element_text(size = 10))
            
            p
          })
          
          # Combine into 3x3 layout
          combined_plot <- wrap_plots(scenario_plots, ncol = 3) + 
            plot_layout(guides = "collect") & 
            theme(legend.position = "bottom")
          
          plot_list[["scenarios_by_test"]] <- combined_plot
          
          if (save_plots) {
            ggsave(file.path(output_dir, "sROC_scenarios_by_test.pdf"),
                   combined_plot, 
                   width = plot_width * 1.2, 
                   height = plot_height * 1.2)
          }
        }
        
        # 3. Individual plots with CI/PI
        if ("individual" %in% plot_type) {
          
          individual_plots <- list()
          
          for (s in 1:n_scenarios) {
            for (t in 1:n_tests) {
              
              # Get data for this scenario
              scenario_data <- scenario_results[[s]]
              
              # Create individual plot with CI/PI using existing function
              plot_results <- R_fn_sROC_plot_NMA(
                stan_model_file_name = paste0(scenario_names[s], " - ", test_names[t]),
                tibbles_provided = TRUE,
                tibbles = list(
                  Se_summaries = scenario_data$Se_summaries,
                  Sp_summaries = scenario_data$Sp_summaries,
                  Fp_summaries = scenario_data$Fp_summaries,
                  Se_pred_summaries = scenario_data$Se_pred_summaries,
                  Sp_pred_summaries = scenario_data$Sp_pred_summaries,
                  Fp_pred_summaries = scenario_data$Fp_pred_summaries
                ),
                n_index_tests = n_tests,
                n_thr = scenario_data$n_thr,
                test_names = test_names,
                base_size = base_size,
                point_size = point_size,
                linewidth = line_size
              )
              
              # Extract the plot with CI and PI (plot 5)
              individual_plot <- plot_results$plot_list[[5]]
              
              individual_plots[[paste0("scenario_", s, "_test_", t)]] <- individual_plot
              
              if (save_plots) {
                ggsave(file.path(output_dir, 
                                 paste0("sROC_individual_S", s, "_T", t, "_",
                                        gsub("[^A-Za-z0-9]", "_", scenario_names[s]), "_",
                                        test_names[t], ".pdf")),
                       individual_plot,
                       width = plot_width * 0.8, 
                       height = plot_height * 0.8)
              }
            }
          }
          
          plot_list[["individual"]] <- individual_plots
        }
        
        return(plot_list)
  
}




 







#' Create comparison heatmap showing effect of covariates
#' @export
R_fn_plot_MR_comparison_heatmap <- function(
                                            all_comparisons_data,
                                            comparison_type = c("Se", "Sp"),
                                            n_comparisons = 20,
                                            base_size = 14,
                                            save_plot = TRUE,
                                            output_dir = "MR_plots",
                                            plot_width = 12,
                                            plot_height = 10
) {
  
        require(ggplot2)
        require(dplyr)
        require(tidyr)
        require(viridis)
        
        comparison_type <- match.arg(comparison_type)
        
        # Select top comparisons by absolute median difference
        if (comparison_type == "Se") {
          top_comparisons <- all_comparisons_data %>%
            group_by(test1_name, test2_name, threshold1, threshold2) %>%
            summarise(
              mean_abs_diff = mean(abs(diff_Se_median), na.rm = TRUE),
              .groups = "drop"
            ) %>%
            arrange(desc(mean_abs_diff)) %>%
            slice_head(n = n_comparisons) %>%
            mutate(
              comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
                                        test2_name, " (", threshold2, ")")
            ) %>%
            pull(comparison_label)
          
          # Prepare data for heatmap
          heatmap_data <- all_comparisons_data %>%
            mutate(
              comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
                                        test2_name, " (", threshold2, ")")
            ) %>%
            filter(comparison_label %in% top_comparisons) %>%
            select(comparison_label, scenario, diff_Se_median) %>%
            pivot_wider(names_from = scenario, values_from = diff_Se_median)
          
          # Convert to matrix for heatmap
          heatmap_matrix <- as.matrix(heatmap_data[, -1])
          rownames(heatmap_matrix) <- heatmap_data$comparison_label
          
        } else {
          top_comparisons <- all_comparisons_data %>%
            group_by(test1_name, test2_name, threshold1, threshold2) %>%
            summarise(
              mean_abs_diff = mean(abs(diff_Sp_median), na.rm = TRUE),
              .groups = "drop"
            ) %>%
            arrange(desc(mean_abs_diff)) %>%
            slice_head(n = n_comparisons) %>%
            mutate(
              comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
                                        test2_name, " (", threshold2, ")")
            ) %>%
            pull(comparison_label)
          
          heatmap_data <- all_comparisons_data %>%
            mutate(
              comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
                                        test2_name, " (", threshold2, ")")
            ) %>%
            filter(comparison_label %in% top_comparisons) %>%
            select(comparison_label, scenario, diff_Sp_median) %>%
            pivot_wider(names_from = scenario, values_from = diff_Sp_median)
          
          heatmap_matrix <- as.matrix(heatmap_data[, -1])
          rownames(heatmap_matrix) <- heatmap_data$comparison_label
        }
        
        # Create heatmap using ggplot2
        heatmap_long <- heatmap_data %>%
          pivot_longer(-comparison_label, names_to = "scenario", values_to = "difference")
        
        p <- ggplot(heatmap_long, 
                    aes(x = scenario, y = comparison_label, fill = difference)) +
          geom_tile() +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                               midpoint = 0,
                               name = paste0("Difference in ", 
                                             ifelse(comparison_type == "Se", "Sensitivity", "Specificity"))) +
          theme_bw(base_size = base_size) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 9),
            axis.title = element_blank(),
            legend.position = "right"
          ) +
          labs(title = paste0("Effect of Meta-regression Covariates on ",
                              ifelse(comparison_type == "Se", "Sensitivity", "Specificity"),
                              " Differences"))
        
        if (save_plot) {
          if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
          }
          
          ggsave(
            file.path(output_dir, paste0("heatmap_", comparison_type, "_comparisons.pdf")),
            p,
            width = plot_width,
            height = plot_height
          )
        }
        
        return(p)
  
}






 




#' Extract sROC data for a single scenario
#' @export
extract_sroc_data_single_scenario <- function(
                                              Se_summaries,
                                              Sp_summaries,
                                              Fp_summaries,
                                              Se_pred_summaries,
                                              Sp_pred_summaries,
                                              Fp_pred_summaries,
                                              scenario_name = "Scenario",
                                              scenario_index = 1
) {
  
        require(dplyr)
        
        # Create a tibble from the summaries
        n_rows <- nrow(Se_summaries)
        
        tibble(
          scenario = scenario_name,
          scenario_index = scenario_index,
          variable = Se_summaries$variable,
          # Parse test and threshold from variable name
          test = as.integer(gsub(".*\\[(\\d+),.*", "\\1", variable)),
          threshold = as.integer(gsub(".*,(\\d+)\\]", "\\1", variable)),
          # Se/Sp/Fp values
          Se_median = Se_summaries$`50%`,
          Se_lower = Se_summaries$`2.5%`,
          Se_upper = Se_summaries$`97.5%`,
          Sp_median = Sp_summaries$`50%`,
          Sp_lower = Sp_summaries$`2.5%`,
          Sp_upper = Sp_summaries$`97.5%`,
          Fp_median = Fp_summaries$`50%`,
          Fp_lower = Fp_summaries$`2.5%`,
          Fp_upper = Fp_summaries$`97.5%`,
          # Predictive intervals
          Se_pred_lower = Se_pred_summaries$`2.5%`,
          Se_pred_upper = Se_pred_summaries$`97.5%`,
          Sp_pred_lower = Sp_pred_summaries$`2.5%`,
          Sp_pred_upper = Sp_pred_summaries$`97.5%`,
          Fp_pred_lower = Fp_pred_summaries$`2.5%`,
          Fp_pred_upper = Fp_pred_summaries$`97.5%`
        ) %>%
          filter(Se_median > -0.5)  # Remove missing values
  
}









#' Combine sROC data from multiple scenarios
#' @export
combine_sroc_data_scenarios <- function(scenario_results) {
  
        require(dplyr)
        require(purrr)
        
        # Extract data from each scenario
        combined_data <- map_df(1:length(scenario_results), function(i) {
          res <- scenario_results[[i]]
          scenario_name <- names(scenario_results)[i]
          
          extract_sroc_data_single_scenario(
            Se_summaries = res$Se_summaries,
            Sp_summaries = res$Sp_summaries,
            Fp_summaries = res$Fp_summaries,
            Se_pred_summaries = res$Se_pred_summaries,
            Sp_pred_summaries = res$Sp_pred_summaries,
            Fp_pred_summaries = res$Fp_pred_summaries,
            scenario_name = scenario_name,
            scenario_index = i
          )
        })
        
        return(combined_data)
  
}




#' Create basic sROC plot
#' @export
create_sroc_plot_base <- function(
                                  data,
                                  x_var = "Fp_median",
                                  y_var = "Se_median",
                                  color_by = NULL,
                                  facet_by = NULL,
                                  title = NULL,
                                  subtitle = NULL,
                                  base_size = 14,
                                  point_size = 2,
                                  line_size = 0.8
) {
  
        require(ggplot2)
        
        p <- ggplot(data, aes_string(x = x_var, y = y_var))
        
        # Add color aesthetic if specified
        if (!is.null(color_by)) {
          p <- p + aes_string(color = color_by, group = color_by)
        }
        
        # Add lines and points
        p <- p +
          geom_line(size = line_size) +
          geom_point(size = point_size) +
          scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          labs(x = "False positive rate (1 - Specificity)",
               y = "Sensitivity") +
          theme_bw(base_size = base_size)
        
        # Add faceting if specified
        if (!is.null(facet_by)) {
          p <- p + facet_wrap(as.formula(paste("~", facet_by)))
        }
        
        # Add title and subtitle
        if (!is.null(title)) {
          p <- p + ggtitle(title, subtitle = subtitle)
        }
        
        return(p)
        
}





#' Create sROC plots with different layouts
#' @export
create_sroc_plots_flexible <- function(
                                      scenario_results,
                                      test_names = NULL,
                                      layout = c("separate", "by_test", "by_scenario", "grid"),
                                      base_size = 14,
                                      point_size = 2,
                                      line_size = 0.8
) {
  
        require(ggplot2)
        require(dplyr)
        require(patchwork)
        
        # Combine all scenario data
        combined_data <- combine_sroc_data_scenarios(scenario_results)
        
        # Add test names if provided
        if (!is.null(test_names)) {
          combined_data <- combined_data %>%
            mutate(test_name = factor(test_names[test], levels = test_names))
        } else {
          combined_data <- combined_data %>%
            mutate(test_name = factor(paste0("Test ", test)))
        }
        
        # Get unique values
        n_scenarios <- length(unique(combined_data$scenario))
        n_tests <- length(unique(combined_data$test))
        scenario_names <- unique(combined_data$scenario)
        
        plots <- list()
        
        # 1. Separate plots for each test-scenario combination
        if ("separate" %in% layout) {
          
          for (s in scenario_names) {
            for (t in unique(combined_data$test)) {
              
              plot_data <- combined_data %>%
                filter(scenario == s, test == t)
              
              test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
              
              p <- create_sroc_plot_base(
                data = plot_data,
                title = paste0(test_name, " - ", s),
                base_size = base_size,
                point_size = point_size,
                line_size = line_size
              )
              
              plots[[paste0(s, "_", test_name)]] <- p
            }
          }
        }
        
        # 2. One plot per test, all scenarios overlaid
        if ("by_test" %in% layout) {
          
          for (t in unique(combined_data$test)) {
            
            plot_data <- combined_data %>%
              filter(test == t)
            
            test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
            
            p <- create_sroc_plot_base(
              data = plot_data,
              color_by = "scenario",
              title = paste0("sROC curves for ", test_name),
              subtitle = "Colored by meta-regression scenario",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size
            )
            
            plots[[paste0("by_test_", test_name)]] <- p
          }
        }
        
        # 3. One plot per scenario, all tests overlaid
        if ("by_scenario" %in% layout) {
          
          for (s in scenario_names) {
            
            plot_data <- combined_data %>%
              filter(scenario == s)
            
            p <- create_sroc_plot_base(
              data = plot_data,
              color_by = "test_name",
              title = s,
              subtitle = "Colored by test",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size
            )
            
            plots[[paste0("by_scenario_", s)]] <- p
          }
        }
        
        # 4. Grid layout
        if ("grid" %in% layout) {
          
          # Grid with tests as panels, scenarios as colors
          p1 <- create_sroc_plot_base(
            data = combined_data,
            color_by = "scenario",
            facet_by = "test_name",
            title = "sROC curves by test",
            subtitle = "Scenarios shown in different colors",
            base_size = base_size,
            point_size = point_size,
            line_size = line_size
          )
          
          plots[["grid_by_test"]] <- p1
          
          # Grid with scenarios as panels, tests as colors
          p2 <- create_sroc_plot_base(
            data = combined_data,
            color_by = "test_name",
            facet_by = "scenario",
            title = "sROC curves by scenario",
            subtitle = "Tests shown in different colors",
            base_size = base_size,
            point_size = point_size,
            line_size = line_size
          )
          
          plots[["grid_by_scenario"]] <- p2
        }
        
        return(plots)
      
}


#' create_prediction_region_smooth
#' @keywords internal
#' @export
create_prediction_region_smooth <- function(df, 
                                            model_name, 
                                            n_points = 100) {
        # 
        # pred_data <- df %>% 
        #   filter(Model == model_name)
        
        pred_data <- df
        
        # Create smooth curves for boundaries
        fp_seq <- seq(0, 1, length.out = n_points)
        
        # Interpolate upper and lower bounds
        upper_se <- approx(pred_data$Fp_pred_lower, pred_data$Se_pred_upper, 
                           xout = fp_seq, rule = 2)$y
        lower_se <- approx(pred_data$Fp_pred_upper, pred_data$Se_pred_lower, 
                           xout = fp_seq, rule = 2)$y
        
        # Create polygon points
        polygon_points <- data.frame(
          x = c(fp_seq, rev(fp_seq)),
          y = c(upper_se, rev(lower_se))#,
          # Model = model_name
        )
        
        return(polygon_points)
  
}




  


#' create_confidence_polygon - Create credible interval polygon function
#' @keywords internal
#' @export
create_confidence_polygon <- function(df, 
                                      model_name = NULL) {
  
        # Create upper curve (high sensitivity, low false positive)
        upper <- df %>% 
          # dplyr::filter(Model == model_name) %>%
          dplyr::select(x = Fp_lower,
                        y = Se_upper) %>%
          dplyr::arrange(x)
        
        # Create lower curve (low sensitivity, high false positive)
        lower <- df %>% 
          # dplyr::filter(Model == model_name) %>%
          dplyr::select( x = Fp_upper, 
                         y = Se_lower) %>%
          dplyr::arrange(desc(x))
        
        # Combine into a closed polygon
        out <- tibble(bind_rows(upper, lower))
        n_rows <- nrow(out)
        out
        # out <- out %>%
        #   dplyr::mutate(Model = rep(model_name, n_rows))
  
}

 


#' create_prediction_polygon -  Create prediction interval polygon function
#' @keywords internal
#' @export
create_prediction_polygon <- function(df, 
                                      model_name = NULL) {
  
  # Create upper curve (high sensitivity, low false positive)
  upper <- df %>% 
    # dplyr::filter(Model == model_name) %>%
    dplyr::select( x = Fp_pred_lower, 
                   y = Se_pred_upper) %>% 
    dplyr::arrange(x)
  
  # Create lower curve (low sensitivity, high false positive)
  lower <- df %>% 
    # dplyr::filter(Model == model_name) %>%
    dplyr::select( x = Fp_pred_upper, 
                   y = Se_pred_lower) %>%
    dplyr::arrange(desc(x))
  
  # Combine into a closed polygon
  out <- tibble(bind_rows(upper, lower))
  n_rows <- nrow(out)
  out
  # out <- out %>% dplyr::mutate(Model = rep(model_name, n_rows))
  
}




  
add_regions_to_sroc <- function(
    plot,
    tibble_NMA,
    show_conf = TRUE,
    show_pred = FALSE,
    conf_region_alpha = 0.15,
    pred_region_alpha = 0.30,
    conf_region_colour = "blue",
    pred_region_colour = "green"
) {

        require(ggplot2)
  
        polygon_Conf_list <- polygon_Pred_list <- list()
        
        for (t in 1:n_index_tests) {
          
          tibble_test_t <- dplyr::filter(tibble_NMA, test == t)
          
          # Confidence region
          tibble_Conf_polygon_test_t <- create_confidence_polygon(tibble_test_t)
          n_rows = nrow(tibble_Conf_polygon_test_t)
          tibble_Conf_polygon_test_t <- tibble_Conf_polygon_test_t %>% 
            dplyr::mutate(test = rep(t, n_rows),
                          test_name = factor(rep(test_names[t], n_rows), levels = test_names))
          polygon_Conf_list[[t]] <- tibble_Conf_polygon_test_t
          
          # Prediction region
          tibble_Pred_polygon_test_t <- create_prediction_region_smooth(tibble_test_t)
          n_rows = nrow(tibble_Pred_polygon_test_t)
          tibble_Pred_polygon_test_t <- tibble_Pred_polygon_test_t %>% 
            dplyr::mutate(test = rep(t, n_rows), 
                          test_name = factor(rep(test_names[t], n_rows), levels = test_names))
          polygon_Pred_list[[t]] <- tibble_Pred_polygon_test_t
          
        }
        
        polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% dplyr::filter(!(is.na(x)))
        polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% dplyr::filter(!(is.na(x)))

        if (show_conf) {
          plot <- plot +
            geom_polygon(data = polygon_Conf_tibble,
                         aes(x = x, y = y),
                         fill = conf_region_colour,
                         alpha = conf_region_alpha)  
        }

        if (show_pred) {
          plot <- plot +
            geom_polygon(data = polygon_Pred_tibble,
                         aes(x = x, y = y),
                         fill = pred_region_colour,
                         alpha = pred_region_alpha)  
        }

        return(plot)

}








#' Create forest plots for specific comparisons across scenarios (GENERAL VERSION)
#' @export
R_fn_plot_forest_comparisons_general <- function(
                                                all_comparisons_data,
                                                selected_comparisons = NULL,
                                                comparison_type = c("Se", "Sp"),
                                                scenario_var = "scenario",  # Column name containing scenario labels
                                                group_var = NULL,           # Optional: column to use for coloring
                                                shape_var = NULL,           # Optional: column to use for shapes
                                                n_auto_select = 6,          # Number of comparisons to auto-select
                                                base_size = 14,
                                                point_size = 3,
                                                error_bar_width = 0.2,
                                                error_bar_alpha = 0.7,
                                                facet_ncol = 2,
                                                facet_scales = "free_x",
                                                plot_title = NULL,
                                                x_label = NULL,
                                                y_label = "Scenario",
                                                save_plot = FALSE,
                                                output_path = NULL
) {
  
        require(ggplot2)
        require(dplyr)
        
        comparison_type <- match.arg(comparison_type)
        
        # Set default labels based on comparison type
        if (is.null(x_label)) {
          x_label <- paste0("Difference in ", 
                            ifelse(comparison_type == "Se", "Sensitivity", "Specificity"))
        }
        
        if (is.null(plot_title)) {
          plot_title <- paste0(ifelse(comparison_type == "Se", "Sensitivity", "Specificity"),
                               " Differences Across Scenarios")
        }
        
        # Create comparison label if it doesn't exist
        if (!"comparison_label" %in% names(all_comparisons_data)) {
          all_comparisons_data <- all_comparisons_data %>%
            mutate(comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
                                             test2_name, " (", threshold2, ")"))
        }
        
        # Auto-select comparisons if not provided
        if (is.null(selected_comparisons)) {
          # Find comparisons with largest variance across scenarios
          variance_summary <- all_comparisons_data %>%
            group_by(comparison_label) %>%
            summarise(
              var_diff = var(get(paste0("diff_", comparison_type, "_median")), na.rm = TRUE),
              mean_abs_diff = mean(abs(get(paste0("diff_", comparison_type, "_median"))), na.rm = TRUE),
              .groups = "drop"
            ) %>%
            filter(mean_abs_diff > 0.05) %>%  # At least 5% difference
            arrange(desc(var_diff)) %>%
            slice_head(n = n_auto_select)
          
          selected_comparisons <- variance_summary$comparison_label
        }
        
        # Filter to selected comparisons
        plot_data <- all_comparisons_data %>%
          filter(comparison_label %in% selected_comparisons)
        
        # Get the correct column names
        x_col <- paste0("diff_", comparison_type, "_median")
        xmin_col <- paste0("diff_", comparison_type, "_lower")
        xmax_col <- paste0("diff_", comparison_type, "_upper")
        
        # Create base plot
        p <- ggplot(plot_data, 
                    aes_string(y = scenario_var, x = x_col)) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          geom_errorbarh(aes_string(xmin = xmin_col, xmax = xmax_col), 
                         height = error_bar_width, 
                         alpha = error_bar_alpha)
        
        # Add points with optional color and shape mappings
        if (!is.null(group_var) && !is.null(shape_var)) {
          p <- p + geom_point(aes_string(color = group_var, shape = shape_var), 
                              size = point_size)
        } else if (!is.null(group_var)) {
          p <- p + geom_point(aes_string(color = group_var), 
                              size = point_size)
        } else if (!is.null(shape_var)) {
          p <- p + geom_point(aes_string(shape = shape_var), 
                              size = point_size)
        } else {
          p <- p + geom_point(size = point_size)
        }
        
        # Add faceting
        p <- p + 
          facet_wrap(~ comparison_label, ncol = facet_ncol, scales = facet_scales) +
          labs(x = x_label,
               y = y_label,
               title = plot_title) +
          theme_bw(base_size = base_size) +
          theme(
            axis.text.y = element_text(size = base_size * 0.8),
            legend.position = "bottom",
            strip.text = element_text(size = base_size * 0.8, face = "bold"),
            panel.spacing = unit(1, "lines")
          )
        
        # Save if requested
        if (save_plot && !is.null(output_path)) {
          ggsave(output_path, p, width = 10, height = 8)
        }
        
        return(p)
  
}








#' Extract Se/Sp vs threshold data for scenarios (GENERAL)
#' @export
extract_threshold_data_scenarios <- function(scenario_results) {
  
        require(dplyr)
        require(purrr)
        
        # Extract data from each scenario
        combined_data <- map_df(1:length(scenario_results), function(i) {
          res <- scenario_results[[i]]
          scenario_name <- names(scenario_results)[i]
          
          extract_sroc_data_single_scenario(
            Se_summaries = res$Se_summaries,
            Sp_summaries = res$Sp_summaries,
            Fp_summaries = res$Fp_summaries,
            Se_pred_summaries = res$Se_pred_summaries,
            Sp_pred_summaries = res$Sp_pred_summaries,
            Fp_pred_summaries = res$Fp_pred_summaries,
            scenario_name = scenario_name,
            scenario_index = i
          )
        })
        
        # Normalize thresholds for each test
        combined_data <- combined_data %>%
          group_by(test, scenario) %>%
          mutate(
            # Find actual threshold range for each test
            min_thr = min(threshold[!is.na(Se_median)]),
            max_thr = max(threshold[!is.na(Se_median)]),
            # Normalize to 0-1
            threshold_norm = (threshold - min_thr) / (max_thr - min_thr),
            # For plotting actual threshold values
            threshold_actual = threshold
          ) %>%
          ungroup()
        
        return(combined_data)
  
}





#' Create base threshold plot (GENERAL)
#' @export
create_threshold_plot_base <- function(
                                        data,
                                        y_var = "Se_median",
                                        x_var = "threshold_norm",
                                        color_by = NULL,
                                        linetype_by = NULL,
                                        shape_by = NULL,
                                        group_by = NULL,
                                        facet_by = NULL,
                                        title = NULL,
                                        subtitle = NULL,
                                        x_label = "Normalized threshold (0-1)",
                                        y_label = NULL,
                                        base_size = 14,
                                        point_size = 2,
                                        line_size = 0.8,
                                        hide_x_axis = TRUE
) {
        
        require(ggplot2)
        
        # Default y label based on variable
        if (is.null(y_label)) {
          y_label <- if (grepl("Se", y_var)) "Sensitivity" else "Specificity"
        }
        
        # Create base plot
        p <- ggplot(data, aes_string(x = x_var, y = y_var))
        
        # Add aesthetics
        aes_list <- list()
        if (!is.null(color_by)) aes_list$color <- color_by
        if (!is.null(linetype_by)) aes_list$linetype <- linetype_by
        if (!is.null(shape_by)) aes_list$shape <- shape_by
        if (!is.null(group_by)) aes_list$group <- group_by
        
        if (length(aes_list) > 0) {
          p <- p + aes_string(!!!aes_list)
        }
        
        # Add geoms
        p <- p +
          geom_line(size = line_size) +
          geom_point(size = point_size) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          labs(x = x_label, y = y_label) +
          theme_bw(base_size = base_size)
        
        # Add faceting if specified
        if (!is.null(facet_by)) {
          p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free_x")
        }
        
        # Hide x-axis if requested
        if (hide_x_axis) {
          p <- p + theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank())
        }
        
        # Add title and subtitle
        if (!is.null(title)) {
          p <- p + ggtitle(title, subtitle = subtitle)
        }
        
        return(p)
  
}






#' Create Se/Sp vs threshold plots with different layouts (GENERAL)
#' @export
create_threshold_plots_flexible <- function(
                                            scenario_results,
                                            test_names = NULL,
                                            plot_type = c("Se", "Sp", "both"),
                                            layout = c("by_test", "by_scenario", "grid"),
                                            base_size = 14,
                                            point_size = 2,
                                            line_size = 0.8,
                                            show_actual_thresholds = FALSE
) {
      
      require(ggplot2)
      require(dplyr)
      require(patchwork)
      
      plot_type <- match.arg(plot_type, several.ok = TRUE)
      layout <- match.arg(layout, several.ok = TRUE)
      
      # Get the data
      combined_data <- extract_threshold_data_scenarios(scenario_results)
      
      # Add test names if provided
      if (!is.null(test_names)) {
        combined_data <- combined_data %>%
          mutate(test_name = factor(test_names[test], levels = test_names))
      } else {
        combined_data <- combined_data %>%
          mutate(test_name = factor(paste0("Test ", test)))
      }
      
      # Determine x variable
      x_var <- ifelse(show_actual_thresholds, "threshold_actual", "threshold_norm")
      x_label <- ifelse(show_actual_thresholds, "Test threshold", "Normalized threshold (0-1)")
      
      plots <- list()
      
      # 1. By test (one plot per test, all scenarios shown)
      if ("by_test" %in% layout) {
        
        for (t in unique(combined_data$test)) {
          
          test_data <- combined_data %>% filter(test == t)
          test_name <- unique(test_data$test_name)
          
          if ("Se" %in% plot_type || "both" %in% plot_type) {
            p_se <- create_threshold_plot_base(
              data = test_data,
              y_var = "Se_median",
              x_var = x_var,
              color_by = "scenario",
              group_by = "scenario",
              title = paste0("Sensitivity vs Threshold - ", test_name),
              x_label = x_label,
              y_label = "Sensitivity",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size,
              hide_x_axis = !show_actual_thresholds
            )
            plots[[paste0("Se_by_test_", test_name)]] <- p_se
          }
          
          if ("Sp" %in% plot_type || "both" %in% plot_type) {
            p_sp <- create_threshold_plot_base(
              data = test_data,
              y_var = "Sp_median",
              x_var = x_var,
              color_by = "scenario",
              group_by = "scenario",
              title = paste0("Specificity vs Threshold - ", test_name),
              x_label = x_label,
              y_label = "Specificity",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size,
              hide_x_axis = !show_actual_thresholds
            )
            plots[[paste0("Sp_by_test_", test_name)]] <- p_sp
          }
        }
      }
      
      # 2. By scenario (one plot per scenario, all tests shown)
      if ("by_scenario" %in% layout) {
        
        for (s in unique(combined_data$scenario)) {
          
          scenario_data <- combined_data %>% filter(scenario == s)
          
          if ("Se" %in% plot_type || "both" %in% plot_type) {
            p_se <- create_threshold_plot_base(
              data = scenario_data,
              y_var = "Se_median",
              x_var = x_var,
              color_by = "test_name",
              group_by = "test_name",
              title = paste0("Sensitivity - ", s),
              x_label = x_label,
              y_label = "Sensitivity",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size,
              hide_x_axis = !show_actual_thresholds
            )
            plots[[paste0("Se_by_scenario_", s)]] <- p_se
          }
          
          if ("Sp" %in% plot_type || "both" %in% plot_type) {
            p_sp <- create_threshold_plot_base(
              data = scenario_data,
              y_var = "Sp_median",
              x_var = x_var,
              color_by = "test_name",
              group_by = "test_name",
              title = paste0("Specificity - ", s),
              x_label = x_label,
              y_label = "Specificity",
              base_size = base_size,
              point_size = point_size,
              line_size = line_size,
              hide_x_axis = !show_actual_thresholds
            )
            plots[[paste0("Sp_by_scenario_", s)]] <- p_sp
          }
        }
      }
      
      # 3. Grid layouts
      if ("grid" %in% layout) {
        
        if ("Se" %in% plot_type || "both" %in% plot_type) {
          p_se_grid <- create_threshold_plot_base(
            data = combined_data,
            y_var = "Se_median",
            x_var = x_var,
            color_by = "scenario",
            group_by = "scenario",
            facet_by = "test_name",
            title = "Sensitivity vs Threshold by Test",
            x_label = x_label,
            y_label = "Sensitivity",
            base_size = base_size,
            point_size = point_size,
            line_size = line_size,
            hide_x_axis = !show_actual_thresholds
          )
          plots[["Se_grid_by_test"]] <- p_se_grid
        }
        
        if ("Sp" %in% plot_type || "both" %in% plot_type) {
          p_sp_grid <- create_threshold_plot_base(
            data = combined_data,
            y_var = "Sp_median",
            x_var = x_var,
            color_by = "scenario",
            group_by = "scenario",
            facet_by = "test_name",
            title = "Specificity vs Threshold by Test",
            x_label = x_label,
            y_label = "Specificity",
            base_size = base_size,
            point_size = point_size,
            line_size = line_size,
            hide_x_axis = !show_actual_thresholds
          )
          plots[["Sp_grid_by_test"]] <- p_sp_grid
        }
      }
      
      return(plots)
  
}












