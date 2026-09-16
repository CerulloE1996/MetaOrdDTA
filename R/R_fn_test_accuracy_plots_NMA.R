







# #' R_fn_sROC_plot_NMA
# #' @keywords internal
# #' @export
# R_fn_sROC_plot_NMA <- function( stan_model_file_name,
#                                 stan_mod_samples = NULL,
#                                 ##
#                                 tibbles_provided = FALSE,
#                                 tibbles = NULL,
#                                 ##
#                                 # covariates,
#                                 # new_cov_data,
#                                 ##
#                                 n_index_tests,
#                                 n_thr,
#                                 ##
#                                 test_names = NULL,
#                                 relevant_thresholds = NULL,
#                                 ##
#                                 df_true = NULL,
#                                 ##
#                                 conf_region_colour = "blue", 
#                                 pred_region_colour = "green",
#                                 ##
#                                 conf_region_alpha = 0.30,
#                                 pred_region_alpha = 0.15,
#                                 ##
#                                 base_size = 16,
#                                 point_size = 3,
#                                 linewidth = 0.5
# ) {
# 
#                 require(ggplot2)
#                 require(dplyr)
#                 require(patchwork)
#                 ##
#                 model_name <- stan_model_file_name
#                 n_thr_max <- max(n_thr)
#                 ##
#                 # Set default test names if not provided
#                 if (is.null(test_names)) {
#                   test_names <- paste0("Test ", 1:n_index_tests)
#                 }
#                 ##
#                 Se_median_array <- Sp_median_array <- Fp_median_array <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 Se_mean_array   <- Sp_mean_array   <- Fp_mean_array   <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 ##
#                 Se_lower_array <- Sp_lower_array <- Fp_lower_array <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 Se_upper_array <- Sp_upper_array <- Fp_upper_array <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 ##
#                 Se_pred_lower_array <- Sp_pred_lower_array <- Fp_pred_lower_array <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 Se_pred_upper_array <- Sp_pred_upper_array <- Fp_pred_upper_array <- array(NA, dim = c(n_index_tests, n_thr_max))
#                 ##
#                 ## Extract summaries as tibbles:
#                 ##
#                 if (tibbles_provided == TRUE) {
#                   
#                     Se <- tibbles$Se_summaries
#                     Sp <- tibbles$Sp_summaries
#                     Fp <- tibbles$Fp_summaries
#                     ##
#                     Se_pred <- tibbles$Se_pred_summaries
#                     Sp_pred <- tibbles$Sp_pred_summaries
#                     Fp_pred <- tibbles$Fp_pred_summaries
#                     
#                 } else {
#                   
#                     Se <- stan_mod_samples$summary(c("Se_baseline"), mean, 
#                                                    quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     Sp <- stan_mod_samples$summary(c("Sp_baseline"), mean, 
#                                                    quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     Fp <- stan_mod_samples$summary(c("Fp_baseline"), mean, 
#                                                    quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     ##
#                     Se_pred <- stan_mod_samples$summary(c("Se_baseline_pred"), mean, 
#                                                         quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     Sp_pred <- stan_mod_samples$summary(c("Sp_baseline_pred"), mean, 
#                                                         quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     Fp_pred <- stan_mod_samples$summary(c("Fp_baseline_pred"), mean, 
#                                                         quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975)))
#                     
#                 }
#                 ##
#                 counter <- 1
#                 for (k in 1:n_thr_max) {
#                   for (t in 1:n_index_tests) {
#                     # test_vec[counter] <- t
#                     ## Posterior medians (of pooled estimates):
#                     Se_median_array[t, k] <- Se$`50%`[counter]
#                     Sp_median_array[t, k] <- Sp$`50%`[counter]
#                     Fp_median_array[t, k] <- Fp$`50%`[counter]
#                     ## Posterior means (of pooled estimates):
#                     Se_mean_array[t, k] <- Se$mean[counter]
#                     Sp_mean_array[t, k] <- Sp$mean[counter]
#                     Fp_mean_array[t, k] <- Fp$mean[counter]
#                     ## Posterior lower 95% (of pooled estimates):
#                     Se_lower_array[t, k] <- Se$`2.5%`[counter]
#                     Sp_lower_array[t, k] <- Sp$`2.5%`[counter]
#                     Fp_lower_array[t, k] <- Fp$`2.5%`[counter]
#                     ## Posterior upper 95% (of pooled estimates):
#                     Se_upper_array[t, k] <- Se$`97.5%`[counter]
#                     Sp_upper_array[t, k] <- Sp$`97.5%`[counter]
#                     Fp_upper_array[t, k] <- Fp$`97.5%`[counter]
#                     ## Posterior lower prediction 95% (of pooled estimates):
#                     Se_pred_lower_array[t, k] <- Se_pred$`2.5%`[counter]
#                     Sp_pred_lower_array[t, k] <- Sp_pred$`2.5%`[counter]
#                     Fp_pred_lower_array[t, k] <- Fp_pred$`2.5%`[counter]
#                     ## Posterior upper prediction 95% (of pooled estimates):
#                     Se_pred_upper_array[t, k] <- Se_pred$`97.5%`[counter]
#                     Sp_pred_upper_array[t, k] <- Sp_pred$`97.5%`[counter]
#                     Fp_pred_upper_array[t, k] <- Fp_pred$`97.5%`[counter]
#                     ##
#                     counter <- counter + 1
#                   }
#                 }
#                 ##
#                 test_vec <- thr_vec <- test_name_vec <- c()
#                 counter <- 1 
#                 for (t in 1:n_index_tests) {
#                   for (k in 1:n_thr_max) {
#                     test_vec[counter] <- t
#                     test_name_vec[counter] <- test_names[t]  # Use actual test names
#                     thr_vec[counter] <- k
#                     
#                     counter <- counter + 1
#                   }
#                 }
#                 ##
#                 n_rows_total <- length(test_vec)
#                 ##
#                 ## Make tibble with test names
#                 tibble_NMA <- tibble(  Model = rep(model_name, n_rows_total),
#                                        ##
#                                        test = factor(test_vec),
#                                        test_name = factor(test_name_vec, levels = test_names),  # Use actual test names
#                                        threshold = thr_vec,
#                                        ##
#                                        Se_median = c(t(Se_median_array)),
#                                        Sp_median = c(t(Sp_median_array)),
#                                        Fp_median = c(t(Fp_median_array)),
#                                        ##
#                                        Se_mean = c(t(Se_mean_array)),
#                                        Sp_mean = c(t(Sp_mean_array)),
#                                        Fp_mean = c(t(Fp_mean_array)),
#                                        ##
#                                        Se_lower = c(t(Se_lower_array)),
#                                        Sp_lower = c(t(Sp_lower_array)),
#                                        Fp_lower = c(t(Fp_lower_array)),
#                                        ##
#                                        Se_upper = c(t(Se_upper_array)),
#                                        Sp_upper = c(t(Sp_upper_array)),
#                                        Fp_upper = c(t(Fp_upper_array)),
#                                        ##
#                                        Se_pred_lower = c(t(Se_pred_lower_array)),
#                                        Sp_pred_lower = c(t(Sp_pred_lower_array)),
#                                        Fp_pred_lower = c(t(Fp_pred_lower_array)),
#                                        ##
#                                        Se_pred_upper = c(t(Se_pred_upper_array)),
#                                        Sp_pred_upper = c(t(Sp_pred_upper_array)),
#                                        Fp_pred_upper = c(t(Fp_pred_upper_array)))
#                 ##
#                 tibble_NMA <- dplyr::filter(tibble_NMA, Se_median != -1)
#                 ##
#                 if (!(is.null(relevant_thresholds))) { 
#                   tibble_NMA <-  filter_tibble_by_thresholds(  tibble_data = tibble_NMA,
#                                                                relevant_thresholds = relevant_thresholds)
#                 }
#                 ##
#                 ## Update polygon creation to use test names
#                 polygon_Conf_list <- polygon_Pred_list <- list()
#                 
#                 for (t in 1:n_index_tests) {
#                   
#                     tibble_test_t <- dplyr::filter(tibble_NMA, test == t)
#                     
#                     # Confidence region
#                     tibble_Conf_polygon_test_t <- create_confidence_polygon(tibble_test_t, 
#                                                                             model_name = model_name)
#                     n_rows = nrow(tibble_Conf_polygon_test_t)
#                     tibble_Conf_polygon_test_t <- tibble_Conf_polygon_test_t %>% 
#                       dplyr::mutate( Model = rep(model_name, n_rows),
#                                      test = rep(t, n_rows),
#                                      test_name = factor(rep(test_names[t], n_rows), levels = test_names))
#                     polygon_Conf_list[[t]] <- tibble_Conf_polygon_test_t
#                     
#                     # Prediction region
#                     tibble_Pred_polygon_test_t <- create_prediction_region_smooth(tibble_test_t,
#                                                                             model_name = model_name)
#                     n_rows = nrow(tibble_Pred_polygon_test_t)
#                     tibble_Pred_polygon_test_t <- tibble_Pred_polygon_test_t %>% 
#                       dplyr::mutate( Model = rep(model_name, n_rows),
#                                      test = rep(t, n_rows), 
#                                     test_name = factor(rep(test_names[t], n_rows), levels = test_names))
#                     polygon_Pred_list[[t]] <- tibble_Pred_polygon_test_t
#                     
#                 }
#                 
#                 polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% dplyr::filter(!(is.na(x)))
#                 polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% dplyr::filter(!(is.na(x)))
#                 
#                   plot_list <- list()
#                   
#                   # Plot 1: All tests on one plot (use test names for color legend)
#                   plot_1 <-  ggplot(tibble_NMA, 
#                                     mapping = aes(x = Fp_median, 
#                                                   y = Se_median,
#                                                   colour = test_name)) + 
#                     geom_line(linewidth = linewidth) + 
#                     geom_point(size = point_size) + 
#                     theme_bw(base_size = base_size) + 
#                     xlab("False positive rate (Fp)") + 
#                     ylab("Sensitivity (Se)") + 
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) + 
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     labs(colour = "Test")
#                   print(plot_1)
#                   plot_list[[1]] <- plot_1
#                   
#                   # Plot 2: Each test on separate panel (use test names for facets)
#                   plot_2 <-  ggplot(tibble_NMA, 
#                                     mapping = aes(x = Fp_median,
#                                                   y = Se_median,
#                                                   colour = Model)) + 
#                     geom_line(linewidth = linewidth) + 
#                     geom_point(size = point_size) + 
#                     theme_bw(base_size = base_size) + 
#                     facet_wrap(~ test_name) +  # Use test_name instead of test_char
#                     xlab("False positive rate (Fp)") + 
#                     ylab("Sensitivity (Se)") + 
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) + 
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) + 
#                     theme(legend.position = "none")
#                   print(plot_2)
#                   plot_list[[2]] <- plot_2
#                   
#                   ##
#                   ## ---- Plot 3 (plot w/ 95% confidence region):
#                   ##
#                   plot_3 <-   plot_2  +
#                     ##
#                     geom_polygon(data = polygon_Conf_tibble,
#                                  aes(x = x, y = y, 
#                                      colour = Model),
#                                  fill = conf_region_colour,
#                                  alpha = conf_region_alpha)  +
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1))
#                   # xlab("False positive rate (Fp)") +
#                   # ylab("Sensitivity (Se)")
#                   # geom_polygon(data = polygon_Pred_tibble, aes(x = x, y = y), fill = pred_region_colour, alpha = 0.25) +
#                   ##
#                   print(plot_3)
#                   plot_list[[3]] <- plot_3
#                   ##
#                   ## ---- Plot 4 (plot w/ 95% prediction region):
#                   ##
#                   plot_4 <-  plot_2 +  
#                     geom_polygon(
#                                 data = polygon_Pred_tibble,
#                                 aes(x = x, y = y, colour = Model),
#                                 fill = pred_region_colour,
#                                 alpha = conf_region_alpha)  +
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1))
#                   # xlab("False positive rate (Fp)") +
#                   # ylab("Sensitivity (Se)")
#                   print(plot_4)
#                   plot_list[[4]] <- plot_4
#                   ##
#                   ## ---- Plot 5 (plot w/ BOTH the 95% confidence + prediction regions):
#                   ##
#                   plot_5 <-  plot_2 +
#                     geom_polygon(
#                                 data = polygon_Conf_tibble,
#                                 aes(x = x, y = y, colour = Model),
#                                 fill = conf_region_colour,
#                                 alpha = conf_region_alpha) +
#                     ##
#                     geom_polygon(
#                       data = polygon_Pred_tibble,
#                       aes(x = x, y = y, colour = Model),
#                       fill = pred_region_colour,
#                       alpha = pred_region_alpha) +
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1))
#                   # xlab("False positive rate (Fp)") +
#                   # ylab("Sensitivity (Se)")
#                   print(plot_5)
#                   plot_list[[5]] <- plot_5
#                   ##
#                   ## ---- Plot 6 (plot w/ BOTH the 95% confidence + prediction regions):
#                   ##
#                   plot_6 <-  plot_1 +
#                     geom_polygon(
#                       data = polygon_Conf_tibble,
#                       aes(x = x, 
#                           y = y,
#                           # colour = test_name,
#                           fill = test_name),
#                       alpha = pred_region_alpha) +
#                     ##
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1))
#                   # xlab("False positive rate (Fp)") +
#                   # ylab("Sensitivity (Se)")
#                   print(plot_6)
#                   plot_list[[6]] <- plot_6
#                   ##
#                   ## ---- Plot 7 (plots of Se and Sp vs. threshold for each test - base plot)
#                   ##
#                   # Normalize thresholds for each test to 0-1 range
#                   tibble_NMA_normalized <- tibble_NMA %>%
#                     group_by(test_name) %>%
#                     mutate(
#                       # Find actual threshold range for each test
#                       min_thr = min(threshold[!is.na(Se_median)]),
#                       max_thr = max(threshold[!is.na(Se_median)]),
#                       # Normalize to 0-1
#                       threshold_norm = (threshold - min_thr) / (max_thr - min_thr),
#                       # Create labels for key points
#                       threshold_label = ifelse(
#                         threshold %in% c(min_thr, max_thr, round((min_thr + max_thr)/2)),
#                         as.character(threshold),
#                         ""
#                       )
#                     ) %>%
#                     ungroup()
#                   ##
#                   plot_Se_vs_threshold_base <-   ggplot( tibble_NMA_normalized, 
#                                                     mapping = aes(x = threshold_norm, 
#                                                                   y = Se_median, 
#                                                                   colour = test_name)) + 
#                     geom_line(linewidth = linewidth) + 
#                     geom_point(size = point_size) + 
#                     theme_bw(base_size = base_size) + 
#                     xlab("Test threshold") + 
#                     ylab("Sensitivity (Se)") + 
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     theme(
#                       axis.text.x = element_blank(),  # Remove tick labels
#                       axis.ticks.x = element_blank()  # Remove tick marks
#                     ) + 
#                     labs(colour = "Test") + 
#                     theme(legend.position = "bottom")
#                   ##
#                   plot_Se_vs_threshold_base
#                   ##
#                   plot_Se_vs_threshold_base_w_CrI <- plot_Se_vs_threshold_base +   
#                                                  geom_ribbon(aes( x = threshold_norm, 
#                                                                   ymin = Se_lower, 
#                                                                   ymax = Se_upper,
#                                                                   fill = test_name),
#                                                              alpha = conf_region_alpha) 
#                   ##
#                   plot_Se_vs_threshold_panel <- plot_Se_vs_threshold_base_w_CrI + 
#                     facet_wrap(~ test_name, 
#                                scales = "free", 
#                                ncol = 1) + 
#                     geom_ribbon(aes( x = threshold_norm, 
#                                      ymin = Se_pred_lower, 
#                                      ymax = Se_pred_upper,
#                                      fill = test_name),
#                                 alpha = pred_region_alpha) + 
#                     theme(legend.position = "none")
#                   plot_Se_vs_threshold_panel
#                   ##
#                   ##
#                   plot_Sp_vs_threshold_base <-   ggplot( tibble_NMA_normalized, 
#                                                          mapping = aes(x = threshold_norm, 
#                                                                        y = Sp_median, 
#                                                                        colour = test_name)) + 
#                     geom_line(linewidth = linewidth) + 
#                     geom_point(size = point_size) + 
#                     theme_bw(base_size = base_size) + 
#                     xlab("Test threshold") + 
#                     ylab("Specificity (Sp)") + 
#                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                                        limits = c(0, 1)) +
#                     theme(
#                       axis.text.x = element_blank(),  # Remove tick labels
#                       axis.ticks.x = element_blank()  # Remove tick marks
#                     ) + 
#                     labs(colour = "Test") + 
#                     theme(legend.position = "bottom")
#                   ##
#                   plot_Sp_vs_threshold_base
#                   ##
#                   plot_Sp_vs_threshold_base_w_CrI <- plot_Sp_vs_threshold_base +   
#                                                   geom_ribbon(aes( x = threshold_norm, 
#                                                                    ymin = Sp_lower, 
#                                                                    ymax = Sp_upper,
#                                                                    fill = test_name),
#                                                                    alpha = conf_region_alpha)
#                   plot_Sp_vs_threshold_base_w_CrI
#                   ##
#                   plot_Sp_vs_threshold_panel <- plot_Sp_vs_threshold_base_w_CrI + 
#                     facet_wrap(~ test_name, 
#                                scales = "free", 
#                                ncol = 1) + 
#                     geom_ribbon(aes( x = threshold_norm, 
#                                      ymin = Sp_pred_lower, 
#                                      ymax = Sp_pred_upper,
#                                      fill = test_name),
#                                 alpha = pred_region_alpha) + 
#                     theme(legend.position = "none")
#                   plot_Sp_vs_threshold_panel
#                   ##
#                   plot_7 <- plot_Se_vs_threshold_base + 
#                     theme(legend.position = "none") + 
#                     xlab(" ") +  
#                     plot_Sp_vs_threshold_base + 
#                     patchwork::plot_layout(ncol = 1)
#                   ##
#                   print(plot_7)
#                   plot_list[[7]] <- plot_7
#                   ##
#                   plot_8 <- plot_Se_vs_threshold_base_w_CrI + 
#                     theme(legend.position = "none") + 
#                     xlab(" ") + 
#                     plot_Sp_vs_threshold_base_w_CrI + 
#                     patchwork::plot_layout(ncol = 1)
#                   ##
#                   print(plot_8)
#                   plot_list[[8]] <- plot_8
#                   ##
#                   plot_9 <- plot_Se_vs_threshold_panel + plot_Sp_vs_threshold_panel
#                   print(plot_9)
#                   plot_list[[9]] <- plot_9
#                   ##
#                 
#                   return(list(  plot_list = plot_list,
#                                 ##
#                                 tibble_NMA = tibble_NMA,
#                                 ##
#                                 polygon_Conf_list   = polygon_Conf_list,
#                                 polygon_Conf_tibble = polygon_Conf_tibble,
#                                 polygon_Pred_list   = polygon_Pred_list,
#                                 polygon_Pred_tibble = polygon_Pred_tibble))
#   
# }
# 





 
# Calculate which test is "higher" on sROC curve
#' @export
identify_sroc_position <- function(diff_data, 
                                   print = FALSE
) {
        
        diff_data_filtered <- diff_data %>%
          mutate(
            # Positive = test1 is higher/left on sROC
            # Negative = test2 is higher/left on sROC
            sroc_position_score = diff_Se_median - diff_Sp_median,
            trade_off_ratio = abs(diff_Se_median) / abs(diff_Sp_median),
            high_value = ((trade_off_ratio > 3) | (trade_off_ratio < 0.33)),
            keep_comparison = high_value) %>% 
          filter(keep_comparison == 1)
        
          #   keep_comparison = abs(sroc_position_score) < 0.05 | abs(sroc_position_score) > 0.15) %>%
          # 
        
        if (print) {
            diff_data %>% print(n = 100000)
            diff_data_filtered %>% print(n = 100000)
        }
        
        return(diff_data_filtered)
  
}



#' find_optimal_threshold_ranges
#' @export
find_optimal_threshold_ranges <- function(diff_data, 
                                          max_extreme_se = 0.25,  # Stop if Se diff > 25%
                                          max_extreme_sp = 0.25,  # Stop if Sp diff > 25%
                                          min_comparisons = 3) {   # Show at least 3 per block
  
        diff_data %>%
          # Sort by thresholds to ensure progression
          arrange(test1, test2, threshold1, threshold2) %>%
          group_by(test1_name, test2_name, threshold1) %>%
          mutate(
            # Calculate "extremeness" - how far from reasonable trade-off
            is_extreme = abs(diff_Se_median) > max_extreme_se & 
              abs(diff_Sp_median) > max_extreme_sp,
            
            # Check if heavily dominated (both metrics strongly favor one test)
            is_dominated = sign(diff_Se_median) == sign(diff_Sp_median) & 
              abs(diff_Se_median) > 0.15 & 
              abs(diff_Sp_median) > 0.15,
            
            # Find cumulative extremeness (stop when we hit 2+ extreme comparisons)
            cumulative_extreme = cumsum(is_extreme | is_dominated),
            
            # Mark where pattern "breaks" (big jump in values)
            se_change = abs(diff_Se_median - lag(diff_Se_median, default = 0)),
            sp_change = abs(diff_Sp_median - lag(diff_Sp_median, default = 0)),
            has_jump = se_change > 0.15 | sp_change > 0.15,
            
            # Determine optimal range
            should_include = cumulative_extreme < 2 | row_number() <= min_comparisons
          ) %>%
          filter(should_include) %>%
          # Get the ranges
          summarise(
            threshold2_min = min(threshold2),
            threshold2_max = max(threshold2),
            n_comparisons = n(),
            se_range = max(diff_Se_median) - min(diff_Se_median),
            sp_range = max(diff_Sp_median) - min(diff_Sp_median),
            .groups = "drop"
          )
  
}








#' plot_se_differences
#' @export
plot_se_differences <- function(tibble, 
                                relevant_thresholds = NULL,
                                filter = FALSE,
                                plot_options = NULL,
                                plot_title = "Sensitivity Differences Between Tests/Thresholds") {
  
        library(forcats)
        ##
        if (is.null(plot_settings$base_size)) plot_settings$base_size <- 14
        ##
        if (!(is.null(relevant_thresholds))) { 
          filter <- TRUE
        }
        ##
        if (filter) { 
          filtered_tibble <- filter_clinically_relevant_comparisons( tibble_comp = tibble,
                                                                       relevant_thresholds = relevant_thresholds) 
        }
        ##
        plot_data <- filtered_tibble %>%
          arrange(test1, threshold1, test2, threshold2) %>%
          mutate(
            # Create a comparison type based on what tests are being compared
            comparison_type = paste0(test1_name, " vs ", test2_name),
            # Create comparison label if it doesn't exist
            comparison_label = if("comparison_label" %in% 
                                  names(.)) comparison_label 
                               else paste0(test1_name, " (", threshold1, ") - ", test2_name, " (", threshold2, ")"),
            order_within_facet = row_number(),
            comparison_label_ordered = fct_reorder(comparison_label, order_within_facet, .desc = TRUE)
          ) %>%
          ungroup()
        
        # if (filter) { 
        #   plot_data <- identify_sroc_position(plot_data) 
        # }
        
        plot <- ggplot(plot_data, aes(y = comparison_label_ordered, 
                                      x = diff_Se_median)) +
                geom_vline(xintercept = 0,
                           linetype = "dashed", 
                           alpha = 0.5) +
                geom_errorbarh(aes(xmin = diff_Se_lower, 
                                   xmax = diff_Se_upper), 
                               height = 0.2,
                               alpha = 0.5) +
                geom_point(size = 3) +
                facet_wrap(~ comparison_type, 
                           drop = FALSE,
                           scales = "free_y", ncol = 1, dir = "v") +
                labs(x = "Difference in Sensitivity", 
                     y = "",
                     title = plot_title) +
                theme_bw()
        
        plot <- plot + theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = plot_options$base_size),
          plot.title = element_text(size = plot_options$base_size)
        )
        
        return(list(plot_data = plot_data,
                    plot = plot))
  
}

#' plot_sp_differences
#' @export
plot_sp_differences <- function(tibble, 
                                relevant_thresholds = NULL,
                                filter = FALSE,
                                plot_options = NULL,
                                plot_title = "Specificity Differences Between Tests/Thresholds") {
  
        library(forcats)
        ##
        if (is.null(plot_options$base_size)) plot_options$base_size <- 14
        ##
        if (!(is.null(relevant_thresholds))) { 
          filter <- TRUE
        }
        ##
        if (filter) { 
          filtered_tibble <- filter_clinically_relevant_comparisons( tibble_comp = tibble,
                                                                     relevant_thresholds = relevant_thresholds) 
        }
        ##
        plot_data <- filtered_tibble %>%
          arrange(test1, threshold1, test2, threshold2) %>%
          mutate(
            # Create a comparison type based on what tests are being compared
            comparison_type = paste0(test1_name, " vs ", test2_name),
            # Create comparison label if it doesn't exist
            comparison_label = if("comparison_label" %in% 
                                  names(.)) comparison_label 
                               else paste0(test1_name, " (", threshold1, ") - ", test2_name, " (", threshold2, ")"),
            order_within_facet = row_number(),
            comparison_label_ordered = fct_reorder(comparison_label, order_within_facet, .desc = TRUE)
          ) %>%
          ungroup()
        
        # if (filter) { 
        #   plot_data <- identify_sroc_position(plot_data) 
        # }
        
        plot <- ggplot(plot_data, aes(y = comparison_label_ordered, 
                                      x = diff_Sp_median)) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          geom_errorbarh(aes(xmin = diff_Sp_lower, 
                             xmax = diff_Sp_upper), 
                         height = 0.2, alpha = 0.5) +
          geom_point(size = 3) +
          facet_wrap(~ comparison_type,
                     drop = FALSE,
                     scales = "free_y", ncol = 1) +
          labs(x = "Difference in Specificity", 
               y = "",
               title = plot_title) +
          theme_bw()
        
        plot <- plot + theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = plot_options$base_size),
          plot.title = element_text(size = plot_options$base_size)
        )
        
        return(list(plot_data = plot_data,
                    plot = plot))
  
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
                                                output_dir = "MR_results",
                                                n_cores = NULL,  # NULL = use all available cores - 1
                                                parallel = TRUE  # Option to run sequentially if needed
) {
        
        require(dplyr)
        require(purrr)
        ##
        # if (parallel) { 
          require(foreach)
          require(doParallel)
        # }
        ##
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
        
        # Set up parallel backend if requested
        if (parallel) {
              
              # Determine number of cores
              if (is.null(n_cores)) {
                n_cores <- parallel::detectCores() - 1  # Leave one core free
              }
              n_cores <- min(n_cores, n_scenarios)  # Don't use more cores than scenarios
              
              # Register parallel backend
              cl <- makeCluster(n_cores)
              registerDoParallel(cl)
              
              cat("Processing", n_scenarios, "meta-regression scenarios in parallel using", n_cores, "cores...\n")
          
        } else {
          
              # Register sequential backend
              registerDoSEQ()
              cat("Processing", n_scenarios, "meta-regression scenarios sequentially...\n")
          
        }
        
        if (parallel) {
          
                # Process scenarios in parallel
                scenario_results <- foreach(
                  i = 1:n_scenarios,
                  .packages = c("dplyr", "purrr"),  # Add any other required packages
                  .export = c("R_fn_compute_MR_complete_new_baseline_analysis_v1"),  # Export the function
                  .errorhandling = "pass",  # Return error messages instead of stopping
                  .verbose = FALSE
                ) %dopar% {
                  
                        tryCatch({
                          cat("\nProcessing Scenario", i, ":", scenario_names[i], "on worker", Sys.getpid(), "\n")
                          
                                # Compute results for this scenario
                                results <-  R_fn_compute_MR_complete_new_baseline_analysis_v1(
                                            debugging = FALSE,
                                            model_prep_obj = model_prep_obj,
                                            model_summary_and_trace_obj = model_summary_and_trace_obj,
                                            baseline_case_nd = baseline_case_nd_list[[i]],
                                            baseline_case_d = baseline_case_d_list[[i]],
                                            use_probit_link = use_probit_link,
                                            test_names = test_names,
                                            compute_comparisons = compute_comparisons)
                                
                                # Add scenario info
                                results$scenario_name <- scenario_names[i]
                                results$scenario_index <- i
                                
                                # Save individual results if requested
                                if (save_results) {
                                  saveRDS(results, 
                                          file.path(output_dir, 
                                                    paste0("scenario_", i, "_", 
                                                           gsub("[^A-Za-z0-9]", "_", scenario_names[i]), 
                                                           ".rds")))
                                }
                                
                                return(results)
                          
                        }, error = function(e) {
                          
                                # Return error info instead of stopping everything
                                return(list(
                                  error = TRUE,
                                  scenario_name = scenario_names[i],
                                  scenario_index = i,
                                  error_message = as.character(e)
                                ))
                          
                        })
                }
                
                # Stop cluster if we created one
                if (parallel) {
                  stopCluster(cl)
                }
                
                # Name the results list
                names(scenario_results) <- scenario_names
                
                # Check for errors
                errors <- sapply(scenario_results, function(x) !is.null(x$error) && x$error)
                if (any(errors)) {
                  cat("\nWarning: The following scenarios encountered errors:\n")
                  for (i in which(errors)) {
                    cat("  -", scenario_names[i], ":", scenario_results[[i]]$error_message, "\n")
                  }
                }
                
                cat("\nCompleted processing", sum(!errors), "out of", n_scenarios, "scenarios successfully.\n")
        
        } else { 
                
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




 



# 
# 
# 
# 
#  
# R_fn_plot_MR_comparison_heatmap <- function(
#                                             all_comparisons_data,
#                                             comparison_type = c("Se", "Sp"),
#                                             n_comparisons = 20,
#                                             base_size = 14,
#                                             save_plot = TRUE,
#                                             output_dir = "MR_plots",
#                                             plot_width = 12,
#                                             plot_height = 10
# ) {
#   
#         require(ggplot2)
#         require(dplyr)
#         require(tidyr)
#         require(viridis)
#         
#         comparison_type <- match.arg(comparison_type)
#         
#         # Select top comparisons by absolute median difference
#         if (comparison_type == "Se") {
#           top_comparisons <- all_comparisons_data %>%
#             group_by(test1_name, test2_name, threshold1, threshold2) %>%
#             summarise(
#               mean_abs_diff = mean(abs(diff_Se_median), na.rm = TRUE),
#               .groups = "drop"
#             ) %>%
#             arrange(desc(mean_abs_diff)) %>%
#             slice_head(n = n_comparisons) %>%
#             mutate(
#               comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
#                                         test2_name, " (", threshold2, ")")
#             ) %>%
#             pull(comparison_label)
#           
#           # Prepare data for heatmap
#           heatmap_data <- all_comparisons_data %>%
#             mutate(
#               comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
#                                         test2_name, " (", threshold2, ")")
#             ) %>%
#             filter(comparison_label %in% top_comparisons) %>%
#             select(comparison_label, scenario, diff_Se_median) %>%
#             pivot_wider(names_from = scenario, values_from = diff_Se_median)
#           
#           # Convert to matrix for heatmap
#           heatmap_matrix <- as.matrix(heatmap_data[, -1])
#           rownames(heatmap_matrix) <- heatmap_data$comparison_label
#           
#         } else {
#           top_comparisons <- all_comparisons_data %>%
#             group_by(test1_name, test2_name, threshold1, threshold2) %>%
#             summarise(
#               mean_abs_diff = mean(abs(diff_Sp_median), na.rm = TRUE),
#               .groups = "drop"
#             ) %>%
#             arrange(desc(mean_abs_diff)) %>%
#             slice_head(n = n_comparisons) %>%
#             mutate(
#               comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
#                                         test2_name, " (", threshold2, ")")
#             ) %>%
#             pull(comparison_label)
#           
#           heatmap_data <- all_comparisons_data %>%
#             mutate(
#               comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
#                                         test2_name, " (", threshold2, ")")
#             ) %>%
#             filter(comparison_label %in% top_comparisons) %>%
#             select(comparison_label, scenario, diff_Sp_median) %>%
#             pivot_wider(names_from = scenario, values_from = diff_Sp_median)
#           
#           heatmap_matrix <- as.matrix(heatmap_data[, -1])
#           rownames(heatmap_matrix) <- heatmap_data$comparison_label
#         }
#         
#         # Create heatmap using ggplot2
#         heatmap_long <- heatmap_data %>%
#           pivot_longer(-comparison_label, names_to = "scenario", values_to = "difference")
#         
#         p <- ggplot(heatmap_long, 
#                     aes(x = scenario, y = comparison_label, fill = difference)) +
#           geom_tile() +
#           scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
#                                midpoint = 0,
#                                name = paste0("Difference in ", 
#                                              ifelse(comparison_type == "Se", "Sensitivity", "Specificity"))) +
#           theme_bw(base_size = base_size) +
#           theme(
#             axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#             axis.text.y = element_text(size = 9),
#             axis.title = element_blank(),
#             legend.position = "right"
#           ) +
#           labs(title = paste0("Effect of Meta-regression Covariates on ",
#                               ifelse(comparison_type == "Se", "Sensitivity", "Specificity"),
#                               " Differences"))
#         
#         if (save_plot) {
#           if (!dir.exists(output_dir)) {
#             dir.create(output_dir, recursive = TRUE)
#           }
#           
#           ggsave(
#             file.path(output_dir, paste0("heatmap_", comparison_type, "_comparisons.pdf")),
#             p,
#             width = plot_width,
#             height = plot_height
#           )
#         }
#         
#         return(p)
#   
# }
# 
# 




 




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




#' @keywords internal
#' @export
clean_grid_axes <- function(patchwork_plot, ncol = NULL, remove_titles = FALSE) {
  
  plot_list <- patchwork_plot$patches$plots
  n_plots <- length(plot_list)
  
  if (is.null(ncol)) {
    if (n_plots == 4) ncol <- 2
    else if (n_plots == 9) ncol <- 3
    else ncol <- ceiling(sqrt(n_plots))
  }
  
  nrow <- ceiling(n_plots / ncol)
  
  for (i in 1:n_plots) {
    row <- ceiling(i / ncol)
    col <- ((i - 1) %% ncol) + 1
    
    # For interior plots, remove tick labels
    if (row < nrow) {  # Not bottom row
      plot_list[[i]] <- plot_list[[i]] + 
        theme(axis.text.x = element_blank())
      
      if (remove_titles) {
        plot_list[[i]] <- plot_list[[i]] + 
          theme(axis.title.x = element_blank())
      }
    }
    
    if (col > 1) {  # Not first column
      plot_list[[i]] <- plot_list[[i]] + 
        theme(axis.text.y = element_blank())
      
      if (remove_titles) {
        plot_list[[i]] <- plot_list[[i]] + 
          theme(axis.title.y = element_blank())
      }
    }
  }
  
  # Make sure we still have all plots
  if (length(plot_list) != n_plots) {
    stop("Lost some plots! Started with ", n_plots, " but now have ", length(plot_list))
  }
  
  wrap_plots(plot_list, ncol = ncol) +
    plot_annotation(
      title = patchwork_plot$plot_annotation$title,
      subtitle = patchwork_plot$plot_annotation$subtitle,
      theme = patchwork_plot$plot_annotation$theme
    )
}




#' filter_tibble_by_thresholds
#' @export
filter_tibble_by_thresholds <- function(tibble_data, 
                                        relevant_thresholds) {
  
        if (is.list(relevant_thresholds)) {
          # Filter by test-specific thresholds
          filtered_data <- tibble_data %>%
            group_by(test) %>%
            filter(threshold %in% relevant_thresholds[[first(test_name)]]) %>%
            ungroup()
        } else {
          # Filter all tests by same thresholds
          filtered_data <- tibble_data %>%
            filter(threshold %in% relevant_thresholds)
        }
        return(filtered_data)
  
}





#' Create sROC plots with different layouts
#' @export
create_sROC_plots_flexible <- function(
                                      combined_data,
                                      ##
                                      test_names = NULL,
                                      relevant_thresholds = NULL,
                                      ##
                                      layout = c("separate", 
                                                 "by_test", "by_scenario",
                                                 "grid", "grid_separate"),
                                      ##
                                      aesthetic_mapping = NULL,  # List with parse_fn, color_var, linetype_var, shape_var
                                      ##
                                      separate_sROC_settings = NULL,
                                      ##
                                      by_test_sROC_settings = NULL,
                                      by_scenario_sROC_settings = NULL,
                                      ##
                                      by_test_grid_sROC_settings = NULL,
                                      by_scenario_grid_sROC_settings = NULL,
                                      ##
                                      by_test_grid_separate_sROC_settings = NULL,
                                      by_scenario_grid_separate_sROC_settings = NULL
) {
  
        require(ggplot2)
        require(dplyr)
        require(patchwork)
        ##
        scenario_levels <- unique(combined_data$scenario) ## Preserve original scenario order from the data
        ##
        # 
        # Add test names if provided, preserving order from relevant_thresholds
        if (!is.null(test_names)) {
          # Create test_name column using numeric test index
          combined_data <- combined_data %>%
            mutate(test_name = test_names[test])
          
          # If relevant_thresholds provided, use its order for factor levels
          if (!is.null(relevant_thresholds)) {
            test_name_order <- names(relevant_thresholds)
          } else {
            test_name_order <- test_names
          }
          
          # Set factor levels in the correct order
          combined_data <- combined_data %>%
            mutate(test_name = factor(test_name, levels = test_name_order))
          
        } else {
          combined_data <- combined_data %>%
            mutate(test_name = factor(paste0("Test ", test)))
        }
        ##
        ## DEBUG: Check data structure:
        ##
        cat("\n=== DATA DEBUG ===\n")
        cat("test_names provided:\n")
        print(test_names)
        cat("\nrelevant_thresholds order:\n")
        print(names(relevant_thresholds))
        cat("\nUnique test values in data:\n")
        print(unique(combined_data$test))
        cat("\ntest -> test_name mapping:\n")
        combined_data %>% 
          distinct(test, test_name) %>% 
          arrange(test) %>% 
          print()
        cat("=================\n\n")
        ##
        # After combining data, apply aesthetic mapping if provided
        if (!is.null(aesthetic_mapping)) {
          # Apply the parsing function if provided
          if (!is.null(aesthetic_mapping$parse_fn)) {
            combined_data <- aesthetic_mapping$parse_fn(combined_data)
          }
        }
        ##
        # Get unique values
        n_scenarios <- length(unique(combined_data$scenario))
        n_tests <- length(unique(combined_data$test))
        scenario_names <- unique(combined_data$scenario)
        ##
        is_single_scenario <- (n_scenarios == 1)
        ##
        ##
        # Filter thresholds using the helper function
        if (!is.null(relevant_thresholds)) {
          combined_data <- filter_tibble_by_thresholds(combined_data, relevant_thresholds)
        }
        ##     
        # separate_sROC_settings
        if (is.null(separate_sROC_settings)) separate_sROC_settings <- list()
        if (is.null(separate_sROC_settings$base_size)) separate_sROC_settings$base_size <- 18
        if (is.null(separate_sROC_settings$point_size)) separate_sROC_settings$point_size <- 2
        if (is.null(separate_sROC_settings$line_size)) separate_sROC_settings$line_size <- 0.8
        if (is.null(separate_sROC_settings$rows_in_legend)) separate_sROC_settings$rows_in_legend <- 3
        if (is.null(separate_sROC_settings$legend_position)) separate_sROC_settings$legend_position <- "none"
        if (is.null(separate_sROC_settings$add_regions)) separate_sROC_settings$add_regions <- TRUE
        if (is.null(separate_sROC_settings$show_conf)) separate_sROC_settings$show_conf <- TRUE
        if (is.null(separate_sROC_settings$show_pred)) separate_sROC_settings$show_pred <- TRUE
        if (is.null(separate_sROC_settings$conf_region_alpha)) separate_sROC_settings$conf_region_alpha <- 0.30
        if (is.null(separate_sROC_settings$pred_region_alpha)) separate_sROC_settings$pred_region_alpha <- 0.15
        if (is.null(separate_sROC_settings$conf_region_colour)) separate_sROC_settings$conf_region_colour <- "blue"
        if (is.null(separate_sROC_settings$pred_region_colour)) separate_sROC_settings$pred_region_colour <- "green"
        
        # by_test_sROC_settings
        if (is.null(by_test_sROC_settings)) by_test_sROC_settings <- list()
        if (is.null(by_test_sROC_settings$base_size)) by_test_sROC_settings$base_size <- 18
        if (is.null(by_test_sROC_settings$point_size)) by_test_sROC_settings$point_size <- 2
        if (is.null(by_test_sROC_settings$line_size)) by_test_sROC_settings$line_size <- 0.8
        if (is.null(by_test_sROC_settings$rows_in_legend)) by_test_sROC_settings$rows_in_legend <- 3
        if (is.null(by_test_sROC_settings$legend_position)) by_test_sROC_settings$legend_position <- "bottom"
        if (is.null(by_test_sROC_settings$add_regions)) by_test_sROC_settings$add_regions <- FALSE
        if (is.null(by_test_sROC_settings$show_conf)) by_test_sROC_settings$show_conf <- FALSE
        if (is.null(by_test_sROC_settings$show_pred)) by_test_sROC_settings$show_pred <- FALSE
        if (is.null(by_test_sROC_settings$conf_region_alpha)) by_test_sROC_settings$conf_region_alpha <- NULL
        if (is.null(by_test_sROC_settings$pred_region_alpha)) by_test_sROC_settings$pred_region_alpha <- NULL
        if (is.null(by_test_sROC_settings$conf_region_colour)) by_test_sROC_settings$conf_region_colour <- NULL
        if (is.null(by_test_sROC_settings$pred_region_colour)) by_test_sROC_settings$pred_region_colour <- NULL
        
        # by_scenario_sROC_settings
        if (is.null(by_scenario_sROC_settings)) by_scenario_sROC_settings <- list()
        if (is.null(by_scenario_sROC_settings$base_size)) by_scenario_sROC_settings$base_size <- 18
        if (is.null(by_scenario_sROC_settings$point_size)) by_scenario_sROC_settings$point_size <- 2
        if (is.null(by_scenario_sROC_settings$line_size)) by_scenario_sROC_settings$line_size <- 0.8
        if (is.null(by_scenario_sROC_settings$rows_in_legend)) by_scenario_sROC_settings$rows_in_legend <- 2
        if (is.null(by_scenario_sROC_settings$legend_position)) by_scenario_sROC_settings$legend_position <- "bottom"
        if (is.null(by_scenario_sROC_settings$add_regions)) by_scenario_sROC_settings$add_regions <- FALSE
        if (is.null(by_scenario_sROC_settings$show_conf)) by_scenario_sROC_settings$show_conf <- FALSE
        if (is.null(by_scenario_sROC_settings$show_pred)) by_scenario_sROC_settings$show_pred <- FALSE
        if (is.null(by_scenario_sROC_settings$conf_region_alpha)) by_scenario_sROC_settings$conf_region_alpha <- NULL
        if (is.null(by_scenario_sROC_settings$pred_region_alpha)) by_scenario_sROC_settings$pred_region_alpha <- NULL
        if (is.null(by_scenario_sROC_settings$conf_region_colour)) by_scenario_sROC_settings$conf_region_colour <- NULL
        if (is.null(by_scenario_sROC_settings$pred_region_colour)) by_scenario_sROC_settings$pred_region_colour <- NULL
        
        # by_test_grid_sROC_settings
        if (is.null(by_test_grid_sROC_settings)) by_test_grid_sROC_settings <- list()
        if (is.null(by_test_grid_sROC_settings$base_size)) by_test_grid_sROC_settings$base_size <- 18
        if (is.null(by_test_grid_sROC_settings$point_size)) by_test_grid_sROC_settings$point_size <- 2
        if (is.null(by_test_grid_sROC_settings$line_size)) by_test_grid_sROC_settings$line_size <- 0.8
        if (is.null(by_test_grid_sROC_settings$rows_in_legend)) by_test_grid_sROC_settings$rows_in_legend <- 3
        if (is.null(by_test_grid_sROC_settings$legend_position)) by_test_grid_sROC_settings$legend_position <- "bottom"
        if (is.null(by_test_grid_sROC_settings$add_regions)) by_test_grid_sROC_settings$add_regions <- FALSE
        if (is.null(by_test_grid_sROC_settings$show_conf)) by_test_grid_sROC_settings$show_conf <- FALSE
        if (is.null(by_test_grid_sROC_settings$show_pred)) by_test_grid_sROC_settings$show_pred <- FALSE
        if (is.null(by_test_grid_sROC_settings$conf_region_alpha)) by_test_grid_sROC_settings$conf_region_alpha <- NULL
        if (is.null(by_test_grid_sROC_settings$pred_region_alpha)) by_test_grid_sROC_settings$pred_region_alpha <- NULL
        if (is.null(by_test_grid_sROC_settings$conf_region_colour)) by_test_grid_sROC_settings$conf_region_colour <- NULL
        if (is.null(by_test_grid_sROC_settings$pred_region_colour)) by_test_grid_sROC_settings$pred_region_colour <- NULL
        if (is.null(by_test_grid_sROC_settings$ncol)) by_test_grid_sROC_settings$ncol <- 3
        if (is.null(by_test_grid_sROC_settings$facet_label_width)) by_test_grid_sROC_settings$facet_label_width <- 25
        
        # by_scenario_grid_sROC_settings
        if (is.null(by_scenario_grid_sROC_settings)) by_scenario_grid_sROC_settings <- list()
        if (is.null(by_scenario_grid_sROC_settings$base_size)) by_scenario_grid_sROC_settings$base_size <- 18
        if (is.null(by_scenario_grid_sROC_settings$point_size)) by_scenario_grid_sROC_settings$point_size <- 2
        if (is.null(by_scenario_grid_sROC_settings$line_size)) by_scenario_grid_sROC_settings$line_size <- 0.8
        if (is.null(by_scenario_grid_sROC_settings$rows_in_legend)) by_scenario_grid_sROC_settings$rows_in_legend <- 1
        if (is.null(by_scenario_grid_sROC_settings$legend_position)) by_scenario_grid_sROC_settings$legend_position <- "bottom"
        if (is.null(by_scenario_grid_sROC_settings$add_regions)) by_scenario_grid_sROC_settings$add_regions <- FALSE
        if (is.null(by_scenario_grid_sROC_settings$show_conf)) by_scenario_grid_sROC_settings$show_conf <- FALSE
        if (is.null(by_scenario_grid_sROC_settings$show_pred)) by_scenario_grid_sROC_settings$show_pred <- FALSE
        if (is.null(by_scenario_grid_sROC_settings$conf_region_alpha)) by_scenario_grid_sROC_settings$conf_region_alpha <- NULL
        if (is.null(by_scenario_grid_sROC_settings$pred_region_alpha)) by_scenario_grid_sROC_settings$pred_region_alpha <- NULL
        if (is.null(by_scenario_grid_sROC_settings$conf_region_colour)) by_scenario_grid_sROC_settings$conf_region_colour <- NULL
        if (is.null(by_scenario_grid_sROC_settings$pred_region_colour)) by_scenario_grid_sROC_settings$pred_region_colour <- NULL
        if (is.null(by_scenario_grid_sROC_settings$ncol)) by_scenario_grid_sROC_settings$ncol <- 2
        if (is.null(by_scenario_grid_sROC_settings$facet_label_width)) by_scenario_grid_sROC_settings$facet_label_width <- 25
        
        # by_test_grid_separate_sROC_settings
        if (is.null(by_test_grid_separate_sROC_settings)) by_test_grid_separate_sROC_settings <- list()
        if (is.null(by_test_grid_separate_sROC_settings$base_size)) by_test_grid_separate_sROC_settings$base_size <- 18
        if (is.null(by_test_grid_separate_sROC_settings$point_size)) by_test_grid_separate_sROC_settings$point_size <- 2
        if (is.null(by_test_grid_separate_sROC_settings$line_size)) by_test_grid_separate_sROC_settings$line_size <- 0.8
        if (is.null(by_test_grid_separate_sROC_settings$rows_in_legend)) by_test_grid_separate_sROC_settings$rows_in_legend <- 3
        if (is.null(by_test_grid_separate_sROC_settings$legend_position)) by_test_grid_separate_sROC_settings$legend_position <- "bottom"
        if (is.null(by_test_grid_separate_sROC_settings$add_regions)) by_test_grid_separate_sROC_settings$add_regions <- TRUE
        if (is.null(by_test_grid_separate_sROC_settings$show_conf)) by_test_grid_separate_sROC_settings$show_conf <- TRUE
        if (is.null(by_test_grid_separate_sROC_settings$show_pred)) by_test_grid_separate_sROC_settings$show_pred <- TRUE
        if (is.null(by_test_grid_separate_sROC_settings$conf_region_alpha)) by_test_grid_separate_sROC_settings$conf_region_alpha <- 0.30
        if (is.null(by_test_grid_separate_sROC_settings$pred_region_alpha)) by_test_grid_separate_sROC_settings$pred_region_alpha <- 0.15
        if (is.null(by_test_grid_separate_sROC_settings$conf_region_colour)) by_test_grid_separate_sROC_settings$conf_region_colour <- "blue"
        if (is.null(by_test_grid_separate_sROC_settings$pred_region_colour)) by_test_grid_separate_sROC_settings$pred_region_colour <- "green"
        if (is.null(by_test_grid_separate_sROC_settings$ncol)) by_test_grid_separate_sROC_settings$ncol <- 3
        
        # by_scenario_grid_separate_sROC_settings
        if (is.null(by_scenario_grid_separate_sROC_settings)) by_scenario_grid_separate_sROC_settings <- list()
        if (is.null(by_scenario_grid_separate_sROC_settings$base_size)) by_scenario_grid_separate_sROC_settings$base_size <- 18
        if (is.null(by_scenario_grid_separate_sROC_settings$point_size)) by_scenario_grid_separate_sROC_settings$point_size <- 2
        if (is.null(by_scenario_grid_separate_sROC_settings$line_size)) by_scenario_grid_separate_sROC_settings$line_size <- 0.8
        if (is.null(by_scenario_grid_separate_sROC_settings$rows_in_legend)) by_scenario_grid_separate_sROC_settings$rows_in_legend <- 1
        if (is.null(by_scenario_grid_separate_sROC_settings$legend_position)) by_scenario_grid_separate_sROC_settings$legend_position <- "bottom"
        if (is.null(by_scenario_grid_separate_sROC_settings$add_regions)) by_scenario_grid_separate_sROC_settings$add_regions <- TRUE
        if (is.null(by_scenario_grid_separate_sROC_settings$show_conf)) by_scenario_grid_separate_sROC_settings$show_conf <- TRUE
        if (is.null(by_scenario_grid_separate_sROC_settings$show_pred)) by_scenario_grid_separate_sROC_settings$show_pred <- TRUE
        if (is.null(by_scenario_grid_separate_sROC_settings$conf_region_alpha)) by_scenario_grid_separate_sROC_settings$conf_region_alpha <- 0.30
        if (is.null(by_scenario_grid_separate_sROC_settings$pred_region_alpha)) by_scenario_grid_separate_sROC_settings$pred_region_alpha <- 0.15
        if (is.null(by_scenario_grid_separate_sROC_settings$conf_region_colour)) by_scenario_grid_separate_sROC_settings$conf_region_colour <- "blue"
        if (is.null(by_scenario_grid_separate_sROC_settings$pred_region_colour)) by_scenario_grid_separate_sROC_settings$pred_region_colour <- "green"
        if (is.null(by_scenario_grid_separate_sROC_settings$ncol)) by_scenario_grid_separate_sROC_settings$ncol <- 2
        
        # Initialize plots with nested structure like threshold function
        plots <- list()
        plots$separate <- list()
        plots$by_test <- list()
        plots$by_scenario <- list()
        plots$grid <- list()
        plots$grid_separate <- list()
        
        # 1. Separate plots for each test-scenario combination
        if ("separate" %in% layout) {

              for (s in scenario_names) {
                
                    # Loop through tests in the correct order
                    test_order <- if (!is.null(relevant_thresholds)) {
                      seq_along(names(relevant_thresholds))
                    } else {
                      unique(combined_data$test)
                    }
                    for (t in test_order) {
                    ## for (t in unique(combined_data$test)) {
    
                      plot_data <- combined_data %>%
                        filter(scenario == s, test == t)
    
                      test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
    
                      p <- create_sroc_plot_base(
                        data = plot_data,
                        title = paste0(test_name, " - ", s),
                        base_size = separate_sROC_settings$base_size,
                        point_size = separate_sROC_settings$point_size,
                        line_size = separate_sROC_settings$line_size
                        ##
                      )  +  theme(
                        axis.text.y = element_text(size = separate_sROC_settings$base_size),
                        legend.position = separate_sROC_settings$legend_position,
                        strip.text = element_text(size = separate_sROC_settings$base_size, face = "bold"),
                        panel.spacing = unit(1, "lines")
                        ##
                      ) + guides(color = guide_legend(nrow = separate_sROC_settings$rows_in_legend))
    
                      # Add regions if requested
                      if (separate_sROC_settings$add_regions) {
                            
                            p <- add_regions_to_sroc(
                              plot = p,
                              tibble_NMA = plot_data,
                              test_names = test_names,
                              show_conf = separate_sROC_settings$show_conf,
                              show_pred = separate_sROC_settings$show_pred,
                              conf_region_alpha = separate_sROC_settings$conf_region_alpha,
                              pred_region_alpha = separate_sROC_settings$pred_region_alpha,
                              conf_region_colour = separate_sROC_settings$conf_region_colour,
                              pred_region_colour = separate_sROC_settings$pred_region_colour
                            )
                        
                      }
    
                      plots$separate[[paste0(s, "_", test_name)]] <- p
                      
                    }
                
              }

        }

        # 2. One plot per test, all scenarios overlaid
        if ("by_test" %in% layout) {
          
                # Loop through tests in the correct order
                test_order <- if (!is.null(relevant_thresholds)) {
                  seq_along(names(relevant_thresholds))
                } else {
                  unique(combined_data$test)
                }
                for (t in test_order) {
              # for (t in unique(combined_data$test)) {

                  plot_data <-  combined_data %>%
                                filter(test == t) %>%
                                mutate(scenario = factor(scenario, levels = scenario_levels)) 

                  test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))

                  # Check if we should use special aesthetic mapping
                  if (!is.null(aesthetic_mapping) &&
                      !is.null(aesthetic_mapping$color_var) &&
                      !is.null(aesthetic_mapping$linetype_var)) {
        
                            # Create plot with custom aesthetics
                            p <- ggplot(plot_data,
                                        aes_string(x = "Fp_median",
                                                   y = "Se_median",
                                                   color = aesthetic_mapping$color_var,
                                                   linetype = aesthetic_mapping$linetype_var,
                                                   group = "scenario")) +
                              geom_line(linewidth = by_test_sROC_settings$line_size) +
                              geom_point(aes_string(shape = aesthetic_mapping$shape_var),
                                         size = by_test_sROC_settings$point_size) +
                              labs(title = paste0("sROC curves for ", test_name),     # ADD THIS
                                   subtitle = "Custom aesthetic mapping",              # ADD THIS  
                                   x = "1 - Specificity",                              # ADD THIS
                                   y = "Sensitivity") +                                # ADD THIS
                              scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                                 limits = c(0, 1)) +
                              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                                 limits = c(0, 1))
                            
                            # Apply custom scales if provided
                            if (!is.null(aesthetic_mapping$color_scale)) {
                              p <- p + aesthetic_mapping$color_scale
                            }
                            if (!is.null(aesthetic_mapping$linetype_scale)) {
                              p <- p + aesthetic_mapping$linetype_scale
                            }
                            if (!is.null(aesthetic_mapping$shape_scale)) {
                              p <- p + aesthetic_mapping$shape_scale
                            }
                    

                  } else {

  
                            p <- create_sroc_plot_base(
                              data = plot_data,
                              color_by = "scenario",
                              title = paste0("sROC curves for ", test_name),
                              subtitle = "Colored by meta-regression scenario",
                              base_size = by_test_sROC_settings$base_size,
                              point_size = by_test_sROC_settings$point_size,
                              line_size = by_test_sROC_settings$line_size)

                  }

                  p <- p +
                  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                     limits = c(0, 1)) +
                  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                     limits = c(0, 1)) +
                   theme(
                          axis.text.y = element_text(size = by_test_sROC_settings$base_size),
                          legend.position = by_test_sROC_settings$legend_position,
                          strip.text = element_text(size = by_test_sROC_settings$base_size, face = "bold"),
                          panel.spacing = unit(1, "lines")
                          ##
                        ) + guides(color = guide_legend(nrow = by_test_sROC_settings$rows_in_legend))


                # Add regions if requested
                  if (by_test_sROC_settings$add_regions) {
                    p <- add_regions_to_sroc(
                      plot = p,
                      tibble_NMA = plot_data,
                      group_by = "scenario",  # ADD THIS LINE - group by scenario, not test!
                      test_names = test_names,
                      scenario_names = unique(plot_data$scenario),  # ADD THIS LINE
                      show_conf = by_test_sROC_settings$show_conf,
                      show_pred = by_test_sROC_settings$show_pred,
                      conf_region_alpha = by_test_sROC_settings$conf_region_alpha,
                      pred_region_alpha = by_test_sROC_settings$pred_region_alpha,
                      conf_region_colour = by_test_sROC_settings$conf_region_colour,
                      pred_region_colour = by_test_sROC_settings$pred_region_colour
                    )
                  }

                  plots$by_test[[paste0("by_test_", test_name)]] <- p
              }

        }
        
        # 3. One plot per scenario, all tests overlaid
        if ("by_scenario" %in% layout) {
    
                for (s in scenario_names) {
    
                          plot_data <- combined_data %>% filter(scenario == s)
            
                          p <- create_sroc_plot_base(
                            data = plot_data,
                            color_by = "test_name",
                            title = s,
                            subtitle = "Colored by test",
                            base_size = by_scenario_sROC_settings$base_size,
                            point_size = by_scenario_sROC_settings$point_size,
                            line_size = by_scenario_sROC_settings$line_size
                            ##
                          ) + theme(
                            axis.text.y = element_text(size = by_scenario_sROC_settings$base_size),
                            legend.position = by_scenario_sROC_settings$legend_position,
                            strip.text = element_text(size = by_scenario_sROC_settings$base_size, face = "bold"),
                            panel.spacing = unit(1, "lines")
                            ##
                          ) + guides(color = guide_legend(nrow = by_scenario_sROC_settings$rows_in_legend))
            
                          # Add regions if requested
                          if (by_scenario_sROC_settings$add_regions) {
                            
                                  p <- add_regions_to_sroc(
                                    plot = p,
                                    tibble_NMA = plot_data,
                                    test_names = test_names,
                                    show_conf = by_scenario_sROC_settings$show_conf,
                                    show_pred = by_scenario_sROC_settings$show_pred,
                                    conf_region_alpha = by_scenario_sROC_settings$conf_region_alpha,
                                    pred_region_alpha = by_scenario_sROC_settings$pred_region_alpha,
                                    conf_region_colour = NULL,
                                    pred_region_colour = NULL
                                  )
                                  ##
                                  # Make fill colors match the line colors
                                  # Get the current color scale values
                                  n_tests <- length(unique(plot_data$test_name))
                                  test_colors <- scales::hue_pal()(n_tests)
                                  ##
                                  p <- p + scale_fill_manual(values = test_colors, 
                                                             guide = "none")  # Hide from legend since color already shows it
                                  
                          }
            
                          plots$by_scenario[[paste0("by_scenario_", s)]] <- p
                      
                }

        }
        
        combined_data <- combined_data %>%
                         mutate(scenario = factor(scenario, levels = scenario_levels))
        
        # 4. Grid layout
        if ("grid" %in% layout) {
          
              # Grid with tests as panels, scenarios as colors
              p1 <- create_sroc_plot_base(
                data = combined_data,
                color_by = "scenario",
                facet_by = "test_name",
                title = "sROC curves by test",
                subtitle = "Scenarios shown in different colors",
                base_size = by_test_grid_sROC_settings$base_size,
                point_size = by_test_grid_sROC_settings$point_size,
                line_size = by_test_grid_sROC_settings$line_size
              ) + 
                facet_wrap(~ test_name,
                           drop = FALSE,
                           ncol = by_test_grid_sROC_settings$ncol, 
                           dir = "v", 
                           labeller = label_wrap_gen(width = by_test_grid_sROC_settings$facet_label_width)) +
                theme(
                  axis.text.y = element_text(size = by_test_grid_sROC_settings$base_size),
                  legend.position = by_test_grid_sROC_settings$legend_position,
                  strip.text = element_text(size = by_test_grid_sROC_settings$base_size, face = "bold"),
                  panel.spacing = unit(1, "lines")
                ) + guides(color = guide_legend(nrow = by_test_grid_sROC_settings$rows_in_legend))
              
              plots$grid[["grid_by_test"]] <- p1
              
              # Grid with scenarios as panels, tests as colors
              p2 <- create_sroc_plot_base(
                data = combined_data,
                color_by = "test_name",
                facet_by = "scenario",
                title = "sROC curves by scenario",
                subtitle = "Tests shown in different colors",
                base_size = by_scenario_grid_sROC_settings$base_size,
                point_size = by_scenario_grid_sROC_settings$point_size,
                line_size = by_scenario_grid_sROC_settings$line_size
              ) + 
                facet_wrap(~ scenario, 
                           drop = FALSE,
                           ncol = by_scenario_grid_sROC_settings$ncol, 
                           dir = "v", 
                           labeller = label_wrap_gen(width = by_scenario_grid_sROC_settings$facet_label_width)) +
                theme(
                  axis.text.y = element_text(size = by_scenario_grid_sROC_settings$base_size),
                  legend.position = by_scenario_grid_sROC_settings$legend_position,
                  strip.text = element_text(size = by_scenario_grid_sROC_settings$base_size * 0.7, face = "bold"),
                  panel.spacing = unit(0.5, "lines")
                ) + guides(color = guide_legend(nrow = by_scenario_grid_sROC_settings$rows_in_legend))
              
              plots$grid[["grid_by_scenario"]] <- p2
          
        }
        
        # 5. Grid-separate layout
        if ("grid_separate" %in% layout) {
          
          # Create separate grids for each test showing all scenarios
          # Loop through tests in the correct order
          test_order <- if (!is.null(relevant_thresholds)) {
            seq_along(names(relevant_thresholds))
          } else {
            unique(combined_data$test)
          }
          for (t in test_order) {
        # for (t in unique(combined_data$test)) {
                
                test_data <- combined_data %>%
                  filter(test == t)
                
                test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                
                # Create individual plots for each scenario for this test
                scenario_plots <- list()
                
                for (s in scenario_names) {
                  
                      scenario_data <- test_data %>%
                        filter(scenario == s)
                      
                      p <- create_sroc_plot_base(
                        data = scenario_data,
                        title = s,
                        base_size = by_test_grid_separate_sROC_settings$base_size,
                        point_size = by_test_grid_separate_sROC_settings$point_size,
                        line_size = by_test_grid_separate_sROC_settings$line_size
                      ) + theme(
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none",
                        plot.title = element_text(size = by_test_grid_separate_sROC_settings$base_size * 0.8)
                      )
                      
                      # Add regions if requested
                      if (by_test_grid_separate_sROC_settings$add_regions) {
                        p <- add_regions_to_sroc(
                          plot = p,
                          tibble_NMA = scenario_data,
                          test_names = test_names,
                          show_conf = by_test_grid_separate_sROC_settings$show_conf,
                          show_pred = by_test_grid_separate_sROC_settings$show_pred,
                          conf_region_alpha = by_test_grid_separate_sROC_settings$conf_region_alpha,
                          pred_region_alpha = by_test_grid_separate_sROC_settings$pred_region_alpha,
                          conf_region_colour = by_test_grid_separate_sROC_settings$conf_region_colour,
                          pred_region_colour = by_test_grid_separate_sROC_settings$pred_region_colour
                        )
                      }
                      
                      scenario_plots[[s]] <- p
                  
                }
                
                # Arrange in grid
                grid_plot <- wrap_plots(scenario_plots, ncol = by_test_grid_separate_sROC_settings$ncol) +
                  plot_annotation(
                    title = paste0("sROC curves for ", test_name, " across all scenarios"),
                    caption = "X-axis: 1 - Specificity | Y-axis: Sensitivity",
                    theme = theme(
                      plot.title = element_text(size = by_test_grid_separate_sROC_settings$base_size + 2, face = "bold"),
                      plot.caption = element_text(size = 14, hjust = 0.5)
                    )
                  )
                
            plots$grid_separate[[paste0("grid_separate_test_", test_name)]] <- grid_plot
            
          }
          
          # Also create separate grids for each scenario showing all tests
          for (s in scenario_names) {
            
            scenario_data <- combined_data %>%
              filter(scenario == s)
            
            # Create individual plots for each test for this scenario
            test_plots <- list()
            
            
            # Loop through tests in the correct order
            test_order <- if (!is.null(relevant_thresholds)) {
              seq_along(names(relevant_thresholds))
            } else {
              unique(combined_data$test)
            }
            
            for (t in test_order) {
            # for (t in unique(combined_data$test)) {
              test_data <- scenario_data %>%
                filter(test == t)
              
              test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
              
              p <- create_sroc_plot_base(
                data = test_data,
                title = test_name,
                base_size = by_scenario_grid_separate_sROC_settings$base_size,
                point_size = by_scenario_grid_separate_sROC_settings$point_size,
                line_size = by_scenario_grid_separate_sROC_settings$line_size
              ) + theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none",
                plot.title = element_text(size = by_scenario_grid_separate_sROC_settings$base_size * 0.8)
              )
              
              # Add regions if requested
              if (by_scenario_grid_separate_sROC_settings$add_regions) {
                p <- add_regions_to_sroc(
                  plot = p,
                  tibble_NMA = test_data,
                  test_names = test_names,
                  show_conf = by_scenario_grid_separate_sROC_settings$show_conf,
                  show_pred = by_scenario_grid_separate_sROC_settings$show_pred,
                  conf_region_alpha = by_scenario_grid_separate_sROC_settings$conf_region_alpha,
                  pred_region_alpha = by_scenario_grid_separate_sROC_settings$pred_region_alpha,
                  conf_region_colour = by_scenario_grid_separate_sROC_settings$conf_region_colour,
                  pred_region_colour = by_scenario_grid_separate_sROC_settings$pred_region_colour
                )
              }
              
              test_plots[[test_name]] <- p
            }
            
            # Arrange in grid
            grid_plot <- wrap_plots(test_plots, ncol = by_scenario_grid_separate_sROC_settings$ncol) +
              plot_annotation(
                title = paste0("sROC curves for scenario: ", s),
                caption = "X-axis: 1 - Specificity | Y-axis: Sensitivity",
                theme = theme(
                  plot.title = element_text(size = by_scenario_grid_separate_sROC_settings$base_size + 2, face = "bold"),
                  plot.caption = element_text(size = 14, hjust = 0.5)
                )
              )
            
            plots$grid_separate[[paste0("grid_separate_scenario_", s)]] <- grid_plot
          }
        }

        
        ##
        ## ---- At the END - rename outputs for BASELINE (i.e., NON-MR) case (single scenario!):
        ##
        ##
        if (n_scenarios == 1) {
              
              # Create cleaner names for non-MR outputs
              clean_plots <- list()
              
              # Main plot with all tests overlaid
              if (!is.null(plots$by_scenario)) {
                clean_plots$all_tests <- plots$by_scenario[[1]]
              }
              
              # Individual test panels
              if (!is.null(plots$separate)) {
                clean_plots$by_test <- plots$separate
              }
              
              # Grid of tests (faceted)
              if (!is.null(plots$grid)) {
                clean_plots$faceted <- plots$grid[[1]]
              }
              
              # Separate panels with CI/PI regions
              if (!is.null(plots$grid_separate)) {
                clean_plots$panels_with_regions <- tail(plots$grid_separate, 1)[[1]]
              }
              
              return(clean_plots)
          
        }
        
        return(plots)
      
}








#' Create Se/Sp vs threshold plots with different layouts
#' @export
create_threshold_plots_flexible <- function(
                                            combined_data,
                                            ##
                                            test_names = NULL,
                                            relevant_thresholds = NULL,
                                            ##
                                            plot_type = c("Se", "Sp", "both"),
                                            layout = c("separate", 
                                                       "by_test", "by_scenario",
                                                       "grid", "grid_separate"),
                                            ##
                                            separate_threshold_settings = NULL,
                                            ##
                                            by_test_threshold_settings = NULL,
                                            by_scenario_threshold_settings = NULL,
                                            ##
                                            by_test_grid_threshold_settings = NULL,
                                            by_scenario_grid_threshold_settings = NULL,
                                            ##
                                            by_test_grid_separate_threshold_settings = NULL,
                                            by_scenario_grid_separate_threshold_settings = NULL,
                                            ##
                                            # Aesthetic mapping options (same as sROC)
                                            aesthetic_mapping = NULL,
                                            ##
                                            show_actual_thresholds = FALSE  # Whether to show actual threshold values
) {
  
        require(ggplot2)
        require(dplyr)
        require(patchwork)
        ##
        scenario_levels <- unique(combined_data$scenario)
        ##
        # Add test names if provided
        if (!is.null(test_names)) {
          combined_data <- combined_data %>%
            mutate(test_name = test_names[test])
        } else {
          combined_data <- combined_data %>%
            mutate(test_name = paste0("Test ", test))
        }
        
        # AFTER creating test_name, set factor levels based on relevant_thresholds order
        if (!is.null(relevant_thresholds) && !is.null(test_names)) {
          # Get the order from relevant_thresholds
          desired_order <- names(relevant_thresholds)
          # Apply as factor levels
          combined_data <- combined_data %>%
            mutate(test_name = factor(test_name, levels = desired_order))
        } else if (!is.null(test_names)) {
          # No relevant_thresholds, just use test_names as-is
          combined_data <- combined_data %>%
            mutate(test_name = factor(test_name, levels = test_names))
        } else {
          combined_data <- combined_data %>%
            mutate(test_name = factor(test_name))
        }
        ##
        # # Add test names if provided
        # if (!is.null(test_names)) {
        #   combined_data <- combined_data %>%
        #     mutate(test_name = test_names[test])
        # } else {
        #   combined_data <- combined_data %>%
        #     mutate(test_name = paste0("Test ", test))
        # }
        # ##
        # # AFTER creating test_name, set factor levels based on relevant_thresholds order
        # if (!is.null(relevant_thresholds) && !is.null(test_names)) {
        #   
        #       # Get the order from relevant_thresholds
        #       desired_order <- names(relevant_thresholds)
        #       # Apply as factor levels
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(test_name, levels = desired_order))
        #   
        # } else if (!is.null(test_names)) {
        #   
        #       # No relevant_thresholds, just use test_names as-is
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(test_name, levels = test_names))
        #   
        # } else {
        #   
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(test_name))
        #   
        # }
        # # Establish test order from relevant_thresholds if provided
        # if (!is.null(relevant_thresholds)) {
        #   
        #       test_order <- names(relevant_thresholds)
        #       
        #       # If test_names provided, reorder it to match relevant_thresholds
        #       if (!is.null(test_names)) {
        #         test_names <- test_names[match(test_order, names(test_names))]
        #       } else {
        #         # If no test_names, use the names from relevant_thresholds
        #         test_names <- test_order
        #         names(test_names) <- test_order
        #       }
        #   
        # }
        # ##
        # # Add test names with CORRECT factor levels
        # if (!is.null(test_names)) {
        #   
        #       # Create the correct level order
        #       test_name_levels <- test_names
        #       
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(test_names[test], levels = test_name_levels))
        #   
        # } else {
        #   
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(paste0("Test ", test)))
        #   
        # }
        # ##
        # # CRITICAL: Also ensure 'test' column respects the order
        # if (!is.null(test_names)) {
        #   combined_data <- combined_data %>%
        #     mutate(test = factor(test, levels = 1:length(test_names)))
        # }
        ##
        plot_type <- match.arg(plot_type, several.ok = TRUE)
        # 
        # # # Use the existing extract function
        # combined_data <- extract_threshold_data_scenarios(all_scenario_results)
        # 
        # Add test names if provided
        # if (!is.null(test_names)) {
        #   combined_data <- combined_data %>%
        #     mutate(test_name = factor(test_names[test], levels = test_names))
        # } else {
        #   combined_data <- combined_data %>%
        #     mutate(test_name = factor(paste0("Test ", test)))
        # }
        # # Add test names if provided, preserving order from relevant_thresholds
        # if (!is.null(test_names)) {
        #   
        #       # Create test_name column using numeric test index
        #       combined_data <- combined_data %>%
        #         mutate(test_name = test_names[test])
        #       
        #       # If relevant_thresholds provided, use its order for factor levels
        #       if (!is.null(relevant_thresholds)) {
        #         test_name_order <- names(relevant_thresholds)
        #       } else {
        #         test_name_order <- test_names
        #       }
        #       
        #       # Set factor levels in the correct order
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(test_name, levels = test_name_order))
        #   
        # } else {
        #   
        #       combined_data <- combined_data %>%
        #         mutate(test_name = factor(paste0("Test ", test)))
        #   
        # }
        # # THRESHOLD FILTERING (same as sROC)
        # if (!is.null(relevant_thresholds)) {
        #   if (is.list(relevant_thresholds)) {
        #     # Filter by test-specific thresholds
        #     combined_data <- combined_data %>%
        #       rowwise() %>%
        #       filter(threshold_actual %in% relevant_thresholds[[test_name]]) %>%
        #       ungroup()
        #   } else {
        #     # Filter all tests by same thresholds
        #     combined_data <- combined_data %>%
        #       filter(threshold_actual %in% relevant_thresholds)
        #   }
        # }
        ##
        # Apply aesthetic mapping if provided
        if (!is.null(aesthetic_mapping) && !is.null(aesthetic_mapping$parse_fn)) {
          combined_data <- aesthetic_mapping$parse_fn(combined_data)
        }
        ##
        # Get unique values
        n_scenarios <- length(unique(combined_data$scenario))
        n_tests <- length(unique(combined_data$test))
        scenario_names <- unique(combined_data$scenario)
        ##
        ##
        ##
        # separate_threshold_settings
        if (is.null(separate_threshold_settings)) separate_threshold_settings <- list()
        if (is.null(separate_threshold_settings$base_size)) separate_threshold_settings$base_size <- 14
        if (is.null(separate_threshold_settings$point_size)) separate_threshold_settings$point_size <- 2
        if (is.null(separate_threshold_settings$line_size)) separate_threshold_settings$line_size <- 0.8
        if (is.null(separate_threshold_settings$rows_in_legend)) separate_threshold_settings$rows_in_legend <- 3
        if (is.null(separate_threshold_settings$legend_position)) separate_threshold_settings$legend_position <- "none"
        if (is.null(separate_threshold_settings$add_regions)) separate_threshold_settings$add_regions <- TRUE
        if (is.null(separate_threshold_settings$show_conf)) separate_threshold_settings$show_conf <- TRUE
        if (is.null(separate_threshold_settings$show_pred)) separate_threshold_settings$show_pred <- TRUE
        if (is.null(separate_threshold_settings$conf_region_alpha)) separate_threshold_settings$conf_region_alpha <- 0.30
        if (is.null(separate_threshold_settings$pred_region_alpha)) separate_threshold_settings$pred_region_alpha <- 0.15
        if (is.null(separate_threshold_settings$conf_region_colour_Se)) separate_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(separate_threshold_settings$conf_region_colour_Sp)) separate_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(separate_threshold_settings$pred_region_colour_Se)) separate_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(separate_threshold_settings$pred_region_colour_Sp)) separate_threshold_settings$pred_region_colour_Sp <- "green"
        
        # by_test_threshold_settings
        if (is.null(by_test_threshold_settings)) by_test_threshold_settings <- list()
        if (is.null(by_test_threshold_settings$base_size)) by_test_threshold_settings$base_size <- 14
        if (is.null(by_test_threshold_settings$point_size)) by_test_threshold_settings$point_size <- 2
        if (is.null(by_test_threshold_settings$line_size)) by_test_threshold_settings$line_size <- 0.8
        if (is.null(by_test_threshold_settings$rows_in_legend)) by_test_threshold_settings$rows_in_legend <- 3
        if (is.null(by_test_threshold_settings$legend_position)) by_test_threshold_settings$legend_position <- "bottom"
        if (is.null(by_test_threshold_settings$add_regions)) by_test_threshold_settings$add_regions <- TRUE
        if (is.null(by_test_threshold_settings$show_conf)) by_test_threshold_settings$show_conf <- TRUE
        if (is.null(by_test_threshold_settings$show_pred)) by_test_threshold_settings$show_pred <- FALSE
        if (is.null(by_test_threshold_settings$conf_region_alpha)) by_test_threshold_settings$conf_region_alpha <- 0.30
        if (is.null(by_test_threshold_settings$pred_region_alpha)) by_test_threshold_settings$pred_region_alpha <- 0.15
        if (is.null(by_test_threshold_settings$conf_region_colour_Se)) by_test_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_test_threshold_settings$conf_region_colour_Sp)) by_test_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_test_threshold_settings$pred_region_colour_Se)) by_test_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_test_threshold_settings$pred_region_colour_Sp)) by_test_threshold_settings$pred_region_colour_Sp <- "green"
        
        # by_scenario_threshold_settings
        if (is.null(by_scenario_threshold_settings)) by_scenario_threshold_settings <- list()
        if (is.null(by_scenario_threshold_settings$base_size)) by_scenario_threshold_settings$base_size <- 14
        if (is.null(by_scenario_threshold_settings$point_size)) by_scenario_threshold_settings$point_size <- 2
        if (is.null(by_scenario_threshold_settings$line_size)) by_scenario_threshold_settings$line_size <- 0.8
        if (is.null(by_scenario_threshold_settings$rows_in_legend)) by_scenario_threshold_settings$rows_in_legend <- 2
        if (is.null(by_scenario_threshold_settings$legend_position)) by_scenario_threshold_settings$legend_position <- "bottom"
        if (is.null(by_scenario_threshold_settings$add_regions)) by_scenario_threshold_settings$add_regions <- TRUE
        if (is.null(by_scenario_threshold_settings$show_conf)) by_scenario_threshold_settings$show_conf <- TRUE
        if (is.null(by_scenario_threshold_settings$show_pred)) by_scenario_threshold_settings$show_pred <- FALSE
        if (is.null(by_scenario_threshold_settings$conf_region_alpha)) by_scenario_threshold_settings$conf_region_alpha <- 0.30
        if (is.null(by_scenario_threshold_settings$pred_region_alpha)) by_scenario_threshold_settings$pred_region_alpha <- 0.15
        if (is.null(by_scenario_threshold_settings$conf_region_colour_Se)) by_scenario_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_scenario_threshold_settings$conf_region_colour_Sp)) by_scenario_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_scenario_threshold_settings$pred_region_colour_Se)) by_scenario_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_scenario_threshold_settings$pred_region_colour_Sp)) by_scenario_threshold_settings$pred_region_colour_Sp <- "green"
        
        # by_test_grid_threshold_settings
        if (is.null(by_test_grid_threshold_settings)) by_test_grid_threshold_settings <- list()
        if (is.null(by_test_grid_threshold_settings$base_size)) by_test_grid_threshold_settings$base_size <- 14
        if (is.null(by_test_grid_threshold_settings$point_size)) by_test_grid_threshold_settings$point_size <- 2
        if (is.null(by_test_grid_threshold_settings$line_size)) by_test_grid_threshold_settings$line_size <- 0.8
        if (is.null(by_test_grid_threshold_settings$rows_in_legend)) by_test_grid_threshold_settings$rows_in_legend <- 3
        if (is.null(by_test_grid_threshold_settings$legend_position)) by_test_grid_threshold_settings$legend_position <- "bottom"
        if (is.null(by_test_grid_threshold_settings$add_regions)) by_test_grid_threshold_settings$add_regions <- FALSE
        if (is.null(by_test_grid_threshold_settings$show_conf)) by_test_grid_threshold_settings$show_conf <- FALSE
        if (is.null(by_test_grid_threshold_settings$show_pred)) by_test_grid_threshold_settings$show_pred <- FALSE
        if (is.null(by_test_grid_threshold_settings$conf_region_alpha)) by_test_grid_threshold_settings$conf_region_alpha <- NULL
        if (is.null(by_test_grid_threshold_settings$pred_region_alpha)) by_test_grid_threshold_settings$pred_region_alpha <- NULL
        if (is.null(by_test_grid_threshold_settings$conf_region_colour_Se)) by_test_grid_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_test_grid_threshold_settings$conf_region_colour_Sp)) by_test_grid_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_test_grid_threshold_settings$pred_region_colour_Se)) by_test_grid_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_test_grid_threshold_settings$pred_region_colour_Sp)) by_test_grid_threshold_settings$pred_region_colour_Sp <- "green"
        if (is.null(by_test_grid_threshold_settings$facet_label_width)) by_test_grid_threshold_settings$facet_label_width <- 25
        
        # by_scenario_grid_threshold_settings
        if (is.null(by_scenario_grid_threshold_settings)) by_scenario_grid_threshold_settings <- list()
        if (is.null(by_scenario_grid_threshold_settings$base_size)) by_scenario_grid_threshold_settings$base_size <- 14
        if (is.null(by_scenario_grid_threshold_settings$point_size)) by_scenario_grid_threshold_settings$point_size <- 2
        if (is.null(by_scenario_grid_threshold_settings$line_size)) by_scenario_grid_threshold_settings$line_size <- 0.8
        if (is.null(by_scenario_grid_threshold_settings$rows_in_legend)) by_scenario_grid_threshold_settings$rows_in_legend <- 1
        if (is.null(by_scenario_grid_threshold_settings$legend_position)) by_scenario_grid_threshold_settings$legend_position <- "bottom"
        if (is.null(by_scenario_grid_threshold_settings$add_regions)) by_scenario_grid_threshold_settings$add_regions <- FALSE
        if (is.null(by_scenario_grid_threshold_settings$show_conf)) by_scenario_grid_threshold_settings$show_conf <- FALSE
        if (is.null(by_scenario_grid_threshold_settings$show_pred)) by_scenario_grid_threshold_settings$show_pred <- FALSE
        if (is.null(by_scenario_grid_threshold_settings$conf_region_alpha)) by_scenario_grid_threshold_settings$conf_region_alpha <- NULL
        if (is.null(by_scenario_grid_threshold_settings$pred_region_alpha)) by_scenario_grid_threshold_settings$pred_region_alpha <- NULL
        if (is.null(by_scenario_grid_threshold_settings$conf_region_colour_Se)) by_scenario_grid_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_scenario_grid_threshold_settings$conf_region_colour_Sp)) by_scenario_grid_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_scenario_grid_threshold_settings$pred_region_colour_Se)) by_scenario_grid_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_scenario_grid_threshold_settings$pred_region_colour_Sp)) by_scenario_grid_threshold_settings$pred_region_colour_Sp <- "green"
        if (is.null(by_scenario_grid_threshold_settings$facet_label_width)) by_scenario_grid_threshold_settings$facet_label_width <- 25
        
        # by_test_grid_separate_threshold_settings
        if (is.null(by_test_grid_separate_threshold_settings)) by_test_grid_separate_threshold_settings <- list()
        if (is.null(by_test_grid_separate_threshold_settings$base_size)) by_test_grid_separate_threshold_settings$base_size <- 14
        if (is.null(by_test_grid_separate_threshold_settings$point_size)) by_test_grid_separate_threshold_settings$point_size <- 2
        if (is.null(by_test_grid_separate_threshold_settings$line_size)) by_test_grid_separate_threshold_settings$line_size <- 0.8
        if (is.null(by_test_grid_separate_threshold_settings$rows_in_legend)) by_test_grid_separate_threshold_settings$rows_in_legend <- 3
        if (is.null(by_test_grid_separate_threshold_settings$legend_position)) by_test_grid_separate_threshold_settings$legend_position <- "none"
        if (is.null(by_test_grid_separate_threshold_settings$add_bands)) by_test_grid_separate_threshold_settings$add_bands <- TRUE
        if (is.null(by_test_grid_separate_threshold_settings$show_conf)) by_test_grid_separate_threshold_settings$show_conf <- TRUE
        if (is.null(by_test_grid_separate_threshold_settings$show_pred)) by_test_grid_separate_threshold_settings$show_pred <- TRUE
        if (is.null(by_test_grid_separate_threshold_settings$conf_region_alpha)) by_test_grid_separate_threshold_settings$conf_region_alpha <- 0.30
        if (is.null(by_test_grid_separate_threshold_settings$pred_region_alpha)) by_test_grid_separate_threshold_settings$pred_region_alpha <- 0.15
        if (is.null(by_test_grid_separate_threshold_settings$conf_region_colour_Se)) by_test_grid_separate_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_test_grid_separate_threshold_settings$conf_region_colour_Sp)) by_test_grid_separate_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_test_grid_separate_threshold_settings$pred_region_colour_Se)) by_test_grid_separate_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_test_grid_separate_threshold_settings$pred_region_colour_Sp)) by_test_grid_separate_threshold_settings$pred_region_colour_Sp <- "green"
        
        # by_scenario_grid_separate_threshold_settings
        if (is.null(by_scenario_grid_separate_threshold_settings)) by_scenario_grid_separate_threshold_settings <- list()
        if (is.null(by_scenario_grid_separate_threshold_settings$base_size)) by_scenario_grid_separate_threshold_settings$base_size <- 14
        if (is.null(by_scenario_grid_separate_threshold_settings$point_size)) by_scenario_grid_separate_threshold_settings$point_size <- 2
        if (is.null(by_scenario_grid_separate_threshold_settings$line_size)) by_scenario_grid_separate_threshold_settings$line_size <- 0.8
        if (is.null(by_scenario_grid_separate_threshold_settings$rows_in_legend)) by_scenario_grid_separate_threshold_settings$rows_in_legend <- 1
        if (is.null(by_scenario_grid_separate_threshold_settings$legend_position)) by_scenario_grid_separate_threshold_settings$legend_position <- "none"
        if (is.null(by_scenario_grid_separate_threshold_settings$add_bands)) by_scenario_grid_separate_threshold_settings$add_bands <- TRUE
        if (is.null(by_scenario_grid_separate_threshold_settings$show_conf)) by_scenario_grid_separate_threshold_settings$show_conf <- TRUE
        if (is.null(by_scenario_grid_separate_threshold_settings$show_pred)) by_scenario_grid_separate_threshold_settings$show_pred <- TRUE
        if (is.null(by_scenario_grid_separate_threshold_settings$conf_region_alpha)) by_scenario_grid_separate_threshold_settings$conf_region_alpha <- 0.30
        if (is.null(by_scenario_grid_separate_threshold_settings$pred_region_alpha)) by_scenario_grid_separate_threshold_settings$pred_region_alpha <- 0.15
        if (is.null(by_scenario_grid_separate_threshold_settings$conf_region_colour_Se)) by_scenario_grid_separate_threshold_settings$conf_region_colour_Se <- "red"
        if (is.null(by_scenario_grid_separate_threshold_settings$conf_region_colour_Sp)) by_scenario_grid_separate_threshold_settings$conf_region_colour_Sp <- "green"
        if (is.null(by_scenario_grid_separate_threshold_settings$pred_region_colour_Se)) by_scenario_grid_separate_threshold_settings$pred_region_colour_Se <- "red"
        if (is.null(by_scenario_grid_separate_threshold_settings$pred_region_colour_Sp)) by_scenario_grid_separate_threshold_settings$pred_region_colour_Sp <- "green"
        
        
        # Determine x variable
        x_var <- ifelse(show_actual_thresholds, "threshold_actual", "threshold_norm")
        x_label <- ifelse(show_actual_thresholds, "Test threshold", "Normalized threshold (0-1)")
 
        
        plots <- list()
        plots$separate <- list()
        plots$by_test <- list()
        plots$by_test <- list()
        plots$by_scenario <- list()
        plots$grid_separate <- list()
        
        
        # p <- add_bands_to_threshold(p = p, 
        #                             data = scenario_data,
        #                             y_var = "Se", 
        #                             settings = by_test_grid_separate_threshold_settings)
        
        # Helper function to add bands (ribbons for threshold plots)
        add_bands_to_threshold <- function(p, 
                                           data, 
                                           y_var, 
                                           settings) {
          
                # Check for either add_bands or add_regions (for compatibility)
                if (!isTRUE(settings$add_bands) && !isTRUE(settings$add_regions)) return(p)
                
                if (y_var == "Se") { 
                  pred_color <- settings$pred_region_colour_Se
                  conf_color <- settings$conf_region_colour_Se
                } else if (y_var == "Sp") { 
                  pred_color <- settings$pred_region_colour_Sp
                  conf_color <- settings$conf_region_colour_Sp
                }
          
                # Determine x variable from the plot
                x_var <- as.character(p$mapping$x)[2]  # Extract the actual variable name
                
                if (isTRUE(settings$show_conf)) {
                  lower_var <- paste0(y_var, "_lower")
                  upper_var <- paste0(y_var, "_upper")
                  
                  # Check if we have one group or multiple
                  n_groups <- length(unique(data$test_name))
                  if (n_groups == 1) {
                    # Single group - use fixed color
                    p <- p + geom_ribbon(
                      data = data,
                      aes_string(x = x_var, ymin = lower_var, ymax = upper_var),
                      fill = conf_color,
                      alpha = settings$conf_region_alpha
                    )
                  } else {
                    # Multiple groups - inherit fill
                    p <- p + geom_ribbon(
                      data = data,
                      aes_string(x = x_var, ymin = lower_var, ymax = upper_var, 
                                 fill = "test_name"),
                      alpha = settings$conf_region_alpha
                    )
                  }
                }
                
                if (isTRUE(settings$show_pred)) {
                  pred_lower <- paste0(y_var, "_pred_lower")
                  pred_upper <- paste0(y_var, "_pred_upper")
                  
                  n_groups <- length(unique(data$test_name))
                  if (n_groups == 1) {
                    p <- p + geom_ribbon(
                      data = data,
                      aes_string(x = x_var, ymin = pred_lower, ymax = pred_upper),
                      fill = pred_color,
                      alpha = settings$pred_region_alpha
                    )
                  } else {
                    p <- p + geom_ribbon(
                      data = data,
                      aes_string(x = x_var, ymin = pred_lower, ymax = pred_upper,
                                 fill = "test_name"),
                      alpha = settings$pred_region_alpha
                    )
                  }
                }
                
                return(p)
          
        }
        
        
        
        # Modified version for combined Se/Sp plots
        add_bands_to_threshold_combined <- function(p, 
                                                    data, 
                                                    settings) {
          
                # Check for either add_bands or add_regions (for compatibility)
                if (!isTRUE(settings$add_bands) && !isTRUE(settings$add_regions)) return(p)
                
                # Determine x variable from the plot
                x_var <- as.character(p$mapping$x)[2]
                
                if (isTRUE(settings$show_conf)) {
                  # Add both Se and Sp bands with matching colors from the lines
                  p <- p + 
                    geom_ribbon(data = data,
                                aes_string(x = x_var, 
                                           ymin = "Se_lower", 
                                           ymax = "Se_upper"),
                                fill = settings$conf_region_colour_Se, 
                                alpha = settings$conf_region_alpha) +
                    geom_ribbon(data = data,
                                aes_string(x = x_var, 
                                           ymin = "Sp_lower", 
                                           ymax = "Sp_upper"),
                                fill = settings$conf_region_colour_Sp, 
                                alpha = settings$conf_region_alpha)
                }
                
                if (isTRUE(settings$show_pred)) {
                  p <- p + 
                    geom_ribbon(data = data,
                                aes_string(x = x_var, 
                                           ymin = "Se_pred_lower", 
                                           ymax = "Se_pred_upper"),
                                fill = settings$pred_region_colour_Se, 
                                alpha = settings$pred_region_alpha) +
                    geom_ribbon(data = data,
                                aes_string(x = x_var, 
                                           ymin = "Sp_pred_lower",
                                           ymax = "Sp_pred_upper"),
                                fill = settings$pred_region_colour_Sp,
                                alpha = settings$pred_region_alpha)
                }
                
                return(p)
          
        }
        
        
        
        # 1. Separate plots for each test-scenario combination
        if ("separate" %in% layout) {
          
              for (s in scenario_names) {
                
                # Loop through tests in the correct order
                test_order <- if (!is.null(relevant_thresholds)) {
                  seq_along(names(relevant_thresholds))
                } else {
                  unique(combined_data$test)
                }
                
                for (t in test_order) {
                # for (t in unique(combined_data$test)) {
                  
                      plot_data <- combined_data %>%
                        filter(scenario == s, test == t)
                      
                      test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                      
                      # Create Se plot
                      if ("Se" %in% plot_type || "both" %in% plot_type) {
                        p_se <- create_threshold_plot_base(
                          data = plot_data,
                          y_var = "Se_median",
                          x_var = x_var,
                          title = paste0("Sensitivity - ", test_name, " - ", s),
                          x_label = x_label,
                          y_label = "Sensitivity",
                          base_size = separate_threshold_settings$base_size,
                          point_size = separate_threshold_settings$point_size,
                          line_size = separate_threshold_settings$line_size,
                          hide_x_axis = !show_actual_thresholds
                        ) + theme(legend.position = separate_threshold_settings$legend_position)
                        
                        p_se <- add_bands_to_threshold(p_se, plot_data, "Se", separate_threshold_settings)
                        plots$separate[[paste0("Se_", s, "_", test_name)]] <- p_se
                      }
                      
                      # Create Sp plot
                      if ("Sp" %in% plot_type || "both" %in% plot_type) {
                        p_sp <- create_threshold_plot_base(
                          data = plot_data,
                          y_var = "Sp_median",
                          x_var = x_var,
                          title = paste0("Specificity - ", test_name, " - ", s),
                          x_label = x_label,
                          y_label = "Specificity",
                          base_size = separate_threshold_settings$base_size,
                          point_size = separate_threshold_settings$point_size,
                          line_size = separate_threshold_settings$line_size,
                          hide_x_axis = !show_actual_thresholds
                        ) + theme(legend.position = separate_threshold_settings$legend_position)
                        
                        p_sp <- add_bands_to_threshold(p_sp, plot_data, "Sp", separate_threshold_settings)
                        plots$separate[[paste0("Sp_", s, "_", test_name)]] <- p_sp
                      }
                      
                      # Combined Se/Sp plot (stacked vertically like plot_7)
                      if ("both" %in% plot_type) {
                        p_both <- plots$separate[[paste0("Se_", s, "_", test_name)]] + 
                          theme(legend.position = "none") + xlab("") +
                          plots$separate[[paste0("Sp_", s, "_", test_name)]] +
                          plot_layout(ncol = 1)
                        
                        plots$separate[[paste0("both_", s, "_", test_name)]] <- p_both
                      }
                      
                      # In the separate layout section, after the existing "both" code:
                      if ("both" %in% plot_type) {
                        # Add this new joined version:
                        p_both_joined <- ggplot(plot_data, aes_string(x = x_var)) +
                          # Se line and points (solid)
                          geom_line(aes_string(y = "Se_median"), 
                                    color = separate_threshold_settings$conf_region_colour_Se,
                                    size = separate_threshold_settings$line_size) +
                          geom_point(aes_string(y = "Se_median"), 
                                     color = separate_threshold_settings$conf_region_colour_Se,
                                     size = separate_threshold_settings$point_size) +
                          # Sp line and points (dashed)
                          geom_line(aes_string(y = "Sp_median"), 
                                    color = separate_threshold_settings$conf_region_colour_Sp,
                                    linetype = "dashed",
                                    size = separate_threshold_settings$line_size) +
                          geom_point(aes_string(y = "Sp_median"), 
                                     color = separate_threshold_settings$conf_region_colour_Sp,
                                     shape = 17,
                                     size = separate_threshold_settings$point_size) +
                          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                             limits = c(0, 1)) +
                          labs(x = x_label, 
                               y = "Se (solid) / Sp (dashed)",
                               title = paste0(test_name, " - ", s)) +
                          theme_bw(base_size = separate_threshold_settings$base_size) +
                          theme(legend.position = separate_threshold_settings$legend_position)
                        
                        if (!show_actual_thresholds) {
                          p_both_joined <- p_both_joined + 
                            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                        }
                        
                        # Add bands using combined function
                        p_both_joined <- add_bands_to_threshold_combined(p_both_joined, plot_data, separate_threshold_settings)
                        
                        plots$separate[[paste0("both_joined_", s, "_", test_name)]] <- p_both_joined
                      }
                  
                }
                
              }
          
        }
        
        # 2. One plot per test, all scenarios overlaid
        if ("by_test" %in% layout) {
          
          # Loop through tests in the correct order
          test_order <- if (!is.null(relevant_thresholds)) {
            seq_along(names(relevant_thresholds))
          } else {
            unique(combined_data$test)
          }
          
          for (t in test_order) {
          # for (t in unique(combined_data$test)) {
            
            plot_data <- combined_data %>%
              filter(test == t) %>%
              mutate(scenario = factor(scenario, levels = scenario_levels))
            
            test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
            
            # Check for aesthetic mapping
            if (!is.null(aesthetic_mapping) && 
                !is.null(aesthetic_mapping$color_var) && 
                !is.null(aesthetic_mapping$linetype_var)) {
              
              # Se plot with custom aesthetics
              if ("Se" %in% plot_type || "both" %in% plot_type) {
                      
                      p_se <- ggplot(plot_data, 
                                     aes_string(x = x_var, y = "Se_median",
                                                color = aesthetic_mapping$color_var,
                                                linetype = aesthetic_mapping$linetype_var,
                                                group = "scenario")) +
                        geom_line(size = by_test_threshold_settings$line_size) +
                        geom_point(aes_string(shape = aesthetic_mapping$shape_var), 
                                   size = by_test_threshold_settings$point_size) +
                        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
                        labs(x = x_label, 
                             y = "Sensitivity", 
                             title = paste0("Sensitivity vs Threshold - ", test_name)) +  # ADD THIS TITLE
                        theme_bw(base_size = by_test_threshold_settings$base_size) +
                        theme(legend.position = by_test_threshold_settings$legend_position)
                      
                      # Apply custom scales
                      if (!is.null(aesthetic_mapping$color_scale)) p_se <- p_se + aesthetic_mapping$color_scale
                      if (!is.null(aesthetic_mapping$linetype_scale)) p_se <- p_se + aesthetic_mapping$linetype_scale
                      if (!is.null(aesthetic_mapping$shape_scale)) p_se <- p_se + aesthetic_mapping$shape_scale
                      
                      if (!show_actual_thresholds) {
                        p_se <- p_se + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                      }
                      
                      plots$by_test[[paste0("Se_by_test_", test_name)]] <- p_se
                
              }
              
              # Sp plot with custom aesthetics (similar structure)
              if ("Sp" %in% plot_type || "both" %in% plot_type) {
                
                      p_sp <- ggplot(plot_data, 
                                     aes_string(x = x_var, y = "Sp_median",
                                                color = aesthetic_mapping$color_var,
                                                linetype = aesthetic_mapping$linetype_var,
                                                group = "scenario")) +
                        geom_line(size = by_test_threshold_settings$line_size) +
                        geom_point(aes_string(shape = aesthetic_mapping$shape_var), 
                                   size = by_test_threshold_settings$point_size) +
                        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
                        labs(x = x_label, 
                             y = "Specificity", 
                             title = paste0("Specificity vs Threshold - ", test_name)) +  # ADD THIS TITLE
                        theme_bw(base_size = by_test_threshold_settings$base_size) +
                        theme(legend.position = by_test_threshold_settings$legend_position)
                      
                      # Apply custom scales
                      if (!is.null(aesthetic_mapping$color_scale)) p_sp <- p_sp + aesthetic_mapping$color_scale
                      if (!is.null(aesthetic_mapping$linetype_scale)) p_sp <- p_sp + aesthetic_mapping$linetype_scale
                      if (!is.null(aesthetic_mapping$shape_scale)) p_sp <- p_sp + aesthetic_mapping$shape_scale
                      
                      if (!show_actual_thresholds) {
                        p_sp <- p_sp + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                      }
                      
                      plots$by_test[[paste0("Sp_by_test_", test_name)]] <- p_sp
                
              }
              
            } else {
              
                  # Use default aesthetics
                  if ("Se" %in% plot_type || "both" %in% plot_type) {
                    p_se <- create_threshold_plot_base(
                      data = plot_data,
                      y_var = "Se_median",
                      x_var = x_var,
                      color_by = "scenario",
                      group_by = "scenario",
                      title = paste0("Sensitivity vs Threshold - ", test_name),
                      subtitle = "Colored by meta-regression scenario",
                      x_label = x_label,
                      y_label = "Sensitivity",
                      base_size = by_test_threshold_settings$base_size,
                      point_size = by_test_threshold_settings$point_size,
                      line_size = by_test_threshold_settings$line_size,
                      hide_x_axis = !show_actual_thresholds
                    ) + 
                      theme(legend.position = by_test_threshold_settings$legend_position) +
                      guides(color = guide_legend(nrow = by_test_threshold_settings$rows_in_legend))
                    
                    # No bands for by_test (too messy with multiple scenarios)
                    plots$by_test[[paste0("Se_by_test_", test_name)]] <- p_se
                  }
                  
                  if ("Sp" %in% plot_type || "both" %in% plot_type) {
                    p_sp <- create_threshold_plot_base(
                      data = plot_data,
                      y_var = "Sp_median",
                      x_var = x_var,
                      color_by = "scenario",
                      group_by = "scenario",
                      title = paste0("Specificity vs Threshold - ", test_name),
                      subtitle = "Colored by meta-regression scenario",
                      x_label = x_label,
                      y_label = "Specificity",
                      base_size = by_test_threshold_settings$base_size,
                      point_size = by_test_threshold_settings$point_size,
                      line_size = by_test_threshold_settings$line_size,
                      hide_x_axis = !show_actual_thresholds
                    ) + 
                      theme(legend.position = by_test_threshold_settings$legend_position) +
                      guides(color = guide_legend(nrow = by_test_threshold_settings$rows_in_legend))
                    
                    plots$by_test[[paste0("Sp_by_test_", test_name)]] <- p_sp
                  }
              
            }
            
              # Combined plot
              if ("both" %in% plot_type && 
                  exists(paste0("Se_by_test_", test_name), where = plots$by_test) &&
                  exists(paste0("Sp_by_test_", test_name), where = plots$by_test)) {
                
                    p_both <- plots$by_test[[paste0("Se_by_test_", test_name)]] + 
                      theme(legend.position = "none") + xlab("") +
                      plots$by_test[[paste0("Sp_by_test_", test_name)]] +
                      plot_layout(ncol = 1, guides = "collect") &
                      theme(legend.position = "bottom")
                    
                    plots$by_test[[paste0("both_by_test_", test_name)]] <- p_both
                    
              }
            
            # Add this for both_joined with custom aesthetics:
            if ("both" %in% plot_type) {
              p_both_joined <- ggplot(plot_data, 
                                      aes_string(x = x_var, 
                                                 color = aesthetic_mapping$color_var,
                                                 linetype = aesthetic_mapping$linetype_var,
                                                 group = "scenario")) +
                # Se lines (solid)
                geom_line(aes_string(y = "Se_median"), 
                          size = by_test_threshold_settings$line_size) +
                geom_point(aes_string(y = "Se_median", shape = aesthetic_mapping$shape_var), 
                           size = by_test_threshold_settings$point_size) +
                # Sp lines (dashed - but override linetype)
                geom_line(aes_string(y = "Sp_median"), 
                          linetype = "dashed",
                          size = by_test_threshold_settings$line_size) +
                geom_point(aes_string(y = "Sp_median", shape = aesthetic_mapping$shape_var), 
                           shape = 17,  # Triangle for Sp
                           size = by_test_threshold_settings$point_size) +
                scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                   limits = c(0, 1)) +
                labs(x = x_label, 
                     y = "Se (solid) / Sp (dashed)",
                     title = paste0("Se & Sp vs Threshold - ", test_name),
                     subtitle = "Se shown as solid lines/circles, Sp as dashed lines/triangles") +
                theme_bw(base_size = by_test_threshold_settings$base_size) +
                theme(legend.position = by_test_threshold_settings$legend_position)
              
              # Apply custom scales
              if (!is.null(aesthetic_mapping$color_scale)) {
                p_both_joined <- p_both_joined + aesthetic_mapping$color_scale
              }
              if (!is.null(aesthetic_mapping$linetype_scale)) {
                p_both_joined <- p_both_joined + aesthetic_mapping$linetype_scale
              }
              if (!is.null(aesthetic_mapping$shape_scale)) {
                p_both_joined <- p_both_joined + aesthetic_mapping$shape_scale
              }
              
              if (!show_actual_thresholds) {
                p_both_joined <- p_both_joined + 
                  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
              }
              
              plots$by_test[[paste0("both_joined_by_test_", test_name)]] <- p_both_joined
            }
            
            
          }
        }
        
        # 3. One plot per scenario, all tests overlaid
        if ("by_scenario" %in% layout) {
          
          for (s in scenario_names) {
            
                plot_data <- combined_data %>%
                  filter(scenario == s)
                
                if ("Se" %in% plot_type || "both" %in% plot_type) {
                  
                      p_se <- create_threshold_plot_base(
                        data = plot_data,
                        y_var = "Se_median",
                        x_var = x_var,
                        color_by = "test_name",
                        group_by = "test_name",
                        title = paste0("Sensitivity - ", s),
                        subtitle = "Colored by test",
                        x_label = x_label,
                        y_label = "Sensitivity",
                        base_size = by_scenario_threshold_settings$base_size,
                        point_size = by_scenario_threshold_settings$point_size,
                        line_size = by_scenario_threshold_settings$line_size,
                        hide_x_axis = !show_actual_thresholds
                      ) + 
                        theme(legend.position = by_scenario_threshold_settings$legend_position) +
                        guides(color = guide_legend(nrow = by_scenario_threshold_settings$rows_in_legend))
                      
                      p_se <- add_bands_to_threshold(p_se, plot_data, "Se", by_scenario_threshold_settings)
                      plots$by_scenario[[paste0("Se_by_scenario_", s)]] <- p_se
                  
                }
                
                if ("Sp" %in% plot_type || "both" %in% plot_type) {
                  
                      p_sp <- create_threshold_plot_base(
                        data = plot_data,
                        y_var = "Sp_median",
                        x_var = x_var,
                        color_by = "test_name",
                        group_by = "test_name",
                        title = paste0("Specificity - ", s),
                        subtitle = "Colored by test",
                        x_label = x_label,
                        y_label = "Specificity",
                        base_size = by_scenario_threshold_settings$base_size,
                        point_size = by_scenario_threshold_settings$point_size,
                        line_size = by_scenario_threshold_settings$line_size,
                        hide_x_axis = !show_actual_thresholds
                      ) + 
                        theme(legend.position = by_scenario_threshold_settings$legend_position) +
                        guides(color = guide_legend(nrow = by_scenario_threshold_settings$rows_in_legend))
                      
                      p_sp <- add_bands_to_threshold(p_sp, plot_data, "Sp", by_scenario_threshold_settings)
                      plots$by_scenario[[paste0("Sp_by_scenario_", s)]] <- p_sp
                }
                
                # In by_scenario section, add this after the Sp plot:
                if ("both" %in% plot_type) {
                  
                    # Create joined version
                    p_both_joined <- ggplot(plot_data, aes_string(x = x_var, group = "test_name")) +
                      # Se lines (solid)
                      geom_line(aes_string(y = "Se_median", color = "test_name"), 
                                size = by_scenario_threshold_settings$line_size) +
                      geom_point(aes_string(y = "Se_median", color = "test_name"), 
                                 size = by_scenario_threshold_settings$point_size) +
                      # Sp lines (dashed)
                      geom_line(aes_string(y = "Sp_median", color = "test_name"), 
                                linetype = "dashed",
                                size = by_scenario_threshold_settings$line_size) +
                      geom_point(aes_string(y = "Sp_median", color = "test_name"), 
                                 shape = 17,
                                 size = by_scenario_threshold_settings$point_size) +
                      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                         limits = c(0, 1)) +
                      labs(x = x_label, 
                           y = "Se (solid) / Sp (dashed)",
                           title = paste0("Se & Sp - ", s),
                           subtitle = "Se shown as solid lines/circles, Sp as dashed lines/triangles") +
                      theme_bw(base_size = by_scenario_threshold_settings$base_size) +
                      theme(legend.position = by_scenario_threshold_settings$legend_position) +
                      guides(color = guide_legend(nrow = by_scenario_threshold_settings$rows_in_legend))
                    
                    if (!show_actual_thresholds) {
                      p_both_joined <- p_both_joined + 
                        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                    }
                    
                    plots$by_scenario[[paste0("both_joined_by_scenario_", s)]] <- p_both_joined
                  
                }
            
          }
          
        }
        
        combined_data <- combined_data %>%
          mutate(scenario = factor(scenario, levels = scenario_levels))
        
        # 4. Grid layout
        if ("grid" %in% layout) {
          
            # Grid with tests as panels, scenarios as colors
            if ("Se" %in% plot_type || "both" %in% plot_type) {
              
                  p_se_grid <- create_threshold_plot_base(
                    data = combined_data,
                    y_var = "Se_median",
                    x_var = x_var,
                    color_by = "scenario",
                    group_by = "scenario",
                    facet_by = "test_name",
                    title = "Sensitivity vs Threshold by Test",
                    subtitle = "Scenarios shown in different colors",
                    x_label = x_label,
                    y_label = "Sensitivity",
                    base_size = by_test_grid_threshold_settings$base_size,
                    point_size = by_test_grid_threshold_settings$point_size,
                    line_size = by_test_grid_threshold_settings$line_size,
                    hide_x_axis = !show_actual_thresholds
                  ) + 
                    theme(
                      legend.position = by_test_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_test_grid_threshold_settings$base_size, face = "bold"),
                      panel.spacing = unit(1, "lines")
                    ) + 
                    guides(color = guide_legend(nrow = by_test_grid_threshold_settings$rows_in_legend))
                  
                  plots$grid[["Se_grid_by_test"]] <- p_se_grid
              
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
                    subtitle = "Scenarios shown in different colors",
                    x_label = x_label,
                    y_label = "Specificity",
                    base_size = by_test_grid_threshold_settings$base_size,
                    point_size = by_test_grid_threshold_settings$point_size,
                    line_size = by_test_grid_threshold_settings$line_size,
                    hide_x_axis = !show_actual_thresholds
                  ) + 
                    theme(
                      legend.position = by_test_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_test_grid_threshold_settings$base_size, face = "bold"),
                      panel.spacing = unit(1, "lines")
                    ) + 
                    guides(color = guide_legend(nrow = by_test_grid_threshold_settings$rows_in_legend))
                  
                  plots$grid[["Sp_grid_by_test"]] <- p_sp_grid
              
            }
            
            # Grid with scenarios as panels, tests as colors
            if ("Se" %in% plot_type || "both" %in% plot_type) {
                  
                  p_se_grid_scenario <- create_threshold_plot_base(
                    data = combined_data,
                    y_var = "Se_median",
                    x_var = x_var,
                    color_by = "test_name",
                    group_by = "test_name",
                    facet_by = "scenario",
                    title = "Sensitivity vs Threshold by Scenario",
                    subtitle = "Tests shown in different colors",
                    x_label = x_label,
                    y_label = "Sensitivity",
                    base_size = by_scenario_grid_threshold_settings$base_size,
                    point_size = by_scenario_grid_threshold_settings$point_size,
                    line_size = by_scenario_grid_threshold_settings$line_size,
                    hide_x_axis = !show_actual_thresholds
                  ) + 
                    theme(
                      legend.position = by_scenario_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_scenario_grid_threshold_settings$base_size * 0.7, face = "bold"),
                      panel.spacing = unit(0.5, "lines")
                    ) + 
                    guides(color = guide_legend(nrow = by_scenario_grid_threshold_settings$rows_in_legend))
                  
                  plots$grid[["Se_grid_by_scenario"]] <- p_se_grid_scenario
              
            }
            
            if ("Sp" %in% plot_type || "both" %in% plot_type) {
              
                  p_sp_grid_scenario <- create_threshold_plot_base(
                    data = combined_data,
                    y_var = "Sp_median",
                    x_var = x_var,
                    color_by = "test_name",
                    group_by = "test_name",
                    facet_by = "scenario",
                    title = "Specificity vs Threshold by Scenario",
                    subtitle = "Tests shown in different colors",
                    x_label = x_label,
                    y_label = "Specificity",
                    base_size = by_scenario_grid_threshold_settings$base_size,
                    point_size = by_scenario_grid_threshold_settings$point_size,
                    line_size = by_scenario_grid_threshold_settings$line_size,
                    hide_x_axis = !show_actual_thresholds
                  ) + 
                    theme(
                      legend.position = by_scenario_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_scenario_grid_threshold_settings$base_size * 0.7, face = "bold"),
                      panel.spacing = unit(0.5, "lines")
                    ) + 
                    guides(color = guide_legend(nrow = by_scenario_grid_threshold_settings$rows_in_legend))
                  
                  plots$grid[["Sp_grid_by_scenario"]] <- p_sp_grid_scenario
              
            }
            
            # Combined Se/Sp grids (like plot_9 in original)
            if ("both" %in% plot_type) {
                
                  # Panel plot with Se and Sp side by side
                  if (exists("Se_grid_by_test", where = plots$grid) && exists("Sp_grid_by_test", where = plots$grid)) {
                    p_both_grid <- plots$grid[["Se_grid_by_test"]] + plots$grid[["Sp_grid_by_test"]] +
                      plot_layout(ncol = 2, guides = "collect") &
                      theme(legend.position = "bottom")
                    
                    plots$grid[["both_grid_by_test"]] <- p_both_grid
                  }
                
                  # Combined plot with Se and Sp on same axes - by scenario
                  both_joined_grid_by_scenario <- ggplot(combined_data, aes_string(x = x_var, group = "test_name")) +
                    # Se lines (solid)
                    geom_line(aes_string(y = "Se_median", color = "test_name"), 
                              size = by_scenario_grid_threshold_settings$line_size) +
                    geom_point(aes_string(y = "Se_median", color = "test_name"), 
                               size = by_scenario_grid_threshold_settings$point_size) +
                    # Sp lines (dashed)
                    geom_line(aes_string(y = "Sp_median", color = "test_name"), 
                              linetype = "dashed",
                              size = by_scenario_grid_threshold_settings$line_size) +
                    geom_point(aes_string(y = "Sp_median", color = "test_name"), 
                               shape = 17,
                               size = by_scenario_grid_threshold_settings$point_size) +
                    facet_wrap(~ scenario, scales = "free_x",
                               dir = "v", 
                               drop = FALSE,
                               labeller = label_wrap_gen(width = by_scenario_grid_threshold_settings$facet_label_width)) +
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    labs(x = x_label, 
                         y = "Se (solid) / Sp (dashed)",
                         title = "Sensitivity and Specificity vs Threshold by Scenario",
                         subtitle = "Se shown as solid lines/circles, Sp as dashed lines/triangles") +
                    theme_bw(base_size = by_scenario_grid_threshold_settings$base_size) +
                    theme(
                      legend.position = by_scenario_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_scenario_grid_threshold_settings$base_size * 0.7, face = "bold"),
                      panel.spacing = unit(0.5, "lines")
                    ) +
                    guides(color = guide_legend(nrow = by_scenario_grid_threshold_settings$rows_in_legend))
                  
                  if (!show_actual_thresholds) {
                    both_joined_grid_by_scenario <- both_joined_grid_by_scenario + 
                      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                  }
                  
                  plots$grid[["both_joined_grid_by_scenario"]] <- both_joined_grid_by_scenario
                
                  # Combined plot with Se and Sp on same axes - by test
                  both_joined_grid_by_test <- ggplot(combined_data, aes_string(x = x_var, group = "scenario")) +
                    # Se lines (solid)
                    geom_line(aes_string(y = "Se_median", color = "scenario"), 
                              size = by_test_grid_threshold_settings$line_size) +
                    geom_point(aes_string(y = "Se_median", color = "scenario"), 
                               size = by_test_grid_threshold_settings$point_size) +
                    # Sp lines (dashed to distinguish)
                    geom_line(aes_string(y = "Sp_median", color = "scenario"), 
                              linetype = "dashed",
                              size = by_test_grid_threshold_settings$line_size) +
                    geom_point(aes_string(y = "Sp_median", color = "scenario"), 
                               shape = 17,  # Triangle for Sp points
                               size = by_test_grid_threshold_settings$point_size) +
                    facet_wrap(~ test_name, scales = "free_x", dir = "v", 
                               drop = FALSE,
                               labeller = label_wrap_gen(width = by_test_grid_threshold_settings$facet_label_width)) +
                    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                       limits = c(0, 1)) +
                    labs(x = x_label, 
                         y = "Se (solid) / Sp (dashed)",
                         title = "Sensitivity and Specificity vs Threshold by Test",
                         subtitle = "Se shown as solid lines/circles, Sp as dashed lines/triangles") +
                    theme_bw(base_size = by_test_grid_threshold_settings$base_size) +
                    theme(
                      legend.position = by_test_grid_threshold_settings$legend_position,
                      strip.text = element_text(size = by_test_grid_threshold_settings$base_size, face = "bold"),
                      panel.spacing = unit(1, "lines")
                    ) +
                    guides(color = guide_legend(nrow = by_test_grid_threshold_settings$rows_in_legend))
                  
                  if (!show_actual_thresholds) {
                    both_joined_grid_by_test <- both_joined_grid_by_test + 
                      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                  }
                  
                  plots$grid[["both_joined_grid_by_test"]] <- both_joined_grid_by_test
              
            }
          
        }
        
        
        # 5. Grid-separate layout
        if ("grid_separate" %in% layout) {
              
              # Create separate grids for each test showing all scenarios
              # Loop through tests in the correct order
              test_order <- if (!is.null(relevant_thresholds)) {
                seq_along(names(relevant_thresholds))
              } else {
                unique(combined_data$test)
              }
              
              for (t in test_order) {
              # for (t in unique(combined_data$test)) {
                
                test_data <- combined_data %>%
                  filter(test == t)
                
                test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                
                # Create individual Se plots for each scenario for this test
                if ("Se" %in% plot_type || "both" %in% plot_type) {
                  
                      se_scenario_plots <- list()
                      
                      for (s in scenario_names) {
                        
                            scenario_data <- test_data %>%
                              filter(scenario == s)
                            
                            p <- create_threshold_plot_base(
                              data = scenario_data,
                              y_var = "Se_median",
                              x_var = x_var,
                              title = s,
                              x_label = x_label,
                              y_label = "Sensitivity",
                              base_size = by_test_grid_separate_threshold_settings$base_size,
                              point_size = by_test_grid_separate_threshold_settings$point_size,
                              line_size = by_test_grid_separate_threshold_settings$line_size,
                              hide_x_axis = !show_actual_thresholds
                            ) + 
                              theme(
                                axis.text.y = element_text(size = by_test_grid_separate_threshold_settings$base_size * 0.8),
                                legend.position = "none",
                                plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size * 0.8)
                              )
                            
                            # # ADD THIS RIGHT HERE:
                            # p <- p + theme(
                            #   axis.text.x = element_blank(),
                            #   axis.text.y = element_blank(),
                            #   axis.title.x = element_blank(),
                            #   axis.title.y = element_blank()
                            # )
                            
                            # Add bands if requested
                            p <- add_bands_to_threshold(p, 
                                                        scenario_data, 
                                                        "Se", 
                                                        by_test_grid_separate_threshold_settings)
                            
                            se_scenario_plots[[s]] <- p
                        
                      }
                      
                      # Arrange in 3x3 grid
                      se_grid_plot <- wrap_plots(se_scenario_plots, ncol = by_test_grid_separate_threshold_settings$ncol) +
                        plot_annotation(
                          title = paste0("Sensitivity vs Threshold for ", test_name, " across all scenarios"),
                          theme = theme(plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size + 2, 
                                                                  face = "bold"))
                        )
                      
                      se_grid_plot <-  se_grid_plot +  plot_annotation(
                        subtitle = se_grid_plot$plot_annotation$subtitle,
                        caption = "X-axis: Normalized threshold (0-1) | Y-axis: Se (solid) / Sp (dashed)",
                        theme = theme(
                          plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size + 2,  # Title size
                                                    face = "bold"),      
                          plot.subtitle = element_text(size = by_test_grid_separate_threshold_settings$base_size),  # Subtitle size
                          plot.caption = element_text(size = by_test_grid_separate_threshold_settings$base_size,  # Caption size & centered
                                                      hjust = 0.5)   
                        )
                      )
                      
                      plots$grid_separate[[paste0("Se_grid_separate_test_", test_name)]] <- se_grid_plot
                  
                }
                
                # Create individual Sp plots for each scenario for this test
                if ("Sp" %in% plot_type || "both" %in% plot_type) {
                  
                      sp_scenario_plots <- list()
                      
                      for (s in scenario_names) {
                        scenario_data <- test_data %>%
                          filter(scenario == s)
                        
                        p <- create_threshold_plot_base(
                          data = scenario_data,
                          y_var = "Sp_median",
                          x_var = x_var,
                          title = s,
                          x_label = x_label,
                          y_label = "Specificity",
                          base_size = by_test_grid_separate_threshold_settings$base_size,
                          point_size = by_test_grid_separate_threshold_settings$point_size,
                          line_size = by_test_grid_separate_threshold_settings$line_size,
                          hide_x_axis = !show_actual_thresholds
                        ) + 
                          theme(
                            axis.text.y = element_text(size = by_test_grid_separate_threshold_settings$base_size * 0.8),
                            legend.position = "none",
                            plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size * 0.8)
                          )
                        
                        # # ADD THIS RIGHT HERE:
                        # p <- p + theme(
                        #   axis.text.x = element_blank(),
                        #   axis.text.y = element_blank(),
                        #   axis.title.x = element_blank(),
                        #   axis.title.y = element_blank()
                        # )
                        
                        # Add bands if requested
                        p <- add_bands_to_threshold(p, 
                                                    scenario_data,
                                                    "Sp", 
                                                    by_test_grid_separate_threshold_settings)
                        
                        sp_scenario_plots[[s]] <- p
                      }
                      
                      # Arrange in 3x3 grid
                      sp_grid_plot <- wrap_plots(sp_scenario_plots, ncol = by_test_grid_separate_threshold_settings$ncol) +
                        plot_annotation(
                          title = paste0("Specificity vs Threshold for ", test_name, " across all scenarios"),
                          theme = theme(plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size + 2, 
                                                                  face = "bold"))
                        )
                      
                      sp_grid_plot <-  sp_grid_plot +  plot_annotation(
                        subtitle = sp_grid_plot$plot_annotation$subtitle,
                        caption = "X-axis: Normalized threshold (0-1) | Y-axis: Se (solid) / Sp (dashed)",
                        theme = theme(
                          plot.title = element_text(size = by_test_grid_separate_threshold_settings$base_size + 2,   # Title size
                                                    face = "bold"),     
                          plot.subtitle = element_text(size = by_test_grid_separate_threshold_settings$base_size),  # Subtitle size                
                          plot.caption = element_text(size = by_test_grid_separate_threshold_settings$base_size,
                                                      hjust = 0.5)      # Caption size & centered
                        )
                      )
                      
                      plots$grid_separate[[paste0("Sp_grid_separate_test_", test_name)]] <- sp_grid_plot
                  
                }
                
              }
              
              # Also create separate grids for each scenario showing all tests
              for (s in scenario_names) {
                
                scenario_data <- combined_data %>%
                  filter(scenario == s)
                
                # Create individual Se plots for each test for this scenario
                if ("Se" %in% plot_type || "both" %in% plot_type) {
                  
                          se_test_plots <- list()
                          
                          # Loop through tests in the correct order
                          test_order <- if (!is.null(relevant_thresholds)) {
                            seq_along(names(relevant_thresholds))
                          } else {
                            unique(combined_data$test)
                          }
                          
                          for (t in test_order) {
                          # for (t in unique(combined_data$test)) {
                                  
                                  test_data <- scenario_data %>% filter(test == t)
                                  
                                  test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                                  
                                  p <- create_threshold_plot_base(
                                    data = test_data,
                                    y_var = "Se_median",
                                    x_var = x_var,
                                    title = test_name,
                                    x_label = x_label,
                                    y_label = "Sensitivity",
                                    base_size = by_scenario_grid_separate_threshold_settings$base_size,
                                    point_size = by_scenario_grid_separate_threshold_settings$point_size,
                                    line_size = by_scenario_grid_separate_threshold_settings$line_size,
                                    hide_x_axis = !show_actual_thresholds
                                  ) + 
                                    theme(
                                      axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                      legend.position = "none",
                                      plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                    )
                                  
                                  # REPLACE the existing axis removal with this complete version:
                                  p <- p + theme(
                                    axis.text.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                    plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                  )
                                  
                                  # Add bands if requested
                                  p <- add_bands_to_threshold(p, test_data, "Se", by_scenario_grid_separate_threshold_settings)
                                  
                                  se_test_plots[[test_name]] <- p
                            
                          }
                          
                          # Arrange in 2x2 grid
                          se_grid_plot <- wrap_plots(se_test_plots, ncol = by_scenario_grid_separate_threshold_settings$ncol) +
                            plot_annotation(
                              title = paste0("Sensitivity vs Threshold for scenario: ", s),
                              theme = theme(plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2, 
                                                                      face = "bold"))
                            )
                          
                          se_grid_plot <-  se_grid_plot +  plot_annotation(
                            subtitle = se_grid_plot$plot_annotation$subtitle,
                            caption = "X-axis: Normalized threshold (0-1) | Y-axis: Se (solid) / Sp (dashed)",
                            theme = theme(
                              plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2,  # Title size
                                                        face = "bold"),     
                              plot.subtitle = element_text(size = by_scenario_grid_separate_threshold_settings$base_size), # Subtitle size
                              plot.caption = element_text(size = by_scenario_grid_separate_threshold_settings$base_size, # Caption size & centered
                                                          hjust = 0.5)      
                            )
                          )
                          
                          plots$grid_separate[[paste0("Se_grid_separate_scenario_", s)]] <- se_grid_plot
                      
                }
                
                # Create individual Sp plots for each test for this scenario
                if ("Sp" %in% plot_type || "both" %in% plot_type) {
                      
                          sp_test_plots <- list()
                          
                          # Loop through tests in the correct order
                          test_order <- if (!is.null(relevant_thresholds)) {
                            seq_along(names(relevant_thresholds))
                          } else {
                            unique(combined_data$test)
                          }
                          
                          for (t in test_order) {
                          # for (t in unique(combined_data$test)) {
                            
                                test_data <- scenario_data %>%
                                  filter(test == t)
                                
                                test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                                
                                p <- create_threshold_plot_base(
                                  data = test_data,
                                  y_var = "Sp_median",
                                  x_var = x_var,
                                  title = test_name,
                                  x_label = x_label,
                                  y_label = "Specificity",
                                  base_size = by_scenario_grid_separate_threshold_settings$base_size,
                                  point_size = by_scenario_grid_separate_threshold_settings$point_size,
                                  line_size = by_scenario_grid_separate_threshold_settings$line_size,
                                  hide_x_axis = !show_actual_thresholds
                                ) + 
                                  theme(
                                    axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                    legend.position = "none",
                                    plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                  )
                                
                                # REPLACE the existing axis removal with this complete version:
                                p <- p + theme(
                                  axis.text.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                  plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                )
                                
                                # Add bands if requested
                                p <- add_bands_to_threshold(p, test_data, "Sp", by_scenario_grid_separate_threshold_settings)
                                
                                sp_test_plots[[test_name]] <- p
                            
                          }
                          
                          # Arrange in grid
                          sp_grid_plot <- wrap_plots(sp_test_plots, ncol = by_scenario_grid_separate_threshold_settings$ncol) +
                            plot_annotation(
                              title = paste0("Specificity vs Threshold for scenario: ", s),
                              theme = theme(plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2, 
                                                                      face = "bold"))
                            )
                          
                          sp_grid_plot <-  sp_grid_plot +  plot_annotation(
                            subtitle = sp_grid_plot$plot_annotation$subtitle,
                            caption = "X-axis: Normalized threshold (0-1) | Y-axis: Se (solid) / Sp (dashed)",
                            theme = theme(
                              plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2, # Title size
                                                        face = "bold"),     
                              plot.subtitle = element_text(size = by_scenario_grid_separate_threshold_settings$base_size), # Subtitle size
                              plot.caption = element_text(size = by_scenario_grid_separate_threshold_settings$base_size,   # Caption size & centered
                                                          hjust = 0.5)     
                            )
                          )
                          
                          plots$grid_separate[[paste0("Sp_grid_separate_scenario_", s)]] <- sp_grid_plot
                  
                }
                
              }
                
                # In the grid_separate layout section, add this for "both" plot type:
                
                # Combined Se/Sp plots for grid_separate
                if ("both" %in% plot_type) {
                  
                  # For each test, create grid of scenarios with Se/Sp combined
                  # Loop through tests in the correct order
                  test_order <- if (!is.null(relevant_thresholds)) {
                    seq_along(names(relevant_thresholds))
                  } else {
                    unique(combined_data$test)
                  }
                  
                  for (t in test_order) {
                  # for (t in unique(combined_data$test)) {
                    
                        test_data <- combined_data %>%
                          filter(test == t)
                        
                        test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                        
                        # Create individual combined plots for each scenario
                        both_scenario_plots <- list()
     
                        for (s in scenario_names) {
                          
                              scenario_data <- test_data %>%
                                filter(scenario == s)
                              
                              p <- ggplot(scenario_data, aes_string(x = x_var)) +
                                # Se line and points (solid)
                                geom_line(aes_string(y = "Se_median"), 
                                          color = "red",
                                          size = by_scenario_grid_separate_threshold_settings$line_size) +
                                geom_point(aes_string(y = "Se_median"), 
                                           color = "red",
                                           size = by_scenario_grid_separate_threshold_settings$point_size) +
                                # Sp line and points (dashed)
                                geom_line(aes_string(y = "Sp_median"), 
                                          color = "darkgreen",
                                          linetype = "dashed",
                                          size = by_scenario_grid_separate_threshold_settings$line_size) +
                                geom_point(aes_string(y = "Sp_median"), 
                                           color = "darkgreen",
                                           shape = 17,
                                           size = by_scenario_grid_separate_threshold_settings$point_size) +
                                scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                                   limits = c(0, 1)) +
                                labs(x = x_label, 
                                     y = "Se/Sp",
                                     title = s) +
                                theme_bw(base_size = by_scenario_grid_separate_threshold_settings$base_size) +
                                theme(
                                  axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                  legend.position = "none",
                                  plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                )
                              
                              # REPLACE the existing axis removal with this complete version:
                              p <- p + theme(
                                axis.text.x = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.x = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                              )
                              
                              if (!show_actual_thresholds) {
                                p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                              }
                              
                              # Use the combined bands function
                              p <- add_bands_to_threshold_combined(p, 
                                                                   data = scenario_data,
                                                                   settings = by_scenario_grid_separate_threshold_settings)
                              
                              both_scenario_plots[[s]] <- p
                          
                        }
                        
                        # Arrange in 3x3 grid with legend
                        legend_plot <- ggplot() + 
                          geom_point(aes(x = 1, y = 1, color = "Sensitivity"), shape = 16, size = 3) +
                          geom_point(aes(x = 1, y = 0.5, color = "Specificity"), shape = 17, size = 3) +
                          geom_line(aes(x = c(0.8, 1.2), y = c(1, 1), color = "Sensitivity"), size = 1) +
                          geom_line(aes(x = c(0.8, 1.2), y = c(0.5, 0.5), color = "Specificity"), 
                                    linetype = "dashed", size = 1) +
                          scale_color_manual(values = c("Sensitivity" = "red", 
                                                        "Specificity" = "darkgreen")) +
                          theme_void() +
                          theme(legend.position = "bottom",
                                legend.title = element_blank()) +
                          guides(color = guide_legend(override.aes = list(
                            shape = c(16, 17),
                            linetype = c("solid", "dashed")
                          )))
                        
                        grid_plot <- wrap_plots(both_scenario_plots, ncol = by_scenario_grid_separate_threshold_settings$ncol) +
                          plot_annotation(
                            title = paste0("Se & Sp vs Threshold for ", test_name, " across all scenarios"),
                            subtitle = "Red = Sensitivity (solid), Green = Specificity (dashed)",
                            theme = theme(plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2, 
                                                                    face = "bold"))
                          )
                        
                        plots$grid_separate[[paste0("both_joined_grid_separate_test_", test_name)]] <- grid_plot
 
                  
                  # For each scenario, create grid of tests with Se/Sp combined
                  for (s in scenario_names) {
                        
                              scenario_data <- combined_data %>%
                                filter(scenario == s)
                              
                              # Create individual combined plots for each test
                              both_test_plots <- list()
                              
                              # Replace the commented section with this:
                              # Loop through tests in the correct order
                              test_order <- if (!is.null(relevant_thresholds)) {
                                seq_along(names(relevant_thresholds))
                              } else {
                                unique(combined_data$test)
                              }
                              
                              for (t in test_order) {
                              # for (t in unique(combined_data$test)) {
                                    
                                    test_data <- scenario_data %>%
                                      filter(test == t)
                                    
                                    test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
                                    
                                    p <- ggplot(test_data, aes_string(x = x_var)) +
                                      # Se line and points
                                      geom_line(aes_string(y = "Se_median"), 
                                                color = by_scenario_grid_separate_threshold_settings$conf_region_colour_Se,
                                                size = by_scenario_grid_separate_threshold_settings$line_size) +
                                      geom_point(aes_string(y = "Se_median"), 
                                                 color = by_scenario_grid_separate_threshold_settings$conf_region_colour_Se,
                                                 size = by_scenario_grid_separate_threshold_settings$point_size) +
                                      # Sp line and points
                                      geom_line(aes_string(y = "Sp_median"), 
                                                color = by_scenario_grid_separate_threshold_settings$conf_region_colour_Sp,
                                                linetype = "dashed",
                                                size = by_scenario_grid_separate_threshold_settings$line_size) +
                                      geom_point(aes_string(y = "Sp_median"), 
                                                 color = by_scenario_grid_separate_threshold_settings$conf_region_colour_Sp,
                                                 shape = 17,
                                                 size = by_scenario_grid_separate_threshold_settings$point_size) +
                                      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                                         limits = c(0, 1)) +
                                      labs(x = x_label, 
                                           y = "Se/Sp",
                                           title = test_name) +
                                      theme_bw(base_size = by_scenario_grid_separate_threshold_settings$base_size) +
                                      theme(
                                        axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),
                                        legend.position = "none",
                                        plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                      )
                                    
                                    # REPLACE the existing axis removal with this complete version:
                                    p <- p + theme(
                                      axis.text.x = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.text.y = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8),  # Keep Y-axis visible!
                                      plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size * 0.8)
                                    )
                                    
                                    if (!show_actual_thresholds) {
                                      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                                    }
                                    
                                    # Use the combined bands function
                                    p <- add_bands_to_threshold_combined(p, test_data, by_scenario_grid_separate_threshold_settings)
                                    
                                    both_test_plots[[test_name]] <- p
                                
                              }
                              
                              # Arrange in 2x2 grid
                              grid_plot <- wrap_plots(both_test_plots, ncol = by_scenario_grid_separate_threshold_settings$ncol) +
                                plot_annotation(
                                  title = paste0("Se & Sp vs Threshold for scenario: ", s),
                                  subtitle = "Red = Sensitivity (solid), Green = Specificity (dashed)",
                                  theme = theme(plot.title = element_text(size = by_scenario_grid_separate_threshold_settings$base_size + 2, 
                                                                          face = "bold"))
                                )
                              
                              plots$grid_separate[[paste0("both_joined_grid_separate_scenario_", s)]] <- grid_plot
                        
                  }
                  
                }
                
              }
          
        }
        
        # AT THE END, add clean names for single scenario:
        if (n_scenarios == 1) {
          
                clean_plots <- list()
                
                # Main combined Se/Sp plot
                if (!is.null(plots$by_scenario)) {
                  clean_plots$all_tests_combined_Se <- plots$by_scenario$`Se_by_scenario_Baseline (no covariates)`
                  clean_plots$all_tests_combined_Sp <- plots$by_scenario$`Sp_by_scenario_Baseline (no covariates)`
                  clean_plots$all_tests_combined_Se_and_Sp <- plots$by_scenario[[paste0("both_joined_by_scenario_", scenario_names[1])]]
                }
                
                # Individual test plots
                if (!is.null(plots$separate)) {
                  clean_plots$by_test <- plots$separate
                }
                
                # Grid plots
                if (!is.null(plots$grid)) {
                  clean_plots$Se_grid <- plots$grid[["Se_grid_by_test"]]
                  clean_plots$Sp_grid <- plots$grid[["Sp_grid_by_test"]]
                  clean_plots$both_grid <- plots$grid[["both_grid_by_test"]]
                }
                
                # Panels with regions
                if (!is.null(plots$grid_separate) && length(plots$grid_separate) > 0) {
                  # clean_plots$panels_with_regions <- plots$grid_separate[[length(plots$grid_separate)]]
                  clean_plots$panels_with_regions_Se <- plots$grid_separate$`Se_grid_separate_scenario_Baseline (no covariates`
                  clean_plots$panels_with_regions_Sp <- plots$grid_separate$`Sp_grid_separate_scenario_Baseline (no covariates`
                  clean_plots$panels_with_regions_Se_and_Sp <- plots$grid_separate$`both_joined_grid_separate_scenario_Baseline (no covariates`
                  
                }
                
                 # clean_plots$by_scenario <- plots$by_scenario
                
                return(clean_plots)
                
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




  

#' create_confidence_polygon
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

 

 



#' @export
add_regions_to_sroc <- function(
    plot,
    tibble_NMA,
    group_by = c("test", "scenario"),  # Choose grouping variable
    test_names = NULL,
    scenario_names = NULL,
    show_conf = TRUE,
    show_pred = FALSE,
    conf_region_alpha = 0.15,
    pred_region_alpha = 0.30,
    conf_region_colour = NULL,
    pred_region_colour = NULL
) {
  
        require(ggplot2)
        
        group_by <- match.arg(group_by)
        
        polygon_Conf_list <- polygon_Pred_list <- list()
        
        if (group_by == "test") {
          # Original behavior - group by test
          unique_tests <- sort(unique(tibble_NMA$test))
          n_groups <- length(unique_tests)
          
          # Set default test names if not provided
          if (is.null(test_names)) {
            test_names <- paste0("Test ", unique_tests)
          }
          
          for (t in unique_tests) {
            
            tibble_group <- dplyr::filter(tibble_NMA, test == t)
            
            if(nrow(tibble_group) > 0) {
              
              # Confidence region
              tibble_Conf_polygon <- create_confidence_polygon(tibble_group)
              n_rows = nrow(tibble_Conf_polygon)
              tibble_Conf_polygon <- tibble_Conf_polygon %>% 
                dplyr::mutate(test = rep(t, n_rows),
                              test_name = factor(rep(test_names[t], n_rows), levels = test_names),
                              group_var = test_name)  # Generic group variable
              polygon_Conf_list[[t]] <- tibble_Conf_polygon
              
              # Prediction region
              tibble_Pred_polygon <- create_prediction_region_smooth(tibble_group)
              n_rows = nrow(tibble_Pred_polygon)
              tibble_Pred_polygon <- tibble_Pred_polygon %>% 
                dplyr::mutate(test = rep(t, n_rows), 
                              test_name = factor(rep(test_names[t], n_rows), levels = test_names),
                              group_var = test_name)  # Generic group variable
              polygon_Pred_list[[t]] <- tibble_Pred_polygon
            }
          }
          
        } else {  # group_by == "scenario"
          # New behavior - group by scenario
          unique_scenarios <- unique(tibble_NMA$scenario)
          n_groups <- length(unique_scenarios)
          
          # Set default scenario names if not provided
          if (is.null(scenario_names)) {
            scenario_names <- unique_scenarios
          }
          
          for (i in 1:length(unique_scenarios)) {
            s <- unique_scenarios[i]
            
            tibble_group <- dplyr::filter(tibble_NMA, scenario == s)
            
            if(nrow(tibble_group) > 0) {
              
              # Confidence region
              tibble_Conf_polygon <- create_confidence_polygon(tibble_group)
              n_rows = nrow(tibble_Conf_polygon)
              tibble_Conf_polygon <- tibble_Conf_polygon %>% 
                dplyr::mutate(scenario = rep(s, n_rows),
                              scenario_index = rep(i, n_rows),
                              group_var = factor(rep(s, n_rows), levels = scenario_names))
              polygon_Conf_list[[i]] <- tibble_Conf_polygon
              
              # Prediction region
              tibble_Pred_polygon <- create_prediction_region_smooth(tibble_group)
              n_rows = nrow(tibble_Pred_polygon)
              tibble_Pred_polygon <- tibble_Pred_polygon %>% 
                dplyr::mutate(scenario = rep(s, n_rows),
                              scenario_index = rep(i, n_rows),
                              group_var = factor(rep(s, n_rows), levels = scenario_names))
              polygon_Pred_list[[i]] <- tibble_Pred_polygon
            }
          }
        }
        
        polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% 
          dplyr::filter(!(is.na(x)))
        polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% 
          dplyr::filter(!(is.na(x)))
        
        # Check if there's only one group (fixes the pink issue)
        single_group <- (n_groups == 1)
        
        if (show_conf && nrow(polygon_Conf_tibble) > 0) {
          if (!is.null(conf_region_colour) || single_group) {
            # Use fixed color if specified OR if only one group
            color_to_use <- ifelse(!is.null(conf_region_colour), 
                                   conf_region_colour, 
                                   "blue")  # Default to blue
            plot <- plot +
              geom_polygon(data = polygon_Conf_tibble,
                           aes(x = x, y = y, group = group_var),
                           fill = color_to_use,
                           alpha = conf_region_alpha,
                           inherit.aes = FALSE)
          } else {
            # Inherit color from the plot (only when multiple groups)
            plot <- plot +
              geom_polygon(data = polygon_Conf_tibble,
                           aes(x = x, y = y, fill = group_var, group = group_var),
                           alpha = conf_region_alpha,
                           inherit.aes = FALSE)
          }
        }
        
        if (show_pred && nrow(polygon_Pred_tibble) > 0) {
          if (!is.null(pred_region_colour) || single_group) {
            # Use fixed color if specified OR if only one group
            color_to_use <- ifelse(!is.null(pred_region_colour), 
                                   pred_region_colour, 
                                   "green")  # Default to green
            plot <- plot +
              geom_polygon(data = polygon_Pred_tibble,
                           aes(x = x, y = y, group = group_var),
                           fill = color_to_use,
                           alpha = pred_region_alpha,
                           inherit.aes = FALSE)
          } else {
            # Inherit color from the plot (only when multiple groups)
            plot <- plot +
              geom_polygon(data = polygon_Pred_tibble,
                           aes(x = x, y = y, fill = group_var, group = group_var),
                           alpha = pred_region_alpha,
                           inherit.aes = FALSE)
          }
        }
        
        return(plot)
  
}






# #' Create forest plots for specific comparisons across scenarios (GENERAL VERSION)
# #' @export
# R_fn_plot_forest_comparisons_general <- function(
#                                                 all_comparisons_data,
#                                                 selected_comparisons = NULL,
#                                                 comparison_type = c("Se", "Sp"),
#                                                 scenario_var = "scenario",  # Column name containing scenario labels
#                                                 group_var = NULL,           # Optional: column to use for coloring
#                                                 shape_var = NULL,           # Optional: column to use for shapes
#                                                 n_auto_select = 6,          # Number of comparisons to auto-select
#                                                 base_size = 14,
#                                                 point_size = 3,
#                                                 error_bar_width = 0.2,
#                                                 error_bar_alpha = 0.7,
#                                                 facet_ncol = 2,
#                                                 facet_scales = "free_x",
#                                                 plot_title = NULL,
#                                                 x_label = NULL,
#                                                 y_label = "Scenario",
#                                                 save_plot = FALSE,
#                                                 output_path = NULL
# ) {
#   
#         require(ggplot2)
#         require(dplyr)
#         
#         comparison_type <- match.arg(comparison_type)
#         
#         # Set default labels based on comparison type
#         if (is.null(x_label)) {
#           x_label <- paste0("Difference in ", 
#                             ifelse(comparison_type == "Se", "Sensitivity", "Specificity"))
#         }
#         
#         if (is.null(plot_title)) {
#           plot_title <- paste0(ifelse(comparison_type == "Se", "Sensitivity", "Specificity"),
#                                " Differences Across Scenarios")
#         }
#         
#         # Create comparison label if it doesn't exist
#         if (!"comparison_label" %in% names(all_comparisons_data)) {
#           all_comparisons_data <- all_comparisons_data %>%
#             mutate(comparison_label = paste0(test1_name, " (", threshold1, ") vs ",
#                                              test2_name, " (", threshold2, ")"))
#         }
#         
#         # Auto-select comparisons if not provided
#         if (is.null(selected_comparisons)) {
#           # Find comparisons with largest variance across scenarios
#           variance_summary <- all_comparisons_data %>%
#             group_by(comparison_label) %>%
#             summarise(
#               var_diff = var(get(paste0("diff_", comparison_type, "_median")), na.rm = TRUE),
#               mean_abs_diff = mean(abs(get(paste0("diff_", comparison_type, "_median"))), na.rm = TRUE),
#               .groups = "drop"
#             ) %>%
#             filter(mean_abs_diff > 0.05) %>%  # At least 5% difference
#             arrange(desc(var_diff)) %>%
#             slice_head(n = n_auto_select)
#           
#           selected_comparisons <- variance_summary$comparison_label
#         }
#         
#         # Filter to selected comparisons
#         plot_data <- all_comparisons_data %>%
#           filter(comparison_label %in% selected_comparisons)
#         
#         # Get the correct column names
#         x_col <- paste0("diff_", comparison_type, "_median")
#         xmin_col <- paste0("diff_", comparison_type, "_lower")
#         xmax_col <- paste0("diff_", comparison_type, "_upper")
#         
#         # Create base plot
#         p <- ggplot(plot_data, 
#                     aes_string(y = scenario_var, x = x_col)) +
#           geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
#           geom_errorbarh(aes_string(xmin = xmin_col, xmax = xmax_col), 
#                          height = error_bar_width, 
#                          alpha = error_bar_alpha)
#         
#         # Add points with optional color and shape mappings
#         if (!is.null(group_var) && !is.null(shape_var)) {
#           p <- p + geom_point(aes_string(color = group_var, shape = shape_var), 
#                               size = point_size)
#         } else if (!is.null(group_var)) {
#           p <- p + geom_point(aes_string(color = group_var), 
#                               size = point_size)
#         } else if (!is.null(shape_var)) {
#           p <- p + geom_point(aes_string(shape = shape_var), 
#                               size = point_size)
#         } else {
#           p <- p + geom_point(size = point_size)
#         }
#         
#         # Add faceting
#         p <- p + 
#           facet_wrap(~ comparison_label, ncol = facet_ncol, scales = facet_scales, dir = "v") +
#           labs(x = x_label,
#                y = y_label,
#                title = plot_title) +
#           theme_bw(base_size = base_size) +
#           theme(
#             axis.text.y = element_text(size = base_size * 0.8),
#             legend.position = "bottom",
#             strip.text = element_text(size = base_size * 0.8, face = "bold"),
#             panel.spacing = unit(1, "lines")
#           )
#         
#         # Save if requested
#         if (save_plot && !is.null(output_path)) {
#           ggsave(output_path, p, width = 10, height = 8)
#         }
#         
#         return(p)
#   
# }


`%||%` <- function(x, y) if (is.null(x)) y else x




#' Plot pairwise differences (Se or Sp) - works for both normal and meta-regression
#' @export
plot_pairwise_differences <- function(  tibble_comp,
                                        comparison_type = c("Se", "Sp", "both"),
                                        ##
                                        test_names = NULL,
                                        relevant_thresholds = NULL,
                                        ##
                                        baseline_plot_settings = NULL,
                                        metareg_plot_settings = NULL
) {
      
      require(ggplot2)
      require(dplyr)
      require(forcats)
      require(patchwork)
      
      comparison_type <- match.arg(comparison_type, several.ok = TRUE)
      ##
      will_combine <- ("both" %in% comparison_type)
      
      # Filter if thresholds provided
      if (!is.null(relevant_thresholds)) {
        tibble_comp <- filter_clinically_relevant_comparisons(
          tibble_comp = tibble_comp,
          relevant_thresholds = relevant_thresholds
        )
      }
      
      # Add test names if provided
      if (!is.null(test_names)) {
        tibble_comp <- tibble_comp %>%
          mutate(
            test1_name = test_names[test1],
            test2_name = test_names[test2]
          )
      }
      
      # Ensure comparison_label exists
      if (!"comparison_label" %in% names(tibble_comp)) {
        tibble_comp <- tibble_comp %>%
          mutate(
            comparison_label = paste0(
              test1_name, " (≥", threshold1, ") vs ",
              test2_name, " (≥", threshold2, ")"
            )
          )
      }
      
      # Check if this is meta-regression (multiple scenarios)
      has_scenarios <- "scenario" %in% names(tibble_comp) && length(unique(tibble_comp$scenario)) > 1
      
      # Set defaults for baseline_plot_settings
      if (is.null(baseline_plot_settings)) baseline_plot_settings <- list()
      if (is.null(baseline_plot_settings$base_size)) baseline_plot_settings$base_size <- 14
      if (is.null(baseline_plot_settings$point_size)) baseline_plot_settings$point_size <- 3
      if (is.null(baseline_plot_settings$error_bar_width)) baseline_plot_settings$error_bar_width <- 0.2
      if (is.null(baseline_plot_settings$error_bar_alpha)) baseline_plot_settings$error_bar_alpha <- 0.5
      if (is.null(baseline_plot_settings$y_axis_text_size_multiplier)) baseline_plot_settings$y_axis_text_size_multiplier <- 1.0
      if (is.null(baseline_plot_settings$show_values)) baseline_plot_settings$show_values <- FALSE
      if (is.null(baseline_plot_settings$value_size_multiplier)) baseline_plot_settings$value_size_multiplier <- 0.8
      if (is.null(baseline_plot_settings$ratio_plot))  baseline_plot_settings$ratio_plot <- 3.0
      if (is.null(baseline_plot_settings$ratio_text))  baseline_plot_settings$ratio_text <- 1.5
      ##
      baseline_plot_settings$y_axis_text_size <- baseline_plot_settings$y_axis_text_size_multiplier * baseline_plot_settings$base_size * 0.8
      
      
      # Set defaults for metareg_plot_settings
      if (is.null(metareg_plot_settings)) metareg_plot_settings <- list()
      if (is.null(metareg_plot_settings$base_size)) metareg_plot_settings$base_size <- 14
      if (is.null(metareg_plot_settings$point_size)) metareg_plot_settings$point_size <- 3
      # if (is.null(metareg_plot_settings$line_size)) metareg_plot_settings$line_size <- 0.8
      if (is.null(metareg_plot_settings$error_bar_width)) metareg_plot_settings$error_bar_width <- 0.2
      if (is.null(metareg_plot_settings$error_bar_alpha)) metareg_plot_settings$error_bar_alpha <- 0.7
      if (is.null(metareg_plot_settings$legend_position)) metareg_plot_settings$legend_position <- "bottom" ## 
      if (is.null(metareg_plot_settings$facet_ncol)) metareg_plot_settings$facet_ncol <- 2 ## 
      if (is.null(metareg_plot_settings$facet_scales)) metareg_plot_settings$facet_scales <- "fixed" ##
      if (is.null(metareg_plot_settings$n_auto_select)) metareg_plot_settings$n_auto_select <- 6 ##
      if (is.null(metareg_plot_settings$min_mean_abs_diff)) metareg_plot_settings$min_mean_abs_diff <- 0.01
      if (is.null(metareg_plot_settings$strip_text_size_multiplier)) metareg_plot_settings$strip_text_size_multiplier <- 0.8 ##
      if (is.null(metareg_plot_settings$show_values)) metareg_plot_settings$show_values <- FALSE
      if (is.null(metareg_plot_settings$value_size_multiplier)) metareg_plot_settings$value_size_multiplier <- 0.8
      if (is.null(metareg_plot_settings$ratio_plot))  metareg_plot_settings$ratio_plot <- 3.0
      if (is.null(metareg_plot_settings$ratio_text))  metareg_plot_settings$ratio_text <- 1.5
      if (is.null(metareg_plot_settings$facet_label_width))  metareg_plot_settings$facet_label_width <- 25
      ##
      metareg_plot_settings$y_axis_text_size <- metareg_plot_settings$y_axis_text_size_multiplier * metareg_plot_settings$base_size * 0.8
      
      
      plots <- list()
    
      if (has_scenarios) {   # CASE 1: Meta-regression (multiple scenarios)
        
                # For meta-regression, create forest plots with comparisons on y-axis, faceted by scenario
                if ("Se" %in% comparison_type || "both" %in% comparison_type) {
                          
                          # Auto-select most variable comparisons
                          se_variance <- tibble_comp %>%
                            group_by(comparison_label) %>%
                            summarise(
                              var_diff = var(diff_Se_median, na.rm = TRUE),
                              mean_abs = mean(abs(diff_Se_median), na.rm = TRUE),
                              .groups = "drop"
                            ) %>%
                            filter(!is.na(var_diff), mean_abs > metareg_plot_settings$min_mean_abs_diff) %>%
                            arrange(desc(var_diff)) %>%
                            slice_head(n = metareg_plot_settings$n_auto_select)
                          
                          # Get scenario order from the data (preserves original order)
                          scenario_levels <- unique(tibble_comp$scenario)
                          
                          plot_data_se <- tibble_comp %>%
                            filter(comparison_label %in% se_variance$comparison_label) %>%
                            mutate(
                              scenario = factor(scenario, levels = scenario_levels),
                              comparison_label = factor(comparison_label, levels = se_variance$comparison_label)
                            )
                          
                          p_se <- ggplot(plot_data_se, 
                                         aes(y = comparison_label,
                                             x = diff_Se_median)) +
                            geom_vline(xintercept = 0, 
                                       linetype = "dashed", 
                                       alpha = 0.5) +
                            geom_errorbarh(aes(xmin = diff_Se_lower,
                                               xmax = diff_Se_upper), 
                                           height = metareg_plot_settings$error_bar_width,
                                           alpha = metareg_plot_settings$error_bar_alpha) +
                            geom_point(size = metareg_plot_settings$point_size) +
                            facet_wrap(~ scenario, 
                                       drop = FALSE,
                                       ncol = metareg_plot_settings$facet_ncol,
                                       scales = metareg_plot_settings$facet_scales,
                                       labeller = label_wrap_gen(width = metareg_plot_settings$facet_label_width)) +
                            labs(x = "Difference in Sensitivity",
                                 y = "",
                                 title = "Sensitivity Differences Across Scenarios") +
                            theme_bw(base_size = metareg_plot_settings$base_size) +
                            theme(
                              legend.position = metareg_plot_settings$legend_position,
                              strip.text = element_text(size = metareg_plot_settings$base_size * metareg_plot_settings$strip_text_size_multiplier),
                              axis.text.y = element_text(size = metareg_plot_settings$base_size * 0.7)
                            )
                          
                          if (metareg_plot_settings$show_values) {
                                
                                # Create text data
                                text_data <- plot_data_se %>%
                                  mutate(
                                    value_label = sprintf("%.1f%% (%.1f, %.1f)", 
                                                          diff_Se_median * 100,
                                                          diff_Se_lower * 100, 
                                                          diff_Se_upper * 100)
                                  )
                                
                                # Text-only plot
                                text_plot <- ggplot(text_data, aes(y = comparison_label, x = 0)) +
                                  geom_text(aes(label = value_label),
                                            hjust = 0, 
                                            size = metareg_plot_settings$base_size * metareg_plot_settings$value_size_multiplier / 3.00) +
                                  facet_wrap(~ scenario, 
                                             drop = FALSE,
                                             ncol = metareg_plot_settings$facet_ncol,
                                             scales = metareg_plot_settings$facet_scales,
                                             labeller = label_wrap_gen(width = metareg_plot_settings$facet_label_width)) + 
                                  theme_void() +
                                  theme(
                                    strip.text = element_blank(),
                                    panel.spacing = theme_get()$panel.spacing
                                  ) +
                                  coord_cartesian(xlim = c(-1, 100), clip = "off")
                                
                                # Combine
                                p_se_combined <- p_se + text_plot + 
                                  plot_layout(widths = c(metareg_plot_settings$ratio_plot, 
                                                         metareg_plot_settings$ratio_text))
                                
                                plots$Se_forest <- p_se_combined
                            
                          } else {
                            
                                plots$Se_forest <- p_se
                            
                          }
                  
                }
                
                if ("Sp" %in% comparison_type || "both" %in% comparison_type) {
                  
                          # Auto-select most variable comparisons for Sp
                          sp_variance <- tibble_comp %>%
                            group_by(comparison_label) %>%
                            summarise(
                              var_diff = var(diff_Sp_median, na.rm = TRUE),
                              mean_abs = mean(abs(diff_Sp_median), na.rm = TRUE),
                              .groups = "drop"
                            ) %>%
                            filter(!is.na(var_diff), mean_abs > metareg_plot_settings$min_mean_abs_diff) %>%
                            arrange(desc(var_diff)) %>%
                            slice_head(n = metareg_plot_settings$n_auto_select)
                          
                          # Get scenario order from the data (preserves original order)
                          scenario_levels <- unique(tibble_comp$scenario)
                          
                          plot_data_sp <- tibble_comp %>%
                            filter(comparison_label %in% sp_variance$comparison_label) %>%
                            mutate(
                              scenario = factor(scenario, levels = scenario_levels),
                              comparison_label = factor(comparison_label, levels = sp_variance$comparison_label)
                            )
                          
                          p_sp <- ggplot(plot_data_sp, 
                                         aes(y = comparison_label, 
                                             x = diff_Sp_median)) +
                            geom_vline(xintercept = 0,
                                       linetype = "dashed",
                                       alpha = 0.5) +
                            geom_errorbarh(aes(xmin = diff_Sp_lower,
                                               xmax = diff_Sp_upper), 
                                           height = metareg_plot_settings$error_bar_width,
                                           alpha = metareg_plot_settings$error_bar_alpha) +
                            geom_point(size = metareg_plot_settings$point_size) +
                            facet_wrap(~ scenario, 
                                       drop = FALSE,
                                       ncol = metareg_plot_settings$facet_ncol, 
                                       scales = metareg_plot_settings$facet_scales,
                                       labeller = label_wrap_gen(width = metareg_plot_settings$facet_label_width)) +
                            labs(x = "Difference in Specificity",
                                 y = "",
                                 title = "Specificity Differences Across Scenarios") +
                            theme_bw(base_size = metareg_plot_settings$base_size) +
                            theme(
                              legend.position = metareg_plot_settings$legend_position,
                              strip.text = element_text(size = metareg_plot_settings$base_size * metareg_plot_settings$strip_text_size_multiplier),
                              axis.text.y = element_text(size = metareg_plot_settings$base_size * 0.7)
                            )
                          
                          if (metareg_plot_settings$show_values) {
                            
                                # Create text data
                                text_data <- plot_data_sp %>%
                                  mutate(
                                    value_label = sprintf("%.1f%% (%.1f, %.1f)", 
                                                          diff_Sp_median * 100,
                                                          diff_Sp_lower * 100, 
                                                          diff_Sp_upper * 100)
                                  )
                                
                                # Text-only plot
                                text_plot <- ggplot(text_data, aes(y = comparison_label, x = 0)) +
                                  geom_text(aes(label = value_label),
                                            hjust = 0, 
                                            size = metareg_plot_settings$base_size * metareg_plot_settings$value_size_multiplier / 3.00) +
                                  facet_wrap(~ scenario, 
                                             drop = FALSE,
                                             ncol = metareg_plot_settings$facet_ncol,
                                             scales = metareg_plot_settings$facet_scales, 
                                             labeller = label_wrap_gen(width = metareg_plot_settings$facet_label_width)) +
                                  theme_void() +
                                  theme(
                                    strip.text = element_blank(),
                                    panel.spacing = theme_get()$panel.spacing
                                  ) +
                                  coord_cartesian(xlim = c(-1, 100), clip = "off")
                                
                                # Combine
                                p_sp_combined <- p_sp + text_plot + 
                                  plot_layout(widths = c(metareg_plot_settings$ratio_plot, 
                                                         metareg_plot_settings$ratio_text))
                                
                                plots$Sp_forest <- p_sp_combined
                            
                          } else {
                            
                                plots$Sp_forest <- p_sp
                            
                          }
                  
                }
            
      } else {   # CASE 2: Normal (single scenario or no scenarios)

                # Add comparison type for faceting
                plot_data <-  tibble_comp %>%
                              mutate(
                                same_test = (test1 == test2),
                                comparison_type_facet = ifelse(same_test, 
                                                               paste0(test1_name, " thresholds"),
                                                               paste0(test1_name, " vs ", test2_name))
                              ) %>%
                              arrange(test1, threshold1, test2, threshold2) %>%
                              mutate(
                                # Convert to factor with levels in the order they appear after arranging
                                comparison_type_facet = factor(comparison_type_facet, 
                                                               levels = unique(comparison_type_facet))
                              ) %>%
                              group_by(comparison_type_facet) %>%
                              mutate(
                                order_within_facet = row_number(),
                                comparison_label_ordered = fct_reorder(comparison_label, 
                                                                       order_within_facet, .desc = TRUE)
                              ) %>%
                              ungroup()
                
                if ("Se" %in% comparison_type || "both" %in% comparison_type) {
                  
                      # # Calculate x position for text if not specified
                      # if (baseline_plot_settings$show_values && is.null(baseline_plot_settings$value_x_position)) {
                      #   x_range <- range(c(plot_data$diff_Se_lower, plot_data$diff_Se_upper), na.rm = TRUE)
                      #   text_x_position <- x_range[2] + diff(x_range) * 0.15  # 15% beyond max
                      # } else {
                      #   text_x_position <- baseline_plot_settings$value_x_position
                      # }
                  
                      p_se <- ggplot(plot_data, 
                                     aes(y = comparison_label_ordered,  ## value_size_multiplier
                                         x = diff_Se_median)) +
                              geom_vline(xintercept = 0,
                                         linetype = "dashed", 
                                         alpha = 0.5) +
                              geom_errorbarh(aes(xmin = diff_Se_lower, xmax = diff_Se_upper), 
                                             height = baseline_plot_settings$error_bar_width, 
                                             alpha =  baseline_plot_settings$error_bar_alpha) +
                              geom_point(size =  baseline_plot_settings$point_size) 
                      
                      p_se <- p_se + 
                              facet_wrap(~ comparison_type_facet,
                                   drop = FALSE,
                                   scales = "free_y", 
                                   ncol = 1) +
                              labs(x = "Difference in Sensitivity",
                                   y = "",
                                   title = "Sensitivity Differences") +
                              theme_bw(base_size =  baseline_plot_settings$base_size) + 
                              theme(
                                axis.text.y = element_text(size = baseline_plot_settings$y_axis_text_size),
                                legend.position = "none"
                              )
                      
                      if (baseline_plot_settings$show_values) {
                        
                              # Create a text-only plot
                              text_data <- plot_data %>%
                                mutate(
                                  value_label = sprintf("%.1f%% (%.1f, %.1f)", 
                                                        diff_Se_median * 100,
                                                        diff_Se_lower * 100, 
                                                        diff_Se_upper * 100)
                                )
                              
                              text_plot <- ggplot(text_data, aes(y = comparison_label_ordered, x = 0)) +
                                geom_text(aes(label = value_label),
                                          hjust = 0, 
                                          size = baseline_plot_settings$base_size * baseline_plot_settings$value_size_multiplier / 3.00) +
                                facet_wrap(~ comparison_type_facet,
                                           drop = FALSE, 
                                           scales = "free_y",
                                           ncol = 1) +
                                theme_void() +
                                theme(
                                  strip.text = element_blank(),
                                  panel.spacing = theme_get()$panel.spacing
                                ) +
                                coord_cartesian(xlim = c(-1, 100), 
                                                clip = "off")   
                                
                                # Combine main plot with text panel
                                p_se_combined <- p_se + text_plot + plot_layout(widths = c(baseline_plot_settings$ratio_plot, 
                                                                                           baseline_plot_settings$ratio_text))
                                
                                plots$Se <- p_se_combined
                        
                      } else { 
                                plots$Se <- p_se
                      }
                }
                
                if ("Sp" %in% comparison_type || "both" %in% comparison_type) {
                  
                      p_sp <- ggplot(plot_data, 
                                     aes(y = comparison_label_ordered,
                                         x = diff_Sp_median)) +
                        geom_vline(xintercept = 0, 
                                   linetype = "dashed", 
                                   alpha = 0.5) +
                        geom_errorbarh(aes(xmin = diff_Sp_lower,
                                           xmax = diff_Sp_upper), 
                                       height =  baseline_plot_settings$error_bar_width, 
                                       alpha =  baseline_plot_settings$error_bar_alpha) +
                        geom_point(size =  baseline_plot_settings$point_size) +
                        facet_wrap(~ comparison_type_facet, 
                                   drop = FALSE,
                                   scales = "free_y", 
                                   ncol = 1) +
                        labs(x = "Difference in Specificity",
                             y = "",
                             title = "Specificity Differences") +
                        theme_bw(base_size =  baseline_plot_settings$base_size) +
                        theme(
                          axis.text.y = element_text(size =  baseline_plot_settings$y_axis_text_size),
                          legend.position = "none"
                        )
                      
                      # If we're combining, remove y-axis elements from Sp
                      if (will_combine) {
                        p_sp <- p_sp + 
                          theme(
                            axis.text.y = element_blank(),
                            axis.ticks.y = element_blank()
                          )
                      }
                      
                      if (baseline_plot_settings$show_values) {
                        
                              # Create a text-only plot
                              text_data <- plot_data %>%
                                mutate(
                                  value_label = sprintf("%.1f%% (%.1f, %.1f)", 
                                                        diff_Sp_median * 100,
                                                        diff_Sp_lower * 100, 
                                                        diff_Sp_upper * 100)
                                )
                              
                              text_plot <- ggplot(text_data, aes(y = comparison_label_ordered, x = 0)) +
                                geom_text(aes(label = value_label),
                                          hjust = 0, 
                                          size = baseline_plot_settings$base_size * baseline_plot_settings$value_size_multiplier / 3.00) +
                                facet_wrap(~ comparison_type_facet,
                                           drop = FALSE, 
                                           scales = "free_y", ncol = 1) +
                                theme_void() +
                                theme(
                                  strip.text = element_blank(),
                                  panel.spacing = theme_get()$panel.spacing
                                ) +
                                coord_cartesian(xlim = c(-1, 100),
                                                clip = "off")   
                              
                              # Combine main plot with text panel
                              p_sp_combined <- p_sp + text_plot + plot_layout(widths = c(baseline_plot_settings$ratio_plot, 
                                                                                         baseline_plot_settings$ratio_text))
                              
                              plots$Sp <- p_sp_combined
                        
                      } else { 
                              plots$Sp <- p_sp
                      }
                      
                }
        
      }
      
      # Combine if both requested
      if ("both" %in% comparison_type && length(plots) == 2) {
        
            if (baseline_plot_settings$show_values) { 
              n_cols <- 4
            } else { 
              n_cols <- 2
            }
            
            # if (has_scenarios) {
            #   plots$combined <- plots$Se_forest + plots$Sp_forest + 
            #                     plot_layout(ncol = n_cols)
            # } else {
              
              plots$combined <- plots$Se + plots$Sp +
                plot_layout(ncol = n_cols)
              
              # Then modify just the second plot in the combined layout
              plots$combined[[2]] <- plots$combined[[2]] + 
                theme(
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()
                )
              
            # }
        
      }
      
      # Return single plot if only one type requested
      if (length(comparison_type) == 1 && comparison_type != "both") {
        return(plots[[1]])
      }
      
      return(plots)
      
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
                                        x_label = "Threshold",
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
        
        # Create base plot with required aesthetics
        p <- ggplot(data, aes_string(x = x_var, y = y_var))
        
        # Add optional aesthetics one by one
        if (!is.null(color_by)) {
          p <- p + aes_string(color = color_by)
        }
        if (!is.null(linetype_by)) {
          p <- p + aes_string(linetype = linetype_by)
        }
        if (!is.null(shape_by)) {
          p <- p + aes_string(shape = shape_by)
        }
        if (!is.null(group_by)) {
          p <- p + aes_string(group = group_by)
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
          p <- p + facet_wrap(dir = "v", as.formula(paste("~", facet_by)), scales = "free_x")
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



 






