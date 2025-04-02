



 




#' R_fn_sROC_plot_NMA
#' @keywords internal
#' @export
R_fn_sROC_plot_NMA <- function( stan_model_file_name,
                                stan_mod_samples,
                                ##
                                df_true = NULL,
                                ##
                                conf_region_colour = "blue", 
                                pred_region_colour = "blue",
                                ##
                                n_index_tests,
                                n_thr
) {
                
                require(ggplot2)
                require(dplyr)
                
                model_name <- stan_model_file_name
  
                n_thr_max <- max(n_thr)
                ##
                Se_median_array <- Sp_median_array <- Fp_median_array <- array(dim = c(n_index_tests, n_thr_max))
                Se_mean_array   <- Sp_mean_array   <- Fp_mean_array   <- array(dim = c(n_index_tests, n_thr_max))
                ##
                Se_lower_array <- Sp_lower_array <- Fp_lower_array <- array(dim = c(n_index_tests, n_thr_max))
                Se_upper_array <- Sp_upper_array <- Fp_upper_array <- array(dim = c(n_index_tests, n_thr_max))
                ##
                Se_pred_lower_array <- Sp_pred_lower_array <- Fp_pred_lower_array <- array(dim = c(n_index_tests, n_thr_max))
                Se_pred_upper_array <- Sp_pred_upper_array <- Fp_pred_upper_array <- array(dim = c(n_index_tests, n_thr_max))
                ##
                
                # tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 100)
                # ##
                # tibble_gq <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 100)
                ##
                Se <- stan_mod_samples$summary(c("Se"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                Sp <- stan_mod_samples$summary(c("Sp"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                Fp <- stan_mod_samples$summary(c("Fp"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                ##
                # Se <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se"))  & (!(stringr::str_detect(parameter, "Se_pred"))))
                # Sp <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp"))  & (!(stringr::str_detect(parameter, "Sp_pred"))))
                # Fp  <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp")) & (!(stringr::str_detect(parameter, "Fp_pred"))))
                ##
                # Se_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_pred")))
                # Sp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_pred")))
                # Fp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_pred")))
                Se_pred <- stan_mod_samples$summary(c("Se_pred"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                Sp_pred <- stan_mod_samples$summary(c("Sp_pred"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                Fp_pred <- stan_mod_samples$summary(c("Fp_pred"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
                
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
                
                test_vec <- thr_vec <- c()
                counter <- 1 
                for (t in 1:n_index_tests) {
                  for (k in 1:n_thr_max) {
                    test_vec[counter] <- t
                    thr_vec[counter] <- k
                    ##
                    counter <- counter + 1
                  }
                }
                
                n_rows_total <- length(test_vec)
                
                ##
                ## ---- Make tibble for ggplot (to make sROC plots, etc):
                ##
                # Se_mean_array
                # c(t(Se_mean_array))
                # ##
                # stan_model_file_name <- "NMA_Xu"
                # model_name <- stan_model_file_name
                ##
                tibble_NMA <- tibble(  Model = rep(model_name, n_rows_total),
                                       ##
                                       test = factor(test_vec),
                                       test_char = paste0("test ", test),
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
                # tibble_NMA
                polygon_Conf_list <- polygon_Pred_list <- list()
                ##
                for (t in 1:n_index_tests) {
                  tibble_test_t <- dplyr::filter(tibble_NMA, test == t)
                  ##
                  ## ---- Tibbles for 95% credible region:
                  ##
                  tibble_Conf_polygon_test_t <- create_confidence_polygon(tibble_test_t, model_name = model_name)
                  n_rows = nrow(tibble_Conf_polygon_test_t)
                  tibble_Conf_polygon_test_t <- tibble_Conf_polygon_test_t %>% dplyr::mutate(test = rep(t, n_rows),
                                                                                             test_char = paste0("test ", test))
                  polygon_Conf_list[[t]] <- tibble_Conf_polygon_test_t
                  ##
                  ## ---- Tibbles for 95% prediction region:
                  ##
                  tibble_Pred_polygon_test_t <- create_prediction_polygon(tibble_test_t, model_name = model_name)
                  n_rows = nrow(tibble_Pred_polygon_test_t)
                  tibble_Pred_polygon_test_t <- tibble_Pred_polygon_test_t %>% dplyr::mutate(test = rep(t, n_rows), 
                                                                                             test_char = paste0("test ", test))
                  polygon_Pred_list[[t]] <- tibble_Pred_polygon_test_t
                  # polygon_Pred_list[[t]] <- create_prediction_polygon(tibble_test_t, model_name = model_name)
                }
                ##
                ## ---- Final conf and prediction regions tibbles:
                ##
                polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% dplyr::filter(!(is.na(x))) %>% print(n = 100)
                polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% dplyr::filter(!(is.na(x))) %>% print(n = 100)
                
                {
                      plot_list <- list()
                      ##
                      ## ---- Plot 1 (all tests on 1 plot):
                      ##
                      plot_1 <-  ggplot(tibble_NMA, mapping = aes(x = Fp_median, y = Se_median, colour = test)) + 
                        geom_line(linewidth = 0.5) + 
                        geom_point(size = 3) + 
                        theme_bw(base_size = 16) + 
                        xlab("False positive rate (Fp)") + 
                        ylab("Sensitivity (Se)")
                      print(plot_1)
                      plot_list[[1]] <- plot_1
                      ##
                      ## ---- Plot 2 (each test on separate panel):
                      ##
                      plot_2 <-  ggplot(tibble_NMA, mapping = aes(x = Fp_median, y = Se_median, colour = Model)) + 
                        geom_line(linewidth = 0.5) + 
                        geom_point(size = 3) + 
                        theme_bw(base_size = 16) + 
                        facet_wrap(~ test_char) + 
                        xlab("False positive rate (Fp)") + 
                        ylab("Sensitivity (Se)")
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
                                     alpha = 0.40)
                      # geom_polygon(data = polygon_Pred_tibble, aes(x = x, y = y), fill = pred_region_colour, alpha = 0.25) + 
                      ##
                      print(plot_3)
                      plot_list[[3]] <- plot_3
                      ##
                      ## ---- Plot 4 (plot w/ 95% prediction region):
                      ##
                      plot_4 <-  plot_2 +  geom_polygon(
                        data = polygon_Pred_tibble,
                        aes(x = x, y = y, colour = Model),
                        fill = pred_region_colour, 
                        alpha = 0.40)
                      print(plot_4)
                      plot_list[[4]] <- plot_4
                      ##
                      ## ---- Plot 5 (plot w/ BOTH the 95% confidence + prediction regions):
                      ##
                      plot_5 <-  plot_2 + 
                        geom_polygon(
                          data = polygon_Conf_tibble,
                          aes(x = x, y = y, colour = Model),
                          fill = pred_region_colour, 
                          alpha = 0.40) + 
                        ##
                        geom_polygon(
                          data = polygon_Pred_tibble,
                          aes(x = x, y = y, colour = Model),
                          fill = pred_region_colour, 
                          alpha = 0.20)
                      print(plot_5)
                      plot_list[[5]] <- plot_5
                }
                
                return(list(plot_list = plot_list,
                            ##
                            tibble_NMA = tibble_NMA,
                            ##
                            polygon_Conf_list   = polygon_Conf_list,
                            polygon_Conf_tibble = polygon_Conf_tibble,
                            polygon_Pred_list   = polygon_Pred_list,
                            polygon_Pred_tibble = polygon_Pred_tibble))
                
              
   

}




 























