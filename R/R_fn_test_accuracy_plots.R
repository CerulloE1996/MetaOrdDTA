



 
#' create_confidence_polygon - Create credible interval polygon function
#' @keywords internal
#' @export
create_confidence_polygon <- function(df, 
                                      model_name) {
          
          # Create upper curve (high sensitivity, low false positive)
          upper <- df %>% 
            dplyr::filter(Model == model_name) %>%
            dplyr::select(x = Fp_lower, y = Se_upper) %>%
            dplyr::arrange(x)
          
          # Create lower curve (low sensitivity, high false positive)
          lower <- df %>% 
            dplyr::filter(Model == model_name) %>%
            dplyr::select(x = Fp_upper, y = Se_lower) %>%
            dplyr::arrange(desc(x))
          
          # Combine into a closed polygon
          out <- tibble(bind_rows(upper, lower))
          n_rows <- nrow(out)
          out <- out %>%
            dplyr::mutate(Model = rep(model_name, n_rows))

}









#' create_prediction_polygon -  Create prediction interval polygon function
#' @keywords internal
#' @export
create_prediction_polygon <- function(df, 
                                      model_name) {
          
          # Create upper curve (high sensitivity, low false positive)
          upper <- df %>% 
            dplyr::filter(Model == model_name) %>%
            dplyr::select(x = Fp_pred_lower, 
                   y = Se_pred_upper) %>% 
            dplyr::arrange(x)
          
          # Create lower curve (low sensitivity, high false positive)
          lower <- df %>% 
            dplyr::filter(Model == model_name) %>%
            dplyr::select(x = Fp_pred_upper, 
                   y = Se_pred_lower) %>%
            dplyr::arrange(desc(x))
          
          # Combine into a closed polygon
          out <- tibble(bind_rows(upper, lower))
          n_rows <- nrow(out)
          out <- out %>% dplyr::mutate(Model = rep(model_name, n_rows))
          
}






# conf_region_colour <- "blue"
# pred_region_colour <- "blue"
# 
# stan_model_file_name = stan_model_file_name
# stan_mod_samples = stan_mod_samples
# df_true = NULL
# conf_region_colour = "blue"
# pred_region_colour = "blue"







#' R_fn_sROC_plot
#' @keywords internal
#' @export
R_fn_sROC_plot <- function( stan_model_file_name,
                            stan_mod_samples,
                            df_true = NULL,
                            conf_region_colour = "blue", 
                            pred_region_colour = "blue"
) {
            require(ggplot2)
     
            model_name <- stan_model_file_name
            # stan_mod_samples$summary(c("Se"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
            try({
              Se <- stan_mod_samples$summary(c("Se"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
              Se_pred <- stan_mod_samples$summary(c("Se_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
              ##
              Sp <- stan_mod_samples$summary(c("Sp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
              Sp_pred <- stan_mod_samples$summary(c("Sp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
              ##
              Fp <- stan_mod_samples$summary(c("Fp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
              Fp_pred <- stan_mod_samples$summary(c("Fp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
            })
            ##
            n_thr <- nrow(Se)
            # if (!(is.null(df_true))) { 
            #   df_fitted <- tibble(cbind(df_true, df_fitted)) ; df_fitted
            # }
            ##
            df_fitted_Se <- tibble(Se_median = Se$`50%`, Se_lower = Se$`2.5%`, Se_upper = Se$`97.5%`, 
                                    Se_pred_lower = Se_pred$`2.5%`, Se_pred_upper = Se_pred$`97.5%`)
            df_fitted_Sp <- tibble(Sp_median = Sp$`50%`, Sp_lower = Sp$`2.5%`, Sp_upper = Sp$`97.5%`, 
                                    Sp_pred_lower = Sp_pred$`2.5%`, Sp_pred_upper = Sp_pred$`97.5%`)
            df_fitted_Fp <- tibble(Fp_median = Fp$`50%`, Fp_lower = Fp$`2.5%`, Fp_upper = Fp$`97.5%`, 
                                    Fp_pred_lower = Fp_pred$`2.5%`, Fp_pred_upper = Fp_pred$`97.5%`)
            ##
            df_fitted <- tibble(cbind(df_fitted_Se, df_fitted_Sp, df_fitted_Fp)) %>%
                         dplyr::mutate(Model = rep(model_name, n_thr)) %>% 
                         print(n = 50)
            ##
            try({  
                if (!(is.null(df_true))) { 
                  ## df_fitted <- tibble(cbind(df_true, df_fitted)) ; df_fitted
                  df_fitted$Se_true <- df_true$Se_true
                  df_fitted$Sp_true <- df_true$Sp_true
                  df_fitted$Fp_true <- df_true$Fp_true
                }
            })
            
            ##
            polygon_Conf <- create_confidence_polygon(df_fitted, model_name = model_name)
            polygon_Pred <- create_prediction_polygon(df_fitted, model_name = model_name)
            # Jones_polygon_Conf   <- create_confidence_polygon(df_Jones,   model_name = "Jones")   ; Jones_polygon_Conf
            # Jones_polygon_Pred   <- create_prediction_polygon(df_Jones,   model_name = "Jones")   ; Jones_polygon_Pred
          
            df_all <- df_fitted
            ##
             ##try({  
            ##   if (!(is.null(df_true))) { 
              ##   df_all <- rbind(df_fitted, df_true)
            ##   } else { 
                 df_all <- df_fitted
            ##   }
            ## })
            
             {
                 plot_list <- list()
                 ##
                 ## --------- Plot 1:
                 plot_1 <-  ggplot(df_all, mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
                            geom_line(linewidth = 0.5) + 
                            geom_point(size = 3) + 
                            theme_bw(base_size = 16)
                 plot_1
                 ##
                 if (!(is.null(df_true))) { 
                     plot_1 <- plot_1 + geom_point(color = "green", size = 3, data = df_all,
                                                   mapping = aes(x = Fp_true, y = Se_true))
                     plot_1 <- plot_1 + geom_line(color = "green",  linewidth = 0.5, data = df_all, 
                                                  mapping = aes(x = Fp_true, y = Se_true))
                 }
                 ##
                 print(plot_1)
                 plot_list[[1]] <- plot_1
                
                ## --------- Plot 2:
                 plot_2 <-    ggplot(df_all, mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
                              geom_line(linewidth = 1.0) + 
                              geom_point(size = 3) + 
                              theme_bw(base_size = 16) +
                              geom_polygon(data = polygon_Conf, aes(x = x, y = y), fill = conf_region_colour, alpha = 0.50) +
                              geom_polygon(data = polygon_Pred, aes(x = x, y = y), fill = pred_region_colour, alpha = 0.25)
                             
                 
                 if (!(is.null(df_true))) { 
                   plot_2 <- plot_2 + geom_point(color = "green", size = 3, data = df_all,
                                                 mapping = aes(x = Fp_true, y = Se_true))
                   plot_2 <- plot_2 + geom_line(color = "green",  linewidth = 0.5, data = df_all, 
                                                mapping = aes(x = Fp_true, y = Se_true))
                 }
                 print(plot_2)
                 plot_list[[2]] <- plot_2
                 
             }
             
             return(list(plot_1 = plot_1,
                         plot_2 = plot_2,
                         plot_list = plot_list,
                         polygon_Conf = polygon_Conf,
                         polygon_Pred = polygon_Pred,
                         df_fitted = df_fitted))
   

}







# Create Box-Cox transformed fits for different lambda values
# Function to create Box-Cox grid values
#' box_cox_grid
#' @keywords internal
#' @export
box_cox_grid <- function(x, 
                         lambda) {
  
  if (abs(lambda) < 0.001) {
    return(log(x))
  } else {
    return((x^lambda - 1)/lambda)
  }
  
}



# 
# individual_obs_tibble <- x_individual
# group_name <- "non-diseased"
# study_index <- 1

# 
# 
# # Check what study_id values exist
# unique_studies <- unique(individual_obs_tibble$study_id)
# print("Available study IDs:")
# print(unique_studies)
# 
# # Check what group values exist
# unique_groups <- unique(individual_obs_tibble$group)
# print("Available group values:")
# print(unique_groups)



## Function for log-logistic density:
#' dloglogistic
#' @keywords internal
#' @export
dloglogistic <- function(x, 
                         alpha, 
                         beta) {
  return((alpha/beta) * (x/beta)^(alpha-1) / (1 + (x/beta)^alpha)^2)
}



## Function for Box-Cox density
#' box_cox_density
#' @keywords internal
#' @export
box_cox_density <- function(x, 
                            lambda,
                            data) {
  # Only process positive values (Box-Cox requires positive values)
  if(any(x <= 0)) return(rep(0, length(x)))
  
  # Transform the data
  transformed_data <- box_cox_transform(data, lambda)
  
  # Fit normal to transformed data
  mean_transformed <- mean(transformed_data)
  sd_transformed <- sd(transformed_data)
  
  # Calculate density with Jacobian adjustment
  density_values <- dnorm(box_cox_transform(x, lambda), 
                          mean = mean_transformed, 
                          sd = sd_transformed) * x^(lambda-1)
  
  return(density_values)
}




{
  
  ## 
  ## lambda <- -0.5:
  ##
  box_cox_neg05 <- function(x) {
    return((x^(-0.5) - 1.0)/(-0.5))
  }
  bc_neg05_jacobian <- function(x) { 
    return(x^(-1.5))
  }
  bc_neg05_density <- function(x) {
    
    unajusted_dens <- dnorm(box_cox_neg05(x), 
                            mean = mean(x), 
                            sd   = sd(x))
    
    return(unajusted_dens * bc_05_jacobian(x))
    
  }
  ## 
  ## lambda <- +0.5:
  ##
  box_cox_05 <- function(x) {
    return((x^0.5 - 1.0)/0.5)
  }
  bc_05_jacobian <- function(x) { 
    return(x^(-0.5))
  }
  bc_05_density <- function(x) {
    
    unajusted_dens <- dnorm(box_cox_05(x), 
                            mean = mean(x), 
                            sd   = sd(x))
    
    return(unajusted_dens * bc_05_jacobian(x))
    
  }
  ## 
  ## lambda > 0:
  ##
  box_cox_pos <- function(x, 
                          lambda) {
    return((x^lambda - 1.0)/lambda)
  }
  bc_pos_jacobian <- function(x,
                              lambda) { 
    return(x^(1.0 - lambda))
  }
  bc_pos_density <- function(x,
                             lambda) {
    
    trans_x <- box_cox_pos(x = x, 
                           lambda = lambda)
    
    unajusted_dens <- dnorm(trans_x, 
                            mean = mean(x), 
                            sd   = sd(x))
    
    Jacobian <- bc_pos_jacobian(x = x,
                                lambda = lambda)
    
    return(unajusted_dens * Jacobian)
    
  }
  
}



#' box_cox_transform
#' @keywords internal
#' @export
box_cox_transform <- function(x, 
                              lambda) {
  
  if(abs(lambda) < 0.001) return(log(x))
  return((x^lambda - 1)/lambda)
  
}



#' bc_density
#' @keywords internal
#' @export
bc_density <- function(x, 
                       lambda, 
                       observations) {
  
  ## Transform original data:
  transformed_data <- box_cox_transform(observations, 
                                        lambda)
  transformed_mean <- mean(transformed_data)
  transformed_sd   <- sd(transformed_data)
  
  ## Apply Jacobian adjustment:
  jacobian <- x^(lambda-1)
  
  ## Calculate density:
  dnorm(box_cox_transform(x, lambda), 
        mean = transformed_mean, 
        sd = transformed_sd) * jacobian
  
}


# study_indexes <- c(1:n_studies)
# 
# ##



##
## Function to fit distributions and create plots:
##
#' plot_distribution_fits
#' @keywords internal
#' @export
plot_distribution_fits <- function(  individual_obs_tibble, 
                                     n_studies,
                                     study_indexes = NULL, 
                                     group_name = "non-diseased") {
            
            if (group_name == "non-diseased") { 
              group_value <- 0.0
            } else { 
              group_value <- 1.0
            }
            
            if (is.null(study_indexes)) { 
              study_indexes <- seq(from = 1, to = n_studies)
            }
            
            individual_obs_tibble
            df_filtered <- dplyr::filter(individual_obs_tibble, 
                                         study_id %in% study_indexes, 
                                         group == as.numeric(group_value)) %>% 
              print(n = 100)
            
            ## Extract test values
            observations <- df_filtered$value
            
            ## Calculate bin width based on data range
            bin_width <- max(1, (max(observations) - min(observations)) / 30)
            
            ## Log-normal (need to handle zeros if present)
            obs_for_log <- observations
            if(any(obs_for_log <= 0)) obs_for_log <- observations + 1 # Shift if needed
            
            log_normal_fit <- try(MASS::fitdistr(obs_for_log, "lognormal"))
            if(!inherits(lognormal_fit, "try-error")) {
              lognormal_meanlog <- lognormal_fit$estimate["meanlog"]
              lognormal_sdlog <- lognormal_fit$estimate["sdlog"]
            } else {
              lognormal_meanlog <- mean(log(obs_for_log))
              lognormal_sdlog <- sd(log(obs_for_log))
            }
            
            ## Log-logistic (need to handle zeros if present)
            # Fit log-logistic parameters
            # For log-logistic, we'll estimate shape and scale parameters
            # using method of moments from log-transformed data
            log_data <- log(obs_for_log)
            mu_log   <- mean(log_data)
            sd_log   <- sd(log_data)
            
            ## Convert to log-logistic parameters (approximation)
            ## Shape parameter (alpha) related to the variability
            alpha <- pi/(sd_log * sqrt(3))
            ## Scale parameter (beta) related to the median
            beta <- exp(mu_log)
            # ##
            # ## Create box-cox densities:
            # ##
            # bc_neg1_density   <- bc_density(lambda = -1.0, x = observations)
            # bc_neg05_density  <- bc_density(lambda = -0.5, x = observations)
            # ## For lambda > 0:
            # bc_05_density  <- bc_density(lambda = 0.5, x = observations)
            # bc_1_density   <- bc_density(lambda = 1.0, x = observations)
            # bc_15_density  <- bc_density(lambda = 1.5, x = observations)
            # bc_2_density   <- bc_density(lambda = 2.0, x = observations)
            # # bc_25_density  <- bc_density(lambda = 2.5, x = observations)
            # bc_3_density   <- bc_density(lambda = 3.0, x = observations)
            ##
            ## Plot title:
            ##
            if (length(study_indexes == n_studies)) { 
              plot_title <- paste("Distribution of observed", group_name, "data for ALL", n_studies, "studies")
            } else { 
              plot_title <- paste("Distribution of observed", group_name, "data for all studies", study_index)
            }
            ##
            ## Create plot:
            ##
            plot <- ggplot(data.frame(value = observations), aes(x = value)) +
              ##
              ## Histogram of OBSERVED data:
              ##
              geom_histogram(aes(y = after_stat(density)), binwidth = bin_width, 
                             fill = "lightblue", color = "darkblue", alpha = 0.8) +
              ##
              ## Add density estimate:
              ##
              geom_density(color = "black", size = 2.0) +
              ##
              ## Add log-normal density (if data is positive):
              ##
              stat_function(fun = function(x) dlnorm(x,
                                                     meanlog = lognormal_meanlog, 
                                                     sdlog   = lognormal_sdlog),
                            color = "green", size = 1.0, linetype = "dashed") + 
              ##
              ## Add log-logistic density:
              ##
              stat_function(fun = function(x) dloglogistic(x, 
                                                           alpha, 
                                                           beta),
                            color = "blue", size = 1.0, linetype = "dashed") + 
              ##
              ## Add Box-Cox densities:
              ##
              stat_function(fun = function(x) bc_density(x, lambda = -1.0, observations = observations), color = "darkblue", size = 1.5) +
              stat_function(fun = function(x) bc_density(x, lambda = -0.5, observations = observations), color = "darkred",  size = 1.5) +
              ##
              stat_function(fun = function(x) bc_density(x, lambda = +0.5,  observations = observations), color = "orange", size = 1.5) +
              stat_function(fun = function(x) bc_density(x, lambda = +1.0,  observations = observations), color = "green",  size = 1.5) +
              stat_function(fun = function(x) bc_density(x, lambda = +1.5,  observations = observations), color = "purple", size = 1.5) +
              stat_function(fun = function(x) bc_density(x, lambda = +2.0,  observations = observations), color = "red",    size = 1.5) +
              ## 
              ## Formatting:
              ##
              labs(title = plot_title,
                   subtitle = paste("n =", length(observations)),
                   x = "Test Value", 
                   y = "Density") +
              theme_minimal() +
              ##
              ## Add legends:
              ##
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.90, label = "Log-normal (dashed line)",   color = "green") +
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.85, label = "Log-logistic (dashed line)", color = "blue") + 
              ## Box-cox:
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.75, label = "Box-Cox (λ=-1.0)", color = "darkblue") +
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.70, label = "Box-Cox (λ=-0.5)", color = "darkred") + 
              ## Box-cox:
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.65, label = "Box-Cox (λ=+0.5)", color = "orange") + 
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.60, label = "Box-Cox (λ=+1.0)", color = "green") + 
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.55, label = "Box-Cox (λ=+1.5)", color = "purple") + 
              annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.50, label = "Box-Cox (λ=+2.0)", color = "red")
            
            print(plot)
            
            return(plot)
  
}



















