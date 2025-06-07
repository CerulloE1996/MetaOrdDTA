


#' SD_approx_probit_to_prob
#' @export
SD_approx_probit_to_prob <- function(SD_probit_scale) { 
  
  SD_prob_scale <- SD_probit_scale * dnorm(qnorm(SD_probit_scale))
  ##
  return(signif(SD_prob_scale, 3))
  
}





#' SD_approx_ID_ord_prob_to_C_probit
#' @export
SD_approx_ID_ord_prob_to_C_probit <- function(mean_C_cutpoint_scale, SD_prob_scale) { 
  
  SD_probit_scale <- (1.0 / dnorm(mean_C_cutpoint_scale)) * SD_prob_scale
  ##
  return(signif(SD_probit_scale, 3))
  
}




#' arrange_missing_values
#' @export
arrange_missing_values <- function(data_list) {
        
        # Create a copy of the input list to avoid modifying the original
        result_list <- list(
          non_diseased = data_list[[1]],
          diseased = data_list[[2]]
        )
        
        # Process each matrix in the list
        for (i in 1:2) {
          matrix_data <- result_list[[i]]
          
          # Process each row
          for (row in 1:nrow(matrix_data)) {
            # Get current row data
            row_data <- matrix_data[row, ]
            
            # Split into regular values and missing values (-1)
            regular_values <- row_data[row_data != -1]
            missing_values <- row_data[row_data == -1]
            
            # Combine regular values first, then missing values
            new_row <- c(regular_values, missing_values)
            
            # Replace the original row
            matrix_data[row, ] <- new_row
          }
          
          result_list[[i]] <- matrix_data
        }
        
        return(result_list)
  
}



#' apply_missingness_pattern
#' @export
apply_missingness_pattern <- function( x_complete, 
                                       x_pattern, 
                                       enforce_consecutive_missingness = FALSE) {
  
  
  # Check if inputs are valid
  if (!is.list(x_complete) || !is.list(x_pattern) || length(x_complete) != 2 || length(x_pattern) != 2) {
    stop("Both inputs must be lists with exactly 2 elements (non-diseased and diseased matrices)")
  }
  
  # Create copies of the complete data to modify
  result <- list(
    non_diseased = x_complete[[1]],
    diseased = x_complete[[2]]
  )
  
  # Get dimensions of the first complete matrix to use as reference
  n_rows_complete <- nrow(x_complete[[1]])
  n_cols_complete <- ncol(x_complete[[1]])
  
  # Calculate the missingness pattern based on the first pattern matrix
  missingness_indices <- list()
  pattern_matrix <- x_pattern[[1]]
  n_rows_pattern <- nrow(pattern_matrix)
  n_cols_pattern <- ncol(pattern_matrix)
  
  # Find missing indices for each column
  for (col in 2:min(n_cols_complete, n_cols_pattern)) {
    # Find rows with -1 in the pattern matrix
    missing_rows_in_pattern <- which(pattern_matrix[, col] == -1)
    
    if (length(missing_rows_in_pattern) > 0) {
      # Calculate proportion of missing values in this column of pattern
      missing_proportion <- length(missing_rows_in_pattern) / n_rows_pattern
      
      # Calculate how many rows should be missing in the complete matrix
      n_missing_in_complete <- round(n_rows_complete * missing_proportion)
      
      if (n_missing_in_complete > 0) {
        # Randomly select rows to make missing
        rows_to_make_missing <- sample(1:n_rows_complete, n_missing_in_complete)
        
        # Store these indices
        missingness_indices[[col]] <- rows_to_make_missing
      } else {
        missingness_indices[[col]] <- integer(0)
      }
    } else {
      missingness_indices[[col]] <- integer(0)
    }
  }
  
  # Apply the same missingness pattern to both populations
  for (pop_idx in 1:2) {
    # For each column in the complete matrix
    for (col in 2:min(n_cols_complete, n_cols_pattern)) {
      if (length(missingness_indices[[col]]) > 0) {
        # Set values to -1 using the same indices
        result[[pop_idx]][missingness_indices[[col]], col] <- -1
      }
    }
    
    # Optionally apply the consecutive missingness rule
    if (enforce_consecutive_missingness) {
      # If there's a -1 in column j, all columns to the right should also be -1
      for (row in 1:nrow(result[[pop_idx]])) {
        missing_found <- FALSE
        for (col in 2:ncol(result[[pop_idx]])) {
          if (missing_found || result[[pop_idx]][row, col] == -1) {
            missing_found <- TRUE
            result[[pop_idx]][row, col] <- -1
          }
        }
      }
    }
  }
  
  return(result)
  
}














#' generate_ordinal_dirichlet_prior
#' @export
generate_ordinal_dirichlet_prior <- function(n_cat, 
                                             max_alpha = 5, 
                                             min_alpha = 1,
                                             left_skew_factor = 0.3) {
        
        # n_cat: Number of categories in the ordinal scale
        # max_alpha: Maximum concentration parameter (for early categories)
        # min_alpha: Minimum concentration parameter (for later categories)
        # left_skew_factor: Controls how quickly values decrease (lower = more left-skewed)
        
        # Generate positions from 0 to 1
        positions <- seq(0, 1, length.out = n_cat)
        
        # Apply a transformation to create the left skew
        # Using a power function where lower values of left_skew_factor create more skew
        transformed_positions <- positions^left_skew_factor
        
        # Normalize to [0,1] range
        transformed_positions <- transformed_positions / max(transformed_positions)
        
        # Map to alpha parameter range
        alpha_values <- max_alpha - (max_alpha - min_alpha) * transformed_positions
        
        # Round to create cleaner values (optional)
        alpha_values <- round(alpha_values, 1)
        
        return(alpha_values)
        
}









#' generate_dual_ordinal_dirichlet_priors
#' @export
generate_dual_ordinal_dirichlet_priors <- function(n_cat, 
                                                   nd_peak = 0.15,          # Where non-diseased peaks (as fraction)
                                                   d_peak = 0.30,           # Where diseased peaks (as fraction)
                                                   nd_max_alpha = 5.0,      # Max alpha for non-diseased
                                                   d_max_alpha = 5.0,       # Max alpha for diseased
                                                   nd_min_alpha = 1.0,      # Min alpha for non-diseased
                                                   d_min_alpha = 1.0,
                                                   smoothness = 0.15) {     # Higher = smoother distribution
  
        # Create sequence for category positions
        positions <- seq(0, 1, length.out = n_cat)
        
        # Generate smoother non-diseased alphas
        nd_spread <- smoothness * 1.2  # Slightly wider for non-diseased
        nd_alpha_values <- nd_min_alpha + (nd_max_alpha - nd_min_alpha) * 
          exp(-((positions - nd_peak)^2) / (2 * nd_spread^2))
        
        # Generate smoother diseased alphas
        d_spread <- smoothness * 1.5   # Even wider for diseased
        d_alpha_values <- d_min_alpha + (d_max_alpha - d_min_alpha) * 
          exp(-((positions - d_peak)^2) / (2 * d_spread^2))
        
        # Round for cleaner values
        nd_alpha_values <- round(nd_alpha_values, 1)
        d_alpha_values <- round(d_alpha_values, 1)
        
        return(list(
          non_diseased = nd_alpha_values,
          diseased = d_alpha_values
        ))
  
}



#' divide_studies
#' @export
divide_studies <- function( n_studies, 
                            n_groups) {
  
        # Calculate base size for each group
        base_size <- floor(n_studies / n_groups)
        
        # Calculate remainder
        remainder <- n_studies %% n_groups
        
        # Create list to store group vectors
        groups <- list()
        
        # Start index
        start_idx <- 1
        
        # For each group
        for (i in 1:n_groups) {
            
            # Determine group size (add 1 to base_size if there's still remainder)
            group_size <- base_size + (if (i <= remainder) 1 else 0)
            
            # Calculate end index
            end_idx <- start_idx + group_size - 1
            
            # Store group vector
            groups[[i]] <- start_idx:end_idx
            
            # Update start index for next group
            start_idx <- end_idx + 1
          
        }
        
        return(groups)
        
}
















#' R_fn_apply_thr_missingness
#' @export
R_fn_apply_thr_missingness_to_test <- function(  data_list_of_mat, 
                                                 thr_combo_vec = NULL, 
                                                 apply_missings_to_these_rows
) {

         data_list_of_mat_original <- data_list_of_mat
          
         for (c in 1:2) {
              data_list_of_mat[[c]] <-  R_fn_apply_thr_missingness_to_mat( data_matrix = data_list_of_mat_original[[c]],
                                                                           thr_combo_vec = thr_combo_vec,
                                                                           apply_missings_to_these_rows = apply_missings_to_these_rows)
         }
          
         return(data_list_of_mat)
  
}





 


#' R_fn_apply_missingness_pattern_PHQ_9
#' @export
R_fn_apply_missingness_pattern_PHQ_9 <- function( x_PHQ9, 
                                                  case_list, 
                                                  seed = 123) {
  
          # set.seed(123)
  
          # Create a copy of the original data
          x_PHQ_missing <- x_PHQ9
          
          # Loop through each case
          for (case_name in names(case_list)) {
                  
                  case <- case_list[[case_name]]
                  studies <- case$studies
                
                  if (case_name == "case_3") {
                    
                          # For case 3, we need to generate random thresholds for each study:
                          for (study_idx in studies) {
                            # Select a random threshold between 8 and 15
                            random_thr <- sample(8:15, 1)
                            
                            # Create mask where everything is -1 except the random threshold
                            for (group in 1:2) { # Loop over diseased and non-diseased groups
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (col != random_thr) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  } else if (case_name == "case_4") {
                    
                          ## For case 4, we randomly select one of the threshold combinations for each study:
                          # Extract all threshold combination vectors
                          thr_combo_vecs <- list()
                          for (i in 1:12) {
                            vec_name <- paste0("thr_combo_vec_", i)
                            if (vec_name %in% names(case)) {
                              thr_combo_vecs[[i]] <- case[[vec_name]]
                            }
                          }
                          
                          # Apply a randomly selected threshold combination to each study
                          for (study_idx in studies) {
                            # Randomly select one of the threshold combinations
                            selected_combo_idx <- sample(1:length(thr_combo_vecs), 1)
                            selected_thresholds <- thr_combo_vecs[[selected_combo_idx]]
                            
                            # Create mask where everything is -1 except the selected thresholds
                            for (group in 1:2) {
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (!(col %in% selected_thresholds)) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  } else if (case_name %in% c("case_1", "case_2", "case_5")) {
                    
                          ## For cases 1, 2, and 5, the threshold pattern is the same for all studies in the case:
                          ## Extract the threshold vector (it's called thr_combo_vec_1 for these cases):
                          thresholds <- case$thr_combo_vec_1
                          
                          ## Apply the missingness pattern to all studies in this case:
                          for (study_idx in studies) {
                            for (group in 1:2) {
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (!(col %in% thresholds)) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  }
            
          }
          
          return(x_PHQ_missing)
  
}








#' convert_to_aggregate_counts
#' @export
convert_to_aggregate_counts <- function(y_list, 
                                        n_studies, 
                                        n_tests_inc_ref,
                                        n_thr) {
  
          ## Add BINARY ref test n_thr (=1) to "n_thr":
          n_thr <- c(1, n_thr)
          ##
          n_tests <- n_tests_inc_ref

          n_total_nd <- n_total_d <- numeric(n_studies)
          x_nd_list <- x_d_list <- list()
          Se_per_study_list <- Sp_per_study_list <- list()
          
          ## Initialise lists:
          n_index_tests <- n_tests - 1
          for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
            
            n_thr_t <- n_thr[t]
            n_cat_t <- n_thr_t + 1
            ##
            x_nd_list[[t - 1]] <- matrix(NA, n_studies, n_cat_t)
            x_d_list[[t - 1]]  <- matrix(NA, n_studies, n_cat_t)
            ##
            Se_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t + 1)
            Sp_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t + 1)
            
          }
          
          for (s in 1:n_studies) {
            
                    study_data <- y_list[[s]]
                    disease_status <- study_data[, 1] # Col 1 is disease status / reference test
                    
                    ## Get n_total_d and n_total_nd:
                    n_total_d[s]  <- sum(disease_status == 1)
                    n_total_nd[s] <- sum(disease_status == 0)
                
            for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
              
                    n_thr_t <- n_thr[t]
                    n_cat_t <- n_thr_t + 1
                    index_test_results   <- study_data[, t] ## - 1 ## minus one so {0, .., 27} for e.g. PHQ-9.
                     
                    print(paste("min = ", min(index_test_results)))
                    print(paste("max = ", max(index_test_results)))
                    
                    # x_d_list[[t - 1]][s, 1]  <- n_total_d[s]
                    # x_nd_list[[t - 1]][s, 1] <- n_total_nd[s]
                    # Get counts for each threshold
                    for (k in 1:(n_thr_t + 1)) {
                          ##
                          ## True-positives (TP's):
                          ##
                          x_d_list[[t - 1]][s, k]  <- sum(index_test_results[disease_status == 1] >= k - 1) # Pr(testing POSITIVE at threshold k)
                          Se_per_study_list[[t - 1]][s, k] <-   x_d_list[[t - 1]][s, k]  / n_total_d[s]
                          ##
                          ## False-positives (FP's):
                          ##
                          x_nd_list[[t - 1]][s, k] <- sum(index_test_results[disease_status == 0] >= k - 1) # Pr(testing POSITIVE at threshold k)
                          Sp_per_study_list[[t - 1]][s, k] <-   x_nd_list[[t - 1]][s, k]  / n_total_d[s]
                          
                    }
                    # for (k in 1:(n_thr_t + 1)) {
                    # 
                    #     Se_per_study_list[[t - 1]][s, k] <-        ( sum(index_test_results[disease_status == 1] >= k - 1) / n_total_d[s])
                    #     Sp_per_study_list[[t - 1]][s, k] <- 1.0 -  ( sum(index_test_results[disease_status == 0] >= k - 1)/ n_total_nd[s])
                    #   
                    # }
            }
            
          }

          
          
          return(list(
            n_total_nd = n_total_nd,
            n_total_d  = n_total_d, 
            x_nd_list = x_nd_list,
            x_d_list  = x_d_list,
            Se_per_study_list = Se_per_study_list,
            Sp_per_study_list = Sp_per_study_list
          ))
  
}








#' create_test_in_study_for_NMA - function to create "indicator_test_in_study" indicator matrix
#' @export
create_test_in_study_for_NMA <- function(   seed,
                                            n_studies, 
                                            n_tests_inc_ref, 
                                            prob_present = 0.75, 
                                            min_index_tests_per_study = 1
) {
  
  set.seed(seed, kind = "Mersenne-Twister")
  
  n_tests <- n_tests_inc_ref
  
  ## At least 1 index tests per study -> at least 2 total tests per study:
  min_tests_per_study <- min_index_tests_per_study + 1
  
  ## Initialize the matrix:
  indicator_index_test_in_study <- matrix(0, nrow = n_studies, ncol = n_tests)
  
  ## Ref test is ALWAYS included: 
  reference_test_prob <- 1.0
  
  ## For each study, decide which tests are present:
  for (s in 1:n_studies) {
    
    ## Give reference test higher probability of being included
    present_probs <- rep(prob_present, n_tests)
    ## Assuming test 1 is ref test:
    present_probs[1] <- reference_test_prob  
    
    ## Randomly determine which tests are present (1 = present, 0 = absent)
    test_presence <- rbinom(n_tests, 1, present_probs)
    
    ## Ensure at least min_tests_per_study tests are present
    while (sum(test_presence) < min_tests_per_study) {
      ## If not enough tests, add one more randomly
      missing_tests <- which(test_presence == 0)
      if (length(missing_tests) > 0) {
        add_test <- sample(missing_tests, 1)
        test_presence[add_test] <- 1
      } else {
        break  # All tests already included
      }
    }
    
    indicator_index_test_in_study[s, ] <- test_presence
    
  }
  
  ## Optional: Ensure each test appears in at least one study:
  for (t in 1:n_tests) {
    
    if (sum(indicator_index_test_in_study[, t]) == 0) {
      # If a test doesn't appear in any study, add it to a random study
      random_study <- sample(1:n_studies, 1)
      indicator_index_test_in_study[random_study, t] <- 1
    }
    
  }
  
  indicator_index_test_in_study <- indicator_index_test_in_study[1:n_studies, 2:n_tests]
  
  ## Count number of (observed / non-missing): 
  n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
  
  return(list( indicator_index_test_in_study = indicator_index_test_in_study,
               n_index_tests_per_study = n_index_tests_per_study))
  
}















#' apply_thr_missingness
#' @keywords internal
#' @export
apply_thr_missingness <- function(  agg_data_cumulative, 
                                    studies_subset_vec,
                                    missing_thr_subset_vec,
                                    missing_indicator = -1) { 
  
  agg_data_cumulative$x_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  agg_data_cumulative$x_non_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  
  return(agg_data_cumulative)
  
}
















#' convert_cumulative_to_category
#' @keywords internal
#' @export
convert_cumulative_to_category <- function(cumulative_matrix,
                                           missing_indicator = -1) {
  
  # Get dimensions and create output matrix with same dimensions plus one column
  n_rows <- nrow(cumulative_matrix)
  n_cols <- ncol(cumulative_matrix)
  category_matrix <- matrix(-1, nrow = n_rows, ncol = n_cols + 1)
  
  for (i in 1:n_rows) {
    # Get current row
    row_data <- cumulative_matrix[i, ]
    
    # First category is the first cumulative count
    category_matrix[i, 1] <- row_data[1]
    
    # Second category is second column minus first column
    if (n_cols >= 2) {
      if (!is.na(row_data[1]) && !is.na(row_data[2])) {
        category_matrix[i, 2] <- row_data[2] - row_data[1]
      }
    }
    
    # Process remaining categories by taking differences
    for (j in 2:(n_cols-1)) {
      if (!is.na(row_data[j]) && !is.na(row_data[j+1])) {
        category_matrix[i, j+1] <- row_data[j+1] - row_data[j]
      }
    }
    
    # Last category is the total minus the last cumulative count
    # (Assuming we want all categories to sum to the total)
    if (!is.na(row_data[n_cols])) {
      total <- row_data[n_cols]
      category_matrix[i, n_cols+1] <- 0  # Default to zero
      
      # If we want a true non-zero last category, we would need the actual total N
      # which isn't directly provided in this data format
    }
  }
  
  return(category_matrix[, 1:n_cols])
  
}

















#' categorical_to_individual
#' @keywords internal
#' @export
categorical_to_individual <- function( categorical_matrix,
                                       binary_disease_indicator,
                                       missing_indicator = -1) {
  
  # Initialize list to store results
  results <- list()
  n_studies <- nrow(categorical_matrix)
  n_categories <- ncol(categorical_matrix)
  
  for (study in 1:n_studies) {
    # Get counts for current study
    study_counts <- categorical_matrix[study, ]
    
    # Initialize vectors for this study
    study_id <- vector()
    category_values <- vector()
    
    # For each category
    for (cat in 1:n_categories) {
      count <- study_counts[cat]
      
      if (count == missing_indicator) {
        # For missing data, create one observation with missing_indicator
        study_id <- c(study_id, study)
        category_values <- c(category_values, missing_indicator)
      } else if (count > 0) {
        # For non-missing data, replicate the category value count times
        study_id <- c(study_id, rep(study, count))
        category_values <- c(category_values, rep(cat, count))
      }
    }
    
    # Create data frame for this study
    study_data <- data.frame(
      study_id = study_id,
      group = binary_disease_indicator,  # 1 for diseased, 0 for non-diseased
      value = category_values
    )
    
    results[[study]] <- study_data
  }
  
  # Combine all studies into one data frame
  final_data <- do.call(rbind, results)
  rownames(final_data) <- NULL  # Reset row names
  
  return(tibble(final_data))
  
}










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
 




# individual_obs_tibble <- x_individual
# group_name <- "non-diseased"
# study_index <- 1
# 
# # 
# # 
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

##
## Function to fit distributions and create plots:
##
#' plot_distribution_fits
#' @keywords internal
#' @export
plot_distribution_fits <- function(individual_obs_tibble, 
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




















