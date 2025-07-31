


 
# 

# Helper to handle NULL values
`%||%` <- function(x, y) if (is.null(x)) y else x



 



 
#' simulate_covariates_preserving_missingness_FIXED
#' @keywords internal
#' @export
simulate_covariates_preserving_missingness_FIXED <- function(original_X, 
                                                             n_studies,
                                                             indicator_index_test_in_study,
                                                             covariates_same_between_classes = NULL,
                                                             add_noise = TRUE,
                                                             noise_sd_prop = 0.25,
                                                             flip_prob = 0.15) {
  
          n_classes <- length(original_X)
          n_tests <- length(original_X[[1]])
          
          # Auto-detect number of covariates
          n_covariates <- ncol(original_X[[1]][[1]])
          
          # If not specified, assume all covariates same across disease
          if (is.null(covariates_same_between_classes)) {
            covariates_same_between_classes <- rep(TRUE, n_covariates)
          }
          
          # Auto-detect covariate types by examining the data
          covariate_types <- character(n_covariates)
          
          for (col in 1:n_covariates) {
            # Collect all non-zero values across all tests/classes
            all_values <- c()
            for (c in 1:n_classes) {
              for (t in 1:n_tests) {
                col_vals <- original_X[[c]][[t]][, col]
                # Get non-missing values (assuming 0 in intercept column means missing)
                if (col == 1) {  # Intercept column
                  non_missing_rows <- which(col_vals == 1)
                } else {
                  non_missing_rows <- which(original_X[[c]][[t]][, 1] == 1)  # Where intercept = 1
                }
                if (length(non_missing_rows) > 0) {
                  all_values <- c(all_values, col_vals[non_missing_rows])
                }
              }
            }
            
            # Determine type
            unique_vals <- unique(all_values)
            if (length(unique_vals) == 0) {
              covariate_types[col] <- "unknown"
            } else if (all(unique_vals %in% c(0, 1))) {
              covariate_types[col] <- "binary"
            } else {
              covariate_types[col] <- "continuous"
            }
          }
          
          cat("Detected covariate types:\n")
          for (i in 1:n_covariates) {
            cat(sprintf("  Column %d (%s): %s\n", 
                        i, 
                        colnames(original_X[[1]][[1]])[i], 
                        covariate_types[i]))
          }
          
          # Identify available rows per test
          available_rows_per_test <- list()
          for (t in 1:n_tests) {
            # Rows where intercept = 1 (indicating non-missing)
            available_rows_per_test[[t]] <- which(original_X[[1]][[t]][, 1] == 1)
          }
          
          # ========== HERE'S THE FIX: Sample ONE row per study ==========
          sampled_row_per_study <- rep(NA, n_studies)
          
          # For each study, sample from rows that are available in ALL its tests
          for (s in 1:n_studies) {
            # Which tests does this study have?
            tests_for_study <- which(indicator_index_test_in_study[s, ])
            
            if (length(tests_for_study) > 0) {
              # Find rows available in ALL tests for this study
              common_available_rows <- available_rows_per_test[[tests_for_study[1]]]
              
              if (length(tests_for_study) > 1) {
                for (t in tests_for_study[-1]) {
                  common_available_rows <- intersect(common_available_rows, available_rows_per_test[[t]])
                }
              }
              
              # Sample ONE row for this study
              if (length(common_available_rows) > 0) {
                sampled_row_per_study[s] <- sample(common_available_rows, 1)
              } else {
                warning(sprintf("No common available rows for study %d across tests %s", 
                                s, paste(tests_for_study, collapse=", ")))
              }
            }
          }
          
          # Initialize output
          simulated_X <- list()
          
          # Fill in data
          for (c in 1:n_classes) {
            simulated_X[[c]] <- list()
            
            for (t in 1:n_tests) {
              # Initialize with zeros
              X_simulated <- matrix(0, nrow = n_studies, ncol = n_covariates)
              colnames(X_simulated) <- colnames(original_X[[c]][[t]])
              
              # Fill in data
              studies_with_this_test <- which(indicator_index_test_in_study[, t])
              
              for (s in studies_with_this_test) {
                if (!is.na(sampled_row_per_study[s])) {
                  # Use the SAME row for this study across ALL tests
                  row_to_copy <- sampled_row_per_study[s]
                  X_simulated[s, ] <- original_X[[c]][[t]][row_to_copy, ]
                  
                  # Add noise only once (when c=1) and only if the row is valid
                  if (add_noise && c == 1 && X_simulated[s, 1] == 1) {
                    
                    # Process each covariate based on its type
                    for (col in 2:n_covariates) {  # Skip intercept
                      
                      if (covariate_types[col] == "continuous") {
                        # Get all non-zero values for this column to calculate SD
                        col_values <- original_X[[c]][[t]][available_rows_per_test[[t]], col]
                        if (length(col_values) > 1 && sd(col_values) > 0) {
                          noise_sd <- noise_sd_prop * sd(col_values)
                          X_simulated[s, col] <- X_simulated[s, col] + rnorm(1, 0, noise_sd)
                        }
                        
                      } else if (covariate_types[col] == "binary") {
                        # Flip with probability
                        if (runif(1) < flip_prob) {
                          X_simulated[s, col] <- 1 - X_simulated[s, col]
                        }
                      }
                      # Unknown types are left as-is
                    }
                  }
                }
              }
              
              simulated_X[[c]][[t]] <- X_simulated
            }
          }
          
          # Copy from class 1 to class 2 for covariates that should be same
          if (!all(covariates_same_between_classes)) {
            for (t in 1:n_tests) {
              studies_with_this_test <- which(indicator_index_test_in_study[, t])
              for (s in studies_with_this_test) {
                for (col in 1:n_covariates) {
                  if (covariates_same_between_classes[col]) {
                    simulated_X[[2]][[t]][s, col] <- simulated_X[[1]][[t]][s, col]
                  }
                }
              }
            }
          } else {
            # All same - just copy everything
            simulated_X[[2]] <- simulated_X[[1]]
          }
          
          return(simulated_X)
  
}










#' apply_test_indicators_to_X
#' @export
apply_test_indicators_to_X <- function(X,
                                       indicator_matrix) {
  
          # X structure: X[[disease_status]][[test]]
          # indicator_matrix: rows = studies, cols = tests
          
          # Get dimensions
          n_disease_states <- length(X)
          n_tests <- length(X[[1]])
          n_studies <- nrow(indicator_matrix)
          
          # Create new X with same structure
          X_new <- X
          
          # For each disease state
          for (d in 1:n_disease_states) {
            # For each test
            for (t in 1:n_tests) {
              # Get the indicator column for this test
              test_available <- indicator_matrix[, t]
              
              # Multiply each row of X by the indicator
              # If test not available in study, all covariates become 0
              for (s in 1:n_studies) {
                X_new[[d]][[t]][s, ] <- X[[d]][[t]][s, ] * test_available[s]
              }
            }
          }
          
          return(X_new)
  
}

# # Usage:
# X_corrected <- apply_test_indicators_to_X(X, indicator_index_test_in_study)
# 
# # Check the result for test 1:
# X_corrected[[1]][[1]][1:10, ]  # Should show 0s where indicator_index_test_in_study[,1] == 0
# 
# # To verify it worked correctly:
# # For test 1, rows 2, 5, 6, 8, 12, etc. should be 0
# # For test 2, rows 5, 13, 15, 16, etc. should be 0
# # etc.
# 
# 



#' R_fn_sim_NMA_data_with_covariates_varying
#' @keywords internal
#' @export
R_fn_sim_NMA_data_with_covariates_varying <- function(  seed = 123,
                                                        n_studies_per_test = c(25, 20, 30, 25, 22),  # Vector of studies per test (including reference)
                                                        N_per_study_mean = 500,
                                                        N_per_study_SD = 500,
                                                        assume_perfect_GS = TRUE,
                                                        true_Mean_prev = 0.2,
                                                        true_SD_probit_prev = 0.25,
                                                        simulate_covariates = TRUE,
                                                        covariate_settings,
                                                        original_cov_data = NULL,
                                                        missing_value_covs = -999,
                                                        # study_overlap_prop = 0.333,  # Proportion of studies that have multiple tests
                                                        require_index_test = TRUE  # If TRUE, every study must have at least one index test
) {
  
          random_C <- TRUE
          hetero_sigma <- FALSE
          ##
          study_overlap_prop <- 0.333 # Proportion of studies that have multiple tests
  
          set.seed(seed, kind = "L'Ecuyer-CMRG")
          
          # Setup parameters
          n_binary_tests <- 1
          n_ordinal_tests <- 4
          n_index_tests <- 4
          n_tests <- n_binary_tests + n_ordinal_tests
          
          # Validate input
          if (length(n_studies_per_test) != n_tests) {
            stop("n_studies_per_test must have length equal to number of tests")
          }
          
          # Calculate total unique studies (with overlap)
          # max_studies <- max(n_studies_per_test)
          # n_total_studies <- round(max_studies * (1 + (1 - study_overlap_prop)))
          ##
          n_total_studies <- n_studies_per_test[1]
          ##
          # Create study-test availability matrix
          indicator_test_in_study <- matrix(FALSE, nrow = n_total_studies, ncol = n_tests)
          
          # First, ensure all studies have the reference test
          indicator_test_in_study[, 1] <- TRUE
          
          # Correct test indices:
          # Column 2 = GAD-2 (Test1 in your output)
          # Column 3 = GAD-7 (Test2 in your output)
          # Column 4 = HADS (Test3 in your output)
          # Column 5 = BAI (Test4 in your output)
         
         
          # ========== Extract real test patterns ==========
          # Get test patterns from original data (only for index tests 2-5)
          real_test_patterns <- matrix(0, nrow = nrow(original_cov_data$X[[1]][[1]]), ncol = n_tests - 1)  # Only index tests
          
          # Extract patterns for index tests only (tests 2-5)
          for (t in 2:n_tests) {
            valid_studies <- which(original_cov_data$X[[1]][[t-1]][, "intercept"] == 1)  # t-1 because X doesn't include ref test
            real_test_patterns[valid_studies, t-1] <- 1
          }
          
          # Get unique patterns and their frequencies
          unique_patterns <- unique(real_test_patterns)
          pattern_counts <- table(apply(real_test_patterns, 1, paste, collapse = ""))
          
          cat("Available test patterns in real data:\n")
          for (i in 1:nrow(unique_patterns)) {
            pattern <- unique_patterns[i, ]
            tests_in_pattern <- which(pattern == 1) + 1  # +1 because these are tests 2-5
            pattern_string <- paste(pattern, collapse = "")
            count <- pattern_counts[pattern_string]
            cat(sprintf("  Pattern %d: Tests %s (n=%d studies)\n", 
                        i, paste(tests_in_pattern, collapse = ","), count))
          }
          
          # Calculate total studies needed
          n_total_studies <- max(n_studies_per_test)
          
          # Initialize indicator matrix
          indicator_test_in_study <- matrix(FALSE, nrow = n_total_studies, ncol = n_tests)
          
          # First ensure reference test for all (ALWAYS 1)
          indicator_test_in_study[, 1] <- TRUE
          
          # Sample from real patterns
          sampled_patterns <- sample(1:nrow(unique_patterns), 
                                     n_total_studies, 
                                     replace = TRUE,
                                     prob = as.numeric(pattern_counts))
          
          # Assign patterns (only for index tests 2-5)
          for (s in 1:n_total_studies) {
            pattern <- unique_patterns[sampled_patterns[s], ]
            indicator_test_in_study[s, 2:n_tests] <- as.logical(pattern)
          }
          
          # Fine-tune to match target numbers (only for index tests)
          for (t in 2:n_tests) {
            current_count <- sum(indicator_test_in_study[, t])
            target_count <- n_studies_per_test[t]
            
            if (current_count < target_count) {
              # Need more studies with this test
              candidates <- which(indicator_test_in_study[, t] == FALSE)
              
              # Check which could have this test based on real patterns
              valid_candidates <- c()
              for (cand in candidates) {
                current_pattern <- indicator_test_in_study[cand, 2:n_tests]  # Only index tests
                test_pattern <- current_pattern
                test_pattern[t-1] <- TRUE  # t-1 because pattern is only for index tests
                pattern_string <- paste(as.numeric(test_pattern), collapse = "")
                if (pattern_string %in% names(pattern_counts)) {
                  valid_candidates <- c(valid_candidates, cand)
                }
              }
              
              n_to_add <- min(length(valid_candidates), target_count - current_count)
              if (n_to_add > 0) {
                add_studies <- sample(valid_candidates, n_to_add)
                indicator_test_in_study[add_studies, t] <- TRUE
              }
              
            } else if (current_count > target_count) {
              # Too many studies with this test
              has_this_test <- which(indicator_test_in_study[, t])
              
              valid_removals <- c()
              for (cand in has_this_test) {
                current_pattern <- indicator_test_in_study[cand, 2:n_tests]  # Only index tests
                test_pattern <- current_pattern
                test_pattern[t-1] <- FALSE  # t-1 because pattern is only for index tests
                pattern_string <- paste(as.numeric(test_pattern), collapse = "")
                if (pattern_string %in% names(pattern_counts)) {
                  valid_removals <- c(valid_removals, cand)
                }
              }
              
              n_to_remove <- min(length(valid_removals), current_count - target_count)
              if (n_to_remove > 0) {
                remove_studies <- sample(valid_removals, n_to_remove)
                indicator_test_in_study[remove_studies, t] <- FALSE
              }
            }
          }
          
          # No need for the filtering section since we're using real patterns
          n_studies <- nrow(indicator_test_in_study)
          ##
          indicator_index_test_in_study <- indicator_test_in_study[, 2:ncol(indicator_test_in_study)]
          
          # Initialize parameter storage
          true_params_list <- list()
          ##
          # # Copy parameter arrays (same as before)
          # if (hetero_sigma) {
          #   ## not doing this as DSM - low CV's for sigma found.
          # } else { 
               if (random_C) {
                     true_params_list <- list()
                     # ##
                     # Extract from full model
                     beta_mu_full_original <- structure(c(-0.736, -1.425, -0.86, -1.686, 0.356, -0.247, 0.233, 
                                                 -0.349, 0.101, 0.132, -0.025, 0.317, 0.244, 0.312, -0.224, 0.025, 
                                                 0.062, -0.015, -0.132, 0.001, -0.139, -0.357, -0.25, 0.004, -0.031, 
                                                 -0.15, 0.444, 0, -0.221, -0.099, 0.053, 0.001, -0.042, -0.088, 
                                                 0.251, -0.19, -0.142, 0.011, 0.029, -0.169, 0.083, 0.155, -0.065, 
                                                 0.415, 0.225, 0.267, -0.189, 0.118, 0.475, 0.328, -0.135, 0.002, 
                                                 -0.348, -0.367, -0.129, 0.001),
                                               dim = c(4L, 2L, 7L))
                     ##
                     ## Set beta_mu initially to full model one:
                     ##
                     true_params_list$beta_mu_array <- beta_mu_full_original
                     
                     # For prev_GAD and ref_test: ensure AT LEAST as large as Model 6
                     beta_mu_model6 <-    structure(c(-1.057, -1.361, -1.131, -1.802, 0.149, -0.028, 0.288, 
                                                      -0.53, 0.114, 0.155, -0.014, 0.278, 0.25, 0.341, -0.334, 0, 0.018, 
                                                      0.117, -0.047, 0.368, 0.182, 0.362, -0.304, 0.135, 0.363, 0.427, 
                                                      -0.187, 0.003, -0.5, -0.501, 0.196, 0), 
                                                    dim = c(4L, 2L, 4L))
                     
                     # Take the LARGER of the two (in absolute value)
                     for (t in 1:n_index_tests) {
                           for (c in 1:2) {
                             
                                 # Prevalence (column 2 in both beta_mu_full_original and beta_mu_model6)
                                 if (abs(beta_mu_full_original[t, c, 2]) < abs(beta_mu_model6[t, c, 2])) {
                                     true_params_list$beta_mu_array[t, c, 2] <- beta_mu_model6[t, c, 2]
                                 }
                                 
                                 # Ref test type (columns 6-7 in beta_mu_full_original, cols 3-4 in beta_mu_model6 so minus 3)
                                 for (ref_test_cols_in_X_FULL in c(6:7)) {
                                     ref_test_cols_in_X_MR_Model_6 <- ref_test_cols_in_X_FULL - 3
                                     if (abs(beta_mu_full_original[t, c, ref_test_cols_in_X_FULL]) < abs(beta_mu_model6[t, c, ref_test_cols_in_X_MR_Model_6 - 1])) {
                                       true_params_list$beta_mu_array[t, c, ref_test_cols_in_X_FULL] <- beta_mu_model6[t, c, ref_test_cols_in_X_MR_Model_6 - 1]
                                     }
                                 }
                             
                           }
                     }
                     
                     # For RoB and study_setting: just make them smaller
                     RoB_cols_in_X_FULL <- 3
                     true_params_list$beta_mu_array[, , RoB_cols_in_X_FULL] <- beta_mu_full_original[, , RoB_cols_in_X_FULL] * 0.5 
                     ##
                     study_setting_cols_in_X_FULL <- c(4:5)
                     true_params_list$beta_mu_array[, , study_setting_cols_in_X_FULL] <- beta_mu_full_original[, , study_setting_cols_in_X_FULL] * 0.5 
                     
                     
                     # Print what's actually in each position
                     {
                         print("Beta position 2 (should be prev):")
                         print(true_params_list$beta_mu_array[1, , 2])
                         
                         print("Beta position 3 (should be RoB):")
                         print(true_params_list$beta_mu_array[1, , 3])
                         
                         print("Beta positions 4-5 (should be study_setting):")
                         print(true_params_list$beta_mu_array[1, , 4:5])
                         
                         print("Beta positions 6-7 (should be ref_test):")
                         print(true_params_list$beta_mu_array[1, , 6:7])
                         ##
                         print(paste("beta_mu_full_original = :"))
                         print(beta_mu_full_original)
                         ##
                         print(paste("final = :"))
                         print(beta_mu_model6)
                         ##
                         print(paste("final beta_mu array to use for simulation dummy example:"))
                         print(true_params_list$beta_mu_array)
                     }
                     
                     
                    } else  { 
          }
          ##
          # true_params_list$beta_sigma_array <- structure(c(0.25667275, 0.280819125, 0.25048075, 0.26706375, 
          #                                                      0.325171375, 0.341804375, 0.243068125, 0.310942), dim = c(4L, 
          #                                                                                                               2L))
          # # Hierarchical parameters for beta_sigma
          # if (hetero_sigma) {
          #     if (random_C) {
          #       if (simulate_covariates) {
          #           true_params_list$log_beta_sigma_MU <- c(-1.491510, -1.575652)
          #           true_params_list$log_beta_sigma_SD <- c(0.1489826, 0.2187959)
          #       } else {
          #           true_params_list$log_beta_sigma_MU <- c(-1.314094, -1.204109)
          #           true_params_list$log_beta_sigma_SD <- c(0.1577823, 0.1811997)
          #       }
          #     } else { 
          #       true_params_list$log_beta_sigma_MU <- c(-1.326167, -1.255996)
          #       true_params_list$log_beta_sigma_SD <- c(0.1571329, 0.1686700)
          #     }
          # } else { 
          #     # if (random_C) {
          #     #   if (simulate_covariates) {
          #     #     true_params_list$log_beta_sigma_MU <- hello
          #     #     true_params_list$log_beta_sigma_SD <- hello
          #     #   } else {
          #     #     true_params_list$log_beta_sigma_MU <- hello
          #     #     true_params_list$log_beta_sigma_SD <- hello
          #     #   }
          #     # } else { 
          #     #   true_params_list$log_beta_sigma_MU <- hello
          #     #   true_params_list$log_beta_sigma_SD <- hello
          #     # }
          # }
          ##
          ## Calculate beta_sigma for each test using hierarchical structure
           if (random_C) {
                   true_params_list$beta_sigma_array <-  structure(c(0.267, 0.267, 0.267, 0.267, 
                                                                     0.391, 0.391, 0.391, 0.391), 
                                                                   dim = c(4L, 2L))
                 } else { 
                   true_params_list$beta_sigma_array <- hello
            }
          ##
            if (random_C) {
                    true_params_list$beta_tau_array <-  structure(c(0.073, 0.351, 0.068, 0.172,
                                                                    0.105, 0.143, 0.114, 0.144), 
                                                                  dim = c(4L, 2L))
                } else { 
                    true_params_list$beta_tau_array <-  hello
                }
          ##
               if (random_C) {
                    true_params_list$beta_corr_val <- 0.415
                } else { 
                    true_params_list$beta_corr_val <-    hello
           }
           ##
               if (random_C) {
                     true_params_list$alpha_array <-  structure(c(18.214, 2.265, 2.781, 1.827, 6.485, 5.754, 2.332, 
                                                                  5.623, 23.117, 3.009, 4.087, 0.765, 10.408, 1.135, 3.624, 1.617, 
                                                                  25.15, 6.161, 6.389, 0.884, 15.344, 3.81, 4.758, 2.737, 11.065, 
                                                                  3.89, 8.47, 1.341, 20.463, 3.207, 5.175, 2.822, 10.275, 4.286, 
                                                                  9.787, 1.761, 24.949, 4.786, 7.273, 7.203, 5.376, 2.462, 11.104, 
                                                                  3.325, 27.298, 7.685, 7.228, 4.552, 6.453, 2.89, 12.076, 3.415, 
                                                                  32.748, 12.982, 10.472, 4.72, NaN, NaN, 15.615, 5.713, 34.359, 
                                                                  13.396, 8.906, 6.833, NaN, NaN, 10.61, 4.546, 33.488, 18.349, 
                                                                  10.483, 5.76, NaN, NaN, 11.258, 5.195, 29.546, 25.279, 11.128, 
                                                                  7.879, NaN, NaN, 9.756, 6.06, 22.498, 22.952, 12.59, 4.162, NaN, 
                                                                  NaN, 9.545, 6.868, 22.204, 17.147, 11.92, 13.28, NaN, NaN, 8.243, 
                                                                  6.5, 17.431, 20.853, 16.992, 8.821, NaN, NaN, 9.223, 5.935, 15.465, 
                                                                  21.568, 11.278, 8.901, NaN, NaN, 9.408, 6.678, 13.921, 16.494, 
                                                                  13.326, 17.395, NaN, NaN, 9.614, 4.885, 8.673, 12.661, 18.311, 
                                                                  21.595, NaN, NaN, 7.137, 4.618, 7.678, 7.361, 14.786, 13.271, 
                                                                  NaN, NaN, 5.122, 5.282, 4.731, 4.961, 15.568, 27.593, NaN, NaN, 
                                                                  7.43, 5.128, 2.884, 4.053, 12.57, 15.298, NaN, NaN, 5.444, 4.431, 
                                                                  2.333, 4.264, 10.523, 13.688, NaN, NaN, 4.894, 3.296, 1.356, 
                                                                  1.631, 12.692, 12.839, NaN, NaN, 18.174, 5.272, 3.752, 8.266, 
                                                                  10.508, 13.034, NaN, NaN, NaN, NaN, NaN, NaN, 15.248, 22.748, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 15.58, 13.582, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 13.013, 17.944, NaN, NaN, NaN, NaN, NaN, NaN, 12.375, 
                                                                  7.97, NaN, NaN, NaN, NaN, NaN, NaN, 10.542, 18.443, NaN, NaN, 
                                                                  NaN, NaN, NaN, NaN, 8.037, 17.813, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                                  9.464, 11.579, NaN, NaN, NaN, NaN, NaN, NaN, 8.21, 5.805, NaN, 
                                                                  NaN, NaN, NaN, NaN, NaN, 13.296, 25.709, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 8.933, 19.653, NaN, NaN, NaN, NaN, NaN, NaN, 4.941, 
                                                                  9.716, NaN, NaN, NaN, NaN, NaN, NaN, 9.029, 9.431, NaN, NaN, 
                                                                  NaN, NaN, NaN, NaN, 7.234, 9.886, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                                  13.192, 18.879, NaN, NaN, NaN, NaN, NaN, NaN, 12.255, 5.032, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 15.95, 13.888, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 5.722, 18.15, NaN, NaN, NaN, NaN, NaN, NaN, 5.62, 9.14, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 4.748, 15.735, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 9.298, 21.268, NaN, NaN, NaN, NaN, NaN, NaN, 5.258, 
                                                                  7.051, NaN, NaN, NaN, NaN, NaN, NaN, 10.441, 6.974, NaN, NaN, 
                                                                  NaN, NaN, NaN, NaN, 13.845, 7.219, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                                  6.177, 11.06, NaN, NaN, NaN, NaN, NaN, NaN, 13.505, 11.435, NaN, 
                                                                  NaN, NaN, NaN, NaN, NaN, 3.7, 9.006, NaN, NaN, NaN, NaN, NaN, 
                                                                  NaN, 3.76, 9.297, NaN, NaN, NaN, NaN, NaN, NaN, 4.607, 5.559, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 4.638, 5.336, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 4.533, 4.65, NaN, NaN, NaN, NaN, NaN, NaN, 9.398, 4.597, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 4.363, 6.409, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 4.867, 6.26, NaN, NaN, NaN, NaN, NaN, NaN, 4.926, 7.413, 
                                                                  NaN, NaN, NaN, NaN, NaN, NaN, 6.646, 4.271, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, 4.694, 4.252, NaN, NaN, NaN, NaN, NaN, NaN, 4.582, 
                                                                  4.185, NaN, NaN, NaN, NaN, NaN, NaN, 4.669, 4.233, NaN, NaN, 
                                                                  NaN, NaN, NaN, NaN, 7.546, 4.002, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                                  7.229, 4.034, NaN, NaN, NaN, NaN, NaN, NaN, 4.06, 8.465, NaN, 
                                                                  NaN, NaN, NaN, NaN, NaN, 23.765, 14.544), dim = c(2L, 4L, 64L
                                                                  ))
               } else {
                        ##
                        true_params_list$C_array <- structure(c(-1.279350625, -1.3989578125, -2.088198125, -1.830766875, 
                                                              -2.6383840625, -2.4476190625, -2.6917728125, -2.700195625, -0.53690403125, 
                                                              -0.87724359375, -1.7010740625, -1.7289721875, -2.1997390625, 
                                                              -2.3096178125, -2.36057625, -2.55331375, 0.14995921875, -0.172943375, 
                                                              -1.3778553125, -1.5881553125, -1.85372375, -1.91208875, -2.113551875, 
                                                              -2.363799375, 0.488701875, 0.22914696875, -1.0946284375, -1.381964375, 
                                                              -1.54641, -1.664443125, -1.9521803125, -2.2277634375, 0.913880625, 
                                                              0.7659923125, -0.84729153125, -1.190455625, -1.2822653125, -1.4573559375, 
                                                              -1.7756765625, -1.9919853125, 1.2725884375, 1.2819115625, -0.6035365, 
                                                              -0.96960525, -1.03322625, -1.241093125, -1.6525090625, -1.8843940625, 
                                                              NaN, NaN, -0.40653953125, -0.786740625, -0.75705071875, -0.9706671875, 
                                                              -1.5033575, -1.7933003125, NaN, NaN, -0.1729134375, -0.53436028125, 
                                                              -0.503098125, -0.74578796875, -1.3987778125, -1.6724121875, NaN, 
                                                              NaN, -0.03479440625, -0.36264440625, -0.2592136875, -0.49945003125, 
                                                              -1.2937346875, -1.592301875, NaN, NaN, 0.09812949375, -0.19949378125, 
                                                              -0.044677265625, -0.21000359375, -1.193284375, -1.49689, NaN, 
                                                              NaN, 0.22053890625, -0.024401946875, 0.1324893125, 0.0400360125, 
                                                              -1.0926575, -1.4511075, NaN, NaN, 0.35103753125, 0.163322, 0.32773090625, 
                                                              0.23503415625, -1.00362478125, -1.327385, NaN, NaN, 0.45956578125, 
                                                              0.34110871875, 0.51901975, 0.46777953125, -0.89677496875, -1.253421875, 
                                                              NaN, NaN, 0.583520875, 0.5160165, 0.72745746875, 0.74244659375, 
                                                              -0.83448165625, -1.183728125, NaN, NaN, 0.71920746875, 0.71396990625, 
                                                              0.96566596875, 1.0065576875, -0.7594595625, -1.0621246875, NaN, 
                                                              NaN, 0.8711014375, 0.87908990625, 1.1495478125, 1.27853875, -0.6710876875, 
                                                              -0.933768875, NaN, NaN, 0.99045515625, 1.0483915625, 1.37616875, 
                                                              1.466186875, -0.60155715625, -0.8625738125, NaN, NaN, 1.095938125, 
                                                              1.2662265625, 1.5809475, 1.6287921875, -0.5313750625, -0.7239008125, 
                                                              NaN, NaN, 1.2589646875, 1.5125759375, 1.78248875, 1.8047171875, 
                                                              -0.4751164375, -0.65194665625, NaN, NaN, 1.3979865625, 1.8229225, 
                                                              2.0815734375, 2.1737578125, -0.4277245625, -0.58976925, NaN, 
                                                              NaN, 1.5629325, 2.225593125, 2.580356875, 2.3963771875, -0.375304875, 
                                                              -0.53298515625, NaN, NaN, NaN, NaN, NaN, NaN, -0.32770896875, 
                                                              -0.4773513125, NaN, NaN, NaN, NaN, NaN, NaN, -0.2637201875, -0.38457209375, 
                                                              NaN, NaN, NaN, NaN, NaN, NaN, -0.19615721875, -0.32852490625, 
                                                              NaN, NaN, NaN, NaN, NaN, NaN, -0.1403855, -0.25724896875, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, -0.09052110625, -0.2257794375, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, -0.049528278125, -0.15645459375, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, -0.0176757059375, -0.09103160625, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, 0.0222713875, -0.046481125, NaN, NaN, 
                                                              NaN, NaN, NaN, NaN, 0.0571145625, -0.024236039375, NaN, NaN, 
                                                              NaN, NaN, NaN, NaN, 0.1167284375, 0.07133686875, NaN, NaN, NaN, 
                                                              NaN, NaN, NaN, 0.15891790625, 0.14724853125, NaN, NaN, NaN, NaN, 
                                                              NaN, NaN, 0.18051071875, 0.18661096875, NaN, NaN, NaN, NaN, NaN, 
                                                              NaN, 0.2183054375, 0.22500384375, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                              0.25002765625, 0.265211625, NaN, NaN, NaN, NaN, NaN, NaN, 0.30858903125, 
                                                              0.34360228125, NaN, NaN, NaN, NaN, NaN, NaN, 0.36809125, 0.36376534375, 
                                                              NaN, NaN, NaN, NaN, NaN, NaN, 0.4438813125, 0.4220055625, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, 0.47129015625, 0.502803625, NaN, NaN, 
                                                              NaN, NaN, NaN, NaN, 0.4991404375, 0.543860625, NaN, NaN, NaN, 
                                                              NaN, NaN, NaN, 0.52186671875, 0.6185569375, NaN, NaN, NaN, NaN, 
                                                              NaN, NaN, 0.57020859375, 0.72254278125, NaN, NaN, NaN, NaN, NaN, 
                                                              NaN, 0.59815553125, 0.7613453125, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                              0.6567521875, 0.80196565625, NaN, NaN, NaN, NaN, NaN, NaN, 0.73996571875, 
                                                              0.8442039375, NaN, NaN, NaN, NaN, NaN, NaN, 0.77952334375, 0.91169921875, 
                                                              NaN, NaN, NaN, NaN, NaN, NaN, 0.86567225, 0.98376228125, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, 0.88929209375, 1.048956875, NaN, NaN, 
                                                              NaN, NaN, NaN, NaN, 0.91333184375, 1.118455625, NaN, NaN, NaN, 
                                                              NaN, NaN, NaN, 0.94803990625, 1.1633446875, NaN, NaN, NaN, NaN, 
                                                              NaN, NaN, 0.9834824375, 1.2099759375, NaN, NaN, NaN, NaN, NaN, 
                                                              NaN, 1.0198434375, 1.2519078125, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                              1.0996959375, 1.29196875, NaN, NaN, NaN, NaN, NaN, NaN, 1.13789625, 
                                                              1.3589278125, NaN, NaN, NaN, NaN, NaN, NaN, 1.1853596875, 1.4316540625, 
                                                              NaN, NaN, NaN, NaN, NaN, NaN, 1.236143125, 1.525458125, NaN, 
                                                              NaN, NaN, NaN, NaN, NaN, 1.3126865625, 1.5870815625, NaN, NaN, 
                                                              NaN, NaN, NaN, NaN, 1.3734203125, 1.6527240625, NaN, NaN, NaN, 
                                                              NaN, NaN, NaN, 1.437664375, 1.7259840625, NaN, NaN, NaN, NaN, 
                                                              NaN, NaN, 1.5082409375, 1.8103725, NaN, NaN, NaN, NaN, NaN, NaN, 
                                                              1.6534078125, 1.9122446875, NaN, NaN, NaN, NaN, NaN, NaN, 1.8428928125, 
                                                              2.037545, NaN, NaN, NaN, NaN, NaN, NaN, 1.9887415625, 2.483548125
                        ), dim = c(2L, 4L, 63L))
               }
          ##
          # Test names for indexing
          test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
          n_thr_per_test <- c(1, 
                              6, 21, 21, 63)
          ##
          # Generate covariate matrix for all studies
          X_mat <- simulate_covariates_preserving_missingness_FIXED(original_X = original_cov_data$X, 
                                                           n_studies = n_studies, 
                                                           indicator_index_test_in_study = indicator_index_test_in_study, 
                                                           covariates_same_between_classes = rep(TRUE, covariate_settings$n_covariates),
                                                           add_noise = FALSE,
                                                           noise_sd_prop = 0.0,
                                                           flip_prob = 0.0)
          
          X_mat <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies( 
                                X = X_mat,
                                indicator_index_test_in_study = indicator_index_test_in_study)
          ##
          ##
          # Generate study sizes
          N_per_study_vec <- round(TruncatedNormal::rtnorm(
            n = n_studies, 
            mu = N_per_study_mean, 
            sd = N_per_study_SD, 
            lb = 100
          ), 0)
          ##
          # Generate prevalence parameters
          true_Mean_probit_prev <- qnorm(true_Mean_prev)
          ##
          true_probit_prev_per_study <- rnorm(
            n = n_studies, 
            mean = true_Mean_probit_prev,
            sd = true_SD_probit_prev
          )
          true_prev_per_study <- pnorm(true_probit_prev_per_study)
          ##
          # Initialize storage
          y_list <- list()
          y_df_list <- list()
          ##
          # Initialize storage arrays
          max_threshold_across_all_tests <- max(n_thr_per_test)
          n_TP_at_each_threshold <- array(-1, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          n_FP_at_each_threshold <- array(-1, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Se_all_tests_all_thresholds <- array(-1, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Sp_all_tests_all_thresholds <- array(-1, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          Fp_all_tests_all_thresholds <- array(-1, dim = c(n_studies, n_tests, max_threshold_across_all_tests + 1))
          ##
          Mean_of_thr_for_all_tests_array_per_study_nd <- array(NA, dim = c(n_studies, n_tests, max_threshold_across_all_tests))
          Mean_of_thr_for_all_tests_array_per_study_d  <- array(NA, dim = c(n_studies, n_tests, max_threshold_across_all_tests))
          ##
          s <- 1
          ##
          # str(X_mat)
          ##
          n_pos_OVERALL <- 0
          n_neg_OVERALL <- 0
          n_total_OVERALL <- 0
          ##
          # Generate data for each study
          for (s in 1:n_studies) {
            
                    N <- N_per_study_vec[s]
                    ##
                    true_prev <- true_prev_per_study[s]
                    d_ind <- sort(rbinom(n = N, size = 1, prob = true_prev))
                    ##
                    n_pos <- sum(d_ind)
                    n_pos_OVERALL <- n_pos_OVERALL + n_pos
                    ##
                    n_neg <- N - sum(d_ind)
                    n_neg_OVERALL <- n_neg_OVERALL + n_neg
                    ##
                    n_total <- n_pos + n_neg
                    n_total_OVERALL <- n_total_OVERALL + n_total
                    ##
                    # Initialize threshold arrays for current study
                    thr_for_all_tests_for_current_study_nd    <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                    thr_for_all_tests_for_current_study_d     <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                    thr_for_all_tests_for_current_study_per_n <- array(NA, dim = c(N, n_tests, max_threshold_across_all_tests))
                    ##
                    # Initialize arrays for this study
                    y <- array(NA, dim = c(N, n_tests))
                    
                    # Get which tests are available in this study
                    available_tests <- which(indicator_test_in_study[s, ]) ; available_tests
                    
                    # If reference test is available, generate it
                    if (1 %in% available_tests) {
                      if (assume_perfect_GS) {
                        y[, 1] <- d_ind
                      } else {
                        # Imperfect reference test
                        ref_Se <- 0.95
                        ref_Sp <- 0.90
                        y[d_ind == 1, 1] <- rbinom(n_pos, 1, ref_Se)
                        y[d_ind == 0, 1] <- rbinom(n_neg, 1, 1 - ref_Sp)
                      }
                    }
                    
                    # Only generate random effects if there are ordinal tests in this study
                    ordinal_tests_in_study <- intersect(available_tests, 2:n_tests)
                    
                    if (length(ordinal_tests_in_study) > 0) {
                      
                          # Generate shared study-level random effects
                          beta_corr <- true_params_list$beta_corr_val
                          beta_L_Omega <- matrix(c(1, 0, beta_corr, sqrt(1 - beta_corr^2)), 2, 2)
                          ## Create test-specific Cholesky factor
                          ##
                          beta_sigma <- true_params_list$beta_sigma_array[1, ]
                          beta_L_Sigma <- diag(beta_sigma) %*% beta_L_Omega
                          ## Generate study-level random effects for this test
                          beta_eta_z_s <- rnorm(n = 2, mean = 0.0, sd = 1.0)
                          beta_eta_s <- as.vector(beta_L_Sigma %*% beta_eta_z_s)
                          
                          # Generate data for each available ordinal test
                          for (t in ordinal_tests_in_study) {
                                  
                                  ## Convert to ordinal categories
                                  n_thr_t <- n_thr_per_test[t]; ##  - 1
                                  n_cat_t <- n_thr_t + 1
                                  ##
                                  test_name <- test_names[t - 1]
                                  ##
                                  # Get test-specific parameters
                                  beta_mu_nd <- true_params_list$beta_mu_array[t - 1, 1, ]
                                  beta_mu_d  <- true_params_list$beta_mu_array[t - 1, 2, ]
                                  ##
                                  beta_tau   <- true_params_list$beta_tau_array[t - 1, ]
                                  ##
                                  alpha_nd   <- true_params_list$alpha_array[1, t - 1, 1:n_cat_t]
                                  alpha_d    <- true_params_list$alpha_array[2, t - 1, 1:n_cat_t]
                                  # length(alpha_nd)
                                  ##
                                  if (random_C) {
                                      thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                        alpha_vec =  alpha_nd,
                                        use_probit_link =  TRUE)
                                      thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]] <- generate_ordered_thresholds_induced_dirichlet(
                                        alpha_vec =  alpha_d,
                                        use_probit_link =  TRUE)
                                  } else { 
                                    thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]] <- true_params_list$C_array[1, t - 1, 1:n_thr_t]
                                    thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]  <- true_params_list$C_array[2, t - 1, 1:n_thr_t]
                                  }
                                  ##
                                  {
                                    ##
                                    Mean_of_thr_for_all_tests_array_per_study_nd[s, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]]
                                    Mean_of_thr_for_all_tests_array_per_study_d[s, t, 1:n_thr_per_test[t]]  <- thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]
                                    ##
                                    # Assign thresholds to diseased and non-diseased subjects
                                    diseased_indices    <- which(d_ind == 1)
                                    nondiseased_indices <- which(d_ind == 0)
                                    ##
                                    for (n in diseased_indices) {
                                      thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_d[t, 1:n_thr_per_test[t]]
                                    }
                                    for (n in nondiseased_indices) {
                                      thr_for_all_tests_for_current_study_per_n[n, t, 1:n_thr_per_test[t]] <- thr_for_all_tests_for_current_study_nd[t, 1:n_thr_per_test[t]]
                                    }
                                  }
                                  
                                  ## Generate test-specific deviations
                                  beta_delta <- c()
                                  beta_delta[1] <- rnorm(1, mean = 0, sd = beta_tau[1])
                                  beta_delta[2] <- rnorm(1, mean = 0, sd = beta_tau[2])
                                  
                                  ## Total random effects
                                  beta_random_nd <- beta_eta_s[1] + beta_delta[1]
                                  beta_random_d  <- beta_eta_s[2] + beta_delta[2]
                                  
                                  ## Calculate study-specific location parameters with covariates
                                  Xbeta_nd <- sum(X_mat[[1]][[t - 1]][s, ] * beta_mu_nd)
                                  Xbeta_d  <- sum(X_mat[[2]][[t - 1]][s, ] * beta_mu_d)
                                  
                                  ## Total location parameters including random effects
                                  location_nd_s <- Xbeta_nd + beta_random_nd
                                  location_d_s  <- Xbeta_d  + beta_random_d
                                  ##
                                  ## Generate latent continuous values
                                  latent_values <- numeric(N)
                                  latent_values[d_ind == 0] <- rnorm(n_neg, mean = location_nd_s, sd = 1)
                                  latent_values[d_ind == 1] <- rnorm(n_pos, mean = location_d_s,  sd = 1)
                                  ##
                                
                                  {
                                      # First category (below first threshold)
                                      threshold_lower <- rep(-9999, N)
                                      threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, 1]
                                      first_category <- (latent_values > threshold_lower) & (latent_values <= threshold_upper)
                                      y[first_category, t] <- 0
                                      
                                      # Middle categories
                                      for (k in 2:(n_thr_t)) {
                                        threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, k - 1]
                                        threshold_upper <- thr_for_all_tests_for_current_study_per_n[1:N, t, k]
                                        category_k <- (latent_values > threshold_lower) & (latent_values <= threshold_upper)
                                        y[category_k, t] <- k - 1
                                      }
                                      
                                      # Last category (above last threshold)
                                      threshold_lower <- thr_for_all_tests_for_current_study_per_n[1:N, t, n_thr_t]
                                      threshold_upper <- rep(9999, N)
                                      last_category <- (latent_values > threshold_lower) & (latent_values <= threshold_upper)
                                      y[last_category, t] <- n_cat_t - 1
                                  }
                                  ##
                                  # Calculate TP and FP at each threshold
                                  # Initialize counts for first threshold (everyone positive)
                                  n_TP_at_each_threshold[s, t, 1] <- n_pos
                                  n_FP_at_each_threshold[s, t, 1] <- n_neg
                                  ##
                                  # Calculate counts for each threshold
                                  for (k in 2:(n_cat_t)) {
                                    df_pos_y <- y[d_ind == 1, t]
                                    df_neg_y <- y[d_ind == 0, t]
                                    ##
                                    n_TP_at_each_threshold[s, t, k] <- sum(df_pos_y >= (k - 1))
                                    n_FP_at_each_threshold[s, t, k] <- sum(df_neg_y >= (k - 1))
                                  }
                                  ##
                                  # Calculate sensitivity and specificity
                                  if (n_pos > 0) {
                                    Se_all_tests_all_thresholds[s, t, 1:n_cat_t] <- n_TP_at_each_threshold[s, t, 1:n_cat_t] / n_pos
                                  }
                                  if (n_neg > 0) {
                                    Fp_all_tests_all_thresholds[s, t, 1:n_cat_t] <- n_FP_at_each_threshold[s, t, 1:n_cat_t] / n_neg
                                    Sp_all_tests_all_thresholds[s, t, 1:n_cat_t] <- 1 - Fp_all_tests_all_thresholds[s, t, 1:n_cat_t]
                                  }
                            
                          }
                          
                    }
                    
                    # Store results
                    y_list[[s]] <- y
                    study_df <- data.frame(
                      Study = s,
                      y,
                      d_ind = d_ind
                    )
                    colnames(study_df)[2:(n_tests + 1)] <- c("Ref_test", test_names)
                    y_df_list[[s]] <- study_df
            
          }
          
          # Combine results
          y_tibble <- dplyr::bind_rows(y_df_list)
          
          # Create summary of which studies have which tests
          study_test_summary <- data.frame(
            study = 1:n_studies,
            indicator_test_in_study
          )
          colnames(study_test_summary)[-1] <- c("Ref_test", test_names)
          
          # Return results
          return(list(
            y_list = y_list,
            y_df_list = y_df_list,
            y_tibble = y_tibble,
            X_mat = X_mat,
            N_per_study_vec = N_per_study_vec,
            true_prev_per_study = true_prev_per_study,
            n_thr_per_test = n_thr_per_test,
            n_thr_per_test_excl_ref = n_thr_per_test[-1],
            n_TP_at_each_threshold = n_TP_at_each_threshold,
            n_FP_at_each_threshold = n_FP_at_each_threshold,
            Se_all_tests_all_thresholds = Se_all_tests_all_thresholds,
            Sp_all_tests_all_thresholds = Sp_all_tests_all_thresholds,
            Fp_all_tests_all_thresholds = Fp_all_tests_all_thresholds,
            test_names = test_names,
            n_studies = n_studies,
            true_params_list = true_params_list,
            indicator_index_test_in_study = indicator_index_test_in_study,
            study_test_summary = study_test_summary,
            n_studies_per_test = colSums(indicator_index_test_in_study)  # Actual studies per test
          ))
  
        

  
}

 