# ============================================================================================================
# Helper Functions for MetaOrdDTA Simulation
# ============================================================================================================






#       
# # Create file name string
# create_file_name_string <- function(scenario) {
#   
#         paste0("test_", scenario$index_test_MA,
#                "_N_", scenario$N_per_study_mean,
#                "_n_studies_", scenario$n_studies,
#                "_param_", scenario$model_parameterisation,
#                "_soft_", scenario$softplus,
#                "_rand_thr_", scenario$random_thresholds,
#                "_alpha_pi_", scenario$alpha_pi)
#   
# }




create_file_name_string <- function( scenario_row, 
                                     reduction_factor
) {

        paste0("test_", scenario_row$index_test_MA,
               "red_fact_", reduction_factor, ## new
               "DGM_", scenario_row$DGM,
               "N_", scenario_row$N_per_study_mean,  # No underscore before N
               "n_studies_", scenario_row$n_studies,  # No underscore before n_studies
               "param_", scenario_row$model_parameterisation,  # No underscore before param
               "soft_", scenario_row$softplus,  # No underscore before soft
               "rand_thr_", scenario_row$random_thresholds,  # No underscore before rand_thr
               "alpha_pi_", scenario_row$alpha_pi)  # No underscore before alpha_pi

}




rename_files_add_red_fact <- function( directory,
                                       reduction_factor,
                                       dry_run = TRUE
) {
  
        # Get all RDS files
        files <- list.files(directory, 
                            pattern = "^seed_.*\\.RDS$", 
                            full.names = FALSE)
        
        cat(sprintf("Found %d files to rename\n", length(files)))
        
        for (file in files) {
          
          ## Extract the seed number and the rest
          if (grepl("^seed_\\d+_test_", file)) {
            
                # Split: "seed_4121_test_BAI..." -> "seed_4121_" + "test_BAI..."
                parts <- strsplit(file, "test_")[[1]]
                seed_part <- parts[1]  # "seed_4121_"
                rest_part <- parts[2]  # "BAIDGM_Xu_random_C..."
                
                # Create new filename with red_fact inserted
                new_file <- paste0(seed_part, 
                                   "test_", 
                                   strsplit(rest_part, "DGM_")[[1]][1],  # "BAI"
                                   "red_fact_", reduction_factor,
                                   "DGM_", 
                                   strsplit(rest_part, "DGM_")[[1]][2])  # Rest after "DGM_"
                
                if (dry_run) {
                  cat(sprintf("Would rename:\n  %s\n  -> %s\n\n", file, new_file))
                } else {
                  file.rename(
                    file.path(directory, file),
                    file.path(directory, new_file)
                  )
                  cat(sprintf("Renamed: %s\n", file))
                }
            
            }
          
        }
        
        if (dry_run) {
          cat("\n*** DRY RUN - No files actually renamed ***\n")
          cat("Set dry_run = FALSE to actually rename files\n")
        }

}

 






get_true_values <- function(DGM, 
                            test
) {
  
        # Initialize lists
        true_DGM_Se <- list()
        true_DGM_Sp <- list()
        true_DGM_Fp <- list()
        
        if (test == "GAD_2" || test == "GAD-2") {
          
              if (DGM == "Xu_random_C") { 
                
                    true_DGM_Se <- c(95.318, 85.87, 63.02, 47.236, 27.83, 15.6)/100
                    true_DGM_Fp <- c(59.589, 32.905, 13.731, 7.66, 3.262, 1.429)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Xu_fixed_C") {
                
                    true_DGM_Se <- c(94.064, 85.133, 63.202, 47.307, 27.227, 13.043)/100
                    true_DGM_Fp <- c(61.364, 32.558, 12.755, 6.989, 2.867, 1.19)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Jones") {
                
                    true_DGM_Se <- c(95.456, 83.1, 63.551, 42.164, 24.192, 12.151)/100
                    true_DGM_Fp <- c(60.355, 31.206, 14.954, 6.836, 3.045, 1.334)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "R&G_random_C") {
                
                    true_DGM_Se <- c(96.176, 85.168, 62.919, 47.441, 30.192, 17.019)/100
                    true_DGM_Fp <- c(60.493, 34.353, 14.575, 7.832, 3.341, 1.277)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              }
          
        } else if (test == "HADS") {
          
              if (DGM == "Xu_random_C") {
                    
                    true_DGM_Se <- c(98.195, 97.838, 96.432, 95.247, 93.466, 90.61, 85.748, 80.685, 
                                     73.559, 63.65, 54.293, 46.918, 38.05, 28.523, 21.017, 15.068, 
                                     11.544, 9.083, 7.075, 4.881, 4.042)/100
                    true_DGM_Fp <- c(92.824, 84.73, 75.338, 65.239, 54.946, 45.399, 35.677, 27.033, 
                                     19.814, 14.339, 10.693, 7.51, 5.303, 3.568, 2.209, 1.449, 0.852, 
                                     0.527, 0.347, 0.211, 0.139)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Xu_fixed_C") {
                    
                    true_DGM_Se <- c(99.491, 99.247, 97.938, 96.351, 94.374, 91.498, 86.541, 81.125, 
                                     73.816, 63.739, 54.087, 46.279, 37.224, 27.342, 19.288, 12.74, 
                                     9.258, 6.843, 4.81, 2.109, 1.208)/100
                    true_DGM_Fp <- c(93.422, 85.736, 76.484, 66.134, 56.019, 46.139, 35.456, 26.526, 
                                     19.198, 13.898, 10.345, 7.254, 4.968, 3.161, 1.807, 1.131, 0.613, 
                                     0.338, 0.18, 0.067, 0.011)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Jones") {
                
                    true_DGM_Se <- c(99.888, 99.556, 98.745, 97.127, 94.344, 90.091, 84.233, 76.827, 
                                     68.134, 58.67, 48.897, 39.432, 30.772, 23.266, 17.008, 12.027, 
                                     8.247, 5.469, 3.518, 2.196, 1.332)/100
                    true_DGM_Fp <- c(92.451, 84.824, 75.383, 64.84, 53.97, 43.523, 34.02, 25.813, 
                                     19.018, 13.642, 9.52, 6.477, 4.302, 2.791, 1.773, 1.102, 0.67, 
                                     0.4, 0.234, 0.134, 0.076)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "R&G_random_C") {
                
                    true_DGM_Se <- c(99.844, 99.415, 98.524, 97.028, 94.781, 91.778, 87.34, 81.7, 
                                     74.824, 67.179, 60.115, 52.512, 44.642, 36.561, 28.316, 22.163, 
                                     16.586, 12.465, 9.371, 6.466, 5.036)/100
                    true_DGM_Fp <- c(92.811, 84.788, 75.295, 65.295, 55.257, 45.983, 36.382, 27.891, 
                                     20.59, 14.814, 10.881, 7.702, 5.227, 3.351, 1.969, 1.209, 0.698, 
                                     0.416, 0.251, 0.134, 0.088)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              }
          
        } else if (test == "BAI") {
          
              if (DGM == "Xu_random_C") {
                
                    true_DGM_Se <- c(98.165, 97.618, 96.739, 95.884, 93.859, 92.615, 91.362, 89.56, 
                                     88.148, 86.218, 85.219, 82.263, 80.332, 78.382, 74.833, 70.772, 
                                     68.394, 63.505, 60.846, 58.597, 56.453, 54.356, 50.815, 48.705, 
                                     46.015, 44.798, 42.21, 39.824, 38.186, 37.321, 33.917, 31.312, 
                                     30.039, 28.747, 27.441, 25.02, 24.331, 22.686, 20.412, 19.329, 
                                     17.543, 15.165, 14.314, 13.438, 12.612, 11.389, 10.173, 9.134, 
                                     8.138, 7.528, 6.925, 6.409, 5.953, 5.274, 4.622, 3.902, 3.488, 
                                     3.095, 2.718, 2.283, 1.905, 1.533, 0.844)/100
                    ##
                    true_DGM_Fp <- c(88.562, 80.336, 72.952, 67.435, 61.549, 56.597, 50.648, 46.35, 
                                     42.087, 38.317, 34.661, 31.583, 27.836, 25.626, 23.301, 20.496, 
                                     18.459, 16.608, 15.197, 14.12, 12.87, 11.969, 10.742, 9.564, 
                                     8.685, 7.905, 7.29, 6.851, 6.338, 5.929, 5.313, 4.914, 4.681, 
                                     4.311, 4.028, 3.567, 3.16, 2.679, 2.508, 2.347, 2.215, 1.981, 
                                     1.853, 1.621, 1.344, 1.222, 0.99, 0.918, 0.856, 0.782, 0.704, 
                                     0.635, 0.519, 0.462, 0.404, 0.348, 0.285, 0.238, 0.198, 0.162, 
                                     0.112, 0.07, 0.05)/100
                    ##
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Xu_fixed_C") {
                    
                    true_DGM_Se <- c(99.011, 98.544, 97.668, 96.81, 94.69, 93.412, 92.18, 90.283, 
                                     88.818, 86.889, 85.876, 82.948, 81.039, 79.021, 75.354, 71.185, 
                                     68.692, 63.746, 60.947, 58.593, 56.38, 54.23, 50.576, 48.365, 
                                     45.592, 44.417, 41.764, 39.272, 37.628, 36.799, 33.303, 30.636, 
                                     29.333, 28.015, 26.716, 24.258, 23.65, 21.887, 19.63, 18.52, 
                                     16.633, 14.203, 13.361, 12.487, 11.664, 10.422, 9.183, 8.126, 
                                     7.118, 6.528, 5.96, 5.484, 5.072, 4.414, 3.776, 3.063, 2.687, 
                                     2.293, 1.937, 1.566, 1.208, 0.877, 0.239)/100
                    ##
                    true_DGM_Fp <- c(88.642, 80.961, 73.593, 68.007, 61.472, 56.666, 50.811, 46.624, 
                                     42.504, 38.592, 34.803, 31.571, 27.857, 25.808, 23.47, 20.864, 
                                     18.944, 17.064, 15.719, 14.574, 13.385, 12.421, 11.146, 9.891, 
                                     8.959, 8.179, 7.577, 7.11, 6.596, 6.142, 5.475, 5.028, 4.807, 
                                     4.438, 4.144, 3.652, 3.214, 2.7, 2.529, 2.369, 2.246, 2, 1.875, 
                                     1.619, 1.313, 1.189, 0.949, 0.891, 0.83, 0.757, 0.681, 0.621, 
                                     0.492, 0.439, 0.381, 0.327, 0.258, 0.216, 0.176, 0.139, 0.086, 
                                     0.045, 0.026)/100
                    ##
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "Jones") {
                
                    true_DGM_Se <- c(99.628, 99.188, 98.586, 97.773, 96.759, 95.532, 94.08, 92.397, 
                                     90.494, 88.395, 86.122, 83.66, 81.031, 78.229, 75.326, 72.289, 
                                     69.176, 66.079, 62.94, 59.78, 56.594, 53.464, 50.387, 47.384, 
                                     44.435, 41.552, 38.779, 36.078, 33.453, 30.998, 28.618, 26.41, 
                                     24.338, 22.332, 20.476, 18.696, 17.061, 15.543, 14.101, 12.777, 
                                     11.57, 10.456, 9.421, 8.488, 7.622, 6.85, 6.135, 5.487, 4.897, 
                                     4.37, 3.886, 3.452, 3.072, 2.722, 2.407, 2.128, 1.877, 1.655, 
                                     1.454, 1.276, 1.123, 0.986, 0.863)/100
                    true_DGM_Fp <- c(88.285, 80.51, 73.28, 66.546, 60.267, 54.535, 49.281, 44.442, 
                                     40.058, 36.056, 32.468, 29.16, 26.149, 23.458, 21.037, 18.841, 
                                     16.901, 15.13, 13.537, 12.115, 10.845, 9.709, 8.687, 7.762, 6.926, 
                                     6.193, 5.539, 4.938, 4.405, 3.94, 3.523, 3.145, 2.811, 2.508, 
                                     2.243, 2.005, 1.792, 1.6, 1.429, 1.275, 1.139, 1.018, 0.908, 
                                     0.811, 0.722, 0.643, 0.574, 0.512, 0.456, 0.407, 0.362, 0.323, 
                                     0.289, 0.258, 0.231, 0.206, 0.183, 0.164, 0.146, 0.13, 0.117, 
                                     0.104, 0.093)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              } else if (DGM == "R&G_random_C") {
                    
                    true_DGM_Se <- c(99.749, 99.3, 98.626, 97.931, 96.839, 95.8, 94.261, 92.873, 
                                     91.173, 89.522, 87.823, 85.811, 83.324, 81.626, 79.435, 76.254, 
                                     74.014, 70.974, 69.209, 67.173, 65.147, 63.493, 60.649, 57.92, 
                                     55.782, 53.744, 51.646, 50.163, 48.354, 47.687, 44.314, 42.413, 
                                     41.389, 39.878, 38.482, 35.722, 35.191, 32.31, 31.091, 30.206, 
                                     28.632, 25.861, 25.908, 23.858, 21.558, 20.51, 18.364, 16.689, 
                                     16.35, -100, 14.202, 14.425, 12.609, 12.413, 10.388, 9.906, 9.384, 
                                     -100, -100, 7.155, -100, 4.563, 3.094)/100
                    true_DGM_Fp <- c(89.437, 81.813, 74.441, 68.949, 62.382, 57.464, 51.616, 47.346, 
                                     42.8, 39.076, 35.753, 32.372, 28.716, 26.552, 23.992, 20.912, 
                                     18.958, 16.64, 15.362, 14.063, 12.869, 11.973, 10.545, 9.325, 
                                     8.471, 7.695, 6.949, 6.505, 5.904, 5.735, 4.836, 4.37, 4.139, 
                                     3.808, 3.52, 2.973, 2.897, 2.404, 2.225, 2.097, 1.876, 1.528, 
                                     1.528, 1.307, 1.075, 0.981, 0.789, 0.664, 0.637, -100, 0.493, 
                                     0.508, 0.394, 0.393, 0.283, 0.26, 0.24, -100, -100, 0.15, -100, 
                                     0.071, 0.038)/100
                    true_DGM_Sp <- 1 - true_DGM_Fp
                
              }
          
          # } else if (test == "GAD_7" || test == "GAD-7") {
          #   
          #   # true_DGM_Se <- c(97.301, 95.945, 94.527, 92.268, 89.168, 84.232, 79.27, 71.388, 
          #   #                  65.608, 59.321, 52.195, 44.815, 37.824, 31.994, 25.598, 21.024, 
          #   #                  17.119, 12.366, 8.378, 4.91, 2.771)/100
          #   # 
          #   # true_DGM_Fp <- c(79.396, 65.052, 52.811, 42.426, 33.381, 26.387, 20.438, 14.645, 
          #   #                  11.63, 9.01, 7.124, 5.612, 4.498, 3.445, 2.525, 1.805, 1.29, 
          #   #                  1.003, 0.659, 0.422, 0.287)/100
          #   # 
          #   # true_DGM_Sp <- c(20.604, 34.948, 47.189, 57.574, 66.619, 73.613, 79.562, 85.355, 
          #   #                  88.37, 90.99, 92.876, 94.388, 95.502, 96.555, 97.475, 98.195, 
          #   #                  98.71, 98.997, 99.341, 99.578, 99.713)/100
          
        } else {
          stop(paste("Unknown test type:", test))
        }
        
        # Return as percentages (multiply by 100) to match your original code
        return(list(
          Se = true_DGM_Se * 100,
          Sp = true_DGM_Sp * 100,
          Fp = true_DGM_Fp * 100
        ))
        
}






# Create all scenarios
create_MetaOrd_scenarios <- function() {
  
        require(dplyr)
        
        # Define all parameter combinations
        scenario_df <- expand.grid(
          index_test_MA = c("GAD_2", "HADS", "BAI"),
          DGM = c("Xu_fixed_C", "Jones", "Xu_random_C", "R&G_random_C"),
          n_studies = c(10, 50),
          model_parameterisation = c("Jones", "Xu", "Bivariate", "R&G"),
          alpha_pi = c("one", "flat", "custom"),
          # box_cox = c(FALSE, TRUE),
          softplus = c(FALSE),
          random_thresholds = c(FALSE, TRUE),
          stringsAsFactors = FALSE
        ) %>%
          tibble() %>%
          # Filter out invalid combinations (from your original code logic)
          filter(
            # !(model_parameterisation == "Jones" & box_cox == FALSE),
            !((DGM == "Jones") & (alpha_pi == "custom")),
            !(DGM == "Xu_fixed_C" & (random_thresholds == TRUE)),
            !((DGM == "Jones") & (random_thresholds == TRUE)),
            !(model_parameterisation == "Bivariate" & random_thresholds == TRUE),
            !(model_parameterisation == "R&G" & random_thresholds == FALSE),
            !(model_parameterisation %in% c("Jones", "Bivariate") & random_thresholds == TRUE)
          ) %>%
          mutate(
            scenario_id = row_number(),
            Dirichlet_random_effects_type = case_when(
              random_thresholds & model_parameterisation == "Xu" ~ "kappa",
              random_thresholds & model_parameterisation == "R&G" ~ "kappa",
              TRUE ~ "fixed"
            ),
            # Add other fixed parameters
            N_per_study_mean = 500,
            missing_thr = TRUE
          )
        
        return(scenario_df)
  
}





# Create run matrix for selected scenarios (like your LC-MVP)
create_run_matrix_for_MetaOrd <- function(scenario_df, 
                                          N_sims_per_scenario
) {
  
        run_matrix <- data.frame()
        
        for (i in 1:nrow(scenario_df)) {
          
              scenario_row <- scenario_df[i, ]
              
              for (sim_id in 1:N_sims_per_scenario) {
                
                  run_matrix <- rbind(run_matrix, data.frame(
                    scenario_id = scenario_row$scenario_id,
                    sim_id = sim_id,
                    index_test_MA = scenario_row$index_test_MA,
                    n_studies = scenario_row$n_studies,
                    model_parameterisation = scenario_row$model_parameterisation,
                    alpha_pi = scenario_row$alpha_pi,
                    box_cox = scenario_row$box_cox,
                    softplus = scenario_row$softplus,
                    random_thresholds = scenario_row$random_thresholds,
                    Dirichlet_random_effects_type = scenario_row$Dirichlet_random_effects_type,
                    stringsAsFactors = FALSE
                  ))
                  
              }
          
        }
        
        return(tibble(run_matrix))
  
}







check_MCSE_status <- function(output_dir, 
                              file_name_string,
                              max_seed_to_check,
                              min_runs_needed = 30,
                              target_MCSE_RMSE,
                              test,
                              DGM,
                              n_studies,
                              min_ESS_thr = NULL,
                              vec_index_inner = NULL
) {
  
          # Set threshold ranges
          if (is.null(vec_index_inner)) {
              if (test == "GAD_2") {
                 vec_index_inner <- 1:6
              } else if (test == "HADS") {
                vec_index_inner <- 3:17
              } else if (test == "BAI") {
                vec_index_inner <- 5:43
              }
          }
  
          if (is.null(min_ESS_thr)) { 
            if (n_studies == 10) { 
               min_ESS_thr <- 4000
              # min_ESS_thr <- 10000
            } else { 
               min_ESS_thr <- 500
              #min_ESS_thr <- 2000
            }
          }
  
          true_vals <- get_true_values( test = test,
                                        DGM = DGM)
          true_DGM_Se <- true_vals$Se
          true_DGM_Sp <- true_vals$Sp
          true_DGM_Fp <- true_vals$Fp
  
          # Load completed results (from your existing code)
          complete_results_list <- list()
          
          for (seed in 1:max_seed_to_check) {
                
                seed_file <- file.path(output_dir,
                                       paste0("seed_", seed, "_", 
                                              file_name_string, "_results_list.RDS"))
                if (file.exists(seed_file)) {
                  try({
                    complete_results_list[[seed]] <- readRDS(seed_file)
                  }, silent = TRUE)
                }
            
          }
          
          # Count valid results
          complete_results_list <- complete_results_list[!sapply(complete_results_list, is.null)]
          n_sims <- length(complete_results_list)
          
          if (n_sims < min_runs_needed) {
            return(list(
              mcse_met = FALSE,
              n_sims = n_sims,
              mcse_value = NA,
              min_sims_needed = NA
            ))
          }
          
          # # Get true values (from your existing code)
          # true_DGM_Se <- true_vals_list_for_each_DGM[[1]]$Se_true_observed_list[[1]]
          # true_DGM_Sp <- true_vals_list_for_each_DGM[[1]]$Sp_true_observed_list[[1]]
          ##
          metrics <-  compute_simulation_metrics( test = test,
                                                  output_dir = output_dir,
                                                  complete_results_list = complete_results_list,
                                                  true_DGM_Se = true_DGM_Se,
                                                  true_DGM_Sp = true_DGM_Sp,
                                                  min_ESS_thr = min_ESS_thr,
                                                  vec_index_inner = vec_index_inner)
          
          metrics$overall
          
          # Calculate minimum sims needed (from your existing code)
          min_sims_needed <-  calc_min_sim_size_RMSE( metrics = metrics,
                                                      target_MCSE_RMSE = target_MCSE_RMSE,
                                                      vec_index_inner = vec_index_inner)
          
          # Get worst MCSE
          tables <- create_metrics_summary_table( metrics = metrics,
                                                  vec_index_inner = vec_index_inner)
          
          # This is from your existing code pattern
          tables$by_threshold$MCSE_RMSE_Se <- as.numeric(tables$by_threshold$MCSE_RMSE_Se)
          tables$by_threshold$MCSE_RMSE_Sp <- as.numeric(tables$by_threshold$MCSE_RMSE_Sp)
          ##
          mean_MCSE_RMSE_Se <- mean(tables$by_threshold$MCSE_RMSE_Se, na.rm = TRUE)
          mean_MCSE_RMSE_Sp <- mean(tables$by_threshold$MCSE_RMSE_Sp, na.rm = TRUE)
          ##
          max_MCSE <- max(mean_MCSE_RMSE_Se, mean_MCSE_RMSE_Sp)
          ##
          # MCSE_RMSE_Se_Sp <- c(mean(tables$by_threshold$MCSE_RMSE_Se),
          #                      mean(tables$by_threshold$MCSE_RMSE_Sp))
          # max_MCSE <- max(MCSE_RMSE_Se_Sp, na.rm = TRUE)
          ##
          mcse_met <- (!is.na(max_MCSE) && max_MCSE <= target_MCSE_RMSE)
          
          # tables$overall
          # tables$overall$Value
          
          return(list(
            mcse_met = mcse_met,
            n_sims = n_sims,
            mcse_value = max_MCSE,
            min_sims_needed = min_sims_needed,
            metrics = metrics,
            tables = tables
          ))
  
}

# metrics$overall


      





# ============================================================================================================
# Process Single Scenario Function - MetaOrdDTA Version
# ============================================================================================================
process_single_scenario_MetaOrd <- function(output_dir,
                                            reduction_factor,
                                            scenario_row,
                                            target_MCSE_RMSE,
                                            max_seed_to_check,
                                            min_runs_needed,
                                            min_ESS_thr,
                                            vec_index_inner
) {
        
        # Extract scenario parameters
        DGM <- scenario_row$DGM
        test <- scenario_row$index_test_MA
        n_studies <- scenario_row$n_studies
        model_parameterisation <- scenario_row$model_parameterisation
        ##
        # if (is.null(vec_index_inner)) {   # Set threshold ranges
        #   if (test == "GAD_2") {
        #     vec_index_inner <- 1:6
        #   } else if (test == "HADS") {
        #     vec_index_inner <- 1:21
        #     # vec_index_inner <- 3:17
        #     # vec_index_inner <- 5:15
        #     # vec_index_inner <- 7:13
        #   } else if (test == "BAI") {
        #     vec_index_inner <- 3:43
        #   }
        # }
        ##
        # min_k <- head(vec_index_inner, 1)
        # max_k <- tail(vec_index_inner, 1)
        ##
        # Get true values
        true_vals <- get_true_values(test = test,
                                     DGM = DGM)
        true_DGM_Se <- true_vals$Se
        true_DGM_Sp <- true_vals$Sp
        
        # Create file name string
        file_name_string <- create_file_name_string( scenario_row = scenario_row,
                                                     reduction_factor = reduction_factor)
        
        # Load results
        complete_results_list <- list()
        for (seed in 1:max_seed_to_check) {
          file_path <- file.path(output_dir, 
                                 paste0("seed_", seed, "_", file_name_string, "_results_list.RDS"))
          if (file.exists(file_path)) {
            try({
              complete_results_list[[seed]] <- readRDS(file_path)
            }, silent = TRUE)
          }
        }
        
        # Remove NULL entries
        complete_results_list <- complete_results_list[!sapply(complete_results_list, is.null)]
        
        if (length(complete_results_list) < min_runs_needed) {
              
              warning(sprintf("Insufficient results for scenario %d (found %d, need %d)", 
                              scenario_row$scenario_id, length(complete_results_list), min_runs_needed))
              return(NULL)
          
        }
        
        # Compute metrics
        metrics <-  compute_simulation_metrics(
                    test = test,
                    complete_results_list = complete_results_list,
                    output_dir = output_dir,
                    true_DGM_Se = true_DGM_Se,
                    true_DGM_Sp = true_DGM_Sp,
                    min_ESS_thr = min_ESS_thr,
                    vec_index_inner = vec_index_inner)
        
        mean_min_ESS <- metrics$mean_min_ESS
     
        # Create summary table
        metrics_tables <- create_metrics_summary_table(
                          metrics = metrics,
                          vec_index_inner = vec_index_inner)
        
        # # Calculate MCSE status (average across thresholds)
        # MCSE_RMSE_Se_Sp <- c(
        #   mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE),
        #   mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
        # )
        # max_MCSE <- max(MCSE_RMSE_Se_Sp, na.rm = TRUE)
        # Take MAX MCSE across thresholds, not mean!
        max_MCSE_Se <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
        max_MCSE_Sp <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
        max_MCSE <- max(max_MCSE_Se, max_MCSE_Sp)
        mcse_met <- (max_MCSE <= target_MCSE_RMSE)
        
        # Calculate minimum simulations needed
        min_sims_needed <-  calc_min_sim_size_RMSE(
                            metrics = metrics,
                            target_MCSE_RMSE = target_MCSE_RMSE,
                            vec_index_inner = vec_index_inner)
        
        return(list(
          mean_min_ESS = mean_min_ESS,
          metrics = metrics,
          metrics_tables = metrics_tables,
          n_sims = metrics$n_sims,
          n_excluded = metrics$overall$n_excluded,
          max_MCSE = max_MCSE,
          mcse_met = mcse_met,
          min_sims_needed = min_sims_needed,
          test = test,
          n_studies = n_studies,
          model_parameterisation = model_parameterisation
        ))
  
}




# ============================================================================================================
# Create Summary Row Function - MetaOrdDTA Version
# ============================================================================================================
create_summary_row_MetaOrd <- function(result,
                                       scenario_row
) {
  
        if (is.null(result)) {
          return(tibble())
        }
        
        metrics <- result$metrics
        
        tibble(
          # Scenario info
          scenario_id = scenario_row$scenario_id,
          test = scenario_row$index_test_MA,
          DGM = scenario_row$DGM,
          n_studies = scenario_row$n_studies,
          model_parameterisation = scenario_row$model_parameterisation,
          alpha_pi = scenario_row$alpha_pi,
          softplus = scenario_row$softplus,
          random_thresholds = scenario_row$random_thresholds,
          Dirichlet_random_effects_type = scenario_row$Dirichlet_random_effects_type,
          N_per_study_mean = scenario_row$N_per_study_mean,
          ##
          # Sample size info
          n_sims = result$n_sims,
          n_excluded = result$n_excluded,
          min_sims_needed = result$min_sims_needed$avg_required,
          max_MCSE = result$max_MCSE,
          mcse_met = result$mcse_met,
          ##
          # Average metrics across thresholds
          ## ---- RMSE:
          avg_RMSE_Se = metrics$overall$avg_RMSE_Se,
          avg_RMSE_Sp = metrics$overall$avg_RMSE_Sp,
          avg_RMSE = metrics$overall$avg_RMSE_Se + metrics$overall$avg_RMSE_Sp,
          ##
          max_RMSE_Se = metrics$overall$max_RMSE_Se,
          max_RMSE_Sp = metrics$overall$max_RMSE_Sp,
          ##
          # MCSE for RMSE
          avg_MCSE_RMSE_Se = metrics$overall$avg_MCSE_RMSE_Se,
          avg_MCSE_RMSE_Sp = metrics$overall$avg_MCSE_RMSE_Sp,
          ##
          max_MCSE_RMSE_Se = metrics$overall$max_MCSE_RMSE_Se,
          max_MCSE_RMSE_Sp = metrics$overall$max_MCSE_RMSE_Sp,
          ##
          # Bias
          avg_abs_bias_Se = metrics$overall$avg_abs_bias_Se,
          avg_abs_bias_Sp = metrics$overall$avg_abs_bias_Sp,
          avg_abs_bias = metrics$overall$avg_abs_bias_Se + metrics$overall$avg_abs_bias_Sp,
          ##
          max_abs_bias_Se = metrics$overall$max_abs_bias_Se,
          max_abs_bias_Sp = metrics$overall$max_abs_bias_Sp,
          ##
          # MCSE for Bias
          avg_MCSE_bias_Se = metrics$overall$avg_MCSE_bias_Se,
          avg_MCSE_bias_Sp = metrics$overall$avg_MCSE_bias_Sp,
          ##
          # Coverage
          avg_coverage_Se = metrics$overall$avg_coverage_prob_Se,
          avg_coverage_Sp = metrics$overall$avg_coverage_prob_Sp,
          ##
          # CI Width
          avg_width_Se = metrics$overall$avg_interval_width_Se,
          avg_width_Sp = metrics$overall$avg_interval_width_Sp,
          ##
          mean_min_ESS = metrics$mean_min_ESS
        )
  
}







# ============================================================================================================
# Process All Scenarios - MetaOrdDTA Version
# ============================================================================================================
process_all_scenarios_MetaOrd <- function(  scenario_df,
                                            output_dir,
                                            test_filter = NULL,
                                            reduction_factor = NULL,
                                            n_studies_filter = NULL,
                                            model_filter = NULL,
                                            target_MCSE_RMSE,
                                            max_seed_to_check,
                                            min_runs_needed = 30,
                                            min_ESS_thr = 1000,
                                            vec_index_inner
) {
        
        require(dplyr)
        
        # if (is.null(vec_index_inner)) {   # Set threshold ranges
        #   if (test == "GAD_2") {
        #     vec_index_inner <- 1:6
        #   } else if (test == "HADS") {
        #     # vec_index_inner <- 3:17
        #     vec_index_inner <- 5:15
        #   } else if (test == "BAI") {
        #     vec_index_inner <- 3:43
        #   }
        # }
  
        if (is.null(reduction_factor)) { 
          reduction_factor <- 1.0
        }
        
        # Filter scenarios if needed
        scenarios_to_process <- scenario_df
        
        if (!is.null(test_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(index_test_MA %in% test_filter)
        }
        
        if (!is.null(n_studies_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(n_studies %in% n_studies_filter)
        }
        
        if (!is.null(model_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(model_parameterisation %in% model_filter)
        }
        
        # Store all results
        all_results <- list()
        summary_table <- tibble()
        
        for (i in 1:nrow(scenarios_to_process)) {
          
                scenario_row <- scenarios_to_process[i, ]
                
                # cat("\n", paste(rep("=", 80), collapse = ""), "\n")
                # cat(sprintf("Processing: %s, %d studies, %s model\n", 
                #             scenario_row$index_test_MA, 
                #             scenario_row$n_studies, 
                #             scenario_row$model_parameterisation))
                # cat(paste(rep("=", 80), collapse = ""), "\n")
                
                # Process this scenario
                result <- process_single_scenario_MetaOrd(
                  reduction_factor = reduction_factor,
                  output_dir = output_dir,
                  scenario_row = scenario_row,
                  target_MCSE_RMSE = target_MCSE_RMSE,
                  max_seed_to_check = max_seed_to_check,
                  min_runs_needed = min_runs_needed,
                  min_ESS_thr = min_ESS_thr,
                  vec_index_inner = vec_index_inner)
                
                mean_min_ESS <- result$mean_min_ESS
                
                # print(paste("mean_min_ESS = ", mean_min_ESS))
           
                # Store results
                if (!is.null(result)) {
                  # if (mean_min_ESS > min_ESS_thr) {
                      result$scenario_info <- scenario_row
                      all_results[[paste0("scenario_", scenario_row$scenario_id)]] <- result
                      
                      # Add to summary table
                      summary_row <- create_summary_row_MetaOrd(result, scenario_row)
                      summary_table <- bind_rows(summary_table, summary_row)
                  # }
                }
                
                
        }
        
        summary_table <- summary_table %>% arrange(test, DGM)
        
        return(list(
          all_results = all_results,
          summary_table = summary_table
        ))
  
}








# ---- MetaOrdDTA MODEL SUMMARY ----
create_MetaOrd_summary <- function( data, 
                                    metric = "RMSE"
) {
  
        group_by_alpha <- FALSE
        
        require(dplyr)
        require(knitr)
        
        cat("\n================================================================================\n")
        cat("MetaOrdDTA MODEL SUMMARY:\n")
        cat("================================================================================\n\n")
        
        # Helper function to find statistically equivalent models based on CI overlap
        find_equivalent_models <- function(df,
                                           metric_col,
                                           mcse_col = NULL,
                                           find_worst = FALSE,
                                           z_score = 1.96) {  # 95% CI by default
          
                  if(find_worst) {
                    df <- df %>% arrange(desc(!!sym(metric_col)))  # Sort descending for worst
                  } else {
                    df <- df %>% arrange(!!sym(metric_col))  # Sort ascending for best
                  }
                  
                  if(nrow(df) == 0) return(df)
                  
                  if(!is.null(mcse_col)) {
                    # Calculate CIs for all models
                    df <- df %>%
                      mutate(
                        lower_CI = !!sym(metric_col) - z_score * !!sym(mcse_col),
                        upper_CI = !!sym(metric_col) + z_score * !!sym(mcse_col)
                      )
                    
                    best_lower <- df$lower_CI[1]
                    best_upper <- df$upper_CI[1]
                    
                    # Check which models have overlapping CIs with the best
                    equivalent <- df %>%
                      filter(lower_CI < best_upper & upper_CI > best_lower) %>%
                      select(-lower_CI, -upper_CI)
                  } else {
                    # Fallback if no MCSE column provided
                    equivalent <- df %>% filter(row_number() == 1)
                  }
                  
                  equivalent <- equivalent %>%
                    mutate(
                      rank = row_number(),
                      significant = ifelse(rank == 1 & nrow(equivalent) == 1, "*", ""),
                      equiv_group = paste0("Tied (n=", nrow(equivalent), ")")
                    )
                  
                  return(equivalent)
          
        }
        
        # Best performing scenarios per test/n_studies
        cat("\n1. BEST PERFORMING MODELS:\n")
        cat("----------------------------------------\n")
        cat("Note: * indicates statistically significant difference (non-overlapping 95% CIs)\n")
        cat("      Multiple models shown if CIs overlap with best\n\n")
        
        # cat("1a. BEST FOR SENSITIVITY (RMSE_Se in %):\n")
        best_Se <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se * 1,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se * 1,
                 Bias_Se_pct = abs(avg_abs_bias_Se) * 1,
                 MCSE_Bias_Se_pct = avg_MCSE_bias_Se * 1)
        
        # Group by test, n_studies, and optionally alpha_pi
        if (group_by_alpha) {
          best_Se <- best_Se %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_Se <- best_Se %>% group_by(test, DGM, n_studies)
        }
        
        if (metric == "RMSE") {
          
                best_Se <- best_Se %>% 
                  group_modify(~ find_equivalent_models(.x, "RMSE_Se_pct", "MCSE_Se_pct")) %>%
                  select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                         RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se, 
                         significant, equiv_group) %>%
                  rename(`Test` = test,
                         `N Studies` = n_studies,
                         `Model` = model_parameterisation,
                         `Alpha_pi` = alpha_pi,
                         `Rand_Thr` = random_thresholds,
                         `RMSE_Se(%)` = RMSE_Se_pct,
                         `MCSE(%)` = MCSE_Se_pct,
                         `|Bias_Se|(%)` = Bias_Se_pct,
                         `Coverage_Se(%)` = avg_coverage_Se,
                         Sig = significant,
                         Group = equiv_group)
          
        } else { 
                
                best_Se <- best_Se %>% 
                  group_modify(~ find_equivalent_models(.x, "Bias_Se_pct", "MCSE_Bias_Se_pct")) %>%
                  select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                         RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se, 
                         significant, equiv_group) %>%
                  rename(`Test` = test,
                         `N Studies` = n_studies,
                         `Model` = model_parameterisation,
                         `Alpha_pi` = alpha_pi,
                         `Rand_Thr` = random_thresholds,
                         `RMSE_Se(%)` = RMSE_Se_pct,
                         `MCSE(%)` = MCSE_Se_pct,
                         `|Bias_Se|(%)` = Bias_Se_pct,
                         `Coverage_Se(%)` = avg_coverage_Se,
                         Sig = significant,
                         Group = equiv_group)
          
        }
        
        # print(kable(best_Se, digits = 2, format = "simple"))
        
        # cat("\n1b. BEST FOR SPECIFICITY (RMSE_Sp in %):\n")
        best_Sp <- data %>%
          mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
                 Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1,
                 MCSE_Bias_Sp_pct = avg_MCSE_bias_Sp * 1)
        
        if (group_by_alpha) {
          best_Sp <- best_Sp %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_Sp <- best_Sp %>% group_by(test, DGM, n_studies)
        }
        
        if (metric == "RMSE") {
          
          best_Sp <- best_Sp %>% 
            group_modify(~ find_equivalent_models(.x, "RMSE_Sp_pct", "MCSE_Sp_pct")) %>%
            select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                   RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp, 
                   significant, equiv_group) %>%
            rename(`Test` = test,
                   `N Studies` = n_studies,
                   `Model` = model_parameterisation,
                   `Alpha_pi` = alpha_pi,
                   `Rand_Thr` = random_thresholds,
                   `RMSE_Sp(%)` = RMSE_Sp_pct,
                   `MCSE(%)` = MCSE_Sp_pct,
                   `|Bias_Sp|(%)` = Bias_Sp_pct,
                   `Coverage_Sp(%)` = avg_coverage_Sp,
                   Sig = significant,
                   Group = equiv_group)
          
        } else { 
          
          best_Sp <- best_Sp %>% 
            group_modify(~ find_equivalent_models(.x, "Bias_Sp_pct", "MCSE_Bias_Sp_pct")) %>%
            select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                   RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp, 
                   significant, equiv_group) %>%
            rename(`Test` = test,
                   `N Studies` = n_studies,
                   `Model` = model_parameterisation,
                   `Alpha_pi` = alpha_pi,
                   `Rand_Thr` = random_thresholds,
                   `RMSE_Sp(%)` = RMSE_Sp_pct,
                   `MCSE(%)` = MCSE_Sp_pct,
                   `|Bias_Sp|(%)` = Bias_Sp_pct,
                   `Coverage_Sp(%)` = avg_coverage_Sp,
                   Sig = significant,
                   Group = equiv_group)
          
        }
        
        # print(kable(best_Sp, digits = 2, format = "simple"))
        
        cat("\n1c. BEST OVERALL (RMSE_Se + RMSE_Sp in %):\n")
        best_overall <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp,
                 MCSE_total_pct = sqrt(MCSE_Se_pct^2 + MCSE_Sp_pct^2),
                 Bias_Se_pct = abs(avg_abs_bias_Se),
                 Bias_Sp_pct = abs(avg_abs_bias_Sp),
                 Bias_total_pct = Bias_Se_pct + Bias_Sp_pct)
        
        if (group_by_alpha) {
          best_overall <- best_overall %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_overall <- best_overall %>% group_by(test, DGM, n_studies)
        }
        
        best_overall <- best_overall %>%
          group_modify(~ find_equivalent_models(.x, "RMSE_total_pct", "MCSE_total_pct")) %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                 RMSE_total_pct, MCSE_total_pct,
                 Bias_total_pct, Bias_Se_pct, Bias_Sp_pct,
                 avg_coverage_Se, avg_coverage_Sp, 
                 significant, equiv_group) %>%
          rename(`Test` = test,
                 `DGM` = DGM,
                 `N Studies` = n_studies,
                 `Model` = model_parameterisation,
                 `Alpha_pi` = alpha_pi,
                 `Rand_Thr` = random_thresholds,
                 `RMSE_Total(%)` = RMSE_total_pct,
                 `MCSE_Total(%)` = MCSE_total_pct,
                 `|Bias_Se|(%)` = Bias_Se_pct,
                 `|Bias_Sp|(%)` = Bias_Sp_pct,
                 `|Bias_Total|(%)` = Bias_total_pct,
                 `Cvg_Se(%)` = avg_coverage_Se,
                 `Cvg_Sp(%)` = avg_coverage_Sp,
                 Sig = significant,
                 Group = equiv_group)
        
        print(kable(best_overall, digits = 2, format = "simple"))
        
        
        # SECTION 2: WORST performing scenarios
        cat("\n\n2. WORST PERFORMING MODELS:\n")
        cat("----------------------------------------\n")
        cat("Note: Shows models NOT in the best group for each scenario\n\n")
        ##
        # # 2a. WORST FOR SENSITIVITY
        # cat("2a. WORST FOR SENSITIVITY (RMSE_Se in %):\n")
        # 
        # # First get the best performers
        # best_Se_names <- data %>%
        #   mutate(RMSE_Se_pct = avg_RMSE_Se,
        #          MCSE_Se_pct = avg_MCSE_RMSE_Se,
        #          Bias_Se_pct = abs(avg_abs_bias_Se),
        #          MCSE_Bias_Se_pct = avg_MCSE_bias_Se)
        # 
        # if (group_by_alpha) {
        #   best_Se_names <- best_Se_names %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   best_Se_names <- best_Se_names %>% group_by(test, DGM, n_studies)
        # }
        # 
        # if (metric == "RMSE") {
        #   best_Se_names <- best_Se_names %>%
        #     group_modify(~ find_equivalent_models(.x, "RMSE_Se_pct", "MCSE_Se_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # } else {
        #   best_Se_names <- best_Se_names %>%
        #     group_modify(~ find_equivalent_models(.x, "Bias_Se_pct", "MCSE_Bias_Se_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # }
        # 
        # # Now get worst as complement
        # worst_Se <- data %>%
        #   mutate(RMSE_Se_pct = avg_RMSE_Se,
        #          MCSE_Se_pct = avg_MCSE_RMSE_Se,
        #          Bias_Se_pct = abs(avg_abs_bias_Se)) %>%
        #   anti_join(best_Se_names, by = c("test", "DGM", "n_studies", "model_parameterisation", "alpha_pi"))
        # 
        # if (group_by_alpha) {
        #   worst_Se <- worst_Se %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   worst_Se <- worst_Se %>% group_by(test, DGM, n_studies)
        # }
        # 
        # worst_Se <- worst_Se %>%
        #   arrange(test, n_studies, desc(RMSE_Se_pct)) %>%
        #   select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
        #          RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se) %>%
        #   rename(`Test` = test,
        #          `N Studies` = n_studies,
        #          `Model` = model_parameterisation,
        #          `Alpha_pi` = alpha_pi,
        #          `Rand_Thr` = random_thresholds,
        #          `RMSE_Se(%)` = RMSE_Se_pct,
        #          `MCSE(%)` = MCSE_Se_pct,
        #          `|Bias_Se|(%)` = Bias_Se_pct,
        #          `Coverage_Se(%)` = avg_coverage_Se)
        ##
        # if(nrow(worst_Se) > 0) {
        #   print(kable(worst_Se, digits = 2, format = "simple"))
        # } else {
        #   cat("(All models tied as best for all scenarios)\n")
        # # }
        # 
        # # 2b. WORST FOR SPECIFICITY
        # cat("\n2b. WORST FOR SPECIFICITY (RMSE_Sp in %):\n")
        # 
        # best_Sp_names <- data %>%
        #   mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
        #          MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
        #          Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1,
        #          MCSE_Bias_Sp_pct = avg_MCSE_bias_Sp * 1)
        # 
        # if (group_by_alpha) {
        #   best_Sp_names <- best_Sp_names %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   best_Sp_names <- best_Sp_names %>% group_by(test, DGM, n_studies)
        # }
        # 
        # if (metric == "RMSE") {
        #   best_Sp_names <- best_Sp_names %>%
        #     group_modify(~ find_equivalent_models(.x, "RMSE_Sp_pct", "MCSE_Sp_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # } else {
        #   best_Sp_names <- best_Sp_names %>%
        #     group_modify(~ find_equivalent_models(.x, "Bias_Sp_pct", "MCSE_Bias_Sp_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # }
        # 
        # worst_Sp <- data %>%
        #   mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
        #          MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
        #          Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1) %>%
        #   anti_join(best_Sp_names, by = c("test",  "DGM", "n_studies", "model_parameterisation", "alpha_pi"))
        # 
        # if (group_by_alpha) {
        #   worst_Sp <- worst_Sp %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   worst_Sp <- worst_Sp %>% group_by(test, DGM, n_studies)
        # }
        # 
        # worst_Sp <- worst_Sp %>%
        #   arrange(test, n_studies, desc(RMSE_Sp_pct)) %>%
        #   select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
        #          RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp) %>%
        #   rename(`Test` = test,
        #          `N Studies` = n_studies,
        #          `Model` = model_parameterisation,
        #          `Alpha_pi` = alpha_pi,
        #          `Rand_Thr` = random_thresholds,
        #          `RMSE_Sp(%)` = RMSE_Sp_pct,
        #          `MCSE(%)` = MCSE_Sp_pct,
        #          `|Bias_Sp|(%)` = Bias_Sp_pct,
        #          `Coverage_Sp(%)` = avg_coverage_Sp)
        
        # if(nrow(worst_Sp) > 0) {
        #   print(kable(worst_Sp, digits = 2, format = "simple"))
        # } else {
        #   cat("(All models tied as best for all scenarios)\n")
        # }
        
        # 2c. WORST OVERALL
        cat("\n2c. WORST OVERALL (RMSE_Se + RMSE_Sp in %):\n")
        
        # Get the list of best model identifiers BEFORE renaming
        best_overall_names <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp,
                 MCSE_total_pct = sqrt(MCSE_Se_pct^2 + MCSE_Sp_pct^2)) %>%
          group_by(test, DGM, n_studies) %>%
          group_modify(~ find_equivalent_models(.x, "RMSE_total_pct", "MCSE_total_pct")) %>%
          ungroup() %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds)
        
        ## ---- Now get ALL models NOT in that list
        worst_overall <- data %>%
          anti_join(best_overall_names, by = c("test", "DGM", "n_studies", 
                                               "model_parameterisation", 
                                               "alpha_pi", 
                                               "random_thresholds")) %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 Bias_Se_pct = abs(avg_abs_bias_Se),
                 Bias_Sp_pct = abs(avg_abs_bias_Sp),
                 Bias_total_pct = Bias_Se_pct + Bias_Sp_pct,
                 MCSE_total_pct = sqrt(avg_MCSE_RMSE_Se^2 + avg_MCSE_RMSE_Sp^2)) %>%
          arrange(desc(RMSE_total_pct)) %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                 RMSE_total_pct, MCSE_total_pct, Bias_total_pct, Bias_Se_pct, Bias_Sp_pct,
                 avg_coverage_Se, avg_coverage_Sp) %>%
          rename(`Test` = test, 
                 `DGM` = DGM,
                 `N Studies` = n_studies, 
                 `Model` = model_parameterisation,
                 `Alpha_pi` = alpha_pi, 
                 `Rand_Thr` = random_thresholds,
                 `RMSE_Total(%)` = RMSE_total_pct, 
                 `MCSE_Total(%)` = MCSE_total_pct,
                 `|Bias_Total|(%)` = Bias_total_pct, 
                 `|Bias_Se|(%)` = Bias_Se_pct,
                 `|Bias_Sp|(%)` = Bias_Sp_pct, 
                 `Cvg_Se(%)` = avg_coverage_Se,
                 `Cvg_Sp(%)` = avg_coverage_Sp)
        
        if(nrow(worst_overall) > 0) {
          print(kable(worst_overall, digits = 2, format = "simple"))
        } else {
          cat("(All models tied as best)\n")
        }
        
}






