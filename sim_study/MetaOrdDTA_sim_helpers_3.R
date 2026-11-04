# ============================================================================================================
# Helper Functions for MetaOrdDTA Simulation
# ============================================================================================================




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
          rest_part <- parts[2]  # "BAIDGM_ord_biv_random_C..."

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

    if (DGM == "ord_biv_random_C") {

#       true_DGM_Se <- c(95.318, 85.87, 63.02, 47.236, 27.83, 15.6)/100
      true_DGM_Se <- c(97.321, 87.998, 65.313, 48.267, 26.922, 12.975)/100
#       true_DGM_Fp <- c(59.589, 32.905, 13.731, 7.66, 3.262, 1.429)/100
      true_DGM_Fp <- c(59.623, 32.744, 13.598, 7.495, 3.036, 1.223)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    } else if (DGM == "ord_biv_fixed_C") {

#       true_DGM_Se <- c(94.064, 85.133, 63.202, 47.307, 27.227, 13.043)/100
      true_DGM_Se <- c(97.202, 88.063, 65.417, 48.687, 27.694, 13.136)/100
#       true_DGM_Fp <- c(61.364, 32.558, 12.755, 6.989, 2.867, 1.19)/100
      true_DGM_Fp <- c(60.442, 32.476, 12.869, 7.087, 2.925, 1.226)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    } else if (DGM == "Jones") {

      true_DGM_Se <- c(95.456, 83.1, 63.551, 42.164, 24.192, 12.151)/100
      true_DGM_Fp <- c(60.355, 31.206, 14.954, 6.836, 3.045, 1.334)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    } else if (DGM == "R&G_random_C") {

#       true_DGM_Se <- c(96.176, 85.168, 62.919, 47.441, 30.192, 17.019)/100
      true_DGM_Se <- c(97.543, 87.703, 66.175, 50.498, 30.885, 15.610)/100
#       true_DGM_Fp <- c(60.493, 34.353, 14.575, 7.832, 3.341, 1.277)/100
      true_DGM_Fp <- c(60.064, 33.161, 14.205, 7.818, 3.179, 1.093)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    }

  } else if (test == "HADS") {

    if (DGM == "ord_biv_random_C") {

#       true_DGM_Se <- c(98.195, 97.838, 96.432, 95.247, 93.466, 90.61, 85.748, 80.685,
#                        73.559, 63.65, 54.293, 46.918, 38.05, 28.523, 21.017, 15.068,
#                        11.544, 9.083, 7.075, 4.881, 4.042)/100
      true_DGM_Se <- c(99.573, 99.308, 97.900, 96.498, 94.533, 91.612, 86.604, 81.294,
                       73.966, 63.702, 53.997, 46.368, 37.072, 27.144, 19.208, 12.921,
                       9.295, 6.744, 4.568, 2.043, 1.161)/100
#       true_DGM_Fp <- c(92.824, 84.73, 75.338, 65.239, 54.946, 45.399, 35.677, 27.033,
#                        19.814, 14.339, 10.693, 7.51, 5.303, 3.568, 2.209, 1.449, 0.852,
#                        0.527, 0.347, 0.211, 0.139)/100
      true_DGM_Fp <- c(93.220, 85.034, 75.604, 65.467, 55.108, 45.524, 35.778, 27.103,
                       19.852, 14.371, 10.726, 7.493, 5.235, 3.463, 2.062, 1.299,
                       0.710, 0.384, 0.205, 0.078, 0.010)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    } else if (DGM == "ord_biv_fixed_C") {

      # true_DGM_Se <- c(99.491, 99.247, 97.938, 96.351, 94.374, 91.498, 86.541, 81.125,
      #                  73.816, 63.739, 54.087, 46.279, 37.224, 27.342, 19.288, 12.74,
      #                  9.258, 6.843, 4.81, 2.109, 1.208)/100
      # true_DGM_Fp <- c(93.422, 85.736, 76.484, 66.134, 56.019, 46.139, 35.456, 26.526,
      #                  19.198, 13.898, 10.345, 7.254, 4.968, 3.161, 1.807, 1.131, 0.613,
      #                  0.338, 0.18, 0.067, 0.011)/100
      ##
#       true_DGM_Se <- c(99.485, 99.206, 97.92, 96.4, 94.475, 91.68, 86.846, 81.529,
#                        74.425, 64.555, 55.086, 47.324, 38.341, 28.45, 20.373, 13.597,
#                        9.9, 7.268, 5.062, 2.249, 1.234)/100
      true_DGM_Se <- c(99.516, 99.305, 98.069, 96.590, 94.738, 92.012, 87.253, 82.025,
                       74.941, 65.059, 55.527, 47.799, 38.705, 28.677, 20.438, 13.633,
                       10.003, 7.470, 5.310, 2.357, 1.378)/100
#       true_DGM_Fp <- c(93.503, 85.892, 76.685, 66.342, 56.266, 46.369, 35.665, 26.722,
#                        19.367, 14.037, 10.447, 7.339, 5.033, 3.215, 1.851, 1.159, 0.632,
#                        0.351, 0.187, 0.071, 0.013)/100
      true_DGM_Fp <- c(93.451, 85.793, 76.565, 66.218, 56.127, 46.221, 35.525, 26.608,
                       19.269, 13.957, 10.402, 7.296, 4.998, 3.190, 1.826, 1.145,
                       0.621, 0.344, 0.185, 0.069, 0.011)/100
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

#       true_DGM_Se <- c(99.844, 99.415, 98.524, 97.028, 94.781, 91.778, 87.34, 81.7,
#                        74.824, 67.179, 60.115, 52.512, 44.642, 36.561, 28.316, 22.163,
#                        16.586, 12.465, 9.371, 6.466, 5.036)/100
      true_DGM_Se <- c(99.913, 99.610, 98.916, 97.663, 95.665, 92.859, 88.566, 82.954,
                       75.995, 68.111, 60.819, 52.898, 44.592, 35.969, 27.272, 20.776,
                       14.925, 10.574, 7.205, 3.803, 1.852)/100
#       true_DGM_Fp <- c(92.811, 84.788, 75.295, 65.295, 55.257, 45.983, 36.382, 27.891,
#                        20.59, 14.814, 10.881, 7.702, 5.227, 3.351, 1.969, 1.209, 0.698,
#                        0.416, 0.251, 0.134, 0.088)/100
      true_DGM_Fp <- c(93.257, 85.169, 75.634, 65.567, 55.448, 46.019, 36.309, 27.736,
                       20.434, 14.656, 10.767, 7.605, 5.137, 3.249, 1.876, 1.129,
                       0.628, 0.350, 0.187, 0.069, 0.024)/100
      true_DGM_Sp <- 1 - true_DGM_Fp

    }

  } else if (test == "BAI") {

    if (DGM == "ord_biv_random_C") {

#       true_DGM_Se <- c(98.165, 97.618, 96.739, 95.884, 93.859, 92.615, 91.362, 89.56,
#                        88.148, 86.218, 85.219, 82.263, 80.332, 78.382, 74.833, 70.772,
#                        68.394, 63.505, 60.846, 58.597, 56.453, 54.356, 50.815, 48.705,
#                        46.015, 44.798, 42.21, 39.824, 38.186, 37.321, 33.917, 31.312,
#                        30.039, 28.747, 27.441, 25.02, 24.331, 22.686, 20.412, 19.329,
#                        17.543, 15.165, 14.314, 13.438, 12.612, 11.389, 10.173, 9.134,
#                        8.138, 7.528, 6.925, 6.409, 5.953, 5.274, 4.622, 3.902, 3.488,
#                        3.095, 2.718, 2.283, 1.905, 1.533, 0.844)/100
      true_DGM_Se <- c(99.269, 98.740, 97.859, 96.997, 95.068, 93.874, 92.660, 91.080,
                       89.741, 87.531, 86.483, 83.520, 81.555, 79.563, 76.073, 71.950,
                       68.930, 62.995, 59.704, 57.151, 54.752, 52.430, 48.505, 46.177,
                       43.223, 41.843, 38.967, 36.231, 34.462, 33.495, 30.305, 27.844,
                       26.629, 25.424, 24.185, 21.956, 21.341, 19.782, 17.735, 16.667,
                       15.020, 12.825, 12.060, 11.259, 10.441, 9.341, 8.229, 7.296,
                       6.358, 5.803, 5.288, 4.843, 4.428, 3.811, 3.256, 2.622,
                       2.254, 1.901, 1.563, 1.234, 0.927, 0.631, 0.125)/100
      ##
#       true_DGM_Fp <- c(88.562, 80.336, 72.952, 67.435, 61.549, 56.597, 50.648, 46.35,
#                        42.087, 38.317, 34.661, 31.583, 27.836, 25.626, 23.301, 20.496,
#                        18.459, 16.608, 15.197, 14.12, 12.87, 11.969, 10.742, 9.564,
#                        8.685, 7.905, 7.29, 6.851, 6.338, 5.929, 5.313, 4.914, 4.681,
#                        4.311, 4.028, 3.567, 3.16, 2.679, 2.508, 2.347, 2.215, 1.981,
#                        1.853, 1.621, 1.344, 1.222, 0.99, 0.918, 0.856, 0.782, 0.704,
#                        0.635, 0.519, 0.462, 0.404, 0.348, 0.285, 0.238, 0.198, 0.162,
#                        0.112, 0.07, 0.05)/100
      true_DGM_Fp <- c(88.703, 79.717, 72.225, 66.131, 59.622, 54.539, 48.608, 44.363,
                       40.112, 36.862, 33.283, 30.263, 26.577, 24.470, 22.187, 19.456,
                       17.480, 15.615, 14.259, 13.175, 12.011, 11.121, 9.919, 8.813,
                       7.957, 7.196, 6.606, 6.168, 5.692, 5.299, 4.713, 4.342,
                       4.139, 3.791, 3.524, 3.095, 2.723, 2.289, 2.139, 1.994,
                       1.871, 1.662, 1.543, 1.338, 1.092, 0.981, 0.775, 0.723,
                       0.672, 0.607, 0.543, 0.489, 0.383, 0.337, 0.289, 0.245,
                       0.192, 0.156, 0.124, 0.096, 0.057, 0.027, 0.014)/100
      ##
      true_DGM_Sp <- 1 - true_DGM_Fp

    } else if (DGM == "ord_biv_fixed_C") {

#       true_DGM_Se <- c(99.011, 98.544, 97.668, 96.81, 94.69, 93.412, 92.18, 90.283,
#                        88.818, 86.889, 85.876, 82.948, 81.039, 79.021, 75.354, 71.185,
#                        68.692, 63.746, 60.947, 58.593, 56.38, 54.23, 50.576, 48.365,
#                        45.592, 44.417, 41.764, 39.272, 37.628, 36.799, 33.303, 30.636,
#                        29.333, 28.015, 26.716, 24.258, 23.65, 21.887, 19.63, 18.52,
#                        16.633, 14.203, 13.361, 12.487, 11.664, 10.422, 9.183, 8.126,
#                        7.118, 6.528, 5.96, 5.484, 5.072, 4.414, 3.776, 3.063, 2.687,
#                        2.293, 1.937, 1.566, 1.208, 0.877, 0.239)/100
      true_DGM_Se <- c(99.196, 98.891, 98.255, 97.677, 96.134, 95.264, 94.410, 93.083,
                       92.136, 90.394, 89.672, 87.234, 85.661, 84.098, 81.021, 77.457,
                       74.851, 69.426, 66.431, 64.107, 61.922, 59.723, 55.937, 53.680,
                       50.710, 49.510, 46.673, 43.970, 42.214, 41.396, 38.012, 35.436,
                       34.211, 32.982, 31.749, 29.287, 28.803, 27.058, 24.710, 23.635,
                       21.607, 18.889, 18.057, 17.207, 16.319, 14.909, 13.469, 12.320,
                       11.163, 10.525, 9.871, 9.376, 8.892, 8.065, 7.234, 6.185,
                       5.677, 5.145, 4.593, 4.018, 3.408, 2.743, 0.875)/100
      ##
#       true_DGM_Fp <- c(88.642, 80.961, 73.593, 68.007, 61.472, 56.666, 50.811, 46.624,
#                        42.504, 38.592, 34.803, 31.571, 27.857, 25.808, 23.47, 20.864,
#                        18.944, 17.064, 15.719, 14.574, 13.385, 12.421, 11.146, 9.891,
#                        8.959, 8.179, 7.577, 7.11, 6.596, 6.142, 5.475, 5.028, 4.807,
#                        4.438, 4.144, 3.652, 3.214, 2.7, 2.529, 2.369, 2.246, 2, 1.875,
#                        1.619, 1.313, 1.189, 0.949, 0.891, 0.83, 0.757, 0.681, 0.621,
#                        0.492, 0.439, 0.381, 0.327, 0.258, 0.216, 0.176, 0.139, 0.086,
#                        0.045, 0.026)/100
      true_DGM_Fp <- c(87.699, 79.661, 71.984, 66.318, 59.705, 54.853, 48.948, 44.742,
                       40.545, 37.295, 33.590, 30.404, 26.797, 24.822, 22.562, 20.027,
                       18.171, 16.409, 15.082, 13.997, 12.876, 11.920, 10.716, 9.533,
                       8.626, 7.885, 7.308, 6.885, 6.393, 5.979, 5.328, 4.908,
                       4.713, 4.373, 4.108, 3.635, 3.208, 2.714, 2.568, 2.423,
                       2.312, 2.079, 1.963, 1.715, 1.406, 1.292, 1.043, 0.997,
                       0.951, 0.888, 0.826, 0.767, 0.627, 0.578, 0.519, 0.463,
                       0.380, 0.333, 0.286, 0.241, 0.161, 0.093, 0.062)/100
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

#       true_DGM_Se <- c(99.749, 99.3, 98.626, 97.931, 96.839, 95.8, 94.261, 92.873,
#                        91.173, 89.522, 87.823, 85.811, 83.324, 81.626, 79.435, 76.254,
#                        74.014, 70.974, 69.209, 67.173, 65.147, 63.493, 60.649, 57.92,
#                        55.782, 53.744, 51.646, 50.163, 48.354, 47.687, 44.314, 42.413,
#                        41.389, 39.878, 38.482, 35.722, 35.191, 32.31, 31.091, 30.206,
#                        28.632, 25.861, 25.908, 23.858, 21.558, 20.51, 18.364, 16.689,
#                        16.35, -100, 14.202, 14.425, 12.609, 12.413, 10.388, 9.906, 9.384,
#                        -100, -100, 7.155, -100, 4.563, 3.094)/100
      true_DGM_Se <- c(99.842, 99.489, 98.923, 98.288, 97.256, 96.251, 94.714, 93.309,
                       91.622, 90.120, 88.295, 86.291, 83.731, 81.993, 79.708, 76.607,
                       74.069, 70.775, 68.445, 66.443, 64.266, 62.275, 59.192, 56.522,
                       53.957, 52.068, 49.758, 47.789, 46.008, 44.688, 41.970, 39.972,
                       38.952, 37.529, 36.222, 33.870, 32.417, 30.193, 28.657, 27.610,
                       26.270, 24.121, 23.238, 21.748, 19.992, 18.673, 16.756, 15.824,
                       14.831, 14.000, 13.196, 12.494, 11.245, 10.344, 9.358, 8.246,
                       7.279, 6.525, 5.738, 4.928, 3.816, 2.627, 1.190)/100
#       true_DGM_Fp <- c(89.437, 81.813, 74.441, 68.949, 62.382, 57.464, 51.616, 47.346,
#                        42.8, 39.076, 35.753, 32.372, 28.716, 26.552, 23.992, 20.912,
#                        18.958, 16.64, 15.362, 14.063, 12.869, 11.973, 10.545, 9.325,
#                        8.471, 7.695, 6.949, 6.505, 5.904, 5.735, 4.836, 4.37, 4.139,
#                        3.808, 3.52, 2.973, 2.897, 2.404, 2.225, 2.097, 1.876, 1.528,
#                        1.528, 1.307, 1.075, 0.981, 0.789, 0.664, 0.637, -100, 0.493,
#                        0.508, 0.394, 0.393, 0.283, 0.26, 0.24, -100, -100, 0.15, -100,
#                        0.071, 0.038)/100
      true_DGM_Fp <- c(89.065, 80.850, 73.265, 67.429, 60.551, 55.476, 49.431, 45.033,
                       40.672, 37.390, 33.963, 30.734, 27.224, 25.148, 22.720, 19.869,
                       17.842, 15.542, 14.104, 12.976, 11.849, 10.899, 9.565, 8.527,
                       7.622, 7.007, 6.309, 5.758, 5.292, 4.965, 4.338, 3.914,
                       3.709, 3.436, 3.197, 2.795, 2.564, 2.233, 2.021, 1.884,
                       1.717, 1.466, 1.370, 1.216, 1.047, 0.928, 0.770, 0.698,
                       0.625, 0.568, 0.514, 0.470, 0.395, 0.345, 0.294, 0.241,
                       0.198, 0.167, 0.137, 0.109, 0.074, 0.043, 0.014)/100
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

  # Return as percentages (multiply by 100) to match the original code
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
        DGM = c("ord_biv_fixed_C", "Jones", "ord_biv_random_C", "R&G_random_C"),
        n_studies = c(10, 50),
        model_parameterisation = c("Jones", "ord_bivariate", "Bivariate", "R&G"),
        alpha_pi = c("one", "flat", "informative", "informative_from_the_literature"),
        # box_cox = c(FALSE, TRUE),
        softplus = c(FALSE),
        random_thresholds = c(FALSE, TRUE),
        stringsAsFactors = FALSE
      ) %>%
        tibble() %>%
        # Filter out invalid combinations (from the original code logic)
        filter(
          # !(model_parameterisation == "Jones" & box_cox == FALSE),
          !((DGM %in% c("Jones", "R&G_random_C")) & (alpha_pi %in% c("informative",
                                                                      "informative_from_the_literature"))),
          !((model_parameterisation %in% c("Jones", "Bivariate")) &
              (alpha_pi %in% c("informative", "informative_from_the_literature"))),
          !((model_parameterisation == "R&G") & (alpha_pi == "informative")),
          !(DGM == "ord_biv_fixed_C" & (random_thresholds == TRUE)),
          !((DGM == "Jones") & (random_thresholds == TRUE)),
          !(model_parameterisation == "Bivariate" & random_thresholds == TRUE),
          !(model_parameterisation == "R&G" & random_thresholds == FALSE),
          !(model_parameterisation %in% c("Jones", "Bivariate") & random_thresholds == TRUE)
        ) %>%
        mutate(
          scenario_id = row_number(),
          Dirichlet_random_effects_type = case_when(
            random_thresholds & model_parameterisation == "ord_bivariate" ~ "kappa",
            random_thresholds & model_parameterisation == "R&G" ~ "kappa",
            TRUE ~ "fixed"
          ),
          # Add other fixed parameters
          N_per_study_mean = 500,
          missing_thr = TRUE
        )

      return(scenario_df)

}





# Create run matrix for selected scenarios (like the LC-MVP)
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
                              ##
                              min_runs_needed = 30,
                              target_MCSE_RMSE,
                              ##
                              test,
                              DGM,
                              n_studies,
                              ##
                              min_ESS_thr = NULL,
                              vec_index_inner,
                              min_studies_per_threshold
) {

      if (is.null(min_ESS_thr)) {

              # if (test == "BAI") {
              #
              #       if (n_studies == 10) {
              #         min_ESS_thr <- 2500
              #       } else {
              #         min_ESS_thr <- 250
              #       }
              #
              # } else {
              #
              #       if (n_studies == 10) {
              #         min_ESS_thr <- 4000
              #       } else {
              #         min_ESS_thr <- 500
              #       }
              #
              # }
              ##
              if (test == "BAI") {
                  if (n_studies == 10) {
                     min_ESS_thr <- 1000
                  } else {
                     min_ESS_thr <- 250
                  }
              } else if (test == "HADS") {
                  if (n_studies == 10) {
                    min_ESS_thr <- 3000
                  } else {
                    min_ESS_thr <- 500
                  }
              } else {  ## GAD-2
                  if (n_studies == 10) {
                    min_ESS_thr <- 3000
                  } else {
                    min_ESS_thr <- 400 ## ???
                  }
              }

      }
      ##
      true_vals <- get_true_values( test = test,
                                    DGM = DGM)
      true_DGM_Se <- true_vals$Se
      true_DGM_Sp <- true_vals$Sp
      true_DGM_Fp <- true_vals$Fp
      ##
      complete_results_list <- list()
      ##
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
      ##
      complete_results_list <- complete_results_list[!sapply(complete_results_list, is.null)]
      n_sims <- length(complete_results_list)
      ##
      if (n_sims < min_runs_needed) {

          return(list(
            mcse_met = FALSE,
            n_sims = n_sims,
            mcse_value = NA,
            min_sims_needed = NA
          ))

      }
      ##
      metrics <-  compute_simulation_metrics( test = test,
                                              output_dir = output_dir,
                                              complete_results_list = complete_results_list,
                                              true_DGM_Se = true_DGM_Se,
                                              true_DGM_Sp = true_DGM_Sp,
                                              min_ESS_thr = min_ESS_thr,
                                              vec_index_inner = vec_index_inner,
                                              min_studies_per_threshold = min_studies_per_threshold)
      ##
      min_sims_needed <-  calc_min_sim_size_RMSE( metrics = metrics,
                                                  target_MCSE_RMSE = target_MCSE_RMSE,
                                                  vec_index_inner = vec_index_inner)
      ##
      metrics_tables <- create_metrics_summary_table( metrics = metrics,
                                                      vec_index_inner = vec_index_inner)
      ##
      max_MCSE_Se <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
      max_MCSE_Sp <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
      max_MCSE <- max(max_MCSE_Se, max_MCSE_Sp)
      ##
      avg_MCSE_Se <- mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
      avg_MCSE_Sp <- mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
      avg_MCSE <- max(avg_MCSE_Se, avg_MCSE_Sp)
      ##
      MCSE_val <- (max_MCSE + avg_MCSE)/2.0
      ##
      mcse_met <- (!is.na(MCSE_val) && MCSE_val <= target_MCSE_RMSE)
      # ##
      # avg_MCSE_Se <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
      # avg_MCSE_Sp <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
      # max_MCSE <- max(avg_MCSE_Se, avg_MCSE_Sp)
      # ##
      # # MCSE_RMSE_Se_Sp <- c(mean(tables$by_threshold$MCSE_RMSE_Se),
      # #                      mean(tables$by_threshold$MCSE_RMSE_Sp))
      # # max_MCSE <- max(MCSE_RMSE_Se_Sp, na.rm = TRUE)
      # ##
      # mcse_met <- (!is.na(max_MCSE) && max_MCSE <= target_MCSE_RMSE)

      # tables$overall
      # tables$overall$Value

      return(list(
        mcse_met = mcse_met,
        n_sims = n_sims,
        MCSE_val = MCSE_val,
        min_sims_needed = min_sims_needed,
        metrics = metrics,
        tables = metrics_tables
      ))

}

# metrics$overall



