


#' arrange_missing_values_paired
#' @keywords internal
#' @export
arrange_missing_values_paired <- function(x_list, 
                                          cutpoint_index_list) {
  
          result_x <- list()
          result_cutpoint <- list()
          
          for (c in 1:2) {
            x_matrix <- x_list[[c]]
            cutpoint_matrix <- cutpoint_index_list[[c]]
            
            for (row in 1:nrow(x_matrix)) {
              # Get the row data
              x_row <- x_matrix[row, ]
              cutpoint_row <- cutpoint_matrix[row, ]
              
              # Find non-missing positions (excluding first column which is total)
              non_missing_pos <- which(x_row != -1)
              missing_pos <- which(x_row == -1)
              
              # Reorder both x and cutpoint_index together
              new_order <- c(non_missing_pos, missing_pos)
              
              x_matrix[row, ] <- x_row[new_order]
              cutpoint_matrix[row, ] <- cutpoint_row[new_order]
            }
            
            result_x[[c]] <- x_matrix
            result_cutpoint[[c]] <- cutpoint_matrix
          }
          
          return(list(x = result_x, cutpoint_index = result_cutpoint))
  
}



#' R_fn_prep_MA_data
#' @keywords internal
#' @export
R_fn_prep_MA_data <- function(  x,
                                ##
                                X, ## covariates - for MA: this is a list of 2 matrices X_nd and X_d OR a matrix (X) if covs shared between 2 classes (for R&G model MUST be shared).   
                                model_parameterisation,
                                ##
                                n_index_tests_per_study,      ## dummy arg (needed for R pkg)
                                indicator_index_test_in_study ## dummy arg (needed for R pkg)
) {
    
 
    # if (is.null(compute_sim_study_metrics)) { 
    #    compute_sim_study_metrics <- 0
    #    vec_index_inner_thr <- 1
    # }
    ##
    x_non_diseased <- x[[1]]
    x_diseased     <- x[[2]]
    ##
    ## Get "n_thr":
    ##
    try({ 
        n_thr <- ncol(x_diseased) - 1
    })
    ##
    n_cat <- n_thr + 1
    ##
    ## Get "n_studies":
    ##
    n_studies <- nrow(x_diseased)
    ##
    ## Get totals:
    ##
    n_total_nd <-  x_non_diseased[, 1]
    n_total_d  <-  x_diseased[, 1]
    ##
    ## Make "stan_data_list":
    ##
    stan_data_list <- list()
    ##
    stan_data_list$n_studies <-  n_studies
    stan_data_list$n_thr     <-  n_thr
    ##
    stan_data_list$n_non_diseased <- n_total_nd
    stan_data_list$n_diseased <-     n_total_d
    ##
    # stan_data_list$n <- list()
    stan_data_list$x <- list()
    stan_data_list$x_with_missings <- list()
    ##
    stan_data_list$cutpoint_index <- list()
    stan_data_list$n_obs_cutpoints <- rep(NA, n_studies)
    # ##
    # int n_covariates_nd;
    # int n_covariates_d;
    # int n_covariates_max;
    # matrix[n_studies, n_covariates_nd] X_nd; // must be user-inoutted - study-level covariates for D-
    # matrix[n_studies, n_covariates_d]  X_d;  // must be user-inoutted - study-level covariates for D+
    # array[2] int n_covs;
    # vector[n_covariates_nd] baseline_case_nd; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
    # vector[n_covariates_d] baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
    # array[2] int n_covs;
    ##
    cov_info_list <- MetaOrdDTA:::R_fn_get_covariate_info_MA( X = X, 
                                                              model_parameterisation = model_parameterisation,
                                                              n_studies = n_studies)
    ##
    stan_data_list <- c(stan_data_list, cov_info_list)
    ##
    for (c in 1:2) {
        stan_data_list$cutpoint_index[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr + 1);
        stan_data_list$x[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr + 1);
    }
    stan_data_list$x_with_missings <- x
    ##
    stan_data_list$x <-  stan_data_list$x_with_missings
    ##
    cutpoint_index <- create_cutpoint_index_matrix(data_matrix = stan_data_list$x[[1]], 
                                                   max_score = stan_data_list$n_thr)
    stan_data_list$cutpoint_index <- list(cutpoint_index, 
                                           cutpoint_index)
    ##
    stan_data_list$n_obs_cutpoints <- rowSums(stan_data_list$x[[1]][, c(-1)] != -1)
    
    # for (s in 1:n_studies) {
    #   
    #         for (c in 1:2) {
    #             
    #             cutpoint_counter <- 0;
    #             stan_data_list$x[[c]][s, 1] <- stan_data_list$x_with_missings[[c]][s, 1]
    #             prev_non_missing_count <- x[[c]][s, 1]
    #             stan_data_list$cutpoint_index[[c]][s, 1] <- 0 ## because 1st col is total counts so no cutpoint to map
    #             
    #             for (k in 1:(n_thr + 1)) {
    #               
    #                       not_missing <- FALSE
    #                       if  (stan_data_list$x_with_missings[[c]][s, k] != -1) { 
    #                         not_missing <- TRUE
    #                       }
    #                       if (k == 1) { ## skip first col. of x as this is just the total counts
    #                         not_missing <- FALSE
    #                         stan_data_list$cutpoint_index[[c]][s, 2] = 1
    #                       } else {
    #                           
    #                           if (not_missing == TRUE)   {
    #                             
    #                                 cutpoint_counter = cutpoint_counter + 1;
    #                                 
    #                                 # if (k != (n_thr + 1)) {
    #                                 # if (k <= n_thr) {
    #                                   stan_data_list$cutpoint_index[[c]][s, cutpoint_counter + 1] = k - 1  ## - 1;
    #                                 # }
    #                                 # }
    #                                 
    #                                 stan_data_list$x[[c]][s, cutpoint_counter + 1] <- stan_data_list$x_with_missings[[c]][s, k]
    #                                 
    #                                 ##  previous_non_missing_index <- k
    #                                 prev_non_missing_count <- stan_data_list$x[[c]][s, cutpoint_counter + 1]
    #                                   
    #                           } else { 
    #                                 ## stan_data_list$x[[c]][s, cutpoint_counter] <- prev_non_missing_count
    #                           }
    #                       }
    #               
    #             }
    # 
    #                 stan_data_list$n_obs_cutpoints[s] <- cutpoint_counter## - 1
    #             
    #         }
    # }
    ##
    # stan_data_list$x <- arrange_missing_values(stan_data_list$x)
    # stan_data_list$cutpoint_index <- arrange_missing_values(stan_data_list$cutpoint_index)
    rearranged <- arrange_missing_values_paired(stan_data_list$x, stan_data_list$cutpoint_index)
    stan_data_list$x <- rearranged$x
    stan_data_list$cutpoint_index <- rearranged$cutpoint_index
    ##
    # After creating your list of matrices, convert to 3D array:
    x_array <- array(NA, dim = c(2, n_studies, n_thr + 1))
    for (c in 1:2) {
      x_array[c, , ] <- as.integer(stan_data_list$x[[c]])
    }
    stan_data_list$x <- x_array
    ##
    # Same for cutpoint_index:
    cutpoint_index_array <- array(NA, dim = c(2, n_studies, n_thr + 1))
    for (c in 1:2) {
      cutpoint_index_array[c, , ] <- as.integer(stan_data_list$cutpoint_index[[c]])
    }
    stan_data_list$cutpoint_index <- cutpoint_index_array
    ##
    ##
    stan_data_list$cts_thr_values    <- seq(from = 1, to = n_thr, by = 1)
    ##
    # stan_data_list$box_cox <-  if_null_then_set_to(box_cox, default_box_cox(cts)[[1]])
    # stan_data_list$softplus <- if_null_then_set_to(softplus, default_softplus()[[1]])
    ##
    stan_data_list$same_cutpoints_between_groups <- 0
    ##
    # stan_data_list$x_cat <- x_cat
    ##
    print(paste("n_obs_cutpoints = ")) ; print(stan_data_list$n_obs_cutpoints)
    print(paste("cutpoint_index = ")) ; print(stan_data_list$cutpoint_index)
    # print(paste("n = ")) ; print(stan_data_list$n[[1]])
    print(paste("x = ")) ; print(stan_data_list$x[[1]])
    print(paste("x_with_missings = ")) ; print(stan_data_list$x_with_missings[[1]])
    # print(paste("x_with_missings = ")) ; print(stan_data_list$x_with_missings[[1]])
    
    n_tests <- 1
    
    # stan_data_list$compute_sim_study_metrics <- compute_sim_study_metrics
    # ##
    # stan_data_list$vec_index_inner_thr <- vec_index_inner_thr
    # stan_data_list$n_inner_thr <- length(vec_index_inner_thr)
    # ##
    # # stan_data_list$vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner_thr)
    # stan_data_list$n_outer_thr <- length(stan_data_list$vec_index_outer_thr)
    
  
  return( list( stan_data_list = stan_data_list,
                n_tests = n_tests,
                n_studies = n_studies,
                n_thr = n_thr,
                n_cat = n_cat))
  
  
}


