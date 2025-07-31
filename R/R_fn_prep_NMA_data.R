


#' expand_test_data_to_all_studies
#' @keywords internal
#' @export
expand_test_data_to_all_studies <- function(x, 
                                            indicator_index_test_in_study, 
                                            study_mappings = NULL,
                                            missing_value = -1) {
            
          n_tests <- length(x)
          n_studies <- nrow(indicator_index_test_in_study)
          
          # Check if x is already in expanded format
          # This means all matrices should have n_studies rows
          is_already_expanded <- TRUE
          for (t in 1:n_tests) {
            # Check both diseased and non-diseased matrices
            for (c in 1:2) {
              if (nrow(x[[t]][[c]]) != n_studies) {
                is_already_expanded <- FALSE
                break
              }
            }
            if (!is_already_expanded) break
          }
          
          # If already expanded, just return x as is
          if (is_already_expanded) {
            n_studies_per_test <- colSums(indicator_index_test_in_study)
            return(list(x_expanded = x, 
                        n_studies_per_test = n_studies_per_test,
                        study_mappings = NULL))  # or you might want to preserve/calculate study_mappings
          }
          
          # Otherwise, proceed with expansion as before
          # Calculate n_studies_per_test from the indicator matrix
          n_studies_per_test <- colSums(indicator_index_test_in_study)
          
          # If study_mappings not provided, create them
          if (is.null(study_mappings)) {
            # For each test, figure out which global study indices correspond to its rows
            study_mappings <- list()
            
            for (t in 1:n_tests) {
              # Get indices of studies that have this test
              studies_with_test <- which(indicator_index_test_in_study[, t] == 1)
              
              # Assumption: rows in x[[t]] correspond to these studies in order
              n_rows_test <- nrow(x[[t]][[1]])
              
              if (length(studies_with_test) != n_rows_test) {
                warning(paste("Test", t, "has", n_rows_test, "rows but indicator shows", 
                              length(studies_with_test), "studies"))
              }
              
              # Map row indices to global study indices
              study_mappings[[t]] <- studies_with_test[1:n_rows_test]
            }
          }
          
          # Create expanded x with all studies
          x_expanded <- list()
          
          for (t in 1:n_tests) {
            x_expanded[[t]] <- list()
            
            # Get dimensions for this test
            n_cols <- ncol(x[[t]][[1]])
            
            for (c in 1:2) {  # diseased/non-diseased
              # Initialize with "missing_value" for all studies
              x_expanded[[t]][[c]] <- matrix(missing_value, nrow = n_studies, ncol = n_cols)
              
              # Fill in data for studies that have this test
              for (row_idx in 1:nrow(x[[t]][[c]])) {
                global_study_idx <- study_mappings[[t]][row_idx]
                x_expanded[[t]][[c]][global_study_idx, ] <- x[[t]][[c]][row_idx, ]
              }
            }
          }
          
          return(list(x_expanded = x_expanded, 
                      n_studies_per_test = n_studies_per_test,
                      study_mappings = study_mappings))
  
}

 





#' R_fn_convert_X_to_stan_array
#' @keywords internal
#' @export
R_fn_convert_X_to_stan_array <- function(X_list) {
  
  n_tests <- length(X_list)
  n_studies <- nrow(X_list[[1]])
  n_covariates <- ncol(X_list[[1]])
  
  # Create 3D array: [n_tests, n_studies, n_covariates]
  X_array <- array(NA, dim = c(n_tests, n_studies, n_covariates))
  
  for (t in 1:n_tests) {
    X_array[t, , ] <- X_list[[t]]
  }
  
  return(X_array)
  
}




## NOTE: X =  covariates - for NMA:
## this is a list of 2 lists X_nd and X_d - where each list in X_nd/X_d corresponds 
## to an index test, and then each element - "X_d[[c]][[t]]]" - is a matrix.
## If covs shared between 2 classes (for R&G model MUST be shared), then X is a list 
## of n_index_tests matrices instead.(i.e. a list of matrices, NOT a list of lists)

#' R_fn_prep_NMA_data
#' @keywords internal
#' @export
R_fn_prep_NMA_data <- function( x,
                                ##
                                intercept_only,
                                cov_data, 
                                ##
                                model_parameterisation,
                                ##
                                n_index_tests_per_study,
                                indicator_index_test_in_study
) {
  
  
  
    outs <- expand_test_data_to_all_studies( x = x,
                                             indicator_index_test_in_study = indicator_index_test_in_study)
    ##
    x <- outs$x_expanded
    n_studies_per_test <- outs$n_studies_per_test
    print(paste("n_studies_per_test = ", n_studies_per_test))
    ##
    if (is.null(n_index_tests_per_study)) {
      stop("if doing NMA, you must input 'n_index_tests_per_study' R vector")
    }
    if (is.null(indicator_index_test_in_study)) {
      stop(paste0("If doing NMA, you must input 'indicator_index_test_in_study' R matrix,", 
                   "\n", 
                  "(where each value is either 0 or 1, and each row in the matrix corresponds to a study,",
                  "\n", 
                  "and each column in the matrix corresponds to an index test)."))
    }
    ##
    ##
    ##
    # x_non_diseased <- x[[1]] ## for NMA this will be another list, where each list element corresponds to test t and is a matrix
    # x_diseased     <- x[[2]]
    ##
    ## ---- Find "n_tests" / "n_index_tests" for NMA:
    ##
    { 
        n_tests <- length(x)
        n_index_tests <- n_tests
        ##
        print(paste("n_tests = ", n_tests))
        print(paste("n_index_tests = ", n_index_tests))
    }
    ##
    ## find n_thr vector for NMA:
    ##
    { 
      n_thr <- c()
      for (t in 1:n_index_tests) {
         n_thr[t] <- ncol(x[[t]][[1]]) - 1
      }
      # ## Add a "1" to the start since BINARY ref test:
      # n_thr <- c(1, n_thr)
      ##
      cat("n_thr = ", n_thr)
 
    }
    ##
    ## ---- Get "n_studies":
    ##
    n_studies <- nrow(indicator_index_test_in_study)
    ## n_studies <- nrow(x[[1]][[1]])
    ##
    ##
    ##
    {
        n_cat <- n_thr + 1
        ##
        ## Make "stan_data_list":
        ##
        stan_data_list <- list()
        ##
        stan_data_list$n_studies <-  n_studies
        stan_data_list$n_index_tests   <-  n_index_tests
        ##
        stan_data_list$n_studies_per_test <- n_studies_per_test
        ##
        stan_data_list$n_thr <- n_thr  
        stan_data_list$n_thr_max <- max(n_thr)
        stan_data_list$n_cat <- n_cat 
        stan_data_list$n_cat_max <- max(n_cat)
        ##
        stan_data_list$n_index_tests_per_study <- n_index_tests_per_study
        ##
        ## Index tests only:
        ##
        stan_data_list$indicator_index_test_in_study <- indicator_index_test_in_study
        ##
        stan_data_list$x <- list()
        stan_data_list$x_with_missings <- list()
        ##
        stan_data_list$cutpoint_index <- list()
        stan_data_list$n_obs_cutpoints <- array(0, dim = c(n_index_tests, n_studies))
    }
    ##
    ## ---- Covariates:
    ##
    cov_info_list <- MetaOrdDTA:::R_fn_get_covariate_info_NMA(  intercept_only = intercept_only,
                                                                cov_data = cov_data, 
                                                                model_parameterisation = model_parameterisation,
                                                                n_index_tests = n_index_tests,
                                                                n_studies = n_studies)
    print(paste("cov_info_list$X = "))
    print(str(cov_info_list$X))
    
    X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies( X = cov_info_list$X,
                                                             indicator_index_test_in_study = indicator_index_test_in_study)
    cov_info_list$X_nd <- X[[1]]
    cov_info_list$X_d  <- X[[2]]
    try({ cov_info_list$X <- NULL })
    ##
    cov_info_list$X_nd <- R_fn_convert_X_to_stan_array(X[[1]])
    cov_info_list$X_d  <- R_fn_convert_X_to_stan_array(X[[2]])
    ##
    cat("X_nd dimensions:", dim(cov_info_list$X_nd), "\n")
    cat("X_d dimensions:", dim(cov_info_list$X_d), "\n")
    ##
    stan_data_list <- c(stan_data_list, cov_info_list)
    ##
    stan_data_list$x <- x
    stan_data_list$x_with_missings <- x
    ##
    cutpoint_index <- list()
    for (t in 1:n_index_tests) {
        cutpoint_index[[t]] <- create_cutpoint_index_matrix( data_matrix = stan_data_list$x[[t]][[1]], 
                                                             max_score = stan_data_list$n_thr[[t]])
    }
    # stan_data_list$cutpoint_index <- list(cutpoint_index, 
    #                                       cutpoint_index)
    for (t in 1:n_index_tests) {
      stan_data_list$cutpoint_index[[t]] <- list(cutpoint_index[[t]], cutpoint_index[[t]])
    }
    ##
    ##
    for (t in 1:n_index_tests) {
      
      stan_data_list$n_obs_cutpoints[t, ] <- rowSums( stan_data_list$x[[t]][[1]][, c(-1)] != -1)
      ##
      rearranged <- arrange_missing_values_paired(stan_data_list$x[[t]], 
                                                  stan_data_list$cutpoint_index[[t]])
      ##
      ##
      stan_data_list$x[[t]] <- rearranged$x
      stan_data_list$cutpoint_index[[t]] <- rearranged$cutpoint_index
        
              
    }
      ##
      ## ---- Other values needed:
      ##
      stan_data_list$cts_thr_values_nd <- list()
      stan_data_list$cts_thr_values_d  <- list()
      stan_data_list$cts_thr_values    <- list()
      ##
      for (t in 1:n_index_tests) {
        # n_thr_t <- stan_data_list$n_thr[t]
        stan_data_list$cts_thr_values_nd[[t]] <- seq(from = 1, to = stan_data_list$n_thr_max, by = 1)
        stan_data_list$cts_thr_values_d[[t]]  <- seq(from = 1, to = stan_data_list$n_thr_max, by = 1)
        stan_data_list$cts_thr_values[[t]]    <- seq(from = 1, to = stan_data_list$n_thr_max, by = 1)
      }
      ##
      ## ---- Get "n_total_cat" and "n_total_cutpoints"
      ##
      {
        stan_data_list$n_total_cat <- 0
        stan_data_list$n_total_cutpoints <- 0
        ##
        for (t in 1:n_index_tests) {
              ##
              n_thr_t <- stan_data_list$n_thr[t]
              n_cat_t <- n_thr_t + 1
              ##
              for (s in 1:n_studies_per_test[t]) {
                    
                    # if (indicator_index_test_in_study[s, t] == 1) {
                    for (k in 1:n_thr_t) {
                      stan_data_list$n_total_cutpoints = stan_data_list$n_total_cutpoints + 1;
                    }
                    for (k in 1:n_cat_t) {
                      stan_data_list$n_total_cat = stan_data_list$n_total_cat + 1;
                    }
                    # }
                
              }
          
        }
      }
      ##
      ## ---- Get "n_total_summary_cat":
      ##
      stan_data_list$n_total_summary_cat <- 0
      ##
      for (t in 1:n_index_tests) {
          for (k in 1:(stan_data_list$n_thr[t] + 1)) {
            stan_data_listn_total_summary_cat = stan_data_list$n_total_summary_cat + 1;
          }
      }
      ##
      ## ---- Print:
      ##
      {
        print(paste("stan_data_list$n_total_cutpoints = ", stan_data_list$n_total_cutpoints))
        print(paste("stan_data_list$n_total_cat = ", stan_data_list$n_total_cat))
        print(paste("stan_data_list$n_total_summary_cat = ", stan_data_list$n_total_summary_cat))
      }
    ####
    #### ---- Convert double-nested lists to 4D arrays:
    ####
      ## Initialize 4D arrays:
      x_with_missings_array <- array(-1, dim = c(n_index_tests, 2, n_studies, stan_data_list$n_thr_max + 1))
      x_array <- array(-1, dim = c(n_index_tests, 2, n_studies, stan_data_list$n_thr_max + 1))
      cutpoint_index_array <- array(-1, dim = c(n_index_tests, 2, n_studies, stan_data_list$n_thr_max + 1))
      ##
      ## Fill 4D arrays w/ list elements:
      ##
      for (t in 1:n_index_tests) {
        for (c in 1:2) {
          # Make sure to check if the data exists
          if (!is.null(stan_data_list$x_with_missings[[t]][[c]])) {
            # Get the actual dimensions of this test's data
            n_cols_t <- ncol(stan_data_list$x_with_missings[[t]][[c]])
            # Only fill the columns that exist for this test
            x_with_missings_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$x_with_missings[[t]][[c]])
          }
          if (!is.null(stan_data_list$x[[t]][[c]])) {
            n_cols_t <- ncol(stan_data_list$x[[t]][[c]])
            x_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$x[[t]][[c]])
          }
          if (!is.null(stan_data_list$cutpoint_index[[t]][[c]])) {
            n_cols_t <- ncol(stan_data_list$cutpoint_index[[t]][[c]])
            cutpoint_index_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$cutpoint_index[[t]][[c]])
          }
        }
      }
      ## Replace list structures with arrays in "stan_data_list":
      stan_data_list$x_with_missings <- x_with_missings_array
      # stan_data_list$n <- n_array
      stan_data_list$x <- x_array
      stan_data_list$cutpoint_index <- cutpoint_index_array
  
  return( list( stan_data_list = stan_data_list,
                n_tests = n_tests,
                n_index_tests = n_index_tests,
                n_studies = n_studies,
                n_studies_per_test = n_studies_per_test,
                n_thr = n_thr,
                n_cat = n_cat))
  
  
}

