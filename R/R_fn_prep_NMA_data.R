

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
                                X, 
                                model_parameterisation,
                                ##
                                n_index_tests_per_study,
                                indicator_index_test_in_study
) {
  
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
    x_non_diseased <- x[[1]] ## for NMA this will be another list, where each list element corresponds to test t and is a matrix
    x_diseased     <- x[[2]]
    ##
    ## ---- Find "n_tests" / "n_index_tests" for NMA:
    ##
    { 
        n_tests <- length(x_non_diseased)
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
         n_thr[t] <- ncol(x_non_diseased[[t]]) - 1
      }
      # ## Add a "1" to the start since BINARY ref test:
      # n_thr <- c(1, n_thr)
      ##
      cat("n_thr = ", n_thr)
 
    }
    ##
    ## ---- Get "n_studies":
    ##
    n_studies <- nrow(x_non_diseased[[1]])
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
        ##
    }
##
    for (c in 1:2) {
          stan_data_list$x[[c]] <- list()
          stan_data_list$cutpoint_index[[c]] <- list()
          ##
          for (t in 1:n_index_tests) {
              stan_data_list$x[[c]][[t]] <- array(-1, dim = c(n_studies, stan_data_list$n_thr_max + 1)) ## "staggered" array.
              stan_data_list$cutpoint_index[[c]][[t]] <- array(-1, dim = c(n_studies, stan_data_list$n_thr_max + 1)) ## "staggered" array.
          }
    }
    # array[n_index_tests] int n_covariates_nd;
    # array[n_index_tests] int n_covariates_d;
    # int n_covariates_max; // this is the max across ALL index tests
    # array[n_index_tests] matrix[n_studies, n_covariates_nd] X_nd;
    # array[n_index_tests] matrix[n_studies, n_covariates_d]  X_d;
    # array[n_index_tests] vector[n_covariates_nd] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
    # array[n_index_tests] vector[n_covariates_d]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
    ##
    ## ---- Covariates:
    ##
    cov_info_list <- MetaOrdDTA:::R_fn_get_covariate_info_NMA(  X = X, 
                                                                model_parameterisation = model_parameterisation,
                                                                n_index_tests = n_index_tests,
                                                                n_studies = n_studies)
    stan_data_list <- c(stan_data_list, cov_info_list)
    ##
    stan_data_list$x_with_missings <- x
    ##
    for (t in 1:n_index_tests) {
      
        n_thr_t <-  stan_data_list$n_thr[t]
        
        for (s in 1:n_studies) {
          
          for (c in 1:2) {
                  
                  cutpoint_counter <- 0;
                  stan_data_list$x[[c]][[t]][s, 1] <- stan_data_list$x_with_missings[[c]][[t]][s, 1]
                  prev_non_missing_count <- x[[c]][[t]][s, 1]
                  stan_data_list$cutpoint_index[[c]][[t]][s, 1] <- 0 ## because 1st col is total counts so no cutpoint to map
                  
                  for (k in 1:(n_thr_t + 1)) {
                    
                            not_missing <- FALSE
                            if  (stan_data_list$x_with_missings[[c]][[t]][s, k] != -1) { 
                              not_missing <- TRUE
                            }
                            if (k == 1) { ## skip first col. of x as this is just the total counts
                              not_missing <- FALSE
                              stan_data_list$cutpoint_index[[c]][[t]][s, 2] = 1
                            } else {
                                    
                                    if (not_missing == TRUE)   {
                                      
                                      cutpoint_counter = cutpoint_counter + 1;
                                      
                                      # if (k != (n_thr + 1)) {
                                      # if (k <= n_thr) {
                                      stan_data_list$cutpoint_index[[c]][[t]][s, cutpoint_counter + 1] = k - 1  ## - 1;
                                      # }
                                      # }
                                      
                                      stan_data_list$x[[c]][[t]][s, cutpoint_counter + 1] <- stan_data_list$x_with_missings[[c]][[t]][s, k]
                                      
                                      ##  previous_non_missing_index <- k
                                      prev_non_missing_count <- stan_data_list$x[[c]][[t]][s, cutpoint_counter + 1]
                                      
                                    } else { 
                                      ## stan_data_list$x[[c]][s, cutpoint_counter] <- prev_non_missing_count
                                    }
                              
                            }
                    
                  }
                  
                  stan_data_list$n_obs_cutpoints[t, s] <- cutpoint_counter## - 1
            
          }
          
        }
        
              
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
            for (s in 1:n_studies) {
                  
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
          if (!is.null(stan_data_list$x_with_missings[[c]][[t]])) {
            # Get the actual dimensions of this test's data
            n_cols_t <- ncol(stan_data_list$x_with_missings[[c]][[t]])
            # Only fill the columns that exist for this test
            x_with_missings_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$x_with_missings[[c]][[t]])
          }
          if (!is.null(stan_data_list$x[[c]][[t]])) {
            n_cols_t <- ncol(stan_data_list$x[[c]][[t]])
            x_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$x[[c]][[t]])
          }
          if (!is.null(stan_data_list$cutpoint_index[[c]][[t]])) {
            n_cols_t <- ncol(stan_data_list$cutpoint_index[[c]][[t]])
            cutpoint_index_array[t, c, , 1:n_cols_t] <- as.matrix(stan_data_list$cutpoint_index[[c]][[t]])
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
                n_thr = n_thr,
                n_cat = n_cat))
  
  
}

