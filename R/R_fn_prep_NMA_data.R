



#' R_fn_prep_NMA_data
#' @keywords internal
#' @export
R_fn_prep_NMA_data <- function( x = NULL,
                                n_index_tests_per_study = NULL,
                                indicator_index_test_in_study = NULL
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
         n_thr[t] <- ncol(x_non_diseased[[t]])
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
    ## ---- Get totals (use last col from first INDEX test):
    ##
    {
        t <- 1
        n_total_nd <-  x_non_diseased[[t]][, n_thr[t]]
        n_total_d  <-  x_diseased[[t]][, n_thr[t]]
        (cat("n_total_nd = ", n_total_nd))
        (cat("n_total_d = ", n_total_d))
    }
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
        stan_data_list$n_non_diseased <- n_total_nd
        stan_data_list$n_diseased <-     n_total_d
        ##
        stan_data_list$n_total_non_diseased <- 
        ##
        stan_data_list$n_index_tests_per_study <- n_index_tests_per_study
        ##
        ## Index tests only:
        ##
        stan_data_list$indicator_index_test_in_study <- indicator_index_test_in_study
        ##
        stan_data_list$n <- list()
        stan_data_list$x <- list()
        stan_data_list$x_with_missings <- list()
        ##
        stan_data_list$cutpoint_index <- list()
        stan_data_list$n_obs_cutpoints <- array(0, dim = c(n_studies, n_index_tests))
        ##
    }
    ##
    for (c in 1:2) {
          
          stan_data_list$x_with_missings[[c]] <- list()
          stan_data_list$x[[c]] <- list()
          stan_data_list$n[[c]] <- list()
          stan_data_list$cutpoint_index[[c]] <- list()
          ##
          for (t in 1:n_index_tests) {
              stan_data_list$x_with_missings[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
              stan_data_list$x[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
              stan_data_list$n[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
              stan_data_list$cutpoint_index[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
          }
      
    }
    ##
    ## Do "x_with_missings" first:
    ##
    for (t in 1:n_index_tests) {
      
          n_thr_t <- stan_data_list$n_thr[t]
          ##
          for (s in 1:n_studies) {
            ##
            for (k in 1:n_thr_t) {
              ##
              for (c in 1:2) {
                if (c == 1) stan_data_list$x_with_missings[[c]][[t]][s, k] <- x_non_diseased[[t]][s, k]
                if (c == 2) stan_data_list$x_with_missings[[c]][[t]][s, k] <- x_diseased[[t]][s, k]
              }
              
            }
            
          }
          
    }
    ##
    ##
    # n_obs_cutpoints <-  matrix(nrow = n_studies, ncol = n_index_tests)
    ##
    for (t in 1:n_index_tests) {
              
              n_thr_t <-  stan_data_list$n_thr[t]
              ##
              for (s in 1:n_studies) {
                
                    for (c in 1:2) {
                            
                        cutpoint_counter = 0;
                        previous_non_missing_index = 1;
                        
                        stan_data_list$cutpoint_index[[c]][[t]][s, 1] <- 1
                        
                        for (k in 2:(n_thr_t + 1)) {
                          
                              if (k == (n_thr_t + 1)) {
                                cond <- TRUE
                              } else {
                                cond <- (stan_data_list$x_with_missings[[c]][[t]][s, k] != -1)
                              }
                              
                              if (cond == TRUE)   {
                                        
                                        cutpoint_counter = cutpoint_counter + 1;
                                        
                                        if (k != (n_thr_t + 1)) {
                                          stan_data_list$cutpoint_index[[c]][[t]][s, cutpoint_counter + 1] = k;
                                        }
                                        ##
                                        stan_data_list$x[[c]][[t]][s, cutpoint_counter] <- stan_data_list$x_with_missings[[c]][[t]][s, previous_non_missing_index]
                                        ##
                                        if (cutpoint_counter > 1) {
                                          stan_data_list$n[[c]][[t]][s, cutpoint_counter - 1] <-   stan_data_list$x[[c]][[t]][s, cutpoint_counter]
                                        }
                                        ##
                                        if (k == (n_thr_t )) {
                                          stan_data_list$n[[1]][[t]][s, n_thr_t] <- stan_data_list$n_non_diseased[s]
                                          stan_data_list$n[[2]][[t]][s, n_thr_t] <- stan_data_list$n_diseased[s]
                                        }
                                        
                                        previous_non_missing_index <- k
                                
                              }
                              ## print(paste("k = ", k))
                        }
                        stan_data_list$n_obs_cutpoints[s, t] <- cutpoint_counter
                        ## print(paste("c = ", c))
                    }
                    ## print(paste("s = ", s))
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
                  ####  for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
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
    {

            if (!(is.array(stan_data_list$x_with_missings))) {

              ## Initialize 4D arrays:
              x_with_missings_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
              n_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
              x_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
              cutpoint_index_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))


              ## Fill 4D arrays w/ list elements:
              for (t in 1:n_index_tests) {
                for (c in 1:2) {
                  # Make sure to check if the data exists
                  if (!is.null(stan_data_list$x_with_missings[[c]][[t]])) {
                    x_with_missings_array[t, c, , ] <- as.matrix(stan_data_list$x_with_missings[[c]][[t]])
                  }

                  if (!is.null(stan_data_list$n[[c]][[t]])) {
                    n_array[t, c, , ] <- as.matrix(stan_data_list$n[[c]][[t]])
                  }

                  if (!is.null(stan_data_list$x[[c]][[t]])) {
                    x_array[t, c, , ] <- as.matrix(stan_data_list$x[[c]][[t]])
                  }

                  if (!is.null(stan_data_list$cutpoint_index[[c]][[t]])) {
                    cutpoint_index_array[t, c, , ] <- as.matrix(stan_data_list$cutpoint_index[[c]][[t]])
                  }
                }
              }

              ## Replace list structures with arrays in "stan_data_list":
              stan_data_list$x_with_missings <- x_with_missings_array
              stan_data_list$n <- n_array
              stan_data_list$x <- x_array
              stan_data_list$cutpoint_index <- cutpoint_index_array

            }

    }
  
  return( list( stan_data_list = stan_data_list,
                n_tests = n_tests,
                n_index_tests = n_index_tests,
                n_studies = n_studies,
                n_thr = n_thr,
                n_cat = n_cat))
  
  
}

