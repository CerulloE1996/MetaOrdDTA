



#' R_fn_prep_MA_data
#' @keywords internal
#' @export
R_fn_prep_MA_data <- function(x = NULL,
                              n_index_tests_per_study = NULL,      ## dummy arg (needed for R pkg)
                              indicator_index_test_in_study = NULL ## dummy arg (needed for R pkg)
) {
    
 
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
    #
    for (c in 1:2) {
        stan_data_list$cutpoint_index[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr);
        # stan_data_list$n[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr);
        # stan_data_list$x[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr);
        # stan_data_list$x_with_missings[[c]] <- matrix(-1, nrow = n_studies, ncol = n_thr);
    }
    stan_data_list$x <- x
    stan_data_list$x_with_missings <- x
    ##
    for (s in 1:n_studies) {
      
            for (c in 1:2) {
                
                cutpoint_counter <- 0;
                previous_non_missing_index <- 1;
                
                # stan_data_list$cutpoint_index[[c]][s, 1] <- 1
                
                # for (k in 2:(n_thr + 1)) {
                  for (k in 2:(n_thr + 1)) {
                  
                      # if (k == (n_thr + 1)) {
                      #   cond <- TRUE
                      # } else {
                        cond <- (stan_data_list$x_with_missings[[c]][s, k] != -1)
                      # }
                      
                      if (cond == TRUE)   {
                        
                            cutpoint_counter = cutpoint_counter + 1;
                            
                            # if (k != (n_thr + 1)) {
                              stan_data_list$cutpoint_index[[c]][s, cutpoint_counter] = k - 1;
                            # }
                            
                            stan_data_list$x[[c]][s, cutpoint_counter] <-   stan_data_list$x_with_missings[[c]][s, previous_non_missing_index]
                            
                            # if (cutpoint_counter > 1) {
                            #   stan_data_list$n[[c]][s, cutpoint_counter] <- stan_data_list$x[[c]][s, cutpoint_counter - 1]
                            # }
                            
                            # if (k == n_thr) {
                            #   stan_data_list$n[[1]][s, n_thr] <- stan_data_list$n_non_diseased[s]
                            #   stan_data_list$n[[2]][s, n_thr] <- stan_data_list$n_diseased[s]
                            # }
                            
                            previous_non_missing_index <- k
                              
                      }
                      ## print(paste("k = ", k))
                  
                }
                stan_data_list$n_obs_cutpoints[s] <- cutpoint_counter
                ## print(paste("c = ", c))
                
            }
                  ## print(paste("s = ", s))
    }
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
  
  return( list( stan_data_list = stan_data_list,
                n_tests = n_tests,
                n_studies = n_studies,
                n_thr = n_thr,
                n_cat = n_cat))
  
  
}


