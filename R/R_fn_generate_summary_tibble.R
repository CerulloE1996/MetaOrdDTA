

# 
# input_array = temp_trace
# iter_dim <- 1
# chain_dim <- 2
# params_dim <-  factorial(3) - (iter_dim + chain_dim)
#                                       

#' generate_summary_tibble
#' @keywords internal
#' @export
generate_summary_tibble <- function(n_threads = NULL,
                                    input_array, 
                                    iter_dim,
                                    chain_dim,
                                    param_names, 
                                    # n_to_compute, 
                                    compute_nested_rhat,
                                    # n_chains, 
                                    n_superchains) {
        
              n_threads <- if_null_then_set_to(n_threads, round(parallel::detectCores() / 2))
              
              
              
              
              ## The correct dim order (for using Bayesplot R pkg) is: Iteration, Chain, Parameter
              params_dim <-  factorial(3) - (iter_dim + chain_dim)
              reformatted_trace <-  make_array_dims_equal_Bayesplot_order(  input_array = input_array,
                                                                iter_dim = iter_dim,
                                                                chain_dim = chain_dim,
                                                                params_dim = params_dim)
              n_iter   <- dim(reformatted_trace)[iter_dim]
              n_chains <- dim(reformatted_trace)[chain_dim]
              n_params <- dim(reformatted_trace)[params_dim]

              ## Initialize summary dataframe:
              summary_df <- data.frame(     parameter = param_names,
                                            mean = NA,
                                            sd = NA,
                                            `2.5%` = NA,
                                            `50%` = NA,
                                            `97.5%` = NA,
                                            n_eff = NA,
                                            Rhat = NA,
                                            n_Rhat = NA,
                                            check.names = FALSE)
              
              # # Effective Sample Size (ESS) and Rhat - using the fast custom RcppParallel fn "BayesMVP::Rcpp_compute_MCMC_diagnostics()"
              posterior_draws_as_std_vec_of_mats <- list()
              ##
              # n_iter <- dim(reformatted_trace)[2]
              mat <- matrix(nrow = n_iter, ncol = n_chains)
              str(mat)
              for (p in 1:n_params) {
                posterior_draws_as_std_vec_of_mats[[p]] <- mat
              }
              ##
              # n_params <- n_params

             #  message(print(paste("str(reformatted_trace) = ")))
             #  message(print(paste(str(reformatted_trace))))
             # # message(print(str(posterior_draws_as_std_vec_of_mats)))

              for (p in 1:n_params) {
                  for (kk in 1:n_chains) {
                    for (i in 1:n_iter) {
                        # all_iters_param_i_chain_kk <- 1:n_iter
                        posterior_draws_as_std_vec_of_mats[[p]][i, kk] <- reformatted_trace[i, kk, p]
                    }
                  }
              }
              ##
              ##
              if (n_params < n_threads) { n_threads = n_params }

              ## Compute summary stats using custom Rcpp/C++ functions:
              outs <-  (BayesMVP:::Rcpp_compute_chain_stats(   posterior_draws_as_std_vec_of_mats,
                                                    stat_type = "mean",
                                                    n_threads = n_threads))
              means_between_chains <- outs$statistics[, 1]
              ##
              outs <-  (BayesMVP:::Rcpp_compute_chain_stats(   posterior_draws_as_std_vec_of_mats,
                                                    stat_type = "sd",
                                                    n_threads = n_threads))
              SDs_between_chains <- outs$statistics[, 1]
              ##
              outs <-  (BayesMVP:::Rcpp_compute_chain_stats(   posterior_draws_as_std_vec_of_mats,
                                                    stat_type = "quantiles",
                                                    n_threads = n_threads))
              quantiles_between_chains <- outs$statistics
          
              #### Compute Effective Sample Size (ESS) and Rhat using custom Rcpp/C++ functions
              outs <-  (BayesMVP:::Rcpp_compute_MCMC_diagnostics(  posterior_draws_as_std_vec_of_mats,
                                                        diagnostic = "split_ESS",
                                                        n_threads = n_threads))
              ess_vec <- outs$diagnostics[, 1]
              # ess_tail_vec <- outs$diagnostics[, 2]

              outs <-  (BayesMVP:::Rcpp_compute_MCMC_diagnostics(  posterior_draws_as_std_vec_of_mats,
                                                        diagnostic = "split_rhat",
                                                        n_threads = n_threads))
              rhat_vec <- outs$diagnostics[, 1]
              # rhat_tail_vec <- outs$diagnostics[, 2]
              
          

              for (p in seq_len(n_params)) {
                
                try({ 

                        #### Get all values for this parameter across iterations and chains
                        param_values <- as.vector(reformatted_trace[ , , p])

                        #### Calculate summary statistics
                        if (p %% (round(n_params/10)) == 0) {
                            message(paste("summary_df$mean[p] = "))
                            print(str(summary_df$mean[p]))
                        }
                        ##
                        ##message(paste("means_between_chains[p] = "))
                        ## print(str(means_between_chains[p]))
                        ##
                        summary_df$mean[p] <- means_between_chains[p]
                        summary_df$sd[p] <- SDs_between_chains[p]
                        try({
                            summary_df[p, c("2.5%", "50%", "97.5%")] <- quantiles_between_chains[p, ]
                            #### summary_df[p, c("2.5%", "50%", "97.5%")] <- quantiles_between_chains[, p]
                        })
                        summary_df$n_eff[p] <- round(ess_vec[p])
                        summary_df$Rhat[p] <- rhat_vec[p]
                })

              }
            
              if (compute_nested_rhat == TRUE) {

                        nested_rhat_vec <- numeric(n_params)
                        superchain_ids <- create_superchain_ids(n_chains = n_chains, n_superchains = n_superchains)

                        for (p in seq_len(n_params)) {
                          try({ 
                            nested_rhat_vec[p] <- posterior::rhat_nested(reformatted_trace[ , , p], superchain_ids = superchain_ids)
                            summary_df$n_Rhat[p] <- nested_rhat_vec[p]
                          })
                        }

              }

              summary_tibble <- tibble::tibble(summary_df)
              print(summary_tibble, n = 100)

              return(summary_tibble)

    
}





 












