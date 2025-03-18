


#' generate_summary_tibble
#' @keywords internal
#' @export
generate_summary_tibble <- function(n_threads = NULL,
                                    trace, 
                                    param_names, 
                                    n_to_compute, 
                                    compute_nested_rhat,
                                    n_chains, 
                                    n_superchains) {
  
        
            n_cores <- parallel::detectCores()
            n_threads <- n_cores

              #### Initialize summary dataframe
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
              mat <- matrix(nrow = n_iter, ncol = n_chains)

              n_params <- n_to_compute
              
              comment(print(str(trace)))
             # comment(print(str(posterior_draws_as_std_vec_of_mats)))

              for (i in 1:n_params) {
                posterior_draws_as_std_vec_of_mats[[i]] <- mat
                for (kk in 1:n_chains) {
                  posterior_draws_as_std_vec_of_mats[[i]][1:n_iter, kk] <- trace[i, 1:n_iter, kk] 
                }
              }
              
              if (n_params < n_threads) { n_threads = n_params }
              
              #### Compute summary stats using custom Rcpp/C++ functions
              outs <-  (BayesMVP:::Rcpp_compute_chain_stats(   posterior_draws_as_std_vec_of_mats,
                                                    stat_type = "mean",
                                                    n_threads = n_threads))
              means_between_chains <- outs$statistics[, 1]

              outs <-  (BayesMVP:::Rcpp_compute_chain_stats(   posterior_draws_as_std_vec_of_mats,
                                                    stat_type = "sd",
                                                    n_threads = n_threads))
              SDs_between_chains <- outs$statistics[, 1]

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
              
              
              for (i in seq_len(n_to_compute)) {
                
                    #### Get all values for this parameter across iterations and chains
                    param_values <- as.vector(trace[i, , ])
                    
                    #### Calculate summary statistics
                    summary_df$mean[i] <- means_between_chains[i]
                    summary_df$sd[i] <- SDs_between_chains[i]
                    try({  
                        summary_df[i, c("2.5%", "50%", "97.5%")] <- quantiles_between_chains[i, ]
                        #### summary_df[i, c("2.5%", "50%", "97.5%")] <- quantiles_between_chains[, i]
                    })
                    summary_df$n_eff[i] <- round(ess_vec[i])
                    summary_df$Rhat[i] <- rhat_vec[i]
                
              }
              
              if (compute_nested_rhat == TRUE) {
                
                    nested_rhat_vec <- numeric(n_to_compute)
                    superchain_ids <- create_superchain_ids(n_chains = n_chains, n_superchains = n_superchains)
                    
                    for (i in seq_len(n_to_compute)) {
                      nested_rhat_vec[i] <- posterior::rhat_nested(trace[i, , ], superchain_ids = superchain_ids)
                      summary_df$n_Rhat[i] <- nested_rhat_vec[i]
                    }
                                    
              }
              
              summary_tibble <- tibble::tibble(summary_df)
              print(summary_tibble, n = 100)
              
              return(summary_tibble)

    
}





 







