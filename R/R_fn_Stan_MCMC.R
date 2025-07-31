


 






# stan_fitted_model <-  Stan_model_sample_output$Stan_mod_samples


#' R_fn_compute_cmdstanr_summary_tibble
#' @keywords internal
#' @export
R_fn_compute_cmdstanr_summary_tibble <- function( stan_fitted_model, 
                                                  variables_vec) {

      cmdstanr_summary_tibble <- stan_fitted_model$summary( variables = variables_vec,
                                                            "mean", "sd", 
                                                            ~ quantile(.x, probs = c(0.025,  0.50, 0.975), na.rm = TRUE), 
                                                            "rhat" , "ess_bulk", "ess_tail") %>% print()
      return(cmdstanr_summary_tibble)


}









#' R_fn_compute_L_and_epsilon_from_Stan_model
#' @keywords internal
#' @export
R_fn_compute_L_and_epsilon_from_stan_model <- function( Min_ESS,
                                                        stan_fitted_model) {
  
              # ##
              ##
              stan_modeL_metadata <- stan_fitted_model$metadata()
              ##
              n_chains <- stan_fitted_model$num_chains()
              iter_warmup <- stan_modeL_metadata$iter_warmup
              iter_sampling <- stan_modeL_metadata$iter_sampling
              ##
              save_warmup <- stan_modeL_metadata$save_warmup
              compute_burnin_HMC_params <-  ifelse(save_warmup == 1, TRUE, FALSE)
              ##
              mean_epsilon_sampling  <-   (mean(c(stan_fitted_model$sampler_diagnostics()[,,5]))) # mean # of Leapfrog steps (sampling only)
              ##
              mean_L_sampling  <-   (mean(c(2^stan_fitted_model$sampler_diagnostics()[,,1] - 1))) # mean # of Leapfrog steps (sampling only)
              ##
              L_mean_per_chain_sampling <-  c()
              L_max_per_chain_sampling  <- c()
              for (i in 1:n_chains) {
                L_mean_per_chain_sampling[i]    <-  mean ( (  2^stan_fitted_model$sampler_diagnostics()[,,1] - 1  )[1:iter_sampling,i,1]) 
                L_max_per_chain_sampling[i]       <-  max ( (  2^stan_fitted_model$sampler_diagnostics()[,,1] - 1  )[1:iter_sampling,i,1])
              }
              ##
              max_of_mean_Ls_per_chain_sampling  <-  max(L_mean_per_chain_sampling)
              stan_total_gradient_evals_sampling <-  mean_L_sampling*n_chains*iter_sampling
              ##
              if (compute_burnin_HMC_params) {
                      mean_L_burnin  <-   (mean(c(2^stan_fitted_model$sampler_diagnostics(inc_warmup = TRUE)[1:iter_warmup,,1] - 1)))  # mean L (total)
                      ##
                      L_means_per_chain_burnin <-  c()
                      L_max_per_chain_burnin   <- c()
                      for (i in 1:n_chains) {
                        L_means_per_chain_burnin[i]    <-  mean( (  2^stan_fitted_model$sampler_diagnostics(inc_warmup = TRUE)[1:iter_warmup,,1] - 1  )[ ,i,1])
                        L_max_per_chain_burnin[i]      <-  max(  (  2^stan_fitted_model$sampler_diagnostics(inc_warmup = TRUE)[1:iter_warmup,,1] - 1  )[ ,i,1]) 
                      }
                      max_of_mean_Ls_per_chain_burnin <-        max(L_means_per_chain_burnin)
                      ##
                      stan_total_gradient_evals_burnin   <- mean_L_burnin*n_chains*iter_warmup 
                      stan_total_gradient_evals          <- stan_total_gradient_evals_sampling + stan_total_gradient_evals_burnin
                      stan_min_ess_per_grad_eval         <- Min_ESS / stan_total_gradient_evals
              }
              ##
              ## Sampling phase:
              ##
              stan_min_ess_per_grad_eval_sampling       <- Min_ESS / stan_total_gradient_evals_sampling ; stan_min_ess_per_grad_eval_sampling
              ##
              ##
              sampling_time_seconds <- max(stan_fitted_model$time()$chains$sampling)
              Min_ESS_per_sec_sampling <- Min_ESS / sampling_time_seconds
              ##
              stan_grad_evals_per_sec_sampling          <- Min_ESS_per_sec_sampling / stan_min_ess_per_grad_eval_sampling
              stan_grad_evals_per_sec_sampling_div_1000 <- stan_grad_evals_per_sec_sampling/ 1000
              ##
              ## Printing:
              ##
              print(paste("Mean L across chains (sampling) = ", mean_L_sampling))
              print(paste("Max L across chains (sampling) = ", max_of_mean_Ls_per_chain_sampling))
              ##
              print(paste("Mean HMC step-size (epsilon) across chains (sampling) = ", signif(mean_epsilon_sampling, 4)))
              ##
              print(paste("ESS/grad (sampling) - adjusted = ", signif(1000 * stan_min_ess_per_grad_eval_sampling, 3)))
              print(paste("grad/sec - adjusted ",  signif(Min_ESS_per_sec_sampling  / (1000 * stan_min_ess_per_grad_eval_sampling), 3)  ))
              print(paste("grad/sec - adjusted (sampling) = ",  signif(stan_grad_evals_per_sec_sampling_div_1000, 3)))
              ##
              if (compute_burnin_HMC_params) {
                print(paste("Mean L across chains (burnin)  = ", mean_L_burnin))
                print(paste("Max L across chains (burnin)  = ", max_of_mean_Ls_per_chain_burnin))
                print(paste("ESS/grad (total) = ", signif(1000 * stan_min_ess_per_grad_eval, 3)))
              }
              ##
              ## Compute path length (tau):
              ##
              mean_tau_sampling <- mean_L_sampling * mean_epsilon_sampling
              ##
              ## Output:
              ##
              out_list_sampling <- list(  mean_epsilon_sampling = mean_epsilon_sampling,
                                          mean_L_sampling = mean_L_sampling,
                                          mean_tau_sampling = mean_tau_sampling,
                                          ##
                                          L_mean_per_chain_sampling = L_mean_per_chain_sampling,
                                          L_max_per_chain_sampling  = L_max_per_chain_sampling,
                                          ##
                                          max_of_mean_Ls_per_chain_sampling = max_of_mean_Ls_per_chain_sampling,
                                          ##
                                          stan_total_gradient_evals_sampling = stan_total_gradient_evals_sampling,
                                          stan_min_ess_per_grad_eval_sampling = stan_min_ess_per_grad_eval_sampling,
                                          stan_grad_evals_per_sec_sampling = stan_grad_evals_per_sec_sampling,
                                          stan_grad_evals_per_sec_sampling_div_1000 = stan_grad_evals_per_sec_sampling_div_1000)
              ##
              out_list_burnin_and_total <- NULL
              if (compute_burnin_HMC_params) {
                  out_list_burnin_and_total <- list(  mean_L_burnin = mean_L_burnin,
                                                      ##
                                                      L_means_per_chain_burnin = L_means_per_chain_burnin,
                                                      L_max_per_chain_burnin = L_max_per_chain_burnin,
                                                      ##
                                                      max_of_mean_Ls_per_chain_burnin = max_of_mean_Ls_per_chain_burnin,
                                                      stan_total_gradient_evals = stan_total_gradient_evals,
                                                      stan_min_ess_per_grad_eval = stan_min_ess_per_grad_eval)
              }
              ##
              return(list( out_list_sampling = out_list_sampling,
                           out_list_burnin_and_total = out_list_burnin_and_total))

}




