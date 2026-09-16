





#' R_fn_sample_stan_model
#' @keywords internal
#' @export
R_fn_sample_stan_model  <-    function(   debugging = FALSE,
                                          ##
                                          stan_data_list,
                                          stan_model_obj,
                                          ##
                                          n_chains,
                                          init_lists_per_chain,
                                          ##
                                          n_superchains = NULL,
                                          seed,
                                          n_burnin,
                                          n_iter,
                                          adapt_delta,
                                          max_treedepth,
                                          metric_shape,
                                          ##
                                          cmdstanr_args
) {
  
          if (debugging) message("Printing from R_fn_sample_stan_model:")
          ##
          final_stan_args <- list(
            seed = seed,
            data = stan_data_list,
            init = init_lists_per_chain,
            chains = n_chains,
            parallel_chains = n_chains,
            iter_sampling = n_iter,
            iter_warmup = n_burnin,
            max_treedepth = max_treedepth,
            adapt_delta = adapt_delta,
            metric = metric_shape
          )
          ##
          final_stan_args <- modifyList(final_stan_args, cmdstanr_args)
          ##
          {
        
                #### stan_model_obj <- cmdstanr::cmdstan_model(Stan_model_file_path) ## , force_recompile = TRUE)
                
                tictoc::tic()
                
                try({  stan_data_list$priors <- NULL }, silent = TRUE)
                
                stan_mod_samples <- do.call(stan_model_obj$sample, final_stan_args)
                ##
                # stan_mod_samples <- stan_model_obj$sample(   seed = seed,
                #                                              data = stan_data_list,
                #                                              init =   init_lists_per_chain, 
                #                                              chains = n_chains,
                #                                              parallel_chains = n_chains, 
                #                                              iter_sampling = n_iter,
                #                                              iter_warmup = n_burnin,
                #                                              max_treedepth = max_treedepth,
                #                                              adapt_delta = adapt_delta,
                #                                              metric = metric_shape,
                #                                              ...)
                
                
                try({
                  print(tictoc::toc(log = TRUE))
                  log.txt <- tictoc::tic.log(format = TRUE)
                  tictoc::tic.clearlog()
                  time_total <- unlist(log.txt)
                  ##
                  extract_numeric_string <-  stringr::str_extract(time_total, "\\d+\\.\\d+")   
                  time_total <- as.numeric(extract_numeric_string)
                })
                
                print(paste("time_total = ",  time_total))
                
                gc(reset = TRUE)
            
          }
          
          
          out_list <- list(  
            stan_mod_samples = stan_mod_samples,
            time_total = time_total)
          
          
          return(out_list)
  
}










#' R_fn_sample_model
#' @keywords internal
#' @export
R_fn_sample_model  <-    function(      debugging = FALSE,
                                        algorithm = "Stan",
                                        ##
                                        # Model_type,
                                        #  model_args_list,
                                        ##
                                        stan_data_list,
                                        stan_model_obj,
                                        ##
                                        n_chains,
                                        init_lists_per_chain,
                                        ##
                                        n_superchains = NULL,
                                        seed,
                                        n_burnin,
                                        n_iter,
                                        adapt_delta,
                                        max_treedepth,
                                        metric_shape,
                                        ##
                                        cmdstanr_args
) { 
  
                if (debugging) message("Printing from R_fn_sample_model:")
                ##
                ## ---- Sample:
                ##
                if (algorithm %in% c("Stan", "stan", "NUTS", "NUTS-HMC")) {
                  
                        stan_sampling_outs_list <- R_fn_sample_stan_model( debugging = debugging,
                                                                           ##
                                                                           stan_data_list = stan_data_list,
                                                                           stan_model_obj = stan_model_obj,
                                                                           ##
                                                                           n_chains = n_chains,
                                                                           init_lists_per_chain = init_lists_per_chain, 
                                                                           ##
                                                                           seed = seed,
                                                                           n_burnin = n_burnin,
                                                                           n_iter = n_iter,
                                                                           adapt_delta = adapt_delta,
                                                                           max_treedepth = max_treedepth,
                                                                           metric_shape = metric_shape,
                                                                           ##
                                                                           cmdstanr_args = cmdstanr_args)
                                                                           
                        
                        return(stan_sampling_outs_list)
                        
                } else { 
                  
                        ## Add BayesMVP sampling code here for more efficient sampling (for future / to-do)
                  
                }
  
  

}

















