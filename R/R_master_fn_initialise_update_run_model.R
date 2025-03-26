








#' initialise_update_run_model
#' @keywords internal
#' @export
initialise_update_run_model <- function(      debugging = FALSE,
                                              ##
                                              x, 
                                              # n_index_tests_per_study = NULL,
                                              indicator_index_test_in_study = NULL,
                                              ##
                                              internal_obj,
                                              ##
                                              basic_model_options,
                                              advanced_model_options,
                                              MCMC_params,
                                              other_advanced_options,
                                              ##
                                              priors,
                                              ##
                                              init_lists_per_chain
) {
            
          if (basic_model_options$network) {
            try({ 
              n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
            }, silent = TRUE)
            print(paste("n_index_tests_per_study = "))
            cat(n_index_tests_per_study)
          } else { 
            n_index_tests_per_study <- NULL
          }
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_1"))
          ##
          ## ---- Update R lists (if any of inputs change):
          ##
          prep_data_and_model_outs <-  prep_data_and_model( debugging = debugging, 
                                                                         ##
                                                                         x = x,
                                                                         # n_index_tests_per_study = n_index_tests_per_study,
                                                                         indicator_index_test_in_study = indicator_index_test_in_study,
                                                                         ##
                                                                         internal_obj = internal_obj,
                                                                         ##
                                                                         basic_model_options = basic_model_options,
                                                                         advanced_model_options = advanced_model_options,
                                                                         MCMC_params = MCMC_params,
                                                                         other_advanced_options = other_advanced_options,
                                                                         ##
                                                                         priors = priors,
                                                                         ##
                                                                         init_lists_per_chain = init_lists_per_chain)
          ##
          internal_obj <- prep_data_and_model_outs$internal_obj
          ##
          basic_model_options <- prep_data_and_model_outs$basic_model_options
          advanced_model_options <- prep_data_and_model_outs$advanced_model_options
          MCMC_params <- prep_data_and_model_outs$MCMC_params
          other_advanced_options <- prep_data_and_model_outs$other_advanced_options
          ##
          init_lists_per_chain <- prep_data_and_model_outs$init_lists_per_chain
          ##
          priors <- prep_data_and_model_outs$priors
          ##
          ## ---- Set MCMC params if not already set:
          ##
          {
              MCMC_params$n_chains <- if_null_then_set_to( MCMC_params$n_chains, round(parallel::detectCores()/2))
              MCMC_params$seed     <- if_null_then_set_to( MCMC_params$seed, 123) 
              ##
              MCMC_params$n_burnin <- if_null_then_set_to( MCMC_params$n_burnin, 500)
              MCMC_params$n_iter   <- if_null_then_set_to( MCMC_params$n_iter, 1000)
              ##
              MCMC_params$n_superchains <- if_null_then_set_to( MCMC_params$n_superchains, MCMC_params$n_chains)
              ##
              MCMC_params$adapt_delta   <- if_null_then_set_to( MCMC_params$adapt_delta, 0.80)
              MCMC_params$max_treedepth <- if_null_then_set_to( MCMC_params$max_treedepth, 10)
              MCMC_params$metric_shape  <- if_null_then_set_to( MCMC_params$metric_shape, "diag_e") 
          }
          
          ##
          ## ---- Prepare / cleanup / update the "stan_data_list" list for the Stan model:
          ##
          try({ 
                {
                      ##
                      if (basic_model_options$cts) { 
                        internal_obj$outs_data$stan_data_list$box_cox <-  advanced_model_options$box_cox
                      }
                      ##
                      internal_obj$outs_data$stan_data_list$softplus  <- advanced_model_options$softplus ## this is needed for ALL models
                      ##
                      internal_obj$outs_data$stan_data_list$n_tests   <- internal_obj$outs_data$n_tests
                      internal_obj$outs_data$stan_data_list$n_studies <- internal_obj$outs_data$n_studies
                      internal_obj$outs_data$stan_data_list$n_thr     <- internal_obj$outs_data$n_thr
                      internal_obj$outs_data$stan_data_list$n_cat     <- internal_obj$outs_data$n_cat
                      ##
                      internal_obj$outs_data$stan_data_list <- c(internal_obj$outs_data$stan_data_list, priors) ## add priors to "stan_data_list"
                }
                ##
          })
          # print(paste("internal_obj$outs_data = "))
          # print(internal_obj$outs_data)
          print(paste("internal_obj$outs_data$stan_data_list = "))
          print(internal_obj$outs_data$stan_data_list)
          ##
          ## ---- Initialise Stan model (shouldn't have to compile again):
          ##
          try({
                  ##
                  internal_obj$outs_stan_model_name$stan_model_file_path <-  system.file(package = "MetaOrdDTA", 
                                                                 file.path("stan_models", "stan_models_MA", 
                                                                 internal_obj$outs_stan_model_name$stan_model_file_name))
                  ##
                  if (debugging == TRUE) {
                        ##
                        ## TEMP / FOR DEBUGGING:
                        ##
                        internal_obj$stan_model_MA_path <- "/home/enzo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/stan_models/stan_models_MA"
                        
                        if (basic_model_options$prior_only) {
                          internal_obj$outs_stan_model_name$stan_model_file_path <- file.path( internal_obj$outs_stan_model_name$stan_model_MA_path,
                                                                          "prior_only",
                                                                          internal_obj$outs_stan_model_name$stan_model_file_name)
                        } else { 
                          internal_obj$outs_stan_model_name$stan_model_file_path <- file.path( internal_obj$outs_stan_model_name$stan_model_MA_path,
                                                                          internal_obj$outs_stan_model_name$stan_model_file_name)
                        }
                  }
                  print("internal_obj$outs_stan_model_name = ")
                  print(internal_obj$outs_stan_model_name)
          })
          ##
          ## ---- Remove any duplicates from "stan_data_list":
          ##
          stan_data_list_temp <- internal_obj$outs_data$stan_data_list
          internal_obj$outs_data$stan_data_list <- R_fn_remove_duplicates_from_list(stan_data_list_temp)
          # internal_obj$outs_data$stan_data_list$
          # ##
          # ##
          # ## BOOKMARK: try and make it so that the following fn ("R_fn_init_stan_inits_external") only needs to re-run
          # ## when a parameter/input is updated by the user.
          # ##
          # try({ 
          #     for (i in 1:length(init_lists_per_chain)) {
          #       
          #             # Remove the unnamed element if it's not needed
          #             if (init_lists_per_chain[[i]][[1]] == 1) { 
          #               init_lists_per_chain[[i]][[1]] <- NULL
          #             }
          #             
          #             # Or give it a name if it's required by your model
          #             # names(init_lists_per_chain[[i]])[1] <- "parameter_name"
          #     }
          # })
          ##
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_3"))
          # ##
          # # try(( init_lists_per_chain <- NULL))
          # # init_lists_per_chain[[1]]
          # init_lists_per_chain[[1]]$inits
          try({
            
                # 
                # internal_obj$outs_stan_compile$stan_model_obj$variables()$parameters
                # 
                # ( is.numeric((internal_obj$outs_data$stan_data_list$n)) == FALSE)
                # 
                # str(internal_obj$outs_data$stan_data_list$n)
                # old_n <- internal_obj$outs_data$stan_data_list$n
                # new_n <- list()
                # indicator_index_test_in_study
                # for (c in 1:2) { 
                #     new_n[[c]] <- list()
                #     ##
                #     for (t in 1:n_tests) { 
                #       try({ 
                #         new_n[[c]][[t]] <- old_n[[c]][[t]]
                #       })
                #     }
                # }
                # 
                # str(new_n)
                # 
                # any(!(is.numeric(new_n[[1]][[4]])))
                # 
                # any(!(is.numeric(new_n[[1]])))
                # 
                # 
                internal_obj$outs_stan_init <- R_fn_init_stan_inits_external(    
                                                                          stan_data_list       = internal_obj$outs_data$stan_data_list,
                                                                          stan_model_file_path = internal_obj$outs_stan_model_name$stan_model_file_path,
                                                                          stan_model_obj       = internal_obj$outs_stan_compile$stan_model_obj,
                                                                          n_chains             = MCMC_params$n_chains,
                                                                          init_lists_per_chain = init_lists_per_chain)
          
          })
          # internal_obj$outs_data$stan_data_list$n[[1]][[1]]
          # init_lists_per_chain[[1]]
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_4"))
          ##
          ## ----------------- Sample Stan model: --------------------------------------------------------
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_5"))
          ##
          try({ internal_obj$stan_data_list$priors <- NULL }, silent = TRUE)
          ##
          try({
            try({
              internal_obj$outs_stan_sampling <- R_fn_sample_stan_model(     stan_data_list = internal_obj$outs_data$stan_data_list,
                                                           stan_model_obj = internal_obj$outs_stan_init$stan_model_obj,
                                                           ##
                                                           n_chains = MCMC_params$n_chains, 
                                                           init_lists_per_chain = init_lists_per_chain,
                                                           ##
                                                           seed = MCMC_params$seed,
                                                           n_burnin = MCMC_params$n_burnin, 
                                                           n_iter = MCMC_params$n_iter, 
                                                           adapt_delta = MCMC_params$adapt_delta, 
                                                           max_treedepth = MCMC_params$max_treedepth,
                                                           metric_shape = MCMC_params$metric_shape)
                    ##
                    # internal_obj$outs_stan_sampling <- outs_stan_sampling
                    # ##
                    # internal_obj$stan_mod_samples <- outs_stan_sampling$stan_mod_samples
                    # internal_obj$time_total       <- outs_stan_sampling$time_total
                    
            })
            ##
            if (debugging) message(paste("initialise_update_run_model", "hello_there_6"))
            ##
          })
          
          try({ internal_obj$stan_data_list$priors <- NULL }, silent = TRUE)
          
          if (debugging) message(paste("initialise_update_run_model", "hello_there_7"))
          
          return(list(  
                        # full_outputs = full_outputs,
                        ##
                        internal_obj = internal_obj,
                        ##
                        basic_model_options = basic_model_options,
                        advanced_model_options = advanced_model_options,
                        MCMC_params = MCMC_params,
                        ##
                        priors = priors,
                        ##
                        init_lists_per_chain = init_lists_per_chain))
        
}


 









