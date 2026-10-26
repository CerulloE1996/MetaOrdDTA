








#' initialise_update_run_model
#' @keywords internal
#' @export
initialise_update_run_model <- function(      debugging = FALSE,
                                              ##
                                              x, 
                                              ##
                                              intercept_only,
                                              cov_data,
                                              ##
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
                                              init_lists_per_chain,
                                              ##
                                              cmdstanr_args,
                                              ##
                                              algorithm,
                                              NicoStan_args
) {
            
          if (basic_model_options$network) {
            try({ 
              n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
            }, silent = TRUE)
            if (debugging) message("n_index_tests_per_study = ", paste(n_index_tests_per_study, collapse = ", "))
          } else { 
            n_index_tests_per_study <- NULL
          }
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_1"))
          ##
          ## ---- Update R lists (if any of inputs change):
          ##
          prep_data_and_model_outs <-  prep_data_and_model(  debugging = debugging, 
                                                             ##
                                                             x = x,
                                                             ##
                                                             intercept_only = intercept_only,
                                                             cov_data = cov_data,
                                                             ##
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
                      internal_obj$outs_data$stan_data_list$compound_symmetry     <- advanced_model_options$compound_symmetry
                      internal_obj$outs_data$stan_data_list$hetero_sigma     <- advanced_model_options$hetero_sigma
                      ##
                      internal_obj$outs_data$stan_data_list <- c(internal_obj$outs_data$stan_data_list, priors) ## add priors to "stan_data_list"
          }
          ##
          ## ---- Initialise Stan model (shouldn't have to compile again):
          ##
          ## The compiler already resolved the correct MA/NMA and prior-only path.
          ## The previous code overwrote every path with the MA directory.
          internal_obj$outs_stan_model_name$stan_model_file_path <-
            internal_obj$outs_stan_compile$stan_model_file_path
          if (debugging) {
            message("Stan model: ", internal_obj$outs_stan_model_name$stan_model_file_path)
          }
          ##
          ## ---- Remove any duplicates from "stan_data_list":
          ##
          stan_data_list_temp <- internal_obj$outs_data$stan_data_list
          internal_obj$outs_data$stan_data_list <- R_fn_remove_duplicates_from_list(stan_data_list_temp)
          ##
          ## ---- prior_only = TRUE: the model's own Stan file with its likelihood switched off (every study's
          ##      K_fold_CV_indicator = 0), so that the draws are the model's real induced prior:
          ##
          if (isTRUE(as.logical(basic_model_options$prior_only))) {
                stan_model_file_path_of_the_fit <- internal_obj$outs_stan_model_name$stan_model_file_path
                model_reads_the_K_fold_CV_indicator <- any(grepl("K_fold_CV_indicator",
                                                                 readLines(stan_model_file_path_of_the_fit)))
                if (!model_reads_the_K_fold_CV_indicator) {
                      stop(paste0("prior_only = TRUE: the Stan model ", stan_model_file_path_of_the_fit,
                                  " has no K_fold_CV_indicator, so its likelihood cannot be switched off."))
                }
                n_studies_of_the_fit <- internal_obj$outs_data$n_studies
                internal_obj$outs_data$stan_data_list$K_fold_CV_indicator <- rep(0, n_studies_of_the_fit)
                message(paste0("prior_only = TRUE: the likelihood of all ", n_studies_of_the_fit,
                               " studies is switched off (K_fold_CV_indicator = 0) in ",
                               basename(stan_model_file_path_of_the_fit)))
          }
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
          #             # Or give it a name if it's required by the model
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
          {
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
                                                stan_functions_directory =
                                                      internal_obj$outs_stan_compile$stan_functions_directory,
                                                n_chains             = MCMC_params$n_chains,
                                                init_lists_per_chain = init_lists_per_chain)
          
          }
          # internal_obj$outs_data$stan_data_list$n[[1]][[1]]
          # init_lists_per_chain[[1]]
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_4"))
          ##
          ## ----------------- Sample Stan model: --------------------------------------------------------
          ##
          if (debugging) message(paste("initialise_update_run_model", "hello_there_5"))
          ##
          ## (algorithm = "NicoStan": the NicoStan version of the model with NicoSIMD's fused AVX functions,
          ##  sampled with NicoStan; see R_fn_sample_model_with_NicoStan.R)
          if (identical(algorithm, "NicoStan")) {
            internal_obj$outs_stan_sampling <- R_fn_sample_model_with_NicoStan(
                                                     stan_data_list = internal_obj$outs_data$stan_data_list,
                                                     stan_model_file_name =
                                                           internal_obj$outs_stan_model_name$stan_model_file_name,
                                                     n_chains = MCMC_params$n_chains,
                                                     init_lists_per_chain = init_lists_per_chain,
                                                     seed = MCMC_params$seed,
                                                     n_burnin = MCMC_params$n_burnin,
                                                     n_iter = MCMC_params$n_iter,
                                                     NicoStan_args = NicoStan_args)
          } else {
          if (is.null(internal_obj$outs_stan_init$stan_model_obj)) {
            stop("The Stan model has not been compiled. Recreate the model with compile = TRUE before sampling.")
          }
          internal_obj$outs_stan_sampling <- R_fn_sample_model(
                                                     stan_data_list = internal_obj$outs_data$stan_data_list,
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
                                                     metric_shape = MCMC_params$metric_shape,
                                                     ##
                                                     cmdstanr_args = cmdstanr_args)
          }

          if (debugging) message(paste("initialise_update_run_model", "hello_there_6"))
          
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
                        init_lists_per_chain = init_lists_per_chain,
                        ##
                        cov_data = cov_data,
                        ##
                        cmdstanr_args = cmdstanr_args))
        
}


 








