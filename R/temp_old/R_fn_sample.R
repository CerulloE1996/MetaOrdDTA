








#' R_fn_sample_model
#' @keywords internal
#' @export
R_fn_sample_model  <-    function(      Model_type,
                                        init_object, ## This may be updated within this function
                                        init_lists_per_chain, ## This may be updated within this function
                                        parallel_method,
                                        Stan_data_list,  ## here as may be updated
                                        model_args_list, ## here as may be updated
                                        y,  ## here as may be updated
                                        N,  ## here as may be updated
                                        ##
                                        manual_tau,
                                        tau_if_manual,
                                        ##
                                        sample_nuisance,  ## here as may be updated
                                        diffusion_HMC,
                                        partitioned_HMC,
                                        vect_type,
                                        Phi_type,
                                        inv_Phi_type,
                                        n_params_main,  ## here as may be updated
                                        n_nuisance,  ## here as may be updated
                                        n_chains_burnin,  ## here as may be updated
                                        n_chains_sampling,
                                        n_superchains,
                                        seed,
                                        n_burnin,
                                        n_iter,
                                        adapt_delta,
                                        LR_main,
                                        LR_us,
                                        n_adapt ,
                                        clip_iter,
                                        gap,
                                        metric_type_main,
                                        metric_shape_main,
                                        metric_type_nuisance,
                                        metric_shape_nuisance,
                                        max_eps_main,
                                        max_eps_us,
                                        max_L,
                                        tau_mult,
                                        ratio_M_us,
                                        ratio_M_main,
                                        interval_width_main ,
                                        interval_width_nuisance,
                                        force_autodiff,
                                        force_PartialLog,
                                        multi_attempts,
                                        n_nuisance_to_track) { 
  
                message("Printing from R_fn_sample_model:")
  
                # bookmark
                if (sample_nuisance == FALSE) {
                  n_nuisance <- 9 # dummy
                }
 
                Model_args_as_Rcpp_List <- init_object$Model_args_as_Rcpp_List
                # print(paste("Model_args_as_Rcpp_List = "))
                # print(str(Model_args_as_Rcpp_List))
                
                bs_model <- init_object$bs_model
                
                Model_args_as_Rcpp_List$n_nuisance <- n_nuisance
                
                if (Model_type != "Stan")  {
                    n_params_main <- Model_args_as_Rcpp_List$n_params_main # <- n_params_main
                    n_nuisance <-    Model_args_as_Rcpp_List$n_nuisance #<- n_nuisance
                } else if  (Model_type == "Stan")   { 
                    Model_args_as_Rcpp_List$n_params_main <- n_params_main
                    Model_args_as_Rcpp_List$n_nuisance <- n_nuisance
                }
                      
                      if (Model_type == "Stan") { 
                        ## n_class <- 0 
                        n_class <- Stan_data_list$n_class
                      } else { 
                        n_class <-  Model_args_as_Rcpp_List$Model_args_ints[2]
                      }
                      
                      print(paste("n_class = ", n_class))

                n_params <- n_nuisance + n_params_main
                print(paste("n_params = ", n_params))
                print(paste("n_params_main = ",  n_params_main))
                print(paste("n_nuisance = ",  n_nuisance))
             
                index_nuisance <- 1:n_nuisance
                index_main <- (1 + n_nuisance):n_params
                
                if (Model_type == "Stan") {
                  Model_args_as_Rcpp_List$model_so_file <-    init_object$model_so_file
                  Model_args_as_Rcpp_List$json_file_path <-   init_object$json_file_path
                } else { 
                  Model_args_as_Rcpp_List$model_so_file <-  "none" #   init_object$dummy_model_so_file
                  Model_args_as_Rcpp_List$json_file_path <- "none" #  init_object$dummy_json_file_path
                }
                
                if (Model_type == "Stan") {
                   cmdstanr::write_stan_json(data = Stan_data_list, file = Model_args_as_Rcpp_List$json_file_path)
                } else { 
                  ####  cmdstanr::write_stan_json(data = Stan_data_list, file = Model_args_as_Rcpp_List$dummy_json_file_path)
                }
                
                print(  Model_args_as_Rcpp_List$model_so_file)
                print(  Model_args_as_Rcpp_List$json_file_path)
                
                RcppParallel::setThreadOptions(numThreads = n_chains_burnin);
                
                init_burnin_object <-                BayesMVP:::init_and_run_burnin( Model_type = Model_type,
                                                                                     init_object = init_object,
                                                                                     sample_nuisance = sample_nuisance,
                                                                                     parallel_method = parallel_method,
                                                                                     Stan_data_list = Stan_data_list,
                                                                                     y = y, ## only used in C++ for manual models! (all data passed via Stan_data_list / JSON strings for .stan models!)
                                                                                     N = N, ## only used in C++ for manual models! (all data passed via Stan_data_list / JSON strings for .stan models!)
                                                                                     n_chains_burnin = n_chains_burnin,
                                                                                     n_params_main = n_params_main,
                                                                                     n_nuisance = n_nuisance,
                                                                                     ##
                                                                                     manual_tau = manual_tau,
                                                                                     tau_if_manual = tau_if_manual,
                                                                                     ##
                                                                                     diffusion_HMC = diffusion_HMC,
                                                                                     metric_type_main = metric_type_main,
                                                                                     metric_shape_main = metric_shape_main,
                                                                                     metric_type_nuisance = metric_type_nuisance,
                                                                                     metric_shape_nuisance = metric_shape_nuisance,
                                                                                     seed = seed,
                                                                                     n_burnin = n_burnin,
                                                                                     adapt_delta = adapt_delta,
                                                                                     LR_main = LR_main,
                                                                                     LR_us = LR_us,
                                                                                     n_adapt = n_adapt,
                                                                                     partitioned_HMC = partitioned_HMC,
                                                                                     clip_iter = clip_iter,
                                                                                     gap = gap,
                                                                                     max_eps_main = max_eps_main,
                                                                                     max_eps_us = max_eps_us,
                                                                                     max_L = max_L,
                                                                                     ratio_M_us = ratio_M_us,
                                                                                     ratio_M_main = ratio_M_main,
                                                                                     interval_width_main = interval_width_main,
                                                                                     interval_width_nuisance = interval_width_nuisance,
                                                                                     tau_mult = tau_mult,
                                                                                     n_nuisance_to_track = min(n_nuisance, 10),
                                                                                     force_autodiff = force_autodiff,
                                                                                     force_PartialLog = force_PartialLog,
                                                                                     multi_attempts = multi_attempts,
                                                                                     Model_args_as_Rcpp_List = Model_args_as_Rcpp_List)
        

                {

                          theta_main_vectors_all_chains_input_from_R <- init_burnin_object$theta_main_vectors_all_chains_input_from_R  # inits stored here
                          theta_us_vectors_all_chains_input_from_R <- init_burnin_object$theta_us_vectors_all_chains_input_from_R  # inits stored here
        
                          Model_args_as_Rcpp_List <-  init_burnin_object$Model_args_as_Rcpp_List
                          EHMC_args_as_Rcpp_List <-   init_burnin_object$EHMC_args_as_Rcpp_List 
                          EHMC_Metric_as_Rcpp_List <- init_burnin_object$EHMC_Metric_as_Rcpp_List
                          EHMC_burnin_as_Rcpp_List <- init_burnin_object$EHMC_burnin_as_Rcpp_List
        
                          time_burnin <- init_burnin_object$time_burnin
                          
                          n_chains_burnin <-  init_burnin_object$n_chains_burnin
                          n_burnin <-  init_burnin_object$n_burnin

                }

                {

                          post_burnin_prep_inits <-  BayesMVP:::R_fn_post_burnin_prep_for_sampling( n_chains_sampling = n_chains_sampling,
                                                                                                    n_superchains = n_superchains,
                                                                                                    n_params_main = n_params_main,
                                                                                                    n_nuisance = n_nuisance,
                                                                                                    theta_main_vectors_all_chains_input_from_R,
                                                                                                    theta_us_vectors_all_chains_input_from_R)
        
                          theta_main_vectors_all_chains_input_from_R <- post_burnin_prep_inits$theta_main_vectors_all_chains_input_from_R
                          theta_us_vectors_all_chains_input_from_R <- post_burnin_prep_inits$theta_us_vectors_all_chains_input_from_R

                }


                {

                          gc(reset = TRUE)
        
                          tictoc::tic("post-burnin timer")
        
                          if (Model_type != "Stan") {
                              Model_args_as_Rcpp_List$model_so_file <- "none"
                              Model_args_as_Rcpp_List$json_file_path <- "none"
                          }
        
        
                          RcppParallel::setThreadOptions(numThreads = n_chains_sampling); #### BOOKMARK
                          # RcppParallel::setThreadOptions(numThreads = parallel::detectCores()); #### BOOKMARK
                          
                          if (parallel_method == "OpenMP") { 
                            fn <- BayesMVP:::Rcpp_fn_OpenMP_EHMC_sampling
                          } else { ###  use RcppParallel
                            fn <- BayesMVP:::Rcpp_fn_RcppParallel_EHMC_sampling
                          }
                          
                          ### Call C++ parallel sampling function
                          result <-       (fn(                          n_threads_R = n_chains_sampling,
                                                                        sample_nuisance_R = sample_nuisance,
                                                                        n_nuisance_to_track = n_nuisance_to_track,
                                                                        seed_R = seed,
                                                                        iter_one_by_one = FALSE,
                                                                        n_iter_R = n_iter,
                                                                        partitioned_HMC_R = partitioned_HMC,
                                                                        Model_type_R = Model_type,
                                                                        force_autodiff_R = force_autodiff,
                                                                        force_PartialLog = force_PartialLog,
                                                                        multi_attempts_R = multi_attempts,
                                                                        theta_main_vectors_all_chains_input_from_R = theta_main_vectors_all_chains_input_from_R, # inits stored here
                                                                        theta_us_vectors_all_chains_input_from_R = theta_us_vectors_all_chains_input_from_R,  # inits stored here
                                                                        y =  y,  ## only used in C++ for manual models! (all data passed via Stan_data_list / JSON strings for .stan models!) 
                                                                        Model_args_as_Rcpp_List =  Model_args_as_Rcpp_List,
                                                                        EHMC_args_as_Rcpp_List =   EHMC_args_as_Rcpp_List,
                                                                        EHMC_Metric_as_Rcpp_List = EHMC_Metric_as_Rcpp_List))
                          
                          try({
                              print(tictoc::toc(log = TRUE))
                              log.txt <- tictoc::tic.log(format = TRUE)
                              tictoc::tic.clearlog()
                              time_sampling <- unlist(log.txt)
                              ##
                              extract_numeric_string <-  stringr::str_extract(time_sampling, "\\d+\\.\\d+")   
                              time_sampling <- as.numeric(extract_numeric_string)
                          })
        
                          print(paste("time_sampling = ",  time_sampling))
                          
                          gc(reset = TRUE)

                }
                
                time_total <- NULL
                try({
                   time_total <- time_sampling + time_burnin
                   print(paste("time_burnin = ",  time_burnin))
                   print(paste("time_total = ",  time_total))
                })

  out_list <- list(LR_main = LR_main, 
                   LR_us = LR_us, 
                   adapt_delta = adapt_delta,
                   n_chains_burnin = n_chains_burnin,
                   n_burnin = n_burnin,
                   metric_type_main = metric_type_main,
                   metric_shape_main = metric_shape_main,
                   metric_type_nuisance = metric_type_nuisance,
                   metric_shape_nuisance = metric_shape_nuisance,
                   diffusion_HMC = diffusion_HMC,
                   partitioned_HMC = partitioned_HMC,
                   n_superchains = n_superchains,
                   interval_width_main = interval_width_main,
                   interval_width_nuisance = interval_width_nuisance,
                   force_autodiff = force_autodiff,
                   force_PartialLog = force_PartialLog,
                   multi_attempts = multi_attempts,
                   time_burnin = time_burnin,
                   time_sampling = time_sampling,
                   time_total = time_total,
                   result = result, 
                   init_burnin_object = init_burnin_object)
  
  
  return(out_list)





}




