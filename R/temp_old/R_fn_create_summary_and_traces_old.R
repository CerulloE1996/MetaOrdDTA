




## Helper R function to remove "log_lik" parameter from trace array 
#' R_fn_remove_log_lik_from_array
#' @keywords internal
#' @export
R_fn_remove_log_lik_from_array <- function(arr) {
  
      format(object.size(arr), units = "MB")
      
      # Find which parameters don't contain "log_lik"
      keep_params <- !grepl("log_lik", dimnames(arr)[[1]])
      
      # Create new array without the log_lik parameters
      arr_filtered <- arr[, , keep_params, drop = FALSE]
      
      return(arr_filtered)
      
}






 
#### ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' create_summary_and_traces
#' @keywords internal
create_summary_and_traces <- function(    package = "MetaOrdDTA",
                                          use_bridgestan = NULL,
                                          use_BayesMVP_for_faster_summaries = NULL,
                                          ##
                                          Stan_model_sample_output = NULL,
                                          stan_param_names_list = NULL,
                                          ##
                                          model_results = NULL,
                                          init_object = NULL,
                                          ##
                                          n_nuisance = NULL,
                                          ##
                                          compute_main_params = TRUE, # excludes nuisance params. and log-lik 
                                          compute_transformed_parameters = TRUE,
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = FALSE, 
                                          ##
                                          compute_nested_rhat = FALSE,
                                          n_superchains = NULL,
                                          save_trace_tibbles = FALSE
) {
  
  require(stringr)
  
  if (use_bridgestan) {
       require(bridgestan)
  }

  ## Start timer: 
  tictoc::tic()
  ##
  if (is.null(use_BayesMVP_for_faster_summaries)) { 
    use_BayesMVP_for_faster_summaries <- check_if_BayesMVP_R_pkg_installed(TRUE, FALSE)
  }
  if (package == "BayesMVP") {
        use_BayesMVP_for_faster_summaries <- TRUE
        if (is.null(n_nuisance)) { 
          stop("ERROR: need 'n_nuisance' if using BayesMVP")
        }
  } else if (package == "MetaOrdDTA") {
        ## get fitted stan model object (R6 class) - needed later in the fn for MetaOrdDTA.
        stan_fitted_model <- Stan_model_sample_output$Stan_mod_samples
        n_nuisance <- 0
  }
  ##
  if (is.null(use_bridgestan)) {
     if (package == "BayesMVP") {
       use_bridgestan <- TRUE
     } else if (package == "MetaOrdDTA") {
       use_bridgestan <- FALSE
     }
  }
  ##
  if (package == "BayesMVP") {
            ## Extract essential model info from "init_object" object
            Model_type <- init_object$Model_type
            sample_nuisance <- init_object$sample_nuisance
            #### y <- init_object$y
            ##
            ## Extract traces:
            ##
            main_trace <- model_results$result[[1]]
            div_trace <- model_results$result[[2]]
            nuisance_trace <- model_results$result[[3]]
            ##
            if (save_log_lik_trace == TRUE) {
                if (Model_type != "Stan") {
                  trace_log_lik_mnl_models <-  model_results$result[[6]]
                } else {
                  trace_log_lik_mnl_models <- NULL
                }
            } else {
                trace_log_lik_mnl_models <- NULL
            }
            ##
            trace_log_lik_mnl_models <- NULL ## for BayesMVP only
            ##
            #### time_burnin <- model_results$time_burnin
            try({
                time_burnin <- model_results$init_burnin_object$time_burnin
                if(is.na(time_burnin)) {
                  time_burnin <- model_results$time_burnin
                }
            })
            try({
                time_sampling <- model_results$time_sampling
            })
            try({
                time_total_wo_summaries <- time_burnin + time_sampling
            })
  }
  ##
  if (package == "MetaOrdDTA") {
        ##
        ## Make "stan_model_outputs" list (extract everything we actually need here):
        ##
        {
          stan_model_outputs <- list()
          stan_model_outputs$diagnostic_summary  <- Stan_model_sample_output$Stan_mod_samples$diagnostic_summary()
          stan_model_outputs$sampler_diagnostics <- Stan_model_sample_output$Stan_mod_samples$sampler_diagnostics()
          stan_model_outputs$draws_array         <- Stan_model_sample_output$Stan_mod_samples$draws()
        }
        ##
        ## Extract Stan traces:
        ##
        stan_draws_array         <- stan_model_outputs$draws_array
        str(stan_draws_array)
        stan_diagnostic_summary  <- stan_model_outputs$diagnostic_summary
        stan_sampler_diagnostics <- stan_model_outputs$sampler_diagnostics
        ##
        # str(stan_model_outputs$draws_array)
        ##
        main_trace_array_w_names <- stan_model_outputs$draws_array
        # str(main_trace_array_w_names)
        n_params <- print(dim(main_trace_array_w_names)[3])
        ##
        ## Extract Stan timing info:
        ##
        try({ 
          {
            # Stan_model_samples$stan_model
            ##
            stan_times <- Stan_model_sample_output$Stan_mod_samples$time()
            time_total_wo_summaries <- stan_times$total
            ##
            stan_times_per_chain_burnin <- stan_times$chains$warmup
            time_burnin <- max(stan_times_per_chain_burnin)
            ##
            stan_times_per_chain_samp   <- stan_times$chains$sampling
            time_sampling <- max(stan_times_per_chain_samp)
          }
        })
  }
  ##
  if (package == "MetaOrdDTA") {
          try({
                n_iter <- Stan_model_sample_output$n_iter
                ##
                ## Other MCMC / HMC info
                ##
                n_chains_burnin   <- n_chains  ## For MetaOrdDTA only
                n_chains_sampling <- n_chains  ## For MetaOrdDTA only
                # init_lists_per_chain <- Stan_model_samples$init_lists_per_chain  ## For MetaOrdDTA only
                ##
                n_burnin <- Stan_model_sample_output$n_burnin
                ##
                ##
                adapt_delta <- Stan_model_sample_output$adapt_delta
                max_treedepth <- Stan_model_sample_output$max_treedepth  ## For MetaOrdDTA only
                metric_shape  <- Stan_model_sample_output$metric_shape  ## For MetaOrdDTA only
                ##
                n_superchains <- Stan_model_sample_output$n_superchains
          })
          ##
          try({
               ##  L_main_during_burnin_vec <- ?????
               ##
               stan_diagnostic_summary <- stan_model_outputs$diagnostic_summary
               n_divs_per_chain      <- stan_diagnostic_summary$num_divergent ## For MetaOrdDTA only
               n_max_trees_per_chain <- stan_diagnostic_summary$num_max_treedepth ## For MetaOrdDTA only
               ebfmi_per_chain       <- stan_diagnostic_summary$ebfmi ## For MetaOrdDTA only
          })
  } else if (package == "BayesMVP") { ## For BayesMVP only:
            try({
                ## Other MCMC / HMC info:
                n_chains_burnin <- model_results$n_chains_burnin  ## for BayesMVP only
                n_burnin <- model_results$n_burnin  ## for BayesMVP only
                ##
                LR_main <- model_results$LR_main  ## for BayesMVP only
                LR_us <- model_results$LR_us  ## for BayesMVP only
                adapt_delta <- model_results$adapt_delta  ## for BayesMVP only
                ##
                metric_type_main <- model_results$metric_type_main  ## for BayesMVP only
                metric_shape_main <- model_results$metric_shape_main  ## for BayesMVP only
                metric_type_nuisance <- model_results$metric_type_nuisance  ## for BayesMVP only
                metric_shape_nuisance <- model_results$metric_shape_nuisance  ## for BayesMVP only
                ##
                diffusion_HMC <- model_results$diffusion_HMC  ## for BayesMVP only
                partitioned_HMC <- model_results$partitioned_HMC  ## for BayesMVP only
                ##
                n_superchains <- model_results$n_superchains  ## for BayesMVP only
            })
            ##
            try({
                interval_width_main <- model_results$interval_width_main  ## for BayesMVP only
                interval_width_nuisance <- model_results$interval_width_nuisance  ## for BayesMVP only
                force_autodiff <- model_results$force_autodiff  ## for BayesMVP only
                force_PartialLog <- model_results$force_PartialLog  ## for BayesMVP only
                multi_attempts <- model_results$multi_attempts  ## for BayesMVP only
                ##
                L_main_during_burnin_vec <- model_results$init_burnin_object$L_main_during_burnin_vec ## for BayesMVP only
                L_us_during_burnin_vec <- model_results$init_burnin_object$L_us_during_burnin_vec ## for BayesMVP only
                L_main_during_burnin <- model_results$init_burnin_object$L_main_during_burnin ## for BayesMVP only
                L_us_during_burnin <- model_results$init_burnin_object$L_us_during_burnin ## for BayesMVP only
            })
            ##
            n_divs <- sum(unlist(div_trace)) ## for BayesMVP only
            pct_divs <- 100 * n_divs / length(unlist(div_trace)) ## for BayesMVP only
            ##
            n_chains <- length(main_trace) ## for BayesMVP only
            ##
            n_chains_burnin <- init_object$n_chains_burnin
            n_iter <-   dim(main_trace[[1]])[2]
            ##
            n_params_main <- dim(main_trace[[1]])[1]
            print(paste("n_params_main = ", n_params_main))
            ##
            if (is.null(n_superchains)) {
              n_superchains <- round(n_chains / n_chains_burnin)
            }
  }
  ##
  try({  
      if (is.null(compute_nested_rhat)) {
        compute_nested_rhat <- FALSE
        if (n_chains > 15) {
          compute_nested_rhat <- TRUE
        }
      }
  })
  ##
  if (package == "BayesMVP") {
          
          #### Stan_model_file_path <- (file.path(pkg_dir, "inst/stan_models/PO_LC_MVP_bin.stan"))  ### TEMP
          try({   
             Stan_model_file_path <- init_object$Stan_model_file_path
          })
          ##
          if (Model_type == "Stan") {
            json_file_path <- init_object$json_file_path
            model_so_file <-  init_object$model_so_file
          } else {
            json_file_path <- init_object$dummy_json_file_path
            model_so_file <-  init_object$dummy_model_so_file
          }
          ##
          cmdstanr::write_stan_json(data = init_object$Stan_data_list, file = json_file_path)
          ##
          print(paste("init_object$json_file_path = "))
          print(init_object$json_file_path)
          ##
          print(paste("init_object$model_so_file = "))
          print(init_object$model_so_file)
          ##
          print(paste("init_object$dummy_json_file_path = "))
          print(init_object$dummy_json_file_path)
          ##
          print(paste("init_object$dummy_model_so_file = "))
          print(init_object$dummy_model_so_file)
          ##
          Sys.setenv(STAN_THREADS = "true")
          ##
          if (use_bridgestan == TRUE) {
              bs_model  <- NULL
              bs_model  <- init_object$bs_model        ## bookmark_2 (previous code)
              if (is.null(bs_model)) { 
                bs_model <- bridgestan::StanModel$new(Stan_model_file_path, data = json_file_path, 1234) # creates .so file
              }
              bs_names  <- bs_model$param_names() ## bookmark_2 (previous code)
              bs_names_inc_tp <-  (bs_model$param_names(include_tp = TRUE))
              bs_names_inc_tp_and_gq <-  (bs_model$param_names(include_tp = TRUE, include_gq = TRUE))
              ##
              pars_names <- bs_names_inc_tp_and_gq
              ##
              n_params  <- length(bs_names)
              n_params_inc_tp <- length(bs_names_inc_tp)
              n_params_inc_tp_and_gq <- length(bs_names_inc_tp_and_gq) 
              
              bs_index <-  1:n_params
              bs_index_inc_tp <-  1:n_params_inc_tp
              bs_index_inc_tp_and_gq <-  1:n_params_inc_tp_and_gq
              
              # Mow find names, indexes and N's of tp and gq ONLY
              index_tp <- setdiff(bs_index_inc_tp, bs_index)
              names_tp <- pars_names[index_tp]
              index_tp_wo_log_lik <- setdiff(index_tp, index_log_lik)
              names_tp_wo_log_lik <- pars_names[index_tp_wo_log_lik]
              ## replace names_tp etc. to be w/o log_lik (as stored in seperate array only if user chooses)
              names_tp <- names_tp_wo_log_lik
              index_tp <- names_tp_wo_log_lik
              ##
              index_gq <- setdiff(bs_index_inc_tp_and_gq, bs_index_inc_tp)
              names_gq <- pars_names[index_gq]
              ##
              n_tp <- length(names_tp)
              n_gq <- length(names_gq)
              ##
              # exclude nuisance params from summary
              index_wo_nuisance <- (n_nuisance + 1):n_par_inc_tp_and_gq
              names_wo_nuisance <- pars_names[index_wo_nuisance]
              n_params_wo_nuisance <- length(names_wo_nuisance) ; n_params_wo_nuisance
              ##
              index_wo_log_lik <-  grep("^log_lik", pars_names, invert = TRUE)
              names_wo_log_lik <-  pars_names[index_wo_log_lik]
          }
  }
  ##
  ## For EITHER R package:
  ##
  try({  
        if (is.null(stan_param_names_list)) { 
          if (package == "BayesMVP") {
                  stan_param_names_list <- init_object$stan_param_names_list
          } else { 
            stop("ERROR: 'stan_param_names_list' cannot be NULL if using MetaOrdDTA")
          }
        } else { 
        
        }
        ##
        grouped_params <- list()
        ##
        grouped_params$bs_names_main <- stan_param_names_list$parameters
        grouped_params$bs_names_tp   <- stan_param_names_list$transformed_parameters
        grouped_params$bs_names_inc_tp        <- c(grouped_params$bs_names_main, stan_param_names_list$transformed_parameters) ## new
        ##
        grouped_params$bs_names_gq   <- stan_param_names_list$generated_quantities
        grouped_params$bs_names_inc_tp_and_gq <- c(grouped_params$bs_names_inc_tp, stan_param_names_list$generated_quantities) ## new
        ##
        grouped_params$pars_names <- grouped_params$bs_names_inc_tp_and_gq
        ##
        print(paste(grouped_params$bs_names_main))
        print(paste(grouped_params$bs_names_tp))
        print(paste(grouped_params$bs_names_gq))
  })
  ##
  ##
  try({  
    
        if (package == "BayesMVP") {
              ##
              ## For BayesMVP:
              ##
              pars_names   <- init_object$init_vals_object$param_names
              bs_names_all <- pars_names
              param_names  <- bs_names_all
        } else { 
              ##
              ## Extract full ** Stan ** 3D draws array:
              ##
              # all_draws_array <- stan_model_outputs$draws_array
               str(stan_draws_array)
              ##
              # length((dimnames(stan_draws_array))$variable)
              exclude_NA <- FALSE
              ##
              stan_param_names_all <- MetaOrdDTA:::find_all_Stan_param_names(  stan_draws_array, 
                                                                  print = FALSE, 
                                                                  exclude_NA = exclude_NA)
              bs_names_all <- (stan_param_names_all)
              length(bs_names_all)
              ##
              param_names <- bs_names_all
              bs_names <- bs_names_all
        }
        ## 
        ## ---- For EITHER BayesMVP or MetaOrdDTA (if NOT using bridgestan):
        ##
        if (use_bridgestan == FALSE) {
            try({ 
                  bs_names_main <- MetaOrdDTA:::filter_param_names_string_batched(   all_param_names_vec = bs_names_all, 
                                                                        param_strings_vec = grouped_params$bs_names_main, 
                                                                        condition = "exact_match",
                                                                        print = FALSE)
                  bs_names_main <- unname(unlist(bs_names_main)) ; bs_names_main
                  n_params_main <- print(length(bs_names_main))
                  ##
                  # bs_names_all[150:200]
                  # grouped_params$bs_names_tp
                  bs_names_tp <- MetaOrdDTA:::filter_param_names_string_batched(  all_param_names_vec = bs_names_all, 
                                                                     param_strings_vec = grouped_params$bs_names_tp, 
                                                                     condition = "exact_match",
                                                                     print = FALSE)
                  bs_names_tp <- unname(unlist(bs_names_tp)) ; bs_names_tp
                  ##
                  bs_names_inc_tp <- c(bs_names_main, bs_names_tp)
                  n_par_inc_tp  <- length(bs_names_inc_tp) 
                  # length(bs_names_inc_tp)
                  ##
                  bs_names_gq <- MetaOrdDTA:::filter_param_names_string_batched(  all_param_names_vec = bs_names_all, 
                                                                     param_strings_vec = grouped_params$bs_names_gq,
                                                                     condition = "exact_match",
                                                                     print = FALSE)
                  bs_names_gq <- unname(unlist(bs_names_gq)) ; bs_names_gq
                  bs_names_inc_tp_and_gq <- c(bs_names_main, bs_names_tp, bs_names_gq)
                  # length(bs_names_inc_tp_and_gq)
                  # bs_names_inc_tp        <- c(bs_names, stan_param_names_list$transformed_parameters) ## new
                  # bs_names_inc_tp_and_gq <- c(bs_names_inc_tp, stan_param_names_list$generated_quantities) ## new
                  #### bs_names_inc_tp        <-  bs_model$param_names(include_tp = TRUE) ## bookmark_2 (previous code)
                  #### bs_names_inc_tp_and_gq <-  bs_model$param_names(include_tp = TRUE, include_gq = TRUE) ## bookmark_2 (previous code)
                  pars_names <- bs_names_inc_tp_and_gq
                  length(pars_names)
                  ####  pars_names <- init_object$param_names
                  ### index_lp  <- grep("^lp__", pars_names, invert = FALSE)
                  ##
                  ##
                  # if (init_object$param_names[1] == "lp__") { 
                  #   pars_names <- pars_names[-1]
                  # }
                  ##
                  # if (index_lp == 1) {  # if Stan model generates __lp variable (not all models will)
                  #     pars_names <- pars_names[-c(1)]
                  # } else { 
                  #     pars_names <- pars_names
                  # }
                  ##
                  #  pars_names <- bs_model$param_names(  include_tp = TRUE, include_gq = TRUE)
                  ##
                  n_par_inc_tp_and_gq <- length(pars_names) 
                  n_par_inc_tp_and_gq
                  ##
                  index_log_lik  <- grep("^log_lik", pars_names, invert = FALSE)
                  names_log_lik  <- pars_names[index_log_lik]
                  ##
                  if (length(index_log_lik) == 0) {  
                    # if (Model_type == "Stan") {  ## if log_lik doesn't exist in Stan model
                    warning("No log_lik parameter found in Stan model. Log_lik will not be computed even if save_log_lik = TRUE")
                    # }
                  }
            })
        }
  })
  ##
  ##
  try({  
      ##
      n_params  <- length(bs_names)
      n_params_inc_tp <- length(bs_names_inc_tp)
      n_params_inc_tp_and_gq <- length(bs_names_inc_tp_and_gq) 
      ##
      bs_index <-  1:n_params
      bs_index_inc_tp <-  1:n_params_inc_tp
      bs_index_inc_tp_and_gq <-  1:n_params_inc_tp_and_gq
      ##
      ## Mow find names, indexes and N's of tp and gq ONLY:
      ##
      index_tp <- setdiff(bs_index_inc_tp, bs_index)
      names_tp <- pars_names[index_tp]
      n_tp <- length(index_tp)
      ##
      index_tp_wo_log_lik <- setdiff(index_tp, index_log_lik)
      names_tp_wo_log_lik <- pars_names[index_tp_wo_log_lik]
      n_tp_wo_log_lik <- length(index_tp_wo_log_lik)
      # ## replace names_tp etc. to be w/o log_lik (as stored in seperate array only if user chooses)
      # names_tp <- names_tp_wo_log_lik
      # index_tp <- names_tp_wo_log_lik
      ##
      index_gq <- setdiff(bs_index_inc_tp_and_gq, bs_index_inc_tp)
      names_gq <- pars_names[index_gq]
      n_gq <- length(names_gq)
      ##
  })
  ##
  try({
      index_wo_log_lik <-  grep("^log_lik", pars_names, invert = TRUE)
      names_wo_log_lik <-  pars_names[index_wo_log_lik]
      length(index_wo_log_lik)
      length(pars_names)
      ##
      if (package == "BayesMVP") {
            ##
            ## setdiff(names_wo_nuisance, names_wo_log_lik)
            ##
            print(head(index_wo_log_lik))
            try({
              index_wo_nuisance <-  grep("^nuisance", pars_names, invert = TRUE)
              names_wo_nuisance <-  pars_names[index_wo_nuisance]
              print(head(index_wo_nuisance, 10))
            })
            ##
            try({ 
              names_wo_nuisance_and_log_lik <-  intersect(names_wo_log_lik, names_wo_nuisance)
              index_wo_nuisance_and_log_lik <-  intersect(index_wo_log_lik, index_wo_nuisance)
              n_params_wo_nuisance_and_log_lik <- length(index_wo_nuisance_and_log_lik)
            })
            ##
            try({ 
              print(paste("n_params_main = ", n_params_main))
              index_params_main <- index_wo_nuisance_and_log_lik[1]:(index_wo_nuisance_and_log_lik[1] + n_params_main - 1)
              print(head(index_params_main))
            })
      } else if (package == "MetaOrdDTA") { 
        
        
      }
  })
   # include_nuisance <- FALSE #### BOOKMARK / DEBUG
   ##
   ##
   if (package == "BayesMVP") {
      if (use_bridgestan == TRUE) {
         ## print(n_params_main)
         pars_indicies_to_track <- 1:n_par_inc_tp_and_gq
         n_params_full <- n_par_inc_tp_and_gq
         all_param_outs_trace_inc_log_lik <- BayesMVP:::fn_compute_param_constrain_from_trace_parallel( unc_params_trace_input_main = main_trace,
                                                                                                        unc_params_trace_input_nuisance = nuisance_trace,
                                                                                                        pars_indicies_to_track = pars_indicies_to_track,
                                                                                                        n_params_full = n_params_full,
                                                                                                        n_nuisance = n_nuisance,
                                                                                                        n_params_main = n_params_main,
                                                                                                        include_nuisance = include_nuisance,
                                                                                                        model_so_file = model_so_file,
                                                                                                        json_file_path = json_file_path)
      } else if (use_bridgestan == FALSE) { 
        
            # trace_params_main_as_posterior_array <- posterior::as_draws_array(trace_params_main)
            # 
            # str(trace_params_main_as_posterior_array)
            # 
            # stan_standalone_gq_outs <- stan_model_obj$generate_quantities(  fitted_params = trace_params_main_as_posterior_array,
            #                                                                 data = Stan_data_list,
            #                                                                 seed = 123,
            #                                                                 output_basename = NULL,
            #                                                                 parallel_chains = parallel::detectCores())
            # 
            # standalone_gq_draws <- stan_standalone_gq_outs$draws()
            # 
            # str(standalone_gq_draws)
            # n_gq
            # 
            # str(trace_params_main)
            # stan_model_obj <- cmdstan_model(outs_init_stan_model$Stan_model_file_path)
            # stan_model_obj
        
        
      }
     
   }
   ##
   ##
   offset <- 0
   ##
   ## trace_params_main <- array(dim = c(n_params_main, n_iter, n_chains))
   ## message(print(str(all_param_outs_trace_inc_log_lik)))
   ##
   trace_params_main <- NULL
   ##
   trace_tp <- NULL
   trace_tp_wo_log_lik <- NULL
   ##
   trace_gq <- NULL
   trace_log_lik <- NULL
   ##
   if (package == "BayesMVP") {
         if (use_bridgestan == TRUE) {
                   ##
                   ## Main params ("parameters" block):
                   ##
                   if (compute_main_params == TRUE) {
                     # kk <- 1
                     # all_param_outs_trace_inc_log_lik[[kk]][index_params_main - offset, 1:n_iter]
                     ## trace_params_main[1:n_params_main, 1:n_iter, kk] <- all_param_outs_trace_inc_log_lik[[kk]][index_params_main - offset, 1:n_iter]
                      for (kk in 1:n_chains) {
                        try({
                           trace_params_main[1:n_params_main, 1:n_iter, kk] <-   all_param_outs_trace_inc_log_lik[[kk]][index_params_main - offset, 1:n_iter]  ##  main_trace_list[[kk]][1:n_params_main, 1:n_iter]##  all_param_outs_trace[[kk]][index_params_main - offset, 1:n_iter] ## BOOKMARK
                        }, silent = TRUE)
                      }
                     
                   }
                   ##
                   ## Transformed parameters:
                   ##
                   try({
                     ##
                     message(print(paste("index_tp_wo_log_lik = ")))
                     message(print(head(index_tp_wo_log_lik)))
                     message(print(length(index_tp_wo_log_lik)))
                     ##
                     message(print(paste("index_tp = ")))
                     message(print(head(index_tp)))
                     message(print(length(index_tp)))
                     ##
                     n_tp_wo_log_lik <- length(index_tp_wo_log_lik)
                     
                     if (compute_transformed_parameters == TRUE) {
                       trace_tp <- array(dim = c(n_tp_wo_log_lik, n_iter, n_chains))
                       for (kk in 1:n_chains) {
                         trace_tp[1:n_tp_wo_log_lik, 1:n_iter, kk] <- all_param_outs_trace_inc_log_lik[[kk]][index_tp_wo_log_lik - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
                       }
                       # str(trace_tp)
                     }
                   })
                   ##
                   ## Generated quantities:
                   ##
                   try({
                     if (compute_generated_quantities == TRUE) {
                       trace_gq <- array(dim = c(n_gq, n_iter, n_chains))
                       for (kk in 1:n_chains) {
                         trace_gq[1:n_gq, 1:n_iter, kk] <- all_param_outs_trace_inc_log_lik[[kk]][index_gq - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
                       }
                     }
                     # str(trace_gq)
                   })
                   str(all_param_outs_trace_inc_log_lik)
                   length(index_gq)
                   str(all_param_outs_trace_inc_log_lik[[kk]])
                   ##
                   try({
                     if (save_log_lik_trace == TRUE) {
        
                           trace_log_lik <- array(dim = c(N, n_iter, n_chains))
        
                           if (Model_type == "Stan") {
        
                                  for (kk in 1:n_chains) {
                                      trace_log_lik[1:N,  1:n_iter, kk] <- all_param_outs_trace_inc_log_lik[[kk]][index_log_lik - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
                                  }
        
                           } else {  ## if built-in / manual model

                                 trace_log_lik <- trace_log_lik_mnl_models

                           }
        
                     }
                     # str(trace_log_lik)
                   })
         } else if (use_bridgestan == FALSE) { 
           
             #### BOOKMARK - Put code to compute traces for BayesMVP WITHOUT usign bridgestan here
           
         }
     
   } else if (package == "MetaOrdDTA") {
     
                 ##
                 ## Main params ("parameters" block):
                 ##
                 if (compute_main_params == TRUE) {
                     str(stan_draws_array)
                     trace_params_main_outs <- MetaOrdDTA:::filter_param_draws_string_batched( all_draws_array = stan_draws_array,
                                                                                      param_strings_vec = grouped_params$bs_names_main,
                                                                                      condition = "exact_match",
                                                                                      exclude_NA = exclude_NA,
                                                                                      print = FALSE)
                     
                     str(trace_params_main_outs)
                     trace_params_main <- trace_params_main_outs$param_draws_filtered_array
                     str(trace_params_main)
                 }
                 ##
                 ## Transformed parameters trace (inc. "log_lik"):
                 ##
                 grouped_params$bs_names_tp
                 ##
                 if (compute_transformed_parameters == TRUE) {
                   str(stan_draws_array)
                   trace_tp_outs <- MetaOrdDTA:::filter_param_draws_string_batched(  all_draws_array = stan_draws_array,
                                                                        param_strings_vec = grouped_params$bs_names_tp,
                                                                        condition = "exact_match",
                                                                        exclude_NA = exclude_NA,
                                                                        print = FALSE)
                   trace_tp_list <- trace_tp_outs$param_draws_filtered_list
                   trace_tp      <- trace_tp_outs$param_draws_filtered_array
                   str(trace_tp_list)
                 }
                 
                 ##
                 grouped_params$bs_names_tp_wo_log_lik <- grouped_params$bs_names_tp[grouped_params$bs_names_tp != "log_lik"]
                 ##
                 ##
                 ## Transformed parameters trace (excl. "log_lik"):
                 ##
                 if (compute_transformed_parameters == TRUE) {
                   str(stan_draws_array)
                   trace_tp_wo_log_lik_outs <- MetaOrdDTA:::filter_param_draws_string_batched( all_draws_array = stan_draws_array,
                                                                                  param_strings_vec = grouped_params$bs_names_tp,
                                                                                  condition = "exact_match",
                                                                                  exclude_NA = exclude_NA,
                                                                                  print = FALSE)
                   trace_tp_wo_log_lik_list <- trace_tp_wo_log_lik_outs$param_draws_filtered_list
                   trace_tp_wo_log_lik      <- trace_tp_wo_log_lik_outs$param_draws_filtered_array
                   # str(trace_tp_wo_log_lik)
                 }
                 ##
                 ## Generated quantities trace:
                 ##
                 if (compute_generated_quantities == TRUE) {
                   str(stan_draws_array)
                   trace_gq_outs <- MetaOrdDTA:::filter_param_draws_string_batched(   all_draws_array = stan_draws_array,
                                                                         param_strings_vec = grouped_params$bs_names_gq,
                                                                         condition = "exact_match",
                                                                         exclude_NA = exclude_NA,
                                                                         print = FALSE)
                   trace_gq_list <- trace_gq_outs$param_draws_filtered_list
                   trace_gq      <- trace_gq_outs$param_draws_filtered_array
                   str(trace_gq)
                 }
                 
                 trace_log_lik <- NULL
                 if (save_log_lik_trace == TRUE) {
                   str(stan_draws_array)
                   trace_log_lik_outs <- MetaOrdDTA:::filter_param_draws_string_batched(    all_draws_array = stan_draws_array,
                                                                               param_strings_vec = "log_lik",
                                                                               condition = "exact_match",
                                                                               exclude_NA = exclude_NA,
                                                                               print = TRUE)
                   trace_log_lik_list <- trace_log_lik_outs$param_draws_filtered_list
                   trace_log_lik      <- trace_log_lik_outs$param_draws_filtered_array
                   str(trace_log_lik)
                 }
                 ##
                 ## Re-format arrays so that dims  = c(param, n_iter, n_chains):
                 ##
                 re_format_3D_array_dims <- function(input_array) { 
                     
                           try({
                             str(input_array)
                             input_array_2 <- base::aperm(input_array, c(1,3,2))
                             str(input_array_2)
                             input_array_3 <- base::aperm(input_array_2, c(2,1,3))
                             str(input_array_3)
                             # trace_params_main <- input_array_3
                             
                             return(input_array_3)
                             
                           })
                   
                 }
                 ##
                 trace_params_main_3 <- NULL
                 trace_tp_3 <- NULL
                 trace_tp_wo_log_lik_3 <- NULL
                 trace_gq_3 <- NULL
                 trace_log_lik_3 <- NULL
                 ##
                 try({
                   trace_params_main_3 <- re_format_3D_array_dims(trace_params_main)
                   trace_params_main <- trace_params_main_3
                 })
                 ##
                 try({
                   trace_tp_3 <- re_format_3D_array_dims(trace_tp)
                   trace_tp <- trace_tp_3
                 })
                 ##
                 try({
                   trace_tp_wo_log_lik_3 <- re_format_3D_array_dims(trace_tp_wo_log_lik)
                   trace_tp_wo_log_lik <- trace_tp_wo_log_lik_3
                 })
                 ##
                 try({
                   trace_gq_3 <- re_format_3D_array_dims(trace_gq)
                   trace_gq <- trace_gq_3
                 })
                 ##
                 try({
                   trace_log_lik_3 <- re_format_3D_array_dims(trace_log_lik)
                   trace_log_lik <- trace_log_lik_3
                 })
   }
   ##
   ##
   n_cores <- parallel::detectCores()
   ##   
   ## --------- MAIN PARAMETERS / "PARAMETERS" BLOCK IN STAN  ----------------------------------------
   ##
     {
          Min_ESS_main <- NULL
          summary_tibble_main_params <- NULL
          ##
          Max_nested_rhat_main <- NULL
          ##
          # n_params_main <- length(names_main)
          # n_params_main
          ##
          # names_main
          # str(trace_params_main)
          ##
          # if (compute_main_params == TRUE) {
          
          if (package == "MetaOrdDTA") { 
                names_main <- bs_names_main
                temp_trace_to_use_for_summaries <- trace_params_main_3
                n_chains <- dim(temp_trace_main_use_for_summaries)[3]
                n_superchains <- MetaOrdDTA:::if_null_then_set_to(n_superchains, n_chains)
          } else if (package == "BayesMVP") { 
            temp_trace_to_use_for_summaries <- trace_params_main_3
                names_main <- head(names_wo_nuisance_and_log_lik, n_params_main)
                n_chains <- dim(temp_trace_main_use_for_summaries)[3]
          }
            
          # trace[i, 1:n_iter, kk] 
          # str(temp_trace_to_use_for_summaries)
          
          if ((package == "BayesMVP") || (use_BayesMVP_for_faster_summaries == TRUE)) {
          
                  # str(temp_trace_to_use_for_summaries)
            
                  summary_tibble_main_params <- BayesMVP:::generate_summary_tibble( n_threads = n_cores,
                                                                                    trace = temp_trace_to_use_for_summaries,
                                                                                    param_names = names_main,
                                                                                    n_to_compute = n_params_main,
                                                                                    compute_nested_rhat = compute_nested_rhat,
                                                                                    n_chains = n_chains, 
                                                                                    n_superchains = n_superchains)
                  
                  Min_ESS_main <-  min(na.rm = TRUE, summary_tibble_main_params$n_eff[1:n_params_main])
                  Max_rhat_main <- max(na.rm = TRUE, summary_tibble_main_params$Rhat[1:n_params_main])
                  ##
                  if (compute_nested_rhat == TRUE) {
                    Max_nested_rhat_main <- max(na.rm = TRUE, summary_tibble_main_params$n_Rhat[1:n_params_main])
                  }
              
          } else if ((package == "MetaOrdDTA") || (use_BayesMVP_for_faster_summaries == FALSE)) { 

                  summary_tibble_main_params <- MetaOrdDTA:::R_fn_compute_cmdstanr_summary_tibble( stan_fitted_model = stan_fitted_model, 
                                                                                         variables_vec =  grouped_params$bs_names_main)
                  Min_ESS_main < - round(min(na.rm = TRUE, summary_tibble_main_params$ess_bulk[1:n_params_main]))
                  Max_rhat_main <- max(na.rm = TRUE, summary_tibble_main_params$rhat[1:n_params_main])
            
          }
          
                        
                        

        # }
        
        Min_ESS <- Min_ESS_main ## bookmark
   }
   ##
   ## --------- TRANSFORMED PARAMETERS -------------------------------------------------------------- 
   ##
   {
     summary_tibble_transformed_parameters <- NULL
     if ((compute_transformed_parameters == TRUE) && (n_tp > 0)) {
       
       ## Remove log-lik trace from tp trace array (since we store it in seperate array called "trace_log_lik" if user chooses to store it)
       ## trace_tp <- R_fn_remove_log_lik_from_array(trace_tp)
       #     length(names_tp_wo_log_lik)
       # n_tp_wo_log_lik
       if ((package == "BayesMVP") || (use_BayesMVP_for_faster_summaries == TRUE)) {
         
                   ## str(trace_tp_wo_log_lik_3)
                   
                   if (package == "MetaOrdDTA") { 
                     temp_trace_to_use_for_summaries <- trace_tp_wo_log_lik_3
                     n_chains <- dim(temp_trace_main_use_for_summaries)[3]
                     n_superchains <- MetaOrdDTA:::if_null_then_set_to(n_superchains, n_chains)
                   } else if (package == "BayesMVP") { 
                     temp_trace_to_use_for_summaries <- trace_tp_wo_log_lik_3
                     n_chains <- dim(temp_trace_main_use_for_summaries)[3]
                   }
         
                   summary_tibble_transformed_parameters <- BayesMVP:::generate_summary_tibble(  n_threads = n_cores,
                                                                                                 trace = temp_trace_to_use_for_summaries,
                                                                                                 param_names = names_tp_wo_log_lik,
                                                                                                 n_to_compute = n_tp_wo_log_lik,
                                                                                                 compute_nested_rhat = compute_nested_rhat,
                                                                                                 n_chains = n_chains, 
                                                                                                 n_superchains = n_superchains)
       
       } else if ((package == "MetaOrdDTA") || (use_BayesMVP_for_faster_summaries == FALSE)) { 

                   summary_tibble_transformed_parameters <- MetaOrdDTA:::R_fn_compute_cmdstanr_summary_tibble( stan_fitted_model = stan_fitted_model, 
                                                                                                  variables_vec =  names_tp_wo_log_lik)
         
       }
       
     }
   }
   ##
   ## --------- GENERATED QUANTITIES ---------------------------------------------------------------- 
   ##
   {
        summary_tibble_generated_quantities <- NULL
        ##
        # str(trace_gq)
        # length(names_gq)
        # n_gq
        ##
        if  ((compute_generated_quantities == TRUE) && (n_gq > 0))  { 
                      
          if ((package == "BayesMVP") || (use_BayesMVP_for_faster_summaries == TRUE)) {
            
                          if (package == "MetaOrdDTA") { 
                              temp_trace_to_use_for_summaries <- trace_gq_3
                              n_chains <- dim(temp_trace_main_use_for_summaries)[3]
                              n_superchains <- MetaOrdDTA:::if_null_then_set_to(n_superchains, n_chains)
                          } else if (package == "BayesMVP") { 
                              temp_trace_to_use_for_summaries <- trace_gq_3
                              n_chains <- dim(temp_trace_main_use_for_summaries)[3]
                          }
                
                          summary_tibble_generated_quantities <- BayesMVP:::generate_summary_tibble(  n_threads = n_cores,
                                                                                                      trace = temp_trace_to_use_for_summaries,
                                                                                                      param_names = names_gq,
                                                                                                      n_to_compute = n_gq,
                                                                                                      compute_nested_rhat = compute_nested_rhat,
                                                                                                      n_chains = n_chains, 
                                                                                                      n_superchains = n_superchains)
                  
              } else if ((package == "MetaOrdDTA") || (use_BayesMVP_for_faster_summaries == FALSE)) {  
    
                          summary_tibble_generated_quantities <- MetaOrdDTA:::R_fn_compute_cmdstanr_summary_tibble( stan_fitted_model = stan_fitted_model, 
                                                                                                                    variables_vec =  names_gq)
                
              }
                            
        }
   }
   ##
   ## -----------------------------  DF / tibble creation (for "posterior" and "bayesplot" R packages)  --------------------------------------------
   ##
   {
         trace_params_main_tibble <- NULL
         trace_transformed_params_tibble <- NULL
         trace_generated_quantities_tibble <- NULL  
         # ##
         trace_params_main_reshaped <- NULL
         trace_tp_reshaped <- NULL
         trace_gq_reshaped <- NULL
         ##
         # 
         # ## re-shape data (so can easily use with "posterior" & "bayesplot" R packages)
         # trace_params_main_reshaped <-  base::aperm(trace_params_main_3, c(2, 3, 1))
         # str(trace_log_lik)
         # str(trace_tp)
         # str(trace_tp_wo_log_lik)
         # if (compute_transformed_parameters == TRUE)   trace_tp_reshaped <- base::aperm(trace_tp_wo_log_lik_3, c(2, 3, 1))
         # if (compute_generated_quantities == TRUE)     trace_gq_reshaped <- base::aperm(trace_gq_3, c(2, 3, 1))
         # 
         # ## add the parameter names to the array
         # dimnames(trace_params_main_reshaped) <- list(iterations = 1:n_iter, 
         #                                             chains = 1:n_chains, 
         #                                             parameters = names_main)
         # ##
         # str(trace_tp_wo_log_lik_3)
         # str(trace_tp_wo_log_lik_2)
         # str(trace_tp_wo_log_lik)
         # str(trace_tp_reshaped)
         # length(names_tp)
         # length(names_tp_wo_log_lik)
         # 
         # reshaped_traces_list <- list()
         # reshaped_traces_list$trace_tp <- base::aperm(trace_, c(2, 3, 1))
         # 
         # trace_tp_reshaped <- base::aperm(trace_tp_wo_log_lik_3, c(2, 3, 1))
         # 
         # if (compute_transformed_parameters == TRUE)  dimnames(trace_tp_reshaped) <- list(iterations = 1:n_iter, 
         #                                                                                   chains = 1:n_chains,
         #                                                                                   parameters = names_tp_wo_log_lik)
         # ##
         # if (compute_generated_quantities == TRUE) dimnames(trace_gq_reshaped) <-    list(iterations = 1:n_iter, 
         #                                                                                 chains = 1:n_chains, 
         #                                                                                 parameters = names_gq)
     
         #### -----------------------------  DF / tibble creation (for "posterior" and "bayesplot" R packages)  ------
         try({  
           {
               ##
               ## convert from array -> to df/tibble format 
               ##
               str(trace_params_main) ## make sure its "dims = (iter, chains, params)"
               str(trace_tp_wo_log_lik) ## make sure its "dims = (iter, chains, params)"
               str(trace_gq) ## make sure its "dims = (iter, chains, params)"
            
               if (save_trace_tibbles == TRUE) {
                   trace_params_main_tibble <- dplyr::tibble(posterior::as_draws_df(trace_params_main))
                   if (compute_transformed_parameters == TRUE)  trace_transformed_params_tibble   <- dplyr::tibble(posterior::as_draws_df(trace_tp_wo_log_lik))
                   if (compute_generated_quantities == TRUE)    trace_generated_quantities_tibble <- dplyr::tibble(posterior::as_draws_df(trace_gq))
               }
           }
         })
         # ####  ----- Make OVERALL (full) draws array   ------------------------------------------------------------
         # {
         #     dims_main <- dim(trace_params_main)[3] ; dims_main
         #     dims_tp_wo_log_lik   <- dim(trace_tp_wo_log_lik)[3] ; dims_tp_wo_log_lik
         #     dims_gq   <- dim(trace_gq)[3] ; dims_gq
         #     ##
         #     dim_total <- dims_main
         #     if (compute_transformed_parameters == TRUE) dim_total <- dim_total + dims_tp_wo_log_lik
         #     if (compute_generated_quantities == TRUE)   dim_total <- dim_total + dims_gq
         #     dim_total
         #     ##
         #     names_total <- names_main
         #     if (compute_transformed_parameters == TRUE)  names_total <- c(names_total, names_tp_wo_log_lik)
         #     if (compute_generated_quantities == TRUE)    names_total <- c(names_total, names_gq)
         #     ##
         #     length(names_main)
         #     length(names_tp_wo_log_lik)
         #     length(names_gq)
         #     length(names_total)
         #     ##
         #     length(names_main) + 
         #     length(names_tp_wo_log_lik) + 
         #     length(names_gq) + 
         #     length(names_log_lik)
         # }
         
         
   }
   ##
   # draws_array <- array(NA, dim = c(n_iter, n_chains, dim_total))
   # draws_array[1:n_iter, 1:n_chains, 1:n_params_main] <- trace_params_main_reshaped # main first 
   # 
   # if (compute_transformed_parameters == TRUE)   {
   #   draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_tp)] <- trace_tp_reshaped
   # }
   # if ((compute_generated_quantities == TRUE) &&  (compute_transformed_parameters == TRUE))   {
   #   draws_array[1:n_iter, 1:n_chains, (n_params_main + n_tp + 1):(n_params_main + n_tp + n_gq)] <- trace_gq_reshaped
   # }
   # 
   # if ((compute_generated_quantities == TRUE) &&  (compute_transformed_parameters == FALSE))   {
   #   draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_gq)] <- trace_gq_reshaped
   # }
   # if ((compute_generated_quantities == FALSE) &&  (compute_transformed_parameters == TRUE))   {
   #   draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_tp)] <- trace_tp_reshaped
   # }
   ##
#    dimnames(draws_array) <- list( iterations = 1:n_iter, 
#                                   chains = 1:n_chains, 
#                                   parameters = names_total)
   
   
     try({
       print(tictoc::toc(log = TRUE))
       log.txt <- tictoc::tic.log(format = TRUE)
       tictoc::tic.clearlog()
       time_summaries <- unlist(log.txt)
       ##
       extract_numeric_string <-  stringr::str_extract(time_summaries, "\\d+\\.\\d+")   
       time_summaries <- as.numeric(extract_numeric_string)
     })
     ##
  {
          time_total <- time_summaries + time_burnin + time_sampling
          ##
          try({ 
              ESS_per_sec_samp <- Min_ESS_main / time_sampling
              ESS_per_sec_total <- Min_ESS_main / time_total
          })
          ##
          if (package == "BayesMVP") {
               EHMC_args_as_Rcpp_List <- model_results$init_burnin_object$EHMC_args_as_Rcpp_List
          }
          ##
          try({
            message(((paste("Max R-hat (parameters block, main only) = ", round(Max_rhat_main, 4)))))
          }, silent = TRUE)
          try({
            message(((paste("Max R-hat (parameters block, main only) = ", round(Max_nested_rhat_main, 4)))))
          }, silent = TRUE)
          try({
            message(((paste("Min ESS (parameters block, main only) = ", round(Min_ESS_main, 0)))))
            message(((paste("Min ESS / sec [samp.] (parameters block, main only) = ", signif(ESS_per_sec_samp, 3)))))
            message(((paste("Min ESS / sec [overall] (parameters block, main only) = ", signif(ESS_per_sec_total, 3)))))
          }, silent = TRUE)
          ##
          ## ESS / grad:
          ##
          if (package == "BayesMVP") {
            
                    try({
                      L_main_during_sampling <- (EHMC_args_as_Rcpp_List$tau_main / EHMC_args_as_Rcpp_List$eps_main)
                      n_grad_evals_sampling_main <- L_main_during_sampling * n_iter * n_chains_sampling
                      Min_ess_per_grad_main_samp <-  Min_ESS_main / n_grad_evals_sampling_main
                    }, silent = TRUE)
                    ##
                    try({
                      L_us_during_sampling <- (EHMC_args_as_Rcpp_List$tau_us / EHMC_args_as_Rcpp_List$eps_us)
                      n_grad_evals_sampling_us <-  L_us_during_sampling  * n_iter * n_chains_sampling
                      Min_ess_per_grad_us_samp <-  Min_ESS_main / n_grad_evals_sampling_us
                    }, silent = TRUE)
                    ##
                    try({ 
                      if (partitioned_HMC == TRUE) { ## i.e. if nuisance are sampled seperately
                        weight_nuisance_grad <- 0.3333333
                        weight_main_grad <- 0.6666667 ## main grad takes ~ 2x as long to compute as nuisance grad
                        Min_ess_per_grad_samp_weighted <- (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) / (weight_nuisance_grad + weight_main_grad)
                      } else if (partitioned_HMC == FALSE) {  # if not partitioned, grad isnt seperate and "main" grad = "all" grad so use weight_main_grad only !!
                        weight_nuisance_grad <- 0.00
                        weight_main_grad <- 1.00
                        Min_ess_per_grad_samp_weighted <- (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) / (weight_nuisance_grad + weight_main_grad)
                      }
                      message(((paste("Min ESS / grad [samp., weighted] (parameters block, main only) = ", signif(1000 *  Min_ess_per_grad_samp_weighted, 3)))))
                    }, silent = TRUE)
                    ##
                    try({ 
                      grad_evals_per_sec <- ESS_per_sec_samp / Min_ess_per_grad_samp_weighted
                      message(((paste("Grad evals / sec [samp.] (parameters block, main only) = ", signif(grad_evals_per_sec, 3)))))
                    }, silent = TRUE)
            
          } else if (package == "MetaOrdDTA") { 
            
                    L_main_during_burnin <- NULL
                    ##
                    Stan_HMC_info_list <- MetaOrdDTA:::R_fn_compute_L_and_epsilon_from_Stan_model( Min_ESS = Min_ESS, 
                                                                                                   stan_fitted_model = stan_fitted_model)
                    ##
                    mean_epsilon_sampling <- Stan_HMC_info_list$out_list_sampling$mean_epsilon_sampling ; mean_epsilon_sampling
                    eps <- mean_epsilon_sampling
                    eps_main <- mean_epsilon_sampling
                    mean_L_sampling <- Stan_HMC_info_list$out_list_sampling$mean_L_sampling ; mean_L_sampling
                    ##
                    mean_tau_sampling <-  Stan_HMC_info_list$out_list_sampling$mean_tau_sampling ; mean_tau_sampling
                    tau <- mean_tau_sampling
                    tau_main <- mean_tau_sampling
                    ##
                    L_mean_per_chain_sampling <- Stan_HMC_info_list$out_list_sampling$L_mean_per_chain_sampling ; L_mean_per_chain_sampling
                    L_max_per_chain_sampling <- Stan_HMC_info_list$out_list_sampling$L_max_per_chain_sampling ; L_max_per_chain_sampling
                    ##
                    max_of_mean_Ls_per_chain_sampling <- Stan_HMC_info_list$out_list_sampling$max_of_mean_Ls_per_chain_sampling ; max_of_mean_Ls_per_chain_sampling
                    ##
                    stan_total_gradient_evals_sampling <- Stan_HMC_info_list$out_list_sampling$stan_total_gradient_evals_sampling
                    stan_min_ess_per_grad_eval_sampling <- Stan_HMC_info_list$out_list_sampling$stan_min_ess_per_grad_eval_sampling
                    stan_grad_evals_per_sec_sampling <- Stan_HMC_info_list$out_list_sampling$stan_grad_evals_per_sec_sampling
                    stan_grad_evals_per_sec_sampling_div_1000 <- Stan_HMC_info_list$out_list_sampling$stan_grad_evals_per_sec_sampling_div_1000
                    ##
                    L_main_during_sampling <- mean_L_sampling
                    ##
                    try({
                      n_grad_evals_sampling_main <-  L_main_during_sampling * n_iter * n_chains_sampling
                      Min_ess_per_grad_main_samp <-  Min_ESS_main / n_grad_evals_sampling_main
                    }, silent = TRUE)
                    ##
                    message(((paste("Min ESS / grad [samp., weighted] (parameters block, main only) = ", signif(1000 *  Min_ess_per_grad_main_samp, 3)))))
                    ##
                    try({ 
                      grad_evals_per_sec <- ESS_per_sec_samp / Min_ess_per_grad_main_samp
                      message(((paste("Grad evals / sec [samp.] (parameters block, main only) = ", signif(grad_evals_per_sec, 3)))))
                    }, silent = TRUE)
          }
          ##
          ##
          try({ 
                sampling_time_to_Min_ESS <- time_sampling
                sampling_time_to_100_ESS <- (100 / Min_ESS_main) * sampling_time_to_Min_ESS
                sampling_time_to_1000_ESS <- (1000 / Min_ESS_main) * sampling_time_to_Min_ESS
                sampling_time_to_10000_ESS <- (10000 / Min_ESS_main) * sampling_time_to_Min_ESS
                ##
                ## w/o summary time
                total_time_to_100_ESS_wo_summaries <-   time_burnin + sampling_time_to_100_ESS
                total_time_to_1000_ESS_wo_summaries <-  time_burnin + sampling_time_to_1000_ESS
                total_time_to_10000_ESS_wo_summaries <- time_burnin + sampling_time_to_10000_ESS
                ##
                ## w/ summary time (note: assuming that time_summaries scales linearly w/ 
                # the min ESS required (i.e. n_iter and/or n_chains) Might be over-estimate)
                summary_time_to_Min_ESS <- time_summaries
                summary_time_to_100_ESS <- (100 / Min_ESS_main) * summary_time_to_Min_ESS
                summary_time_to_1000_ESS <- (1000 / Min_ESS_main) * summary_time_to_Min_ESS
                summary_time_to_10000_ESS <- (10000 / Min_ESS_main) * summary_time_to_Min_ESS
                ##
                total_time_to_100_ESS_with_summaries <-   time_burnin + sampling_time_to_100_ESS   + summary_time_to_100_ESS
                total_time_to_1000_ESS_with_summaries <-  time_burnin + sampling_time_to_1000_ESS  + summary_time_to_1000_ESS
                total_time_to_10000_ESS_with_summaries <- time_burnin + sampling_time_to_10000_ESS + summary_time_to_10000_ESS
          })
          ##
          try({ 
              n_divs <- sum(n_divs_per_chain)
              total_transitions <- n_chains * n_iter
              pct_divs <- 100*(n_divs / total_transitions)
          })
          ##
          if (n_divs > 0) { 
            warning("divergences detected!")
            try({ 
              message(print(paste("Number of divergences = ", n_divs)))
              message(print(paste("% divergences = ", pct_divs)))
            })
          } else { 
            message("No divergences detected!")
          }
          
    }
    ##
    ## ---------  Make R lists to output  --------------------------------------------------------------------
    ##
    {
      summary_tibbles   <- NULL
      traces_as_arrays  <- NULL
      traces_as_tibbles <- NULL
      ##
      ## ---- List to store summary tibbles / DF's:
      ##
      summary_tibbles <- list( summary_tibble_main_params = summary_tibble_main_params,
                               summary_tibble_transformed_parameters = summary_tibble_transformed_parameters,
                               summary_tibble_generated_quantities = summary_tibble_generated_quantities)
      ##
      ## ---- List to store traces (as 3D arrays):
      ##
      ##
      traces_as_arrays <- list( 
                               # # draws_array = draws_array,
                               # trace_params_main = trace_params_main_reshaped,
                               # trace_transformed_params = trace_tp_reshaped,
                               # trace_generated_quantities = trace_gq_reshaped
                               trace_params_main = trace_params_main,
                               trace_tp = trace_tp,
                               trace_tp_wo_log_lik = trace_tp_wo_log_lik,
                               trace_gq = trace_gq,
                               trace_log_lik = trace_log_lik)
      ##
      ## ---- List to store traces (as tibbles/DF's):
      ##
      traces_as_tibbles <- list( trace_params_main_tibble = trace_params_main_tibble,
                                 trace_transformed_params_tibble = trace_transformed_params_tibble,
                                 trace_generated_quantities_tibble = trace_generated_quantities_tibble)
    }
    if (package == "BayesMVP") {
      
              ##
              ##  ---- List to store HMC info ("HMC_info" list):
              ##
              HMC_info <- list(   tau_main = EHMC_args_as_Rcpp_List$tau_main,
                                  eps_main = EHMC_args_as_Rcpp_List$eps_main,
                                  L_main_during_burnin = L_main_during_burnin,
                                  L_main_during_sampling = L_main_during_sampling,
                                  ##
                                  tau_us = EHMC_args_as_Rcpp_List$tau_us,
                                  eps_us = EHMC_args_as_Rcpp_List$eps_us,
                                  L_us_during_burnin = L_us_during_burnin,
                                  L_us_during_sampling = L_us_during_sampling,
                                  ##
                                  n_chains_sampling = n_chains_sampling,
                                  n_chains_burnin = n_chains_burnin,
                                  n_superchains = n_superchains,
                                  n_iter = n_iter,
                                  n_burnin = n_burnin,
                                  ##
                                  adapt_delta = adapt_delta,
                                  LR_main = LR_main,
                                  LR_us = LR_us,
                                  ##
                                  metric_type_main = metric_type_main,
                                  metric_shape_main = metric_shape_main,
                                  metric_type_nuisance = metric_type_nuisance,
                                  metric_shape_nuisance = metric_shape_nuisance,
                                  ##
                                  diffusion_HMC = diffusion_HMC,
                                  partitioned_HMC = partitioned_HMC,
                                  ##
                                  interval_width_main = interval_width_main,
                                  interval_width_nuisance = interval_width_nuisance,
                                  ##
                                  force_autodiff = force_autodiff,
                                  force_PartialLog = force_PartialLog,
                                  multi_attempts = multi_attempts)
              ##
              ## ---- List to store efficiency information ("efficiency_info" list)
              ##
              efficiency_info <- list(              Max_rhat_main = Max_rhat_main,
                                                    Max_nested_rhat_main = Max_nested_rhat_main,
                                                    ##
                                                    Min_ESS_main = Min_ESS_main, 
                                                    ESS_per_sec_samp = ESS_per_sec_samp, 
                                                    ESS_per_sec_total = ESS_per_sec_total,
                                                    ##
                                                    time_burnin = time_burnin, 
                                                    time_sampling = time_sampling, 
                                                    time_summaries = time_summaries,
                                                    time_total_wo_summaries = time_total_wo_summaries, 
                                                    time_total = time_total,
                                                    ##
                                                    Min_ess_per_grad_samp_weighted = Min_ess_per_grad_samp_weighted,
                                                    grad_evals_per_sec = grad_evals_per_sec,
                                                    ##
                                                    sampling_time_to_Min_ESS = sampling_time_to_Min_ESS,
                                                    sampling_time_to_100_ESS = sampling_time_to_100_ESS,
                                                    sampling_time_to_1000_ESS = sampling_time_to_1000_ESS,
                                                    sampling_time_to_10000_ESS = sampling_time_to_10000_ESS,
                                                    ##
                                                    total_time_to_100_ESS_wo_summaries = total_time_to_100_ESS_wo_summaries,
                                                    total_time_to_1000_ESS_wo_summaries = total_time_to_1000_ESS_wo_summaries,
                                                    total_time_to_10000_ESS_wo_summaries = total_time_to_10000_ESS_wo_summaries,
                                                    ##
                                                    total_time_to_100_ESS_with_summaries = total_time_to_100_ESS_with_summaries,
                                                    total_time_to_1000_ESS_with_summaries = total_time_to_1000_ESS_with_summaries,
                                                    total_time_to_10000_ESS_with_summaries = total_time_to_10000_ESS_with_summaries)
              
    } else if (package == "MetaOrdDTA") {
              
              ##
              ## List to store HMC info:
              ##
              HMC_info <- list(   Stan_HMC_info_list = Stan_HMC_info_list,
                                  ##
                                  eps = eps,
                                  ##
                                  L_main_during_burnin = L_main_during_burnin,
                                  L_main_during_sampling = L_main_during_sampling,
                                  ##
                                  tau = tau,
                                  ##
                                  n_chains_sampling = n_chains_sampling,
                                  n_chains_burnin = n_chains_burnin,
                                  n_chains = n_chains,
                                  ##
                                  n_superchains = n_superchains,
                                  ##
                                  n_iter = n_iter,
                                  n_burnin = n_burnin,
                                  ##
                                  adapt_delta = adapt_delta,
                                  max_treedepth = max_treedepth,
                                  metric_shape = metric_shape)
              ##
              ## list to store efficiency information:
              ##
              efficiency_info <- list(              n_iter = n_iter,
                                                    ##
                                                    Max_rhat_main = Max_rhat_main,
                                                    Max_nested_rhat_main = Max_nested_rhat_main,
                                                    ##
                                                    Min_ESS_main = Min_ESS_main, 
                                                    ESS_per_sec_samp = ESS_per_sec_samp, 
                                                    ESS_per_sec_total = ESS_per_sec_total,
                                                    ##
                                                    time_burnin = time_burnin, 
                                                    time_sampling = time_sampling, 
                                                    time_summaries = time_summaries,
                                                    time_total_wo_summaries = time_total_wo_summaries, 
                                                    time_total = time_total,
                                                    ##
                                                    Min_ess_per_grad_main_samp = Min_ess_per_grad_main_samp,
                                                    grad_evals_per_sec = grad_evals_per_sec,
                                                    ##
                                                    sampling_time_to_Min_ESS = sampling_time_to_Min_ESS,
                                                    sampling_time_to_100_ESS = sampling_time_to_100_ESS,
                                                    sampling_time_to_1000_ESS = sampling_time_to_1000_ESS,
                                                    sampling_time_to_10000_ESS = sampling_time_to_10000_ESS,
                                                    ##
                                                    total_time_to_100_ESS_wo_summaries = total_time_to_100_ESS_wo_summaries,
                                                    total_time_to_1000_ESS_wo_summaries = total_time_to_1000_ESS_wo_summaries,
                                                    total_time_to_10000_ESS_wo_summaries = total_time_to_10000_ESS_wo_summaries,
                                                    ##
                                                    total_time_to_100_ESS_with_summaries = total_time_to_100_ESS_with_summaries,
                                                    total_time_to_1000_ESS_with_summaries = total_time_to_1000_ESS_with_summaries,
                                                    total_time_to_10000_ESS_with_summaries = total_time_to_10000_ESS_with_summaries
                                                    )
    }
   
   


    
    ##
    ## ------------ final lists to output:
    ##
     if (package == "BayesMVP") {
          if (save_log_lik_trace == FALSE)   trace_log_lik <- NULL
          if (save_nuisance_trace == FALSE)  nuisance_trace <- NULL
          ##
          if (save_nuisance_trace == TRUE) {
            ## do nothing
          } else {
            nuisance_trace <- NULL
          }
      ##
     } else if (package == "MetaOrdDTA") { 
       nuisance_trace <- NULL
     }
     ##
     ## ---- "HMC_diagnostic_info" list (nested inside the "HMC_info" list made above):
     ##
     n_divs_per_chain      <- NULL
     n_max_trees_per_chain <- NULL
     ebfmi_per_chain       <- NULL
     divergences           <- NULL
     HMC_diagnostic_info   <- NULL
     ##
     try({ 

        {
            n_divs_per_chain      <- stan_diagnostic_summary$num_divergent ## For MetaOrdDTA only
            n_max_trees_per_chain <- stan_diagnostic_summary$num_max_treedepth ## For MetaOrdDTA only
            ebfmi_per_chain       <- stan_diagnostic_summary$ebfmi ## For MetaOrdDTA only
            ##
            divergences <- list( n_divs = n_divs,
                                 pct_divs = pct_divs,
                                 n_divs_per_chain = n_divs_per_chain)
        }
        ##
        HMC_diagnostic_info <-  list(divergences = divergences,
                                     n_max_trees_per_chain = n_max_trees_per_chain,
                                     ebfmi_per_chain = ebfmi_per_chain)
        ##
        ## Add it to the "HMC_info" list:
        ##
        HMC_info <- HMC_info$HMC_diagnostic_info
     })
    ##
    ## ---- "summaries" list:
    ##
    summaries <- list(summary_tibbles = summary_tibbles)
    ##
    ## ---- "traces" list:
    ##
    # traces_as_arrays <- list( 
    #   trace_params_main = trace_params_main,
    #   trace_tp = trace_tp,
    #   trace_tp_wo_log_lik = trace_tp_wo_log_lik,
    #   trace_gq = trace_gq,
    #   trace_log_lik = trace_log_lik)
    ##
    traces <- list(traces_as_arrays = traces_as_arrays,
                   traces_as_tibbles = traces_as_tibbles,
                   # trace_log_lik = trace_log_lik,
                   nuisance_trace = nuisance_trace)
    ##
    ## ---------- list to output: --------
    ##
    output_list <- list(  HMC_info = HMC_info,
                          efficiency_info = efficiency_info,
                          summaries = summaries,
                          traces = traces)
    
    # ### store output (optional)
    # if (store_outputs == TRUE) {
    #   if (is.null(store_outputs_dir)) { 
    #     store_outputs_dir <- getwd() # store in users wd if no dir specified
    #   }
    #   saveRDS(object = output_list, file = paste("BayesMVP_seed_", seed, 
    #                                              "Model_type_", Model_type, 
    #                                              "N_", N, 
    #                                              ))
    #   
    # }
                                          
    
  ### output 
  return(output_list)
  
  
  
  
}














