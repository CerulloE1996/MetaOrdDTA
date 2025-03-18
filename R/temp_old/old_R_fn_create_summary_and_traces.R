





 

## Helper R function to remove "log_lik" parameter from trace array 
#' R_fn_remove_log_lik_from_array
#' @keywords internal
R_fn_remove_log_lik_from_array <- function(arr) {
  
  format(object.size(arr), units = "MB")
  
  # Find which parameters don't contain "log_lik"
  keep_params <- !grepl("log_lik", dimnames(arr)[[3]])
  
  # Create new array without the log_lik parameters
  arr_filtered <- arr[, , keep_params, drop = FALSE]
  
  return(arr_filtered)
  
}





#### ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' create_summary_and_traces
#' @keywords internal
create_summary_and_traces <- function(    model_results,
                                          init_object,
                                          n_nuisance, 
                                          compute_main_params = TRUE, # excludes nuisance params. and log-lik 
                                          compute_transformed_parameters = TRUE,
                                          compute_generated_quantities = FALSE,
                                          save_log_lik_trace = TRUE, 
                                          save_nuisance_trace = FALSE,
                                          compute_nested_rhat = NULL,
                                          n_superchains = NULL,
                                          save_trace_tibbles = FALSE
) {
  

  
  ## Extract essential model info from "init_object" object
  Model_type <- init_object$Model_type
  y <- init_object$y
  sample_nuisance <- init_object$sample_nuisance
  #### N <- nrow(y) 
  
  ## Extract traces 
  main_trace <- model_results$result[[1]]
  div_trace <- model_results$result[[2]]
  nuisance_trace <- model_results$result[[3]]
  
  if (save_log_lik_trace == TRUE) {
      if (Model_type != "Stan") {
        log_lik_trace_mnl_models <-  model_results$result[[6]] 
      } else { 
        log_lik_trace_mnl_models <- NULL
      }
  } else { 
        log_lik_trace_mnl_models <- NULL
  }
  
  
  time_burnin <- model_results$time_burnin
  time_sampling <- model_results$time_sampling
  time_total_wo_summaries <- time_burnin + time_sampling
  
  #### Other MCMC / HMC info
  n_chains_burnin <- model_results$n_chains_burnin
  n_burnin <- model_results$n_burnin
  
  LR_main <- model_results$LR_main
  LR_us <- model_results$LR_us
  adapt_delta <- model_results$adapt_delta
  
  metric_type_main <- model_results$metric_type_main
  metric_shape_main <- model_results$metric_shape_main
  metric_type_nuisance <- model_results$metric_type_nuisance
  metric_shape_nuisance <- model_results$metric_shape_nuisance
  
  diffusion_HMC <- model_results$diffusion_HMC
  partitioned_HMC <- model_results$partitioned_HMC
  
  n_superchains <- model_results$n_superchains
  interval_width_main <- model_results$interval_width_main
  interval_width_nuisance <- model_results$interval_width_nuisance
  force_autodiff <- model_results$force_autodiff
  force_PartialLog <- model_results$force_PartialLog
  multi_attempts <- model_results$multi_attempts
  
  {

  tictoc::tic()
    
  if (sample_nuisance == TRUE) {
      n_nuisance <- n_nuisance
      n_nuisance_tracked <- dim(nuisance_trace[[1]])[1]
  } else { 
      n_nuisance <- 0
      n_nuisance_tracked <- 0
  }
  
  n_divs <- sum(unlist(div_trace))
  pct_divs <- 100 * n_divs / length(unlist(div_trace))
  
  n_chains <- length(main_trace)
  
  if (is.null(compute_nested_rhat)) {
    if (n_chains > 15) {
      compute_nested_rhat <- TRUE
    }
  }
  
  n_chains_burnin <- init_object$n_chains_burnin
  n_iter <-   dim(main_trace[[1]])[2]
  
  n_params_main <- dim(main_trace[[1]])[1]
  
  print(paste("n_params_main = ", n_params_main))
  
  if (is.null(n_superchains)) {
    n_superchains <- round(n_chains / n_chains_burnin)
  }
  
  # Stan_model_file_path <- (file.path(pkg_dir, "inst/stan_models/PO_LC_MVP_bin.stan"))  ### TEMP
  Stan_model_file_path <- init_object$Stan_model_file_path
  
  if (Model_type == "Stan") {
    json_file_path <- init_object$json_file_path
    model_so_file <- init_object$model_so_file
  } else { 
    json_file_path <- init_object$dummy_json_file_path
    model_so_file <- init_object$dummy_model_so_file
  }
  
  
  
  Sys.setenv(STAN_THREADS = "true")
  bs_model <- StanModel$new(Stan_model_file_path, data = json_file_path, 1234) # creates .so fil
  bs_names  <-  (bs_model$param_names())
  
  comment(print(paste("bs_names - head = ", head(bs_names))))
  comment(print(paste("bs_names - tail = ", tail(bs_names))))
  
  bs_names_inc_tp <-  (bs_model$param_names(include_tp = TRUE))
  bs_names_inc_tp_and_gq <-  (bs_model$param_names(include_tp = TRUE, include_gq = TRUE))
  
  
  
  ## pars_names <- bs_names_inc_tp_and_gq
  pars_names <- init_object$param_names
  
  if (model_obj$init_object$param_names[1] == "lp__") { 
    pars_names <- pars_names[-1]
  }
  
  
  comment(print(paste("pars_names - head = ", head(pars_names))))
  comment(print(paste("pars_names - tail = ", tail(pars_names))))
  
  index_lp  <- grep("^lp__", pars_names, invert = FALSE)
  
  # if (index_lp == 1) {  # if Stan model generates __lp variable (not all models will)
  #     pars_names <- pars_names[-c(1)]
  # } else { 
  #     pars_names <- pars_names
  # }
  # 
  

 #  pars_names <- bs_model$param_names(  include_tp = TRUE, include_gq = TRUE)
 #   pars_names <- init_object$init_vals_object$param_names
  n_par_inc_tp_and_gq <- length(pars_names) 
  try({ 
  names_nuisance_tracked <- head(pars_names, n_nuisance_tracked)
  })
  
  index_log_lik  <- grep("^log_lik", pars_names, invert = FALSE)
  names_log_lik <- pars_names[index_log_lik]
  if (length(index_log_lik) == 0) {  
    if (Model_type == "Stan") {  ## if log_lik doesn't exist in Stan model
      warning("No log_lik parameter found in Stan model. Log_lik will not be computed even if save_log_lik = TRUE")
    }
  } 
  
  n_params  <- length(bs_names)
  n_params_inc_tp <- length(bs_names_inc_tp)
  n_params_inc_tp_and_gq <- length(bs_names_inc_tp_and_gq) 
  
  bs_index <-  1:n_params
  bs_index_inc_tp <-  1:n_params_inc_tp
  bs_index_inc_tp_and_gq <-  1:n_params_inc_tp_and_gq
  
  # Mow find names, indexes and N's of tp and gq ONLY
  index_tp <- setdiff(bs_index_inc_tp, bs_index)
  index_gq <- setdiff(bs_index_inc_tp_and_gq, bs_index_inc_tp)
  
  names_tp <- pars_names[index_tp]
  names_gq <- pars_names[index_gq]
 
  n_tp <- length(names_tp)
  n_gq <- length(names_gq)
    
  # exclude nuisance params from summary
  index_wo_nuisance <- (n_nuisance + 1):n_par_inc_tp_and_gq
  names_wo_nuisance <- pars_names[index_wo_nuisance]
  n_params_wo_nuisance <- length(names_wo_nuisance) ; n_params_wo_nuisance
  
  index_wo_log_lik <-  grep("^log_lik", pars_names, invert = TRUE)
  names_wo_log_lik <-  pars_names[index_wo_log_lik]
  
  # setdiff(names_wo_nuisance, names_wo_log_lik)
  
  names_wo_nuisance_and_log_lik <-  intersect(names_wo_log_lik, names_wo_nuisance)
  index_wo_nuisance_and_log_lik <-  intersect(index_wo_log_lik, index_wo_nuisance)
  n_params_wo_nuisance_and_log_lik <- length(index_wo_nuisance_and_log_lik)
  
  index_params_main <- index_wo_nuisance_and_log_lik[1]:( index_wo_nuisance_and_log_lik[1] + n_params_main - 1)

  if (sample_nuisance == TRUE) {
    
       try({  
          if (n_nuisance_tracked == n_nuisance) {
            include_nuisance <- TRUE
          } else { 
            include_nuisance <- FALSE
            warning("assumed all nuisance params = 0 as not all nuisance params were tracked during sampling. Hence some outputs (e.g. log_lik) won't be correct.")
          }
       })
    
  } else { 
    
       include_nuisance <- FALSE
     
  }
  
  
  # if (include_nuisance == TRUE)  {
  #   include_log_lik <- TRUE
  # } else { 
  #   include_log_lik <- FALSE
  #   warning("NOTE: can only compute log_lik a posteriori if ALL nuisance are tracked. \n Note it's also possible to compute log_lik during sampling by specifying: save_log_lik_trace = TRUE")
  # }
  
   
  ## print(n_params_main)
   pars_indicies_to_track <- 1:n_par_inc_tp_and_gq
   n_params_full <- n_par_inc_tp_and_gq
   all_param_outs_trace <-    (BayesMVP:::fn_compute_param_constrain_from_trace_parallel(    unc_params_trace_input_main = main_trace,
                                                                                  unc_params_trace_input_nuisance = nuisance_trace, # doesnt matter as not tracking nuisance
                                                                                  pars_indicies_to_track = pars_indicies_to_track,
                                                                                  n_params_full = n_params_full,
                                                                                  n_nuisance = n_nuisance,
                                                                                  n_params_main = n_params_main,
                                                                                  include_nuisance = include_nuisance,
                                                                                  model_so_file = model_so_file,
                                                                                  json_file_path = json_file_path))
  
  
  
  # n_subset_trace <-  n_params_wo_nuisance_and_log_lik
  # indicies_subset <-  index_wo_nuisance_and_log_lik
  # for (kk in 1:n_chains) {
  #   params_subset_trace[[kk]] <- matrix(nrow = n_params_wo_nuisance_and_log_lik, ncol = n_iter)
  #   params_subset_trace[[kk]][, 1:n_iter] <- 
  # }
  ##  subset_out_trace <- array(dim = c(n_params_wo_nuisance_and_log_lik, n_iter, n_chains))

  ### str(all_param_outs_trace[[1]])
   
   comment(print(paste("index_params_main = ", index_params_main)))
   comment(print(paste("index_tp (head) = ", head(index_tp))))
   comment(print(paste("index_gq (head) = ", head(index_gq))))
   
   # if (index_lp == 1) { 
   #   offset <- 1
   # } else { 
   #   offset <- 0
   # }
   
   offset <- 0
   
   comment(print(paste("offset = ", offset)))
   

   comment(print(paste("n_chains = ", n_chains)))
   comment(print(paste("n_iter = ", n_iter)))
   comment(print(paste("n_params_main = ", n_params_main)))
   
   comment(print(paste("length(all_param_outs_trace) = ", length(all_param_outs_trace))))
   comment(print(paste("length(index_params_main) = ", length(index_params_main))))
   
   comment(print(paste("index_params_main = ", index_params_main)))
   
   trace_params_main <- array(dim = c(n_params_main, n_iter, n_chains))
   
   comment(print(str(all_param_outs_trace)))
  
   if (compute_main_params == TRUE) {
     kk <- 1
     all_param_outs_trace[[kk]][index_params_main - offset, 1:n_iter] 
    ##   trace_params_main[1:n_params_main, 1:n_iter, kk] <- all_param_outs_trace[[kk]][index_params_main - offset, 1:n_iter] 
      for (kk in 1:n_chains) {
        try({ 
           trace_params_main[1:n_params_main, 1:n_iter, kk] <- all_param_outs_trace[[kk]][index_params_main - offset, 1:n_iter] 
        }, silent = TRUE)
      }
   }
     
   try({ 
     trace_tp <- NULL
     if (compute_transformed_parameters == TRUE) { 
       trace_tp <- array(dim = c(n_tp, n_iter, n_chains))
       for (kk in 1:n_chains) {
         trace_tp[1:n_tp,  1:n_iter, kk] <- all_param_outs_trace[[kk]][index_tp - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
       }
     }
   })
 
   try({ 
     trace_gq <- NULL
     if (compute_generated_quantities == TRUE) { 
       trace_gq <- array(dim = c(n_gq, n_iter, n_chains))
       for (kk in 1:n_chains) {
         trace_gq[1:n_gq,  1:n_iter, kk] <- all_param_outs_trace[[kk]][index_gq - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
       }
     }
   })
   
   try({ 
     log_lik_trace <- NULL
     if (save_log_lik_trace == TRUE) {
       
           log_lik_trace <- array(dim = c(N, n_iter, n_chains))
       
           if (Model_type == "Stan") {
             
                  for (kk in 1:n_chains) {
                      log_lik_trace[1:N,  1:n_iter, kk] <- all_param_outs_trace[[kk]][index_log_lik - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
                  }
             
           } else {  ## if built-in / manual model
             
                 log_lik_trace <- log_lik_trace_mnl_models
               
           }
       
     }
   })
   
   try({ 
      nuisance_trace <- NULL
      if (save_nuisance_trace == TRUE) {
        index_nuisance <- 1:n_nuisance_tracked
        nuisance_trace <- array(dim = c(n_nuisance_tracked, n_iter, n_chains))
        for (kk in 1:n_chains) {
          nuisance_trace[1:n_nuisance_tracked,  1:n_iter, kk] <- all_param_outs_trace[[kk]][index_nuisance - offset, 1:n_iter] #  params_subset_trace[[kk]][param, 1:n_iter]
        }
      }
   })

   
   ## Remove log-lik trace from tp trace array (since we store it in seperate array called "log_lik_trace" if user chooses to store it)
   trace_tp <- R_fn_remove_log_lik_from_array(trace_tp)
   
   n_cores <- parallel::detectCores()
   
  ### --------- MAIN PARAMETERS / "PARAMETERS" BLOCK IN STAN  ----------------------------------------
  Min_ESS_main <- NULL
  summary_tibble_main_params <- NULL
  names_main <- head(names_wo_nuisance_and_log_lik, n_params_main)
  
  if (compute_main_params == TRUE) { 
    
        summary_tibble_main_params <- BayesMVP:::generate_summary_tibble(  n_threads = n_cores,
                                                                trace = trace_params_main,
                                                                param_names = names_main,
                                                                n_to_compute = n_params_main,
                                                                compute_nested_rhat = compute_nested_rhat,
                                                                n_chains = n_chains, 
                                                                n_superchains = n_superchains)
                  
                  
        Min_ESS_main <- min(na.rm = TRUE, summary_tibble_main_params$n_eff[1:n_params_main])
              
  }
 

  ### --------- GENERATED QUANTITIES ---------------------------------------------------------------- 
  summary_tibble_generated_quantities <- NULL
  
  if  ((compute_generated_quantities == TRUE) && (n_gq > 0))  { 
            
        summary_tibble_generated_quantities <- BayesMVP:::generate_summary_tibble(   n_threads = n_cores,
                                                                          trace = trace_gq,
                                                                          param_names = names_gq,
                                                                          n_to_compute = n_gq,
                                                                          compute_nested_rhat = compute_nested_rhat,
                                                                          n_chains = n_chains, 
                                                                          n_superchains = n_superchains)
                  
    
  }
  
  
  ### --------- TRANSFORMED PARAMETERS -------------------------------------------------------------- 
  summary_tibble_transformed_parameters <- NULL
  if ((compute_transformed_parameters == TRUE) && (n_tp > 0)) {
        
        summary_tibble_transformed_parameters <- BayesMVP:::generate_summary_tibble(     n_threads = n_cores,
                                                                              trace = trace_tp,
                                                                              param_names = names_tp,
                                                                              n_to_compute = n_tp,
                                                                              compute_nested_rhat = compute_nested_rhat,
                                                                              n_chains = n_chains, 
                                                                              n_superchains = n_superchains)
    
  }
  
   
  #### -----------------------------  DF / tibble creation (for "posterior" and "bayesplot" R packages)  --------------------------------------------
  trace_params_main_tibble <- NULL ;   trace_params_main_reshaped <- NULL
  trace_transformed_params_tibble <- NULL ;   trace_tp_reshaped <- NULL
  trace_generated_quantities_tibble <- NULL ; trace_gq_reshaped <- NULL

  ## re-shape data (so can easily use with "posterior" & "bayesplot" R packages)
  trace_params_main_reshaped <-  base::aperm(trace_params_main, c(2, 3, 1))
  if (compute_transformed_parameters == TRUE)   trace_tp_reshaped <- base::aperm(trace_tp, c(2, 3, 1))
  if (compute_generated_quantities == TRUE)     trace_gq_reshaped <- base::aperm(trace_gq, c(2, 3, 1))
  
  ## add the parameter names to the array
  dimnames(trace_params_main_reshaped) <- list(iterations = 1:n_iter, 
                                               chains = 1:n_chains, 
                                               parameters = names_main)
  if (compute_transformed_parameters == TRUE)  dimnames(trace_tp_reshaped) <- list(iterations = 1:n_iter, 
                                                                                   chains = 1:n_chains,
                                                                                   parameters = names_tp)
  if (compute_generated_quantities == TRUE) dimnames(trace_gq_reshaped) <-    list(iterations = 1:n_iter, 
                                                                                   chains = 1:n_chains, 
                                                                                   parameters = names_gq)
  
   ## then convert from array -> to df/tibble format  
  if (save_trace_tibbles == TRUE) {
       trace_params_main_tibble <- dplyr::tibble(posterior::as_draws_df(trace_params_main_reshaped))
       if (compute_transformed_parameters == TRUE)  trace_transformed_params_tibble <- dplyr::tibble(posterior::as_draws_df(trace_tp_reshaped))
       if (compute_generated_quantities == TRUE) trace_generated_quantities_tibble <-  dplyr::tibble(posterior::as_draws_df(trace_gq_reshaped))
  }
  
   ##  ----- Make OVERALL (full) draws array   --------------- 
   dims_main <- dim(trace_params_main_reshaped)[3]
   dims_tp <- dim(trace_tp_reshaped)[3]
   dims_gq <- dim(trace_gq_reshaped)[3]
   
   dim_total <- dims_main
   if (compute_transformed_parameters == TRUE) dim_total <- dim_total + dims_tp
   if (compute_generated_quantities == TRUE)   dim_total <- dim_total + dims_gq
   
   
   names_total <- names_main
   if (compute_transformed_parameters == TRUE)  names_total <- c(names_total, names_tp)
   if (compute_generated_quantities == TRUE)    names_total <- c(names_total, names_gq)
   
   draws_array <- array(NA, dim = c(n_iter, n_chains, dim_total))
   draws_array[1:n_iter, 1:n_chains, 1:n_params_main] <- trace_params_main_reshaped # main first 
   
   if (compute_transformed_parameters == TRUE)   {
     draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_tp)] <- trace_tp_reshaped
   }
   if ((compute_generated_quantities == TRUE) &&  (compute_transformed_parameters == TRUE))   {
     draws_array[1:n_iter, 1:n_chains, (n_params_main + n_tp + 1):(n_params_main + n_tp + n_gq)] <- trace_gq_reshaped
   }
   
   if ((compute_generated_quantities == TRUE) &&  (compute_transformed_parameters == FALSE))   {
     draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_gq)] <- trace_gq_reshaped
   }
   if ((compute_generated_quantities == FALSE) &&  (compute_transformed_parameters == TRUE))   {
     draws_array[1:n_iter, 1:n_chains, (n_params_main + 1):(n_params_main + n_tp)] <- trace_tp_reshaped
   }

   
   dimnames(draws_array) <- list( iterations = 1:n_iter, 
                                  chains = 1:n_chains, 
                                  parameters = names_total)
   
  

  
  
  
  
  try({
    {
      print(tictoc::toc(log = TRUE))
      log.txt <- tictoc::tic.log(format = TRUE)
      tictoc::tic.clearlog()
      time_summaries <- unlist(log.txt)
    }
  })
  
  try({
    time_summaries <- as.numeric(substr(start = 0, stop = 100, strsplit(strsplit(time_summaries, "[:]")[[1]], "[s]")  [[1]][1]))
  })
  
  
  
}

  
    time_total <- time_summaries + time_burnin + time_sampling

    ESS_per_sec_samp <- Min_ESS_main / time_sampling
    ESS_per_sec_total <- Min_ESS_main / time_total
    
    EHMC_args_as_Rcpp_List <- model_results$init_burnin_object$EHMC_args_as_Rcpp_List
    
    try({
      comment(print((paste("Min ESS (parameters block, main only) = ", round(Min_ESS_main, 0)))))
      comment(print((paste("Min ESS / sec [samp.] (parameters block, main only) = ", signif(ESS_per_sec_samp, 3)))))
      comment(print((paste("Min ESS / sec [overall] (parameters block, main only) = ", signif(ESS_per_sec_total, 3)))))
    })
    try({ 
      n_grad_evals_sampling_main <-  (EHMC_args_as_Rcpp_List$tau_main / EHMC_args_as_Rcpp_List$eps_main)  * n_iter * n_chains_sampling
      Min_ess_per_grad_main_samp <-  Min_ESS_main / n_grad_evals_sampling_main
    })
    try({
      n_grad_evals_sampling_us <-  (EHMC_args_as_Rcpp_List$tau_us / EHMC_args_as_Rcpp_List$eps_us)  * n_iter * n_chains_sampling
      Min_ess_per_grad_us_samp <-  Min_ESS_main / n_grad_evals_sampling_us
    })
    try({ 
      
      if (partitioned_HMC == TRUE) {
          weight_nuisance_grad <- 0.3333333
          weight_main_grad <- 0.6666667 ## main grad takes ~ 2x as long to compute as nuisance grad 
          Min_ess_per_grad_samp_weighted <- (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) / (weight_nuisance_grad + weight_main_grad)
      } else {  # if not partitioned, grad isnt seperate and "main" grad = "all" grad so use weight_main_grad only !!
          weight_nuisance_grad <- 0.00
          weight_main_grad <- 1.00
          Min_ess_per_grad_samp_weighted <- (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) / (weight_nuisance_grad + weight_main_grad)
      }
      
      comment(print((paste("Min ESS / grad [samp., weighted] (parameters block, main only) = ", signif(1000 *  Min_ess_per_grad_samp_weighted, 3)))))
    })
    try({ 
      grad_evals_per_sec <- ESS_per_sec_samp / Min_ess_per_grad_samp_weighted
      comment(print((paste("Grad evals / sec [samp.] (parameters block, main only) = ", signif(grad_evals_per_sec, 3)))))
    })
    try({ 
          sampling_time_to_Min_ESS <- time_sampling
          sampling_time_to_100_ESS <- (100 / Min_ESS_main) * sampling_time_to_Min_ESS
          sampling_time_to_1000_ESS <- (1000 / Min_ESS_main) * sampling_time_to_Min_ESS
          sampling_time_to_10000_ESS <- (10000 / Min_ESS_main) * sampling_time_to_Min_ESS
          
          ## w/o summary time
          total_time_to_100_ESS_wo_summaries <-   time_burnin + sampling_time_to_100_ESS
          total_time_to_1000_ESS_wo_summaries <-  time_burnin + sampling_time_to_1000_ESS
          total_time_to_10000_ESS_wo_summaries <- time_burnin + sampling_time_to_10000_ESS
          
          ## w/ summary time (note: assuming that time_summaries scales linearly w/ 
          # the min ESS required (i.e. n_iter and/or n_chains) Might be over-estimate)
          summary_time_to_Min_ESS <- time_summaries
          summary_time_to_100_ESS <- (100 / Min_ESS_main) * summary_time_to_Min_ESS
          summary_time_to_1000_ESS <- (1000 / Min_ESS_main) * summary_time_to_Min_ESS
          summary_time_to_10000_ESS <- (10000 / Min_ESS_main) * summary_time_to_Min_ESS
          
          total_time_to_100_ESS_with_summaries <-   time_burnin + sampling_time_to_100_ESS   + summary_time_to_100_ESS
          total_time_to_1000_ESS_with_summaries <-  time_burnin + sampling_time_to_1000_ESS  + summary_time_to_1000_ESS
          total_time_to_10000_ESS_with_summaries <- time_burnin + sampling_time_to_10000_ESS + summary_time_to_10000_ESS
    })
    
    if (n_divs > 0) { 
      warning("divergences detected!")
      try({ 
        comment(print(paste("Number of divergences = ", n_divs)))
        comment(print(paste("% divergences = ", pct_divs)))
      })
    } else { 
      message("No divergences detected!")
    }
    
    
    ### ---------  Make R lists to output  --------------------------------------
    ## list to store summary tibbles / DF's
    summary_tibbles <- list( summary_tibble_main_params = summary_tibble_main_params,
                             summary_tibble_transformed_parameters = summary_tibble_transformed_parameters,
                             summary_tibble_generated_quantities = summary_tibble_generated_quantities)
    
    ## list to store traces (as 3D arrays)

    traces_as_arrays <- list(draws_array = draws_array,
                             trace_params_main = trace_params_main_reshaped,
                             trace_transformed_params = trace_tp_reshaped,
                             trace_generated_quantities = trace_gq_reshaped)
    
    ## list to store traces (as tibbles/DF's)
    traces_as_tibbles <- list( trace_params_main_tibble = trace_params_main_tibble,
                               trace_transformed_params_tibble = trace_transformed_params_tibble,
                               trace_generated_quantities_tibble = trace_generated_quantities_tibble)
    
    HMC_info <- list( tau_main = EHMC_args_as_Rcpp_List$tau_main,
                      eps_main = EHMC_args_as_Rcpp_List$eps_main,
                      tau_us = EHMC_args_as_Rcpp_List$tau_us,
                      eps_us = EHMC_args_as_Rcpp_List$eps_us,
                      n_chains_sampling = n_chains_sampling,
                      n_chains_burnin = n_chains_burnin,
                      n_iter = n_iter,
                      n_burnin = n_burnin,
                      LR_main = LR_main,
                      LR_us = LR_us,
                      adapt_delta = adapt_delta,
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
                      multi_attempts = multi_attempts)
    
    ## list to store efficiency information
    efficiency_info <- list(              n_iter = n_iter,
                                          Min_ESS_main = Min_ESS_main, 
                                          ESS_per_sec_samp = ESS_per_sec_samp, 
                                          ESS_per_sec_total = ESS_per_sec_total,
                                          time_burnin = time_burnin, 
                                          time_sampling = time_sampling, 
                                          time_summaries = time_summaries,
                                          time_total_wo_summaries = time_total_wo_summaries, 
                                          time_total = time_total, 
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
    
    
    ### final lists to output 
    if (save_log_lik_trace == FALSE) log_lik_trace <- NULL
    if (save_nuisance_trace == FALSE)  nuisance_trace <- NULL
    traces <- list(traces_as_arrays = traces_as_arrays,
                   traces_as_tibbles = traces_as_tibbles,
                   log_lik_trace = log_lik_trace,
                   nuisance_trace = nuisance_trace)
    
    divergences <- list( n_divs = n_divs, 
                         pct_divs = pct_divs)
    
    
    summaries <- list(summary_tibbles = summary_tibbles,
                      divergences = divergences,
                      efficiency_info = efficiency_info, 
                      HMC_info = HMC_info)
    
    
    output_list <- list(  summaries = summaries, ### summary info (incl. efficiency info + divergences)
                          traces = traces       ### trace arrays + trace tibbles
                          ## all_param_outs_trace = all_param_outs_trace
    )
    
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



# Now call it with your actual data:
# summary_table <- create_stan_summary(your_trace_vector, pars_names_wo_nuisance)




