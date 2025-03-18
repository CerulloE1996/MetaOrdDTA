




#' init_EHMC_Metric_as_Rcpp_List
#' @keywords internal
#' @export
init_EHMC_Metric_as_Rcpp_List   <- function(metric_shape) {
        
     EHMC_Metric_as_Rcpp_List <- list(metric_shape = metric_shape)
     return(EHMC_Metric_as_Rcpp_List)
   
}




#' init_EHMC_args_as_Rcpp_List
#' @keywords internal
#' @export
init_EHMC_args_as_Rcpp_List   <- function(diffusion_HMC) {

      tau <- 1
      eps <- 1

      EHMC_args_as_Rcpp_List <- list( tau = tau,
                                      eps = eps)
        
      return(EHMC_args_as_Rcpp_List)

}



#' init_EHMC_burnin_as_Rcpp_List
#' @keywords internal
#' @export
init_EHMC_burnin_as_Rcpp_List   <- function( adapt_delta,
                                             max_treedepth
                                             ) {
  
      adapt_delta <- adapt_delta
      max_treedepth <- max_treedepth
      
      # ------------ put in the list to put into C++
      EHMC_burnin_as_Rcpp_List <- list(   ### for main params
                                          adapt_delta = adapt_delta,
                                          max_treedepth = max_treedepth)
      
      
      return(EHMC_burnin_as_Rcpp_List)
      
}




#' init_and_run_burnin
#' @keywords internal
#' @export
init_and_run_burnin   <- function(  Model_type,
                                    init_object,
                                    Stan_data_list,
                                    # model_args_list,
                                    y,
                                    N,
                                    parallel_method,
                                    ##
                                    manual_tau,
                                    tau_if_manual,
                                    ##
                                    n_chains_burnin,
                                    sample_nuisance,
                                    Phi_type,
                                    n_params_main,
                                    n_nuisance,
                                    seed,
                                    n_burnin,
                                    adapt_delta,
                                    LR_main,
                                    LR_us,
                                    n_adapt,
                                    partitioned_HMC,
                                    diffusion_HMC,
                                    clip_iter,
                                    gap ,
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
                                    interval_width_main,
                                    interval_width_nuisance,
                                    force_autodiff,
                                    force_PartialLog,
                                    multi_attempts,
                                    n_nuisance_to_track,
                                    Model_args_as_Rcpp_List) {
  
  
  n_params <- n_params_main + n_nuisance
  index_us <- 1:n_nuisance
  index_main <- (n_nuisance + 1):n_nuisance
 
  
  {  # ---------------------------------------------------------------- list for EHMC params / EHMC struct in C++
    
        EHMC_Metric_as_Rcpp_List <- init_EHMC_Metric_as_Rcpp_List(   n_params_main = n_params_main, 
                                                                     n_nuisance = n_nuisance, 
                                                                     metric_shape_main = metric_shape_main)  
    
  }
  
  
  
  { # ----------------------------------------------------------------- list for EHMC params / EHMC struct in C++  #
    
        EHMC_args_as_Rcpp_List <- init_EHMC_args_as_Rcpp_List(diffusion_HMC = diffusion_HMC) 
    
  }

 
  {  # ----------------------------------------------------------------- list for EHMC params / EHMC struct in C++
      
        EHMC_burnin_as_Rcpp_List <- init_EHMC_burnin_as_Rcpp_List( n_params_main = n_params_main,
                                                                   n_nuisance = n_nuisance,
                                                                   adapt_delta = adapt_delta,
                                                                   LR_main = LR_main,
                                                                   LR_us = LR_us)
    
  }
  

  theta_main_vectors_all_chains_input_from_R <- init_object$theta_main_vectors_all_chains_input_from_R
  theta_us_vectors_all_chains_input_from_R <-   init_object$theta_us_vectors_all_chains_input_from_R
  
  # print(theta_main_vectors_all_chains_input_from_R)
  # print(theta_us_vectors_all_chains_input_from_R)
  
  ###  theta_vec_mean <- rowMeans(rbind(theta_us_vectors_all_chains_input_from_R, theta_main_vectors_all_chains_input_from_R))
  
 
 # print(theta_main_vectors_all_chains_input_from_R)
  

  
  EHMC_args_as_Rcpp_List$diffusion_HMC <- diffusion_HMC
  
  

 ##  Model_args_as_Rcpp_List$Model_args_strings[1, 2] <-  Phi_type
  
  results_burnin <-  BayesMVP:::R_fn_EHMC_SNAPER_ADAM_burnin( Model_type = Model_type,
                                                              sample_nuisance = sample_nuisance,
                                                              parallel_method = parallel_method,
                                                              n_params_main = n_params_main,
                                                              n_nuisance = n_nuisance,
                                                              Stan_data_list = Stan_data_list,
                                                              # model_args_list = model_args_list,
                                                              y = y,
                                                              N = N,
                                                              ##
                                                              manual_tau = manual_tau, ##
                                                              tau_if_manual = tau_if_manual, ##
                                                              ##
                                                              n_chains_burnin = n_chains_burnin,
                                                              seed = seed,
                                                              n_burnin = n_burnin,
                                                              LR_main = LR_main,
                                                              LR_us = LR_us,
                                                              n_adapt = n_adapt,
                                                              partitioned_HMC = partitioned_HMC,
                                                              clip_iter = clip_iter,
                                                              gap = gap,
                                                              metric_type_main = metric_type_main,
                                                              metric_shape_main = metric_shape_main,
                                                              metric_type_nuisance = metric_type_nuisance,
                                                              metric_shape_nuisance = metric_shape_nuisance,
                                                              max_eps_main = max_eps_main,
                                                              max_eps_us = max_eps_us,
                                                              max_L = max_L,
                                                              tau_mult = tau_mult,
                                                              ratio_M_us = ratio_M_us,
                                                              ratio_M_main = ratio_M_main,
                                                              interval_width_main = interval_width_main,
                                                              interval_width_nuisance = interval_width_nuisance,
                                                              force_autodiff = force_autodiff,
                                                              force_PartialLog = force_PartialLog,
                                                              multi_attempts = multi_attempts,
                                                              theta_main_vectors_all_chains_input_from_R = theta_main_vectors_all_chains_input_from_R,
                                                              theta_us_vectors_all_chains_input_from_R = theta_us_vectors_all_chains_input_from_R,
                                                              n_nuisance_to_track = n_nuisance_to_track,
                                                              Model_args_as_Rcpp_List = Model_args_as_Rcpp_List,
                                                              EHMC_args_as_Rcpp_List = EHMC_args_as_Rcpp_List,
                                                              EHMC_Metric_as_Rcpp_List = EHMC_Metric_as_Rcpp_List,
                                                              EHMC_burnin_as_Rcpp_List = EHMC_burnin_as_Rcpp_List)
  
  
  ### update Rcpp lists to pass onto sampling / post-burnin function 
  time_burnin <- results_burnin$time_burnin
  EHMC_args_as_Rcpp_List <- results_burnin$EHMC_args_as_Rcpp_List
  EHMC_Metric_as_Rcpp_List <- results_burnin$EHMC_Metric_as_Rcpp_List
  EHMC_burnin_as_Rcpp_List <- results_burnin$EHMC_burnin_as_Rcpp_List
  
  ## for init values for sampling phase
  theta_main_vectors_all_chains_input_from_R <- results_burnin$theta_main_vectors_all_chains_input_from_R
  theta_us_vectors_all_chains_input_from_R <- results_burnin$theta_us_vectors_all_chains_input_from_R
  
  return(list(time_burnin = time_burnin,
              Model_args_as_Rcpp_List = Model_args_as_Rcpp_List,
              EHMC_Metric_as_Rcpp_List = EHMC_Metric_as_Rcpp_List,
              EHMC_args_as_Rcpp_List = EHMC_args_as_Rcpp_List,
              EHMC_burnin_as_Rcpp_List = EHMC_burnin_as_Rcpp_List,
              theta_main_vectors_all_chains_input_from_R = theta_main_vectors_all_chains_input_from_R,
              theta_us_vectors_all_chains_input_from_R = theta_us_vectors_all_chains_input_from_R))
  
}






















