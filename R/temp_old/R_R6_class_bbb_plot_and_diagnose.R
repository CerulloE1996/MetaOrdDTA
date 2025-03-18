
require(R6)

#' @title MetaOrd "plot and diagnose" Class
#' @description Class for handling MetaOrd model summaries, samples (trace) and diagnostics.
#'  
#' 
#' @details 
#' This class is called via the main "MetaOrd_model" R6 class. 
#' 
#' @section Model Types:
#' * meta_ord: Meta-analysis model for a single ordinal test (e.g., questionnaires). Uses model proposed by Cerullo et al. 
#'   Supports missing thresholds in some or all studies.
#'   
#' * meta_cts: Meta-analysis model for a single continuous test (e.g., biomarkers). Uses the model proposed by Jones et al. 
#'   Supports missing thresholds in some or all studies.
#'   
#' * network_meta_ord: Meta-analysis model for multiple ordinal tests (e.g., questionnaires). Uses the model proposed by Cerullo et al. 
#'   and NMA model proposed by Nyaga et al. Supports missing thresholds in some or all studies.
#'   
#' * network_meta_cts: Meta-analysis model for multiple continuous tests (e.g., biomarkers). Uses the model proposed by Jones et al. and 
#'   Nyaga et al. Supports missing thresholds in some or all studies.
#' 
#' @section Class Relationship:
#' This class is not meant to be instantiated/called directly by users. Instead, it is created 
#' and returned by the MetaOrd_model$summary() method. It provides methods for:
#' * Extracting MCMC diagnostics (divergences, ESS, etc.)
#' * Creating trace and density plots
#' * Computing efficiency metrics
#' * Accessing posterior draws in various formats
#' 
#' 
#' @section Common Usage Patterns:
#' \describe{
#'   \item{Basic Diagnostics}{get_divergences() -> get_efficiency_metrics()}
#'   \item{Parameter Summaries}{get_summary_main() -> get_summary_transformed()}
#'   \item{Visualization}{plot_traces() -> plot_densities()}
#'   \item{Custom Analysis}{get_posterior_draws() -> Your analysis}
#' }
#' 
#' 
#' 
#' 
#' @export
MetaOrd_plot_and_diagnose <- R6Class("MetaOrd_plot_and_diagnose",
                             
                               public = list(
                                 
                                    ##
                                    ## ---- Must define the fields here first (since this is a seperate R6 class):
                                    ##
                                    model_samples_obj = NULL,
                                    model_summary_and_trace_obj = NULL,
                                    ##
                                    #' @field model_samples_obj The model samples object, obtained from the main "MetaOrd_model" R6 class. 
                                    #' @field model_summary_and_trace_obj The summary + trace object, obtained from the main "MetaOrd_model" R6 class.
                                    ##
                                    #' @param model_samples_obj The model samples object, obtained from the main "MetaOrd_model" R6 class. 
                                    #' @param model_summary_and_trace_obj The summary + trace object, obtained from the main "MetaOrd_model" R6 class.
                                    ##
                                    initialize = function( model_samples_obj,
                                                           model_summary_and_trace_obj) {
                                              ##
                                              ## ---- Basic checks:
                                              ##
                                              if (is.null(model_samples_obj)) {
                                                stop("No model samples object ('model_samples_obj') available")
                                              }
                                              if (is.null(model_summary_and_trace_obj)) {
                                                stop("No model summary + trace object ('model_summary_and_trace_obj') available")
                                              }
                                              ##
                                              ## ---- Store models directly, not nesting them:
                                              ##
                                              self$model_samples_obj <- model_samples_obj
                                              self$model_summary_and_trace_obj <- model_summary_and_trace_obj
                                              ##
                                              ## ---- Print info about what we stored to help with debugging:
                                              ##
                                              cat("MetaOrd_plot_and_diagnose initialized with:\n")
                                              cat("- model_samples_obj class:", class(model_samples_obj), "\n")
                                              print(str(model_samples_obj))
                                              cat("- model_summary_and_trace_obj class:", class(model_summary_and_trace_obj), "\n")
                                              print(str(model_summary_and_trace_obj))
                                              
                                              # invisible(self)
                                          
                                    },
                                    ##
                                    ## ---- Convenience methods for summaries (inc. divergences, R-hat, ESS, posterior means, etc):
                                    ##
                                    ## ---- Get divergence info:
                                    ##
                                    #'@return Returns self$model_summary_and_trace_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                                    #'self$model_summary_and_trace_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                                    get_divergences = function() {
                                         return(self$model_summary_and_trace_obj$HMC_info$HMC_diagnostic_info$divergences)
                                    },
                                    ##
                                    ## ---- HMC diagnostic info:
                                    ##
                                    #'@return Returns self$model_summary_and_trace_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                                    #'self$model_summary_and_trace_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                                    get_HMC_diagnostic_info = function() { 
                                         return(self$model_summary_and_trace_obj$HMC_info$HMC_diagnostic_info) ## Extract the "HMC_info" list
                                    },
                                    ##
                                    ## --------- Convenience method for HMC algorithm metrics:
                                    ##
                                    #' @return Returns a list containing HMC algorithm info (including timings and some efficiency metrics)
                                    get_HMC_info = function() { 
                                      return(self$model_summary_and_trace_obj$HMC_info) ## Extract the "HMC_info" list
                                    },
                                    ##
                                    ## ---- Summary tibble methods:
                                    ##
                                    #'@return Returns main parameter summaries as a tibble.
                                    ##
                                    get_summary_main = function() {  #### --------------------------------------------------- WORKING
                                         return(self$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_main_params)
                                    },
                                    #'@return Returns transformed parameter summaries as a tibble.
                                    ##
                                    get_summary_transformed = function() {
                                         return(self$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters)
                                    },
                                    #'@return Returns generated quantities as a tibble.
                                    ##
                                    get_summary_generated_quantities = function() {  #### --------------------------------------------------- WORKING
                                         return(self$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_generated_quantities)
                                    },
                                    ##
                                    ## ---- Posterior draws / 3D traces methods:
                                    ##
                                    #'Convenience method for ALL traces
                                    #'@return Returns all traces including as tibbles (if enabled) and as 3D array's (iterations, chains, parameters). 
                                    ##
                                    get_all_traces_list = function() {  #### --------------------------------------------------- WORKING
                                          return(self$model_summary_and_trace_obj$traces)
                                    },
                                    #'Convenience method for log_lik trace
                                    #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                                    get_trace_main = function() { #### --------------------------------------------------- WORKING
                                          return(self$model_summary_and_trace_obj$traces$traces_as_arrays$trace_params_main)
                                    },
                                    #'Convenience method for log_lik trace
                                    #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                                    get_trace_transformed = function() {  #### --------------------------------------------------- WORKING
                                          return(self$model_summary_and_trace_obj$traces$traces_as_arrays$trace_tp_wo_log_lik)
                                    },
                                    #'Convenience method for log_lik trace
                                    #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                                    get_trace_generated_quantities = function() { #### --------------------------------------------------- WORKING
                                          return(self$model_summary_and_trace_obj$traces$traces_as_arrays$trace_gq)
                                    },
                                    #'Convenience method for log_lik trace
                                    #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                                    get_trace_log_lik = function() {  #### --------------------------------------------------- WORKING
                                          return(self$model_summary_and_trace_obj$traces$traces_as_arrays$trace_log_lik)
                                    },
                                    #'Convenience method for posterior draws (as tibbles)
                                    #'@return Returns trace as tibbles. The trace gets returned as 3 seperate tibbles (one tibble for the 
                                    #'main parameters, one for the transformed parameters, and one for generated quantities). 
                                    ##
                                    get_traces_as_tibbles = function() {  #### --------------------------------------------------- WORKING
                                        list(
                                           trace_as_tibble_main_params          =  self$model_summary_and_trace_obj$traces$traces_as_tibbles$trace_params_main_tibble,
                                           trace_as_tibble_transformed_params   =  self$model_summary_and_trace_obj$traces$traces_as_tibbles$trace_transformed_params_tibble,
                                           trace_as_tibble_generated_quantities =  self$model_summary_and_trace_obj$traces$traces_as_tibbles$trace_generated_quantities_tibble
                                        )
                                    },
                                    ##
                                    ## --------- Convenience method for efficiency metrics:
                                    ##
                                    #' @return Returns a list containing efficiency metrics including:
                                    #' \itemize{
                                    #'   \item time_burnin: Time spent in burnin phase
                                    #'   \item time_sampling: Time spent sampling
                                    #'   \item time_total_MCMC: Total time spent doing MCMC/HMC sampling, excluding time to compute summaries and diagnostics.
                                    #'   \item time_total_inc_summaries: Total time spent doing MCMC/HMC sampling, including time to compute summaries and diagnostics.
                                    #'   \item Min_ESS_main: The minimum ESS of the main model parameters. 
                                    #'   \item Min_ESS_per_sec_sampling: The minimum ESS per second for the sampling phase. 
                                    #'   \item Min_ESS_per_sec_total: The minimum ESS per second for the total model run time, including any time spent computing
                                    #'   summaries and diagnostics. 
                                    #'   \item Min_ESS_per_grad_sampling:  The minimum ESS per gradient evaluation for the sampling phase. 
                                    #'   \item grad_evals_per_sec: The number of gradient evaluations performed per second. 
                                    #'   \item est_time_to_100_ESS_sampling: The estimated sampling time to reach a minimum ESS of 100.
                                    #'   \item est_time_to_1000_ESS_sampling: The estimated sampling time to reach a minimum ESS of 1000.
                                    #'   \item est_time_to_10000_ESS_sampling: The estimated sampling time to reach a minimum ESS of 10,000.
                                    #'   \item est_time_to_100_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                                    #'   to reach a minimum ESS of 100.
                                    #'   \item est_time_to_1000_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                                    #'   to reach a minimum ESS of 1000.
                                    #'   \item est_time_to_10000_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                                    #'   to reach a minimum ESS of 10,000.
                                    #'   \item est_time_to_100_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                                    #'   diagnostics) to reach a minimum ESS of 100.
                                    #'   \item est_time_to_1000_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                                    #'   diagnostics) to reach a minimum ESS of 1000.
                                    #'   \item est_time_to_10000_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                                    #'   diagnostics) to reach a minimum ESS of 10,000.
                                    #' }
                                    ##
                                    get_efficiency_metrics = function() {
                                                return(self$model_summary_and_trace_obj$efficiency_info) ## Extract the "efficiency_info" list
                                    },
                                    ##
                                    ## --------- Convenience fn to compute "time to X ESS" using the "$time_to_target_ESS()" method
                                    ##
                                    #'@param target_ESS The target ESS. 
                                    #'@return Returns a list called "target_ESS_times" which contains the estimated sampling time to reach the target ESS 
                                    #'("sampling_time_to_target_ESS"), the estimated total time excluding the time it takes to compute model summaries 
                                    #'("total_time_to_target_ESS_wo_summaries"), and the total estimated time to reach the 
                                    #'target ESS ("total_time_to_target_ESS_with_summaries").
                                    time_to_target_ESS = function(target_ESS) {
                                          
                                              if (is.null(self$model_summary_and_trace_obj)) {
                                                stop("No summary object available")   
                                              }
                                          
                                              ## get required efficiency info first
                                              Min_ESS_main   <- self$model_summary_and_trace_obj$efficiency_info$Min_ESS_main
                                              time_burnin    <- self$model_summary_and_trace_obj$efficiency_info$time_burnin
                                              time_sampling  <- self$model_summary_and_trace_obj$efficiency_info$time_sampling
                                              time_summaries <- self$model_summary_and_trace_obj$efficiency_info$time_summaries
                                              
                                              sampling_time_to_Min_ESS <- time_sampling
                                              
                                              ## est. sampling time to target_ESS
                                              sampling_time_to_target_ESS <- (target_ESS / Min_ESS_main) * sampling_time_to_Min_ESS
                                              
                                              ## total w/o summary time
                                              total_time_to_target_ESS_wo_summaries <-   time_burnin + sampling_time_to_target_ESS
                                              
                                              ## total w/ summary time
                                              summary_time_to_Min_ESS <- time_summaries
                                              summary_time_to_target_ESS <- (target_ESS / Min_ESS_main) * summary_time_to_Min_ESS
                                              
                                              total_time_to_target_ESS_with_summaries <-   time_burnin + sampling_time_to_target_ESS   + summary_time_to_target_ESS
                                              
                                              target_ESS_times <- list(sampling_time_to_target_ESS = sampling_time_to_target_ESS,
                                                                       total_time_to_target_ESS_wo_summaries = total_time_to_target_ESS_wo_summaries,
                                                                       total_time_to_target_ESS_with_summaries = total_time_to_target_ESS_with_summaries)
                                              
                                              return(target_ESS_times)
                                      
                                      
                                    },
                                    ##
                                    ## --------- Convenience fn to compute "n_iter to X ESS" using the "$iter_to_target_ESS()" method
                                    ##
                                    #'@param target_ESS The target ESS. 
                                    #'@return Returns a list called "target_ESS_iter" which contains the estimated number of sampling iterations (n_iter) to 
                                    #'reach the target ESS ("sampling_iter_to_target_ESS").
                                    iter_to_target_ESS = function(target_ESS) {
                                        
                                              if (is.null(self$model_summary_and_trace_obj)) {
                                                stop("No summary object available")   
                                              }
                                              
                                              ## get required efficiency info first:
                                              Min_ESS_main   <- self$model_summary_and_trace_obj$efficiency_info$Min_ESS_main
                                              n_iter         <- self$model_summary_and_trace_obj$efficiency_info$n_iter
                                              n_iter_to_min_ESS <- n_iter
                                              
                                              ## est. sampling time to target_ESS
                                              sampling_iter_to_target_ESS <- (target_ESS / Min_ESS_main) * n_iter_to_min_ESS
                                              
                                              
                                              target_ESS_iter <- list(sampling_iter_to_target_ESS = sampling_iter_to_target_ESS)
                                              
                                              return(target_ESS_iter)
                                        
                                    },
                                    ##
                                    ## --------- MCMC trace plots method
                                    ##
                                    #'@param debugging Whether debugging mode is turned on (\code{TRUE}) or off (\code{FALSE}). The default is FALSE.
                                    #'@param params The parameters to generate trace plots for. This is a character vector - e.g. to plot trace plots for 
                                    #'beta: params = c("beta"). 
                                    #'The default is NULL and the trace plots for all model parameters (which have a trace array) will be plotted. 
                                    #'@param batch_size The number of trace plots to display per panel. Default is 9. 
                                    #'@return If no parameters specified (i.e. params = NULL), then this will return an object containing the trace plots for all
                                    #'model parameters  which have a trace array.
                                    plot_traces = function(debugging = debugging,
                                                           params = NULL, 
                                                           batch_size = 9
                                                           ) {
                                
                                              if (is.null(self$model_summary_and_trace_obj)) {
                                                stop("No summary object available")   
                                              }
                                              
                                              if (!is.null(params) && !is.character(params)) {
                                                stop("params must be NULL or a character vector")
                                              }
                                
                                              if (!is.numeric(batch_size) || batch_size < 1) {
                                                stop("batch_size must be a positive integer")
                                              }
                                              
                                              if (debugging) {
                                                  cat("Summary object structure:\n")
                                                  str(self$model_summary_and_trace_obj)
                                              }
                                              
                                              # Get the draws array
                                              draws_array <- self$model_summary_and_trace_obj$traces$traces_as_arrays$draws_array 
                                              
                                              if (is.null(draws_array)) {
                                                stop("draws_array is NULL in model_summary_and_trace_obj")
                                              }
                                
                                              if (is.null(params)) { # plot all params
                                                
                                                     bayesplot::mcmc_trace(draws_array)
                                                
                                              } else {   # plot specific params using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot 
                                              
                                                   MetaOrdDTA:::plot_multiple_params_batched( draws_array = draws_array, 
                                                                                              param_string = param_string,
                                                                                              condition = "exact_match", 
                                                                                              plot_type = "trace",
                                                                                              batch_size = batch_size,
                                                                                              print = FALSE)
                                              }
                                      
                                    },
                                    ##
                                    ## --------- MCMC density plots method
                                    ##
                                    #'@param params The parameters to generate posterior density plots for. This is a character vector - e.g. to plot the 
                                    #'posterior density plots for beta: params = c("beta"). 
                                    #'The default is NULL and the posterior density plots for all model parameters (which have a trace array) will be plotted. 
                                    #'@param batch_size The number of posterior density plots to display per panel. Default is 9. 
                                    #'@return If no parameters specified (i.e. params = NULL), then this will return an object containing the posterior density
                                    #'plots for all model parameters which have a trace array.
                                    plot_densities = function(debugging = debugging,
                                                              params = NULL, 
                                                              batch_size =  9
                                                              ) {
                                      
                                              if (is.null(self$model_summary_and_trace_obj)) {
                                                stop("No summary object available")  
                                              }
                              
                                              if (!is.null(params) && !is.character(params)) {
                                                stop("params must be NULL or a character vector")
                                              }
                                              
                                              if (!is.numeric(batch_size) || batch_size < 1) {
                                                stop("batch_size must be a positive integer")
                                              }
                              
                                              if (debugging) {
                                                cat("Summary object structure:\n")
                                                str(self$model_summary_and_trace_obj)
                                              }
                              
                                              # Get the draws array
                                              draws_array <- self$model_summary_and_trace_obj$traces$traces_as_arrays$draws_array
                                              
                                              if (is.null(params)) { # plot all params
                                                
                                                    bayesplot::mcmc_dens(draws_array)
                                                
                                              } else {   # plot specific params using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot
                                                    
                                                    MetaOrdDTA:::plot_multiple_params_batched(    draws_array = draws_array, 
                                                                                                  param_prefixes = params, 
                                                                                                  plot_type = "density", 
                                                                                                  batch_size = batch_size)
                                              }
                                                        
                                    },
                                    ##
                                    ## --------- Test accuracy-specific plots method: sROC plot (w/ confidence + prediction regions/curves):
                                    ##
                                    #'@param df_true ... 
                                    #'@param conf_region_colour ...
                                    #'@param pred_region_colour ...
                                    #'@return sROC plot(s).
                                    #'
                                    plot_sROC = function( df_true = NULL,
                                                          conf_region_colour = "blue", 
                                                          pred_region_colour = "blue",
                                                          Stan_model_file_name = self$model_samples_obj$Stan_model_file_name,
                                                          Stan_mod_samples     = self$model_samples_obj$full_outputs$Stan_model_sample_output$Stan_mod_samples
                                                          ) {
                                      
                                              if (is.null(model_samples_obj)) {
                                                stop("No model samples object ('model_samples_obj') available")
                                              }
                                              if (is.null(model_summary_and_trace_obj)) {
                                                stop("No model summary + trace object ('model_summary_and_trace_obj') available")
                                              }
                                              ##
                                              if (debugging) {
                                                cat("Summary object structure:\n")
                                                str(self$model_summary_and_trace_obj)
                                              }
                                              
                                              sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( Stan_model_file_name = Stan_model_file_name,
                                                                                             Stan_mod_samples = Stan_mod_samples,
                                                                                             df_true = df_true,
                                                                                             conf_region_colour = conf_region_colour,
                                                                                             pred_region_colour = pred_region_colour)
                                              ##
                                              df_fitted <- sROC_plot_outs$df_fitted
                                              plot_1    <- sROC_plot_outs$plot_1
                                              plot_2    <- sROC_plot_outs$plot_2
                                    
                                              out_list <- list(plot_1 = plot_1, 
                                                               plot_2 = plot_2, 
                                                               df_fitted = df_fitted)
                                              
                                              return(out_list)
                                    }
                                    

                
                        
                      )
)




























































