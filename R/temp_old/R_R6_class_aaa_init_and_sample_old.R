
require(R6)

#' @title MetaOrdDTA Model Class
#' @description R6 Class for MetaOrdDTA model initialization and sampling.
#'
#' @details 
#' This class handles:
#' * Model initialization and compilation
#' * MCMC sampling with NUTS-HMC (using Stan).
#' * Parameter updates and diagnostics
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
#' @section EXAMPLE - Typical workflow:
#'  \preformatted{
#'  #### NOTE: Please see the ".R" example files for full details & examples, 
#'  #### as this is NOT a complete working example, but is provided just to show 
#'  #### what order you should call things in. 
#'  ####  -----------  Compile + initialise the model using "MetaOrd_model$new(...)"  ---------------- 
#'    require(MetaOrdDTA)
#'    ## NOTE: The "model_args_list" argument is only needed for BUILT-IN (not Stan) models.
#'    ##
#'    ## --- INITIALISE MODEL:  --------------------------------------------------------------- 
#'    ##
#'    model_obj <- MetaOrdDTA::MetaOrd_model$new(   
#'                    Model_type = "meta_ord",
#'                    x = x, ## data
#'                    model_args_list = model_args_list,
#'                    n_chains = n_chains
#'                    init_lists_per_chain = init_lists_per_chain)
#'    ##
#'    ## --- SAMPLE MODEL: --------------------------------------------------------------------
#'    ##
#'    model_samples <-  model_obj$sample(
#'                    seed = 123,
#'                    n_burnin = 500,
#'                    n_iter = 1000,
#'                    x = x,
#'                    n_chains = 16,
#'                    init_lists_per_chain = init_lists_per_chain,
#'                    model_args_list = model_args_list)
#'    ##
#'    ## --- EXTRACT MODEL RESULTS SUMMARY + DIAGNOSTICS ---------------------------------------  
#'    ##         
#'    model_fit <- model_samples$summary(
#'                    save_log_lik_trace = FALSE, 
#'                    compute_nested_rhat = FALSE,
#'                    compute_transformed_parameters = FALSE) #'    
#'    
#'    ## extract # divergences + % of sampling iterations which have divergences
#'    model_fit$get_divergences()
#'    
#'    #### --- TRACE PLOTS  ----------------------------------------------------------------------
#'    ## trace_plots_all <- model_samples$plot_traces() # if want the trace for all parameters 
#'    trace_plots <- model_fit$plot_traces(
#'                    params = c("beta", "Omega", "p"), 
#'                    batch_size = 12)
#'    
#'    ## you can extract parameters by doing: "trace$param_name()". 
#'    ## For example:
#'    ## display each panel for beta and Omega ("batch_size" controls the # of plots per panel)
#'    trace_plots$beta[[1]] # 1st (and only) panel
#'    
#'    
#'    #### --- POSTERIOR DENSITY PLOTS -----------------------------------------------------------
#'    ## density_plots_all <- model_samples$plot_densities() # if want the densities 
#'    ## for all parameters.
#'    ## Let's plot the densities for: sensitivity, specificity, and prevalence 
#'    density_plots <- model_fit$plot_densities(
#'                      params = c("Se_bin", "Sp_bin", "p"), 
#'                      batch_size = 12)
#'    
#'    ## you can extract parameters by doing: "trace$param_name()". 
#'    ## For example:
#'    ## display each panel for beta and Omega ("batch_size" controls the # of plots per panel)
#'    density_plots$Se[[1]] # Se - 1st (and only) panel
#'    density_plots$Sp[[1]] # Sp - 1st (and only) panel
#'    density_plots$p[[1]] # p (prevelance) - 1st (and only) panel
#'    
#'    
#'    ###### --- EXTRACT PARAMETER SUMMARY TIBBLES: ----------------------------------------------
#'    ## The "model_summary" object (created using the "$summary()" method) contains
#'    ##  many useful objects. 
#'    ## For example:
#'    require(dplyr)
#'    ## nice summary tibble for main parameters, includes ESS/Rhat, etc:
#'    model_fit$get_summary_main() %>% print(n = 50) 
#'    ## nice summary tibble for transformed parameters, includes ESS/Rhat, etc:
#'    model_fit$get_summary_transformed() %>% print(n = 150) 
#'    ## nice summary tibble for generated quantities, includes 
#'    ## ESS/Rhat, etc (for LC-MVP this includes Se/Sp/prevalence):
#'    model_fit$get_summary_generated_quantities () %>% print(n = 150) 
#'    
#' }
#' 
#' 
#'
#' @export
MetaOrd_model <- R6Class("MetaOrd_model",
                    
            public = list(
                  ##
                  ## -------------------------  Store all "important" parameters as class members / define them first:
                  ##
                  #'@field model_prep_obj Model initialization object.
                  model_prep_obj = NULL,
                  ##
                  #'@field model_samples_obj Model samples object.
                  model_samples_obj = NULL, 
                  ##
                  #'@field model_summary_and_trace_obj Model summary + trace object.
                  model_summary_and_trace_obj = NULL,
                  ##
                  ## ----- "x" (data):
                  #'@field x The dataset. This should be an R list with 2 items, where the first element in the list is a matrix with n_rows = n_studies and n_col = n_thr, where "n_studies" 
                  #' is equal to the total number of studies in the meta-analysis, and "n_thr" is the total number of observed thresholds across all studies. 
                  #' The entries in this matrix should be equal to the total number of individuals who score less than or equal to threshold k (where k corresponds to the row number of the 
                  #' matrix \code{x[[1]]}, since each row represents a threshold). 
                  #' Similarly, the second element of the R list (i.e., \code{x[[2]]}) is the same except that it's for the diseased individuals.
                  #' For example, entry \code{x[[1]][2, 3] would correspond to the total number of non-diseased individuals in study 2 who score less than or equal to threshold 3.}
                  x = NULL,
                  ## ----- "test_type":
                  #'@field test_type Type of test (and hence the type of model to fit) - can be either "continuous" (or simply "cts") or "ordinal" (or simply "ord"). The former uses the
                  #' model from Jones et al and assumes, that the test is continuous (e.g., biokarkers) and the latter assumes the data is ordinal (e.g., questionnaires). 
                  #' The latter is generally more flexible so if you are unsure which model type to use, we suggest using "ordinal". 
                  #' Furthermore, if performing network-meta-analysis (NMA) and you have a mix of continuous and ordinal tests, we also suggest using "ordinal". 
                  #' The default is "ordinal".
                  test_type = NULL, 
                  ## ----- "network":
                  #' @field network Whether to perform meta-analysis (MA) or network-meta-analysis (NMA). The default will depend on the structure of the data (i.e. on \code{x}).
                  network = NULL,
                  ## ----- "cts":
                  #' @field cts Whether to use continuous model (i.e., Jones et al-based) or ordinal. Default is ordinal.
                  cts = NULL,
                  ## ----- "prior_only":
                  #' @field prior_only Whether to run prior-only (for prior predictive checking) or not. Default is FALSE. 
                  prior_only = NULL,
                  ## ----- "priors":
                  #' @field priors  An R list of priors. Default will depend on model.
                  priors = NULL,
                  ## ----- "model_parameterisation":
                  #' @field model_parameterisation Advanced option.
                  model_parameterisation = NULL,
                  ## ----- "random_thresholds":
                  #' @field random_thresholds Advanced option.
                  random_thresholds = NULL,
                  ## ----- "Dirichlet_random_effects_type":
                  #' @field Dirichlet_random_effects_type Advanced option.
                  Dirichlet_random_effects_type = NULL,
                  ## ----- "box_cox":
                  #' @field box_cox Advanced option.
                  box_cox = NULL,
                  ## ----- "softplus":
                  #' @field softplus Advanced option.
                  softplus = NULL,
                  ## ---- "n_chains":
                  #'@field n_chains The total number of burn-in chains. The default is Min(8, n_cores), where n_cores is the number of cores on the CPU. 
                  n_chains = NULL,
                  ## ---- "init_lists_per_chain":
                  #'@field init_lists_per_chain A list of dimension n_chains, where each element of the list is another list which contains the #
                  #' initial values for the chain. Note that this is the same format for initial values that Stan (both rstan and cmdstanr) uses. 
                  init_lists_per_chain = NULL,
                  ##
                  ## ---- Other class members:
                  ##
                  #'@field model_args_list  List containing model-specific arguments.
                  model_args_list = NULL,
                  ##
                  #'@field Stan_data_list
                  Stan_data_list = NULL,
                  #'@field Stan_model_file_name
                  Stan_model_file_name = NULL,
                  #'@field stan_model_obj
                  stan_model_obj = NULL,
                  #'@field stan_functions_dir
                  stan_functions_dir = NULL,
                  ##
                  #'@field n_studies
                  n_studies = NULL,
                  #'@field n_cat
                  n_cat = NULL,
                  #'@field n_thr
                  n_thr = NULL,
                  ## ---------- constructor - initialize using the prep_data_and_model fn (this wraps the prep_data_and_model function with $new()) - store all important parameters:
                  #'@description Create a new "MetaOrd_model" object.
                  #'@param x The dataset. See class documentation for details.
                  ##
                  #'@param test_type Type of test. See class documentation for details.
                  #'@param network See class documentation for details.
                  #'@param cts Whether to use a continuous ("cts") model or whether to use an ordinal ("ord") model.
                  ##
                  #'@param prior_only See class documentation for details.
                  ##
                  #'@param priors See class documentation for details.
                  ##
                  #'@param model_parameterisation See class documentation for details.
                  #'@param random_thresholds See class documentation for details.
                  #'@param Dirichlet_random_effects_type See class documentation for details.
                  ##
                  #'@param box_cox See class documentation for details.
                  ##
                  #'@param softplus See class documentation for details.
                  ##
                  #'@param n_chains Number of chains used for burnin. See class documentation for details.
                  #'@param init_lists_per_chain List of initial values for each chain. See class documentation for details.
                  #'@return Returns self$model_initial_obj, an object generated from the "MetaOrdDTA::prep_data_and_model" function.
                  ##
                  initialize = function(  x,
                                          ##
                                          test_type = NULL,
                                          network = NULL,
                                          cts = NULL,
                                          ##
                                          prior_only = NULL,
                                          ##
                                          priors = NULL,
                                          ##
                                          ## "Advanced" options for ordinal models ONLY:
                                          ##
                                          model_parameterisation = NULL,
                                          random_thresholds = NULL,
                                          Dirichlet_random_effects_type = NULL,
                                          ##
                                          ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                          ##
                                          box_cox = NULL,
                                          ##
                                          ## "Advanced" options for ALL models:
                                          ##
                                          softplus = NULL,
                                          ##
                                          n_chains = NULL,
                                          init_lists_per_chain = NULL
                                          ) {
                    
                            ##
                            ## Store important parameters as class members:
                            ##
                            self$x <- x
                            ## 
                            ##
                            self$test_type  <- test_type
                            self$network    <- network
                            self$cts <- cts
                            ##
                            self$prior_only <- prior_only
                            ##
                            self$priors <- priors
                            ##
                            self$model_parameterisation <- model_parameterisation
                            self$random_thresholds <- random_thresholds
                            self$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
                            ##
                            self$box_cox <- box_cox
                            ##
                            self$softplus <- softplus
                            ##
                            self$n_chains <- if_null_then_set_to(n_chains, round(parallel::detectCores()/2))
                            self$init_lists_per_chain <- init_lists_per_chain
                            ##
                            ## Other class members:
                            ##
                            self$model_args_list <- NULL
                            ##
                            self$Stan_data_list <- NULL
                            self$Stan_model_file_name <- NULL
                            self$stan_model_obj <- NULL
                            self$stan_functions_dir <- NULL
                            ##
                            self$n_studies <- NULL
                            self$n_cat <- NULL
                            self$n_thr <- NULL
                            ##
                            ## -----------  call initialising fn's: ------------------------------------------------------------------------------------------------------------
                            ##
                            self$model_prep_obj <-                MetaOrdDTA:::prep_data_and_model(   x = self$x,
                                                                                    ##
                                                                                    test_type = self$test_type,
                                                                                    network = self$network,
                                                                                    cts = self$cts,
                                                                                    ##
                                                                                    prior_only = self$prior_only,
                                                                                    ##
                                                                                    model_args_list = self$model_args_list,
                                                                                    ##
                                                                                    priors = self$priors,
                                                                                    ##
                                                                                    ## "Advanced" options for ordinal models ONLY:
                                                                                    ##
                                                                                    model_parameterisation = self$model_parameterisation,
                                                                                    random_thresholds = self$random_thresholds,
                                                                                    Dirichlet_random_effects_type = self$Dirichlet_random_effects_type,
                                                                                    ##
                                                                                    ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                                                                    ##
                                                                                    box_cox = self$box_cox,
                                                                                    ##
                                                                                    ## "Advanced" options for ALL models:
                                                                                    ##
                                                                                    softplus = self$softplus,
                                                                                    ##
                                                                                    ## Other stuff:
                                                                                    ##
                                                                                    init_lists_per_chain = self$init_lists_per_chain,
                                                                                    n_chains = self$n_chains)
          
                            ##
                            ## Make and/or update params:
                            ##
                            {
                                self$model_args_list <- self$model_prep_obj$model_args_list
                                ##
                                self$test_type <- self$model_prep_obj$test_type
                                self$network   <- self$model_prep_obj$network
                                self$cts       <- self$model_prep_obj$cts
                                ##
                                self$prior_only <- self$model_prep_obj$prior_only
                                ##
                                self$priors <- self$model_prep_obj$priors
                                ##
                                ## "Advanced" options for ordinal models ONLY:
                                ##
                                self$model_parameterisation        <- self$model_prep_obj$model_parameterisation
                                self$random_thresholds             <- self$model_prep_obj$random_thresholds
                                self$Dirichlet_random_effects_type <- self$model_prep_obj$Dirichlet_random_effects_type
                                ##
                                ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                ##
                                self$box_cox  <- self$model_prep_obj$box_cox
                                ##
                                ## "Advanced" options for ALL models:
                                ##
                                self$softplus <- self$model_prep_obj$softplus
                                ##
                                ## Other stuff:
                                ##
                                self$init_lists_per_chain <- self$model_prep_obj$init_lists_per_chain
                                self$n_chains             <- self$model_prep_obj$n_chains
                                ##
                                self$Stan_data_list       <- self$model_prep_obj$Stan_data_list
                                self$Stan_model_file_name <- self$model_prep_obj$Stan_model_file_name
                                self$stan_model_obj       <- self$model_prep_obj$stan_model_obj
                                self$stan_functions_dir   <- self$model_prep_obj$stan_functions_dir
                                ##
                                self$n_studies <- self$model_prep_obj$n_studies
                                self$n_cat     <- self$n_cat
                                self$n_thr     <- self$model_prep_obj$n_thr
                            }
                            
                            return(self)
                            
                            
                   },
                  ##  
                  ## --------  wrap the sample fn -------------------------------------------------------------------------------------------------------------------------------
                  ##
                  #'@description
                  #'Sample from the model
                  #'@param init_lists_per_chain List of initial values for each chain. See class documentation for details.
                  #'@param model_args_list List of model arguments. See class documentation for details.
                  #'@param Stan_data_list List of Stan data (optional). See class documentation for details.
                  #'@param x The dataset. See class documentation for details.
                  #'@param n_chains Number of chains used for sampling (i.e., post-burnin). Note that this must match the length of "init_lists_per_chain".
                  #'@param n_superchains  Number of superchains for nested R-hat (nR-hat; see Margossian et al, 2023) computation. 
                  #'@param seed Random MCMC seed.
                  #'@param n_burnin The number of burnin iterations. Default is 500. 
                  #'@param n_iter The number of sampling (i.e., post-burnin) iterations. Default is 1000. 
                  #'@param adapt_delta The Metropolis-Hastings target acceptance rate. Default is 0.80. If there are divergences, sometimes increasing this can help, at the cost 
                  #'of efficiency (since this will decrease the step-size (epsilon) and increase the number of leapfrog steps (L) per iteration).
                  #'@param max_treedepth Proportional to the maximum number of leapfrog steps. The default is 10.
                  #'@param metric_shape The shape of the metric to use for the main parameters, which is adapted during the burnin period. Can be either \code{"diag"} or 
                  #'\code{"dense"}. The default is \code{"diag"} (i.e. a diagonal metric).
                  #'@return Returns self invisibly, allowing for method chaining of this classes (MetaOrd_model) methods. E.g.: model$sample(...)$summary(...). 
                  ##
                  sample = function(  x = self$x, ## data is allowed to be updated without re-initialising (as long as "network" and "cts" are the same)
                                      ##
                                      prior_only = self$prior_only,
                                      ##
                                      priors = self$priors,
                                      ##
                                      ## "Advanced" options for ordinal models ONLY:
                                      ##
                                      model_parameterisation        = self$model_parameterisation,
                                      random_thresholds             = self$random_thresholds,
                                      Dirichlet_random_effects_type = self$Dirichlet_random_effects_type,
                                      ##
                                      ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                      ##
                                      box_cox = self$box_cox,
                                      ##
                                      ## "Advanced" options for ALL models:
                                      ##
                                      softplus = self$softplus,
                                      ##
                                      ## Other stuff:
                                      ##
                                      init_lists_per_chain = self$init_lists_per_chain,
                                      n_chains = self$n_chains,
                                      ##
                                      ## MCMC / Stan params:
                                      ##
                                      n_superchains = NULL,
                                      seed          = NULL,
                                      n_burnin      = 500,
                                      n_iter        = 1000,
                                      adapt_delta   = 0.80,
                                      max_treedepth = 10,
                                      metric_shape  = "diag_e"
                                      ) {
                              ##
                              ## ---- Validate initialization:
                              ##
                              if (is.null(self$model_prep_obj)) {
                                stop("Model was not properly initialize")
                              }
      
                              seed <- if_null_then_set_to(seed, 123) 
                              n_superchains <- if_null_then_set_to(n_superchains, n_chains)
      
                              ##
                              ## ---- Extract params which can't be updated (i.e. params which require "$new()" again):
                              ##
                              test_type <- self$test_type
                              network   <- self$network
                              cts       <- self$cts
                              
                              ##
                              ## ---- Update params (if any changed):
                              ##
                              params_same <- 1
                              ##
                              {
                                  ##
                                  ## Update class members if new values provided:
                                  ##
                                  if (!identical(self$x, x))  { self$x <- x ; params_same <- 0 }
                                  ##
                                  if (!identical(self$prior_only, prior_only)) { self$prior_only <- prior_only ; params_same <- 0 }
                                  if (!identical(self$priors, priors))         { self$priors     <- priors     ; params_same <- 0 }
                                  ##
                                  if (!identical(self$model_parameterisation, model_parameterisation))               { self$model_parameterisation        <- model_parameterisation        ; params_same <- 0 }
                                  if (!identical(self$random_thresholds, random_thresholds))                         { self$random_thresholds             <- random_thresholds             ; params_same <- 0 }
                                  if (!identical(self$Dirichlet_random_effects_type, Dirichlet_random_effects_type)) { self$Dirichlet_random_effects_type <- Dirichlet_random_effects_type ; params_same <- 0 }
                                  ##
                                  if (!identical(self$box_cox, box_cox))   { self$box_cox  <- box_cox   ; params_same <- 0 }
                                  ##
                                  if (!identical(self$softplus, softplus)) { self$softplus <- softplus ; params_same <- 0 }
                                  ##
                                  if (!identical(self$priors, priors))     { self$priors   <- priors     ; params_same <- 0 }
                                  ##
                                  #### if (!identical(self$model_args_list, model_args_list))  { self$model_args_list <- model_args_list ; params_same <- 0 }
                                  ##
                                  if (!identical(self$n_chains, n_chains))                          { self$n_chains <- n_chains ; params_same <- 0 }
                                  if (!identical(self$init_lists_per_chain, init_lists_per_chain))  { self$init_lists_per_chain <- init_lists_per_chain ; params_same <- 0 }
                                  
                              }
                          
                              # then update model if any of needed parameters changed
                              if (params_same == 0) {
                                  ##
                                  ## ---------- call "update_model" fn  -----------------------------------------------------------------------------------------
                                  ##
                                  #     self$model_prep_obj <-   MetaOrdDTA:::update_model(    test_type = test_type,
                                  #                                                       model_prep_obj = model_prep_obj,
                                  #                                                       y = self$y,
                                  #                                                       n_chains = self$n_chains
                                  #                                                       init_lists_per_chain = self$init_lists_per_chain,
                                  #                                                       model_args_list = self$model_args_list,
                                  #                                                       Stan_data_list =  self$Stan_data_list)
  
                              }
                              #
                              if (!is.null(adapt_delta) && (adapt_delta <= 0 || adapt_delta >= 1)) {
                                stop("adapt_delta must be between 0 and 1")
                              }
                              
                              ##
                              ## -----------  call "initialise_update_run_model" fn ----------------------------------------------------------------------------
                              ##
                              self$model_samples_obj <-       MetaOrdDTA:::initialise_update_run_model(      debugging = FALSE,
                                                                                     ##
                                                                                     x = self$x,
                                                                                     ##
                                                                                     test_type = self$test_type,
                                                                                     network = self$network,
                                                                                     cts = self$cts,
                                                                                     ##
                                                                                     prior_only = self$prior_only,
                                                                                     ##
                                                                                     model_args_list = self$model_args_list, ## BOOKMARK /  DELETE
                                                                                     ##
                                                                                     priors = self$priors,
                                                                                     ##
                                                                                     ## "Advanced" options for ordinal models ONLY:
                                                                                     ##
                                                                                     model_parameterisation = self$model_parameterisation,
                                                                                     random_thresholds = self$random_thresholds,
                                                                                     Dirichlet_random_effects_type = self$Dirichlet_random_effects_type,
                                                                                     ##
                                                                                     ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                                                                     ##
                                                                                     box_cox = self$box_cox,
                                                                                     ##
                                                                                     ## "Advanced" options for ALL models:
                                                                                     ##
                                                                                     softplus = self$softplus,
                                                                                     ##
                                                                                     ## Other stuff:
                                                                                     ##
                                                                                     init_lists_per_chain = self$init_lists_per_chain,
                                                                                     n_chains = self$n_chains,
                                                                                     ##
                                                                                     ## MCMC:
                                                                                     ##
                                                                                     seed = seed,
                                                                                     n_iter = n_iter,
                                                                                     n_burnin = n_burnin,
                                                                                     n_superchains = n_superchains,
                                                                                     adapt_delta = adapt_delta,
                                                                                     max_treedepth = max_treedepth,
                                                                                     metric_shape = metric_shape)
         
                                                                                          
                                                                                   
                              return(self)
                        
                                      
                    
                  },
                  ##        
                  ## --------  wrap the create_summary_and_traces fn + call the "BayesMVP::MVP_class_plot" R6 class  ----------------------------------------------------------
                  ##
                  #'@description
                  #'Create and compute summary statistics, traces and model diagnostics. 
                  #'@param compute_main_params Whether to compute the main parameter summaries. Default is TRUE.
                  #'@param compute_transformed_parameters Whether to compute transformed parameter summaries. Default is TRUE. For Stan models, this will be for 
                  #'all of the parameters defined in the "transformed parameters" block, EXCEPT for the (transformed) nuisance parameters and log_lik (see "save_log_lik_trace"
                  #'for more information on log_lik). 
                  #'@param compute_generated_quantities Whether to compute the summaries for generated quantities. Default is TRUE.
                  #'@param save_log_lik_trace Whether to save the log-likelihood (log_lik) trace. Default is FALSE. For Stan models, this will only work 
                  #'if there is a "log_lik" parameter defined in the "transformed parameters" model block.
                  #'@param compute_nested_rhat Whether to compute the nested rhat diagnostic (nR-hat) (Margossian et al, 2023). This is useful when
                  #'running many (usually short) chains. Also see "n_superchains" argument. 
                  #'@param n_superchains The number of superchains to use for the computation of nR-hat. 
                  #'Only relevant if \code{compute_nested_rhat = TRUE}.
                  #'@param save_trace_tibbles Whether to save the trace as tibble dataframes as well as 3D arrays. Default is FALSE. 
                  #'@return Returns a new MetaOrd_plot_and_diagnose object (from the "MetaOrd_plot_and_diagnose" R6 class) for creating MCMC diagnostics and plots.
                  ##
                  summary = function(       compute_main_params = TRUE,
                                            compute_transformed_parameters = TRUE,
                                            compute_generated_quantities = TRUE,
                                            save_log_lik_trace = TRUE,
                                            compute_nested_rhat = NULL,
                                            n_superchains = NULL,
                                            save_trace_tibbles = FALSE
                                            ) {
                    
                        ## validate initialization:
                        if (is.null(self$model_prep_obj)) {
                          stop("Model was not properly initialize")
                        }
                    
                        ##
                        # initial_object <- self$initial_object
                        ##
                        Stan_model_sample_output <- self$model_samples_obj$full_outputs$Stan_model_sample_output
                        stan_param_names_list    <- self$model_samples_obj$stan_param_names_list
                        ##
                        ## create model fit object (includes model summary tables + traces + divergence info) by calling "BayesMVP::create_summary_and_traces" -----------------
                        ##
                        self$model_summary_and_trace_obj <-           MetaOrdDTA:::create_summary_and_traces(    package = "MetaOrdDTA",
                                                                                                      use_bridgestan = FALSE,
                                                                                                      use_BayesMVP_for_faster_summaries = MetaOrdDTA:::check_if_BayesMVP_R_pkg_installed(TRUE, FALSE),
                                                                                                      ##
                                                                                                      Stan_model_sample_output = Stan_model_sample_output,
                                                                                                      stan_param_names_list = stan_param_names_list,
                                                                                                      ##
                                                                                                      model_results  = NULL,
                                                                                                      init_object = NULL,
                                                                                                      ##
                                                                                                      n_nuisance = 0,
                                                                                                      ##
                                                                                                      compute_main_params = compute_main_params,
                                                                                                      compute_transformed_parameters = compute_transformed_parameters,
                                                                                                      compute_generated_quantities = compute_generated_quantities,
                                                                                                      ##
                                                                                                      save_log_lik_trace = save_log_lik_trace,
                                                                                                      ##
                                                                                                      compute_nested_rhat = compute_nested_rhat,
                                                                                                      n_superchains = n_superchains,
                                                                                                      save_trace_tibbles = save_trace_tibbles)
                        
                        ##
                        ## return the plotting class instance with the summary:
                        ##
                        MetaOrd_class_plot_object <- MetaOrdDTA::MetaOrd_plot_and_diagnose$new( model_samples_obj = self$model_samples_obj,
                                                                                                model_summary_and_trace_obj = self$model_summary_and_trace_obj)
                                                                                      
                        
                        return(MetaOrd_class_plot_object)
                    
                  }
                          
            )
)









