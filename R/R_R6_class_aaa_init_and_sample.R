
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
                    debugging = NULL,
                    ##
                    ## ----  Core / internal objects (NOT user-inputted)
                    ##
                    outs_data = list(
                      stan_data_list = NULL,
                      n_tests = NULL,
                      n_studies = NULL,
                      n_thr = NULL,
                      n_cat = NULL
                    ),
                    ##
                    outs_stan_model_name = list( 
                      stan_model_file_name = NULL
                    ),
                    ##
                    outs_stan_compile = list( 
                        stan_model_obj = NULL, 
                        # stan_model_file_name = NULL,
                        stan_model_file_path = NULL,
                        ##
                        pkg_root_directory = NULL,
                        stan_models_directory = NULL,
                        stan_functions_directory = NULL,
                        ##
                        stan_MA_directory = NULL,
                        stan_MA_prior_directory = NULL,
                        ##
                        stan_NMA_directory = NULL,
                        stan_NMA_directory = NULL
                    ),
                    ##
                    outs_stan_init = list( 
                      inits_unconstrained_vec_per_chain = NULL,
                      stan_param_names_list = NULL,
                      stan_param_names_main = NULL,
                      stan_init_pseudo_sampling_outs = NULL,
                      stan_model_obj = NULL,
                      json_file_path = NULL,
                      stan_model_file_path = NULL
                    ),
                    ##
                    outs_stan_sampling = list( 
                      stan_mod_samples = NULL,
                      time_total = NULL
                    ),
                    ##
                    ## ---- Main "internal_obj" list:
                    ##
                    internal_obj = list(
                        outs_data = NULL,
                        ##
                        outs_stan_model_name = NULL,
                        outs_stan_compile = NULL,
                        outs_stan_init = NULL,
                        outs_stan_sampling = NULL,
                        ##
                        HMC_info = NULL, ## new
                        efficiency_info = NULL, ## new
                        summaries = NULL, ## new
                        traces = NULL ## new
                    ),
                    ##
                    ## ---- Data (user-inputted - MANDATORY argument)
                    ##
                    x = NULL,
                    ##
                    indicator_index_test_in_study = NULL,
                    ##
                    ## ---- Basic modelling options (user-inputted element-by-element)
                    ##
                    basic_model_options = list(
                      network = NULL,
                      cts = NULL,
                      prior_only = NULL
                    ),
                    #' @field network Whether to perform meta-analysis (MA) or network-meta-analysis (NMA). The default will depend on the structure of the data (i.e. on \code{x}).
                    #' @field cts Whether to use continuous model (i.e., Jones et al-based) or ordinal. Default is ordinal. Type of test (and hence the type of model to fit) -
                    #' can be either "continuous" (or simply "cts") or "ordinal" (or simply "ord"). The former uses the
                    #' model from Jones et al and assumes, that the test is continuous (e.g., biokarkers) and the latter assumes the data is ordinal (e.g., questionnaires). 
                    #' The latter is generally more flexible so if you are unsure which model type to use, we suggest using "ordinal". 
                    #' Furthermore, if performing network-meta-analysis (NMA) and you have a mix of continuous and ordinal tests, we also suggest using "ordinal". 
                    #' The default is "ordinal".
                    #' @field prior_only Whether to run prior-only (for prior predictive checking) or not. Default is FALSE. 
                    ##
                    ## ---- "advanced" user-inputted modelling options  (OPTIONAL - user-inputted element-by-element)
                    ##
                    advanced_model_options = list(
                      model_parameterisation = NULL,
                      random_thresholds = NULL,
                      Dirichlet_random_effects_type = NULL,
                      box_cox = NULL,
                      softplus = NULL
                    ),
                    ##
                    ## ---- priors  (OPTIONAL - user-inputted element-by-element)
                    ##
                    priors = NULL,
                    #' @field priors  An R list of priors. Default will depend on model.
                    ##
                    ## ---- initial values for each chain (OPTIONAL - user-inputted element-by-element)
                    ##
                    init_lists_per_chain = NULL,
                    #'@field init_lists_per_chain A list of dimension n_chains, where each element of the list is another list which contains the #
                    #' initial values for the chain. Note that this is the same format for initial values that Stan (both rstan and cmdstanr) uses. 
                    ##
                    ## ----  (OPTIONAL - user-inputted element-by-element)
                    ##
                    MCMC_params = list(
                      seed = NULL,
                      n_superchains = NULL,
                      n_chains = NULL,
                      n_iter = NULL,
                      n_burnin = NULL,
                      adapt_delta = NULL,
                      max_treedepth = NULL,
                      metric_shape = NULL
                    ),
                    ##
                    ## ---- Other class members:
                    ##
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
                    initialize = function(  debugging = FALSE,
                                            ##
                                            x = self$x,
                                            ##
                                            indicator_index_test_in_study = self$indicator_index_test_in_study,
                                            ##
                                            network = NULL,
                                            cts = NULL,
                                            prior_only = NULL,
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
                                            priors = NULL,
                                            ##
                                            n_chains = NULL,
                                            init_lists_per_chain = self$init_lists_per_chain
                                            ) {
                      
                              # message(("aaa_class_hello_1"))
                              ##
                              ## ---- Store important parameters as class members:
                              ##
                              self$debugging <- debugging
                              ##
                              ## ---- Data:
                              ##
                              self$x <- x
                              ##
                              self$indicator_index_test_in_study <- indicator_index_test_in_study
                              ##
                              ## ---- "basic_model_options" list:
                              ##
                              self$basic_model_options$network <- network
                              self$basic_model_options$cts <- cts
                              self$basic_model_options$prior_only <- prior_only
                              ##
                              ## ---- "advanced_model_options" list:
                              ##
                              self$advanced_model_options$model_parameterisation        <- model_parameterisation
                              self$advanced_model_options$random_thresholds             <- random_thresholds
                              self$advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
                              self$advanced_model_options$box_cox  <- box_cox
                              self$advanced_model_options$softplus <- softplus
                              ##
                              ## ---- priors:
                              ##
                              self$priors <- priors
                              ##
                              ## ---- initial values (per chain):
                              ##
                              self$init_lists_per_chain <- init_lists_per_chain
                              ##
                              ## ---- MCMC params (only n_chains can be inputted for initialisation):
                              ##
                              self$MCMC_params$n_chains <- n_chains
                              ##
                              # message(("aaa_class_hello_2"))
                              ##
                              ## -----------  call initialising fn's: ------------------------------------------------------------------------------------------------------------
                              ##
                              local_prep_data_and_model_outs <-  MetaOrdDTA:::prep_data_and_model(  
                                                                                    debugging = self$debugging,
                                                                                    ##
                                                                                    x = self$x,
                                                                                    indicator_index_test_in_study = self$indicator_index_test_in_study,
                                                                                    ##
                                                                                    internal_obj = self$internal_obj,
                                                                                    ##
                                                                                    basic_model_options = self$basic_model_options,
                                                                                    advanced_model_options = self$advanced_model_options,
                                                                                    MCMC_params = self$MCMC_params,
                                                                                    ##
                                                                                    priors = self$priors,
                                                                                    init_lists_per_chain = self$init_lists_per_chain)
                              
                              
                              ## output of MetaOrdDTA:::prep_data_and_model:
                              # return(list(  full_outputs = full_outputs,
                              #               ##
                              #               internal_obj = internal_obj,
                              #               basic_model_options = basic_model_options,
                              #               advanced_model_options = advanced_model_options,
                              #               MCMC_params = MCMC_params,
                              #               ##
                              #               priors = priors,
                              #               ##
                              #               init_lists_per_chain = init_lists_per_chain))
                                                                                                 
                              
                              ##
                              ## ---- Update R grouped lists:
                              ##
                              self$internal_obj           <- local_prep_data_and_model_outs$internal_obj
                              ##
                              self$basic_model_options    <- local_prep_data_and_model_outs$basic_model_options
                              self$advanced_model_options <- local_prep_data_and_model_outs$advanced_model_options
                              self$MCMC_params            <- local_prep_data_and_model_outs$MCMC_params
                              ##
                              self$priors <- local_prep_data_and_model_outs$priors
                              ##
                              self$init_lists_per_chain <- local_prep_data_and_model_outs$init_lists_per_chain
                              
                              return(self)
                              
                              
                     },
                    ##  
                    ## --------  wrap the sample fn -------------------------------------------------------------------------------------------------------------------------------
                    ##
                    #'@description
                    #'Sample from the model
                    #'@param init_lists_per_chain List of initial values for each chain. See class documentation for details.
                    #'@param model_args_list List of model arguments. See class documentation for details.
                    #'@param stan_data_list List of Stan data (optional). See class documentation for details.
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
                    sample = function(  debugging = self$debugging,
                                        ##
                                        x = self$x, ## data is allowed to be updated without re-initialising (as long as "network" and "cts" are the same)
                                        indicator_index_test_in_study = self$indicator_index_test_in_study,
                                        ##
                                        cts = self$basic_model_options$cts,
                                        network = self$basic_model_options$network,
                                        prior_only = self$basic_model_options$prior_only,
                                        ##
                                        model_parameterisation        = self$advanced_model_options$model_parameterisation,
                                        random_thresholds             = self$advanced_model_options$random_thresholds,
                                        Dirichlet_random_effects_type = self$advanced_model_options$Dirichlet_random_effects_type,
                                        box_cox  = self$advanced_model_options$box_cox,
                                        softplus = self$advanced_model_options$softplus,
                                        ##
                                        priors = self$priors,
                                        ##
                                        init_lists_per_chain = self$init_lists_per_chain,
                                        ##
                                        ## MCMC / Stan params:
                                        ##
                                        seed          = NULL,
                                        n_superchains = NULL,
                                        n_chains = self$MCMC_params$n_chains,
                                        n_burnin      = 500,
                                        n_iter        = 1000,
                                        adapt_delta   = 0.80,
                                        max_treedepth = 10,
                                        metric_shape  = "diag_e"
                                        ) {
                      
                                ##
                                ## ---- Validate initialization:
                                ##
                                if (is.null(self$internal_obj)) {
                                  stop("Model was not properly initialize")
                                }
                                ##
                                # seed <- if_null_then_set_to(seed, 123) 
                                # n_superchains <- if_null_then_set_to(n_superchains, n_chains)
                                # ##
                                # ## ---- Extract params which can't be updated (i.e. params which require "$new()" again):
                                # ##
                                # network   <- self$network
                                # cts       <- self$cts
                                ##
                                ## ---- Update MCMC params:
                                ##
                                self$MCMC_params$n_superchains <- n_superchains
                                self$MCMC_params$n_chains <- n_chains
                                self$MCMC_params$seed <- seed
                                self$MCMC_params$n_burnin <- n_burnin
                                self$MCMC_params$n_iter <- n_iter
                                self$MCMC_params$adapt_delta <- adapt_delta
                                self$MCMC_params$max_treedepth <- max_treedepth
                                self$MCMC_params$metric_shape <- metric_shape
                      
                                ##
                                ## ---- Update params (if any changed):
                                ##
                                params_same <- 1
                                ##
                                {
                                    ##
                                    ## Update class members ** only if ** new values provided:
                                    ##
                                    if (!identical(self$x, x))  { self$x <- x ; params_same <- 0 }
                                    if (!identical(self$indicator_index_test_in_study, indicator_index_test_in_study))  { self$indicator_index_test_in_study <- indicator_index_test_in_study ; params_same <- 0 }
                                    ##
                                    if (!identical(self$basic_model_options$network, network)) { self$basic_model_options$network <- network ; params_same <- 0 }
                                    if (!identical(self$basic_model_options$cts, cts)) { self$basic_model_options$cts <- cts ; params_same <- 0 }
                                    if (!identical(self$basic_model_options$prior_only, prior_only)) { self$basic_model_options$prior_only <- prior_only ; params_same <- 0 }
                                    ##
                                    if (!identical(self$advanced_model_options$model_parameterisation, model_parameterisation))               { self$advanced_model_options$model_parameterisation        <- model_parameterisation        ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$random_thresholds, random_thresholds))                         { self$advanced_model_options$random_thresholds             <- random_thresholds             ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$Dirichlet_random_effects_type, Dirichlet_random_effects_type)) { self$advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$box_cox, box_cox))   { self$advanced_model_options$box_cox  <- box_cox  ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$softplus, softplus)) { self$advanced_model_options$softplus <- softplus ; params_same <- 0 }
                                    ##
                                    if (!identical(self$priors, priors))         { self$priors     <- priors     ; params_same <- 0 }
                                    ##
                                    if (!identical(self$init_lists_per_chain, init_lists_per_chain))  { self$init_lists_per_chain <- init_lists_per_chain ; params_same <- 0 }
                                    ##
                                    if (!identical(self$MCMC_params$n_chains, n_chains)) { self$MCMC_params$n_chains <- n_chains ; params_same <- 0 }
                                    
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
                                    #                                                       stan_data_list =  self$stan_data_list)
    
                                }
                                #
                                if (!is.null(adapt_delta) && (adapt_delta <= 0 || adapt_delta >= 1)) {
                                  stop("adapt_delta must be between 0 and 1")
                                }
                                
                                ##
                                ## -----------  call "initialise_update_run_model" fn ----------------------------------------------------------------------------
                                ##
                                local_model_samples_obj <-       MetaOrdDTA:::initialise_update_run_model(     
                                                                                               debugging = self$debugging,
                                                                                               ##
                                                                                               x = self$x,
                                                                                               indicator_index_test_in_study = self$indicator_index_test_in_study,
                                                                                               ##
                                                                                               internal_obj = self$internal_obj,
                                                                                               ##
                                                                                               basic_model_options = self$basic_model_options,
                                                                                               advanced_model_options = self$advanced_model_options,
                                                                                               MCMC_params = self$MCMC_params,
                                                                                               ##
                                                                                               priors = self$priors,
                                                                                               ##
                                                                                               init_lists_per_chain = self$init_lists_per_chain)

                                self$internal_obj <- local_model_samples_obj$internal_obj
                                ##
                                self$basic_model_options    <- local_model_samples_obj$basic_model_options
                                self$advanced_model_options <- local_model_samples_obj$advanced_model_options
                                self$MCMC_params            <- local_model_samples_obj$MCMC_params
                                ##
                                self$priors <- local_model_samples_obj$priors
                                ##
                                self$init_lists_per_chain <- local_model_samples_obj$init_lists_per_chain
                                                                                     
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
                    summary = function(       debugging = self$debugging,
                                              ##
                                              compute_main_params = TRUE,
                                              compute_transformed_parameters = TRUE,
                                              compute_generated_quantities = TRUE,
                                              save_log_lik_trace = TRUE,
                                              compute_nested_rhat = NULL,
                                              n_superchains = NULL,
                                              ##
                                              save_trace_tibbles = FALSE,
                                              ##
                                              use_BayesMVP_for_faster_summaries = NULL
                                              ) {
                        
                            ## validate initialization:
                            if (is.null(self$internal_obj)) {
                              stop("Model was not properly initialize")
                            }
                            ##
                            print(paste("internal_obj = "))
                            print(str(self$internal_obj))
                            print(self$internal_obj)
                            ##
                            use_BayesMVP_for_faster_summaries <- if_null_then_set_to(use_BayesMVP_for_faster_summaries, 
                                                                                     MetaOrdDTA:::check_if_BayesMVP_R_pkg_installed(TRUE, FALSE))
                            ##
                            # initial_object <- self$initial_object
                            ##
                            # stan_model_sample_output <- self$model_samples_obj$full_outputs$stan_model_sample_output ## add "stan_model_sample_output" to "internal_obj" ???
                            # stan_param_names_list    <- self$internal_obj$stan_param_names_list
                            ##
                            ## create model fit object (includes model summary tables + traces + divergence info) by calling "BayesMVP::create_summary_and_traces" -----------------
                            ##
                            local_model_summary_and_trace_obj <- MetaOrdDTA:::create_summary_and_traces(  
                                                                                  #### debugging = self$debugging,
                                                                                  ##
                                                                                  package = "MetaOrdDTA",
                                                                                  use_bridgestan = FALSE,
                                                                                  use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries,
                                                                                  ##
                                                                                  # stan_model_sample_output = stan_model_sample_output,
                                                                                  # stan_param_names_list = stan_param_names_list,
                                                                                  internal_obj = self$internal_obj,
                                                                                  MCMC_params  = self$MCMC_params,
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
                            
                            self$internal_obj$HMC_info <- local_model_summary_and_trace_obj$HMC_info
                            self$internal_obj$efficiency_info <- local_model_summary_and_trace_obj$efficiency_info
                            self$internal_obj$summaries <- local_model_summary_and_trace_obj$summaries
                            self$internal_obj$traces <- local_model_summary_and_trace_obj$traces
                            
                            return(self)
                      
                   },
                   ##
                   ## ------------- plots + MCMC & efficiency diagnostics:   ------------------------------------------------------------------------------------------------------
                   ##
                   ##
                   ## ---- Get divergence info:
                   ##
                   #'@return Returns self$internal_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                   #'self$internal_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                   get_divergences = function() {
                     return(self$internal_obj$HMC_info$HMC_diagnostic_info$divergences)
                   },
                   ##
                   ## ---- HMC diagnostic info:
                   ##
                   #'@return Returns self$internal_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                   #'self$internal_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                   get_HMC_diagnostic_info = function() { 
                     return(self$internal_obj$HMC_info$HMC_diagnostic_info) ## Extract the "HMC_info" list
                   },
                   ##
                   ## --------- Convenience method for HMC algorithm metrics:
                   ##
                   #' @return Returns a list containing HMC algorithm info (including timings and some efficiency metrics)
                   get_HMC_info = function() { 
                     return(self$internal_obj$HMC_info) ## Extract the "HMC_info" list
                   },
                   ##
                   ## ---- Summary tibble methods:
                   ##
                   #'@return Returns main parameter summaries as a tibble.
                   ##
                   get_summary_main = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_main_params)
                   },
                   #'@return Returns transformed parameter summaries as a tibble.
                   ##
                   get_summary_transformed = function() {
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters)
                   },
                   #'@return Returns generated quantities as a tibble.
                   ##
                   get_summary_generated_quantities = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_generated_quantities)
                   },
                   ##
                   ## ---- Posterior draws / 3D traces methods:
                   ##
                   #'Convenience method for ALL traces
                   #'@return Returns all traces including as tibbles (if enabled) and as 3D array's (iterations, chains, parameters). 
                   ##
                   get_all_traces_list = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_main = function() { #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_params_main)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_transformed = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_tp_wo_log_lik)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_generated_quantities = function() { #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_gq)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_log_lik = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_log_lik)
                   },
                   #'Convenience method for posterior draws (as tibbles)
                   #'@return Returns trace as tibbles. The trace gets returned as 3 seperate tibbles (one tibble for the 
                   #'main parameters, one for the transformed parameters, and one for generated quantities). 
                   ##
                   get_traces_as_tibbles = function() {  #### --------------------------------------------------- WORKING
                     list(
                       trace_as_tibble_main_params          =  self$internal_obj$traces$traces_as_tibbles$trace_params_main_tibble,
                       trace_as_tibble_transformed_params   =  self$internal_obj$traces$traces_as_tibbles$trace_transformed_params_tibble,
                       trace_as_tibble_generated_quantities =  self$internal_obj$traces$traces_as_tibbles$trace_generated_quantities_tibble
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
                     return(self$internal_obj$efficiency_info) ## Extract the "efficiency_info" list
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
                     
                     if (is.null(self$internal_obj)) {
                       stop("No summary object available")   
                     }
                     
                     ## get required efficiency info first
                     Min_ESS_main   <- self$internal_obj$efficiency_info$Min_ESS_main
                     time_burnin    <- self$internal_obj$efficiency_info$time_burnin
                     time_sampling  <- self$internal_obj$efficiency_info$time_sampling
                     time_summaries <- self$internal_obj$efficiency_info$time_summaries
                     
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
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")   
                               }
                               
                               ## get required efficiency info first:
                               Min_ESS_main   <- self$internal_obj$efficiency_info$Min_ESS_main
                               n_iter         <- self$internal_obj$efficiency_info$n_iter
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
                   #'@param param_string The parameters to generate trace plots for. This is a character vector - e.g. to plot trace plots for 
                   #'beta: param_string = c("beta"). 
                   #'The default is NULL and the trace plots for all model parameters (which have a trace array) will be plotted. 
                   #'@param batch_size The number of trace plots to display per panel. Default is 9. 
                   #'@return If no parameters specified (i.e. param_string = NULL), then this will return an object containing the trace plots for all
                   #'model parameters  which have a trace array.
                   plot_traces = function( debugging = FALSE,
                                           param_string = NULL, 
                                           condition = "exact_match",
                                           batch_size = 9
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")  
                               }
                               if (!is.null(param_string) && !is.character(param_string)) {
                                 stop("param_string must be NULL or a character vector")
                               }
                               if (!is.numeric(batch_size) || batch_size < 1) {
                                 stop("batch_size must be a positive integer")
                               }
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               # ##
                               # ## ---- Format the array for Bayesplot R package:
                               # ##
                               # draws_array <- MetaOrdDTA:::format_named_array_for_bayesplot(self$internal_obj$traces$traces_as_arrays$draws_array)
                               ##
                               if (is.null(param_string)) { # plot all param_string
                                 
                                     bayesplot::mcmc_trace(draws_array)
                                 
                               } else {   # plot specific param_string using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot 
                                 
                                     MetaOrdDTA:::plot_param_group_batched(       debugging = debugging,
                                                                                  draws_array = self$internal_obj$traces$traces_as_arrays$draws_array,
                                                                                  plot_type = "trace",
                                                                                  param_string = param_string,
                                                                                  condition = condition,
                                                                                  batch_size = batch_size)
                                                                        
                               }
                               
                   },
                   ##
                   ## --------- MCMC density plots method
                   ##
                   #'@param debugging Whether debugging mode is turned on (\code{TRUE}) or off (\code{FALSE}). The default is FALSE.
                   #'@param param_string The parameters to generate posterior density plots for. This is a character vector - e.g. to plot the 
                   #'posterior density plots for beta: param_string = c("beta"). 
                   #'The default is NULL and the posterior density plots for all model parameters (which have a trace array) will be plotted. 
                   #'@param batch_size The number of posterior density plots to display per panel. Default is 9. 
                   #'@return If no parameters specified (i.e. param_string = NULL), then this will return an object containing the posterior density
                   #'plots for all model parameters which have a trace array.
                   plot_densities = function( debugging = FALSE,
                                              param_string = NULL, 
                                              condition = "exact_match",
                                              batch_size =  9
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")  
                               }
                               if (!is.null(param_string) && !is.character(param_string)) {
                                 stop("param_string must be NULL or a character vector")
                               }
                               if (!is.numeric(batch_size) || batch_size < 1) {
                                 stop("batch_size must be a positive integer")
                               }
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               # ##
                               # ## ---- Format the array for Bayesplot R package:
                               # ##
                               # draws_array <- MetaOrdDTA:::format_named_array_for_bayesplot(self$internal_obj$traces$traces_as_arrays$draws_array)
                               ##
                               if (is.null(param_string)) { # plot all param_string
                                 
                                     bayesplot::mcmc_dens(draws_array)
                                 
                               } else {   # plot specific param_string using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot
                                 
                                     MetaOrdDTA:::plot_param_group_batched(    plot_type = "density",
                                                                               debugging = debugging,
                                                                               draws_array = self$internal_obj$traces$traces_as_arrays$draws_array,
                                                                               param_string = param_string,
                                                                               condition = condition,
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
                   plot_sROC = function( debugging = FALSE,
                                         df_true = NULL,
                                         conf_region_colour = "blue", 
                                         pred_region_colour = "blue"
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model samples object ('internal_obj') available")
                               }
                               if (is.null(self$internal_obj)) {
                                 stop("No model summary + trace object ('internal_obj') available")
                               }
                               ##
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               ##
                               sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( stan_model_file_name = self$internal_obj$outs_stan_model_name$stan_model_file_name,
                                                                              stan_mod_samples = self$internal_obj$outs_stan_sampling$stan_mod_samples,
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













