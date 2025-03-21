

 
 
##
rm(list = ls())
# ##
# .rs.restartR()  # In RStudio
# # 



#  
{
  ## First remove any possible package fragments:
  ## Find user_pkg_install_dir:
  user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
  print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
  ##
  ## Find pkg_install_path + pkg_temp_install_path:
  pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
  pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
  ##
  ## Remove any (possible) MetaOrdDTA package fragments:
  remove.packages("MetaOrdDTA")
  unlink(pkg_install_path, recursive = TRUE, force = TRUE)
  unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
}
# 

os <- .Platform$OS.type


if (os == "unix") { 
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") { 
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")



# if (.Platform$OS.type == "windows") {
#   user_root_dir <- Sys.getenv("PWD")
# } else {
#   user_root_dir <- Sys.getenv("HOME")
# }



# #### -------- INNER pkg stuff:
# ## Only works if do this first (on Linux)!! :
# source("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/MetaOrdDTA/src/R_script_load_OMP_Linux.R")
# ## Document INNER package:
# devtools::clean_dll(local_INNER_pkg_dir)
# Rcpp::compileAttributes(local_INNER_pkg_dir)
# ### devtools::document(local_INNER_pkg_dir)
# roxygen2::roxygenize(local_INNER_pkg_dir)



 
#### -------- ACTUAL (LOCAL) INSTALL:
## Document:
devtools::clean_dll(local_pkg_dir)
roxygen2::roxygenize(local_pkg_dir)
##
## Install (outer pkg):
##
devtools::clean_dll(local_pkg_dir)
devtools::install(local_pkg_dir,
                  upgrade = "never",
                  quick = TRUE)
##
## May need to restart R:
##
.rs.restartR()  # In RStudio

# ?devtools::install





##  -----------  Compile + initialise the model using "MVP_model$new(...)"  ----------------------------------
##
## ----  Prepare the data to use: ----------------------------------------------------------------------------
##
## Read in data ("x" - aggregative cumulative counts w/ "missing thresholds" in each study marked as "-1"):
##
# x <- list(x_nd, x_d)
# saveRDS(object = x, file = file.path(local_pkg_dir, "MetaOrdDTA_example_data_1.RDS"))
# x <- readRDS(file.path(local_pkg_dir, "inst", "examples", "MetaOrdDTA_example_data_1.RDS"))

x <- list(x_nd, x_d)

x_full <- x

x_subset <- list()
N_subset <- n_studies
for (c in 1:2) {
  x_subset[[c]] <- x[[c]][1:N_subset, ]
}
 


##
## ----  Initialise / select model: --------------------------------------------------
##
n_studies <- nrow(x_subset[[1]])
n_thr <- ncol(x_subset[[1]])
##
#### n_chains <- ifelse(parallel::detectCores() > 32, 64, 16)
n_chains <- 8
##
 

#### ----
{
  model_parameterisation = "Jones"
  box_cox <- TRUE
  cts <- TRUE
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
  ##
  inits_for_Xu_fixed_thr_model <- list(beta_mu = c(-1, +1),
                                       beta_SD = c(0.01, 0.01),
                                       beta_z = array(0.01, dim = c(2, n_studies)),
                                       C_raw_vec = rep(-2.0, n_thr),
                                       C = seq(from = -2, to = 2, length = n_thr))
  init_lists_per_chain <- replicate(n_chains, inits_for_Xu_fixed_thr_model, simplify = FALSE)
}


#### ----
{
  model_parameterisation = "Xu"
  ##
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "SD"
  ##
  priors <- list(prior_alpha = rep(1.0, n_thr + 1), 
                 ##
                 prior_beta_mu_mean = c(0.0, 0.0),
                 prior_beta_mu_SD = c(1.0, 1.0),
                 prior_beta_SD_mean = c(0.0, 0.0),
                 prior_beta_sd_SD = c(0.5, 0.5),
                 prior_beta_corr_LKJ = 2.0)
  ##
  inits_for_Xu_fixed_thr_model <- list(beta_mu = c(-1, +1),
                                       beta_SD = c(0.01, 0.01),
                                       beta_z = array(0.01, dim = c(2, n_studies)),
                                       C_raw_vec = rep(-2.0, n_thr),
                                       C = seq(from = -2, to = 2, length = n_thr))
  init_lists_per_chain <- replicate(n_chains, inits_for_Xu_fixed_thr_model, simplify = FALSE)
}


#### ----
{
      model_parameterisation = "Xu"
      ##
      box_cox <- FALSE
      cts <- FALSE
      ##
      random_thresholds <-  TRUE
   #   Dirichlet_random_effects_type <- "SD"
     #  Dirichlet_random_effects_type <- "kappa"
       Dirichlet_random_effects_type <- "alpha"
      ##
      priors <- list(  prior_alpha = rep(1.0, n_thr + 1), 
                       ##
                       prior_beta_mu_mean = c(0.0, 0.0),
                       prior_beta_mu_SD = c(1.0, 1.0),
                       prior_beta_SD_mean = c(0.0, 0.0),
                       prior_beta_sd_SD = c(0.5, 0.5),
                       ##
                       prior_beta_corr_LKJ = 2.0, 
                       ##
                       prior_dirichlet_cat_means_alpha = rep(1.0, n_thr + 1),
                       ##
                       prior_dirichlet_cat_SDs_mean    = rep(0.0, n_thr + 1),
                       prior_dirichlet_cat_SDs_SD      = rep(0.025, n_thr + 1),
                       ##
                       kappa_lb = 1.0,
                       prior_kappa_mean = rep(0, 1),
                       prior_kappa_SD = rep(500, 1),
                       ##
                       alpha_lb = 1.0,
                       prior_alpha_mean = rep(0, 1),
                       prior_alpha_SD = rep(50, 1)
                       )
      ## inits:
      cutpoint_vec <- seq(from = -3.0, to = 3.0, length = n_thr)
      C_array <- array(dim = c(n_studies, n_thr))
      for (s in 1:n_studies) { 
        C_array[s, ] <- cutpoint_vec
      }
      ##
      n_sets_of_C <- 1
      C_raw <- list()
      for (c in 1:n_sets_of_C) {
         C_raw_mat <- array(-2.0, dim = c(n_studies, n_thr))
         C_raw[[c]] <- C_raw_mat
      }
      ##
      # rand_simplex <- gtools::rdirichlet(n = 1, alpha = rep(10, n_cat))
      phi_raw <- generate_inits_for_raw_simplex_vec(n_thr = n_thr, seed = seed, width = 0.01)
      ##
      dirichlet_cat_means_phi <- rep( 1/(n_thr + 1), n_thr + 1)
      inits_for_Xu_rand_thr_model <- list( beta_mu = c(-1, +1),
                                           beta_SD = c(0.01, 0.01),
                                           beta_z = array(0.01, dim = c(2, n_studies)),
                                           ##
                                           # dirichlet_cat_means_phi = list(c(rand_simplex), c(rand_simplex)),
                                           C_raw_vec = rep(-2.0, n_thr),
                                           C_array = C_array, 
                                           C_raw = C_raw_mat,
                                           ##
                                           alpha = (rep(10.0, n_thr + 1)),
                                           ##
                                           dirichlet_cat_means_phi = list(dirichlet_cat_means_phi),
                                           kappa = rep(1000, 1),
                                           ##
                                           phi_raw = phi_raw
                                           )
      init_lists_per_chain <- replicate(n_chains, inits_for_Xu_rand_thr_model, simplify = FALSE)
      
      dirichlet_cat_means_phi*1000
      
 

      
      
}
  

MetaOrdDTA:::R_fn_compile_stan_model_basic_given_file_name( stan_model_file_name = "DTA_NMA_Nyaga_Jones.stan",
                                               cts = TRUE,
                                               network = TRUE,
                                               prior_only = FALSE,
                                               force_recompile = FALSE,
                                               debugging = TRUE)

inits_for_Xu_rand_thr_model$alpha
##
## Softplus or exp:
##
softplus <- TRUE
network  <- TRUE 

model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  x = x_subset, 
                                                  ##
                                                  # priors = priors,
                                                  ##
                                                  n_chains = n_chains,
                                                  ##
                                                  cts = cts,
                                                  network = network,
                                                  prior_only = FALSE,
                                                  ##
                                                  softplus = softplus,
                                                  ##
                                                  box_cox = box_cox,
                                                  ##
                                                  model_parameterisation = model_parameterisation,
                                                  random_thresholds = random_thresholds,
                                                  Dirichlet_random_effects_type = Dirichlet_random_effects_type, ## only used if random cutpoints
                                                  ##
                                                  init_lists_per_chain = init_lists_per_chain)

 # model_prep_obj$internal_obj$outs_data$stan_data_list


##
## ----  Sample model: ----------------------------------------------------------------
##
model_samples_obj <-  model_prep_obj$sample(   n_burnin = 500,
                                               n_iter = 500,
                                               adapt_delta = 0.80, 
                                               max_treedepth = 10,
                                               metric_shape = "diag_e",
                                               ##
                                               priors = priors,
                                               ##
                                               n_chains = n_chains,
                                               ##
                                               init_lists_per_chain = init_lists_per_chain)


##
## ----  Summarise + output results: -------------------------------------------------
##
model_summary_and_trace_obj <- model_samples_obj$summary(
                                          compute_main_params = TRUE,
                                          compute_transformed_parameters = TRUE, 
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = TRUE,
                                          ##
                                          use_BayesMVP_for_faster_summaries = TRUE)

model_summary_and_trace_obj$get_summary_main() %>% print(n = 1000)

##
## ----  Plots: ----------------------------------------------------------------------
##
# ##
# ## MCMC trace plots:
# ##
# model_summary_and_trace_obj$plot_traces(param_string = c("Se", "Sp"))
# ##
# ## Densities:
# ##
# model_summary_and_trace_obj$plot_densities(param_string = c("Se", "Sp"))
##
## sROC plots:
##
# model_summary_and_trace_obj$plot_sROC()
##
df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, 
                  Sp_true = true_Sp_OVERALL_weighted/100, 
                  Fp_true = (100 - true_Sp_OVERALL_weighted)/100)


df_true
###
model_summary_and_trace_obj$plot_sROC(df_true = df_true)


tibble_gq <- model_summary_and_trace_obj$get_summary_generated_quantities()
tibble_Se <- filter(tibble_gq, str_detect(parameter, "Se"))
tibble_Sp <- filter(tibble_gq, str_detect(parameter, "Sp"))



if (model_parameterisation == "Jones") {

    Se_Jones <- head(tibble_Se$mean, n_thr)
    Sp_Jones <- head(tibble_Sp$mean, n_thr)

} else { 
  
    Se_Cerullo <- head(tibble_Se$mean, n_thr)
    Sp_Cerullo <- head(tibble_Sp$mean, n_thr)
    
}


mean(abs(100*df_true$Se_true - 100*Se_Jones))
max(abs(100*df_true$Se_true  - 100*Se_Jones))
sum(abs(100*df_true$Se_true  - 100*Se_Jones))
##
mean(abs(100*df_true$Sp_true - 100*Sp_Jones))
max(abs(100*df_true$Sp_true  - 100*Sp_Jones))
sum(abs(100*df_true$Sp_true  - 100*Sp_Jones))
####
mean(abs(100*df_true$Se_true - 100*Se_Cerullo))
max(abs(100*df_true$Se_true  - 100*Se_Cerullo))
sum(abs(100*df_true$Se_true  - 100*Se_Cerullo))
##
mean(abs(100*df_true$Sp_true - 100*Sp_Cerullo))
max(abs(100*df_true$Sp_true  - 100*Sp_Cerullo))
sum(abs(100*df_true$Sp_true  - 100*Sp_Cerullo))
####
####
####
mean(abs(100*df_true$Se_true - 100*Se_Cerullo)) - mean(abs(100*df_true$Se_true - 100*Se_Jones))
max(abs(100*df_true$Se_true  - 100*Se_Cerullo)) - max(abs(100*df_true$Se_true  - 100*Se_Jones))
sum(abs(100*df_true$Se_true  - 100*Se_Cerullo)) - sum(abs(100*df_true$Se_true  - 100*Se_Jones))
##
mean(abs(100*df_true$Sp_true - 100*Sp_Cerullo)) - mean(abs(100*df_true$Sp_true - 100*Sp_Jones))
max(abs(100*df_true$Sp_true  - 100*Sp_Cerullo)) - max(abs(100*df_true$Sp_true  - 100*Sp_Jones))
sum(abs(100*df_true$Sp_true  - 100*Sp_Cerullo)) - sum(abs(100*df_true$Sp_true  - 100*Sp_Jones))


# 
# 
# model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_mu_mean
# model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_mu_SD
# model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_SD_mean
# model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_SD_SD
# model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_corr_LKJ
# 
# model_samples_obj$advanced_model_options
# model_samples_obj$basic_model_options
# 
# model_samples_obj$priors
# 
# model_samples_obj$outs_data$stan_data_list
# 
# model_samples_obj$outs_data$stan_data_list
# 
# 
# inits_for_Xu_fixed_thr_model
# 
# 
# debugging <- FALSE
# ##
# n_iter = 500
# n_burnin = 500
# ##
# priors <- NULL
# ##
# x = x
# n_chains = n_chains
# ##
# cts = FALSE
# network = FALSE
# prior_only = FALSE
# ##
# softplus = TRUE
# ##
# ##
# model_parameterisation = "Xu"
# random_thresholds = FALSE
# Dirichlet_random_effects_type = "none"
# box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
# ##
# init_lists_per_chain = inits_for_Xu_fixed_thr_model
# 
# 
# basic_model_options$network <- network
# basic_model_options$cts <- cts
# basic_model_options$prior_only <- prior_only
# ##
# advanced_model_options$model_parameterisation <- model_parameterisation
# advanced_model_options$random_thresholds <- random_thresholds
# advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
# advanced_model_options$box_cox <- box_cox
# advanced_model_options$softplus <- softplus
# ##
# MCMC_params$n_chains <- n_chains
# 
# 
# 
# {
# basic_model_options = list(
#   network = NULL,
#   cts = NULL,
#   prior_only = NULL
# )
# ####
# advanced_model_options = list(
#   model_parameterisation = NULL,
#   random_thresholds = NULL,
#   Dirichlet_random_effects_type = NULL,
#   box_cox = NULL,
#   softplus = NULL
# )
# ####
# priors = NULL
# ####
# # init_lists_per_chain = NULL,
# ####
# MCMC_params = list(
#   seed = NULL,
#   n_superchains = NULL,
#   n_chains = NULL,
#   n_iter = NULL,
#   n_burnin = NULL,
#   adapt_delta = NULL,
#   max_treedepth = NULL,
#   metric_shape = NULL
# )
# }
# 
# 
# 
# 
# 
# outs_data = list(
#   stan_data_list = NULL,
#   n_tests = NULL,
#   n_studies = NULL,
#   n_thr = NULL,
#   n_cat = NULL
# )
# ##
# outs_stan_model_name = list( 
#   stan_model_file_name = NULL
# )
# ##
# outs_stan_compile = list( 
#   stan_model_obj = NULL, 
#   stan_model_file_name = NULL,
#   stan_model_file_path = NULL,
#   ##
#   pkg_root_directory = NULL,
#   stan_models_directory = NULL,
#   stan_functions_directory = NULL,
#   ##
#   stan_MA_directory = NULL,
#   stan_MA_prior_directory = NULL,
#   ##
#   stan_NMA_directory = NULL,
#   stan_NMA_directory = NULL
# )
# ##
# outs_stan_init = list( 
#   inits_unconstrained_vec_per_chain = NULL,
#   stan_param_names_list = NULL,
#   stan_param_names_main = NULL,
#   stan_init_pseudo_sampling_outs = NULL,
#   stan_model_obj = NULL,
#   json_file_path = NULL,
#   stan_model_file_path = NULL
# )
# ##
# outs_stan_sampling = list( 
#   stan_mod_samples = NULL,
#   time_total = NULL
# )
# ##
# ## ---- Main "internal_obj" list:
# ##
# internal_obj = list(
#   outs_data = NULL,
#   outs_stan_model_name = NULL,
#   outs_stan_compile = NULL,
#   outs_stan_init = NULL,
#   outs_stan_sampling = NULL,
#   ##
#   HMC_info = NULL, ## new
#   efficiency_info = NULL, ## new
#   summaries = NULL, ## new
#   traces = NULL ## new
# )
# 
# 
# 
# 
# 
# 
#                                                 
# # 
# debugging = FALSE
# ##
# internal_obj <- model_prep_obj$internal_obj
# ##
# basic_model_options <- model_prep_obj$basic_model_options
# advanced_model_options <- model_prep_obj$advanced_model_options
# MCMC_params <- model_prep_obj$MCMC_params
# # ##
# # init_lists_per_chain <- model_prep_obj$init_lists_per_chain
# # ##
# # priors <- model_prep_obj$priors
# # ##
# # # model_prep_obj$internal_obj
# # # model_prep_obj$internal_obj$outs_data$Stan_data_list
# # # ##
# # # model_prep_obj$basic_model_options
# # # model_prep_obj$advanced_model_options
# # # model_prep_obj$MCMC_params
# # # ##
# # # model_prep_obj$init_lists_per_chain
# # # 
# # # model_prep_obj$priors
# # # 
# # #  
# # # model_prep_obj$Stan_data_list$cutpoint_index

##
## ----  Sample model:
##
model_samples_obj <-  model_prep_obj$sample(   n_burnin = 500,
                                               n_iter = 500,
                                               ##
                                               n_chains = n_chains,
                                               ##
                                               init_lists_per_chain = init_lists_per_chain)
##
# str(model_samples_obj)
##
## ----  Summarise + output results: --------------------------------------------------
# ##
# # 
# # # initial_object <- self$initial_object
# # ##
# # Stan_model_sample_output <- model_samples_obj$model_samples_obj$full_outputs$Stan_model_sample_output ; Stan_model_sample_output
# # # stan_param_names_list    <- model_samples_obj$model_samples_obj$internal_obj$$model_samples_obj$stan_param_names_list ; stan_param_names_list
# # # 
# internal_obj    <- model_samples_obj$internal_obj ## $model_samples_obj$internal_obj
# MCMC_params     <- model_samples_obj$MCMC_params
# # internal_obj
# internal_obj
# model_samples_obj$internal_obj$Stan_data_list
# ##
# model_samples_obj$internal_obj$Stan_model_file_name
# model_samples_obj$internal_obj$Stan_model_file_path
# model_samples_obj$internal_obj$json_file_path ##### 
# model_samples_obj$internal_obj$stan_functions_dir
# ##
# model_samples_obj$internal_obj$inits_unconstrained_vec_per_chain
# model_samples_obj$internal_obj$Stan_init_pseudo_sampling_outs
# ##
# model_samples_obj$internal_obj$stan_param_names_list
# model_samples_obj$internal_obj$stan_param_names_list
# model_samples_obj$internal_obj$stan_model_obj
# ##
# model_samples_obj$internal_obj$Stan_mod_samples
# model_samples_obj$internal_obj$time_total
 
# # 
# # model_samples_obj$model_samples_obj$internal_obj$stan_param_names_list
# # 
# #   
# #   
# # ##
# ## ----
##
model_summary_and_trace_obj <- model_samples_obj$summary(
                                          compute_main_params = TRUE,
                                          compute_transformed_parameters = TRUE, 
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = TRUE,
                                          ##
                                          use_BayesMVP_for_faster_summaries = FALSE)

# # 
# # str(model_samples_obj$internal_obj)
# # 
# # stan_param_names_list <- internal_obj$stan_param_names_list
# # stan_param_names_list$parameters
# # 
# model_summary_and_trace_obj <- create_summary_and_traces( 
#                                     package = "MetaOrdDTA",
#                                     ##
#                                     use_BayesMVP_for_faster_summaries = TRUE,
#                                     use_bridgestan = FALSE,
#                                     ##
#                                     compute_main_params = TRUE,
#                                     compute_transformed_parameters = TRUE,
#                                     compute_generated_quantities = TRUE,
#                                     ##
#                                     internal_obj = internal_obj,
#                                     MCMC_params = MCMC_params,
#                                     ##
#                                     save_log_lik_trace = TRUE,
#                                     ##
#                                     compute_nested_rhat = FALSE,
#                                     ##
#                                     save_trace_tibbles = FALSE)
# 
# str(model_summary_and_trace_obj)
# 
# 
# 
# package = "MetaOrdDTA"
# ##
# use_BayesMVP_for_faster_summaries = TRUE
# use_bridgestan = FALSE
# ##
# compute_main_params = TRUE
# compute_transformed_parameters = TRUE
# compute_generated_quantities = TRUE
# ##
# internal_obj = internal_obj
# MCMC_params = MCMC_params
# ##
# save_log_lik_trace = TRUE
# ##
# compute_nested_rhat = FALSE
# ##
# save_trace_tibbles = FALSE
# 
# 




  
# str(model_summary_and_trace_obj$traces$traces_as_arrays$draws_array)
# 
# draws_array <- model_summary_and_trace_obj$traces$traces_as_arrays$draws_array
# 
# # 
# # package = "MetaOrdDTA"
# # use_BayesMVP_for_faster_summaries = TRUE
# # use_bridgestan = FALSE
# # compute_main_params = TRUE
# # compute_transformed_parameters = TRUE
# # compute_generated_quantities = TRUE
# # ##
# # save_log_lik_trace = TRUE,
# # Stan_model_sample_output = model_samples_obj$model_samples_obj$full_outputs$Stan_model_sample_output
# # stan_param_names_list = model_samples_obj$model_samples_obj$stan_param_names_list
# # ##
# # compute_nested_rhat = FALSE
# # 
# 


# model_summary_and_trace_obj
##
###extract # divergences + % of sampling iterations which have divergences:

# model_summary_and_trace_obj$internal_obj$HMC_info
# model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_main_params
# model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters
# model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
##
# 
# 
# get_test <- model_summary_and_trace_obj$get_HMC_info() ; str(get_test)
# get_test
# 
# get_test <- model_summary_and_trace_obj$get_HMC_diagnostic_info() ; str(get_test)
# get_test
# # get_test$divergences
# get_test <- model_summary_and_trace_obj$get_divergences() ; str(get_test)
# get_test
# 
# 
# # get_test <- model_summary_and_trace_obj$time_to_target_ESS(target_ESS = 1000) ; str(get_test)
# 
# get_test <- model_summary_and_trace_obj$get_divergences() ; str(get_test)
# 
# 
# # get_test <- model_summary_and_trace_obj$get_summary_main() ; (get_test)
# # get_test <- model_summary_and_trace_obj$get_summary_transformed() ; (get_test) ## not working?
# # get_test <- model_summary_and_trace_obj$get_summary_generated_quantities() ; (get_test)
# 
# 
# model_summary_and_trace_obj$

# 
# 
# model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
# 
# 
# model_summary_and_trace_obj$model_summary_and_trace_obj$HMC_diagnostic_info$divergences
# 
# str(model_summary_and_trace_obj)
# model_summary_and_trace_obj$model_summary_and_trace_obj
# 
# 
# model_summary_and_trace_obj$get_efficiency_metrics() 
# 
# model_summary_and_trace_obj$time_to_target_ESS(5000)
# 
# model_summary_and_trace_obj$get_HMC_info()
# model_summary_and_trace_obj$get_HMC_diagnostic_info()
# 


##
## MCMC trace plots:
##
model_summary_and_trace_obj$plot_traces(param_string = c("Se", "Sp"))
##
## Densities:
##
model_summary_and_trace_obj$plot_densities(param_string = c("Se", "Sp"))
##
## sROC plots:
##
model_summary_and_trace_obj$plot_sROC()
##
df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, 
                  Sp_true = true_Sp_OVERALL_weighted/100, 
                  Fp_true = (100 - true_Sp_OVERALL_weighted)/100)


df_true
###
model_summary_and_trace_obj$plot_sROC(df_true = df_true)


 



# ##
# ## ---- Get HMC info:
# ##
# model_HMC_info <- model_summary_and_trace_obj$get_HMC_info()
# model_HMC_info
# ##
# str(model_summary_and_trace_obj$get_all_traces_list())
# get_all_traces_list <- model_summary_and_trace_obj$get_all_traces_list()
# get_all_traces_list$traces_as_arrays$trace_params_main
# 
# get_all_traces_list <- list()
# 
# ##
# trace_tp <- get_all_traces_list$traces_as_arrays$trace_tp
# trace_gq <- get_all_traces_list$traces_as_arrays$trace_gq
# 
# ## 3-D array: An array with dimensions ⁠Iteration, Chain, Parameter⁠ in that order.
# str(trace_params_main)
# dimnames(trace_params_main[,,3])
# 
# names <- dimnames(trace_gq[,,3])
# head(names$variable, 100)




draws_array <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$draws_array
draws_array_2 <- format_named_array_for_bayesplot(draws_array)

str(draws_array)


# trace_params_main_reformatted <- format_named_array_for_bayesplot(trace_params_main)
# trace_tp_reformatted <- format_named_array_for_bayesplot(trace_tp)
# trace_gq_reformatted <- format_named_array_for_bayesplot(trace_gq)
# trace_list <- list(trace_params_main_reformatted, 
#                    trace_tp_reformatted,
#                    trace_gq_reformatted)

# ############
# draws_array <-  MetaOrdDTA:::combine_draws_arrays(trace_list)
# str(trace_list)
# # str(draws_array)
# 
# names$variable

draws_array = draws_array_2
param_string = c("Se")
condition = "containing"
print = TRUE
plot_type = "density"
batch_size = 9


## ⁠Iteration, Chain, Parameter⁠ in that order.


outs <- MetaOrdDTA:::plot_param_group_batched(  draws_array = draws_array,
                                   plot_type = "density",
                                   param_string = c("Se", "Sp"),
                                   condition = "exact_match",
                                   batch_size = 9,
                                   debugging = TRUE)

outs[[1]][[1]]

outs[[1]][[1]]
 
outs[[2]][[1]]
##
plot_param_group_batched

str(outs)
outs
all_param_names_vec = names$variable
param_strings_vec = c("Se")
single_param_string =  c("Se", "Sp")
  
condition = "containing"
print = F
batch_size = 9

plot_type = "density"

str(draws_array)

filter_param_names_string_batched(
                         all_param_names_vec = names$variable,
                         param_strings_vec = c("Se"),
                         condition = "containing",
                         print = TRUE)##,



traces <- model_summary_and_trace_obj$get_all_traces()
str(traces)
##
traces_as_arrays <- traces$traces_as_arrays
str(traces_as_arrays)
##
log_lik_trace <- traces$log_lik_trace
str(log_lik_trace)
##
model_summary_and_trace_obj$get_HMC_info()



model_summary_and_trace_obj$model_summary_and_trace_obj$efficiency_info
model_summary_and_trace_obj$get_efficiency_metrics()
all_traces <- model_summary_and_trace_obj$get_all_traces()
str(all_traces)


model_samples_obj$internal_obj$traces
model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays$trace_params_main



str(model_summary_and_trace_obj$model_summary_and_trace_obj$summaries)
##
model_summary_and_trace_obj$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_main_params
model_summary_and_trace_obj$get_summary_main()
##
model_summary_and_trace_obj$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters
model_summary_and_trace_obj$get_summary_transformed()
##
model_summary_and_trace_obj$model_summary_and_trace_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
model_summary_and_trace_obj$get_summary_generated_quantities()


model_summary_and_trace_obj$plot_sROC()



##
##
## ----  MCMC plots - trace plots: ----------------------------------------------------------------------------
##
## E.g., Let's plot the densities for: summary sensitivity and specificity:
trace_plots <- model_summary_and_trace_obj$plot_traces(params = c("Se", "Sp", "p"), 
                                         batch_size = 12)
# you can extract parameters by doing: "trace$param_name()". 
# For example:
# display each panel for beta and Omega ("batch_size" controls the # of plots per panel. Default = 9)
trace_plots$Se[[1]] # Se - 1st (and only) panel
trace_plots$Sp[[1]] # Sp - 1st (and only) panel
##
## ---- MCMC plots - density plots: ----------------------------------------------------------------------------
##
## E.g., Let's plot the densities for: summary sensitivity and specificity:
density_plots <- model_summary_and_trace_obj$plot_densities(params = c("Se", "Sp", "p"), 
                                              batch_size = 12)
density_plots$Se[[1]] # Se - 1st (and only) panel
density_plots$Sp[[1]] # Sp - 1st (and only) panel
##
## ---- Test-accuracy specific plots: --------------------------------------------------------------------------
## ---- (e.g., sROC curves, confidence and prediction regions, etc):
##
## ---------------------BOOKMARK - PUT SROC AND OTHER TEST-ACCURACY PLOTS HERE ---------------------------------
##
## ----------------- Test accuracy-specific plots (e.g. sROC curve / confidence + prediction regions) ----------------
##
##
## sROC plot of results:
##
Stan_model_file_name <- model_samples_obj$Stan_model_file_name ; Stan_model_file_name
 ## model_samples_obj$model_summary_and_trace_obj ## $full_outputs$Stan_model_sample_output$Stan_mod_samples

 Stan_mod_samples     <- model_summary_and_trace_obj$model_samples_obj$full_outputs$Stan_model_sample_output$Stan_mod_samples
 Stan_model_file_path     <- model_summary_and_trace_obj$model_samples_obj$Stan_model_file_path##full_outputs$Stan_model_sample_output$Stan_mod ##$model_samples_obj$full_outputs$Stan_model_sample_output$Stan_mod_samples
 
 
  str(Stan_mod_samples)
  ##

model_summary_and_trace_obj

sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( Stan_model_file_name = Stan_model_file_name,
                                               Stan_mod_samples = Stan_mod_samples,
                                               df_true = NULL,
                                               conf_region_colour = "blue",
                                               pred_region_colour = "blue")
##
df_fitted <- sROC_plot_outs$df_fitted
plot_1 <- sROC_plot_outs$plot_1
plot_2 <- sROC_plot_outs$plot_2
df_fitted
plot_1
plot_2
plot_list
##
print(plot_2)
##
## sROC plot of results (+ true values over-layed) 
## NOTE: only use this if you know the true values for Se/Sp (e.g., simulation study):
##
df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, 
                  Sp_true = true_Sp_OVERALL_weighted/100, 
                  Fp_true = (100 - true_Sp_OVERALL_weighted)/100)


df_true
###
sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( Stan_model_file_name = Stan_model_file_name,
                                               Stan_mod_samples = Stan_mod_samples,
                                               df_true = df_true,
                                               conf_region_colour = "blue",
                                               pred_region_colour = "blue")
##
df_fitted <- sROC_plot_outs$df_fitted
plot_1 <- sROC_plot_outs$plot_1
plot_2 <- sROC_plot_outs$plot_2
plot_list
##
print(plot_2)
##
##
##
##
##
##
##
## ---- Other features: ----------------------------------------------------------------------------------------
##
## The "model_summary" object (created using the "$summary()" method) contains many useful objects. 
## For example:
require(dplyr)
## nice summary tibble for main parameters, includes ESS/Rhat, etc:
model_summary$get_summary_main() %>% print(n = 50) 
## nice summary tibble for transformed parameters, includes ESS/Rhat, etc:
model_summary$get_summary_transformed() %>% print(n = 150) 
## nice summary tibble for generated quantities, includes ESS/Rhat, etc:
model_summary$get_summary_generated_quantities () %>% print(n = 150) 
##
## users can also easily use the "posterior" R package to compute their own statistics. 
## For example:
# let's say we want to compute something not included in the default
# "$summary()" method of BayesMVP, such as tail-ESS.
# We can just use the posterior R package to compute this:
require(posterior)  
## first extract the trace array object (note: already in a posterior-compatible format!)
posterior_draws <- model_summary$get_posterior_draws()
# then compute tail-ESS using posterior::ess_tail:
posterior::ess_tail(posterior_draws[,,"Se_bin[1]"])

## You can also get the traces as tibbles (stored in seperate tibbles for main params, 
## transformed params, and generates quantities) using the "$get_posterior_draws_as_tibbles()" method:
tibble_traces  <- model_summary$get_posterior_draws_as_tibbles()
tibble_trace_main <- tibble_traces$trace_as_tibble_main_params
tibble_trace_transformed_params <- tibble_traces$trace_as_tibble_transformed_params
tibble_trace_generated_quantities <- tibble_traces$trace_as_tibble_generated_quantities

## You can also easily extract model run time / efficiency information using the "$get_efficiency_metrics()" method:
model_efficiency_metrics <- model_summary$get_efficiency_metrics()
time_burnin <- model_efficiency_metrics$time_burnin  ; time_burnin
time_sampling <- model_efficiency_metrics$time_sampling ; time_sampling
time_total_MCMC <- model_efficiency_metrics$time_total_MCMC  ; time_total_MCMC
## Note that the following ("time_total_inc_summaries") ncludes time to compute R-hat, etc:
time_total_inc_summaries <- model_efficiency_metrics$time_total_inc_summaries ; time_total_inc_summaries #

## We can also extract some more specific efficiency info, again using the "$get_efficiency_metrics()" method:
Min_ESS_main_params <- model_efficiency_metrics$Min_ESS_main   ; Min_ESS_main_params
Min_ESS_per_sec_sampling <- model_efficiency_metrics$Min_ESS_per_sec_samp ; Min_ESS_per_sec_sampling
Min_ESS_per_sec_overall <- model_efficiency_metrics$Min_ESS_per_sec_total ; Min_ESS_per_sec_overall

Min_ESS_per_grad <- model_efficiency_metrics$Min_ESS_per_grad_sampling ; Min_ESS_per_grad
grad_evals_per_sec <- model_efficiency_metrics$grad_evals_per_sec ; grad_evals_per_sec

## extract the "time to X ESS" - these are very useful for knowing how long to
##  run your model for. 
est_time_to_100_ESS <- model_efficiency_metrics$est_time_to_100_ESS_inc_summaries ; est_time_to_100_ESS
est_time_to_1000_ESS <- model_efficiency_metrics$est_time_to_1000_ESS_inc_summaries ; est_time_to_1000_ESS
est_time_to_10000_ESS <- model_efficiency_metrics$est_time_to_10000_ESS_inc_summaries; est_time_to_10000_ESS

##  You can also use the "model_samples$time_to_ESS()" method to estimate 
## "time to X ESS" for general X:
##  For example let's say we determined our target (min) ESS to be ~5000:
est_time_5000_ESS <- model_summary$time_to_target_ESS(target_ESS = 5000) ; est_time_5000_ESS
est_time_5000_ESS

### You can also extract the log_lik trace (note: for Stan models this will only work )
log_lik_trace <- model_summary$get_log_lik_trace()
str(log_lik_trace) # will be NULL unless you specify  "save_log_lik_trace = TRUE" in the "$summary()" method
## can then use log_lik_trace e.g. to compute LOO-IC using the loo package 

















# ##
# stan_param_names_list <- model_samples$result$stan_param_names_list
# Stan_model_sample_output <- model_samples$result$full_outputs$Stan_model_sample_output
# 
# outs <- create_summary_and_traces( package = "MetaOrdDTA",
#                                    use_bridgestan = NULL,
#                                    save_log_lik_trace = TRUE,
#                                    ##
#                                    Stan_model_sample_output = Stan_model_sample_output,
#                                    stan_param_names_list = stan_param_names_list,
#                                    ##
#                                    compute_transformed_parameters = FALSE,
#                                    compute_generated_quantities = TRUE)
# 
# {
#   package = "MetaOrdDTA"
#   use_bridgestan = NULL
#   save_log_lik_trace = TRUE
#   ##
#   Stan_model_sample_output = Stan_model_sample_output
#   stan_param_names_list = stan_param_names_list
#   ##
#   compute_transformed_parameters = FALSE
#   compute_generated_quantities = TRUE
# }
#                            
  
# MetaOrdDTA:::initialise_update_run_model(x = x,
#                                          ##
#                                          n_chains = 8, 
#                                          test_type = "ord",
#                                          ##
#                                          network = FALSE,
#                                          cts = FALSE,
#                                          ##
#                                          prior_only = FALSE,
#                                          ##
#                                          model_args_list = NULL,
#                                          ##
#                                          priors = NULL,
#                                          ##
#                                          ## "Advanced" options for ordinal models ONLY:
#                                          ##
#                                          model_parameterisation = NULL,
#                                          random_thresholds = NULL,
#                                          Dirichlet_random_effects_type = NULL,
#                                          ##
#                                          ## "Advanced" options for cts (i.e., Jones) models ONLY:
#                                          ##
#                                          box_cox = NULL,
#                                          ##
#                                          ## "Advanced" options for ALL models:
#                                          ##
#                                          softplus = NULL,
#                                          ##
#                                          ##
#                                          ##
#                                          seed = 123,
#                                          n_burnin = 500,
#                                          n_iter = 500,
#                                          adapt_delta = 0.80,
#                                          max_treedepth = 10,
#                                          metric_shape = "diag_e",
#                                          n_superchains = 8)
  
# ##
# ## * MIGHT * have to restart session:
# #### rstudioapi::restartSession()
# ##
# ##
# {
#   install_success <- FALSE
#   
#   try({  
#     # ## Sometimes only works if do this first (on Linux) :
#     # source("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/MetaOrdDTA/src/R_script_load_OMP_Linux.R")
#     # ##
#     ## Install (inner pkg):
#     require(MetaOrdDTA)
#     MetaOrdDTA::install_BayesMVP() ## CUSTOM_FLAGS = CUSTOM_FLAGS)
#     require(MetaOrdDTA)
#     install_success <- TRUE
#   })
#   
#   try({  
#     beepr::beep("ping")
#   })
#   
#   rstudioapi::restartSession()
# }
# 
# 
# 
# 
# 
# 

# ##### Using GitHub:
# ## Install (outer pkg):
# BayesMVP_repo_link <- "https://github.com/CerulloE1996/MetaOrdDTA"
# remotes::install_github(repo = BayesMVP_repo_link, upgrade = "never")
# ###
# ## Install (inner pkg):
# require(MetaOrdDTA)
# MetaOrdDTA::install_BayesMVP()
# require(MetaOrdDTA)






# CUSTOM_FLAGS <- list(  "CCACHE_PATH = ccache",
#                        "CXX_COMPILER = g++",
#                        "CPP_COMPILER = gcc", 
#                        "CXX_STD = CXX17",
#                        "CPU_BASE_FLAGS = -O3 -march=native -mtune=native",
#                        "FMA_FLAGS = -mfma",
#                        "AVX_FLAGS = -mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl")
#                        # "OMP_FLAGS",
#                        # "OMP_LIB_PATH",
#                        # "OMP_LIB_FLAGS",
#                        # "PKG_CPPFLAGS",
#                        # "PKG_CXXFLAGS",
#                        # "CPPFLAGS",
#                        # "CXXFLAGS",
#                        # "PKG_LIBS")
# 
# 




 
#  
# 
# 
# 
# 
# 
# ### USER INSTALLATION PROCESS: --------------------------------------------------------------------------------------------
# 
# 
# #### ----------------- install cmdstanr first:
# 
# #### Install the cmdstanr "outer" R package:
# remotes::install_github("stan-dev/cmdstanr", force = TRUE)
# #### Load cmdstanr R outer package:
# require(cmdstanr) 
# #### Then, install the latest version of CmdStan:
# install_cmdstan(cores = parallel::detectCores(),
#                 overwrite = TRUE,
#                 cpp_options = list("STAN_MODEL_LDFLAGS" = "-shared",   "CXXFLAGS" = "-fPIC"))
# 
# 
# 
# 
# 
# 
# 
# #### ----------------- Then install bridgestan:
# remotes::install_github("https://github.com/roualdes/bridgestan", subdir="R")
# #### Load bridgestan:
# require(bridgestan)
# 
# 
# 
# 
# 
# #### ----------------- Then install BayesMVP:
# ## First remove any possible package fragments:
# ## Find user_pkg_install_dir:
# user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
# print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
# ##
# ## Find pkg_install_path + pkg_temp_install_path:
# pkg_install_path <- file.path(user_pkg_install_dir, "BayesMVP")
# pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-BayesMVP") 
# ##
# ## Remove any (possible) BayesMVP package fragments:
# remove.packages("BayesMVP")
# unlink(pkg_install_path, recursive = TRUE, force = TRUE)
# unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
# 
# 
# ## Install OUTER R package from GitHub:
# remotes::install_github("https://github.com/CerulloE1996/BayesMVP", force = TRUE, upgrade = "never")
# ## Then restart R session:
# rstudioapi::restartSession()
# ## Then install INNTER (i.e. the "real") package:
# require(BayesMVP)
# BayesMVP::install_BayesMVP()
# require(BayesMVP)
# 
# 
# ## Restart session:
# rstudioapi::restartSession()
# ## Install INNER R package:
# BayesMVP:::install_BayesMVP()
# ## Restart session:
# rstudioapi::restartSession()
# 
# 
# ###### other / testing
# 
# setwd(local_INNER_pkg_dir)
# Rcpp::sourceCpp(file.path(local_INNER_pkg_dir, "src", "cpu_check.cpp"))
# checkCPUFeatures()
# 
# 
# 








