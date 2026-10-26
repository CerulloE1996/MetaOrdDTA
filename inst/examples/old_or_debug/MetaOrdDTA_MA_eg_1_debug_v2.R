

 

os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")




# # 
# # ##
rm(list = ls())
# ##
# .rs.restartR()  # In RStudio
# #


os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")




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


#
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
#
# 






{
    
  # ##
  # ## ---- Load NMA data:
  # ##
  setwd(local_pkg_dir)
  ##
  #### seed <- 123
  ##
  seed <- 1
  # seed <- 2
  # seed <- 3
  # seed <- 4
  # seed <- 5
  # seed <- 6
  # seed <- 7
  # seed <- 8
  # seed <- 9
  # seed <- 10
  ##
  n_studies <- 5
  # n_studies <- 10
  ## n_studies <- 25
  ## n_studies <- 50
  ## n_studies <- 100
  #### n_studies <- 250
  # ##
  source(file.path(getwd(), "inst", "examples", "NMA_missing_thr_prep_data.R"))
  ##
  data <- readRDS(file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))
  ##
  x_NMA <- data$x_NMA
  ##
  x_PHQ9 <- data$x_PHQ9
  ## x_PHQ9_w_missing_thr <- data$x_PHQ9_w_missing_thr
  ##
  x_MA <- list()
  for (c in 1:2) {
    x_MA[[c]] <- x_NMA[[c]][[index_test_chosen_index - 1]]
  }
  ##
  ## Select which data ("x") to use:
  ##
  str(x_NMA)
  x <- x_MA
  x <- x_PHQ9
  ## x <- x_PHQ9_w_missing_thr
  ##
  n_thr <- ncol(x[[1]]) - 1
  ##
  x
  n_thr
  
  
  ##
  n_studies <- nrow(x[[1]])
  n_studies

}

##
## ----  Initialise / select model: --------------------------------------------------------------------------- NMA:
##
network <- FALSE
##
softplus <- FALSE




#### ----
{
  model_parameterisation = "Jones"
  box_cox <- TRUE
  cts <- TRUE
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
  # ##
  # inits <- list(beta_mu = c(-1, +1),
  #                                      beta_SD = c(0.01, 0.01),
  #                                      beta_z = array(0.01, dim = c(2, n_studies)),
  #                                      C_raw_vec = rep(-2.0, n_thr),
  #                                      C = seq(from = -2, to = 2, length = n_thr))
  # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
}

# {
#   model_parameterisation = "ord_bivariate"
#   box_cox <- FALSE
#   cts <- FALSE
#   random_thresholds <-  FALSE
#   Dirichlet_random_effects_type <- "fixed"
#   # ##
#   # inits <- list(beta_mu = c(-1, +1),
#   #                                      beta_SD = c(0.01, 0.01),
#   #                                      beta_z = array(0.01, dim = c(2, n_studies)),
#   #                                      C_raw_vec = rep(-2.0, n_thr),
#   #                                      C = seq(from = -2, to = 2, length = n_thr))
#   # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
# }
# 
# {
#   model_parameterisation = "R&G"
#   box_cox <- FALSE
#   cts <- FALSE
#   random_thresholds <-  FALSE
#   Dirichlet_random_effects_type <- "fixed"
#   # ##
#   # inits <- list(beta_mu = c(-1, +1),
#   #                                      beta_SD = c(0.01, 0.01),
#   #                                      beta_z = array(0.01, dim = c(2, n_studies)),
#   #                                      C_raw_vec = rep(-2.0, n_thr),
#   #                                      C = seq(from = -2, to = 2, length = n_thr))
#   # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
# }
# 
# {
#   model_parameterisation = "ord_bivariate"
#   box_cox <- FALSE
#   cts <- FALSE
#   random_thresholds <-  TRUE
#   Dirichlet_random_effects_type <- "alpha"
#   # ##
#   # inits <- list(beta_mu = c(-1, +1),
#   #                                      beta_SD = c(0.01, 0.01),
#   #                                      beta_z = array(0.01, dim = c(2, n_studies)),
#   #                                      C_raw_vec = rep(-2.0, n_thr),
#   #                                      C = seq(from = -2, to = 2, length = n_thr))
#   # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
# }
# 
# # {
# #   model_parameterisation = "R&G"
# #   box_cox <- FALSE
# #   cts <- FALSE
# #   random_thresholds <-  TRUE
# #   Dirichlet_random_effects_type <- "alpha"
# #   # ##
# #   # inits <- list(beta_mu = c(-1, +1),
# #   #                                      beta_SD = c(0.01, 0.01),
# #   #                                      beta_z = array(0.01, dim = c(2, n_studies)),
# #   #                                      C_raw_vec = rep(-2.0, n_thr),
# #   #                                      C = seq(from = -2, to = 2, length = n_thr))
# #   # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
# # }
# # 





n_chains <-  8 ##  parallel::detectCores() / 2
# internal_obj$outs_data$stan_data_list$x
# internal_obj$outs_data$stan_data_list$n
# 
# any(is.na(unlist(internal_obj$outs_data$stan_data_list$n)))
init_lists_per_chain <- NULL
##
##
vec_index_inner_thr <- 5:20
n_inner_thr <- length(vec_index_inner_thr)
vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner_thr)
n_outer_thr <- length(vec_index_outer_thr)
##
##
##
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  
                  debugging = TRUE,
                  ##
                  x = x, 
                  indicator_index_test_in_study = NULL,
                  # n_index_tests_per_study = NULL,
                  ##
                  # priors = priors,
                  ##
                  n_chains = n_chains,
                  ##
                  cts = cts,
                  network = network,
                  ##
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
                  init_lists_per_chain = init_lists_per_chain,
                  ##
                  advanced_compile = FALSE, ## default is standard/"basic" compile
                  ##
                  ## force_recompile = TRUE,
                  ##
                  set_custom_CXX_CPP_flags = TRUE, 
                  CCACHE_PATH = "/usr/bin/ccache", 
                  custom_cpp_user_header_file_path = NULL, ## if wish to include custom C++ files
                  CXX_COMPILER_PATH = "override g++",
                  CPP_COMPILER_PATH = "override gcc",
                  MATH_FLAGS = "-fno-math-errno  -fno-signed-zeros -fno-trapping-math",
                  FMA_FLAGS = "-mfma",
                  AVX_FLAGS = "-mavx -mavx2",
                  THREAD_FLAGS = "-pthread -DSTAN_THREADS",
                  ##
                  compute_sim_study_metrics = 1,
                  vec_index_inner_thr = vec_index_inner_thr
                  )


{
  
  # model_prep_obj$internal_obj$outs_data$stan_data_list
  init_lists_per_chain <- model_prep_obj$init_lists_per_chain
  priors <- model_prep_obj$priors
  ##
  n_studies <- nrow(x[[1]])
  ##  n_index_tests <- 1# length(x[[1]])
  
  stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list

}




# # ##
# model_prep_obj$internal_obj$outs_data$stan_data_list$compute_sim_study_metrics <- 1
# model_prep_obj$internal_obj$outs_data$stan_data_list$vec_index_inner_thr <- vec_index_inner_thr
# model_prep_obj$internal_obj$outs_data$stan_data_list$vec_index_outer_thr <- vec_index_outer_thr
# model_prep_obj$internal_obj$outs_data$stan_data_list$n_inner_thr <- n_inner_thr
# model_prep_obj$internal_obj$outs_data$stan_data_list$n_outer_thr <- n_outer_thr
# ##
model_prep_obj$internal_obj$outs_data$stan_data_list$n_inner_thr

# stan_data_list$x_with_missings
# stan_data_list$x
# stan_data_list$cutpoint_index
# stan_data_list$n_obs_cutpoints
# ##
# stan_data_list$cts_thr_values


# x_wo_totals <- x[[1]][1:n_studies, 2:ncol(x[[1]])]
# x_wo_totals
# which(x_wo_totals[1, ] != -1)
# 
# x[[1]]
#  
# stan_data_list$cts_thr_values
# stan_data_list$x[[2]]
# x[[2]]
# stan_data_list$n[[2]]
# stan_data_list$n_diseased
# stan_data_list$n_obs_cutpoints
# stan_data_list$cutpoint_index
# # 
# # ##
# str(init_lists_per_chain)
# 
# model_prep_obj$internal_obj$outs_data$stan_data_list$box_cox = FALSE

 

init_lists_per_chain <- vector(length = n_chains, "list")
priors <- NULL


for (kk in 1:n_chains) {
  init_lists_per_chain[[kk]]$C_raw_vec <- list( rep(-3.0, n_thr), rep(-3.0, n_thr))
 # init_lists_per_chain[[kk]]$C <- list( seq(from = -2.0, to = 2.0, length = n_thr), 
  #                                      seq(from = -2.0, to = 2.0, length = n_thr))
}

if (model_parameterisation == "ord_bivariate") {
      
          for (kk in 1:n_chains) {
                
                mat <- matrix(-3.0, nrow = n_studies, ncol = n_thr)
                init_lists_per_chain[[kk]]$C_raw <-  list(mat, mat)
                init_lists_per_chain[[kk]]$alpha <-  list( rep(5.0, n_thr + 1), rep(5.0, n_thr + 1))
                ##
                init_lists_per_chain[[kk]]$beta_mu <- c(-1.01, +1.01)
                init_lists_per_chain[[kk]]$beta_z <- array(0.001, dim = c(n_studies, 2))
                init_lists_per_chain[[kk]]$beta_SD <- rep(0.001, 2)
                init_lists_per_chain[[kk]]$beta_corr <- -0.001
                ##
                p_vec <- rep(1.0/(n_thr + 1), n_thr + 1)
                init_lists_per_chain[[kk]]$p_vec <- list(p_vec, p_vec)
                ##
                dirichlet_phi <- rep(1/(n_thr +  1), n_thr + 1)
                init_lists_per_chain[[kk]]$dirichlet_phi <- list(dirichlet_phi, dirichlet_phi)
                kappa <- 50
                init_lists_per_chain[[kk]]$kappa <- c(kappa, kappa) 
                ##
                dirichlet_phi_raw <- MetaOrdDTA:::generate_inits_for_raw_simplex_vec(n_thr = n_thr, seed = 123)
                init_lists_per_chain[[kk]]$dirichlet_phi_raw <- list(dirichlet_phi_raw, dirichlet_phi_raw)
            
          }
            
            # init_lists_per_chain = NULL
            # 
            # priors$prior_raw_scale_mu_SD <- 0.25
            # priors$prior_raw_scale_SD_SD <- 0.125
            # 
            # priors$prior_beta_mu_SD <- 1
            
            priors$prior_alpha_SD
            
            prior_alpha_mean <- 1.0
            prior_alpha_SD   <- 5.0
            alpha_lb <- 1
            
            priors$prior_alpha_mean <- list( rep(prior_alpha_mean, n_thr + 1),  rep(prior_alpha_mean, n_thr + 1))
            priors$prior_alpha_SD   <- list( rep(prior_alpha_SD, n_thr + 1), rep(prior_alpha_SD, n_thr + 1))
            priors$alpha_lb <- alpha_lb
            
            priors$prior_kappa_mean <- 0
            priors$prior_kappa_SD <- 500
            priors$kappa_lb <- 1
            ##
            prior_dirichlet_phi <- rep(1.0, n_thr + 1)
            prior_dirichlet_phi[] <- 1.0
            priors$prior_dirichlet_phi <- list(prior_dirichlet_phi, prior_dirichlet_phi)
            
            ##
            priors$prior_beta_mu_mean <- rep(0.0, 2)
            priors$prior_beta_mu_SD   <- rep(2.5, 2)
            priors$prior_beta_SD_mean <- rep(0.0, 2)
            priors$prior_beta_SD_SD   <- rep(0.5, 2)
            ##
            
            ##
            # priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1) )
            alpha_priors <- generate_dual_ordinal_dirichlet_priors(n_cat = n_thr + 1,
                                                   nd_peak = 0.15,
                                                   d_peak = 0.30,
                                                   nd_max_alpha = 2.5,
                                                   d_max_alpha = 2.5,
                                                   nd_min_alpha = 1.0,
                                                   d_min_alpha = 1.0,
                                                   smoothness = 0.15)
            ##
            priors$prior_dirichlet_alpha <- list( alpha_priors$non_diseased, alpha_priors$diseased)

            
      
} else if (model_parameterisation == "R&G") { 
      
          for (kk in 1:n_chains) {
                mat <- matrix(-3.0, nrow = n_studies, ncol = n_thr)
                init_lists_per_chain[[kk]]$C_raw <- mat
                init_lists_per_chain[[kk]]$alpha <- rep(10.0, n_thr + 1)
            
          }
         ##  mat <- matrix(-3.0, nrow = n_studies, ncol = n_thr)
          for (kk in 1:n_chains) {
             init_lists_per_chain[[kk]]$C_raw_vec <-  rep(-3.0, n_thr)
          }
          
          prior_alpha_mean <- 1
          prior_alpha_SD   <- 25
          alpha_lb <- 1
          ##
          priors$prior_alpha_mean <- rep(prior_alpha_mean, n_thr + 1)
          priors$prior_alpha_SD   <- rep(prior_alpha_SD, n_thr + 1)
          priors$alpha_lb <- alpha_lb

          
          priors$prior_beta_mu_mean <- 0.0
          priors$prior_beta_mu_SD   <- 5.0
          priors$prior_beta_SD_mean <- 0.0
          priors$prior_beta_SD_SD   <- 1.0 ## 0.25 ## on log-normal scale so should be small !!
          ##
          priors$prior_raw_scale_mu_mean <- 0.0
          priors$prior_raw_scale_mu_SD   <- 1.0
          priors$prior_raw_scale_SD_mean <- 0.0
          priors$prior_raw_scale_SD_SD   <- 1.0 ## 0.25 ## on log-normal scale so should be small !!
          
          ##
          ## priors$prior_dirichlet_alpha <-  rep(2.0, n_thr + 1)
          priors$prior_dirichlet_alpha <- generate_ordinal_dirichlet_prior(n_cat = n_thr + 1,
                                                                           max_alpha = 2.5, 
                                                                           min_alpha = 1.0,
                                                                           left_skew_factor = 0.5)
  
} else if (model_parameterisation == "Jones") { 
  
        priors$prior_beta_mu_mean <- rep(0.0, 2)
        priors$prior_beta_mu_SD   <- rep(5.0, 2)
        priors$prior_beta_SD_mean <- rep(0.0, 2)
        priors$prior_beta_SD_SD   <- rep(1.0, 2)
        ##
        priors$prior_raw_scale_mu_mean <- rep(0.0, 2)
        priors$prior_raw_scale_mu_SD   <- rep(1.0, 2)
        priors$prior_raw_scale_SD_mean <- rep(0.0, 2)
        priors$prior_raw_scale_SD_SD   <- rep(1.0, 2)
        ##
        priors$prior_boxcox_lambda_mean <- 0.00
        priors$prior_boxcox_lambda_mean <- 1.00
        for (kk in 1:n_chains) {
          
          # mat <- matrix(-3.0, nrow = n_studies, ncol = n_thr)
          # init_lists_per_chain[[kk]]$C_raw <-  list(mat, mat)
          # init_lists_per_chain[[kk]]$alpha <-  list( rep(10.0, n_thr + 1), rep(10.0, n_thr + 1))
          ##
          init_lists_per_chain[[kk]]$beta_mu <- c(-1.0, +1.0)
          init_lists_per_chain[[kk]]$beta_SD <- rep(0.001, 2)
          init_lists_per_chain[[kk]]$beta_z <- array(0.001, dim = c(n_studies, 2))
          init_lists_per_chain[[kk]]$beta_corr <- 0.001
          ##
          init_lists_per_chain[[kk]]$raw_scale_mu <- c(-1.0, +1.0)
          init_lists_per_chain[[kk]]$raw_scale_SD <- rep(0.001, 2)
          init_lists_per_chain[[kk]]$raw_scale_z <- array(0.001, dim = c(n_studies, 2))
          init_lists_per_chain[[kk]]$raw_scale_corr <- 0.001
          ##
          init_lists_per_chain[[kk]]$lambda <- rep(0.001, 2)

          
        }
        
  
}

{

    priors$beta_corr_lb <- -1.0
    priors$beta_corr_ub <- +1.0 ## 0.0
    
    init_lists_per_chain[[1]]
    
    model_parameterisation
    random_thresholds
    softplus
    
    try({  
    for (c in 1:2) {
        priors$prior_dirichlet_phi[[c]][] <- 2.0
    }
    })
    priors$prior_dirichlet_phi
    
    priors$prior_kappa_SD <- 50
    priors$kappa_lb <- 0
    
    priors$prior_dirichlet_alpha
    #### init_lists_per_chain  = NULL

}

priors

##
## ----  Sample model: ----------------------------------------------------------------
##
model_samples_obj <-  model_prep_obj$sample(   
                             n_burnin = 500,
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


# evals / sec [samp.] (parameters block, main only) =  33100
# [1] "Number of divergences =  69"
# Number of divergences =  69
# [1] "% divergences =  0.0165865384615385"
# % divergences =  0.0165865384615385
# List of 3
# 
# 0.0165865384615385 * 2 * 100
# 
# internal_obj <- model_samples_obj$internal_obj
# MCMC_params <- model_samples_obj$MCMC_params

##
## ----  Summarise + output results: -------------------------------------------------
##

# compute_main_params = TRUE
# compute_transformed_parameters = FALSE
# compute_generated_quantities = TRUE
# save_log_lik_trace = TRUE
# use_BayesMVP_for_faster_summaries = TRUE
# 
# package = "MetaOrdDTA"
# use_bridgestan = FALSE
# use_BayesMVP_for_faster_summaries = TRUE
# ##
# internal_obj = internal_obj
# MCMC_params  = MCMC_params
# ##
# model_results  = NULL
# init_object = NULL
# ##
# n_nuisance = 0


model_summary_and_trace_obj <- model_samples_obj$summary(
                                          compute_main_params = TRUE,
                                          compute_transformed_parameters = FALSE, 
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = TRUE,
                                          ##
                                          use_BayesMVP_for_faster_summaries = TRUE)

tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 200)
# tibble_tp <- model_summary_and_trace_obj$get_summary_transformed() %>% print(n = 100)
tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)
##
tibble_all <- rbind(tibble_main, tibble_gq)
# tibble_all <- rbind(tibble_main, tibble_tp, tibble_gq)
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_mean_overall_measure"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp_mean_overall_measure"))) %>% print(n = 1000)
##
# 
# 
# +try({  
#   
#     # tibble_main %>% print(n = 1000)
#     # tibble_gq   %>% print(n = 1000)
#     # 
#     # #### dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "LRpos")))
#     dplyr::filter(tibble_tp, (stringr::str_detect(parameter, "C"))) %>% print(n = 1000)
#     # dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "biv_equiv_C"))) %>% print(n = 1000)
#     # 
#     # dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "biv_equiv_C"))) %>% print(n = 1000)
#     # 
#     # stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
#     # Se <- stan_mod_samples$summary(c("Se"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
#     # 
#     # 
#     # stan_mod_samples$summary(c("C_mu"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
#     # 
#     # stan_mod_samples$summary(c("alpha"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
#     # 
#     # stan_mod_samples$summary(c("beta_SD"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
#     # 
#     # 
#     

#     ##
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "alpha"))) %>% print(n = 1000)
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
#     
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dev"))) %>% print(n = 1000)
#     
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "C"))) %>% print(n = 1000)
#     
#     
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_phi"))) %>% print(n = 1000)
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "prob_ord_mu"))) %>% print(n = 1000)
#     
#     
#     dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_phi"))) %>% print(n = 1000)
# 
# })


## 
## ----  Plots: ----------------------------------------------------------------------
##
# model_summary_and_trace_obj$plot_sROC()

# model_summary_and_trace_obj$plot_traces(param_string = "alpha")


{
        n_thr_max <- max(model_prep_obj$internal_obj$outs_data$n_thr)
        ##
        Se_median_array <- Sp_median_array <- Fp_median_array    <- array(dim = c(1, n_thr_max))
        Se_mean_array <- Sp_mean_array <- Fp_mean_array    <- array(dim = c(1, n_thr_max))
        ##
        Se_lower_array <- Sp_lower_array <- Fp_lower_array <- array(dim = c(1, n_thr_max))
        Se_upper_array <- Sp_upper_array <- Fp_upper_array <- array(dim = c(1, n_thr_max))
        ##
        Se_pred_lower_array <- Sp_pred_lower_array <- Fp_pred_lower_array <- array(dim = c(1, n_thr_max))
        Se_pred_upper_array <- Sp_pred_upper_array <- Fp_pred_upper_array <- array(dim = c(1, n_thr_max))
        
        # str(Se_array)
        
        Se <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se")) &
                                     (!(stringr::str_detect(parameter, "Se_pred"))))
        Sp <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp")) &
                              (!(stringr::str_detect(parameter, "Sp_pred"))))
        Fp  <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp")) &
                              (!(stringr::str_detect(parameter, "Fp_pred"))))
        ##
        Se_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_pred")))
        Sp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_pred")))
        Fp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_pred")))
        
        counter <- 1
        for (k in 1:n_thr_max) {
          for (t in 1:1) {
              # test_vec[counter] <- t
              ## Posterior medians (of pooled estimates):
              Se_median_array[t, k] <- 100 * Se$`50%`[counter]
              Sp_median_array[t, k] <- 100 * Sp$`50%`[counter]
              Fp_median_array[t, k] <- 100 * Fp$`50%`[counter]
              ## Posterior means (of pooled estimates):
              Se_mean_array[t, k] <- 100 * Se$mean[counter]
              Sp_mean_array[t, k] <- 100 * Sp$mean[counter]
              Fp_mean_array[t, k] <- 100 * Fp$mean[counter]
              ## Posterior lower 95% (of pooled estimates):
              Se_lower_array[t, k] <- 100 * Se$`2.5%`[counter]
              Sp_lower_array[t, k] <- 100 * Sp$`2.5%`[counter]
              Fp_lower_array[t, k] <- 100 * Fp$`2.5%`[counter]
              ## Posterior upper 95% (of pooled estimates):
              Se_upper_array[t, k] <- 100 * Se$`97.5%`[counter]
              Sp_upper_array[t, k] <- 100 * Sp$`97.5%`[counter]
              Fp_upper_array[t, k] <- 100 * Fp$`97.5%`[counter]
              ## Posterior lower prediction 95% (of pooled estimates):
              Se_pred_lower_array[t, k] <- 100 * Se_pred$`2.5%`[counter]
              Sp_pred_lower_array[t, k] <- 100 * Sp_pred$`2.5%`[counter]
              Fp_pred_lower_array[t, k] <- 100 * Fp_pred$`2.5%`[counter]
              ## Posterior upper prediction 95% (of pooled estimates):
              Se_pred_upper_array[t, k] <- 100 * Se_pred$`97.5%`[counter]
              Sp_pred_upper_array[t, k] <- 100 * Sp_pred$`97.5%`[counter]
              Fp_pred_upper_array[t, k] <- 100 * Fp_pred$`97.5%`[counter]
              ##
              counter <- counter + 1
          }
        }
}
 
# Se %>% print(n = 100) # Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
# Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
# Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
# Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()



  
  

  ##------------------------
  {
    t <- index_test_chosen_index - 1
    
    vec_index_inner <-  5:20
    vec_index_inner_2 <- vec_index_inner  - 0
  }
  
  
  {
    Se_abs_diffs_sim_DGM <- abs(1 * true_DGM_Se[vec_index_inner] - 1 * Se_mean_array[1, 1:n_thr][vec_index_inner_2])
    Sp_abs_diffs_sim_DGM  <- abs(1 * true_DGM_Sp[vec_index_inner] - 1 * Sp_mean_array[1, 1:n_thr][vec_index_inner_2])
    ##
    Se_mean_of_diffs_sim_DGM  <- mean(Se_abs_diffs_sim_DGM )
    Sp_mean_of_diffs_sim_DGM  <- mean(Sp_abs_diffs_sim_DGM )
    ##
    Se_sum_of_diffs_sim_DGM  <- sum(Se_abs_diffs_sim_DGM )
    Sp_sum_of_diffs_sim_DGM  <- sum(Sp_abs_diffs_sim_DGM )
    ##
    Se_max_of_diffs_sim_DGM  <- max(Se_abs_diffs_sim_DGM )
    Sp_max_of_diffs_sim_DGM  <- max(Sp_abs_diffs_sim_DGM )
    ##
    Se_median_of_diffs_sim_DGM  <- median(Se_abs_diffs_sim_DGM )
    Sp_median_of_diffs_sim_DGM  <- median(Sp_abs_diffs_sim_DGM )
  }
  
  
  # Sp_mean_array
  # true_DGM_Sp
  
  
  
  {
    message(paste("Se_mean_of_diffs_sim (using DGM true values) = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_mean_of_diffs_sim_DGM, 3))
    ##
    # message(paste("Se_median_of_diffs_sim (using DGM true values) = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Se_median_of_diffs_sim_DGM , 3))
    ##
    message(paste("Se_max_of_diffs_sim  (using DGM true values) =  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_max_of_diffs_sim_DGM , 3))
    ##
    ##
    ##
    message(paste("Sp_mean_of_diffs_sim   (using DGM true values) = = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_mean_of_diffs_sim_DGM, 3))
    ##
    # message(paste("Sp_median_of_diffs_sim   (using DGM true values) = = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Sp_median_of_diffs_sim_DGM , 3))
    ##
    message(paste("Sp_max_of_diffs_sim   (using DGM true values) =  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_max_of_diffs_sim_DGM , 3))
  }
  ##
  
  
  
  {
    Se_abs_diffs_sim_EMP <- abs(1 * true_Se_OVERALL_weighted[[t]][vec_index_inner] - 1 * Se_mean_array[1, 1:n_thr][vec_index_inner_2])
    Sp_abs_diffs_sim_EMP <- abs(1 * true_Sp_OVERALL_weighted[[t]][vec_index_inner] - 1 * Sp_mean_array[1, 1:n_thr][vec_index_inner_2])
    ##
    Se_mean_of_diffs_sim_EMP <- mean(Se_abs_diffs_sim_EMP)
    Sp_mean_of_diffs_sim_EMP <- mean(Sp_abs_diffs_sim_EMP)
    ##
    Se_sum_of_diffs_sim_EMP <- sum(Se_abs_diffs_sim_EMP)
    Sp_sum_of_diffs_sim_EMP <- sum(Sp_abs_diffs_sim_EMP)
    ##
    Se_max_of_diffs_sim_EMP <- max(Se_abs_diffs_sim_EMP)
    Sp_max_of_diffs_sim_EMP <- max(Sp_abs_diffs_sim_EMP)
    ##
    Se_median_of_diffs_sim_EMP <- median(Se_abs_diffs_sim_EMP)
    Sp_median_of_diffs_sim_EMP <- median(Sp_abs_diffs_sim_EMP)
  }
  
  
  
  
  
  {
    message(paste("Se_mean_of_diffs_sim (using Empirical true values) = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_mean_of_diffs_sim_EMP, 3))
    ##
    # message(paste("Se_median_of_diffs_sim (using Empirical true values) = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Se_median_of_diffs_sim_EMP, 3))
    ##
    message(paste("Se_max_of_diffs_sim  (using Empirical true values) =  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_max_of_diffs_sim_EMP, 3))
    ##
    ##
    ##
    message(paste("Sp_mean_of_diffs_sim   (using Empirical true values) = = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_mean_of_diffs_sim_EMP, 3))
    ##
    # message(paste("Sp_median_of_diffs_sim   (using Empirical true values) = = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Sp_median_of_diffs_sim_EMP, 3))
    ##
    message(paste("Sp_max_of_diffs_sim   (using Empirical true values) =  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_max_of_diffs_sim_EMP, 3))
  }
  
  
  
  
  {
    vec_index_other_1 <- setdiff(1:n_thr, vec_index_inner)
    vec_index_other_2 <- vec_index_other_1  - 0
  }
  
  
  {
    Se_abs_diffs_sim_DGM_ALL <- abs(1 * true_DGM_Se[vec_index_other_1] - 1 * Se_mean_array[1, 1:n_thr][vec_index_other_2])
    Sp_abs_diffs_sim_DGM_ALL <- abs(1 * true_DGM_Sp[vec_index_other_1] - 1 * Sp_mean_array[1, 1:n_thr][vec_index_other_2])
    ##
    Se_mean_of_diffs_sim_DGM_ALL <- mean(Se_abs_diffs_sim_DGM_ALL)
    Sp_mean_of_diffs_sim_DGM_ALL <- mean(Sp_abs_diffs_sim_DGM_ALL)
    ##
    Se_sum_of_diffs_sim_DGM_ALL <- sum(Se_abs_diffs_sim_DGM_ALL)
    Sp_sum_of_diffs_sim_DGM_ALL <- sum(Sp_abs_diffs_sim_DGM_ALL)
    ##
    Se_max_of_diffs_sim_DGM_ALL <- max(Se_abs_diffs_sim_DGM_ALL)
    Sp_max_of_diffs_sim_DGM_ALL <- max(Sp_abs_diffs_sim_DGM_ALL)
    ##
    Se_median_of_diffs_sim_DGM_ALL <- median(Se_abs_diffs_sim_DGM_ALL)
    Sp_median_of_diffs_sim_DGM_ALL <- median(Sp_abs_diffs_sim_DGM_ALL)
  }
  
  
  # Sp_mean_array
  # true_DGM_Sp
  
  
  
  {
    message(paste("Se_mean_of_diffs_sim (using DGM true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_mean_of_diffs_sim_DGM_ALL, 3))
    ##
    # message(paste("Se_median_of_diffs_sim (using DGM true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Se_median_of_diffs_sim_DGM_ALL, 3))
    ##
    message(paste("Se_max_of_diffs_sim  (using DGM true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_max_of_diffs_sim_DGM_ALL, 3))
    ##
    ##
    ##
    message(paste("Sp_mean_of_diffs_sim   (using DGM true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_mean_of_diffs_sim_DGM_ALL, 3))
    # ##
    # message(paste("Sp_median_of_diffs_sim   (using DGM true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Sp_median_of_diffs_sim_DGM_ALL, 3))
    ##
    message(paste("Sp_max_of_diffs_sim   (using DGM true values) [ALL thr]  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_max_of_diffs_sim_DGM_ALL, 3))
  }
  ##
  
  
  
  {
    Se_abs_diffs_sim <- abs(1 * true_Se_OVERALL_weighted[[t]][vec_index_other_1] - 1 * Se_mean_array[1, 1:n_thr][vec_index_other_2])
    Sp_abs_diffs_sim <- abs(1 * true_Sp_OVERALL_weighted[[t]][vec_index_other_1] - 1 * Sp_mean_array[1, 1:n_thr][vec_index_other_2])
    ##
    Se_mean_of_diffs_sim_EMP_ALL <- mean(Se_abs_diffs_sim)
    Sp_mean_of_diffs_sim_EMP_ALL <- mean(Sp_abs_diffs_sim)
    ##
    Se_sum_of_diffs_sim_EMP_ALL <- sum(Se_abs_diffs_sim)
    Sp_sum_of_diffs_sim_EMP_ALL <- sum(Sp_abs_diffs_sim)
    ##
    Se_max_of_diffs_sim_EMP_ALL <- max(Se_abs_diffs_sim)
    Sp_max_of_diffs_sim_EMP_ALL <- max(Sp_abs_diffs_sim)
    ##
    Se_median_of_diffs_sim_EMP_ALL <- median(Se_abs_diffs_sim)
    Sp_median_of_diffs_sim_EMP_ALL <- median(Sp_abs_diffs_sim)
  }
  
  
  
  
  
  {
    message(paste("Se_mean_of_diffs_sim (using Empirical true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_mean_of_diffs_sim_EMP_ALL, 3))
    ##
    # message(paste("Se_median_of_diffs_sim (using Empirical true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Se_median_of_diffs_sim_EMP_ALL, 3))
    ##
    message(paste("Se_max_of_diffs_sim  (using Empirical true values) [ALL thr]  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Se_max_of_diffs_sim_EMP_ALL, 3))
    ##
    ##
    ##
    message(paste("Sp_mean_of_diffs_sim   (using Empirical true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_mean_of_diffs_sim_EMP_ALL, 3))
    ##
    # message(paste("Sp_median_of_diffs_sim   (using Empirical true values) [ALL thr] = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(Sp_median_of_diffs_sim_EMP_ALL, 3))
    ##
    message(paste("Sp_max_of_diffs_sim   (using Empirical true values) [ALL thr]  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(Sp_max_of_diffs_sim_EMP_ALL, 3))
  }
  
  
  
  
  
  ##
  ## ---- COMBINED scores:
  ##
  {
    message(paste("Se_mean_of_diffs_sim (using DGM true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(0.25*Se_mean_of_diffs_sim_DGM + 0.75*Se_mean_of_diffs_sim_DGM_ALL, 3))
    ##
    # message(paste("Se_median_of_diffs_sim (using DGM true values) = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(0.25*Se_median_of_diffs_sim_DGM + 0.75*Se_median_of_diffs_sim_DGM_ALL, 3))
    ##
    message(paste("Se_max_of_diffs_sim  (using DGM true values) [combined weighted score]  = ")) ## , Se_sum_of_diffs_sim))
    print(signif(0.25*Se_max_of_diffs_sim_DGM + 0.75*Se_max_of_diffs_sim_DGM_ALL, 3))
    ##
    ##
    ##
    message(paste("Sp_mean_of_diffs_sim   (using DGM true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
    print(signif(0.25*Sp_mean_of_diffs_sim_DGM + 0.75*Sp_mean_of_diffs_sim_DGM_ALL, 3))
    ##
    # message(paste("Sp_median_of_diffs_sim   (using DGM true values) = = ")) ## , Se_sum_of_diffs_sim))
    # print(signif(0.25*Sp_median_of_diffs_sim_DGM + 0.75*Sp_median_of_diffs_sim_DGM_ALL, 3))
    ##
    message(paste("Sp_max_of_diffs_sim   (using DGM true values) [combined weighted score]   = ")) ## , Se_sum_of_diffs_sim))
    print(signif(0.25*Sp_max_of_diffs_sim_DGM + 0.75*Sp_max_of_diffs_sim_DGM_ALL, 3))
  }
  
  
  # {
  #   message(paste("Se_mean_of_diffs_sim (using Empirical true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
  #   print(signif(0.25*Se_mean_of_diffs_sim_EMP + 0.75*Se_mean_of_diffs_sim_EMP_ALL, 3))
  #   ##
  #   # message(paste("Se_median_of_diffs_sim (using Empirical true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
  #   # print(signif(0.25*Se_median_of_diffs_sim_EMP + 0.75*Se_median_of_diffs_sim_EMP_ALL, 3))
  #   ##
  #   message(paste("Se_max_of_diffs_sim  (using Empirical true values) [combined weighted score]  = ")) ## , Se_sum_of_diffs_sim))
  #   print(signif(0.25*Se_max_of_diffs_sim_EMP + 0.75*Se_max_of_diffs_sim_EMP_ALL, 3))
  #   ##
  #   ##
  #   ##
  #   message(paste("Sp_mean_of_diffs_sim   (using Empirical true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
  #   print(signif(0.25*Sp_mean_of_diffs_sim_EMP + 0.75*Sp_mean_of_diffs_sim_EMP_ALL, 3))
  #   ##
  #   # message(paste("Sp_median_of_diffs_sim   (using Empirical true values) [combined weighted score] = ")) ## , Se_sum_of_diffs_sim))
  #   # print(signif(0.25*Sp_median_of_diffs_sim_EMP + 0.75*Sp_median_of_diffs_sim_EMP_ALL, 3))
  #   ##
  #   message(paste("Sp_max_of_diffs_sim   (using Empirical true values) [combined weighted score]  = ")) ## , Se_sum_of_diffs_sim))
  #   print(signif(0.25*Sp_max_of_diffs_sim_EMP + 0.75*Sp_max_of_diffs_sim_EMP_ALL, 3))
  # }
  # 
  # 
  # 
  
  
  
  
  

{
    # x[[1]][, 10]
    print(paste("seed = ", seed))
    print(paste("model_parameterisation = ", model_parameterisation))
    print(paste("box_cox = ", box_cox))
    print(paste("softplus = ", softplus))
    print(paste("n_studies = ", n_studies))
    print(paste("bs_het_C = ", bs_het_C))
}


  
  
  
  
 t
 true_Sp_OVERALL_weighted[[t]]
 true_DGM_Sp
 Sp_mean_array


 t
 true_Se_OVERALL_weighted[[t]]
 true_DGM_Se
 Se_mean_array



 dplyr::filter(tibble_all, (stringr::str_detect(parameter, "surv_prob"))) %>% print(n = 1000)
 
 
 
 
 stan_data_list$x_with_missings
 stan_data_list$x
 stan_data_list$cutpoint_index
 stan_data_list$n_obs_cutpoints
 ##
 stan_data_list$cts_thr_values
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 



# 
# ####
# ####
# ####
# ppc_outs <- induced_Dirichlet_ppc_plot(method = "alpha",
#                            N = 5000,
#                            n_cat = n_thr + 1,
#                            other_args_list = list(use_log_alpha = FALSE,
#                                                   prior_alpha_mean = prior_alpha_mean,
#                                                   prior_alpha_sd = prior_alpha_SD,
#                                                   alpha_lb = alpha_lb,
#                                                   alpha_ub = NULL))
# 
# 
# ppc_outs <- induced_Dirichlet_ppc_plot(method = "kappa",
#                                        N = 5000,
#                                        n_cat = n_thr + 1,
#                                        other_args_list = list(use_log_kappa = FALSE,
#                                                               prior_kappa_mean = priors$prior_kappa_mean,
#                                                               prior_kappa_sd = priors$prior_kappa_SD,
#                                                               kappa_lb = priors$kappa_lb,
#                                                               kappa_ub = NULL, 
#                                                               prior_dirichlet_cat_means_alpha = priors$prior_dirichlet_phi[[1]]))
# 
# 
# 
# ## ppc_outs$log_alpha
# 
# 
# 

# 
# 
# # 
# # 
# 
# 
# true_Se_OVERALL_weighted[[4]] - 100 * Se_mean_array[1, 1:n_thr[4]]
# 
# true_Se_OVERALL_weighted[[1]] - 100 * Se_mean_array[1, 1:n_thr[1]]
# true_Se_OVERALL_weighted[[2]] - 100 * Se_mean_array[2, 1:n_thr[2]]
# true_Se_OVERALL_weighted[[3]] - 100 * Se_mean_array[3, 1:n_thr[3]]
# true_Se_OVERALL_weighted[[4]] - 100 * Se_mean_array[4, 1:n_thr[4]]
# #
# #
# # true_Sp_OVERALL_weighted[[1]] - 100 * Sp_mean_array[1, 1:n_thr[1]]
# # true_Sp_OVERALL_weighted[[2]] - 100 * Sp_mean_array[2, 1:n_thr[2]]
# # true_Sp_OVERALL_weighted[[3]] - 100 * Sp_mean_array[3, 1:n_thr[3]]
# # true_Sp_OVERALL_weighted[[4]] - 100 * Sp_mean_array[4, 1:n_thr[4]]
# #
# #
# 
# 
# 
# 
# ##
# ## ----  Plots: ----------------------------------------------------------------------
# ##
# # ##
# # ## MCMC trace plots:
# # ##
# # model_summary_and_trace_obj$plot_traces(param_string = c("Se", "Sp"))
# # ##
# # ## Densities:
# # ##
# # model_summary_and_trace_obj$plot_densities(param_string = c("Se", "Sp"))
# ##
# ## sROC plots:
# ##
# # model_summary_and_trace_obj$plot_sROC()
# ##
# df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100,
#                   Sp_true = true_Sp_OVERALL_weighted/100,
#                   Fp_true = (100 - true_Sp_OVERALL_weighted)/100)
# 
# 
# df_true
# ###
# model_summary_and_trace_obj$plot_sROC(df_true = df_true)
# 
# 
# tibble_gq <- model_summary_and_trace_obj$get_summary_generated_quantities()
# tibble_Se <- filter(tibble_gq, str_detect(parameter, "Se"))
# tibble_Sp <- filter(tibble_gq, str_detect(parameter, "Sp"))
# 
# 
# 
# if (model_parameterisation == "Jones") {
# 
#     Se_Jones <- head(tibble_Se$mean, n_thr)
#     Sp_Jones <- head(tibble_Sp$mean, n_thr)
# 
# } else {
# 
#     Se_Cerullo <- head(tibble_Se$mean, n_thr)
#     Sp_Cerullo <- head(tibble_Sp$mean, n_thr)
# 
# }
# 
# 
# mean(abs(100*df_true$Se_true - 100*Se_Jones))
# max(abs(100*df_true$Se_true  - 100*Se_Jones))
# sum(abs(100*df_true$Se_true  - 100*Se_Jones))
# ##
# mean(abs(100*df_true$Sp_true - 100*Sp_Jones))
# max(abs(100*df_true$Sp_true  - 100*Sp_Jones))
# sum(abs(100*df_true$Sp_true  - 100*Sp_Jones))
# ####
# mean(abs(100*df_true$Se_true - 100*Se_Cerullo))
# max(abs(100*df_true$Se_true  - 100*Se_Cerullo))
# sum(abs(100*df_true$Se_true  - 100*Se_Cerullo))
# ##
# mean(abs(100*df_true$Sp_true - 100*Sp_Cerullo))
# max(abs(100*df_true$Sp_true  - 100*Sp_Cerullo))
# sum(abs(100*df_true$Sp_true  - 100*Sp_Cerullo))
# ####
# ####
# ####
# mean(abs(100*df_true$Se_true - 100*Se_Cerullo)) - mean(abs(100*df_true$Se_true - 100*Se_Jones))
# max(abs(100*df_true$Se_true  - 100*Se_Cerullo)) - max(abs(100*df_true$Se_true  - 100*Se_Jones))
# sum(abs(100*df_true$Se_true  - 100*Se_Cerullo)) - sum(abs(100*df_true$Se_true  - 100*Se_Jones))
# ##
# mean(abs(100*df_true$Sp_true - 100*Sp_Cerullo)) - mean(abs(100*df_true$Sp_true - 100*Sp_Jones))
# max(abs(100*df_true$Sp_true  - 100*Sp_Cerullo)) - max(abs(100*df_true$Sp_true  - 100*Sp_Jones))
# sum(abs(100*df_true$Sp_true  - 100*Sp_Cerullo)) - sum(abs(100*df_true$Sp_true  - 100*Sp_Jones))
# 
# 
# #
# #
# # model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_mu_mean
# # model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_mu_SD
# # model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_SD_mean
# # model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_SD_SD
# # model_samples_obj$internal_obj$outs_data$stan_data_list$prior_beta_corr_LKJ
# #
# # model_samples_obj$advanced_model_options
# # model_samples_obj$basic_model_options
# #
# # model_samples_obj$priors
# #
# # model_samples_obj$outs_data$stan_data_list
# #
# # model_samples_obj$outs_data$stan_data_list
# #
# #
# # inits_for_ord_biv_fixed_thr_model
# #
# 
# 
# 
# #
# #
# {
# 
# 
# 
#         basic_model_options = list(
#           network = NULL,
#           cts = NULL,
#           prior_only = NULL
#         )
#         ####
#         advanced_model_options = list(
#           model_parameterisation = NULL,
#           random_thresholds = NULL,
#           Dirichlet_random_effects_type = NULL,
#           box_cox = NULL,
#           softplus = NULL
#         )
#         ####
#         other_advanced_options = list(
#           advanced_compile = NULL, ## default is false (aka a "basic" compile)
#           ##
#           force_recompile = NULL,
#           quiet = NULL,
#           compile = NULL,
#           ##
#           set_custom_CXX_CPP_flags = NULL,
#           CCACHE_PATH = NULL,
#           custom_cpp_user_header_file_path = NULL,
#           CXX_COMPILER_PATH = NULL,
#           CPP_COMPILER_PATH = NULL,
#           MATH_FLAGS = NULL,
#           FMA_FLAGS = NULL,
#           AVX_FLAGS = NULL,
#           THREAD_FLAGS = NULL
#         )
#         ##
#         ####
#         priors = NULL
#         ####
#         # init_lists_per_chain = NULL,
#         ####
#         MCMC_params = list(
#           seed = NULL,
#           n_superchains = NULL,
#           n_chains = NULL,
#           n_iter = NULL,
#           n_burnin = NULL,
#           adapt_delta = NULL,
#           max_treedepth = NULL,
#           metric_shape = NULL
#         )
# 
#         ##
#         ## ---- Main "internal_obj" list:
#         ##
#         internal_obj = list(
#           outs_data = NULL,
#           outs_stan_model_name = NULL,
#           outs_stan_compile = NULL,
#           outs_stan_init = NULL,
#           outs_stan_sampling = NULL,
#           ##
#           HMC_info = NULL, ## new
#           efficiency_info = NULL, ## new
#           summaries = NULL, ## new
#           traces = NULL ## new
#         )
# 
# 
# 
# 
#       debugging <- FALSE
#       ##
#       n_iter = 500
#       n_burnin = 500
#       ##
#       priors <- NULL
#       ##
#       x = x
#       n_chains = n_chains
#       ##
#       cts = FALSE
#       network = FALSE
#       ##
#       prior_only = FALSE
#       ##
#       softplus = TRUE
#       ##
#       ##
#       model_parameterisation = "Gatsonis"
#       random_thresholds = FALSE
#       Dirichlet_random_effects_type = "none"
#       box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
#       ##
#       init_lists_per_chain = NULL
# 
# 
#       basic_model_options$network <- network
#       basic_model_options$cts <- cts
#       basic_model_options$prior_only <- prior_only
#       ##
#       advanced_model_options$model_parameterisation <- model_parameterisation
#       advanced_model_options$random_thresholds <- random_thresholds
#       advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
#       advanced_model_options$box_cox <- box_cox
#       advanced_model_options$softplus <- softplus
#       ##
#       MCMC_params$n_chains <- n_chains
# 
# 
# 
# 
# 
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
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
