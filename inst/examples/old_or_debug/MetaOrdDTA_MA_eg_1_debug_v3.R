

 

os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")


# 
# 
# # # # # # 
# # # # ##
# rm(list = ls())
# # ##
# # .rs.restartR()  # In RStudio
# # #
# 
# 
# os <- .Platform$OS.type
# 
# 
# if (os == "unix") {
#   user_root_dir <- Sys.getenv("PWD")
# } else if (os == "windows") {
#   user_root_dir <- Sys.getenv("USERPROFILE")
# }
# local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
# #
# 
# 
#  {
#   ## First remove any possible package fragments:
#   ## Find user_pkg_install_dir:
#   user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
#   print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
#   ##
#   ## Find pkg_install_path + pkg_temp_install_path:
#   pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
#   pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
#   ##
#   ## Remove any (possible) MetaOrdDTA package fragments:
#   remove.packages("MetaOrdDTA")
#   unlink(pkg_install_path, recursive = TRUE, force = TRUE)
#   unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
# }
# 
# #
# 
# 
# #
# #### -------- ACTUAL (LOCAL) INSTALL:
# ## Document:
# devtools::clean_dll(local_pkg_dir)
# roxygen2::roxygenize(local_pkg_dir)
# ##
# ## Install (outer pkg):
# ##
# devtools::clean_dll(local_pkg_dir)
# devtools::install(local_pkg_dir,
#                   upgrade = "never",
#                   quick = TRUE)
# ##
# ## May need to restart R:
# ##
# .rs.restartR()  # In RStudio
# 
# # ?devtools::install
# # #
# # # #
# # # 

 


missing_thr <- TRUE
## missing_thr <- FALSE

n_chains <-  8 ##  parallel::detectCores() / 2
n_iter   <- round(500/(n_chains/4)) ; n_iter
n_burnin <- 500
 

##
n_studies <- 10
# n_studies <- 25
## n_studies <- 50
## n_studies <- 100


# ##
# n_sims <- 10
# ##

{
  if (n_studies == 10)   n_sims <- 10 ## need to estimate
  if (n_studies == 25)   n_sims <- 150 ## based on MCSE of 0.125% and 10 sims.
  if (n_studies == 50)   n_sims <- 50 ## based on MCSE of 0.125% and 10 sims - got 32 min. for Jones.
  if (n_studies == 100)  n_sims <- 10  ## ??? check
}



n_sims <- 5


{
  

    
    output_dir = MetaOrdDTA::metaorddta_output_path("simulation_results")
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Initialize storage for results
    # Set up summary storage
    summary_df <- data.frame(
      seed = integer(n_sims),
      n_studies = integer(n_sims),
      model_parameterisation = character(n_sims),
      box_cox = logical(n_sims),
      softplus = logical(n_sims),
      random_thresholds = logical(n_sims),
      ##
      mean_Se_abs_diffs_inner = numeric(n_sims),
      mean_Sp_abs_diffs_inner = numeric(n_sims),
      ##
      mean_Se_abs_diffs_outer = numeric(n_sims),
      mean_Sp_abs_diffs_outer = numeric(n_sims),
      ##
      mean_Se_abs_diffs_weighted = numeric(n_sims),
      mean_Sp_abs_diffs_weighted = numeric(n_sims),
      ##
      Se_mean_8  = numeric(n_sims),
      Se_mean_9  = numeric(n_sims),
      Se_mean_10 = numeric(n_sims),
      Se_mean_11 = numeric(n_sims),
      Se_mean_12 = numeric(n_sims),
      Se_mean_13 = numeric(n_sims),
      Se_mean_14 = numeric(n_sims),
      Se_mean_15 = numeric(n_sims),
      Se_mean_16 = numeric(n_sims),
      Se_mean_17 = numeric(n_sims),
      Se_mean_18 = numeric(n_sims),
      ##
      Sp_mean_8  = numeric(n_sims),
      Sp_mean_9  = numeric(n_sims),
      Sp_mean_10 = numeric(n_sims),
      Sp_mean_11 = numeric(n_sims),
      Sp_mean_12 = numeric(n_sims),
      Sp_mean_13 = numeric(n_sims),
      Sp_mean_14 = numeric(n_sims),
      Sp_mean_15 = numeric(n_sims),
      Sp_mean_16 = numeric(n_sims),
      Sp_mean_17 = numeric(n_sims),
      Sp_mean_18 = numeric(n_sims),
      ##
      stringsAsFactors = FALSE
    )
    
    
    network <- FALSE
    ##
    softplus <- FALSE
    # softplus <- TRUE


}



file_name <- "DTA_MA_Xu_FIXEDthr_daig_J.stan"
##
# use_custom_file <- TRUE
use_custom_file <- FALSE



# #### ----
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

{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
  # ##
  # inits <- list(beta_mu = c(-1, +1),
  #                                      beta_SD = c(0.01, 0.01),
  #                                      beta_z = array(0.01, dim = c(2, n_studies)),
  #                                      C_raw_vec = rep(-2.0, n_thr),
  #                                      C = seq(from = -2, to = 2, length = n_thr))
  # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
}
# 
{
  model_parameterisation = "R&G"
  box_cox <- FALSE
  cts <- FALSE
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
  # ##
  # inits <- list(beta_mu = c(-1, +1),
  #                                      beta_SD = c(0.01, 0.01),
  #                                      beta_z = array(0.01, dim = c(2, n_studies)),
  #                                      C_raw_vec = rep(-2.0, n_thr),
  #                                      C = seq(from = -2, to = 2, length = n_thr))
  # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
}
# 
{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  random_thresholds <-  TRUE
  Dirichlet_random_effects_type <- "alpha"
  # ##
  # inits <- list(beta_mu = c(-1, +1),
  #                                      beta_SD = c(0.01, 0.01),
  #                                      beta_z = array(0.01, dim = c(2, n_studies)),
  #                                      C_raw_vec = rep(-2.0, n_thr),
  #                                      C = seq(from = -2, to = 2, length = n_thr))
  # init_lists_per_chain <- replicate(n_chains, inits, simplify = FALSE)
}
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




sim_index <- 1




{
  
  
# Run simulations with different seeds
results_list <- c()
for (sim_index in 1:n_sims) {
  
  require(MetaOrdDTA)
  
  seed <- sim_index
  cat(sprintf("\n\n--- Running simulation %d/%d with seed = %d ---\n\n", sim_index, n_sims, seed))
  
  # Set seed for reproducibility
  set.seed(seed)
  

    {
          # ##
          # ## ---- Load NMA data:
          # ##
          setwd(local_pkg_dir)
          ##
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
          # x <- x_MA
          if (missing_thr == TRUE) x <- x_PHQ9_w_missing_thr
          else                     x <- x_PHQ9
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

# internal_obj$outs_data$stan_data_list$x
# internal_obj$outs_data$stan_data_list$n
# 
# any(is.na(unlist(internal_obj$outs_data$stan_data_list$n)))
init_lists_per_chain <- NULL
##
##
if (n_thr == 27) {
    vec_index_inner_thr <- 5:20
    n_inner_thr <- length(vec_index_inner_thr)
    vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner_thr)
    n_outer_thr <- length(vec_index_outer_thr)
    compute_sim_study_metrics <- 1
} else { 
    vec_index_inner_thr <- NULL
    compute_sim_study_metrics <- NULL
}
##
##
if (use_custom_file == TRUE) custom_file_name <- file_name
else  custom_file_name <- NULL
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
                  compute_sim_study_metrics = compute_sim_study_metrics,
                  vec_index_inner_thr = vec_index_inner_thr,
                  ##
                  custom_file_name = custom_file_name
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


 

init_lists_per_chain <- vector(length = n_chains, "list")
priors <- NULL


for (kk in 1:n_chains) {
  init_lists_per_chain[[kk]]$C_raw_vec <- list( rep(-3.0, n_thr), rep(-3.0, n_thr))
 # init_lists_per_chain[[kk]]$C <- list( seq(from = -2.0, to = 2.0, length = n_thr), 
  #                                      seq(from = -2.0, to = 2.0, length = n_thr))
}

if (model_parameterisation == "Xu") {
      
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
            priors$prior_beta_mu_SD   <- rep(1.0, 2)
            priors$prior_beta_SD_mean <- rep(0.0, 2)
            priors$prior_beta_SD_SD   <- rep(0.5, 2)
            ##
            
            ##
            # priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1) )
            alpha_priors <- generate_dual_ordinal_dirichlet_priors(n_cat = n_thr + 1,
                                                   nd_peak = 0.15,
                                                   d_peak = 0.30,
                                                   nd_max_alpha = 10,
                                                   d_max_alpha = 10,
                                                   nd_min_alpha = 1,
                                                   d_min_alpha = 1,
                                                   smoothness = 0.15)
            # ##
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
          priors$prior_beta_mu_SD   <- 1.0
          priors$prior_beta_SD_mean <- 0.0
          priors$prior_beta_SD_SD   <- 0.5
          ##
          priors$prior_raw_scale_mu_mean <- 0.0
          priors$prior_raw_scale_mu_SD   <- 1.0
          priors$prior_raw_scale_SD_mean <- 0.0
          priors$prior_raw_scale_SD_SD   <- 0.5 ## 0.25 ## on log-normal scale so should be small !!
          
          ##
          ## priors$prior_dirichlet_alpha <-  rep(2.0, n_thr + 1)
          priors$prior_dirichlet_alpha <- generate_ordinal_dirichlet_prior(n_cat = n_thr + 1,
                                                                           max_alpha = 5, 
                                                                           min_alpha = 1,
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
        priors$prior_raw_scale_SD_SD   <- rep(0.5, 2)
        ##
        priors$prior_boxcox_lambda_mean <- 0.00
        priors$prior_boxcox_lambda_SD   <- 1.00
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

model_prep_obj$outs_stan_model_name

##
## ----  Sample model: ----------------------------------------------------------------
##
model_samples_obj <-  model_prep_obj$sample(   
                             n_burnin = n_burnin,
                             n_iter   = n_iter,
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
# Se_mean_overall_measure <- dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_mean_overall_measure"))) %>% print(n = 1000)
# Sp_mean_overall_measure <- dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp_mean_overall_measure"))) %>% print(n = 1000)
# ##
# Se_mean_overall_measure_mean <- Se_mean_overall_measure$mean*100
# Sp_mean_overall_measure_mean <- Sp_mean_overall_measure$mean*100
##

tibble_gq %>% print(n = 100)

dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp"))) %>% print(n = 1000)


model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
##
## 
## ----  Plots: ----------------------------------------------------------------------
##
# model_summary_and_trace_obj$plot_sROC()

# model_summary_and_trace_obj$plot_traces(param_string = "alpha")


{
        n_thr_max <- max(model_prep_obj$internal_obj$outs_data$n_thr)
        ##
        Se_median_vec <- Sp_median_vec <- Fp_median_vec <- c()
        Se_mean_vec   <- Sp_mean_vec   <- Fp_mean_vec   <- c()
        ##
        Se_lower_vec <- Sp_lower_vec <- Fp_lower_vec <- c()
        Se_upper_vec <- Sp_upper_vec <- Fp_upper_vec <- c()
        ##
        Se_pred_lower_vec <- Sp_pred_lower_vec <- Fp_pred_lower_vec <- c()
        Se_pred_upper_vec <- Sp_pred_upper_vec <- Fp_pred_upper_vec <- c()
        
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
              Se_median_vec[k] <- 100 * Se$`50%`[counter]
              Sp_median_vec[k] <- 100 * Sp$`50%`[counter]
              Fp_median_vec[k] <- 100 * Fp$`50%`[counter]
              ## Posterior means (of pooled estimates):
              Se_mean_vec[k] <- 100 * Se$mean[counter]
              Sp_mean_vec[k] <- 100 * Sp$mean[counter]
              Fp_mean_vec[k] <- 100 * Fp$mean[counter]
              ## Posterior lower 95% (of pooled estimates):
              Se_lower_vec[k] <- 100 * Se$`2.5%`[counter]
              Sp_lower_vec[k] <- 100 * Sp$`2.5%`[counter]
              Fp_lower_vec[k] <- 100 * Fp$`2.5%`[counter]
              ## Posterior upper 95% (of pooled estimates):
              Se_upper_vec[k] <- 100 * Se$`97.5%`[counter]
              Sp_upper_vec[k] <- 100 * Sp$`97.5%`[counter]
              Fp_upper_vec[k] <- 100 * Fp$`97.5%`[counter]
              ## Posterior lower prediction 95% (of pooled estimates):
              Se_pred_lower_vec[k] <- 100 * Se_pred$`2.5%`[counter]
              Sp_pred_lower_vec[k] <- 100 * Sp_pred$`2.5%`[counter]
              Fp_pred_lower_vec[k] <- 100 * Fp_pred$`2.5%`[counter]
              ## Posterior upper prediction 95% (of pooled estimates):
              Se_pred_upper_vec[k] <- 100 * Se_pred$`97.5%`[counter]
              Sp_pred_upper_vec[k] <- 100 * Sp_pred$`97.5%`[counter]
              Fp_pred_upper_vec[k] <- 100 * Fp_pred$`97.5%`[counter]
              ##
              counter <- counter + 1
          }
        }
}
 
# Se %>% print(n = 100) # Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
# Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
# Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
# Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()


mean(abs(Se_mean_vec[5:20] - true_DGM_Se[5:20]))
mean(abs(Sp_mean_vec[5:20] - true_DGM_Sp[5:20]))
  
  

  ##------------------------
  {
     vec_index_all         <-  1:n_thr
     ##
     if (n_thr == 27) {
        vec_index_inner <-  8:18 ## most clinically important
        vec_index_outer <- setdiff(vec_index_all, vec_index_inner) ## less clinically important
        ##
        ## ---- inner thr:
        ##
        Se_abs_diffs_inner <- abs(true_DGM_Se[vec_index_inner] - Se_mean_vec[vec_index_inner])
        Sp_abs_diffs_inner <- abs(true_DGM_Sp[vec_index_inner] - Sp_mean_vec[vec_index_inner])
        mean_Se_abs_diffs_inner <- mean(Se_abs_diffs_inner)
        mean_Sp_abs_diffs_inner <- mean(Sp_abs_diffs_inner)
        ##
        ## ---- Outer thr:
        ##
        Se_abs_diffs_outer <- abs(true_DGM_Se[vec_index_outer] - Se_mean_vec[vec_index_outer])
        Sp_abs_diffs_outer <- abs(true_DGM_Sp[vec_index_outer] - Sp_mean_vec[vec_index_outer])
        mean_Se_abs_diffs_outer <- mean(Se_abs_diffs_outer)
        mean_Sp_abs_diffs_outer <- mean(Sp_abs_diffs_outer)
        ##
        ## ---- weighted score:
        ##
        weight <- 0.80
        mean_Se_abs_diffs_weighted <- weight*mean_Se_abs_diffs_inner + (1.0 - weight)*mean_Se_abs_diffs_outer
        mean_Sp_abs_diffs_weighted <- weight*mean_Sp_abs_diffs_inner + (1.0 - weight)*mean_Sp_abs_diffs_outer
        ##
        ## ---- Thr 11:
        ##
        Se_abs_diffs_inner <- abs(true_DGM_Se[vec_index_inner] - Se_mean_vec[11])
        Sp_abs_diffs_inner <- abs(true_DGM_Sp[vec_index_inner] - Sp_mean_vec[11])
        mean_Se_abs_diffs_inner <- mean(Se_abs_diffs_inner)
        mean_Sp_abs_diffs_inner <- mean(Sp_abs_diffs_inner)
        ##
     } else { 
       ### BOOKMARK
     }
    
        # vec_index_other_2 <- vec_index_other_1  - 0
        ##
        # TRUE_DGM_Se_mean_overall_measure <- 0.75*mean(true_DGM_Se[vec_index_inner]) + 0.25*mean(true_DGM_Se[vec_index_outer])
        # TRUE_DGM_Sp_mean_overall_measure <- 0.75*mean(true_DGM_Sp[vec_index_inner]) + 0.25*mean(true_DGM_Sp[vec_index_outer])
        ##
        # message(paste("Se - weighted measure across thresholds = "))
        # print(signif(Se_mean_overall_measure_mean), 3)
        ##
        # message(paste("Sp - weighted measure across thresholds = "))
        # print(signif(Sp_mean_overall_measure_mean), 3)
        message(paste("-------------------------------------------------------- ")) ## , Se_sum_of_diffs_sim))
        ##
        ## ---- Diffs from true values (inner):
        ##
        message(paste("Se_mean_of_diffs_sim (using DGM true values) - inner thr. = ")) ## , Se_sum_of_diffs_sim))
        print(signif(mean_Se_abs_diffs_inner), 3)
        ##
        message(paste("Sp_mean_of_diffs_sim (using DGM true values) - inner thr. = ")) ## , Se_sum_of_diffs_sim))
        print(signif(mean_Sp_abs_diffs_inner), 3)
        message(paste("-------------------------------------------------------- ")) ## , Se_sum_of_diffs_sim))
        ##
        ## ---- Diffs from true values (weighted):
        ##
        message(paste("Se_mean_of_diffs_sim (using DGM true values) - weighted score = ")) ## , Se_sum_of_diffs_sim))
        print(signif(mean_Se_abs_diffs_weighted), 3)
        ##
        message(paste("Sp_mean_of_diffs_sim (using DGM true values) - weighted score = ")) ## , Se_sum_of_diffs_sim))
        print(signif(mean_Sp_abs_diffs_weighted), 3)
        message(paste("-------------------------------------------------------- ")) ## , Se_sum_of_diffs_sim))
        {
          # x[[1]][, 10]
          print(paste("seed = ", seed))
          print(paste("model_parameterisation = ", model_parameterisation))
          print(paste("box_cox = ", box_cox))
          print(paste("softplus = ", softplus))
          print(paste("n_studies = ", n_studies))
          # print(paste("bs_het_C = ", bs_het_C))
        }
  }


 

    
    ##
    ## ---- Store results for this simulation:
    ##
    results_list[[sim_index]] <- list(
      seed = seed,
      n_studies = n_studies,
      model_parameterisation = model_parameterisation,
      box_cox = box_cox,
      softplus = softplus,
      random_thresholds = random_thresholds,
      ##
      true_DGM_Se = true_DGM_Se,
      true_DGM_Sp = true_DGM_Sp,
      ##
      Se_mean_vec = Se_mean_vec,
      Sp_mean_vec = Sp_mean_vec,
      ##
      mean_Se_abs_diffs_inner = mean_Se_abs_diffs_inner,
      mean_Sp_abs_diffs_inner = mean_Sp_abs_diffs_inner,
      ##
      mean_Se_abs_diffs_outer = mean_Se_abs_diffs_outer,
      mean_Sp_abs_diffs_outer = mean_Sp_abs_diffs_outer,
      ##
      mean_Se_abs_diffs_weighted = mean_Se_abs_diffs_weighted,
      mean_Sp_abs_diffs_weighted = mean_Sp_abs_diffs_weighted,
      ##
      Se_mean_8  = Se_mean_vec[8],
      Se_mean_9  = Se_mean_vec[9],
      Se_mean_10 = Se_mean_vec[10],
      Se_mean_11 = Se_mean_vec[11],
      Se_mean_12 = Se_mean_vec[12],
      Se_mean_13 = Se_mean_vec[13],
      Se_mean_14 = Se_mean_vec[14],
      Se_mean_15 = Se_mean_vec[15],
      Se_mean_16 = Se_mean_vec[16],
      Se_mean_17 = Se_mean_vec[17],
      Se_mean_18 = Se_mean_vec[18],
      ##
      Sp_mean_8  = Sp_mean_vec[8],
      Sp_mean_9  = Sp_mean_vec[9],
      Sp_mean_10 = Sp_mean_vec[10],
      Sp_mean_11 = Sp_mean_vec[11],
      Sp_mean_12 = Sp_mean_vec[12],
      Sp_mean_13 = Sp_mean_vec[13],
      Sp_mean_14 = Sp_mean_vec[14],
      Sp_mean_15 = Sp_mean_vec[15],
      Sp_mean_16 = Sp_mean_vec[16],
      Sp_mean_17 = Sp_mean_vec[17],
      Sp_mean_18 = Sp_mean_vec[18]
     ##  model_summary = model_summary_and_trace_obj  # Optional: save full summary object
    )
    ##
    ## ---- Update summary dataframe:
    ##
    create_summary_df <- function(n_sims, min_k, max_k) {
      # Create vector of column names
      base_cols <- c(
        "seed", "n_studies", "model_parameterisation", "box_cox", 
        "softplus", "random_thresholds",
        "mean_Se_abs_diffs_inner", "mean_Sp_abs_diffs_inner",
        "mean_Se_abs_diffs_outer", "mean_Sp_abs_diffs_outer",
        "mean_Se_abs_diffs_weighted", "mean_Sp_abs_diffs_weighted"
      )
      
      # Create empty vectors for threshold-specific column names
      se_cols <- character(max_k - min_k + 1)
      sp_cols <- character(max_k - min_k + 1)
      
      # Fill the threshold column names
      for (i in 1:length(se_cols)) {
        k <- min_k + i - 1
        se_cols[i] <- paste0("Se_mean_", k)
        sp_cols[i] <- paste0("Sp_mean_", k)
      }
      
      # Combine all column names
      all_cols <- c(base_cols, se_cols, sp_cols)
      
      # Create empty data frame with the right structure
      df <- data.frame(matrix(NA, nrow = n_sims, ncol = length(all_cols)))
      names(df) <- all_cols
      
      return(df)
    }
    
    # Usage example:
    summary_df <- create_summary_df(n_sims = n_sims, min_k = 5, max_k = 18)
    ##
    # Create base values for the row
    row_data <- list(
      seed = seed,
      n_studies = n_studies,
      model_parameterisation = model_parameterisation,
      box_cox = box_cox,
      softplus = softplus,
      random_thresholds = random_thresholds,
      ##
      mean_Se_abs_diffs_inner = mean_Se_abs_diffs_inner,
      mean_Sp_abs_diffs_inner = mean_Sp_abs_diffs_inner,
      ##
      mean_Se_abs_diffs_outer = mean_Se_abs_diffs_outer,
      mean_Sp_abs_diffs_outer = mean_Sp_abs_diffs_outer,
      ##
      mean_Se_abs_diffs_weighted = mean_Se_abs_diffs_weighted,
      mean_Sp_abs_diffs_weighted = mean_Sp_abs_diffs_weighted
    )
    
    # Add threshold-specific columns in a loop
    for (k in 5:18) {
      # Add sensitivity for this threshold
      row_data[[paste0("Se_mean_", k)]] <- Se_mean_vec[k]
      
      # Add specificity for this threshold
      row_data[[paste0("Sp_mean_", k)]] <- Sp_mean_vec[k]
    }
    
    # Assign all data to the row
    summary_df[sim_index, ] <- row_data
    
    ##
    ## ---- Save results for this simulation:
    ##
    file_name_string <- paste0( "n_studies_", n_studies, 
                                "param_", model_parameterisation, 
                                "soft_", softplus,
                                "rand_thr_", random_thresholds
                              )
    
    sim_result_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string))
    saveRDS(results_list[[sim_index]], sim_result_file)
    ##
    ## ---- Report progress:
    ##
    cat(sprintf("\n--- Simulation %d/%d completed ---\n", sim_index, n_sims))
    # cat(sprintf("Se difference from true value: %.3f\n", Se_diff))
    # cat(sprintf("Sp difference from true value: %.3f\n", Sp_diff))
  
  
  } ### end of simulations (N_{sim}) loop !!

  # Calculate summary statistics across ALL simulations
  summary_stats <- list(
    n_sims = n_sims,
    n_studies = n_studies,
    model_parameterisation = model_parameterisation,
    box_cox = box_cox,
    softplus = softplus,
    random_thresholds = random_thresholds,
    # ##
    # mean_Se = mean(summary_df$Se_mean_overall_measure),
    # median_Se = median(summary_df$Se_mean_overall_measure),
    # max_Se = max(summary_df$Se_mean_overall_measure),
    # ##
    # mean_Sp = mean(summary_df$Sp_mean_overall_measure),
    # median_Sp = median(summary_df$Sp_mean_overall_measure),
    # max_Sp = max(summary_df$Sp_mean_overall_measure),
    # ##
    # Se_sd = sd(summary_df$Se_mean_overall_measure),
    # Sp_sd = sd(summary_df$Sp_mean_overall_measure),
    ##
    mean_Se_diffs_inner   = mean(summary_df$mean_Se_abs_diffs_inner),
    median_Se_diffs_inner = median(summary_df$mean_Se_abs_diffs_inner),
    max_Se_diffs_inner    = max(summary_df$mean_Se_abs_diffs_inner),
    ##
    mean_Sp_diffs_inner   = mean(summary_df$mean_Sp_abs_diffs_inner),
    median_Sp_diffs_inner = median(summary_df$mean_Sp_abs_diffs_inner),
    max_Sp_diffs_inner    = max(summary_df$mean_Sp_abs_diffs_inner),
    ##
    Se_diffs_sd_inner = sd(summary_df$mean_Se_abs_diffs_inner),
    Sp_diffs_sd_inner = sd(summary_df$mean_Sp_abs_diffs_inner),
    ##
    ## ----
    ##
    mean_Se_diffs_outer   = mean(summary_df$mean_Se_abs_diffs_outer),
    median_Se_diffs_outer = median(summary_df$mean_Se_abs_diffs_outer),
    max_Se_diffs_outer    = max(summary_df$mean_Se_abs_diffs_outer),
    ##
    mean_Sp_diffs_outer   = mean(summary_df$mean_Sp_abs_diffs_outer),
    median_Sp_diffs_outer = median(summary_df$mean_Sp_abs_diffs_outer),
    max_Sp_diffs_outer    = max(summary_df$mean_Sp_abs_diffs_outer),
    ##
    Se_diffs_sd_outer = sd(summary_df$mean_Se_abs_diffs_outer),
    Sp_diffs_sd_outer = sd(summary_df$mean_Sp_abs_diffs_outer),
    ##
    ## ----
    ##
    mean_Se_diffs_weighted   = mean(summary_df$mean_Se_abs_diffs_weighted),
    median_Se_diffs_weighted = median(summary_df$mean_Se_abs_diffs_weighted),
    max_Se_diffs_weighted    = max(summary_df$mean_Se_abs_diffs_weighted),
    ##
    mean_Sp_diffs_weighted   = mean(summary_df$mean_Sp_abs_diffs_weighted),
    median_Sp_diffs_weighted = median(summary_df$mean_Sp_abs_diffs_weighted),
    max_Sp_diffs_weighted    = max(summary_df$mean_Sp_abs_diffs_weighted),
    ##
    Se_diffs_sd_weighted = sd(summary_df$mean_Se_abs_diffs_weighted),
    Sp_diffs_sd_weighted = sd(summary_df$mean_Sp_abs_diffs_weighted),
    ##
    summary_df = summary_df
  )
  
  
  # Loop through thresholds to add statistics
  for (k in 5:18) {
    # Sensitivity statistics
    se_col <- paste0("Se_mean_", k)
    if (se_col %in% names(summary_df)) {
      summary_stats[[paste0("mean_Se_", k)]] <- mean(summary_df[[se_col]])
      summary_stats[[paste0("median_Se_", k)]] <- median(summary_df[[se_col]])
      summary_stats[[paste0("max_Se_", k)]] <- max(summary_df[[se_col]])
      summary_stats[[paste0("SD_Se_", k)]] <- sd(summary_df[[se_col]])
    }
    
    # Specificity statistics
    sp_col <- paste0("Sp_mean_", k)
    if (sp_col %in% names(summary_df)) {
      summary_stats[[paste0("mean_Sp_", k)]] <- mean(summary_df[[sp_col]])
      summary_stats[[paste0("median_Sp_", k)]] <- median(summary_df[[sp_col]])
      summary_stats[[paste0("max_Sp_", k)]] <- max(summary_df[[sp_col]])
      summary_stats[[paste0("SD_Sp_", k)]] <- sd(summary_df[[sp_col]])
    }
  }
  
  
  # Save overall summary
  summary_file <- file.path(getwd(), output_dir, paste0(file_name_string, ".RDS"))
  saveRDS(summary_stats, summary_file)
  
  # Also save as CSV for easy viewing
  summary_csv <- file.path(getwd(), output_dir, paste0(file_name_string, ".csv"))
  write.csv(summary_df, summary_csv, row.names = FALSE)
  

  
  print(summary_stats)
  
  
  }

# Print summary:
print_se_sp_stats <- function(k, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir) {
  cat(paste("--------------------------------------------------"))
  cat(sprintf("\nSe (@thr = %d):\n", k))
  cat(sprintf("  Mean: %s\n", formatC(summary_stats[[paste0("mean_Se_", k)]], digits=3, format="g")))
  cat(sprintf("  SD: %s\n", formatC(summary_stats[[paste0("SD_Se_", k)]], digits=3, format="g")))
  cat(sprintf("  Bias: %s\n", formatC(summary_stats[[paste0("mean_Se_", k)]] - true_DGM_Se[k], digits=3, format="g")))
  ##
  cat(sprintf("\nSp (@thr = %d):\n", k))
  cat(sprintf("  Mean: %s\n", formatC(summary_stats[[paste0("mean_Sp_", k)]], digits=3, format="g")))
  cat(sprintf("  SD: %s\n", formatC(summary_stats[[paste0("SD_Sp_", k)]], digits=3, format="g")))
  cat(sprintf("  Bias: %s\n", formatC(summary_stats[[paste0("mean_Sp_", k)]] - true_DGM_Sp[k], digits=3, format="g")))
  ##
  cat(sprintf("\nResults saved to: %s\n", output_dir))
}
##
{
      cat("\n\n========= SIMULATION STUDY SUMMARY =========\n")
      ##
      cat(sprintf("Number of simulations: %d\n", n_sims))
      cat(sprintf("Number of studies: %d\n", n_studies))
      cat(sprintf("Model parameterisation: %s\n", model_parameterisation))
      cat(sprintf("Box-Cox: %s\n", box_cox))
      cat(sprintf("Softplus: %s\n", softplus))
      cat(sprintf("Random thresholds: %s\n", random_thresholds))
      ##
      ## ---- Se / Sp:
      ##
      print_se_sp_stats(5, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(6, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(7, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(8, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(9, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(10, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(11, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(12, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(13, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(14, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(15, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(16, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(17, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      print_se_sp_stats(18, summary_stats, true_DGM_Se, true_DGM_Sp, output_dir)
      # ##
      # ## ---- Diffs:
      # ##
      # cat(paste("--------------------------------------------------"))
      # cat("\nSensitivity abs. differences from true values:\n")
      # cat(sprintf("  Mean: %s\n", formatC(summary_stats$mean_Se_diff, digits=3, format="g")))
      # # cat(sprintf("  Median: %s\n", formatC(summary_stats$median_Se_diff, digits=3, format="g")))
      # cat(sprintf("  Maximum: %s\n", formatC(summary_stats$max_Se_diff, digits=3, format="g")))
      # cat(sprintf("  SD: %s\n", formatC(summary_stats$Se_diff_sd, digits=3, format="g")))
      # ##
      # cat("\nSpecificity abs. differences from true values:\n")
      # cat(sprintf("  Mean: %s\n", formatC(summary_stats$mean_Sp_diff, digits=3, format="g")))
      # # cat(sprintf("  Median: %s\n", formatC(summary_stats$median_Sp_diff, digits=3, format="g")))
      # cat(sprintf("  Maximum: %s\n", formatC(summary_stats$max_Sp_diff, digits=3, format="g")))
      # cat(sprintf("  SD: %s\n", formatC(summary_stats$Sp_diff_sd, digits=3, format="g")))
      # cat("\nResults saved to: %s\n", output_dir)
}

  
 
  for (k in 8:17) {
    summary_stats[[paste0("Bias_Se_", k)]] <- summary_stats[[paste0("mean_Se_", k)]] - true_DGM_Se[k]
  }
  
  ##
  ## ---- 
  ##
  compute_average_bias <- function(summary_stats, 
                                   true_DGM_Se, 
                                   true_DGM_Sp, 
                                   min_k = 8,
                                   max_k = 17) {
    
    # Initialize sums
    total_bias_Se <- 0
    total_bias_Sp <- 0
    total_SD_Se <- 0
    total_SD_Sp <- 0
    count <- 0
    
    length_vec <- max_k - min_k + 1
    ##
    bias_Se_vec <- c()
    bias_Sp_vec <- c()
    SD_Se_vec <- c()
    SD_Sp_vec <- c()
    
    # Compute bias for each k and add to total
    for (k in min_k:max_k) {
      
          ##
          ## ---- Se:
          ##
          bias_Se_vec[k - min_k + 1] <- abs(summary_stats[[paste0("mean_Se_", k)]] - true_DGM_Se[k])
          summary_stats[[paste0("Bias_Se_", k)]] <- bias_Se_vec[k - min_k + 1]
          ##
          SD_Se_vec[k - min_k + 1] <- summary_stats[[paste0("SD_Se_", k)]]
          ##
          ## ---- Sp
          ##
          bias_Sp_vec[k - min_k + 1] <- abs(summary_stats[[paste0("mean_Sp_", k)]] - true_DGM_Sp[k])
          summary_stats[[paste0("Bias_Sp_", k)]] <- bias_Sp_vec[k - min_k + 1]
          ##
          SD_Sp_vec[k - min_k + 1] <- summary_stats[[paste0("SD_Sp_", k)]]
          ##
          
    }
    
    # Calculate average biases
    mean_bias_Se <- mean(bias_Se_vec, na.rm = TRUE)
    max_bias_Se <- max(bias_Se_vec, na.rm = TRUE)
    mean_SD_Se   <- mean(SD_Se_vec, na.rm = TRUE)
    max_SD_Se   <- max(SD_Se_vec, na.rm = TRUE)
    ##
    mean_bias_Sp <- mean(bias_Sp_vec, na.rm = TRUE)
    max_bias_Sp <- max(bias_Sp_vec, na.rm = TRUE)
    mean_SD_Sp   <- mean(SD_Sp_vec, na.rm = TRUE)
    max_SD_Sp   <- max(SD_Sp_vec, na.rm = TRUE)
    
    # Return results
    return(list(
          avg_bias_Se = avg_bias_Se,
          max_bias_Se = max_bias_Se,
          avg_SD_Se = avg_SD_Se,
          max_SD_Se = max_SD_Se,
          ##
          avg_bias_Sp = avg_bias_Sp,
          max_bias_Sp = max_bias_Sp,
          avg_SD_Sp = avg_SD_Sp,
          max_SD_Sp = max_SD_Sp,
          ##
          summary_stats = summary_stats  # Return updated summary_stats with Bias values
    ))
  }
  



  {
        ##
        ## ---- Compute average bias for thresholds 8 through 17
        ##
        {
            result <- compute_average_bias(summary_stats, true_DGM_Se, true_DGM_Sp, min_k = 8, max_k = 18)
            ##
            message(paste("--------------------------------------------------"))
            ##
            cat(sprintf("Mean Se bias: %s\n", formatC(result$avg_bias_Se, digits=3, format="g")))
            cat(sprintf("Max Se bias: %s\n", formatC(result$max_bias_Se, digits=3, format="g")))
            cat(sprintf("Mean Se SD: %s\n", formatC(result$avg_SD_Se, digits=3, format="g")))
            cat(sprintf("Max Se SD: %s\n", formatC(result$max_SD_Se, digits=3, format="g")))
            ##
            message(paste("--------------------------------------------------"))
            ##
            cat(sprintf("Mean Sp bias: %s\n", formatC(result$avg_bias_Sp, digits=3, format="g")))
            cat(sprintf("Max Sp bias: %s\n", formatC(result$max_bias_Sp, digits=3, format="g")))
            cat(sprintf("Mean Sp SD: %s\n", formatC(result$avg_SD_Sp, digits=3, format="g")))
            cat(sprintf("Max Sp SD: %s\n", formatC(result$max_SD_Sp, digits=3, format="g")))
            ##
            # The function also returns the updated summary_stats with all the Bias values
            summary_stats <- result$summary_stats
            ##
            message(paste("--------------------------------------------------"))
            ##
            print(model_parameterisation)
        }
        
        t
        
        ##
        ## ----
        ##
        calc_min_sim_size <- function(summary_stats, 
                                      min_MCSE_pct,
                                      min_k = 8, 
                                      max_k = 18) {
          
                min_N_sim <- c()
                index <- 1
                
                for (k in min_k:max_k) {
                  # Calculate minimum N required for this threshold
                  sd_key <- paste0("SD_Se_", k)
                  if (sd_key %in% names(summary_stats)) {
                    min_N_sim[index] <- ceiling(1 + 0.5*((summary_stats[[sd_key]]/min_MCSE_pct)^2))
                    index <- index + 1
                  } else {
                    warning(paste("Key", sd_key, "not found in summary_stats"))
                  }
                }
                
                # Return both the vector and the maximum value
                return(list(
                 ## min_N_sim = min_N_sim,
                  min_required = min(min_N_sim, na.rm = TRUE),
                  avg_required = mean(min_N_sim, na.rm = TRUE),
                  max_required = max(min_N_sim, na.rm = TRUE)
                ))
          
        }
        ##
        {
          message("if min_MCSE_pct = 0.25%:")
          print(calc_min_sim_size(summary_stats = summary_stats, 
                            min_MCSE_pct = 0.25, 
                            min_k = 8,
                            max_k = 17))
          ##
          message("if min_MCSE_pct = 0.125%:")
          print(calc_min_sim_size(summary_stats = summary_stats, 
                            min_MCSE_pct = 0.125, 
                            min_k = 8,
                            max_k = 17))
          ##
          message("if min_MCSE_pct = 0.05%:")
          print(calc_min_sim_size(summary_stats = summary_stats, 
                            min_MCSE_pct = 0.05, 
                            min_k = 8,
                            max_k = 17))
        }
  }
    
  
 # true_Sp_OVERALL_weighted[[t]]
 # true_DGM_Sp
 # Sp_mean_array
 # 
 # 
 # t
 # true_Se_OVERALL_weighted[[t]]
 # true_DGM_Se
 # Se_mean_array
 # 
 # 
 # 
 # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "surv_prob"))) %>% print(n = 1000)
 # 
 # 
 # 
 # 
 # stan_data_list$x_with_missings
 # stan_data_list$x
 # stan_data_list$cutpoint_index
 # stan_data_list$n_obs_cutpoints
 # ##
 # stan_data_list$cts_thr_values
 # 
 # 
 
 
 
 
 
 
  
 
 
 
 