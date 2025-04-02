

 
 
##
rm(list = ls())
# ##
# .rs.restartR()  # In RStudio
# # 



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




# ##
# ## ---- Load NMA data:
# ##
# setwd(local_pkg_dir)
source(file.path(getwd(), "inst", "examples", "NMA_missing_thr_prep_data.R"))
data <- readRDS(file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))
x_NMA <- data$x_NMA

x_MA <- list()
for (c in 1:2) {
  x_MA[[c]] <- x_NMA[[c]][[4]]
}
x <- x_MA
n_thr <- ncol(x[[1]])

# x <- x_NMA
# indicator_index_test_in_study <- data$indicator_index_test_in_study
# ##
# n_studies <- nrow(x_NMA[[1]][[1]])
# n_index_tests <- 4##  length(x_NMA[[1]])
# 
# subset <- TRUE
# n_subset <- 5
# 
# if (subset) {
#   indicator_index_test_in_study <- indicator_index_test_in_study[1:n_subset,1:4]
#     for (c in 1:2) {
#       for (t in 1:n_index_tests) {
#         x_NMA[[c]][[t]] <- x_NMA[[c]][[t]][1:n_subset, ]
#       }
#     }
#     n_studies <- nrow(x_NMA[[1]][[1]])
#     n_index_tests <- length(x_NMA[[1]])
#     for (c in 1:2) {
#        x_NMA[[c]][[n_index_tests + 1]] <- NULL
#     }
# }
# x <- x_NMA


#  -----------  Compile + initialise the model using "MVP_model$new(...)"  ----------------------------------
#
# ----  Prepare the data to use: ----------------------------------------------------------------------------
#
# Read in data ("x" - aggregative cumulative counts w/ "missing thresholds" in each study marked as "-1"):
#
# x <- list(x_nd, x_d)
# saveRDS(object = x, file = file.path(local_pkg_dir, "MetaOrdDTA_example_data_1.RDS"))
## x <- readRDS(file.path(local_pkg_dir, "inst", "examples", "MetaOrdDTA_example_data_1.RDS"))

# x_full <- x
# 
# x_subset <- list()
# N_subset <- n_studies
# for (c in 1:2) {
#   x_subset[[c]] <- x[[c]][1:N_subset, ]
# }

n_studies <- nrow(x[[1]])
n_thr <- ncol(x[[1]])

##
## ----  Initialise / select model: --------------------------------------------------------------------------- NMA:
##
network <- FALSE
##
softplus <- TRUE




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
# 

{
  model_parameterisation = "R&G"
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






n_chains <-  4 ##  parallel::detectCores() / 2
# internal_obj$outs_data$stan_data_list$x
# internal_obj$outs_data$stan_data_list$n
# 
# any(is.na(unlist(internal_obj$outs_data$stan_data_list$n)))
init_lists_per_chain <- NULL

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
                  THREAD_FLAGS = "-pthread -DSTAN_THREADS"
                  )

{
  
  # model_prep_obj$internal_obj$outs_data$stan_data_list
  init_lists_per_chain <- model_prep_obj$init_lists_per_chain
  priors <- model_prep_obj$priors
  ##
  n_studies <- nrow(x[[1]])
  ##  n_index_tests <- 1# length(x[[1]])
  
  stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
    
  stan_data_list$cutpoint_index

}
 
stan_data_list$x[[2]] 
stan_data_list$n[[1]][1, ]

##
str(init_lists_per_chain)
# 
# model_prep_obj$internal_obj$outs_data$stan_data_list$box_cox = FALSE



# n_thr_random = n_thr * n_studies
# n_total_C_if_random <- sum(n_thr_random)
# 
# 
# init_lists_per_chain[[1]]$C_raw_vec <-  rep(-2.0, n_total_C_if_random)
# 
# priors$n_total_C_if_random <- n_total_C_if_random

# 
# n_thr <- model_prep_obj$internal_obj$outs_data$n_thr
# n_cat <- n_thr + 1
# n_total_pooled_cat  <- sum(n_cat)
# n_thr_random = n_thr * n_studies
# n_total_C_if_random <- sum(n_thr_random)
# ##
# for (kk in 1:n_chains) {
#   init_lists_per_chain[[kk]]$C_raw_vec <- if_null_then_set_to(init_lists_per_chain[[kk]]$C_raw_vec, rep(-2.0, n_total_C_if_random))
# }

##  init_lists_per_chain = model_samples_obj$init_lists_per_chain


# for (kk in 1:n_chains) {
#   init_lists_per_chain[[kk]]$C_raw_vec <- C_raw # rep(-2, n_thr)
#   init_lists_per_chain[[kk]]$raw_scale <- rep(0.001, n_studies)
#   init_lists_per_chain[[kk]]$beta_mu <- 0.01
# init_lists_per_chain[[kk]]$C <- NULL # seq(from = -2.0, to = 2.0, length = n_thr)
# }

init_lists_per_chain <- vector(length = n_chains, "list")
priors <- NULL


for (kk in 1:n_chains) {
  init_lists_per_chain[[kk]]$C_raw_vec <- list( rep(-3.0, n_thr), rep(-3.0, n_thr))
 # init_lists_per_chain[[kk]]$C <- list( seq(from = -2.0, to = 2.0, length = n_thr), 
  #                                      seq(from = -2.0, to = 2.0, length = n_thr))
}

if (model_parameterisation == "Xu") {
      
          for (kk in 1:n_chains) {
                
                mat <- matrix(-5.0, nrow = n_studies, ncol = n_thr)
                init_lists_per_chain[[kk]]$C_raw <-  list(mat, mat)
                init_lists_per_chain[[kk]]$alpha <-  list( rep(10.0, n_thr + 1), rep(10.0, n_thr + 1))
                ##
                init_lists_per_chain[[kk]]$beta_mu <- c(-1.0, +1.0)
                init_lists_per_chain[[kk]]$beta_z <- array(0.001, dim = c(n_studies, 2))
                init_lists_per_chain[[kk]]$beta_SD <- rep(0.001, 2)
                init_lists_per_chain[[kk]]$beta_corr <- 0.001
                ##
                # C <- seq(from = -2, to = +2, length = n_thr)
                # init_lists_per_chain[[kk]]$C <- list(C, C)
            
          }
        
          priors$prior_beta_mu_SD <- rep(1.5, 2)
          priors$prior_beta_SD_SD <- rep(1.0, 2)

} else if (model_parameterisation == "R&G") { 
  
} else if (model_parameterisation == "Jones") { 
  
        priors$prior_beta_mu_SD <- rep(1.5, 2)
        priors$prior_beta_SD_SD <- rep(5, 2)
        ##
        priors$prior_raw_scale_mu_SD <- rep(1.5, 2)
        priors$prior_raw_scale_SD_SD <- rep(5, 2)
  
}


priors



# x <- rep(1, 100)
# y <- c(rep(1, 20), rep(0, 5), rep(1, 20), rep(0, 5), rep(1, 20), rep(0, 5), rep(1, 20), rep(0, 25))
# x_latent <- y_latent <- c()
# ##
# for (i in 1:length(x)) { 
#   if (x[i] == 1) { 
#     x_latent[i] <-  LaplacesDemon:::rtrunc(spec = "norm", a = 0, b = Inf,  n = 1 , mean = 0, sd = 1)
#   } else { 
#     x_latent[i] <-  LaplacesDemon:::rtrunc(spec = "norm", a = -Inf, b = 0, n = 1 , mean = 0, sd = 1)
#   }
#   if (y[i] == 1) { 
#     y_latent[i] <-  LaplacesDemon:::rtrunc(spec = "norm", a = 0, b = Inf,  n = 1 , mean = 0, sd = 1)
#   } else { 
#     y_latent[i] <-  LaplacesDemon:::rtrunc(spec = "norm", a = -Inf, b = 0, n = 1 , mean = 0, sd = 1)
#   }
# }
# 
#  
#  
# 
# 
# cor(x_latent, y_latent)



model_parameterisation
random_thresholds


{
  
     # init_lists_per_chain = NULL
    # 
    # priors$prior_raw_scale_mu_SD <- 0.25
    # priors$prior_raw_scale_SD_SD <- 0.125
    # 
    # priors$prior_beta_mu_SD <- 1
    
    priors$prior_alpha_SD
    
    prior_alpha_mean <- 5
    prior_alpha_SD   <- 2.5
    alpha_lb <- 1
    
    priors$prior_alpha_mean <- list( rep(prior_alpha_mean, n_thr + 1),  rep(prior_alpha_mean, n_thr + 1))
    priors$prior_alpha_SD   <- list( rep(prior_alpha_SD, n_thr + 1), rep(prior_alpha_SD, n_thr + 1))
    priors$alpha_lb <- alpha_lb


}


#### init_lists_per_chain  = NULL

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

tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 100)
tibble_tp <- model_summary_and_trace_obj$get_summary_transformed() %>% print(n = 100)
tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)

tibble_main

# tibble_main %>% print(n = 1000)
# tibble_gq   %>% print(n = 1000)
# 
# #### dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "LRpos")))
# dplyr::filter(tibble_tp, (stringr::str_detect(parameter, "C"))) %>% print(n = 1000)
# dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "biv_equiv_C"))) %>% print(n = 1000)
# 
# dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "biv_equiv_C"))) %>% print(n = 1000)
# 
# stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
# Se <- stan_mod_samples$summary(c("Se"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
# 
# 
# stan_mod_samples$summary(c("C_mu"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
# 
# stan_mod_samples$summary(c("alpha"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
# 
# stan_mod_samples$summary(c("beta_SD"), mean, quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
# 
# 
# dplyr::filter(tibble_main, (stringr::str_detect(parameter, "alpha"))) %>% print(n = 1000)

## 
## ----  Plots: ----------------------------------------------------------------------
##

model_summary_and_trace_obj$plot_sROC()

# model_summary_and_trace_obj$plot_traces(param_string = "alpha")


{
        n_thr_max <- max(model_prep_obj$internal_obj$outs_data$n_thr)
        ##
        Se_median_array <- Sp_median_array <- Fp_median_array    <- array(dim = c(n_index_tests, n_thr_max))
        Se_mean_array <- Sp_mean_array <- Fp_mean_array    <- array(dim = c(n_index_tests, n_thr_max))
        ##
        Se_lower_array <- Sp_lower_array <- Fp_lower_array <- array(dim = c(n_index_tests, n_thr_max))
        Se_upper_array <- Sp_upper_array <- Fp_upper_array <- array(dim = c(n_index_tests, n_thr_max))
        ##
        Se_pred_lower_array <- Sp_pred_lower_array <- Fp_pred_lower_array <- array(dim = c(n_index_tests, n_thr_max))
        Se_pred_upper_array <- Sp_pred_upper_array <- Fp_pred_upper_array <- array(dim = c(n_index_tests, n_thr_max))
        
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
          for (t in 1:n_index_tests) {
              # test_vec[counter] <- t
              ## Posterior medians (of pooled estimates):
              Se_median_array[t, k] <- Se$`50%`[counter]
              Sp_median_array[t, k] <- Sp$`50%`[counter]
              Fp_median_array[t, k] <- Fp$`50%`[counter]
              ## Posterior means (of pooled estimates):
              Se_mean_array[t, k] <- Se$mean[counter]
              Sp_mean_array[t, k] <- Sp$mean[counter]
              Fp_mean_array[t, k] <- Fp$mean[counter]
              ## Posterior lower 95% (of pooled estimates):
              Se_lower_array[t, k] <- Se$`2.5%`[counter]
              Sp_lower_array[t, k] <- Sp$`2.5%`[counter]
              Fp_lower_array[t, k] <- Fp$`2.5%`[counter]
              ## Posterior upper 95% (of pooled estimates):
              Se_upper_array[t, k] <- Se$`97.5%`[counter]
              Sp_upper_array[t, k] <- Sp$`97.5%`[counter]
              Fp_upper_array[t, k] <- Fp$`97.5%`[counter]
              ## Posterior lower prediction 95% (of pooled estimates):
              Se_pred_lower_array[t, k] <- Se_pred$`2.5%`[counter]
              Sp_pred_lower_array[t, k] <- Sp_pred$`2.5%`[counter]
              Fp_pred_lower_array[t, k] <- Fp_pred$`2.5%`[counter]
              ## Posterior upper prediction 95% (of pooled estimates):
              Se_pred_upper_array[t, k] <- Se_pred$`97.5%`[counter]
              Sp_pred_upper_array[t, k] <- Sp_pred$`97.5%`[counter]
              Fp_pred_upper_array[t, k] <- Fp_pred$`97.5%`[counter]
              ##
              counter <- counter + 1
          }
        }
}
 



# Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
# Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
# Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
# Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()



n_thr <- 27
t <- 4

vec_index_1 <- 1:27
vec_index_2 <- vec_index_1 - 0

{
    Se_abs_diffs_sim <- abs(true_Se_OVERALL_weighted[[t]][vec_index_1] - 100 * Se_mean_array[1, 1:n_thr][vec_index_2])
    Sp_abs_diffs_sim <- abs(true_Sp_OVERALL_weighted[[t]][vec_index_1] - 100 * Sp_mean_array[1, 1:n_thr][vec_index_2])
    ##
    Se_mean_of_diffs_sim <- mean(Se_abs_diffs_sim)
    Sp_mean_of_diffs_sim <- mean(Sp_abs_diffs_sim)
    ##
    Se_sum_of_diffs_sim <- sum(Se_abs_diffs_sim)
    Sp_sum_of_diffs_sim <- sum(Sp_abs_diffs_sim)
    ##
    Se_max_of_diffs_sim <- max(Se_abs_diffs_sim)
    Sp_max_of_diffs_sim <- max(Sp_abs_diffs_sim)
}



{
      message(paste("Se_mean_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_mean_of_diffs_sim)
      ##
      message(paste("Se_max_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_max_of_diffs_sim)
      ##
      message(paste("Se_sum_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_sum_of_diffs_sim)

}
##
{
      message(paste("Sp_mean_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_mean_of_diffs_sim)
      ##
      message(paste("Sp_max_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_max_of_diffs_sim)
      ##
      message(paste("Sp_sum_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_sum_of_diffs_sim)

}




ppc_outs <- induced_Dirichlet_ppc_plot(method = "alpha", 
                           N = 5000,
                           n_cat = n_thr + 1,
                           other_args_list = list(use_log_alpha = FALSE, 
                                                  prior_alpha_mean = prior_alpha_mean, 
                                                  prior_alpha_sd = prior_alpha_SD, 
                                                  alpha_lb = alpha_lb,
                                                  alpha_ub = NULL))


## ppc_outs$log_alpha










true_Se_OVERALL_weighted[[4]] - 100 * Se_mean_array[1, 1:n_thr[4]]
  
true_Se_OVERALL_weighted[[1]] - 100 * Se_mean_array[1, 1:n_thr[1]]
true_Se_OVERALL_weighted[[2]] - 100 * Se_mean_array[2, 1:n_thr[2]]
true_Se_OVERALL_weighted[[3]] - 100 * Se_mean_array[3, 1:n_thr[3]]
true_Se_OVERALL_weighted[[4]] - 100 * Se_mean_array[4, 1:n_thr[4]]
# 
# 
# true_Sp_OVERALL_weighted[[1]] - 100 * Sp_mean_array[1, 1:n_thr[1]] 
# true_Sp_OVERALL_weighted[[2]] - 100 * Sp_mean_array[2, 1:n_thr[2]] 
# true_Sp_OVERALL_weighted[[3]] - 100 * Sp_mean_array[3, 1:n_thr[3]] 
# true_Sp_OVERALL_weighted[[4]] - 100 * Sp_mean_array[4, 1:n_thr[4]] 
# 
# 




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
# 
{
  
  
        
        basic_model_options = list(
          network = NULL,
          cts = NULL,
          prior_only = NULL
        )
        ####
        advanced_model_options = list(
          model_parameterisation = NULL,
          random_thresholds = NULL,
          Dirichlet_random_effects_type = NULL,
          box_cox = NULL,
          softplus = NULL
        )
        ####
        other_advanced_options = list(
          advanced_compile = NULL, ## default is false (aka a "basic" compile)
          ##
          force_recompile = NULL,
          quiet = NULL,
          compile = NULL,
          ##
          set_custom_CXX_CPP_flags = NULL,
          CCACHE_PATH = NULL,
          custom_cpp_user_header_file_path = NULL,
          CXX_COMPILER_PATH = NULL,
          CPP_COMPILER_PATH = NULL,
          MATH_FLAGS = NULL,
          FMA_FLAGS = NULL,
          AVX_FLAGS = NULL,
          THREAD_FLAGS = NULL
        )
        ##
        ####
        priors = NULL
        ####
        # init_lists_per_chain = NULL,
        ####
        MCMC_params = list(
          seed = NULL,
          n_superchains = NULL,
          n_chains = NULL,
          n_iter = NULL,
          n_burnin = NULL,
          adapt_delta = NULL,
          max_treedepth = NULL,
          metric_shape = NULL
        )
        
        ##
        ## ---- Main "internal_obj" list:
        ##
        internal_obj = list(
          outs_data = NULL,
          outs_stan_model_name = NULL,
          outs_stan_compile = NULL,
          outs_stan_init = NULL,
          outs_stan_sampling = NULL,
          ##
          HMC_info = NULL, ## new
          efficiency_info = NULL, ## new
          summaries = NULL, ## new
          traces = NULL ## new
        )
        
        
        
      
      debugging <- FALSE
      ##
      n_iter = 500
      n_burnin = 500
      ##
      priors <- NULL
      ##
      x = x
      n_chains = n_chains
      ##
      cts = FALSE
      network = FALSE
      ##
      prior_only = FALSE
      ##
      softplus = TRUE
      ##
      ##
      model_parameterisation = "Gatsonis"
      random_thresholds = FALSE
      Dirichlet_random_effects_type = "none"
      box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
      ##
      init_lists_per_chain = NULL
      
      
      basic_model_options$network <- network
      basic_model_options$cts <- cts
      basic_model_options$prior_only <- prior_only
      ##
      advanced_model_options$model_parameterisation <- model_parameterisation
      advanced_model_options$random_thresholds <- random_thresholds
      advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
      advanced_model_options$box_cox <- box_cox
      advanced_model_options$softplus <- softplus
      ##
      MCMC_params$n_chains <- n_chains
      
      



}





outs_data = list(
  stan_data_list = NULL,
  n_tests = NULL,
  n_studies = NULL,
  n_thr = NULL,
  n_cat = NULL
)
##
outs_stan_model_name = list(
  stan_model_file_name = NULL
)
##
outs_stan_compile = list(
  stan_model_obj = NULL,
  stan_model_file_name = NULL,
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
)
##
outs_stan_init = list(
  inits_unconstrained_vec_per_chain = NULL,
  stan_param_names_list = NULL,
  stan_param_names_main = NULL,
  stan_init_pseudo_sampling_outs = NULL,
  stan_model_obj = NULL,
  json_file_path = NULL,
  stan_model_file_path = NULL
)
##
outs_stan_sampling = list(
  stan_mod_samples = NULL,
  time_total = NULL
)











