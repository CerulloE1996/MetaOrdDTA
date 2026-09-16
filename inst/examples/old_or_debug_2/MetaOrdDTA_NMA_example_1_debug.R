

 
 
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




##
## ---- Load NMA data:
##
setwd(local_pkg_dir)
data <- readRDS("data_example_1_NMA_list.RDS")
x_NMA <- data$x_NMA
x <- x_NMA
indicator_index_test_in_study <- data$indicator_index_test_in_study
##
n_studies <- nrow(x_NMA[[1]][[1]])
n_index_tests <- 4 ##  length(x_NMA[[1]])

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


##  -----------  Compile + initialise the model using "MVP_model$new(...)"  ----------------------------------
##
## ----  Prepare the data to use: ----------------------------------------------------------------------------
##
## Read in data ("x" - aggregative cumulative counts w/ "missing thresholds" in each study marked as "-1"):
##
# x <- list(x_nd, x_d)
# saveRDS(object = x, file = file.path(local_pkg_dir, "MetaOrdDTA_example_data_1.RDS"))
# x <- readRDS(file.path(local_pkg_dir, "inst", "examples", "MetaOrdDTA_example_data_1.RDS"))
# 
# x <- list(x_nd, x_d)
# 
# x_full <- x
# 
# x_subset <- list()
# N_subset <- n_studies
# for (c in 1:2) {
#   x_subset[[c]] <- x[[c]][1:N_subset, ]
# }
#  
# n_studies <- nrow(x_subset[[1]])
# n_thr <- ncol(x_subset[[1]])

##
## ----  Initialise / select model: --------------------------------------------------------------------------- NMA:
##
network <- TRUE
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
  model_parameterisation = "Gatsonis"
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





{
  model_parameterisation = "Gatsonis"
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





n_chains <- 8 #### parallel::detectCores() / 2
# internal_obj$outs_data$stan_data_list$x
# internal_obj$outs_data$stan_data_list$n
# 
# any(is.na(unlist(internal_obj$outs_data$stan_data_list$n)))
init_lists_per_chain <- NULL

model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  debugging = TRUE,
                                                  ##
                                                  x = x_NMA, 
                                                  indicator_index_test_in_study = indicator_index_test_in_study,
                                                  ## n_index_tests_per_study = n_index_tests_per_study,
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
                                                  advanced_compile = TRUE, ## default is standard/"basic" compile
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

# model_prep_obj$internal_obj$outs_data$stan_data_list
init_lists_per_chain <- model_prep_obj$init_lists_per_chain
priors <- model_prep_obj$priors
##
n_studies <- nrow(x[[1]][[1]])
n_index_tests <- length(x[[1]])

init_lists_per_chain

# if (model_parameterisation == "Jones") {
#     priors$prior_beta_mu_SD      <- array(5.0, dim = c(n_index_tests, 2))
#     priors$prior_raw_scale_mu_SD <- array(5.0, dim = c(n_index_tests, 2))
#     ##
#     priors$prior_beta_sigma_SD      <- rep(1.0, 2)
#     priors$prior_raw_scale_sigma_SD <- rep(1.0, 2)
#     ##
#     priors$prior_beta_tau_SD      <- array(1.0, dim = c(n_index_tests, 2))
#     priors$prior_raw_scale_tau_SD <- array(1.0, dim = c(n_index_tests, 2))
#     ##
#     ##
#     ##
#     ##
#     for (kk in 1:n_chains) {
#       init_lists_per_chain[[kk]]$beta_eta <- array(0.001, dim = c(n_studies, 2))
#       init_lists_per_chain[[kk]]$raw_scale_eta <- array(0.001, dim = c(n_studies, 2))
#       ##
#       init_lists_per_chain[[kk]]$beta_delta <- list()
#       init_lists_per_chain[[kk]]$raw_scale_delta <- list()
#       for (t in 1:n_index_tests) {
#         init_lists_per_chain[[kk]]$beta_delta[[t]] <- array(0.001, dim = c(n_studies, 2))
#         init_lists_per_chain[[kk]]$raw_scale_delta[[t]] <- array(0.001, dim = c(n_studies, 2))
#       }
#     }
#     
# }


# 
# if (model_parameterisation == "Xu") {
#   
#     # n_total_C_if_fixed <- 0
#     # for (t in 1:n_index_tests) { 
#     #   # for (k in 1:n_thr[t]) {
#     #      n_total_C_if_fixed <- n_total_C_if_fixed +  ncol(x[[1]][[t]])
#     #   # }
#     # }
#     # priors$n_total_C_if_fixed <- n_total_C_if_fixed  
#     # ##
#     # priors$prior_beta_mu_mean     <- array(0.0, dim = c(n_index_tests, 2))
#     # ##
#     # priors$prior_beta_mu_SD       <- array(1.0, dim = c(n_index_tests, 2))
#     # ##
#     # priors$prior_beta_sigma_SD    <- rep(1.0, 2)
#     # ##
#     # priors$prior_beta_tau_SD      <- array(1.0, dim = c(n_index_tests, 2))
#     # ##
#     # priors$prior_dirichlet_alpha <- array(1.0, dim = c(n_index_tests, 28))
#     ##
#     ##
#     ##
#     for (kk in 1:n_chains) {
#       init_lists_per_chain[[kk]]$C_raw_vec <- rep(-2.0, n_total_C_if_fixed)
#       ##
#       init_lists_per_chain[[kk]]$beta_mu <- array(c(rep(-1.0, n_index_tests), rep(+1.0, n_index_tests)),
#                                                   dim = c(n_index_tests, 2))
#       ##
#       init_lists_per_chain[[kk]]$beta_sigma <- rep(0.001, 2)
#       init_lists_per_chain[[kk]]$beta_tau <- array(0.001, dim = c(n_index_tests, 2))
#       ##
#       init_lists_per_chain[[kk]]$beta_eta <- array(0.001, dim = c(n_studies, 2))
#       ##
#       init_lists_per_chain[[kk]]$beta_delta <- list()
#       for (t in 1:n_index_tests) {
#         init_lists_per_chain[[kk]]$beta_delta[[t]] <- array(0.001, dim = c(n_studies, 2))
#       }
#     }
#   
# }


##
str(init_lists_per_chain)
# 
# model_prep_obj$internal_obj$outs_data$stan_data_list$box_cox = FALSE

priors


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

tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 100)
##
tibble_gq <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 100)


tibble_main %>% print(n = 1000)
tibble_gq   %>% print(n = 1000)

dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "LRpos")))


 

##
## ----  Plots: ----------------------------------------------------------------------
##

sroc_plots_NMA <- model_summary_and_trace_obj$plot_sROC()

sroc_plots_NMA

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

# test_vec <- thr_vec <- c()
# counter <- 1 
# for (t in 1:n_index_tests) {
#   for (k in 1:n_thr_max) {
#       test_vec[counter] <- t
#       thr_vec[counter] <- k
#       ##
#       counter <- counter + 1
#   }
# }
# 
# n_rows_total <- length(test_vec)
# 
# ##
# ## Make tibble for ggplot (to make sROC plots, etc):
# ##
# Se_mean_array
# c(t(Se_mean_array))
# 
# 
# stan_model_file_name <- "NMA_Xu"
# model_name <- stan_model_file_name
# 
# tibble_NMA <- tibble(  Model = rep(model_name, n_rows_total),
#                        ##
#                        test = factor(test_vec),
#                        test_char = paste0("test ", test),
#                        threshold = thr_vec,
#                        ##
#                        Se_median = c(t(Se_median_array)),
#                        Sp_median = c(t(Sp_median_array)),
#                        Fp_median = c(t(Fp_median_array)),
#                        ##
#                        Se_mean = c(t(Se_mean_array)),
#                        Sp_mean = c(t(Sp_mean_array)),
#                        Fp_mean = c(t(Fp_mean_array)),
#                        ##
#                        Se_lower = c(t(Se_lower_array)),
#                        Sp_lower = c(t(Sp_lower_array)),
#                        Fp_lower = c(t(Fp_lower_array)),
#                        ##
#                        Se_upper = c(t(Se_upper_array)),
#                        Sp_upper = c(t(Sp_upper_array)),
#                        Fp_upper = c(t(Fp_upper_array)),
#                        ##
#                        Se_pred_lower = c(t(Se_pred_lower_array)),
#                        Sp_pred_lower = c(t(Sp_pred_lower_array)),
#                        Fp_pred_lower = c(t(Fp_pred_lower_array)),
#                        ##
#                        Se_pred_upper = c(t(Se_pred_upper_array)),
#                        Sp_pred_upper = c(t(Sp_pred_upper_array)),
#                        Fp_pred_upper = c(t(Fp_pred_upper_array)))
# 
# 
# 
# 
# tibble_NMA
# polygon_Conf_list <- polygon_Pred_list <- list()
# for (t in 1:n_index_tests) {
#     tibble_test_t <- dplyr::filter(tibble_NMA, test == t)
#     ##
#     ## Tibbles for 95% credible region:
#     ##
#     tibble_Conf_polygon_test_t <- create_confidence_polygon(tibble_test_t, model_name = model_name)
#     n_rows = nrow(tibble_Conf_polygon_test_t)
#     tibble_Conf_polygon_test_t <- tibble_Conf_polygon_test_t %>% dplyr::mutate(test = rep(t, n_rows),
#                                                                                test_char = paste0("test ", test))
#     polygon_Conf_list[[t]] <- tibble_Conf_polygon_test_t
#     ##
#     ## Tibbles for 95% prediction region:
#     ##
#     tibble_Pred_polygon_test_t <- create_prediction_polygon(tibble_test_t, model_name = model_name)
#     n_rows = nrow(tibble_Pred_polygon_test_t)
#     tibble_Pred_polygon_test_t <- tibble_Pred_polygon_test_t %>% dplyr::mutate(test = rep(t, n_rows), 
#                                                                                test_char = paste0("test ", test))
#     polygon_Pred_list[[t]] <- tibble_Pred_polygon_test_t
#     # polygon_Pred_list[[t]] <- create_prediction_polygon(tibble_test_t, model_name = model_name)
# }
# ##
# ## ---- Final conf and prediction regions tibbles:
# ##
# polygon_Conf_tibble <- tibble(data.table::rbindlist(polygon_Conf_list)) %>% dplyr::filter(!(is.na(x))) %>% print(n = 100)
# polygon_Pred_tibble <- tibble(data.table::rbindlist(polygon_Pred_list)) %>% dplyr::filter(!(is.na(x))) %>% print(n = 100)
# 
# ##
# ## ---- Plot 1 (all tests on 1 plot):
# ##
# plot_1 <-  ggplot(tibble_NMA, mapping = aes(x = Fp_median, y = Se_median, colour = test)) + 
#   geom_line(linewidth = 0.5) + 
#   geom_point(size = 3) + 
#   theme_bw(base_size = 16) + 
#   xlab("False positive rate (Fp)") + 
#   ylab("Sensitivity (Se)")
# plot_1
# 
# ##
# ## ---- Plot 2 (each test on separate panel):
# ##
# plot_2 <-  ggplot(tibble_NMA, mapping = aes(x = Fp_median, y = Se_median, colour = Model)) + 
#   geom_line(linewidth = 0.5) + 
#   geom_point(size = 3) + 
#   theme_bw(base_size = 16) + 
#   facet_wrap(~ test_char) + 
#   xlab("False positive rate (Fp)") + 
#   ylab("Sensitivity (Se)")
# plot_2
# 
# 
# ##
# ## ---- Plot 3 (plot w/ 95% confidence region):
# ##
# conf_region_colour <- pred_region_colour <- "blue"
# ##
# plot_3 <-   plot_2  + 
#   ##
#   geom_polygon(data = polygon_Conf_tibble,
#                aes(x = x, y = y, colour = Model),
#                fill = conf_region_colour, 
#                alpha = 0.40)
#   # geom_polygon(data = polygon_Pred_tibble, aes(x = x, y = y), fill = pred_region_colour, alpha = 0.25) + 
#   ##
# plot_3
# 
# 
# 
# ##
# ## ---- Plot 4 (plot w/ 95% prediction region):
# ##
# plot_4 <-  plot_2 +  geom_polygon(
#                data = polygon_Pred_tibble,
#                aes(x = x, y = y, colour = Model),
#                fill = pred_region_colour, 
#                alpha = 0.40)
# plot_4
# 
# 
# ##
# ## ---- Plot 5 (plot w/ BOTH the 95% confidence + prediction regions):
# ##
# plot_5 <-  plot_2 + 
#   geom_polygon(
#       data = polygon_Conf_tibble,
#       aes(x = x, y = y, colour = Model),
#       fill = pred_region_colour, 
#       alpha = 0.40) + 
#   ##
#   geom_polygon(
#       data = polygon_Pred_tibble,
#       aes(x = x, y = y, colour = Model),
#       fill = pred_region_colour, 
#       alpha = 0.20)
# plot_5
# 
# 
# 
# 
# 

 

#####
#####
#####


Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()

for (t in 1:n_index_tests) {
    Se_abs_diffs_sim[[t]] <- abs(true_Se_OVERALL_weighted[[t]] - 100 * Se_mean_array[t, 1:n_thr[t]])
    Sp_abs_diffs_sim[[t]] <- abs(true_Sp_OVERALL_weighted[[t]] - 100 * Sp_mean_array[t, 1:n_thr[t]])
    ##
    Se_mean_of_diffs_sim[t] <- mean(Se_abs_diffs_sim[[t]])
    Sp_mean_of_diffs_sim[t] <- mean(Sp_abs_diffs_sim[[t]])
    ##
    Se_sum_of_diffs_sim[t] <- sum(Se_abs_diffs_sim[[t]])
    Sp_sum_of_diffs_sim[t] <- sum(Sp_abs_diffs_sim[[t]])
    ##
    Se_max_of_diffs_sim[t] <- max(Se_abs_diffs_sim[[t]])
    Sp_max_of_diffs_sim[t] <- max(Sp_abs_diffs_sim[[t]])
}



{
      message(paste("Se_mean_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_mean_of_diffs_sim)
      ##
      message(paste("Se_sum_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_sum_of_diffs_sim)
      ##
      message(paste("Se_max_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Se_max_of_diffs_sim)
}
##
{
      message(paste("Sp_mean_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_mean_of_diffs_sim)
      ##
      message(paste("Sp_sum_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_sum_of_diffs_sim)
      ##
      message(paste("Sp_max_of_diffs_sim = ")) ## , Se_sum_of_diffs_sim))
      print(Sp_max_of_diffs_sim)
}



true_Se_OVERALL_weighted[[1]] - 100 * Se_mean_array[1, 1:n_thr[1]] 
true_Se_OVERALL_weighted[[2]] - 100 * Se_mean_array[2, 1:n_thr[2]] 
true_Se_OVERALL_weighted[[3]] - 100 * Se_mean_array[3, 1:n_thr[3]] 
true_Se_OVERALL_weighted[[4]] - 100 * Se_mean_array[4, 1:n_thr[4]] 


true_Sp_OVERALL_weighted[[1]] - 100 * Sp_mean_array[1, 1:n_thr[1]] 
true_Sp_OVERALL_weighted[[2]] - 100 * Sp_mean_array[2, 1:n_thr[2]] 
true_Sp_OVERALL_weighted[[3]] - 100 * Sp_mean_array[3, 1:n_thr[3]] 
true_Sp_OVERALL_weighted[[4]] - 100 * Sp_mean_array[4, 1:n_thr[4]] 






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
      cts = TRUE
      network = TRUE
      ##
      prior_only = FALSE
      ##
      softplus = TRUE
      ##
      ##
      model_parameterisation = "Jones"
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

X <- NULL

x <- x_NMA

MetaOrdDTA:::prep_data_and_model(  
  debugging = self$debugging,
  ##
  x = self$x,
  ##
  X = self$X, ## covariates
  ##
  indicator_index_test_in_study = self$indicator_index_test_in_study,
  ##
  internal_obj = self$internal_obj,
  ##
  basic_model_options    = self$basic_model_options,
  advanced_model_options = self$advanced_model_options,
  MCMC_params            = self$MCMC_params,
  other_advanced_options = self$other_advanced_options,
  ##
  priors = self$priors,
  init_lists_per_chain = self$init_lists_per_chain,
  ##
  compute_sim_study_metrics = self$compute_sim_study_metrics,
  vec_index_inner_thr       = self$vec_index_inner_thr)










