


 

## Set wd:
setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")
setwd("/home/enzo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")




require(BayesMVP)
require(TruncatedNormal)
 

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")
source("missing_thr_prep_Stan_data.R")
source("missing_thr_prior_pred_check.R")
source("R_fn_compile_Stan_model.R")


options(scipen = 999999999999)



## Model args:
model_args_list <- list()


#### To overwrite some model options:
{
    
    ## Model options for JONES model:
    ##
    Stan_data_list$use_box_cox <- 1
    ##  
    ## Model options for "Cerullo" model:
    ##
    Stan_data_list$estimate_scales <- 0
    #### Stan_data_list$same_cutpoints_between_groups <- 1
    ##
    Stan_data_list$prior_only <- 0
    ##
    Stan_data_list$use_empirical_cutpoint_means <- 0
    ##
    
      method <- "sigma"
      # method <- "kappa"
    
    if (method == "sigma") { 
      
            Stan_data_list$prior_dirichlet_cat_means_alpha <- rep(1.00, n_cat)
            ##
            Stan_data_list$prior_dirichlet_cat_SDs_mean <- rep(0.00, n_cat)
            # prior_dirichlet_cat_SDs_SD   <- rep(0.50, n_cat) ## example of when "vague" prior is very NOT vague and bad !!!
            Stan_data_list$prior_dirichlet_cat_SDs_SD   <- rep(0.025, n_cat) ## example of when "vague" prior is very NOT vague and bad !!!
            # ##
            # Stan_data_list$prior_dirichlet_cat_means_alpha <- prior_dirichlet_cat_means_alpha
            # Stan_data_list$prior_dirichlet_cat_SDs_mean <- prior_dirichlet_cat_SDs_mean
            # Stan_data_list$prior_dirichlet_cat_SDs_SD <- prior_dirichlet_cat_SDs_SD
            # ##
            other_args_list <- list( prior_dirichlet_cat_means_alpha = prior_dirichlet_cat_means_alpha,
                                     prior_dirichlet_cat_SDs_mean = prior_dirichlet_cat_SDs_mean,
                                     prior_dirichlet_cat_SDs_SD = prior_dirichlet_cat_SDs_SD)
            ##
            inf_dir_samples <- induced_Dirichlet_ppc_plot(   method = "sigma",
                                                             N = 5000,
                                                             n_cat = n_cat,
                                                             other_args_list = other_args_list)
      
      
    } else if (method == "kappa") {
  
      
            prior_kappa_mean <- 0
            prior_kappa_SD <-   200
            prior_dirichlet_cat_means_alpha <- rep(1, n_cat)
            ##
            Stan_data_list$prior_kappa_mean <-  rep(prior_kappa_mean, 1)
            Stan_data_list$prior_kappa_SD <-    rep(prior_kappa_SD, 1)
            Stan_data_list$prior_dirichlet_cat_means_alpha <- prior_dirichlet_cat_means_alpha
            ##
            use_log_kappa <- TRUE
            kappa_lb <- 1
            Stan_data_list$log_alpha_lb <- log(kappa_lb)
            Stan_data_list$log_alpha_ub <- +Inf
            print(paste("kappa_lb = ")) ; print(exp( Stan_data_list$log_kappa_lb ))
            print(paste("kappa_ub = ")) ; print(exp( Stan_data_list$log_kappa_ub ))
            ##
              other_args_list <- list(  use_log_kappa = use_log_kappa,
                                        log_kappa_lb = log_kappa_lb,
                                        log_kappa_ub = log_kappa_ub,
                                        prior_log_kappa_mean = prior_log_kappa_mean,
                                        prior_log_kappa_sd   = prior_log_kappa_sd)
            ##
            inf_dir_samples <- induced_Dirichlet_ppc_plot(   method = "kappa",
                                                             N = 5000,
                                                             n_cat = n_cat,
                                                             use_log_alpha = FALSE,
                                                             use_log_kappa = FALSE,
                                                             other_args_list = other_args_list)
    
    } 
    
    kappa <- inf_dir_samples$kappa
    alpha <- inf_dir_samples$alpha
 
    print(mean(inf_dir_samples$alpha))
 
    # dev.off()
    # plot(density(inf_dir_samples$alpha), xlim = c(0, 100))
    print(mean(inf_dir_samples$alpha))
    print(quantile(inf_dir_samples$alpha, c(0.025, 0.50, 0.975)))

    # # plot(density(exp(inf_dir_samples$log_alpha)))
    # print(mean(exp(inf_dir_samples$log_alpha)))
    # print(quantile(exp(inf_dir_samples$log_alpha), c(0.025, 0.50, 0.975)))
    
}






{
  
  prior_dirichlet_cat_SDs_mean <- rep(0.00, n_cat)
  # prior_dirichlet_cat_SDs_SD   <- rep(0.50, n_cat) ## example of when "vague" prior is very NOT vague and bad !!!
  prior_dirichlet_cat_SDs_SD   <- rep(0.01, n_cat) ## example of when "vague" prior is very NOT vague and bad !!!
  prior_dirichlet_cat_means_alpha <- rep(1.00, n_cat)
  ##
  Stan_data_list$prior_dirichlet_cat_means_alpha <- prior_dirichlet_cat_means_alpha
  Stan_data_list$prior_dirichlet_cat_SDs_mean <- prior_dirichlet_cat_SDs_mean
  Stan_data_list$prior_dirichlet_cat_SDs_SD <- prior_dirichlet_cat_SDs_SD
  ##
  Stan_data_list$kappa_lb <- 1
  
}

## Stan_data_list$use_box_cox <- 0
Stan_data_list$use_box_cox <- 1




# 
# 
# ##  | ------   Initial values  -------------------------------------------------------------------------
# {
#       n_sets_of_cutpoints <- 1
#   
#       Stan_init_list <- list()
#       ##
#       ##
#       Stan_init_list$beta_mu <- rep(0.00001, 2)
#       Stan_init_list$beta_SD <- rep(0.01, 2)
#       Stan_init_list$beta_z <-  array(0.0, dim = c(2, n_studies))
#       ##
#       # Stan_init_list$kappa <- rep(10, n_sets_of_cutpoints)
#       # Stan_init_list$phi <-  array(dim = c(n_sets_of_cutpoints, n_cat), data = rep(1/n_cat, n_cat))
#       ##
#       Stan_init_list$log_scale_mu <- rep(0.00001, 2)
#       Stan_init_list$log_scale_SD <- rep(0.01, 2)
#       Stan_init_list$log_scale_z <-  array(0.00001, dim = c(2, n_studies))
#  
#  
#   
# }
# 
# 





##
# model_args_list$Dirichlet_random_effects_type <- "kappa"
model_args_list$Dirichlet_random_effects_type <- "SD"

ord_model_parameterisation

random_thresholds


# install.packages("sonify")
# 
# require(sonify)
#  
# sonify::sonify(y = c(1,2,3,4))
# 
# ?sonify::sonify
# 
# 
# # Create a simple function to play sounds using aplay (Linux)
# play_sound <- function(frequency = 100, 
#                        duration = 3) {
#   system(sprintf("play -n synth %s sine %s", duration, frequency))
# }
# 
# # Usage
# play_sound(440, 1)  # Play A4 for 1 second
# play_sound(523, 1)  # Play C5 for 1 second
# 
# 
# play_sound(500, 2.5)  # Play C5 for 1 second
# 




##
## Jones-based models ("Jones"):
##
Model <- "DTA_MA_Jones.stan"
##
## Xu-based models "Cerullo-Xu":
##
Model <- "DTA_MA_Xu_FIXEDthr.stan"
Model <- "DTA_MA_Xu_RANDthr_kappa.stan"
Model <- "DTA_MA_Xu_RANDthr_SD.stan"
##
## HSROC/Gatsonis-based models "Cerullo-Gatsonis":
##
Model <- "DTA_MA_Gat_FIXEDthr.stan"
Model <- "DTA_MA_Gat_RANDthr_SD.stan"


##
## - | ----------  Set "user_input_model_args_list":  --------------------------------------------------------------------
##
{
      user_input_model_args_list <- list()
       
      # file <- "/home/enzo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/inst/stan_models/LC_MVP_bin_PartialLog_v5.stan"
      # mod <- cmdstanr::cmdstan_model(file)
       
                  if (Model == "DTA_MA_Jones.stan") {
      
                             user_input_model_args_list$ord_model_parameterisation    <- "Jones"
                             user_input_model_args_list$random_thresholds             <- FALSE
                             user_input_model_args_list$cutpoint_dirichlet_param      <- "none"
                          
                  } else if (Model == "DTA_MA_Xu_FIXEDthr.stan") {
      
                             user_input_model_args_list$ord_model_parameterisation <- "Xu"
                             user_input_model_args_list$random_thresholds          <- FALSE
                             user_input_model_args_list$cutpoint_dirichlet_param   <- "fixed"
      
                  } else if (Model == "DTA_MA_Xu_RANDthr_kappa.stan") {
                   
                                  
                             user_input_model_args_list$ord_model_parameterisation <- "Xu"
                             user_input_model_args_list$random_thresholds          <- TRUE
                             user_input_model_args_list$cutpoint_dirichlet_param   <- "kappa"
                               
                    
                  } else if (Model == "DTA_MA_Xu_RANDthr_SD.stan") {
      
                             user_input_model_args_list$ord_model_parameterisation <- "Xu"
                             user_input_model_args_list$random_thresholds          <- TRUE
                             user_input_model_args_list$cutpoint_dirichlet_param   <- "SD"
      
                   } else if (Model == "DTA_MA_Gat_FIXEDthr.stan") {
          
                             user_input_model_args_list$ord_model_parameterisation <- "Gatsonis"
                             user_input_model_args_list$random_thresholds          <- FALSE
                             user_input_model_args_list$cutpoint_dirichlet_param   <- "fixed"
                    
                   } else if (Model == "DTA_MA_Gat_RANDthr_SD.stan") {
      
                             user_input_model_args_list$ord_model_parameterisation <- "Gatsonis"
                             user_input_model_args_list$random_thresholds          <- TRUE
                             user_input_model_args_list$cutpoint_dirichlet_param   <- "SD"
      
                   }
      
          {
              print(  user_input_model_args_list$ord_model_parameterisation )
              print(  user_input_model_args_list$random_thresholds )
              print(  user_input_model_args_list$cutpoint_dirichlet_param )
          }
}

##
## Test R package for "meta_cts" (i.e., "Jones" model):
##
## Data:
##
x <- list(x_nd = x_nd, x_d = x_d)
##
## Model_type <- "meta_cts"
Model_type <- "meta_ord"
##
## Reset model_args_list:
##
model_args_list <- user_input_model_args_list



self <- list()


init_lists_per_chain <- NULL
##
self$result  <-         MetaOrdDTA:::initialise_and_update_and_run_model(   Model_type = Model_type,
                                                            ##
                                                            # initial_object = self$initial_object,
                                                            ##
                                                            x = x,
                                                            model_args_list = model_args_list,
                                                            ##
                                                            n_chains = n_chains,
                                                            ##
                                                            init_lists_per_chain = NULL, ##  init_lists_per_chain,
                                                            # inits = inits,
                                                            ##
                                                            priors = priors,
                                                            prior_only = prior_only,
                                                            ##
                                                            seed = seed,
                                                            n_iter = n_iter,
                                                            n_burnin = n_burnin,
                                                            n_superchains = n_superchains,
                                                            adapt_delta = adapt_delta,
                                                            max_treedepth = max_treedepth,
                                                            metric_shape = metric_shape)
##
full_outputs <- self$result$full_outputs
##
outs_init_stan_model     <- full_outputs$outs_init_stan_model
Stan_model_sample_output <- full_outputs$Stan_model_sample_output
##
init_lists_per_chain  <- self$result$init_lists_per_chain
Stan_data_list        <- self$result$Stan_data_list
stan_param_names_list <- self$result$stan_param_names_list
Stan_model_file_path  <- self$result$Stan_model_file_path ## BOOKMARK



 
# 
  draws_array <- self$summary_object$traces$traces_as_arrays$draws_array 
# 
# if (is.null(draws_array)) {
#   stop("draws_array is NULL in summary_object")
# }
# 
# if (is.null(params)) { # plot all params
#   
#   bayesplot::mcmc_trace(draws_array)
#   
# } else {   # plot specific params using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot 
#   
#   MetaOrdDTA:::plot_multiple_params_batched( draws_array = self$summary_object$traces$traces_as_arrays$draws_array , 
#                                              param_prefixes = params, 
#                                              plot_type = "trace", 
#                                              batch_size = batch_size)
  
##
## -----------------  Make model summary and traces + other outputs: NEED TO FIX ----------------------------------------------------
##
  
  package = "MetaOrdDTA"
  Stan_model_sample_output = Stan_model_sample_output
  stan_param_names_list = stan_param_names_list
  save_log_lik_trace = TRUE
  ##
  use_bridgestan = NULL
  ##
  model_results = NULL
  init_object = NULL
  ##
  n_nuisance = NULL
  ##
  compute_main_params = TRUE
  compute_transformed_parameters = TRUE
  compute_generated_quantities = TRUE
  save_log_lik_trace = FALSE
  compute_nested_rhat = FALSE
  n_superchains = NULL
  save_trace_tibbles = FALSE
  
  
# rm(summary_and_traces_outs)
self$model_fit_object <- create_summary_and_traces(  package = "MetaOrdDTA", 
                                                       Stan_model_sample_output = Stan_model_sample_output,
                                                       stan_param_names_list = stan_param_names_list, 
                                                       save_log_lik_trace = TRUE)

##
str(self$model_fit_object)
{
    summaries <- self$model_fit_object$summaries
    traces    <- self$model_fit_object$traces
    ##
    summary_tibbles            <- summaries$summary_tibbles
    summary_tibble_main_params <- summary_tibbles$summary_tibble_main_params
    summary_tibble_transformed_parameters <- summary_tibbles$summary_tibble_transformed_parameters
    summary_tibble_generated_quantities <- summary_tibbles$summary_tibble_generated_quantities
    ##
    n_max_trees_per_chain <- summaries$n_max_trees_per_chain
    ebfmi_per_chain       <- summaries$ebfmi_per_chain
    efficiency_info       <- summaries$efficiency_info
    HMC_info              <- summaries$HMC_info
    ##
    traces_as_arrays  <- traces$traces_as_arrays
    traces_as_tibbles <- traces$traces_as_tibbles
    log_lik_trace     <- traces$log_lik_trace
}
##
## ----------------- MCMC plots (tracee + density plots) --------------------------------------------------------------
##
## param_names <- find_matching_param_names( draws_array = model_outputs$draws_array, param_prefix = "beta")
##
draws_array <- self$model_fit_object$traces$traces_as_arrays$draws_array 
##
plot_param_group_batched(draws_array = draws_array, 
                         param_string = "beta",
                         condition = "exact_match", 
                         plot_type = "trace",
                         batch_size = 12,
                         print = FALSE)
##
plot_multiple_params_batched(draws_array = draws_array, 
                         param_string = c("beta_mu", "beta_SD"),
                         condition = "exact_match", 
                         plot_type = "trace",
                         batch_size = 12,
                         print = FALSE)

## ----------------- Test accuracy-specific plots (e.g. sROC curve / confidence + prediction regions) ----------------
##
##
## sROC plot of results:
##
sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( Stan_model_file_name = Stan_model_file_name,
                                  Stan_mod_samples = Stan_mod_samples,
                                  df_fitted = df_fitted,
                                  df_true = NULL,
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
## sROC plot of results (+ true values over-layed):
##
df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, 
                  Sp_true = true_Sp_OVERALL_weighted/100, 
                  Fp_true = (100 - true_Sp_OVERALL_weighted)/100)


df_true
###
sROC_plot_outs <- MetaOrdDTA:::R_fn_sROC_plot( Stan_model_file_name = Stan_model_file_name,
                                  Stan_mod_samples = Stan_mod_samples,
                                  df_fitted = df_fitted,
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





{
        ## Predicted data:
        #### dev    <- Stan_mod_sample$summary(c("dev"))  ####    %>% print(n = 25)
        dev_nd <- Stan_mod_sample$summary(c("dev_nd"))  #### %>% print(n = 25)
        dev_d  <- Stan_mod_sample$summary(c("dev_d"))   #### %>% print(n = 25)
        dev_nd_mat <- round( array(-999999, dim = c(n_studies, n_thr)), 3)
        dev_d_mat  <- round( array(-999999, dim = c(n_studies, n_thr)), 3)
        
        
        dev_nd_medians <- ifelse(dev_nd$median == -1, 0, dev_nd$median)
        dev_nd_means   <- ifelse(dev_nd$mean == -1, 0, dev_nd$mean)
        dev_d_medians <- ifelse(dev_d$median == -1, 0, dev_d$median)
        dev_d_means   <- ifelse(dev_d$mean == -1, 0, dev_d$mean)
        
        # sum(dev_nd_medians, na.rm = TRUE)
        # sum(dev_nd_means,   na.rm = TRUE)
        # sum(dev_d_medians, na.rm = TRUE)
        # sum(dev_d_means,   na.rm = TRUE)
  
        x_hat_nd <- Stan_mod_sample$summary(c("x_hat_nd")) ##  %>% print(n = 10)
        x_hat_d  <- Stan_mod_sample$summary(c("x_hat_d")) ##   %>% print(n = 10)
        x_hat_nd_mat <- round( array(-1, dim = c(n_studies, n_thr)), 3)
        x_hat_d_mat  <- round( array(-1, dim = c(n_studies, n_thr)), 3)
        
        x_nd <-  Stan_data_list$x[[1]]
        x_d  <-  Stan_data_list$x[[2]]
        
        x_nd <- ifelse(x_nd == -1, 0, x_nd)
        x_d  <- ifelse(x_d ==  -1, 0, x_d)
        
        counter <- 0 
        for (k in 1:n_thr) {
          for (s in 1:n_studies) {
                  counter <- counter + 1
                  x_hat_nd_mat[s, k] <-  (x_hat_nd$median[counter])
                  x_hat_d_mat[s, k]  <-  (x_hat_d$median[counter])
                  dev_nd_mat[s, k]   <-  (dev_nd$median[counter])
                  dev_d_mat[s, k]    <-  (dev_d$median[counter])
          }
        }
        
        dev_nd_per_study <- rowSums(dev_nd_mat, na.rm = TRUE)
        dev_d_per_study  <- rowSums(dev_d_mat, na.rm = TRUE)
        
        x_hat_nd_mat <- ifelse(x_hat_nd_mat == -1, 0, x_hat_nd_mat)
        x_hat_d_mat  <- ifelse(x_hat_d_mat == -1, 0, x_hat_d_mat)
        
        abs_diff_mtx_nd <-  abs(x_hat_nd_mat - x_nd)
        abs_diff_mtx_d  <-  abs(x_hat_d_mat - x_d)
        
        ## Overall differences in total cell counts:
        message(paste("Overall % difference in total cell counts (D-) = "))
        print(round(100 * sum( abs_diff_mtx_nd ) / sum( x_hat_nd_mat[x_hat_nd_mat != -1] ), 2))
        message(paste("Overall % difference in total cell counts (D+) = "))
        print(round(100 * sum( abs_diff_mtx_d )  / sum( x_hat_d_mat[x_hat_d_mat != -1] ), 2))
        
        rowSums(abs_diff_mtx_nd)
        rowSums(x_hat_d_mat)
        
        ## Study-specific differences in total cell counts:
        ## message(paste("Study-specific % difference in total cell counts (D-) = "))
        ##  print( 100 * abs_diff_mtx_nd / x_hat_nd_mat )
        ## message(paste("Study-specific % difference in total cell counts (D+) = "))
        ## print( 100 * abs_diff_mtx_d  / x_hat_d_mat )
        
        ## Deviance:
        message(paste("Overall Deviance in non-diseased group (D-) = "))
        print( sum(dev_nd_medians, na.rm = TRUE) )
        message(paste("Overall Deviance in non-diseased group (D+) = "))
        print( sum(dev_d_medians,  na.rm = TRUE) )
        
        try({  
        {
            message(paste("Se (diffs) = "))
            ##
            est_Se_OVERALL <- Stan_mod_sample$summary(c("Se")) 
            est_Se_OVERALL_mean <- est_Se_OVERALL$mean*100
            est_Se_OVERALL_median <- est_Se_OVERALL$median*100
            ##
            print( abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3)) )
            message(paste("Se (SUM of abs. diffs) = "))
            print( sum(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MEAN of abs. diffs) = "))
            print( mean(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MEDIAN of abs. diffs) = "))
            print( median(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MAX of abs. diffs) = "))
            print( max(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            
            message(paste("Sp (diffs) = "))
            ##
            est_Sp_OVERALL <- Stan_mod_sample$summary(c("Sp")) ; round(est_Sp_OVERALL$mean, 3)
            est_Sp_OVERALL_mean <- est_Sp_OVERALL$mean*100
            est_Sp_OVERALL_median <- est_Sp_OVERALL$median*100
            ##
            print( abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3)) )
            message(paste("Sp (SUM of abs. diffs) = "))
            print( sum(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MEAN of abs. diffs) = "))
            print( mean(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MEDIAN of abs. diffs) = "))
            print( median(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MAX of abs. diffs) = "))
            print( max(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
        }
        })
        
        
        # ## Model_type <- 
        # ## Model_type <- "Cerullo_FIXED_cutpoints"
        # ## Model_type <- "Cerullo_RANDOM_HOMOG_cutpoints"
        # Model_type <- "Cerullo_RANDOM_cutpoints"

        
        if (Model_type == "Jones") { 
              colour <- "red"
              par(mfrow = c(2, 1))
              plot(  log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)    ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)    ;   abline(h = 0)
        } else if (Model_type == "Cerullo_FIXED_cutpoints") { 
              colour <- "orange"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        } else if (Model_type == "Cerullo_RANDOM_HOMOG_cutpoints") { 
              colour <- "blue"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        } else if (Model_type == "Cerullo_RANDOM_cutpoints") { 
              colour <- "green"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        }
        
        message(paste("Model_type = ", Model_type))

        
       #  # dev.off()
       #  plot( x = 100.0 - est_Sp_OVERALL_mean, y = est_Se_OVERALL_mean, col = "blue")
       #  lines(x = 100.0 - true_Sp_OVERALL_weighted,     y = true_Se_OVERALL_weighted,     col = "green")
       #  # 
       # 
       # 
       #  true_Sp_OVERALL_weighted
       #  ##
       #  true_Se_OVERALL_weighted
       #  try({ 
       #      Stan_mod_sample$summary(c("Se_MU"))  %>% print(n = 100)
       #      Stan_mod_sample$summary(c("Se_MED"))  %>% print(n = 100)
       #  })
       #  try({ 
       #      Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
       #      Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)
       #  })
       #  
       # ## Stan_mod_sample$summary(c("C"))  %>% print(n = 100)
       #  
       #  try({ 
       #    Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
       #  })
       #  try({ 
       #    Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("C_MED")) %>% print(n = 100)
       #  })
       #  try({ 
       #    Se_MU <- Stan_mod_sample$summary(c("Se_MU")) ##%>% print(n = 100)
       #    Se_MU <- Se_MU$mean
       #    true_Se_OVERALL_weighted - Se_MU*100
       #    # ##
       #    # Se_MED <- Stan_mod_sample$summary(c("Se_MED")) ##%>% print(n = 100)
       #    # Se_MED <- Se_MED$mean
       #    # true_Se_OVERALL_weighted - Se_MED*100
       #    # ##
       #    Se_EMP <- Stan_mod_sample$summary(c("Se_EMP")) ##%>% print(n = 100)
       #    Se_EMP <- Se_EMP$mean
       #    true_Se_OVERALL_weighted - Se_EMP*100
       #    # ##
       #    # Se_SIM_MED <- Stan_mod_sample$summary(c("Se_SIM_MED")) ##%>% print(n = 100)
       #    # Se_SIM_MED <- Se_SIM_MED$mean
       #    # true_Se_OVERALL_weighted - Se_SIM_MED*100
       #    # ##
       #    # Se_SIM_MU <- Stan_mod_sample$summary(c("Se_SIM_MU")) ##%>% print(n = 100)
       #    # Se_SIM_MU <- Se_SIM_MU$mean
       #    # true_Se_OVERALL_weighted - Se_SIM_MU*100
       #    ##
       #    true_Se_OVERALL_weighted
       #    true_Se_OVERALL_weighted[4:7]
       #  })
       #  ##
       #  try({ 
       #    Stan_mod_sample$summary(c("Sp_MU")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("Sp_MED")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("Sp_EMP")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM_MED")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM_MU")) %>% print(n = 100)
       #    true_Sp_OVERALL_weighted
       #    true_Sp_OVERALL_weighted[4:7]
       #  })
       #  ##
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("C_SD")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("C_empirical")) %>% print(n = 100)
       #  # })
       #  ##
       #  # try({
       #  #   Stan_mod_sample$summary(c("unc_C_normal_SD_sq"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("raw_MU"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("raw_SD"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("log_increment_SD"))  %>% print(n = 100)
       #    ##
       #  try({
       #    Stan_mod_sample$summary(c("category_means"))  %>% print(n = 100)
       #    Stan_mod_sample$summary(c("category_SDs"))  %>% print(n = 100)
       #    Stan_mod_sample$summary(c("kappa"))  %>% print(n = 100)
       #  })
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_mu_empirical"))  %>% print(n = 100)
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_mu_medians"))  %>% print(n = 100)
       #  # })
       #  # ##
       #  # # Stan_mod_sample$summary(c("prob_ord_mu"))  %>% print(n = 100)
       #  # # Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
       #  # ##
       #  # try({  
       #  # Stan_mod_sample$summary(c("alpha")) %>% print(n = 100)
       #  # alpha <- Stan_mod_sample$summary(c("alpha"),  quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975)))  %>% print(n = 100)
       #  # # alpha_raw <- Stan_mod_sample$summary(c("alpha_raw"))  %>% print(n = 100)
       #  # # 100*c(alpha$mean)
       #  # })
       #  # ##
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("kappa")) %>% print(n = 100)
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("unc_C_normal_SD")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("unc_C_normal_MU")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("unc_C_normal_MED")) %>% print(n = 100)
       #  # })
       #  # 
       #  # print(mean(inf_dir_samples$alpha))
       #  # print(quantile(inf_dir_samples$alpha, c(0.025, 0.50, 0.975)))
       #  # 
       #  # 
       #  # Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
       #  try({
       #    Stan_mod_sample$summary(c("location_MU_nd")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("location_MU_d")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("scale_MU_nd")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("scale_MU_d")) %>% print(n = 100)
       #  })
       #    ## 
       #  try({
       #    Stan_mod_sample$summary(c("expected_p_ord")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("dirichlet_cat_means_phi")) %>% print(n = 100)
       #  })
          ##
        try({
          Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
          ##
          Stan_mod_sample$summary(c("raw_scale_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("raw_scale_SD")) %>% print(n = 100)
          ##
          Stan_mod_sample$summary(c("lambda")) %>% print(n = 100)
          
        })
        try({
        Stan_mod_sample$summary(c("Se")) %>% print(n = 100)
        Stan_mod_sample$summary(c("Sp")) %>% print(n = 100)
        })
        try({
          Se <- Stan_mod_sample$summary(c("Se"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Se_pred <- Stan_mod_sample$summary(c("Se_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Sp <- Stan_mod_sample$summary(c("Sp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Sp_pred <- Stan_mod_sample$summary(c("Sp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Fp <- Stan_mod_sample$summary(c("Fp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Fp_pred <- Stan_mod_sample$summary(c("Fp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
        })
          
        
        
}
# Stan_mod_sample$summary(c("dirichlet_cat_SDs_sigma")) %>% print(n = 100)


# Stan_data_list$prior_beta_mu_mean
# Stan_data_list$prior_beta_mu_SD
# Stan_data_list$prior_beta_SD_mean
# Stan_data_list$prior_beta_SD_SD


df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, 
                  Sp_true = true_Sp_OVERALL_weighted/100, 
                  Fp_true = (100 - true_Sp_OVERALL_weighted)/100)


if (Model_type == "Cerullo_Gat_RANDOM_cutpoints") {
  
        ###
        df_Cerullo_Se <- tibble(Se_median = Se$`50%`, Se_lower = Se$`2.5%`, Se_upper = Se$`97.5%`, Se_pred_lower = Se_pred$`2.5%`, Se_pred_upper = Se_pred$`97.5%`)
        df_Cerullo_Sp <- tibble(Sp_median = Sp$`50%`, Sp_lower = Sp$`2.5%`, Sp_upper = Sp$`97.5%`, Sp_pred_lower = Sp_pred$`2.5%`, Sp_pred_upper = Sp_pred$`97.5%`)
        df_Cerullo_Fp <- tibble(Fp_median = Fp$`50%`, Fp_lower = Fp$`2.5%`, Fp_upper = Fp$`97.5%`, Fp_pred_lower = Fp_pred$`2.5%`, Fp_pred_upper = Fp_pred$`97.5%`)
        df_Cerullo <- tibble(cbind(df_Cerullo_Se, df_Cerullo_Sp, df_Cerullo_Fp)) %>%
          dplyr::mutate(Model = rep("Cerullo", n_thr))
        ##
        df_Cerullo <- tibble(cbind(df_true, df_Cerullo)) ; df_Cerullo
        
} else if (Model_type == "Jones") { 
  
        ###
        df_Jones_Se <- tibble(Se_median = Se$`50%`, Se_lower = Se$`2.5%`, Se_upper = Se$`97.5%`, Se_pred_lower = Se_pred$`2.5%`, Se_pred_upper = Se_pred$`97.5%`)
        df_Jones_Sp <- tibble(Sp_median = Sp$`50%`, Sp_lower = Sp$`2.5%`, Sp_upper = Sp$`97.5%`, Sp_pred_lower = Sp_pred$`2.5%`, Sp_pred_upper = Sp_pred$`97.5%`)
        df_Jones_Fp <- tibble(Fp_median = Fp$`50%`, Fp_lower = Fp$`2.5%`, Fp_upper = Fp$`97.5%`, Fp_pred_lower = Fp_pred$`2.5%`, Fp_pred_upper = Fp_pred$`97.5%`)
        df_Jones <- tibble(cbind(df_Jones_Se, df_Jones_Sp, df_Jones_Fp)) %>%
          dplyr::mutate(Model = rep("Jones", n_thr))
        ##
        df_Jones <- tibble(cbind(df_true, df_Jones)) ; df_Jones
        
}







Cerullo_polygon_Conf <- create_confidence_polygon(df = df_Cerullo, model_name = "Cerullo") ; Cerullo_polygon_Conf
Jones_polygon_Conf   <- create_confidence_polygon(df = df_Jones,   model_name = "Jones")   ; Jones_polygon_Conf
##

Cerullo_polygon_Pred <- create_prediction_polygon(df = df_Cerullo, model_name = "Cerullo") ; Cerullo_polygon_Pred
Jones_polygon_Pred   <- create_prediction_polygon(df = df_Jones,   model_name = "Jones")   ; Jones_polygon_Pred

df_Cerullo$Se_median
df_Cerullo$Se_true

## df_all <- df_Cerullo
##
df_all <- rbind(df_Jones, df_Cerullo)

df_Jones$Se_median - 
df_Cerullo$Se_median

## --------- Plot 1:
ggplot(data = df_all, 
       mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
  geom_line(size = 0.5) + 
  geom_point(size = 3) + 
  geom_point(color = "green", size = 3, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) +
  geom_line(color = "green",  size = 0.5, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
  theme_bw(base_size = 16)



## --------- Plot 2:
ggplot(data = df_all, 
       mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
  # geom_line(size = 0.5) + 
  # geom_point(size = 3) + 
  # geom_point(color = "green", size = 3, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) +
  # geom_line(color = "green",  size = 0.5, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
  theme_bw(base_size = 16) +
  geom_polygon(data = Cerullo_polygon_Conf, aes(x = x, y = y), fill = "red",  alpha = 0.75) +
  geom_polygon(data = Jones_polygon_Conf,   aes(x = x, y = y), fill = "cyan", alpha = 0.75) +
  geom_polygon(data = Cerullo_polygon_Pred, aes(x = x, y = y), fill = "red",  alpha = 0.25) + 
  geom_polygon(data = Jones_polygon_Pred,   aes(x = x, y = y), fill = "cyan", alpha = 0.25)
  

















x_individual_nd <- categorical_to_individual( Stan_data_list$x_cat[[1]], 
                                              binary_disease_indicator = 0, 
                                              missing_indicator = -1)

x_individual_d <- categorical_to_individual(  Stan_data_list$x_cat[[2]], 
                                              binary_disease_indicator = 1, 
                                              missing_indicator = -1)


x_individual_nd ## %>% print(n = 10000)
x_individual_d

x_individual <- tibble(rbind(x_individual_nd, x_individual_d))



require(MASS)
require(gridExtra) 


individual_obs_tibble <- x_individual_nd
group_name <- "non-diseased"
study_index <- 1


plot_distribution_fits(individual_obs_tibble = x_individual_nd,
                       group_name = "non-diseased")
                       



# Function to compare diseased and non-diseased distributions
compare_distributions <- function(x_cat_nd, x_cat_d, study_index = 1) {
  p1 <- plot_distribution_fits(x_cat_nd, study_index, "Non-diseased")
  p2 <- plot_distribution_fits(x_cat_d, study_index, "Diseased")
  
  grid.arrange(p1, p2, ncol = 1)
}

# Example usage for a specific study
# Assuming x_cat_nd is your non-diseased data and x_cat_d is your diseased data
# compare_distributions(x_cat[[1]], x_cat[[2]], 1)

# To examine multiple studies
examine_multiple_studies <- function(x_cat_nd, x_cat_d, study_indices = 1:3) {
  for (i in study_indices) {
    print(paste("Examining Study", i))
    print(compare_distributions(x_cat_nd, x_cat_d, i))
  }
}

# Example usage to examine first 3 studies
# examine_multiple_studies(x_cat[[1]], x_cat[[2]], 1:3)



































# Cerullo_polygon$Fp_lower
# Cerullo_polygon$Fp_upper
# 
# Cerullo_polygon$Se_lower
# Cerullo_polygon$Se_upper
# 
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Se_lower, ymax = Se_upper, fill = Model)) + 
#   # geom_ribbon(data = df_all, mapping  = aes(xmin = Fp_lower, xmax = Fp_upper, fill = Model))
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Se_pred_lower, ymax = Se_upper, fill = Model)) + 
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Fp_pred_lower, ymax = Fp_pred_upper, fill = Model))
# 
# 
#  geom_ribbon()
# 
# 
#   # geom_point(color = "red",   size = 4, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
#   # geom_line(color = "red",    size = 1, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
#   # theme_bw(base_size = 16)
#   
#   
# 
# 
# mu <- 0.50
# sigma <- 1.10
# 
# log_X <- rnorm(n = 100000, mean = mu, sd = sigma)
# 
# X <- exp(log_X)
# 
# plot(density(X), xlim  = c(0, 10))
#      
# 
# mu_X <- exp(mu + 0.5*sigma^2) ; mu_X ; mean(X)
# med_X <- exp(mu) ; med_X ; median(X)
# sd_X <- sqrt((exp(sigma^2) - 1) * exp(2*mu + sigma^2)) ; sd_X ; sd(X)
# 
# quantile_X_lower <- exp(mu + sigma*qnorm(0.16)) ; quantile_X_lower
# quantile_X_upper <- exp(mu + sigma*qnorm(1 - 0.16)) ; quantile_X_upper
# 
# exp(mu + sigma*qnorm(0.50))
# 
# sd_from_q_X <- (log(quantile_X_upper) - log(quantile_X_lower))/2
# 
# sd_from_q_X <- sqrt((quantile_X_upper-med_X)*(med_X-quantile_X_lower)) ; sd_from_q_X
# 
# 
# sd_X_from_quantile <- 
# 
# plot(density(rnorm(n = 10000, mean = med_X, sd = sd_X)))
# 
# true_Se_OVERALL[4:7]
# 
# 
# Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
# 
# 
# softmax <- function(x) { 
#    return(log(1 + exp(x)))
# }
# softmax_inv <- function(x) { 
#    return(log(- 1 + exp(x)))
# }
# 
# log1p_exp <- softmax
# 
# 
# 
# median_softmax_X_from_X_normal <- function( unc_C_normal_MU, 
#                                             unc_C_normal_SD) {
#   
#         if (unc_C_normal_MU > 5.0) {
#           
#                     unc_C_normal_MED = unc_C_normal_MU
#               
#         } else if (unc_C_normal_MU < -5.0) { 
#           
#                     unc_C_normal_MED = exp(unc_C_normal_MU);
#               
#         } else { 
#        
#               
#                     unc_C_normal_SD_sq <- unc_C_normal_SD * unc_C_normal_SD;
#                     approx_mean <- log1p_exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD_sq);
#                     
#                     lognormal_adjustment <- log(0.5 * unc_C_normal_SD_sq)
#                     
#                     unc_C_normal_MED <- approx_mean / exp(lognormal_adjustment);
#                     
#                     # base_median <- softmax(unc_C_normal_MU);
#                     # 
#                     # correction <- 0
#                     # if (unc_C_normal_MU >= 0) {
#                     #   
#                     #      ## // For μ ≥ 0, correction is smaller:
#                     #      correction <- 0.2 * exp(-0.5 * unc_C_normal_MU) * unc_C_normal_SD^2;
#                     #   
#                     # } else { 
#                     #   
#                     #      ## // For μ < 0, correction increases
#                     #      correction <- 0.3 * exp(0.2 * abs(unc_C_normal_MU)) * unc_C_normal_SD^2;
#                     #   
#                     # }
#                     # 
#                     # ## // Apply correction (add because softplus underestimates median for moderate μ)
#                     # unc_C_normal_MED <- base_median - 0.5 * correction;
#             
#         }
#   
#         return_list <- list(lognormal_adjustment = lognormal_adjustment,
#                             approx_mean = approx_mean,
#                             unc_C_normal_MED = unc_C_normal_MED)
#         
#         return(return_list)
#   
#   
# }
# 
# 
# 
# unc_C_normal_MU <- 1.5
# unc_C_normal_SD <- 2.5
# unc_C_normal <- rnorm(n = 100000, mean = unc_C_normal_MU, sd = unc_C_normal_SD)
# X <- unc_C_normal
# X_softmax <- softmax(unc_C_normal)
# 
# mean(X_softmax)
# median(X_softmax)
# 
# median_softmax_X_from_X_normal(unc_C_normal_MU = unc_C_normal_MU,
#                                unc_C_normal_SD = unc_C_normal_SD)
# 
# mean(X_softmax)
# median(X_softmax)
# 
# 
# 
# plot(density(X_softmax))
# 
# 
# 
# 
# 
# Y <- exp(unc_C_normal) ## Y is log-normal
# mean(Y)
# median(Y)
# 
# 
# mean_Y <- exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD^2 ) ; mean_Y
# median_Y <- exp(unc_C_normal_MU) ; median_Y
# 
# plot(density(Y))
# 
# plot(density(log(1 + exp(unc_C_normal))))
# 
# Y <- 1 + exp(unc_C_normal)
# mean(Y)
# median(Y)
# Y_m_one <- Y - 1 ## this is log-normal
# mean(Y_m_one)
# median(Y_m_one)
# 
# exp(unc_C_normal_MU) ; median(Y_m_one)
# exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2) ; mean(Y_m_one)
# 
# 
# 1 + exp(unc_C_normal_MU) ; median(Y)
# 1 + exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2) ; mean(Y)
# 
# log(1 + exp(unc_C_normal_MU)) ; median(log(1 + exp(unc_C_normal)))
# 
# log(1 + exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2)) ;   mean(log(1 + exp(unc_C_normal)))
# 
# mean(unc_C_normal)
# median(unc_C_normal)
# 
# 
# 
# mean_Ym1 <- exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD^2 ) ; mean_Ym1
# median_Ym1 <- exp(unc_C_normal_MU) ; median_Ym1
# 
# mean_Y <- 1 + mean_Ym1 ; mean_Y ; mean(Y)
# median_Y <- 1 + median_Ym1 ; median_Y ; median(Y)
# 
# log1p_exp_X <- log1p_exp(X)
# median(log1p_exp_X)
# mean(log1p_exp_X)
# plot(density(log1p_exp_X))
# ##
# plot(density(exp(X)), xlim = c(0, 10000))
# 
# 
# 
# soft_C_normal_MED <- log1p_exp(unc_C_normal_MU) ; soft_C_normal_MED
# median(log1p_exp_X)
# mean(log1p_exp_X)
# median(X)
# 
# 
# 
# log(median_Y)
# log(mean_Y)
# 
# 
# 
# 
#   correction <- -0.5 * unc_C_normal_SD^2;  correction  ##  // Approximate correction
#   
# log(exp(median(X_softmax)) - 1) - correction
# exp(log(exp(median(X_softmax)) - 1) - correction)
#   
# approx_MU <- log(exp(soft_C_normal_MED) - 1) + correction; approx_MU
# 
# 
# 
# 
# 
# # 
# # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # 
# # 
# # 
# # 
# # alpha_array <- array(dim = c(2, n_cat))
# # counter <- 1
# # for (k in 1:n_cat) { 
# # for (c in 1:2) {
# #     alpha_array[c, k] <- alpha$`50%`[counter]
# #     counter <- counter + 1
# #   }
# # }
# #  
# # 
# # p_ord_mu_nd <- alpha_array[1, ] / sum(alpha_array[1,])
# # p_ord_mu_d  <- alpha_array[2, ] / sum(alpha_array[2,])
# # 
# # 
# # 
# # 
# # est_Sp_OVERALL_mean
# # true_Sp_OVERALL
# # 
# # 100 - est_Se_OVERALL_mean
# # true_Se_OVERALL
# # 
# # 
# # 100 - est_Se_OVERALL_mean
# # true_Se_OVERALL
# # 
# # 
# # Stan_mod_sample$summary(c("lambda"))  %>% print(n = 100)
# # 
# # 
# # 
# 
# 
# 
# 
# # 10*c(alpha_raw$mean)
# # 
# # mean(c(alpha_raw$mean))
# # 1/11
# # 
# # 
# # ## Summary estimates:
# # Stan_mod_sample$summary(c("Se"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("Sp"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)
# # # ##
# # Stan_mod_sample$summary(c("log_scale_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("log_scale_SD"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("scale"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # ## Stan_mod_sample$summary(c("log_alpha"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("phi"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("kappa"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("C_nd"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("C_d"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("C_mu_empirical"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("lambda"))  %>% print(n = 100)
# # ##
# # # Stan_mod_sample$summary(c("phi_d"))  %>% print(n = 100)
# # # Stan_mod_sample$summary(c("phi_nd"))  %>% print(n = 100)
# # 
# # # ## Study-specific estimates:
# # # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # # Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)
# # 
# # Stan_mod_sample$summary(c("cumul_prob"))  %>% print(n = 100)
# # 
# # 
# # 
# # 
# # Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("phi"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("prob_ord_mu"))  %>% print(n = 100)
# # ##
# # df_alpha <- Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # alpha <- df_alpha$mean
# # alpha_nd <- alpha_d <- c()
# # ##
# # counter <- 0 
# # for (k in 1:n_thr) {
# #   for (c in 1:2) {
# #       counter <- counter + 1
# #       if (c == 1) alpha_nd[k] <- alpha[counter]
# #       if (c == 2) alpha_d[k] <- alpha[counter + 1]
# #   }
# # }
# #   
# #   
# # prob_ord_mu_sim_nd   <- array(NA, dim = c(n_thr, 1000))
# # prob_ord_mu_sim_d    <- array(NA, dim = c(n_thr, 1000))
# # C_mu_sim_nd   <- array(NA, dim = c(n_thr, 1000))
# # C_mu_sim_d    <- array(NA, dim = c(n_thr, 1000))
# # anchor_nd <- 0.0
# # anchor_d  <- 0.0
# # 
# # for (i in 1:1000) {
# #   
# #       prob_ord_mu_sim_nd[,i] <-  c(MCMCpack::rdirichlet(n = 1, alpha_nd))
# #       prob_ord_mu_sim_d[,i]  <-  c(MCMCpack::rdirichlet(n = 1, alpha_d))
# #       
# #       C_mu_sim_nd[1, i] =   anchor_nd - qlogis( plogis(anchor_nd - (-Inf)) -  prob_ord_mu_sim_nd[1 ,i]  )
# #       C_mu_sim_d[1, i]  =   anchor_d  - qlogis( plogis(anchor_d  - (-Inf)) -  prob_ord_mu_sim_d[1, i]  )
# #       for (k in 2:n_thr) {
# #         C_mu_sim_nd[k, i] =  anchor_nd - qlogis( plogis(anchor_nd - C_mu_sim_nd[k - 1, i]) - prob_ord_mu_sim_nd[k, i] );
# #         C_mu_sim_d[k, i]  =  anchor_d  - qlogis( plogis(anchor_d  - C_mu_sim_d[k - 1, i])  - prob_ord_mu_sim_d[k, i] );
# #       }
# #       
# #   
# # }
# # 
# # C_mu_nd <- rowMeans(C_mu_sim_nd)
# # C_mu_d  <- rowMeans(C_mu_sim_d)
# # 
# # 
# # plogis( beta_mu_nd -  C_mu_nd[1:9] )
# # 100.0 - true_Sp_OVERALL
# # 
# # plogis( beta_mu_d  -  C_mu_d[1:9] )
# # true_Se_OVERALL
# # 
# # 
# # ##
# # C_mu <- Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
# # C_mu_means <- C_mu$mean
# # C_mu <- array(dim = c(2, n_thr))
# # counter <- 0 
# # for (k in 1:n_thr) {
# #   for (c in 1:2) {
# #      counter <- counter + 1
# #      C_mu[c, k] <- C_mu_means[counter]
# #   }
# # }
# # C_mu
# # C_mu <- reorder_decreasing(C_mu)
# # ##
# # C_nd_inc <- Stan_mod_sample$summary(c("C_nd_inc"))  %>% print(n = 100)
# # C_d_inc <- Stan_mod_sample$summary(c("C_d_inc"))  %>% print(n = 100)
# # ##
# # beta_mu <- Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # ##
# # plogis( beta_mu_nd - C_nd_inc )
# # 100.0 - true_Sp_OVERALL
# # ##
# #  plogis( - beta_mu_d + C_d_inc$mean)
# # true_Se_OVERALL
# # ##
# # beta_mu <- Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # beta_mu_nd <- beta_mu$mean[1]
# # beta_mu_d  <- beta_mu$mean[2]
# # ##
# # plogis( beta_mu_nd - C_mu[1, ] )
# # 100.0 - true_Sp_OVERALL
# # ##
# # 1.0 - plogis( C_mu[2, ] - beta_mu_d)
# # true_Se_OVERALL
# # ##
# # Stan_mod_sample$summary(c("cutpoints_nd"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("cutpoints_d"))  %>% print(n = 100)
# # 
# # 
# # 
# # prob_cumul_mu_nd*100
# # true_Sp_OVERALL
# # 
# # 
# # prob_cumul_mu_d*100
# # true_Se_OVERALL
# # 
# # 
# # 
# # Stan_mod_sample$summary(c("beta"))  %>% print(n = 100)
# # 
# # 
# # 
# #  
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# ## prior for softplus raw_scale:
# samps <- rnorm(n = 10000, mu = 0.5, sd = 1.0)
# soft_samps <- log(1 + exp(samps))
# 
# round(quantile(soft_samps, probs = c(0.025, 0.50, 0.975)), 2)
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
