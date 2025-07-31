

 
 
# # # # # # # ##
rm(list = ls())


{


os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
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


}













#### ----
{
  model_parameterisation = "Jones"
  ##
  softplus <- FALSE
  ##
  box_cox <- TRUE
  cts <- TRUE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
}


{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
}


{
  model_parameterisation = "R&G"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
}


{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  TRUE
  ##"
  Dirichlet_random_effects_type <- "kappa"
}


{
  model_parameterisation = "R&G"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  TRUE
  ##"
  Dirichlet_random_effects_type <- "kappa"
}








{





os <- .Platform$OS.type


if (os == "unix") {
  user_root_dir <- Sys.getenv("PWD")
} else if (os == "windows") {
  user_root_dir <- Sys.getenv("USERPROFILE")
}
local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")


setwd(local_pkg_dir)

{
  require(MetaOrdDTA)
  require(RcppParallel)
  require(BayesMVP)
}

source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))



{
    
    n_chains <-  8 ##  parallel::detectCores() / 2
    ##
    if (random_thresholds)    n_iter   <- ceiling(1000/(n_chains/4))
    if (!(random_thresholds)) n_iter   <- ceiling(1000/(n_chains/4))
    n_iter
    ##
    n_burnin <- 500
    ##
   #   max_treedepth <- 9
    max_treedepth <- 10
    ##
    adapt_delta <- 0.65
    ##  adapt_delta <- 0.80

}


##
## n_studies <- 10
## n_studies <- 50
########  n_studies <- 5
## n_studies <- 25
######## n_studies <- 100
##
## N_per_study_mean <- 500
######## N_per_study_mean <- 5000
######## N_per_study_mean <- 10000
##
 N_per_study_mean <- 500  ; n_studies <- 10
# N_per_study_mean <- 500  ; n_studies <- 50
 
##
N_per_study_SD <- N_per_study_mean
## N_per_study_SD <- 1 ## N_per_study_mean



missing_thr <- TRUE
# missing_thr <- FALSE






## alpha_pi <- "flat"
alpha_pi <- "weak_info"


{
  
  true_Mean_prev <- 0.25 ;  true_SD_probit_prev <- 0.25
  ##
  scale_ALL_bs_SDs_by <- 0.50
  ##
   bs_het_C_nd_prob_scale <- 0.001      ; bs_het_C_d_prob_scale <- 0.001
  # bs_het_C_nd_prob_scale <- 0.03932276 ; bs_het_C_d_prob_scale <- 0.07712564 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
  ##
  bivariate_locations_nd_bs_het_SD <- 0.45 ; bivariate_locations_d_bs_het_SD <- 0.60 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
  ######## bivariate_locations_nd_bs_het_SD <- 0.25 ; bivariate_locations_d_bs_het_SD <- 0.40  ## ---- made up / arbitrary (ish)
  
}



index_test_chosen_index <- 1 + 1

# if (index_test_MA == "GAD_2") { 
#         
#         index_test_chosen_index <- 1 + 1
#         n_thr <- 6
#         ##
#         vec_index_inner <- 1:n_thr ## NULL
#         compute_sim_study_metrics <- NULL
#         
#         if (bs_het_C_nd_prob_scale < 0.02) {
#          
#             if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 400  ## based on MCSE(Bias) of 0.25% and 30 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
#             ##
#             if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 600 ## based on MCSE(Bias) of 0.25% and 30 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
#          
#         } else { 
#           
#             if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
#             ##
#             if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
#             
#         }
#   
# } else if (index_test_MA == "PHQ_9") { 
#   
#           index_test_chosen_index <- 4 + 1
#           n_thr <- 27
#           ##
#           vec_index_inner <- 5:18
#           n_inner_thr <- length(vec_index_inner)
#           vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
#           n_outer_thr <- length(vec_index_outer_thr)
#           # compute_sim_study_metrics <- 1
#           ##
#           if (bs_het_C_nd_prob_scale < 0.02) {
#             {
#               if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 10 sims.
#               if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
#               ##
#               if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 250 ## based on MCSE(Bias) of 0.25% and 10 sims.
#               if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 30  ## based on MCSE(Bias) of 0.25% and 10 sims.
#             }
#           } else { 
#             if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
#             ##
#             if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
#             if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
#           }
#           
# }


 

## 
if (n_studies < 50) { 
    n_sims <- 100
} else { 
    n_sims <- 25
}

##



{
    print(paste("--------------------------------------------------------------"))
    print(paste("n_sims = ", n_sims))
    print(paste("n_studies = ", n_studies))
    print(paste("N_per_study_mean = ", N_per_study_mean))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("bivariate_locations_nd_bs_het_SD = ", bivariate_locations_nd_bs_het_SD))
    print(paste("bivariate_locations_d_bs_het_SD = ", bivariate_locations_d_bs_het_SD))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("bs_het_C_nd_prob_scale = ", bs_het_C_nd_prob_scale))
    print(paste("bs_het_C_d_prob_scale = ", bs_het_C_d_prob_scale))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("model_parameterisation = ", model_parameterisation))
    print(paste("random_thresholds = ", random_thresholds))
    print(paste("Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))
    print(paste("alpha_pi = ", alpha_pi))
    ##
    print(paste("--------------------------------------------------------------"))
}






{
    output_dir = "simulation_results"
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # # Initialize storage for results
    # # Set up summary storage
    # summary_df <- create_summary_df(n_sims = n_sims,
    #                                 min_k = vec_index_inner[1],
    #                                 max_k = tail(vec_index_inner, 1))
    
    network <- TRUE
    ##
    softplus <- FALSE
    # softplus <- TRUE
}



# file_name <- "DTA_MA_Xu_FIXEDthr_daig_J.stan"
##
# use_custom_file <- TRUE
use_custom_file <- FALSE




}




# 
# source(file.path("inst", "examples", "MetaOrdDTA_MA_eg_1_debug_v4.R"))
# 
# min_ESS
# # priors
# alpha_pi
# 
# 
# 











