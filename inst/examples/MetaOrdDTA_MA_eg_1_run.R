

 
 
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
                      quick = TRUE,
                      dependencies = FALSE)
    ##
    ## May need to restart R:
    ##
    .rs.restartR()  # In RStudio
    
    # ?devtools::install


}
#
# # 
# 



  

start_sim_index <- 1



  # 
#### ----
{
  model_parameterisation = "Bivariate"
  ##
  softplus <- FALSE
  ##
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "none"
}


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
  
  seed <- 123

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
    n_iter   <- ceiling(1000/(n_chains/4)) ## rand-thr actually have better ESS/sec 
    n_iter
    ##
    n_burnin <- 500
    ##
    ##  max_treedepth <- 9
     max_treedepth <- 10
    # 2^max_treedepth
      ##  adapt_delta <- 0.80
      adapt_delta <- 0.65

}



  # N_per_study_mean <- 500  ; n_studies <- 10
  N_per_study_mean <- 500  ; n_studies <- 50
##
N_per_study_SD <- N_per_study_mean
## N_per_study_SD <- 1 ## N_per_study_mean




##
   ##     index_test_MA <- "GAD_2" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "GAD_2" ;  missing_thr <- FALSE   ## FOR TESTING
##
       # index_test_MA <- "HADS" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "HADS" ;  missing_thr <- FALSE   ## FOR TESTING
##
     index_test_MA <- "GAD_7" ;  missing_thr <- TRUE  ## SIM STUDY #2
######## index_test_MA <- "GAD_7" ;  missing_thr <- FALSE   ## FOR TESTING
  ##
 ##       index_test_MA <- "BAI"   ;  missing_thr <- TRUE  ## SIM STUDY #2


       alpha_pi <- "one"
  ##' alpha_pi <- "flat"
## alpha_pi <- "weak_info"


{
  true_Mean_prev <- 0.25 ;  true_SD_probit_prev <- 0.25
  ##
   ############# scale_ALL_bs_SDs_by <- 0.5
   ###   scale_ALL_bs_SDs_by <- 1.0
  ##
 ###   bs_het_C_nd_prob_scale <- 0.03932276 ; bs_het_C_d_prob_scale <- 0.07712564 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
}


if (index_test_MA == "GAD_2") {
        ##
        ## ---- based on "kappa" Xu model fitted to real GAD-2 data:
        ##
        # bivariate_locations_nd_bs_het_SD <- 0.343
        # bivariate_locations_d_bs_het_SD  <- 0.553 
        ##
        index_test_chosen_index <- 1 + 1
        n_thr <- 6
        ##
        vec_index_inner <- 1:n_thr ## NULL
        compute_sim_study_metrics <- NULL
 
          
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 1000 ## based on MCSE(Bias) of 0.25% and 30 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 150  ## based on MCSE(Bias) of 0.25% and 30 sims.
 
          
} else if (index_test_MA == "HADS") {
          ##
          ## ---- based on "kappa" Xu model fitted to real HADS data
          ##
          # bivariate_locations_nd_bs_het_SD <- 0.235
          # bivariate_locations_d_bs_het_SD  <- 0.396
          index_test_chosen_index <- 3 + 1
          n_thr <- 21
          ##
          vec_index_inner <- 3:17
          n_inner_thr <- length(vec_index_inner)
          vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
          n_outer_thr <- length(vec_index_outer_thr)
          # compute_sim_study_metrics <- 1
 
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
          
} else if (index_test_MA == "BAI") {
  
          ##
          ## ---- based on "kappa" Xu model fitted to real HADS data
          ##
          # bivariate_locations_nd_bs_het_SD <- 0.235
          # bivariate_locations_d_bs_het_SD  <- 0.396
          index_test_chosen_index <- 2 + 1
          n_thr <- 63
          ##
          vec_index_inner <- 3:43
          n_inner_thr <- length(vec_index_inner)
          vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
          n_outer_thr <- length(vec_index_outer_thr)
          # compute_sim_study_metrics <- 1
          
          if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
          if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
          ##
          if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
          if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
  
  
} else if (index_test_MA == "GAD_7") {
          ##
          ## ---- based on "kappa" Xu model fitted to real GAD-7 data:
          ##
          # bivariate_locations_nd_bs_het_SD <- 0.488
          # bivariate_locations_d_bs_het_SD  <- 0.669 
          index_test_chosen_index <- 4 + 1
          n_thr <- 21
          ##
          vec_index_inner <- 3:17
          n_inner_thr <- length(vec_index_inner)
          vec_index_outer_thr <- setdiff(1:n_thr, vec_index_inner)
          n_outer_thr <- length(vec_index_outer_thr)
          # compute_sim_study_metrics <- 1
          ##
            if ((n_studies == 10) && (N_per_study_mean == 500))  n_sims <- 2000 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 500))  n_sims <- 200  ## based on MCSE(Bias) of 0.25% and 10 sims.
            ##
            if ((n_studies == 10) && (N_per_study_mean == 2500)) n_sims <- 500 ## based on MCSE(Bias) of 0.25% and 10 sims.
            if ((n_studies == 50) && (N_per_study_mean == 2500)) n_sims <- 100  ## based on MCSE(Bias) of 0.25% and 10 sims.
  
}




 

## 
if (n_studies < 50) { 
    n_sims <- 250
} else { 
    n_sims <- 50
}

##



{
    print(paste("--------------------------------------------------------------"))
    print(paste("n_sims = ", n_sims))
    print(paste("n_studies = ", n_studies))
    print(paste("N_per_study_mean = ", N_per_study_mean))
    ##
    # print(paste("--------------------------------------------------------------"))
    # print(paste("bivariate_locations_nd_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD))
    # print(paste("bivariate_locations_d_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD))
    # ##
    # print(paste("--------------------------------------------------------------"))
    # if (bs_het_C_nd_prob_scale > 0.01) { 
    #     cutpoint_SD <- "significant het"
    # } else { 
    #     cutpoint_SD <- "No het"
    # }
    # print(paste("cutpoint_SD = ", cutpoint_SD))
    ##
    # print(paste("--------------------------------------------------------------"))
    # print(paste("model_parameterisation = ", model_parameterisation))
    # print(paste("--------------------------------------------------------------"))
    # print(paste("random_thresholds = ", random_thresholds))
    print(paste("--------------------------------------------------------------"))
    print(paste("index_test_MA = ", index_test_MA))
    print(paste("--------------------------------------------------------------"))
    # print(paste("Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))
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
  #                                 max_k = tail(vec_index_inner, 1),
  #                                 start_index = start_sim_index)
  
  network <- FALSE
  ##
  softplus <- FALSE
  # softplus <- TRUE
  
  
}



# file_name <- "DTA_MA_Xu_FIXEDthr_daig_J.stan"
##
# use_custom_file <- TRUE
use_custom_file <- FALSE



# source(file.path(local_pkg_dir, "inst", "examples", "NMA_missing_thr_prep_data.R"))


}




# if (n_studies == 50) {
#       if (index_test_MA == "GAD_2") { 
#         
#               sim_config <- list(
#                 index_test_MA = "GAD_2",
#                 vec_index_inner = 1:6,
#                 N_per_study_mean = 500,
#                 n_studies = 50,
#                 alpha_pi = "flat",
#                 n_chains = n_chains,
#                 n_iter = ceiling(1000/(n_chains/4)),
#                 n_burnin = 500,
#                 adapt_delta = 0.65,
#                 max_treedepth = 10,
#                 network = FALSE,
#                 use_custom_file = FALSE,
#                 file_name = NULL,
#                 missing_thr = TRUE
#               )
#               ##
#               # Define your model configurations
#               model_configs <- list(
#                 list(model_parameterisation = "Bivariate",
#                      start_sim_index = 16,
#                      n_sims = 450,
#                      box_cox = FALSE,  
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Jones",    
#                      start_sim_index = 16,
#                      n_sims = 500,
#                      box_cox = TRUE,  
#                      softplus = FALSE, 
#                      cts = TRUE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Xu", 
#                      start_sim_index = 16,
#                      n_sims = 400,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "R&G",   
#                      start_sim_index = 16,
#                      n_sims = 550,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "Xu",   
#                      start_sim_index = 16,
#                      n_sims = 350,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa"),
#                 ##
#                 list(model_parameterisation = "R&G", 
#                      start_sim_index = 16,
#                      n_sims = 400,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa")
#               )
#               
#               
#       } else if (index_test_MA == "HADS") { 
#         
#               sim_config <- list(
#                 index_test_MA = "HADS",
#                 vec_index_inner = 3:17,
#                 N_per_study_mean = 500,
#                 n_studies = 50,
#                 alpha_pi = "flat",
#                 n_chains = n_chains,
#                 n_iter = ceiling(1000/(n_chains/4)),
#                 n_burnin = 500,
#                 adapt_delta = 0.65,
#                 max_treedepth = 10,
#                 network = FALSE,
#                 use_custom_file = FALSE,
#                 file_name = NULL,
#                 missing_thr = TRUE
#               )
#               ##
#               # Define your model configurations
#               model_configs <- list(
#                 list(model_parameterisation = "Bivariate",
#                      start_sim_index = 16,
#                      n_sims = 200,
#                      box_cox = FALSE,  
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Jones",    
#                      start_sim_index = 16,
#                      n_sims = 150,
#                      box_cox = TRUE,  
#                      softplus = FALSE, 
#                      cts = TRUE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Xu", 
#                      start_sim_index = 16,
#                      n_sims = 100,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "R&G",   
#                      start_sim_index = 16,
#                      n_sims = 200,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "Xu",   
#                      start_sim_index = 16,
#                      n_sims = 100,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa"),
#                 ##
#                 list(model_parameterisation = "R&G", 
#                      start_sim_index = 16,
#                      n_sims = 150,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa")
#               )
#         
#         
#       } else if (index_test_MA == "GAD_7") { 
#               
#               sim_config <- list(
#                 index_test_MA = "GAD_7",
#                 vec_index_inner = 3:17,
#                 N_per_study_mean = 500,
#                 n_studies = 50,
#                 alpha_pi = "flat",
#                 n_chains = n_chains,
#                 n_iter = ceiling(1000/(n_chains/4)),
#                 n_burnin = 500,
#                 adapt_delta = 0.65,
#                 max_treedepth = 10,
#                 network = FALSE,
#                 use_custom_file = FALSE,
#                 file_name = NULL,
#                 missing_thr = TRUE
#               )
#               ##
#               # Define your model configurations
#               model_configs <- list(
#                 list(model_parameterisation = "Bivariate",
#                      start_sim_index = 16,
#                      n_sims = 450,
#                      box_cox = FALSE,  
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Jones",    
#                      start_sim_index = 16,
#                      n_sims = 300,
#                      box_cox = TRUE,  
#                      softplus = FALSE, 
#                      cts = TRUE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "none"),
#                 ##
#                 list(model_parameterisation = "Xu", 
#                      start_sim_index = 16,
#                      n_sims = 400,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "R&G",   
#                      start_sim_index = 16,
#                      n_sims = 450,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = FALSE,
#                      Dirichlet_random_effects_type = "fixed"),
#                 ##
#                 list(model_parameterisation = "Xu",   
#                      start_sim_index = 16,
#                      n_sims = 350,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa"),
#                 ##
#                 list(model_parameterisation = "R&G", 
#                      start_sim_index = 16,
#                      n_sims = 350,
#                      box_cox = FALSE, 
#                      softplus = FALSE, 
#                      cts = FALSE, 
#                      random_thresholds = TRUE,
#                      Dirichlet_random_effects_type = "kappa")
#               )
#         
#       }
# } else if (n_studies == 10) { 
#   
#   if (index_test_MA == "GAD_2") { 
#         
#         sim_config <- list(
#           index_test_MA = "GAD_2",
#           vec_index_inner = 1:6,
#           N_per_study_mean = 500,
#           n_studies = 10,
#           alpha_pi = "flat",
#           n_chains = n_chains,
#           n_iter = ceiling(1000/(n_chains/4)),
#           n_burnin = 500,
#           adapt_delta = 0.65,
#           max_treedepth = 10,
#           network = FALSE,
#           use_custom_file = FALSE,
#           file_name = NULL,
#           missing_thr = TRUE
#         )
#         ##
#         # Define your model configurations
#         model_configs <- list(
#           list(model_parameterisation = "Bivariate",
#                start_sim_index = 51,
#                n_sims = 1600,
#                box_cox = FALSE,  
#                softplus = FALSE, 
#                cts = FALSE, 
#                random_thresholds = FALSE,
#                Dirichlet_random_effects_type = "none"),
#           ##
#           list(model_parameterisation = "Jones",    
#                start_sim_index = 51,
#                n_sims = 1600,
#                box_cox = TRUE,  
#                softplus = FALSE, 
#                cts = TRUE, 
#                random_thresholds = FALSE,
#                Dirichlet_random_effects_type = "none"),
#           ##
#           list(model_parameterisation = "Xu", 
#                start_sim_index = 51,
#                n_sims = 1500,
#                box_cox = FALSE, 
#                softplus = FALSE, 
#                cts = FALSE, 
#                random_thresholds = FALSE,
#                Dirichlet_random_effects_type = "fixed"),
#           ##
#           list(model_parameterisation = "R&G",   
#                start_sim_index = 51,
#                n_sims = 2100,
#                box_cox = FALSE, 
#                softplus = FALSE, 
#                cts = FALSE, 
#                random_thresholds = FALSE,
#                Dirichlet_random_effects_type = "fixed"),
#           ##
#           list(model_parameterisation = "Xu",   
#                start_sim_index = 51,
#                n_sims = 1400,
#                box_cox = FALSE, 
#                softplus = FALSE, 
#                cts = FALSE, 
#                random_thresholds = TRUE,
#                Dirichlet_random_effects_type = "kappa"),
#           ##
#           list(model_parameterisation = "R&G", 
#                start_sim_index = 51,
#                n_sims = 1800,
#                box_cox = FALSE, 
#                softplus = FALSE, 
#                cts = FALSE, 
#                random_thresholds = TRUE,
#                Dirichlet_random_effects_type = "kappa")
#         )
#         
#     
#   } else if (index_test_MA == "HADS") { 
#     
#     sim_config <- list(
#       index_test_MA = "HADS",
#       vec_index_inner = 3:17,
#       N_per_study_mean = 500,
#       n_studies = 10,
#       alpha_pi = "flat",
#       n_chains = n_chains,
#       n_iter = ceiling(1000/(n_chains/4)),
#       n_burnin = 500,
#       adapt_delta = 0.65,
#       max_treedepth = 10,
#       network = FALSE,
#       use_custom_file = FALSE,
#       file_name = NULL,
#       missing_thr = TRUE
#     )
#     ##
#     # Define your model configurations
#     model_configs <- list(
#       # ##
#       # list(model_parameterisation = "Bivariate",
#       #      start_sim_index = 1,
#       #      n_sims = 700,
#       #      box_cox = FALSE,  
#       #      softplus = FALSE, 
#       #      cts = FALSE, 
#       #      random_thresholds = FALSE,
#       #      Dirichlet_random_effects_type = "none"),
#       # ##
#       # list(model_parameterisation = "Jones",    
#       #      start_sim_index = 1,
#       #      n_sims = 500,
#       #      box_cox = TRUE,  
#       #      softplus = FALSE, 
#       #      cts = TRUE, 
#       #      random_thresholds = FALSE,
#       #      Dirichlet_random_effects_type = "none"),
#       # ##
#       # list(model_parameterisation = "Xu", 
#       #      start_sim_index = 1,
#       #      n_sims = 400,
#       #      box_cox = FALSE, 
#       #      softplus = FALSE, 
#       #      cts = FALSE, 
#       #      random_thresholds = FALSE,
#       #      Dirichlet_random_effects_type = "fixed"),
#       # ##
#       list(model_parameterisation = "R&G",   
#            start_sim_index = 1,
#            n_sims = 650,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = FALSE,
#            Dirichlet_random_effects_type = "fixed"),
#       # ##
#       # list(model_parameterisation = "Xu",   
#       #      start_sim_index = 1,
#       #      n_sims = 400,
#       #      box_cox = FALSE, 
#       #      softplus = FALSE, 
#       #      cts = FALSE, 
#       #      random_thresholds = TRUE,
#       #      Dirichlet_random_effects_type = "kappa"),
#       # ##
#       list(model_parameterisation = "R&G", 
#            start_sim_index = 1,
#            n_sims = 600,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = TRUE,
#            Dirichlet_random_effects_type = "kappa")
#     )
#     
#     
#   } else if (index_test_MA == "GAD_7") { 
#     
#     sim_config <- list(
#       index_test_MA = "GAD_7",
#       vec_index_inner = 3:17,
#       N_per_study_mean = 500,
#       n_studies = 10,
#       alpha_pi = "flat",
#       n_chains = n_chains,
#       n_iter = ceiling(1000/(n_chains/4)),
#       n_burnin = 500,
#       adapt_delta = 0.65,
#       max_treedepth = 10,
#       network = FALSE,
#       use_custom_file = FALSE,
#       file_name = NULL,
#       missing_thr = TRUE
#     )
#     ##
#     # Define your model configurations
#     model_configs <- list(
#       ##
#       list(model_parameterisation = "Bivariate",
#            start_sim_index = 51,
#            n_sims = 2200,
#            box_cox = FALSE,  
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = FALSE,
#            Dirichlet_random_effects_type = "none"),
#       ##
#       list(model_parameterisation = "Jones",    
#            start_sim_index = 51,
#            n_sims = 1300,
#            box_cox = TRUE,  
#            softplus = FALSE, 
#            cts = TRUE, 
#            random_thresholds = FALSE,
#            Dirichlet_random_effects_type = "none"),
#       ##
#       list(model_parameterisation = "Xu", 
#            start_sim_index = 51,
#            n_sims = 1300,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = FALSE,
#            Dirichlet_random_effects_type = "fixed"),
#       ##
#       list(model_parameterisation = "R&G",   
#            start_sim_index = 51,
#            n_sims = 2500,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = FALSE,
#            Dirichlet_random_effects_type = "fixed"),
#       ##
#       list(model_parameterisation = "Xu",   
#            start_sim_index = 51,
#            n_sims = 1300,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = TRUE,
#            Dirichlet_random_effects_type = "kappa"),
#       ##
#       list(model_parameterisation = "R&G", 
#            start_sim_index = 51,
#            n_sims = 1900,
#            box_cox = FALSE, 
#            softplus = FALSE, 
#            cts = FALSE, 
#            random_thresholds = TRUE,
#            Dirichlet_random_effects_type = "kappa")
#     )
#     
#   }
# }
# 






# model_configs
# sim_config

# 
# sim_configs <- list( 
#   list( index_test_MA = "GAD_2",
#         vec_index_inner = 1:6,
#         N_per_study_mean = 500,
#         n_studies = 10,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   # ##
#   # list( index_test_MA = "GAD_2",
#   #       vec_index_inner = 1:6,
#   #       N_per_study_mean = 500,
#   #       n_studies = 10,
#   #       alpha_pi = "one",
#   #       n_chains = n_chains,
#   #       n_iter = ceiling(1000/(n_chains/4)),
#   #       n_burnin = 500,
#   #       adapt_delta = 0.65,
#   #       max_treedepth = 10,
#   #       network = FALSE,
#   #       use_custom_file = FALSE,
#   #       file_name = NULL,
#   #       missing_thr = TRUE),
#   #
#   # list( index_test_MA = "GAD_2",
#   #       vec_index_inner = 1:6,
#   #       N_per_study_mean = 500,
#   #       n_studies = 50,
#   #       alpha_pi = "one",
#   #       n_chains = n_chains,
#   #       n_iter = ceiling(1000/(n_chains/4)),
#   #       n_burnin = 500,
#   #       adapt_delta = 0.65,
#   #       max_treedepth = 10,
#   #       network = FALSE,
#   #       use_custom_file = FALSE,
#   #       file_name = NULL,
#   #       missing_thr = TRUE),
#   ##
#   ##
#   ##
#   # list( index_test_MA = "GAD_2",
#   #       vec_index_inner = 1:6,
#   #       N_per_study_mean = 500,
#   #       n_studies = 50,
#   #       alpha_pi = "one",
#   #       n_chains = n_chains,
#   #       n_iter = ceiling(1000/(n_chains/4)),
#   #       n_burnin = 500,
#   #       adapt_delta = 0.65,
#   #       max_treedepth = 10,
#   #       network = FALSE,
#   #       use_custom_file = FALSE,
#   #       file_name = NULL,
#   #       missing_thr = TRUE),
#   ##
#   list( index_test_MA = "GAD_2",
#         vec_index_inner = 1:6,
#         N_per_study_mean = 500,
#         n_studies = 10,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   list( index_test_MA = "GAD_2",
#         vec_index_inner = 1:6,
#         N_per_study_mean = 500,
#         n_studies = 10,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   ## ---- 
#   ##
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 10,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   #
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 10,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 50,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 50,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 50,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE),
#   ##
#   list( index_test_MA = "HADS",
#         vec_index_inner = 3:17,
#         N_per_study_mean = 500,
#         n_studies = 50,
#         alpha_pi = "one",
#         n_chains = n_chains,
#         n_iter = ceiling(1000/(n_chains/4)),
#         n_burnin = 500,
#         adapt_delta = 0.65,
#         max_treedepth = 10,
#         network = FALSE,
#         use_custom_file = FALSE,
#         file_name = NULL,
#         missing_thr = TRUE)
#   )
# 
# 
# 

 


# 
# overall_sim_data %>% filter(DGM == "GAD-2") %>% filter(n_studies == 10) %>% print(n = 100)
# overall_sim_data %>% filter(DGM == "GAD-2") %>% filter(n_studies == 50) %>% print(n = 100)
# ##
# overall_sim_data %>% filter(DGM == "HADS")  %>% filter(n_studies == 10) %>% print(n = 100)
# overall_sim_data %>% filter(DGM == "HADS")  %>% filter(n_studies == 50) %>% print(n = 100)
# ##
# overall_sim_data %>% filter(DGM == "BAI")  %>% filter(n_studies == 10) %>% print(n = 100)
# overall_sim_data %>% filter(DGM == "BAI")  %>% filter(n_studies == 50) %>% print(n = 100)
# ##
# # overall_sim_data %>% filter(DGM == "GAD-7") %>% filter(n_studies == 10) %>% print(n = 100)
# # overall_sim_data %>% filter(DGM == "GAD-7") %>% filter(n_studies == 50) %>% print(n = 100)

 
 







n_chains <- 4

start_sim_index <- 6000

# 2141 - 1600 # 550
# 2029 - 1578 # 450
# 1973 - 1299 # 675
# ##
# 440 - 286 # 150
# 410 - 345 # 70
# ##
# 754 - 500 # 250
# 648 - 449 # 200
# ##


sim_configs <- list( 
  # ##
  # ## ---- GAD_2, n_studies = 10:
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 10,
  #       model_parameterisation = "Jones", 
  #       alpha_pi = "flat",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 550,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 10,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 450,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 10,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 675,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa"),
  ##
  ## ---- GAD_2, n_studies = 50:
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 50,
  #       model_parameterisation = "Jones", 
  #       alpha_pi = "flat",
  #       start_sim_index = 3000,
  #       n_sims = hello,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 50,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 150,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "GAD_2",
  #       vec_index_inner = 1:6,
  #       n_studies = 50,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 70,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa"),
  # ##
  # ## ---- HADS, n_studies = 10:
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 10,
  #       model_parameterisation = "Jones", 
  #       alpha_pi = "flat",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 250,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 10,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = 3000,
  #       n_sims = hello,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 10,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 200,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa")
  ##
  ## ---- HADS, n_studies = 50:
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 50,
  #       model_parameterisation = "Jones", 
  #       alpha_pi = "flat",
  #       start_sim_index = 3000,
  #       n_sims = hello,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 50,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = 3000,
  #       n_sims = hello,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "HADS",
  #       vec_index_inner = 3:17,
  #       n_studies = 50,
  #       model_parameterisation = "Xu", 
  #       alpha_pi = "one",
  #       start_sim_index = 3000,
  #       n_sims = hello,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa"),
  # # ##
  # # ## ---- BAI, n_studies = 10:
  # # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 10,
  #       model_parameterisation = "Bivariate",
  #       alpha_pi = "flat",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 245,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "none"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 10,
  #       model_parameterisation = "Jones",
  #       alpha_pi = "flat",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 315,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "none"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 10,
  #       model_parameterisation = "Xu",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 217,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  list( N_per_study_mean = 500,
        n_chains = n_chains,
        n_iter = ceiling(1000/(n_chains/4)),
        n_burnin = 500,
        adapt_delta = 0.65,
        max_treedepth = 10,
        network = FALSE,
        use_custom_file = FALSE,
        file_name = NULL,
        missing_thr = TRUE,
        ##
        index_test_MA = "BAI",
        vec_index_inner = 3:43,
        n_studies = 10,
        model_parameterisation = "Xu",
        alpha_pi = "one",
        start_sim_index = start_sim_index,
        n_sims =  start_sim_index +  400,
        box_cox = FALSE,
        softplus = FALSE,
        cts = FALSE,
        random_thresholds = TRUE,
        Dirichlet_random_effects_type = "kappa"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 10,
  #       model_parameterisation = "R&G",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 600,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 10,
  #       model_parameterisation = "R&G",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims =  start_sim_index + 600,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa"),
  ##
  ## ---- BAI, n_studies = 50:
  ##
  list( N_per_study_mean = 500,
        n_chains = n_chains,
        n_iter = ceiling(1000/(n_chains/4)),
        n_burnin = 500,
        adapt_delta = 0.65,
        max_treedepth = 10,
        network = FALSE,
        use_custom_file = FALSE,
        file_name = NULL,
        missing_thr = TRUE,
        ##
        index_test_MA = "BAI",
        vec_index_inner = 3:43,
        n_studies = 50,
        model_parameterisation = "Bivariate",
        alpha_pi = "flat",
        start_sim_index = start_sim_index,
        n_sims =  start_sim_index + 150,
        box_cox = FALSE,
        softplus = FALSE,
        cts = FALSE,
        random_thresholds = FALSE,
        Dirichlet_random_effects_type = "none"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 50,
  #       model_parameterisation = "Jones",
  #       alpha_pi = "flat",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 63,
  #       box_cox = TRUE,
  #       softplus = FALSE,
  #       cts = TRUE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "none"),
  ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 50,
  #       model_parameterisation = "Xu",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 86,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # #
  list( N_per_study_mean = 500,
        n_chains = n_chains,
        n_iter = ceiling(1000/(n_chains/4)),
        n_burnin = 500,
        adapt_delta = 0.65,
        max_treedepth = 10,
        network = FALSE,
        use_custom_file = FALSE,
        file_name = NULL,
        missing_thr = TRUE,
        ##
        index_test_MA = "BAI",
        vec_index_inner = 3:43,
        n_studies = 50,
        model_parameterisation = "Xu",
        alpha_pi = "one",
        start_sim_index = start_sim_index,
        n_sims =  start_sim_index + 100,
        box_cox = FALSE,
        softplus = FALSE,
        cts = FALSE,
        random_thresholds = TRUE,
        Dirichlet_random_effects_type = "kappa")
  # #
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 50,
  #       model_parameterisation = "R&G",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims = start_sim_index + 150,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = FALSE,
  #       Dirichlet_random_effects_type = "fixed"),
  # ##
  # list( N_per_study_mean = 500,
  #       n_chains = n_chains,
  #       n_iter = ceiling(1000/(n_chains/4)),
  #       n_burnin = 500,
  #       adapt_delta = 0.65,
  #       max_treedepth = 10,
  #       network = FALSE,
  #       use_custom_file = FALSE,
  #       file_name = NULL,
  #       missing_thr = TRUE,
  #       ##
  #       index_test_MA = "BAI",
  #       vec_index_inner = 3:43,
  #       n_studies = 50,
  #       model_parameterisation = "R&G",
  #       alpha_pi = "one",
  #       start_sim_index = start_sim_index,
  #       n_sims =  start_sim_index + 120,
  #       box_cox = FALSE,
  #       softplus = FALSE,
  #       cts = FALSE,
  #       random_thresholds = TRUE,
  #       Dirichlet_random_effects_type = "kappa")
)
  


# complete_results_list[[1]]


# length(model_configs)
length(sim_configs)


# n_sims <- 15
source(file.path("inst", "examples", "MetaOrdDTA_MA_eg_1_debug_v4.R"))




    
    
    {
            seed <- 123
            
            source(file.path(local_pkg_dir, "inst", "examples", "MetaOrdDTA_MA_eg_1_debug_v4.R"))
            source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"))
            
            # First, check what's actually being saved
            ls()  # See all objects in current environmentexists("index_test_MA")  # Check if it exists
            print(index_test_MA)  # See its value
            
            # Make sure ALL required variables are defined before saving
            if (!exists("COMMON_SEED")) COMMON_SEED <- 123
            # if (!exists("index_test_MA")) stop("index_test_MA not defined!")
            
            current_wd <- getwd()
            
            # # Save workspace with explicit check
             save.image(file = "simulation_environment.RData")
            
            # Verify what was saved
            saved_objects <- load("simulation_environment.RData", envir = new.env())
            # print(paste("Saved objects include index_test_MA:", "index_test_MA" %in% saved_objects))

            
            library(doParallel)
            
            cl <- makeCluster(length(sim_configs), outfile = "")
            registerDoParallel(cl)
            
            
            if (exists("ii")) rm(ii)
            if (exists("n_sims")) cat(sprintf("n_sims = %d\n", n_sims))
            if (exists("start_sim_index")) cat(sprintf("start_sim_index = %d\n", start_sim_index))
            
            
            # Define output_dir in global environment
            output_dir <- normalizePath("simulation_results", mustWork = FALSE)
            
            # Create it if it doesn't exist
            if (!dir.exists(output_dir)) {
              dir.create(output_dir, recursive = TRUE)
            }
            
            
            # Export variables but NOT 'c' or 'result'
            # Export EVERYTHING from current environment
            clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)

            
        
        
          n_configs <- length(sim_configs)
          n_configs
        
          
          # sim_config
          
          ii <- 1
    }





results <- foreach(ii = 1:n_configs,
                   .combine = 'list',
                   .errorhandling = "pass") %dopar% {
                     
                     # Use a different variable name to avoid conflicts
                     sim_result <- list(worker_id = ii, 
                                        status = "starting", 
                                        start_time = Sys.time())
                     
                     tryCatch({
                       # # # Load environment
                       #   load("simulation_environment.RData")
                       #   setwd(current_wd)
                 
                       # Load packages
                       library(MetaOrdDTA)
                       library(cmdstanr)
                       library(posterior)
                       library(dplyr)
                       
                       config <- sim_configs[[ii]]
                       
                       # Set variables
                       seed <- COMMON_SEED
                       sim_index <- ii
                       set.seed(seed)
                       options(mc.cores = n_chains)
                       
                       # Extract config values
                       start_sim_index <- config$start_sim_index
                       n_sims <- config$n_sims
                       ##
                       model_parameterisation <- config$model_parameterisation
                       box_cox <- config$box_cox
                       softplus <- config$softplus
                       random_thresholds <- ifelse(is.null(config$random_thresholds), FALSE, config$random_thresholds)
                       Dirichlet_random_effects_type <- config$Dirichlet_random_effects_type
                       cts <- config$cts
                       
                       
                       cat(sprintf("\n===== Worker %d: Running %s =====\n", ii, model_parameterisation))
                       
                       # Run simulation
                       source(file.path(local_pkg_dir, "inst", "examples", "MetaOrdDTA_MA_eg_1_debug_v4.R"))
                       source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"))
                       ##
                       output <- run_single_simulation(start_sim_index = start_sim_index,
                                                       n_sims = n_sims,
                                                       model_parameterisation = model_parameterisation,
                                                       softplus = softplus,
                                                       box_cox = box_cox,
                                                       cts = cts,
                                                       random_thresholds = random_thresholds,
                                                       Dirichlet_random_effects_type = Dirichlet_random_effects_type,
                                                       local_pkg_dir = local_pkg_dir,
                                                       ##
                                                       sim_config = sim_configs[[ii]],
                                                       model_config = sim_configs[[ii]])
                       
                       # Update our result object
                       sim_result$status <- "completed"
                       sim_result$model <- model_parameterisation
                       sim_result$config <- config
                       sim_result$end_time <- Sys.time()
                       sim_result$duration_mins <- as.numeric(difftime(sim_result$end_time, sim_result$start_time, units = "mins"))
                       
                       cat(sprintf("\n===== Worker %d: Completed in %.1f minutes =====\n", ii, sim_result$duration_mins))
                       
                     }, error = function(e) {
                       cat(sprintf("\nWorker %d ERROR: %s\n", ii, e$message))
                       
                       sim_result$status <- "error"
                       sim_result$error_message <- e$message
                       sim_result$error_call <- toString(e$call)
                     })
                     
                     # Return the result
                     sim_result
                   }

stopCluster(cl)

# Check results
cat("\n\n========== FINAL RESULTS ==========\n")
for (i in seq_along(results)) {
  r <- results[[i]]
  if (is.list(r) && !is.null(r$worker_id)) {
    cat(sprintf("\nWorker %d: %s\n", r$worker_id, r$status))
    if (!is.null(r$model)) {
      cat(sprintf("  Model: %s\n", r$model))
    }
    if (!is.null(r$duration_mins)) {
      cat(sprintf("  Duration: %.1f minutes\n", r$duration_mins))
    }
    if (!is.null(r$error_message)) {
      cat(sprintf("  Error: %s\n", r$error_message))
    }
  }
}





# 
# {
#   # ##
#   # ## ---- Load NMA data:
#   # ##
#   # setwd(local_pkg_dir)
#   # setwd("/")
#   ## print(paste(getwd()))
#   ##
#   try({ 
#     rm(true_DGM_Se)
#     rm(true_DGM_Sp)
#   })
#   source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA.R"), local = TRUE)
#   source(file.path(local_pkg_dir, "inst", "examples", "NMA_missing_thr_prep_data.R"), local = TRUE)
#   source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"), , local = TRUE)
#   ##
#   # data <- readRDS(file.path(path, "inst", "examples", "data_example_1_NMA_list.RDS"))
#   ##
#   # x_NMA <- data$x_NMA
#   # ##
#   # x_MA <- list()
#   # for (c in 1:2) {
#   #   x_MA[[c]] <- x_NMA[[c]][[index_test_chosen_index - 1]]
#   # }
#   # ##
#   # ## Select which data ("x") to use:
#   # ##
#   # str(x_NMA)
#   ##
#   ## ---- Select x to use:
#   ##
#   if (index_test_MA == "GAD_7") {
#     if (missing_thr == TRUE) x <- x_GAD7_w_missing_thr
#     else                     x <- x_GAD7
#   } else if (index_test_MA == "GAD_2") {
#     if (missing_thr == TRUE) x <- x_GAD2_w_missing_thr
#     else                     x <- x_GAD2
#   } else if (index_test_MA == "HADS") {
#     if (missing_thr == TRUE) x <- x_HADS_w_missing_thr
#     else                     x <- x_HADS
#   }
#   ##
#   n_thr <- ncol(x[[1]]) - 1
#   ##
#   n_studies <- nrow(x[[1]])
#   n_studies
#   ##
#   # x <- arrange_missing_values(x)
# }
# # 




dput(true_DGM_Se)
dput(true_DGM_Sp)



# 
# 
# 
# start_sim_index = start_sim_index
# n_sims = 3
# model_parameterisation = model_parameterisation
# softplus = softplus
# box_cox = box_cox
# cts = cts
# random_thresholds = random_thresholds
# Dirichlet_random_effects_type = Dirichlet_random_effects_type
# local_pkg_dir =local_pkg_dir
# 
# 
# index_test_MA
# 
# 
# 
# 
# 
# # Run simulation
# output <- run_single_simulation(start_sim_index = start_sim_index,
#                                 n_sims = n_sims,
#                                 model_parameterisation = model_parameterisation,
#                                 softplus = softplus,
#                                 box_cox = box_cox,
#                                 cts = cts,
#                                 random_thresholds = random_thresholds,
#                                 Dirichlet_random_effects_type = Dirichlet_random_effects_type,
#                                 local_pkg_dir =local_pkg_dir,
#                                 ##
#                                 sim_config = sim_config)
# 
# 
# output
# 
# 











dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)


dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline"))) %>% print(n = 1000)



quantile(exp(rnorm(n = 10000, mean = log(2*100), sd = 1)), c(0.025, 0.50, 0.975))
quantile(exp(rnorm(n = 10000, mean = log(100), sd = 1)), c(0.025, 0.50, 0.975))

# Generate samples from the prior
n_samples <- 10000

# Sample log(kappa) from student_t(df=3, location=log(50), scale=1)
log_kappa_samples <- rt(n_samples, df = 5) * 1 + log(200)
log_kappa_samples <- rt(n_samples, df = 5) * 1 + log(50)
# Note: rt() generates standard t-distribution, so we scale and shift

# Transform to get kappa
kappa_samples <- exp(log_kappa_samples)

# Quick summary
summary(kappa_samples)
quantile(kappa_samples, probs = c(0.025, 0.5, 0.975))

# readRDS("simulation_results/test_GAD_7N_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flatscale_het_1het_nd_0.488het_d_0.669cut_SD_1_complete_results_list.RDS")
# 
# readRDS("simulation_results/test_GAD_7N_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_flatscale_het_1het_nd_0.488het_d_0.669cut_SD_1.RDS")
# 
# 



index_test_MA














require(stringr)

 
#  
param <- "Jones"
param <- "Xu"
# param <- "R&G"
##
# n_studies <- 10 ; index_test_MA <- "GAD_2" ; model_parameterisation <- param ; random_thresholds <- FALSE ##  
n_studies <- 10 ; index_test_MA <- "GAD_2" ; model_parameterisation <- param ; random_thresholds <- TRUE  ##   
##
# n_studies <- 50 ; index_test_MA <- "GAD_2" ; model_parameterisation <- param ; random_thresholds <- FALSE ##   
n_studies <- 50 ; index_test_MA <- "GAD_2" ; model_parameterisation <- param ; random_thresholds <- TRUE  ##    
##
##
##
# n_studies <- 10 ; index_test_MA <- "HADS" ; model_parameterisation <- param ; random_thresholds <- FALSE ##   
n_studies <- 10 ; index_test_MA <- "HADS" ; model_parameterisation <- param ; random_thresholds <- TRUE  ##   
##
# n_studies <- 50 ; index_test_MA <- "HADS" ; model_parameterisation <- param ; random_thresholds <- FALSE ##  
n_studies <- 50 ; index_test_MA <- "HADS" ; model_parameterisation <- param ; random_thresholds <- TRUE  ##  
##
##
##
# n_studies <- 10 ; index_test_MA <- "BAI" ; model_parameterisation <- param ; random_thresholds <- FALSE ##   
n_studies <- 10 ; index_test_MA <- "BAI" ; model_parameterisation <- param ; random_thresholds <- TRUE  ##   
##
# n_studies <- 50 ; index_test_MA <- "BAI" ; model_parameterisation <- param ; random_thresholds <- FALSE ##  
n_studies <- 50 ; index_test_MA <- "BAI" ; model_parameterisation <- param ; random_thresholds <- TRUE  #
##
#
# n_studies <- 10 ; index_test_MA <- "GAD_7" ; model_parameterisation <- param ; random_thresholds <- FALSE
# n_studies <- 10 ; index_test_MA <- "GAD_7" ; model_parameterisation <- param ; random_thresholds <- TRUE #'###
# #
# n_studies <- 50 ; index_test_MA <- "GAD_7" ; model_parameterisation <- param ; random_thresholds <- FALSE
# n_studies <- 50 ; index_test_MA <- "GAD_7" ; model_parameterisation <- param ; random_thresholds <- TRUE
##
##
##
##
##
# alpha_pi <- "flat"
alpha_pi <- "one"
##
##
file_name_base <- paste0(    "test_", index_test_MA,
                        "N_", 500,
                        "n_studies_", n_studies,
                        ##
                        "param_", model_parameterisation,
                        "soft_", FALSE,
                        "rand_thr_", random_thresholds,
                        "alpha_pi_", alpha_pi)


file_name <- paste0(file_name_base, "_complete_results_list.RDS")

# 
{
    counter <- 0
    for (seed in 1:5000) {

      try({

        # ##
        # file_name_base <- paste0(    "test_", index_test_MA,
        #                              "N_", 500,
        #                              "n_studies_", n_studies,
        #                              ##
        #                              "param_", model_parameterisation,
        #                              "soft_", FALSE,
        #                              "rand_thr_", 1,
        #                              "alpha_pi_", alpha_pi)


        file_name <- paste0(file_name_base, "_complete_results_list.RDS")

        file_name_specific_run <- paste0("seed_", seed, "_",
                                         file_name_base,
                                         "_results_list.RDS")
        ##
        results_list <- readRDS(file.path("simulation_results", file_name_specific_run))
        counter <- counter + 1

        # file_name_base <- paste0(    "test_", index_test_MA,
        #                              "N_", 500,
        #                              "n_studies_", n_studies,
        #                              ##
        #                              "param_", model_parameterisation,
        #                              "soft_", FALSE,
        #                              "rand_thr_", TRUE,
        #                              "alpha_pi_", alpha_pi)
        #
        # file_name_specific_run <- paste0("seed_", seed, "_",
        #                                  file_name_base,
        #                                  "_results_list.RDS")
        # ##
        # # results_list <- readRDS(file.path("simulation_results", file_name_specific_run))
        #
        # saveRDS(object = results_list, file = file.path(getwd(), "simulation_results", file_name_specific_run))



      })

    }

    print(counter)
}
# 


complete_results_list <- readRDS(file.path(getwd(), "simulation_results", file_name))
# complete_results_list
length(complete_results_list)

# complete_results_list[[length(complete_results_list)]]$true_DGM_Se
# complete_results_list[[length(complete_results_list)]][[7]] %>% print(n=1000)


random_thresholds

min_ESS_vec <- c()
for (n in 1:length(complete_results_list)) {
    tibble_Se_Sp <- MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE, 
                                                  tibble = complete_results_list[[n]][[7]],
                                                  param_strings_vec = c("Se_baseline", "Sp_baseline"),
                                                  condition = "exact_match") 
    
    min_ESS_vec[n] <- (min(tibble_Se_Sp$n_eff, na.rm = TRUE))
}
min_ESS_vec <- (min_ESS_vec[!is.infinite(min_ESS_vec)])
quantile(min_ESS_vec, probs = c(0.05, 0.10, 0.25, 0.50))
##
n_studies
index_test_MA
##
## ---- For GAD-2, n_studies = 10, random-C model:
##
# 5%    10%    25%    50% 
# 917.7 1010.0 1198.0 1450.0 
##
## ---- For GAD-2, n_studies = 50, random-C model:
##
# 5%     10%     25%     50% 
# 700.75  765.50  887.25 1030.50 
##
## ---- For HADS, n_studies = 10, random-C model:
##
# 5%    10%    25%    50% 
# 578.20 699.70 842.75 998.00 
##
## ---- For HADS, n_studies = 50, random-C model:
##
# 5%   10%   25%   50% 
# 303.2 383.2 432.5 490.5 
##
##
## ---- For BAI, n_studies = 10, random-C model:
##
# 5%    10%    25%    50% 
# 233.2  399.2  841.0 1206.0 
##
## ---- For BAI, n_studies = 50, random-C model:
##
# 5%    10%    25%    50% 
# 109.00 115.10 146.25 178.00 




##
compute_equivalent_ESS(N = 50, 
                       N_new = 10, 
                       ESS_observed = 183.616)
# Original: N = 50 studies, ESS = 184
# New:      N = 10 studies
# 
# Equivalent ESS for same MC error:
# Conservative (SD ∝ 1/N^0.33): ESS ≈ 531
# Typical      (SD ∝ 1/N^0.50): ESS ≈ 918
# Optimistic   (SD ∝ 1/N^1.00): ESS ≈ 4590
# 
# Range: 531 - 4590



# complete_results_list[[15]]

# csv_file <- read.csv(file.path("simulation_results", csv_name))
# csv_file


str(complete_results_list)
 
if (stringr::str_detect(string = file_name, pattern = "GAD_2")) { 
  
      index_test_MA <- "GAD_2"
      min_k <- 1
      max_k <- 6
      ##
      true_DGM_Se_GAD_2 <- c(95.563, 85.862, 63.12, 47.106, 27.594, 15.446)
      true_DGM_Sp_GAD_2 <- c(40.389, 67.085, 86.262, 92.329, 96.69, 98.577)
      ##
      true_DGM_Se <- true_DGM_Se_GAD_2
      true_DGM_Sp <- true_DGM_Sp_GAD_2
  
} else if (stringr::str_detect(string = file_name, pattern = "HADS")) { 
  
      index_test_MA <- "HADS"
      min_k <- 3
      max_k <- 17
      ##
      true_DGM_Se_HADS <- c(98.2, 97.834, 96.429, 95.211, 93.447, 90.584, 85.738, 80.568, 
                       73.458, 63.496, 54.196, 46.979, 38.097, 28.518, 20.962, 14.992, 
                       11.453, 8.952, 6.883, 4.428, 3.444)
      true_DGM_Sp_HADS <- c(7.145, 15.157, 24.563, 34.619, 44.948, 54.601, 64.196, 72.843, 
                       80.069, 85.533, 89.248, 92.425, 94.684, 96.402, 97.788, 98.56, 
                       99.17, 99.509, 99.687, 99.845, 99.927)
      ##
      true_DGM_Se <- true_DGM_Se_HADS
      true_DGM_Sp <- true_DGM_Sp_HADS
  
} else if (stringr::str_detect(string = file_name, pattern = "GAD_7")) {
  
      index_test_MA <- "GAD_7"
      min_k <- 3
      max_k <- 17
      ##
      true_DGM_Se_GAD_7 <-  c(97.301, 95.945, 94.527, 92.268, 89.168, 84.232, 79.27, 71.388, 
                        65.608, 59.321, 52.195, 44.815, 37.824, 31.994, 25.598, 21.024, 
                        17.119, 12.366, 8.378, 4.91, 2.771)
      true_DGM_Sp_GAD_7 <-  c(20.604, 34.948, 47.189, 57.574, 66.619, 73.613, 79.562, 85.355, 
                        88.37, 90.99, 92.876, 94.388, 95.502, 96.555, 97.475, 98.195, 
                        98.71, 98.997, 99.341, 99.578, 99.713)
      ##
      true_DGM_Se <- true_DGM_Se_GAD_7
      true_DGM_Sp <- true_DGM_Sp_GAD_7
  
} else if (stringr::str_detect(string = file_name, pattern = "BAI")) {
  
      index_test_MA <- "BAI"
      min_k <- 3
      max_k <- 43
      ##
      true_DGM_Se_BAI <-   c(98.121, 97.614, 96.772, 95.936, 93.841, 92.638, 91.416, 89.622, 
                         88.139, 86.268, 85.294, 82.32, 80.428, 78.427, 74.974, 70.799, 
                         68.462, 63.491, 60.883, 58.539, 56.383, 54.35, 50.759, 48.56, 
                         45.896, 44.695, 42.095, 39.667, 37.971, 37.413, 33.802, 31.229, 
                         29.892, 28.623, 27.3, 24.923, 24.46, 22.646, 20.613, 19.532, 
                         17.637, 15.099, 14.384, 13.499, 12.561, 11.41, 10.189, 8.997, 
                         8.138, NA, 6.857, 6.437, 5.878, 5.314, 4.552, 3.904, 3.548, 
                         NA, NA, 2.332, NA, 1.587, 0.86)
      true_DGM_Sp_BAI <-  c(11.205, 19.163, 26.865, 32.213, 37.886, 42.833, 48.774, 53.075, 
                        58.456, 61.69, 64.889, 68.052, 71.944, 74.159, 76.549, 79.443, 
                        81.474, 83.406, 84.623, 85.742, 86.979, 87.926, 89.172, 90.354, 
                        91.239, 92.089, 92.653, 93.142, 93.623, 93.994, 94.74, 95.086, 
                        95.304, 95.665, 95.942, 96.469, 96.753, 97.327, 97.462, 97.594, 
                        97.737, 98.028, 98.057, 98.35, 98.671, 98.763, 99.002, 99.136, 
                        99.132, NA, 99.304, 99.341, 99.498, 99.506, 99.614, 99.639, 
                        99.692, NA, NA, 99.822, NA, 99.925, 99.945)
      ##
      true_DGM_Se <- true_DGM_Se_BAI
      true_DGM_Sp <- true_DGM_Sp_BAI
      
      
  
}


index_test_MA



complete_results_list[[1000]]$

# Usage example:
metrics <- compute_simulation_metrics(complete_results_list = complete_results_list,
                                      true_DGM_Se = true_DGM_Se,
                                      true_DGM_Sp = true_DGM_Sp,
                                      min_k = min_k,
                                      max_k = max_k)

metrics$overall


complete_results_list[[1]]

{
  
        {
            cat(paste("\n--------------------------------------------------\n"))
            message("if target_MCSE_RMSE = 0.15%:")
            print(calc_min_sim_size_RMSE(
              metrics = metrics,
              target_MCSE_RMSE = 0.15, 
              min_k = min_k,
              max_k = max_k
            ))
            cat(paste("\n--------------------------------------------------\n"))
        }
        
        {
          cat(paste("\n--------------------------------------------------\n"))
          message("if target_MCSE_RMSE = 0.125%:")
          print(calc_min_sim_size_RMSE(
            metrics = metrics,
            target_MCSE_RMSE = 0.125, 
            min_k = min_k,
            max_k = max_k
          ))
          cat(paste("\n--------------------------------------------------\n"))
        }
        
        length(complete_results_list)
        
        
        
        # Create nice tables
        tables <- create_metrics_summary_table(metrics = metrics,
                                               min_k = min_k,
                                               max_k = max_k)
        
        print(tables$overall)
        ##
        cat(paste("\n--------------------------------------------------\n"))
        print(paste("n_sims = ", length(complete_results_list)))
        print(paste("index_test_MA = ", index_test_MA))
        print(paste("n_studies = ", n_studies))
        print(paste("random_thresholds = ", random_thresholds))
        print(paste("model_parameterisation = ", model_parameterisation))
        ##
        cat(paste("\n--------------------------------------------------\n"))
  
}








# Now you can use calc_min_sim_size

 




dput(sim_results$true_DGM_Se[[4]]*100)
dput(sim_results$true_DGM_Sp[[4]]*100)


complete_results_list[[1]][[7]] %>% print(n = 100)
 

{
  
   print(paste("model_parameterisation = ", complete_results_list[[1]]$model_parameterisation))
   print(paste("random_thresholds = ", complete_results_list[[1]]$random_thresholds))
  
  {
    # Print results
    cat("\n========== SIMULATION METRICS SUMMARY ==========\n")
    cat(sprintf("Number of simulations: %d\n", metrics$n_sims))
    cat(sprintf("Thresholds analyzed: %d to %d\n\n", vec_index_inner[1], tail(vec_index_inner, 1)))
    
    cat("--- Overall Performance (averaged across thresholds) ---\n")
    print(tables$overall, row.names = FALSE)
    
    cat("\n--- Performance by Threshold ---\n")
    print(tables$by_threshold, row.names = FALSE)
  }
  
  
   

   
  

# {
#   
#   # cat(paste("\n--------------------------------------------------\n"))
#   # message("if min_MCSE_pct = 0.5%:")
#   # print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec,
#   #                         min_MCSE_pct = 0.50,
#   #                         min_k = vec_index_inner[1],
#   #                         max_k = tail(vec_index_inner, 1)))
#   # cat(paste("\n--------------------------------------------------\n"))
#   cat(paste("\n--------------------------------------------------\n"))
#   message("if min_MCSE_pct = 0.25%:")
#   print(calc_min_sim_size(
#     SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
#     min_MCSE_pct = 0.25, 
#     min_k = vec_index_inner[1],
#     max_k = tail(vec_index_inner, 1)
#   ))
#   
#   cat(paste("\n--------------------------------------------------\n"))
#   ##
#   ##
#   ##
#   cat(paste("\n--------------------------------------------------\n"))
#   message("if min_MCSE_pct = 0.125%:")
#   print(calc_min_sim_size(
#     SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
#     min_MCSE_pct = 0.125, 
#     min_k = vec_index_inner[1],
#     max_k = tail(vec_index_inner, 1)
#   ))
#   cat(paste("\n--------------------------------------------------\n"))
#   ##
#   ##
#   ##
#   cat(paste("\n--------------------------------------------------\n"))
#   message("if min_MCSE_pct = 0.1%:")
#   print(calc_min_sim_size(
#     SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
#     min_MCSE_pct = 0.10, 
#     min_k = vec_index_inner[1],
#     max_k = tail(vec_index_inner, 1)
#   ))
#   cat(paste("\n--------------------------------------------------\n"))
#   ##
#   ##
#   ##
#   cat(paste("\n--------------------------------------------------\n"))
#   message("if min_MCSE_pct = 0.05%:")
#   print(calc_min_sim_size(
#     SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
#     min_MCSE_pct = 0.05, 
#     min_k = vec_index_inner[1],
#     max_k = tail(vec_index_inner, 1)
#   ))
#   
#   cat(paste("\n--------------------------------------------------\n"))
# }


# ##
{
  # cat(paste("\n--------------------------------------------------\n"))
  # message("if target_MCSE_RMSE = 0.25%:")
  # print(calc_min_sim_size_RMSE(
  #   metrics = metrics,
  #   target_MCSE_RMSE = 0.25, 
  #   min_k = vec_index_inner[1],
  #   max_k = tail(vec_index_inner, 1)
  # ))
  # 
  # cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.125%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.125, 
    min_k = min_k,
    max_k = max_k
  ))
  cat(paste("\n--------------------------------------------------\n"))
  ##
  ##
  ##
  # cat(paste("\n--------------------------------------------------\n"))
  # message("if target_MCSE_RMSE = 0.1%:")
  # print(calc_min_sim_size_RMSE(
  #   metrics = metrics,
  #   target_MCSE_RMSE = 0.10, 
  #   min_k = vec_index_inner[1],
  #   max_k = tail(vec_index_inner, 1)
  # ))
  # cat(paste("\n--------------------------------------------------\n"))
  # ##
  # ##
  # ##
  cat(paste("\n--------------------------------------------------\n"))
  message("if target_MCSE_RMSE = 0.05%:")
  print(calc_min_sim_size_RMSE(
    metrics = metrics,
    target_MCSE_RMSE = 0.05,
    min_k = min_k,
    max_k = max_k)
  )
  # 
  # cat(paste("\n--------------------------------------------------\n"))
}


}

#






 


# 
# 
# library(tidyverse)
# library(patchwork)
# 
# 
# 
# 
# # Usage example:
# metrics <- compute_simulation_metrics(complete_results_list = complete_results_list,
#                                       true_DGM_Se = complete_results_list[[1]]$true_DGM_Se, 
#                                       true_DGM_Sp = complete_results_list[[1]]$true_DGM_Sp,
#                                       min_k = 3,
#                                       max_k = 17)
# 
# 


{
##
## ----------- For GAD-2 ------------------------------------------------------------------------------------------- :
##
# ##
# ## ---- Bivariate:
##
file_name <- "test_GAD_2N_500n_studies_50param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_GAD_2N_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
# ---- Xu, Random C:
#
file_name <- "test_GAD_2N_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
# ---- R&G, Random C:
#
file_name <- "test_GAD_2N_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_50studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
# ##
## ----------- For HADS ---------------------------------------------------------------------------------------------- :
##
## ---- Bivariate:
#
file_name <- "test_HADSN_500n_studies_50param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_HADSN_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_HADSN_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_HADSN_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Random C:
##
file_name <- "test_HADSN_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
## ---- R&G, Random C:
##
file_name <- "test_HADSN_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_50studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ----------- For BAI ---------------------------------------------------------------------------------------------- :
##
##
## ---- Bivariate:
##
file_name <- "test_BAIN_500n_studies_50param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_BAIN_500n_studies_50param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_BAIN_500n_studies_50param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_BAIN_500n_studies_50param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Random C:
##
file_name <- "test_BAIN_500n_studies_50param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Random C:
##
file_name <- "test_BAIN_500n_studies_50param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_50studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
##
## ------------------------------ n_studies == 10:-------------------------------------------------------------------------------------
##
##
## ----------- For GAD-2 ------------------------------------------------------------------------------------------- :
##
##
## ---- Bivariate:
##
file_name <- "test_GAD_2N_500n_studies_10param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_GAD_2N_500n_studies_10param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_10param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_GAD_2N_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
# ---- Xu, Random C:
#
file_name <- "test_GAD_2N_500n_studies_10param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
# ---- R&G, Random C:
#
file_name <- "test_GAD_2N_500n_studies_10param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_GAD_2_10studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
# ##
## ----------- For HADS ---------------------------------------------------------------------------------------------- :
##
## ---- Bivariate:
#
file_name <- "test_HADSN_500n_studies_10param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_HADSN_500n_studies_10param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_HADSN_500n_studies_10param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_HADSN_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Random C:
##
file_name <- "test_HADSN_500n_studies_10param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
#
## ---- R&G, Random C:
##
file_name <- "test_HADSN_500n_studies_10param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_HADS_10studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ----------- For BAI ---------------------------------------------------------------------------------------------- :
##
##
## ---- Bivariate:
##
file_name <- "test_BAIN_500n_studies_10param_Bivariatesoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_Bivariate <- readRDS(file.path("simulation_results", file_name))
##
## ---- Jones:
##
file_name <- "test_BAIN_500n_studies_10param_Jonessoft_FALSErand_thr_FALSEalpha_pi_flat_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_Jones <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Fixed:
##
file_name <- "test_BAIN_500n_studies_10param_Xusoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_Xu_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Fixed:
##
file_name <- "test_BAIN_500n_studies_10param_R&Gsoft_FALSErand_thr_FALSEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_RG_fixed_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- Xu, Random C:
##
file_name <- "test_BAIN_500n_studies_10param_Xusoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_Xu_random_alpha2 <- readRDS(file.path("simulation_results", file_name))
##
## ---- R&G, Random C:
##
file_name <- "test_BAIN_500n_studies_10param_R&Gsoft_FALSErand_thr_TRUEalpha_pi_one_complete_results_list.RDS"
complete_results_list_BAI_10studies_500N_RG_random_alpha2 <- readRDS(file.path("simulation_results", file_name))

}


# ##
# ## ---- GAD-2:
# ## 
# model_parameterisation <- "Jones"  ; alpha_pi <- "flat" ;  n_studies <- 10 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"     ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"     ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "GAD_2" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "GAD_2" ; random_thresholds <- TRUE
# ##
# model_parameterisation <- "Jones" ; alpha_pi <- "flat" ;  n_studies <- 50 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "GAD_2" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "GAD_2" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "GAD_2" ; random_thresholds <- TRUE
# ## 
# ## ---- HADS:
# ##
# model_parameterisation <- "Jones" ; alpha_pi <- "flat" ;  n_studies <- 10 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "HADS" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "HADS" ; random_thresholds <- TRUE
# ##
# model_parameterisation <- "Jones" ; alpha_pi <- "flat" ;  n_studies <- 50 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "HADS" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "HADS" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "HADS" ; random_thresholds <- TRUE
## 
## ---- BAI:
##
# model_parameterisation <- "Bivariate" ; alpha_pi <- "flat" ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "Jones" ; alpha_pi <- "flat" ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- FALSE
model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 10 ; test <- "BAI" ; random_thresholds <- TRUE
##
model_parameterisation <- "Bivariate" ; alpha_pi <- "flat" ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "Jones" ; alpha_pi <- "flat" ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- FALSE
model_parameterisation <- "Xu"    ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- TRUE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- FALSE
# model_parameterisation <- "R&G"   ; alpha_pi <- "one"  ;  n_studies <- 50 ; test <- "BAI" ; random_thresholds <- TRUE
 

# library(foreach)
# library(doParallel)  # or doMC on Linux/Mac
# 
# # Set up parallel backend (adjust cores as needed)
# n_cores <- 48
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# 
# # Run the foreach loop
# datasets_sim_list <- foreach(seed = 1:8000,
#                              .packages = c("MetaOrdDTA"),  # Add any required packages here
#                              .errorhandling = 'pass') %dopar% {
# 
#                                set.seed(seed)
# 
#                                require(MetaOrdDTA)
# 
#                                n_studies <- 50
#                                N_per_study_mean <- 500
#                                N_per_study_SD <- N_per_study_mean
#                                ##
#                                true_Mean_prev <- 0.25 ;  true_SD_probit_prev <- 0.25
#                                ##
#                                index_test_chosen_index <- 3
# 
#                                # Source the files
#                                source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA.R"), local = TRUE)
#                                source(file.path(local_pkg_dir, "inst", "examples", "NMA_missing_thr_prep_data.R"), local = TRUE)
#                                source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"), local = TRUE)
# 
#                                # Return a list with all the data for this seed
#                                list(
#                                  x_GAD2_w_missing_thr = x_GAD2_w_missing_thr,
#                                  x_GAD2 = x_GAD2,
#                                  x_HADS_w_missing_thr = x_HADS_w_missing_thr,
#                                  x_HADS = x_HADS,
#                                  x_BAI_w_missing_thr = x_BAI_w_missing_thr,
#                                  x_BAI = x_BAI,
#                                  x_GAD7_w_missing_thr = x_GAD7_w_missing_thr,
#                                  x_GAD7 = x_GAD7
#                                )
#                              }
# 
# # Don't forget to stop the cluster
# stopCluster(cl)
# 
# 
# str(datasets_sim_list)
# length(datasets_sim_list)
# 
# datasets_sim_list[[1]]
# 
# file_name_string <- paste0( "sim_data_list_",
#                             "N_", 500, "_",
#                             "n_studies_", n_studies)
# ##
# data_list_file <- file.path(output_dir,
#                                paste0(file_name_string, ".RDS"))
# ##
# saveRDS(datasets_sim_list, data_list_file)
# 



{
      ##
      true_DGM_Se_GAD_2 <- c(95.563, 85.862, 63.12, 47.106, 27.594, 15.446)
      true_DGM_Sp_GAD_2 <- c(40.389, 67.085, 86.262, 92.329, 96.69, 98.577)
      ##
      true_DGM_Se <- true_DGM_Se_GAD_2
      true_DGM_Sp <- true_DGM_Sp_GAD_2
     
      ##
      true_DGM_Se_HADS <- c(98.2, 97.834, 96.429, 95.211, 93.447, 90.584, 85.738, 80.568, 
                            73.458, 63.496, 54.196, 46.979, 38.097, 28.518, 20.962, 14.992, 
                            11.453, 8.952, 6.883, 4.428, 3.444)
      true_DGM_Sp_HADS <- c(7.145, 15.157, 24.563, 34.619, 44.948, 54.601, 64.196, 72.843, 
                            80.069, 85.533, 89.248, 92.425, 94.684, 96.402, 97.788, 98.56, 
                            99.17, 99.509, 99.687, 99.845, 99.927)
      ##
      true_DGM_Se <- true_DGM_Se_HADS
      true_DGM_Sp <- true_DGM_Sp_HADS
     
     
      true_DGM_Se_BAI <-   c(98.121, 97.614, 96.772, 95.936, 93.841, 92.638, 91.416, 89.622, 
                             88.139, 86.268, 85.294, 82.32, 80.428, 78.427, 74.974, 70.799, 
                             68.462, 63.491, 60.883, 58.539, 56.383, 54.35, 50.759, 48.56, 
                             45.896, 44.695, 42.095, 39.667, 37.971, 37.413, 33.802, 31.229, 
                             29.892, 28.623, 27.3, 24.923, 24.46, 22.646, 20.613, 19.532, 
                             17.637, 15.099, 14.384, 13.499, 12.561, 11.41, 10.189, 8.997, 
                             8.138, NA, 6.857, 6.437, 5.878, 5.314, 4.552, 3.904, 3.548, 
                             NA, NA, 2.332, NA, 1.587, 0.86)
      true_DGM_Sp_BAI <-  c(11.205, 19.163, 26.865, 32.213, 37.886, 42.833, 48.774, 53.075, 
                            58.456, 61.69, 64.889, 68.052, 71.944, 74.159, 76.549, 79.443, 
                            81.474, 83.406, 84.623, 85.742, 86.979, 87.926, 89.172, 90.354, 
                            91.239, 92.089, 92.653, 93.142, 93.623, 93.994, 94.74, 95.086, 
                            95.304, 95.665, 95.942, 96.469, 96.753, 97.327, 97.462, 97.594, 
                            97.737, 98.028, 98.057, 98.35, 98.671, 98.763, 99.002, 99.136, 
                            99.132, NA, 99.304, 99.341, 99.498, 99.506, 99.614, 99.639, 
                            99.692, NA, NA, 99.822, NA, 99.925, 99.945)
      ##
      true_DGM_Se <- true_DGM_Se_BAI
      true_DGM_Sp <- true_DGM_Sp_BAI
}





##
total_target_sims <- 7000
##


file_name_string <- paste0( "test_", test,
                            "N_", 500,
                            "n_studies_", n_studies,
                            ##
                            "param_", model_parameterisation,
                            "soft_", FALSE,
                            "rand_thr_", random_thresholds,
                            "alpha_pi_", alpha_pi)

# 
# Also rebuild the complete results_list from individual seed files
complete_results_list <- vector("list", total_target_sims)
##
for (seed_idx in 1:total_target_sims) {
  try({  
    # Use the same filename pattern as when you save individual seeds
    seed_file <- file.path(output_dir,
                           paste0("seed_", seed_idx, "_", file_name_string, "_results_list.RDS"))
    ##
    if (file.exists(seed_file)) {
      
      complete_results_list[[seed_idx]] <- readRDS(seed_file)
      
      # Set seed for reproducibility
      seed <- seed_idx
      set.seed(seed)
      
      
      
    } else {
      cat(sprintf("Warning: seed file not found: %s\n", seed_file))
    }
  })
}

complete_results_list <- Filter(Negate(is.null), complete_results_list)

length(complete_results_list)

# complete_results_list[[100]]
# 
# complete_results_list[[1000]]


# Save the complete results list
try({
  results_list_file <- file.path(output_dir, 
                                 paste0(file_name_string, "_complete_results_list.RDS"))
  saveRDS(complete_results_list, results_list_file)
  cat(sprintf("Saved complete results list with %d non-NULL entries to: %s\n",
              sum(!sapply(complete_results_list, is.null)), results_list_file))
})






##
## ---- Make list of lists for complete results:
##
complete_results_list_GAD_2 <- list(  complete_results_list_GAD_2_10studies_500N_Bivariate,
                                      complete_results_list_GAD_2_10studies_500N_Jones,
                                      complete_results_list_GAD_2_10studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_GAD_2_10studies_500N_RG_fixed_alpha2,
                                      complete_results_list_GAD_2_10studies_500N_Xu_random_alpha2,
                                      complete_results_list_GAD_2_10studies_500N_RG_random_alpha2,
                                      ##
                                      complete_results_list_GAD_2_50studies_500N_Bivariate,
                                      complete_results_list_GAD_2_50studies_500N_Jones,
                                      complete_results_list_GAD_2_50studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_GAD_2_50studies_500N_RG_fixed_alpha2,
                                      complete_results_list_GAD_2_50studies_500N_Xu_random_alpha2,
                                      complete_results_list_GAD_2_50studies_500N_RG_random_alpha2)
##
##
##
complete_results_list_HADS <-  list(  complete_results_list_HADS_10studies_500N_Bivariate,
                                      complete_results_list_HADS_10studies_500N_Jones,
                                      complete_results_list_HADS_10studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_HADS_10studies_500N_RG_fixed_alpha2,
                                      complete_results_list_HADS_10studies_500N_Xu_random_alpha2,
                                      complete_results_list_HADS_10studies_500N_RG_random_alpha2,
                                      ##
                                      complete_results_list_HADS_50studies_500N_Bivariate,
                                      complete_results_list_HADS_50studies_500N_Jones,
                                      complete_results_list_HADS_50studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_HADS_50studies_500N_RG_fixed_alpha2,
                                      complete_results_list_HADS_50studies_500N_Xu_random_alpha2,
                                      complete_results_list_HADS_50studies_500N_RG_random_alpha2)
##
##
##
complete_results_list_BAI <- list(    complete_results_list_BAI_10studies_500N_Bivariate,
                                      complete_results_list_BAI_10studies_500N_Jones,
                                      complete_results_list_BAI_10studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_BAI_10studies_500N_RG_fixed_alpha2,
                                      complete_results_list_BAI_10studies_500N_Xu_random_alpha2,
                                      complete_results_list_BAI_10studies_500N_RG_random_alpha2,
                                      ##
                                      complete_results_list_BAI_50studies_500N_Bivariate,
                                      complete_results_list_BAI_50studies_500N_Jones,
                                      complete_results_list_BAI_50studies_500N_Xu_fixed_alpha2,
                                      complete_results_list_BAI_50studies_500N_RG_fixed_alpha2,
                                      complete_results_list_BAI_50studies_500N_Xu_random_alpha2,
                                      complete_results_list_BAI_50studies_500N_RG_random_alpha2)
##
{
      complete_results_list <- c(complete_results_list_GAD_2, 
                                 complete_results_list_HADS, 
                                 complete_results_list_BAI)
      
                                            
      n_DGMs <- 3
      n_models <- 6
      n_total <- 2*n_DGMs*n_models
      ##
      DGM_vec <- c(rep("GAD-2", length(complete_results_list_GAD_2)),
                   rep("HADS", length(complete_results_list_HADS)),
                   rep("BAI", length(complete_results_list_BAI)))
      length(DGM_vec)
      ##
      model_name_vec <- c("bivariate", 
                          "Jones", 
                          "O-bivariate, fixed C", 
                          "O-HSROC, fixed C", 
                          "O-bivariate, rand C", 
                          "O-HSROC, rand C")
      model_name_vec <- rep(model_name_vec, 2)
      model_name_vec <- rep(model_name_vec, n_DGMs)
      length(model_name_vec)
      ##
      n_studies_vec <- c(rep(10, n_models), rep(50, n_models))
      n_studies_vec <- rep(n_studies_vec, n_DGMs)
      length(n_studies_vec)
      ##
      N_vec <- rep(500, n_total)
      length(N_vec)
      ## 
      min_k_vec_GAD_2  <- rep(1, 2*n_models)
      min_k_vec_HADS   <- rep(3, 2*n_models)
      min_k_vec_GAD_7  <- rep(3, 2*n_models)
      min_k_vec <- c(min_k_vec_GAD_2, 
                     min_k_vec_HADS,
                     min_k_vec_GAD_7)
      length(min_k_vec)
      ##
      max_k_vec_GAD_2  <- rep(6, 2*n_models)
      max_k_vec_HADS   <- rep(17, 2*n_models)
      max_k_vec_BAI    <- rep(43, 2*n_models)
      max_k_vec <- c(max_k_vec_GAD_2, 
                     max_k_vec_HADS,
                     max_k_vec_BAI)
      length(max_k_vec)
      ##
      true_DGM_Se_list <- true_DGM_Sp_list <- list()
      ##
      true_DGM_Se_list <- c( rep(list(true_DGM_Se_GAD_2), length(complete_results_list_GAD_2)),
                             rep(list(true_DGM_Se_HADS),  length(complete_results_list_HADS)),
                             rep(list(true_DGM_Se_BAI), length(complete_results_list_BAI)))
      ##
      true_DGM_Sp_list <- c( rep(list(true_DGM_Sp_GAD_2), length(complete_results_list_GAD_2)),
                             rep(list(true_DGM_Sp_HADS),  length(complete_results_list_HADS)),
                             rep(list(true_DGM_Sp_BAI), length(complete_results_list_BAI)))
}
     


          
# complete_results_list = complete_results_list[[1]]
# test = DGM_vec[1]
# true_DGM_Se = true_DGM_Se_list[[1]]
# true_DGM_Sp = true_DGM_Sp_list[[1]]
# model_name = model_name_vec[1]
# n_studies = n_studies_vec[1]
# N = N_vec[1]
# min_k = min_k_vec[1]
# max_k = max_k_vec[1]

require(purrr)
require(dplyr)
require(ggplot2)
require(tidyr)
##
first_row <- convert_metrics_to_tibble(  complete_results_list = complete_results_list[[1]], 
                                         output_dir = output_dir,
                            test = DGM_vec[1],
                            true_DGM_Se = true_DGM_Se_list[[1]],
                            true_DGM_Sp = true_DGM_Sp_list[[1]],
                            model_name = model_name_vec[1],
                            n_studies = n_studies_vec[1],
                            N = N_vec[1],
                            min_k = min_k_vec[1], 
                            max_k = max_k_vec[1])
overall_sim_data <- first_row

##
first_rows_thr <- get_threshold_data(  complete_results_list = complete_results_list[[1]], 
                                         output_dir = output_dir,
                                         test = DGM_vec[1],
                                         true_DGM_Se = true_DGM_Se_list[[1]],
                                         true_DGM_Sp = true_DGM_Sp_list[[1]],
                                         model_name = model_name_vec[1],
                                         n_studies = n_studies_vec[1],
                                         N = N_vec[1],
                                         min_k = min_k_vec[1], 
                                         max_k = max_k_vec[1])
threshold_sim_data <- first_rows_thr


for (i in 2:n_total) {
    row_i <-  convert_metrics_to_tibble(  complete_results_list = complete_results_list[[i]], 
                                          output_dir = output_dir,
                                          test = DGM_vec[i],
                                          true_DGM_Se = true_DGM_Se_list[[i]],
                                          true_DGM_Sp = true_DGM_Sp_list[[i]],
                                          model_name = model_name_vec[i],
                                          n_studies = n_studies_vec[i],
                                          N = N_vec[i],
                                          min_k = min_k_vec[i], 
                                          max_k = max_k_vec[i])
    ##
    overall_sim_data <- bind_rows(overall_sim_data, row_i)
    ##
    ##
    ##
    rows_i_thr <-  get_threshold_data(  complete_results_list = complete_results_list[[i]], 
                                          output_dir = output_dir,
                                          test = DGM_vec[i],
                                          true_DGM_Se = true_DGM_Se_list[[i]],
                                          true_DGM_Sp = true_DGM_Sp_list[[i]],
                                          model_name = model_name_vec[i],
                                          n_studies = n_studies_vec[i],
                                          N = N_vec[i],
                                          min_k = min_k_vec[i], 
                                          max_k = max_k_vec[i])
    ##
    threshold_sim_data <- bind_rows(threshold_sim_data, rows_i_thr)
}


overall_sim_data
threshold_sim_data


{
  print(paste("------------------------------- GAD-2 ------------------------------"))
  print(signif(mean(true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(sd(true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(max(true_mean_GAD_2_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  ##
  print(signif(mean(true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(sd(true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(max(true_mean_GAD_2_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  
  
  print(paste("------------------------------- HADS ------------------------------"))
  print(signif(mean(true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(sd(true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(max(true_mean_HADS_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  ##
  print(signif(mean(true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(sd(true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(max(true_mean_HADS_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  
  
  print(paste("------------------------------- BAI ------------------------------"))
  print(signif(mean(true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(sd(true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  print(signif(max(true_mean_BAI_params_list$Ord_bivariate$C_nd_SD_prob_scale), 3))
  ##
  print(signif(mean(true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(sd(true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
  print(signif(max(true_mean_BAI_params_list$Ord_bivariate$C_d_SD_prob_scale), 3))
}




options(pillar.sigfig = 4)

GAD_2_data_10_studies <- overall_sim_data %>% filter(DGM == "GAD-2", n_studies == 10) %>% select(-n_sims_0.125, -n_sims_0.10) %>% print(n = 100)
GAD_2_data_50_studies <- overall_sim_data %>% filter(DGM == "GAD-2", n_studies == 50) %>% select(-n_sims_0.125, -n_sims_0.10) %>% print(n = 100)
##
HADS_data_10_studies <- overall_sim_data %>% filter(DGM == "HADS", n_studies == 10)  %>% select(-n_sims_0.125, -n_sims_0.10) %>% print(n = 100)
HADS_data_50_studies <- overall_sim_data %>% filter(DGM == "HADS", n_studies == 50)  %>% select(-n_sims_0.125, -n_sims_0.10) %>% print(n = 100)
##
# overall_sim_data %>% filter(DGM == "HADS") %>% print(n = 100)
# overall_sim_data %>% filter(DGM == "GAD-7") %>% print(n = 100)
BAI_data_10_studies <- overall_sim_data %>% filter(DGM == "BAI", n_studies == 10) %>% print(n = 100)
BAI_data_50_studies <- overall_sim_data %>% filter(DGM == "BAI", n_studies == 50) %>% print(n = 100)
##

overall_sim_data %>% print(n=100)


data_to_use <- BAI_data_50_studies
##
data_to_use$n_sims
##
round(data_to_use$avg_RMSE_Se, 3)
round(data_to_use$avg_RMSE_Sp, 3)
round(data_to_use$max_RMSE_Se, 3)
round(data_to_use$max_RMSE_Sp, 3)
##
round(data_to_use$avg_abs_bias_Se, 3)
round(data_to_use$avg_abs_bias_Sp, 3)
##
round(data_to_use$avg_coverage_Se*100, 1)
round(data_to_use$avg_coverage_Sp*100, 1)
##
round(data_to_use$avg_width_Se, 1)
round(data_to_use$avg_width_Sp, 1)
##
round(data_to_use$avg_MCSE_RMSE_Se, 4)
round(data_to_use$avg_MCSE_RMSE_Sp, 4)
round(data_to_use$max_MCSE_RMSE_Se, 4)
round(data_to_use$max_MCSE_RMSE_Sp, 4)
##
round(data_to_use$avg_MCSE_bias_Se, 4)
round(data_to_use$avg_MCSE_bias_Sp, 4)




 
# complete_results_list = complete_results_list[[9]]
# 
# complete_results_list$true_DGM_Se
# 
# metrics_list <- compute_simulation_metrics(complete_results_list = complete_results_list,
#                                            true_DGM_Se = complete_results_list[[1]]$true_DGM_Se, 
#                                            true_DGM_Sp = complete_results_list[[1]]$true_DGM_Sp,
#                                            min_k = min_k,
#                                            max_k = max_k)
# tibble(
#   DGM = DGM,
#   Model = model_name,
#   n_studies = n_studies,
#   N = N,
#   # Overall metrics:
#   avg_RMSE_Se = metrics_list$overall$avg_RMSE_Se,
#   avg_RMSE_Sp = metrics_list$overall$avg_RMSE_Sp,
#   avg_MCSE_RMSE_Se = metrics_list$overall$avg_MCSE_RMSE_Se,
#   avg_MCSE_RMSE_Sp = metrics_list$overall$avg_MCSE_RMSE_Sp,
#   avg_abs_bias_Se = metrics_list$overall$avg_abs_bias_Se,
#   avg_abs_bias_Sp = metrics_list$overall$avg_abs_bias_Sp,
#   avg_MCSE_bias_Se = metrics_list$overall$avg_MCSE_bias_Se,
#   avg_MCSE_bias_Sp = metrics_list$overall$avg_MCSE_bias_Sp,
#   avg_coverage_Se = metrics_list$overall$avg_coverage_prob_Se,
#   avg_coverage_Sp = metrics_list$overall$avg_coverage_prob_Sp,
#   avg_width_Se = metrics_list$overall$avg_interval_width_Se,
#   avg_width_Sp = metrics_list$overall$avg_interval_width_Sp
# )
# 



overall_sim_data


# # Assuming you have multiple models' results
# # Example: combine all model results
# all_results <- bind_rows(
#   convert_metrics_to_tibble(metrics_jones, "Jones", 50, 500),
#   convert_metrics_to_tibble(metrics_obiv_fc, "O-bivariate-FC", 50, 500),
#   convert_metrics_to_tibble(metrics_ohsroc_fc, "O-HSROC-FC", 50, 500),
#   convert_metrics_to_tibble(metrics_obiv_rc, "O-bivariate-RC", 50, 500),
#   convert_metrics_to_tibble(metrics_ohsroc_rc, "O-HSROC-RC", 50, 500)
# # )
# 
# 
# 
# 
# # Combine threshold data for all models
# threshold_data <- bind_rows(
#   get_threshold_data(metrics_jones, "Jones", true_DGM_Se, true_DGM_Sp),
#   get_threshold_data(metrics_obiv_fc, "O-bivariate-FC", true_DGM_Se, true_DGM_Sp),
#   # ... etc for all models
# )
# 

# threshold_sim_data %>% print(width = 1000000)
# 
# threshold_sim_data$true_Se
# threshold_sim_data$Sp_median
# 
# 
# threshold_sim_data %>%
#   select(DGM, Model, n_studies, threshold, true_Se, true_Sp, Se_median, Sp_median) %>%
#   head(20) %>%
#   print()







##
## ---- Create all plots for each DGM ---------------------------------------------------
##
##
require(patchwork)
## ---- Create plots for each DGM:
##
# Get unique DGMs
DGM_names <- unique(overall_sim_data$DGM)

# Create lists to store plots
sROC_plots <- list()
rmse_plots <- list()
bias_plots <- list()
coverage_plots <- list()
##
base_size <- 16

for (dgm in DGM_names) {
  
  
        ##
        ## RMSE plots:
        ##
        plot_sROC <- create_sROC_plot(   data = threshold_sim_data,
                                         DGM_name = dgm,
                                         base_size = base_size)
        ##
        sROC_plots[[dgm]] <- plot_sROC
        ##
        ## RMSE plots:
        ##
        plot_se <- create_RMSE_plot(data = overall_sim_data, 
                                    DGM_name = dgm, 
                                    metric_type = "Se",
                                    base_size = base_size)
        ##
        plot_sp <- create_RMSE_plot(data = overall_sim_data, 
                                    DGM_name = dgm, 
                                    metric_type = "Sp",
                                    base_size = base_size)
        ##
        rmse_plots[[dgm]] <- plot_se + plot_sp + 
          plot_layout(ncol = 2) +
          plot_annotation(title = paste("DGM:", dgm))
        ##
        ## Bias plots:
        ##
        bias_se <- create_bias_plot(data = overall_sim_data, 
                                    DGM_name = dgm, 
                                    metric_type = "Se",
                                    base_size = base_size)
        ##
        bias_sp <- create_bias_plot(data = overall_sim_data,
                                    DGM_name = dgm, 
                                    metric_type = "Sp",
                                    base_size = base_size)
        ##
        bias_plots[[dgm]] <- bias_se + bias_sp + 
          plot_layout(ncol = 2) +
          plot_annotation(title = paste("DGM:", dgm))
        ##
        ## Coverage plots:
        ##
        coverage_se <- create_coverage_plot(data = overall_sim_data, 
                                            DGM_name = dgm, 
                                            metric_type = "Se",
                                            base_size = base_size)
        ##
        coverage_sp <- create_coverage_plot(data = overall_sim_data, 
                                            DGM_name = dgm, 
                                            metric_type = "Sp",
                                            base_size = base_size)
        ##
        coverage_plots[[dgm]] <- coverage_se + coverage_sp + 
          plot_layout(ncol = 2) +
          plot_annotation(title = paste("DGM:", dgm))
}

 

##
## ---- Plots for GAD-2:
##
# sROC_plots[["GAD-2"]]
rmse_plots[["GAD-2"]]
bias_plots[["GAD-2"]]
coverage_plots[["GAD-2"]]
##
# ggsave("Sim_study_sROC_GAD_2.png", sROC_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_RMSE_GAD_2.png", rmse_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Bias_GAD_2.png", bias_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_GAD_2.png", coverage_plots[["GAD-2"]], width = 21, height = 9, dpi = 400)

##
## ---- Plots for HADS:
##
# sROC_plots[["HADS"]]
rmse_plots[["HADS"]]
bias_plots[["HADS"]]
coverage_plots[["HADS"]]
##
# ggsave("Sim_study_sROC_HADS.png", sROC_plots[["HADS"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_RMSE_HADS.png", rmse_plots[["HADS"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Bias_HADS.png", bias_plots[["HADS"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_HADS.png", coverage_plots[["HADS"]], width = 21, height = 9, dpi = 400)

##
## ---- Plots for BAI:j
##
# sROC_plots[["BAI"]]
rmse_plots[["BAI"]]
bias_plots[["BAI"]]
coverage_plots[["BAI"]]
##
# ggsave("Sim_study_sROC_BAI.png", sROC_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_RMSE_BAI.png", rmse_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Bias_BAI.png", bias_plots[["BAI"]], width = 21, height = 9, dpi = 400)
ggsave("Sim_study_Coverage_BAI.png", coverage_plots[["BAI"]], width = 21, height = 9, dpi = 400)



 
graphics.off()   # Close all devices
.rs.restartR()   # Restart R session in RStudio

dev.off()  # Close current device

X11()  # For Linux


 library(ggplot2)
library(patchwork)
# ## --- Plot 4: Coverage probability -------------------------------------------------------------------
# ##
# plot_coverage_Se <- overall_sim_data %>%
#   pivot_longer(cols = c(avg_coverage_Se, avg_coverage_Sp),
#                names_to = "metric",
#                values_to = "coverage") %>%
#   mutate(
#     metric = factor(metric,
#                     levels = c("avg_coverage_Se"),
#                     labels = c("Se"))
#   ) %>%
#   ggplot(aes(x = Model, y = coverage, color = Model)) +
#   geom_point(size = 4) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      limits = c(0.7, 1)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(y = "Coverage Probability",
#        title = "95% CI Coverage Across Models",
#        subtitle = "Red line = Nominal 95% coverage") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   facet_wrap(~ DGM + n_studies, scales = "free_y", nrow = 3)
# ##
# plot_coverage_Se
# ##
# plot_coverage_Sp <- overall_sim_data %>%
#   pivot_longer(cols = c(avg_coverage_Se, avg_coverage_Sp),
#                names_to = "metric",
#                values_to = "coverage") %>%
#   mutate(
#     metric = factor(metric,
#                     levels = c("avg_coverage_Sp"),
#                     labels = c("Sp"))
#   ) %>%
#   ggplot(aes(x = Model, y = coverage, color = Model)) +
#   geom_point(size = 4) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      limits = c(0.7, 1)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(y = "Coverage Probability",
#        title = "95% CI Coverage Across Models",
#        subtitle = "Red line = Nominal 95% coverage") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   facet_wrap(~ DGM + n_studies, scales = "free_y", nrow = 3)
# ##
# plot_coverage_Sp
# ##
# ##
# require(patchwork)
# plot_coverage_Se + plot_coverage_Sp
# 
# 
# 
# true_mean_GAD_2_params_list$Ord_bivariate$beta_SD <- c(0.341, 0.559)
# true_mean_GAD_2_params_list$Ord_bivariate$beta_corr <- 0.172
# ##
# true_mean_HADS_params_list$Ord_bivariate$beta_SD <- c(0.235,  0.392)
# true_mean_HADS_params_list$Ord_bivariate$beta_corr <- 0.302
# ##
# true_mean_BAI_params_list$Ord_bivariate$beta_SD <- c(0.257,  0.368)
# true_mean_BAI_params_list$Ord_bivariate$beta_corr <- 0.005574
# 

# ##
# ## --- Combine plots
# ##
# # Using patchwork to combine
# combined_plot <- (plot_rmse | plot_bias) / 
#   plot_coverage /
#   plot_accuracy +
#   plot_layout(heights = c(1, 1, 2))
# 
# # Save
# ggsave("simulation_results_combined.png", combined_plot, 
#        width = 12, height = 14, dpi = 300)
# 















