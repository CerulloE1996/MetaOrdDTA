


# ============================================================================================================
# MetaOrdDTA Simulation Study Configuration
# ============================================================================================================

# Define simulation parameters
get_simulation_config <- function() {
  
  list(
    # Sample sizes
    n_studies_vec = c(10, 50),
    N_per_study_mean = 500,
    
    # Model configurations
    model_parameterisations = c("Jones", "Xu", "Bivariate", "R&G")
    
    # # Prior configurations
    # alpha_pi_options = c("one", "flat", "custom"),
    # 
    # # Computational settings
    # n_chains = 4,
    # n_burnin = 500,
    # n_iter_base = 1000,
    # adapt_delta = 0.65,
    # max_treedepth = 10,
    # 
    # # Threshold settings by test
    # threshold_configs = list(
    #   GAD_2 = list(vec_index_inner = 1:6, n_thr = 6),
    #   HADS = list(vec_index_inner = 3:17, n_thr = 20),
    #   BAI = list(vec_index_inner = 3:43, n_thr = 63),
    #   GAD_7 = list(vec_index_inner = 1:10, n_thr = 10)
    # ),
    # 
    # # MCSE targets
    # target_MCSE_RMSE = 0.25/100,
    # min_ESS_threshold = 500,
    # min_runs_needed = 30
  )
  
}



# Create scenario matrix (similar to the LC-MVP approach)
create_scenario_matrix <- function(config) {
      
      expand.grid(
        index_test_MA = names(config$threshold_configs),
        n_studies = config$n_studies_vec,
        model_parameterisation = config$model_parameterisations,
        alpha_pi = config$alpha_pi_options,
        box_cox = c(FALSE, TRUE),
        softplus = c(FALSE),
        random_thresholds = c(FALSE, TRUE),
        stringsAsFactors = FALSE
      ) %>%
        # Filter out invalid combinations
        filter(
          !(model_parameterisation == "Jones" & !box_cox),
          !(model_parameterisation == "Jones" & random_thresholds),
          !(model_parameterisation == "Bivariate" & random_thresholds),
          !(random_thresholds & !(model_parameterisation %in% c("Xu", "R&G")))
        ) %>%
        mutate(
          scenario_id = row_number(),
          Dirichlet_random_effects_type = case_when(
            random_thresholds & model_parameterisation == "Xu" ~ "kappa",
            random_thresholds & model_parameterisation == "R&G" ~ "kappa",
            TRUE ~ "fixed"
          )
        )
  
}
















