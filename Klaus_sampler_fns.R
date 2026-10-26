##
## ---- Klaus_sampler_fns.R ---------------------------------------------------------------------------------------------------------
##
## ---- The NicoStan settings of the MetaOrdDTA sim study (run_MULTIPLE_MetaOrd_simulation.R lines 650-746, copied as they are;
##      auto_chunking = TRUE added), and the check of a NicoStan fit against a cmdstanr fit of the same model by MCSE z-scores.
##
##
## ---- NicoStan_args of the sim study ------------------------------------------------------------------------------------------------
##
fn_NicoStan_args_of_the_sim_study <- function( n_burnin,
                                               nuisance_params,
                                               interest_only_params ) {

      NicoStan_args <- list(
                                            ## (every $sample setting that changes sampling, written out)
                                            adapt_delta                          = 0.80,  ## the PS7 value 
                                            ###
                                            # metric_shape_main                    = "dense", ## --------
                                            metric_shape_main                    = "diag",    ## --------
                                            ##
                                            nuisance_params                      = nuisance_params,
                                            sample_nuisance                      = !is.null(nuisance_params),
                                            diffusion_HMC                        = TRUE,
                                            diffusion_HMC_integrator             = "kick_flow_kick",  ## signature
                                            partitioned_HMC                      = FALSE,  ## R_fn_sample.R:1606
                                            ##
                                            burnin_algorithm                     = "ESJD", ## --------
                                            ##
                                            n_adapt                              = n_burnin - round(n_burnin / 10),
                                            clip_iter                            = NULL, ## 20,
                                            clip_iter_tau                        = NULL, ## 50,
                                            n_refresh                            = round(n_burnin / 20),
                                            parallel_method                      = "RcppParallel",
                                            ##
                                            n_superchains                        = 4,
                                            ##
                                            learning_rate                        = 0.05,
                                            # learning_rate                        = 0.025,
                                            ##
                                            tau_mult                             = 1.60,  ## R_fn_sample.R:1355
                                            manual_tau                           = FALSE,  ## R_fn_sample.R:1367
                                            tau_if_manual                        = 3.0,  ## R_fn_sample.R:1368
                                            ##
                                            tau_weight_by_p_jump                 = TRUE,  ## R_fn_sample.R:1421
                                            tau_ramp                             = "staged",
                                            eps_reinit_at_ChEES_handover         = FALSE, ## ----
                                            ##
                                            eps_acceptance_mean                  = "harmonic",
                                            ##
                                            eps_reinit_after_pre_burnin          = TRUE,  ## R_fn_sample.R:1425 ## --------
                                            ##
                                            tau_shrink_on_divergence             = FALSE,  ## R_fn_sample.R:1426
                                            tau_shrink_factor                    = 0.95,  ## R_fn_sample.R:1427
                                            tau_shrink_min_divergent_chains      = 2,  ## R_fn_sample.R:1428
                                            ##
                                            share_tau_ii_across_chains_in_burnin = FALSE,  ## R_fn_sample.R:1429
                                            ##
                                            tau_initial                          = "adaptive", ## pi ## --------
                                            tau_adaptation_block                 = "interest_only",  ## --------
                                            interest_only                        = interest_only_params,
                                            ##
                                            tau_sampling_scale                   = "none",  ## R_fn_sample.R:1446
                                            tau_jitter_burnin                    = "uniform",  ## signature
                                            tau_jitter_sampling_type             = "uniform",
                                            ##
                                            randomize_tau_burnin                 = TRUE,  ## signature
                                            randomize_tau_sampling               = TRUE,  ## signature
                                            ##
                                            theta_hat_us_rule                    = "running_mean_frozen",
                                            ##
                                            metric_type_main                     = "Empirical",
                                            ratio_M_main                         = 0.90,  ## R_fn_sample.R:1686
                                            interval_width_main                  = 1,  ## R_fn_sample.R:1689
                                            ##
                                            metric_type_nuisance                 = "Empirical",
                                            metric_shape_nuisance                = "diag",  ## R_fn_sample.R:1697
                                            ratio_M_nuisance                     = 0.90,  ## R_fn_sample.R:1700
                                            interval_width_nuisance              = 1,  ## R_fn_sample.R:1703
                                            ##
                                            max_eps_main                         = 1.0,  ## R_fn_sample.R:1708
                                            max_eps_nuisance                     = 1.0,  ## R_fn_sample.R:1709
                                            max_tau_main                         = 50,
                                            max_tau_nuisance                     = 50,
                                            max_L                                = 1024,  ## R_fn_sample.R:1764
                                            ##
                                            ## (M_decay_*: helper fns 433 type, 436 scale, 440 power)
                                            M_decay_type                         = "inverse",
                                            M_decay_power                        = 0.5,
                                            M_decay_scale                        = (n_burnin - round(n_burnin / 10)) / 5,
                                            ##
                                            beta1_adam                           = 0,  ## R_fn_sample.R:1620
                                            beta2_adam                           = 0.95,  ## R_fn_sample.R:1621
                                            eps_adam                             = 1e-8,  ## R_fn_sample.R:1622
                                            ##
                                            force_autodiff                       = FALSE,  ## R_fn_sample.R:1624
                                            force_PartialLog                     = FALSE,  ## R_fn_sample.R:1625
                                            multi_attempts                       = TRUE,  ## R_fn_sample.R:1626
                                            ##
                                            metric_estimator                     = "pooled",  ## signature
                                            metric_pooled_offdiagonal_shrinkage  = "adaptive",  ## signature
                                            pre_burnin_n_iter                    = NULL,  ## signature
                                            pre_burnin_L                         = NULL,  ## signature
                                            ##
                                            burnin_schedule                      = "automatic",  ## signature
                                            esjd_jump_power                      = 2,  ## signature
                                            tau_cost_exponent                    = 1,  ## signature
                                            ##
                                            metric_update_rescales_eps_and_tau   = TRUE, ## signature ## ----
                                            ##
                                            pooled_metric_warm_up                = 8, ## says how many draws it waits for before its first estimate replaces the starting identity metric.
                                            metric_start_iter                    = round(n_burnin/5), ## says WHEN the metric starts collecting draws.
                                            ##
                                            auto_chunking                        = TRUE)
      return(NicoStan_args)

}
##
## ---- z-scores of a NicoStan fit against a cmdstanr fit of the same model ------------------------------------------------------------
##
## z = (mean_NicoStan - mean_cmdstanr) / sqrt(MCSE_NicoStan^2 + MCSE_cmdstanr^2), MCSE = sd / sqrt(n_eff); when both samplers
## sample the same posterior, z ~ N(0, 1) over the parameters.
##
fn_MCSE_z_scores_of_a_NicoStan_fit_against_a_cmdstanr_fit <- function( tibble_NicoStan,
                                                                       tibble_cmdstanr ) {

      tibble_both <- dplyr::inner_join(
                       dplyr::select(tibble_NicoStan, parameter, mean_NicoStan = mean, sd_NicoStan = sd,
                                     n_eff_NicoStan = n_eff),
                       dplyr::select(tibble_cmdstanr, parameter, mean_cmdstanr = mean, sd_cmdstanr = sd,
                                     n_eff_cmdstanr = n_eff),
                       by = "parameter")
      ##
      tibble_both <- dplyr::mutate(tibble_both,
                                   MCSE_NicoStan = sd_NicoStan / sqrt(n_eff_NicoStan),
                                   MCSE_cmdstanr = sd_cmdstanr / sqrt(n_eff_cmdstanr),
                                   z = (mean_NicoStan - mean_cmdstanr) / sqrt(MCSE_NicoStan^2 + MCSE_cmdstanr^2))
      ##
      ## (constants, e.g. the -1 fillers of unused cutpoints, have sd 0 in both fits: no z)
      tibble_both <- dplyr::filter(tibble_both, is.finite(z), sd_NicoStan > 0, sd_cmdstanr > 0)
      return(tibble_both)

}
##
## ---- Does a set of z-scores follow N(0, 1)? ----------------------------------------------------------------------------------------
##
fn_check_z_scores_follow_N_0_1 <- function( z,
                                            label ) {

      n_z <- length(z)
      tibble_check <- tibble::tibble(
              set_of_parameters            = label,
              n_parameters                 = n_z,
              mean_of_z                    = mean(z),
              SD_of_z                      = sd(z),
              percent_abs_z_above_1.96     = 100 * mean(abs(z) > qnorm(0.975)),
              percent_abs_z_above_2.58     = 100 * mean(abs(z) > qnorm(0.995)),
              largest_abs_z                = max(abs(z)),
              ## (the largest |z| that n_z draws of N(0, 1) exceed with probability 5%)
              largest_abs_z_allowed        = qnorm(1 - 0.025 / n_z),
              p_value_of_Kolmogorov_Smirnov = suppressWarnings(stats::ks.test(z, "pnorm")$p.value))
      tibble_check$passed <- (tibble_check$largest_abs_z <= tibble_check$largest_abs_z_allowed) &&
                             (tibble_check$p_value_of_Kolmogorov_Smirnov > 0.01)
      return(tibble_check)

}
