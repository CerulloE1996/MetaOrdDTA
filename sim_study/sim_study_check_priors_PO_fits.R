##
## ==== sim_study_check_priors_PO_fits.R ===============================================
##
## ---- Prior-only fits of the sim study's models with NicoStan, at the sim runner's own settings (its $sample()
##      call, its nuisance blocks and its parameters of interest), one test per process (the test is the first
##      argument): the median and both 95% limits of Se and Sp at every cutoff, with ESS and divergences.
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      ##
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
      checks_dir <-  file.path(sim_study_dir, "prior_derivation_tools_checks")
      alpha_sets_dir <-  file.path(sim_study_dir, "prior_only_checks_of_the_induced_Dirichlet_alpha_sets")
      outputs_dir <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study")
      test <-  commandArgs(trailingOnly = TRUE)[1]
      name_of_the_model <-  commandArgs(trailingOnly = TRUE)[2]
      n_studies_of_the_data <-  as.numeric(commandArgs(trailingOnly = TRUE)[3])
      stopifnot(n_studies_of_the_data %in% c(10, 50))
      ## the sim study's chains (MetaOrdDTA_run_sims.R line 734) and NicoStan burn-in (run_MULTIPLE_... line 332)
      n_chains <-  4
      n_burnin <-  500
      ## the runner's n_iter for 4 chains at 10 and 50 studies (run_MULTIPLE_MetaOrd_simulation.R lines 285-310)
      n_iter_of_each_model <-  if (n_studies_of_the_data == 10) {
            c(ord_biv_FC = 3500, Bivariate = 10000, ord_biv_RC = 1500, ordHSROC_RC = 1500)
      } else {
            c(ord_biv_FC = 2000, Bivariate = 10000, ord_biv_RC = 1000, ordHSROC_RC = 1000)
      }
      custom_file_of_each_model <-  c(ord_biv_FC = "DTA_MA_ord_biv_FIXEDthr_NicoStan_fused_AVX.stan",
                                      Bivariate = "DTA_MA_strat_bivariate_NicoStan_fused_AVX.stan",
                                      ord_biv_RC = "DTA_MA_ord_biv_RANDthr_kappa_NicoStan_fused_AVX.stan",
                                      ordHSROC_RC = "DTA_MA_Gat_RANDthr_kappa_NicoStan_fused_AVX.stan")
      parameterisation_of_each_model <-  c(ord_biv_FC = "ord_bivariate",
                                           Bivariate = "Bivariate",
                                           ord_biv_RC = "ord_bivariate",
                                           ordHSROC_RC = "R&G")
      random_thresholds_of_each_model <-  c(ord_biv_FC = FALSE,
                                            Bivariate = FALSE,
                                            ord_biv_RC = TRUE,
                                            ordHSROC_RC = TRUE)
      Dirichlet_type_of_each_model <-  c(ord_biv_FC = "fixed",
                                         Bivariate = "fixed",
                                         ord_biv_RC = "kappa",
                                         ordHSROC_RC = "kappa")
      seed <-  1
}
##
## ---- the sim runner of the local-HPC (the current sim study), the same file on both machines:
file_of_the_runner <-  if (Sys.info()[["user"]] == "enzocerullo") {
      file.path(sim_study_dir, "run_MULTIPLE_MetaOrd_simulation.R")
} else {
      "/tmp/run_MULTIPLE_MetaOrd_simulation_of_the_local_HPC.R"
}
source(file_of_the_runner)
fn_find_the_call_with_the_argument <-  function(expression_to_search,
                                                name_of_the_argument) {

      if (is.call(expression_to_search)) {
            if (name_of_the_argument %in% names(expression_to_search)) return(expression_to_search)
            for (part in as.list(expression_to_search)) {
                  call_found <-  fn_find_the_call_with_the_argument(part, name_of_the_argument)
                  if (!is.null(call_found)) return(call_found)
            }
      }
      return(NULL)

}
fn_find_the_assignment_of <-  function(expression_to_search,
                                       name_of_the_variable) {

      if (is.call(expression_to_search)) {
            if (identical(expression_to_search[[1]], as.name("<-")) &&
                identical(expression_to_search[[2]], as.name(name_of_the_variable))) return(expression_to_search[[3]])
            for (part in as.list(expression_to_search)) {
                  if (missing(part)) next
                  found <-  fn_find_the_assignment_of(part, name_of_the_variable)
                  if (!is.null(found)) return(found)
            }
      }
      return(NULL)

}
sample_call_of_the_runner <-  fn_find_the_call_with_the_argument(body(run_multiple_sims), "NicoStan_args")
interest_only_of_the_runner <-  fn_find_the_assignment_of(body(run_multiple_sims), "interest_only_params")
stopifnot(!is.null(sample_call_of_the_runner), !is.null(interest_only_of_the_runner))
##
CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
AVX_FLAGS_of_this_machine <-  if (CPU_has_AVX512) {
      "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
} else {
      "-march=native -mtune=native -mavx2"
}
NicoSIMD_user_header <-  normalizePath(file.path(local_pkg_dir, "..", "NicoStan", "inst", "NicoSIMD", "inst",
                                                 "include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp"))
##
files <-  list.files(file.path(outputs_dir, "simulation_results"))
files <-  files[grepl(paste0("_test_", test, "red_fact_"), files, fixed = TRUE) &
                grepl("DGM_ord_biv_fixed_C", files, fixed = TRUE) &
                endsWith(files, paste0("N_500n_studies_", n_studies_of_the_data,
                                       "param_ord_bivariatesoft_FALSErand_thr_FALSE",
                                       "alpha_pi_one_results_list.RDS"))]
x <-  readRDS(file.path(outputs_dir, "simulation_results",
                        files[which.min(as.numeric(sub("^seed_([0-9]+)_.*$", "\\1", files)))]))$x
n_studies <-  nrow(x[[1]])
##
## ---- The priors of each cell: flat-1 from the saved reference fit, informative from the derived alphas:
##
reference_of_flat_1 <-  readRDS(file.path(checks_dir,
                                          paste0("reference_prior_draws_of_the_fused_ord_biv_fixed_model_at_alpha_one_",
                                                 test, "_16_chains_lower_bound_minus_15.rds")))
derived <-  readRDS(file.path(alpha_sets_dir,
                              paste0("informative_from_the_literature_by_Stan_sampling_", test, ".rds")))
results_of_each_cell <-  list()
## (the ordinal-HSROC random-cutpoint model has no informative prior derived yet: flat-1 only)
priors_to_check <-  if (name_of_the_model == "ordHSROC_RC") "flat_1" else c("flat_1", "informative")
for (prior in priors_to_check) {

      ## (the Bivariate model has no prior-only switch: every count of the data set set to missing)
      x_of_the_fit <-  if (name_of_the_model == "Bivariate") lapply(x, function(m) { m[, -1] <-  -1; m }) else x
      model_prep_obj <-  MetaOrdDTA::MetaOrd_model$new(debugging = FALSE,
                                                       x = x_of_the_fit,
                                                       indicator_index_test_in_study = NULL,
                                                       intercept_only = TRUE,
                                                       cov_data = NULL,
                                                       n_chains = n_chains,
                                                       network = FALSE,
                                                       compound_symmetry = NULL,
                                                       hetero_sigma = NULL,
                                                       prior_only = FALSE,
                                                       model_parameterisation =
                                                             parameterisation_of_each_model[[name_of_the_model]],
                                                       random_thresholds =
                                                             random_thresholds_of_each_model[[name_of_the_model]],
                                                       Dirichlet_random_effects_type =
                                                             Dirichlet_type_of_each_model[[name_of_the_model]],
                                                       custom_file_name = custom_file_of_each_model[[name_of_the_model]],
                                                       advanced_compile = TRUE,
                                                       force_recompile = FALSE,
                                                       quiet = TRUE,
                                                       set_custom_CXX_CPP_flags = TRUE,
                                                       custom_cpp_user_header_file_path = NicoSIMD_user_header,
                                                       FMA_FLAGS = "-mfma",
                                                       AVX_FLAGS = AVX_FLAGS_of_this_machine)
      priors <-  if (prior == "flat_1" && name_of_the_model == "ord_biv_FC") reference_of_flat_1$priors else model_prep_obj$priors
      if (name_of_the_model == "Bivariate") {
            derived_Bivariate <-  readRDS(file.path(checks_dir,
                                                    paste0("Bivariate_", prior, "_prior_by_Stan_sampling_", test,
                                                           "_at_the_lower_bound_minus_15.rds")))$fit_of_each_quantity
            priors$prior_mu_mean <-  list(derived_Bivariate$Se$mu_mean, derived_Bivariate$Sp$mu_mean)
            priors$prior_mu_SD <-  list(derived_Bivariate$Se$mu_SD, derived_Bivariate$Sp$mu_SD)
      }
      init_lists_per_chain <-  MetaOrdDTA::resize_init_list(init_lists_per_chain = model_prep_obj$init_lists_per_chain,
                                                            n_chains_new = n_chains)
      if (name_of_the_model == "ord_biv_RC") {
            ## (the runner's random-cutpoint priors and kappa starts, run_MULTIPLE_MetaOrd_simulation.R lines 352-360)
            priors$prior_beta_SD_SD <-  c(0.5, 0.5)
            priors$prior_kappa_mean <-  c(log(100), log(100))
            priors$prior_kappa_SD <-  c(1.0, 1.0)
            priors$prior_kappa_dof <-  c(5, 3)
            for (kk in 1:n_chains) init_lists_per_chain[[kk]]$kappa <-  c(2000, 500)
            if (prior == "flat_1") {
                  priors$prior_beta_mu_mean <-  reference_of_flat_1$priors$prior_beta_mu_mean
                  priors$prior_beta_mu_SD <-  reference_of_flat_1$priors$prior_beta_mu_SD
            }
      }
      if (name_of_the_model == "ordHSROC_RC") {
            ## (the runner's ordinal-HSROC random-cutpoint flat-1 priors and kappa starts,
            ##  run_MULTIPLE_MetaOrd_simulation.R lines 237-243, 362-369 and 389-391)
            priors$prior_dirichlet_alpha <-  rep(1.0, length(priors$prior_dirichlet_alpha))
            priors$prior_kappa_mean <-  log(100)
            priors$prior_kappa_SD <-  1.0
            priors$prior_kappa_dof <-  5
            priors$kappa_lb <-  0.001
            priors$prior_beta_SD_SD <-  0.5
            priors$prior_raw_scale_SD_SD <-  0.5
            for (kk in 1:n_chains) init_lists_per_chain[[kk]]$kappa <-  1000
            ##
            ## ---- location from the reference's own prior; raw scale from its Stan derivation:
            ##
            derived_flat_1 <-  readRDS(file.path(checks_dir,
                                                 paste0("derived_prior_ord_HSROC_FIXEDthr_flat_1_", test,
                                                        ".rds")))$medians
            priors$prior_beta_mu_mean[] <-  reference_of_flat_1$priors$prior_beta_mu_mean[1, 1]
            priors$prior_beta_mu_SD[] <-  reference_of_flat_1$priors$prior_beta_mu_SD[1, 1]
            priors$prior_raw_scale_mu_mean[] <-  derived_flat_1[["prior_raw_scale_mu_intercept_mean"]]
            priors$prior_raw_scale_mu_SD[] <-  derived_flat_1[["prior_raw_scale_mu_SD"]]
      }
      if (prior == "informative" && name_of_the_model %in% c("ord_biv_FC", "ord_biv_RC")) {
            priors$prior_dirichlet_alpha <-  derived$alpha
            priors$prior_beta_mu_mean[1, 1] <-  derived$location_mean[1]
            priors$prior_beta_mu_mean[2, 1] <-  derived$location_mean[2]
            priors$prior_beta_mu_SD[1, 1] <-  derived$location_SD[1]
            priors$prior_beta_mu_SD[2, 1] <-  derived$location_SD[2]
            priors$prior_beta_SD_SD <-  c(0.5, 0.5)
            priors$use_probit_link <-  1
      }
      priors$K_fold_CV_indicator <-  rep(0, n_studies)
      ##
      ## ---- the runner's NicoStan settings, nuisance blocks and parameters of interest:
      ##
      model_parameterisation <-  parameterisation_of_each_model[[name_of_the_model]]
      random_thresholds <-  random_thresholds_of_each_model[[name_of_the_model]]
      nuisance_params <-  MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                               name_of_the_model = sub("[.]stan$", "", custom_file_of_each_model[[name_of_the_model]]),
                               model_parameterisation = model_parameterisation,
                               random_thresholds = random_thresholds)
      ## (nico_reduce_sum does not take the Bivariate's array[n_thr, n_studies] vector[2] z yet: z sampled as main)
      if (name_of_the_model == "Bivariate") nuisance_params <-  NULL
      interest_only_params <-  eval(interest_only_of_the_runner)
      NicoStan_args <-  eval(sample_call_of_the_runner$NicoStan_args,
                             envir = list2env(list(Stan_file_of_the_fused_version = NULL,
                                                   sampler_settings_to_test = list())))
      message(paste0("\033[36m", test, ", ", prior, ": NicoStan, ", n_chains, " chains x (", n_burnin, " + ",
                     n_iter_of_each_model[[name_of_the_model]], ") = ",
                     n_chains * (n_burnin + n_iter_of_each_model[[name_of_the_model]]), " iterations; nuisance ",
                     paste(nuisance_params, collapse = ", "), "; interest ",
                     paste(interest_only_params, collapse = ", "), "\033[0m"))
      seconds_at_the_start <-  Sys.time()
      model_samples_obj <-  model_prep_obj$sample(seed = seed,
                                                  n_burnin = n_burnin,
                                                  n_iter = n_iter_of_each_model[[name_of_the_model]],
                                                  adapt_delta = eval(sample_call_of_the_runner$adapt_delta),
                                                  max_treedepth = eval(sample_call_of_the_runner$max_treedepth),
                                                  metric_shape = eval(sample_call_of_the_runner$metric_shape),
                                                  priors = priors,
                                                  init_lists_per_chain = init_lists_per_chain,
                                                  n_chains = n_chains,
                                                  algorithm = "NicoStan",
                                                  NicoStan_args = NicoStan_args)
      seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), seconds_at_the_start, units = "secs"))
      outs <-  model_prep_obj$internal_obj$outs_stan_sampling
      print(names(outs))
      draws_array <-  outs$draws_array
      names_of_the_draws <-  dimnames(draws_array)[[3]]
      n_thr <-  sum(startsWith(names_of_the_draws, "Se_baseline["))
      Se <-  matrix(draws_array[, , paste0("Se_baseline[", seq_len(n_thr), "]")], ncol = n_thr)
      Sp <-  matrix(draws_array[, , paste0("Sp_baseline[", seq_len(n_thr), "]")], ncol = n_thr)
      ESS <-  sapply(c(paste0("Se_baseline[", seq_len(n_thr), "]"), paste0("Sp_baseline[", seq_len(n_thr), "]")),
                     function(name) {

                           c(bulk = posterior::ess_bulk(draws_array[, , name]),
                             tail = posterior::ess_tail(draws_array[, , name]))

                     })
      message(paste0("\033[32m", test, ", ", prior, ": NicoStan prior-only fit in ", round(seconds_of_the_fit, 1),
                     " s; min ESS of Se and Sp bulk ", round(min(ESS["bulk", ])), ", tail ",
                     round(min(ESS["tail", ])), "\033[0m"))
      results_of_each_cell[[prior]] <-  list(Se = Se,
                                            Sp = Sp,
                                            ESS = ESS,
                                            divergences = outs$divergences,
                                            seconds_of_the_fit = seconds_of_the_fit,
                                            priors = priors,
                                            NicoStan_args = NicoStan_args)
      saveRDS(results_of_each_cell,
              file.path(checks_dir, paste0("NicoStan_prior_only_fits_of_the_", name_of_the_model, "_model_", test,
                                           if (n_studies_of_the_data == 50) "_50_studies" else "", ".rds")))

}






















