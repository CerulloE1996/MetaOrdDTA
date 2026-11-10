##
## ==== run_one_cmdstanr_PO_fit_alpha_1_Gat_prior_SD_intercept_given.R
##
## ---- One cmdstanr prior-only fit of the R&G random-kappa model and one test with the runner's alpha = (1, ..., 1)
##      prior, except that the SD of the intercept prior (beta_mu) and the SD of the log-scale prior (raw_scale_mu)
##      are given ("runner" keeps the runner's own value); everything else as in
##      run_one_cmdstanr_PO_fit_alpha_1_runner_settings.R.
##
## ---- Usage: Rscript <this file> rg_rand <test> <n_chains> <n_burnin> <n_iter> <output_dir> <seed>
##                     <SD of the intercept | runner> <SD of the log scale | runner>
##                     <mean of the intercept | runner> <mean of the log scale | runner> <number of studies>
##
{
      command_line_arguments <-  commandArgs(trailingOnly = TRUE)
      stopifnot(length(command_line_arguments) == 12)
      model_of_this_fit <-  command_line_arguments[1]
      test_of_this_fit <-  command_line_arguments[2]
      n_chains <-  as.numeric(command_line_arguments[3])
      n_burnin <-  as.numeric(command_line_arguments[4])
      n_iter <-  as.numeric(command_line_arguments[5])
      output_dir <-  command_line_arguments[6]
      MCMC_seed <-  as.numeric(command_line_arguments[7])
      SD_of_the_intercept_given <-  command_line_arguments[8]
      SD_of_the_log_scale_given <-  command_line_arguments[9]
      mean_of_the_intercept_given <-  command_line_arguments[10]
      mean_of_the_log_scale_given <-  command_line_arguments[11]
      n_studies_of_the_data_set <-  command_line_arguments[12]
      stopifnot(model_of_this_fit == "rg_rand",
                test_of_this_fit %in% c("HADS", "GAD_2", "BAI"))
}
local_pkg_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA"
sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
##
suppressPackageStartupMessages({
      require(NicoStan)
      require(MetaOrdDTA)
      require(cmdstanr)
      require(dplyr)
})
cmdstanr::set_cmdstan_path()
##
## ---- The NicoStan settings, adapt_delta and max_treedepth of the runner's own $sample() call:
##
source(file.path(sim_study_dir, "run_MULTIPLE_MetaOrd_simulation.R"))
fn_find_the_call_with_the_argument <- function(expression_to_search,
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
sample_call_of_the_runner <-  fn_find_the_call_with_the_argument(body(run_multiple_sims), "NicoStan_args")
stopifnot(!is.null(sample_call_of_the_runner))
##
## ---- The model the runner fits (its fused AVX version) and its nuisance blocks, as the runner names them:
##
model_parameterisation <-  if (model_of_this_fit == "rg_rand") "R&G" else "ord_bivariate"
random_thresholds <-  model_of_this_fit != "xu_fixed"
Dirichlet_random_effects_type <-  if (random_thresholds) "kappa" else "fixed"
name_of_the_model_family <-  if (model_parameterisation == "ord_bivariate") "DTA_MA_ord_biv_" else "DTA_MA_Gat_"
custom_file_name <-  paste0(name_of_the_model_family,
                            if (random_thresholds) "RANDthr_kappa" else "FIXEDthr",
                            "_NicoStan_fused_AVX.stan")
nuisance_params <-  MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
      name_of_the_model      = sub("[.]stan$", "", custom_file_name),
      model_parameterisation = model_parameterisation,
      random_thresholds      = random_thresholds)
environment_of_the_runner_values <-  list2env(list(n_burnin        = n_burnin,
                                                   nuisance_params = nuisance_params,
                                                   Stan_file_of_the_fused_version = NULL,
                                                   sampler_settings_to_test = list()))
NicoStan_args <-  eval(sample_call_of_the_runner$NicoStan_args, envir = environment_of_the_runner_values)
adapt_delta <-  eval(sample_call_of_the_runner$adapt_delta)
max_treedepth <-  eval(sample_call_of_the_runner$max_treedepth)
##
## ---- A saved 10-study data set of this test (the likelihood is off, so only its sizes matter):
##
folder_of_the_saved_results <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study",
                                          "simulation_results")
end_of_the_file_names <-  paste0("N_500n_studies_", n_studies_of_the_data_set,
                                 "param_ord_bivariatesoft_FALSErand_thr_FALSEalpha_pi_one_results_list.RDS")
files_of_the_test <-  list.files(folder_of_the_saved_results)
files_of_the_test <-  files_of_the_test[grepl(paste0("_test_", test_of_this_fit, "red_fact_"), files_of_the_test,
                                              fixed = TRUE) &
                                        grepl("DGM_ord_biv_fixed_C", files_of_the_test, fixed = TRUE) &
                                        endsWith(files_of_the_test, end_of_the_file_names)]
stopifnot(length(files_of_the_test) > 0)
seeds_of_the_files <-  as.numeric(sub("^seed_([0-9]+)_.*$", "\\1", files_of_the_test))
saved_results_file <-  file.path(folder_of_the_saved_results, files_of_the_test[which.min(seeds_of_the_files)])
x <-  readRDS(saved_results_file)$x
n_studies_actual <-  nrow(x[[1]])
n_thr <-  ncol(x[[1]]) - 1
##
## ---- The model, prepared as the runner prepares it:
##
CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
AVX_FLAGS_of_this_machine <-  if (CPU_has_AVX512) {
      "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
} else {
      "-march=native -mtune=native -mavx2"
}
NicoSIMD_user_header <-  file.path(dirname(local_pkg_dir), "NicoStan", "inst", "NicoSIMD", "inst", "include",
                                   "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp")
model_prep_obj <-  MetaOrdDTA::MetaOrd_model$new(debugging = FALSE,
                                                 x = x,
                                                 indicator_index_test_in_study = NULL,
                                                 intercept_only = TRUE,
                                                 cov_data = NULL,
                                                 n_chains = n_chains,
                                                 network = FALSE,
                                                 compound_symmetry = NULL,
                                                 hetero_sigma = NULL,
                                                 prior_only = FALSE,
                                                 model_parameterisation = model_parameterisation,
                                                 random_thresholds = random_thresholds,
                                                 Dirichlet_random_effects_type = Dirichlet_random_effects_type,
                                                 custom_file_name = custom_file_name,
                                                 advanced_compile = TRUE,
                                                 force_recompile = FALSE,
                                                 quiet = TRUE,
                                                 set_custom_CXX_CPP_flags = TRUE,
                                                 custom_cpp_user_header_file_path = NicoSIMD_user_header,
                                                 FMA_FLAGS = "-mfma",
                                                 AVX_FLAGS = AVX_FLAGS_of_this_machine)
##
## ---- The alpha = (1, ..., 1) priors, as the runner's alpha_pi == "one" branches set them:
##
priors <-  model_prep_obj$priors
priors$prior_dirichlet_alpha <-  list(rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
if (model_parameterisation == "R&G") {
      priors$prior_dirichlet_alpha <-  rep(1.0, n_thr + 1)
      priors$prior_kappa_mean <-  log(100)
      priors$prior_kappa_SD <-  1.0
      priors$kappa_lb <-  0.001
} else {
      priors$prior_dirichlet_alpha <-  list(rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
      priors$prior_kappa_mean <-  c(log(100), log(100))
      priors$prior_kappa_SD <-  c(1.0, 1.0)
      priors$kappa_lb <-  0.001
}
priors$use_probit_link <-  1
priors$K_fold_CV_indicator <-  rep(1, n_studies_actual)
if (model_parameterisation == "R&G") {
      priors$prior_raw_scale_mu_SD <-  1
      priors$prior_raw_scale_SD_SD <-  1
}
init_lists_per_chain <-  MetaOrdDTA::resize_init_list(init_lists_per_chain = model_prep_obj$init_lists_per_chain,
                                                      n_chains_new = n_chains)
if (random_thresholds) {
      if (model_parameterisation == "ord_bivariate") {
            priors$prior_beta_SD_SD <-  c(0.5, 0.5)
            priors$prior_kappa_mean <-  c(log(100), log(100))
            priors$prior_kappa_SD <-  c(1.0, 1.0)
            priors$prior_kappa_dof <-  c(5, 3)
            for (kk in 1:n_chains) init_lists_per_chain[[kk]]$kappa <-  c(2000, 500)
      } else {
            priors$prior_beta_SD_SD <-  c(0.5, 0.5)
            priors$prior_kappa_mean <-  log(100)
            priors$prior_kappa_SD <-  1.0
            priors$prior_kappa_dof <-  5
            for (kk in 1:n_chains) init_lists_per_chain[[kk]]$kappa <-  1000
      }
}
if (model_parameterisation == "ord_bivariate") {
      priors$prior_beta_SD_SD <-  c(0.5, 0.5)
} else {
      priors$prior_beta_SD_SD <-  0.5
      priors$prior_raw_scale_SD_SD <-  0.5
}
##
## ---- The SD of the intercept prior and of the log-scale prior given (the runner's own values printed):
##
message(paste0("  the runner's R&G priors: beta_mu ~ N(", priors$prior_beta_mu_mean[1], ", ", priors$prior_beta_mu_SD[1],
               "), raw_scale_mu ~ N(", priors$prior_raw_scale_mu_mean[1], ", ", priors$prior_raw_scale_mu_SD[1], ")"))
if (SD_of_the_intercept_given != "runner") priors$prior_beta_mu_SD[1] <-  as.numeric(SD_of_the_intercept_given)
if (SD_of_the_log_scale_given != "runner") priors$prior_raw_scale_mu_SD[1] <-  as.numeric(SD_of_the_log_scale_given)
if (mean_of_the_intercept_given != "runner") priors$prior_beta_mu_mean[1] <-  as.numeric(mean_of_the_intercept_given)
if (mean_of_the_log_scale_given != "runner") priors$prior_raw_scale_mu_mean[1] <-  as.numeric(mean_of_the_log_scale_given)
message(paste0("  the priors of this fit: beta_mu ~ N(", priors$prior_beta_mu_mean[1], ", ", priors$prior_beta_mu_SD[1],
               "), raw_scale_mu ~ N(", priors$prior_raw_scale_mu_mean[1], ", ", priors$prior_raw_scale_mu_SD[1], "), alpha ",
               paste(unique(priors$prior_dirichlet_alpha), collapse = " ")))
##
## ---- Every setting printed with its source:
##
message(paste0("Prior-only fit of ", model_of_this_fit, " (", custom_file_name, ") on ", test_of_this_fit, ": ",
               n_studies_actual, " studies, ", n_thr, " thresholds; data set ", basename(saved_results_file)))
message(paste0("  chains ", n_chains, ", burn-in ", n_burnin, ", sampling iterations ", n_iter, ", seed ",
               MCMC_seed, "; from the runner's $sample() call: adapt_delta ", adapt_delta, ", max_treedepth ",
               max_treedepth))
message(paste0("  nuisance blocks (the runner's fn_nuisance_params_of_a_MetaOrdDTA_model): ",
               paste(nuisance_params, collapse = ", ")))
message(paste0("  priors: alpha ", paste(unique(unlist(priors$prior_dirichlet_alpha)), collapse = " "),
               "; kappa mean ", paste(round(priors$prior_kappa_mean, 3), collapse = ", "), ", SD ",
               paste(priors$prior_kappa_SD, collapse = ", "), ", dof ", paste(priors$prior_kappa_dof, collapse = ", "),
               "; beta_SD_SD ", paste(priors$prior_beta_SD_SD, collapse = ", ")))
for (name_of_the_setting in names(NicoStan_args)) {
      if (name_of_the_setting == "nuisance_params") next
      message(paste0("  runner's NicoStan_args$", name_of_the_setting, " = ",
                     paste(NicoStan_args[[name_of_the_setting]], collapse = ", ")))
}
##
## ---- The fit:
##
name_of_the_test_shown <-  c(HADS = "HADS-A", GAD_2 = "GAD-2", BAI = "BAI")[[test_of_this_fit]]
cutoffs_of_interest_of_the_test <-  list(HADS = c(8, 11), GAD_2 = c(3), BAI = c(8, 16, 26))[[test_of_this_fit]]
check <-  model_prep_obj$run_prior_predictive_check_of_Se_and_Sp(
      priors = priors,
      init_lists_per_chain = init_lists_per_chain,
      test_names = name_of_the_test_shown,
      cutoffs_of_interest_of_each_test = list(cutoffs_of_interest_of_the_test),
      lowest_possible_score_of_each_test = 0,
      n_chains = n_chains,
      n_burnin = n_burnin,
      n_iter = n_iter,
      seed = MCMC_seed,
      algorithm = "Stan",
      NicoStan_args = NicoStan_args,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      probs = c(0.025, 0.975),
      show_what_the_alphas_alone_say = FALSE)
##
result <-  list(model = model_of_this_fit,
                test = test_of_this_fit,
                table_of_the_priors = as.data.frame(check$summary),
                min_ESS = check$min_ESS,
                min_ESS_of_each_type = check$min_ESS_of_each_type,
                n_divergent_transitions = check$n_divergent_transitions,
                n_draws = prod(dim(check$draws_array)[1:2]),
                seconds_of_the_fit = check$seconds_of_the_fit,
                NicoStan_args = NicoStan_args,
                n_chains = n_chains,
                n_burnin = n_burnin,
                n_iter = n_iter,
                seed = MCMC_seed,
                draws_of_Se_Sp_and_kappa = check$draws_array[, , grepl("^(Se_baseline|Sp_baseline|kappa)",
                                                                       dimnames(check$draws_array)[[3]]),
                                                             drop = FALSE])
output_file <-  file.path(output_dir,
                          paste0(model_of_this_fit, "_alpha_one_", test_of_this_fit, "_cmdstanr_SD_of_the_intercept_",
                                 SD_of_the_intercept_given, "_SD_of_the_log_scale_", SD_of_the_log_scale_given, "_means_",
                                 mean_of_the_intercept_given, "_", mean_of_the_log_scale_given, "_n_studies_",
                                 n_studies_of_the_data_set, "_seed_", MCMC_seed, ".rds"))
saveRDS(result, output_file)
utils::write.csv(result$table_of_the_priors, sub("\\.rds$", ".csv", output_file), row.names = FALSE)
message(paste0("DONE ", model_of_this_fit, " ", test_of_this_fit, ": ",
               formatC(result$seconds_of_the_fit, format = "f", digits = 1), " s; min ESS bulk / tail / SD ",
               round(result$min_ESS_of_each_type$min_bulk_ESS), " / ", round(result$min_ESS_of_each_type$min_tail_ESS),
               " / ", round(result$min_ESS_of_each_type$min_SD_ESS), "; divergent ", result$n_divergent_transitions,
               " of ", result$n_draws))






















