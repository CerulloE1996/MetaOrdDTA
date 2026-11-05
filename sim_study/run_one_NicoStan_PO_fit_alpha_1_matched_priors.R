##
## ==== run_one_NicoStan_PO_fit_alpha_1_matched_priors.R ==================================================
##
## ---- One NicoStan prior-only fit of one model and one test with the alpha = (1, ..., 1) priors of the sim study
##      (run_MULTIPLE_MetaOrd_simulation.R, its alpha_pi == "one" branches), the model the runner fits (its fused AVX
##      file; the pop-median file for the random-cutpoint models; DTA_MA_Jones.stan for Jones), its NicoStan settings,
##      adapt_delta and max_treedepth taken from the runner's own $sample() call, sampled with the likelihood switched
##      off (run_prior_predictive_check_of_Se_and_Sp), so that the draws of Se_baseline and Sp_baseline are the
##      induced prior. The hyperparameters given replace the runner's own ("runner" keeps the runner's value).
##
## ---- Usage: Rscript <this file> <model> <test> <n_chains> <n_burnin> <n_iter> <output_dir> <seed>
##                     <SD 1 | runner> <SD 2 | runner> <SD 3 | runner> [<mean 1 | runner> <mean 2 | runner>
##                     [<mean 3 | runner>]]
##      model: xu_fixed (the reference), xu_rand_popmed, rg_rand_popmed or jones; test: HADS, GAD_2 or BAI
##      rg_rand_popmed: SD 1 = SD of the intercept prior (prior_beta_mu_SD), SD 2 = SD of the log-scale prior
##                      (prior_raw_scale_mu_SD), SD 3 unused
##      jones:          SD 1 = prior_beta_mu_SD (both classes), SD 2 = prior_raw_scale_mu_SD (both classes),
##                      SD 3 = prior_boxcox_lambda_SD (both classes); mean 1 = prior_beta_mu_mean (both classes),
##                      mean 2 = prior_raw_scale_mu_mean (both classes), mean 3 = prior_boxcox_lambda_mean (both
##                      classes); the means are optional
##      xu_fixed, xu_rand_popmed: the three SDs unused ("runner")
##
{
      command_line_arguments <-  commandArgs(trailingOnly = TRUE)
      stopifnot(length(command_line_arguments) %in% c(10, 12, 13))
      model_of_this_fit <-  command_line_arguments[1]
      test_of_this_fit <-  command_line_arguments[2]
      n_chains <-  as.numeric(command_line_arguments[3])
      n_burnin <-  as.numeric(command_line_arguments[4])
      n_iter <-  as.numeric(command_line_arguments[5])
      output_dir <-  command_line_arguments[6]
      MCMC_seed <-  as.numeric(command_line_arguments[7])
      SD_1_given <-  command_line_arguments[8]
      SD_2_given <-  command_line_arguments[9]
      SD_3_given <-  command_line_arguments[10]
      mean_1_given <-  if (length(command_line_arguments) >= 12) command_line_arguments[11] else "runner"
      mean_2_given <-  if (length(command_line_arguments) >= 12) command_line_arguments[12] else "runner"
      mean_3_given <-  if (length(command_line_arguments) == 13) command_line_arguments[13] else "runner"
      stopifnot(model_of_this_fit %in% c("xu_fixed", "xu_rand_popmed", "rg_rand_popmed", "jones"),
                test_of_this_fit %in% c("HADS", "GAD_2", "BAI"))
      ## (the kappa prior of the random-cutpoint models: "current" = the runner's log(300), SD 0.5, dof 1000;
      ##  "paper" = the paper's log(100), SD 1.0, dof c(5, 3) / 5; "wide" = log(150), SD 1.25, dof 1000;
      ##  chosen by the environment variable KAPPA_PRIOR)
      kappa_prior_of_this_fit <-  Sys.getenv("KAPPA_PRIOR", unset = "current")
      stopifnot(kappa_prior_of_this_fit %in% c("current", "paper", "wide"))
}
##
## ---- The machine: the local HPC / laptop or ALICE (the frozen library of the sim first in .libPaths()):
##
on_ALICE <-  dir.exists("/home/e/ec376")
if (on_ALICE) {
      ## (NICOSTAN_FROZEN_LIBRARY names another installed copy of the sim's packages, e.g. the one with the
      ##  initial step-size search fix, for a seed whose burn-in collapsed under the frozen one)
      frozen_library_of_this_sim <-  normalizePath(Sys.getenv("NICOSTAN_FROZEN_LIBRARY",
                                                              unset = "/home/e/ec376/R_library_4.3.1"))
      local_pkg_dir <-  "/home/e/ec376/R_packages/MetaOrdDTA"
      folder_of_the_saved_results <-  "/home/e/ec376/MetaOrdDTA_outputs/sim_study/simulation_results"
} else {
      frozen_library_of_this_sim <-  "/home/enzocerullo/R/frozen_library_of_the_MetaOrdDTA_sim"
      local_pkg_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA"
      folder_of_the_saved_results <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study",
                                                "simulation_results")
}
.libPaths(c(frozen_library_of_this_sim, .libPaths()))
Sys.setenv(R_LIBS = frozen_library_of_this_sim)
sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
##
suppressPackageStartupMessages({
      require(NicoStan)
      require(MetaOrdDTA)
      require(cmdstanr)
      require(dplyr)
})
stopifnot(startsWith(find.package("NicoStan"), frozen_library_of_this_sim),
          startsWith(find.package("MetaOrdDTA"), frozen_library_of_this_sim))
##
## ---- The NicoStan settings, adapt_delta and max_treedepth of the runner's own $sample() call:
##
source(file.path(sim_study_dir, "run_MULTIPLE_MetaOrd_simulation.R"))
fn_find_the_call_with_the_argument <- function(expression_to_search,
                                               name_of_the_argument) {

      if (is.call(expression_to_search)) {
            if (name_of_the_argument %in% names(expression_to_search)) return(expression_to_search)
            for (part in as.list(expression_to_search)) {
                  ## (an empty argument, as in x[, 1], cannot be passed on: skipped)
                  call_found <-  tryCatch(fn_find_the_call_with_the_argument(part, name_of_the_argument),
                                          error = function(e) NULL)
                  if (!is.null(call_found)) return(call_found)
            }
      }
      return(NULL)

}
sample_call_of_the_runner <-  fn_find_the_call_with_the_argument(body(run_multiple_sims), "NicoStan_args")
stopifnot(!is.null(sample_call_of_the_runner))
##
## ---- The model the runner fits and its nuisance blocks, as the runner names them (run_MULTIPLE 142-171, 576-579,
##      MetaOrdDTA_run_sims.R, "The pop-median model, its user header and its Se / Sp prior fields"):
##
model_parameterisation <-  switch(model_of_this_fit,
                                  xu_fixed       = "ord_bivariate",
                                  xu_rand_popmed = "ord_bivariate",
                                  rg_rand_popmed = "R&G",
                                  jones          = "Jones")
random_thresholds <-  model_of_this_fit %in% c("xu_rand_popmed", "rg_rand_popmed")
Dirichlet_random_effects_type <-  if (random_thresholds) "kappa" else "fixed"
is_a_pop_median_model <-  random_thresholds
use_custom_file <-  FALSE
custom_file_name <-  NULL
if (model_parameterisation %in% c("ord_bivariate", "R&G")) {
      use_custom_file  <- TRUE
      name_of_the_model_family <- if (model_parameterisation == "ord_bivariate") "DTA_MA_ord_biv_" else
                                  "DTA_MA_ord_HSROC_"
      custom_file_name <- paste0(name_of_the_model_family,
                                 if (random_thresholds) "RANDthr_kappa" else "FIXEDthr",
                                 "_NicoStan_fused_AVX.stan")
}
if (is_a_pop_median_model) {
      custom_file_name <-  paste0(name_of_the_model_family,
                                  "RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan")
      prior_fields <-  list(
            use_Se_Sp_prior          = 0,
            Se_Sp_prior_type         = 4,
            n_chosen_max             = 1,
            n_chosen                 = c(0L, 0L),
            chosen_cutoffs           = matrix(1L, 2, 1),
            prior_Se_Sp_dirichlet_A  = list(rep(0, 2), rep(0, 2)),
            prior_Se_Sp_beta         = list(matrix(0, 1, 2), matrix(0, 1, 2)),
            prior_Se_Sp_link_mean_SD = list(matrix(0, 1, 2), matrix(0, 1, 2)))
}
## (the nuisance blocks come from MetaOrdDTA itself, never written here: one source of truth)
nuisance_params <-  MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
      name_of_the_model      = if (model_parameterisation == "Jones") "DTA_MA_Jones" else
                               sub("[.]stan$", "", custom_file_name),
      model_parameterisation = model_parameterisation,
      random_thresholds      = random_thresholds)
## (run_MULTIPLE 582-595)
interest_only_params <- if (model_parameterisation == "ord_bivariate") {
  
                                  if (random_thresholds) c("beta_mu", "dirichlet_phi")
                                  else c("beta_mu", "C_raw_vec")
  
                        } else if (model_parameterisation == "R&G") {
                          
                                  if (random_thresholds) c("beta_mu", "raw_scale_mu", "dirichlet_phi")
                                  else c("beta_mu", "raw_scale_mu", "C_raw_vec")
                          
                        } else {
                              
                                  NULL
                          
                        }
environment_of_the_runner_values <-  list2env(list(n_burnin             = n_burnin,
                                                   nuisance_params      = nuisance_params,
                                                   interest_only_params = interest_only_params))
NicoStan_args <-  eval(sample_call_of_the_runner$NicoStan_args, envir = environment_of_the_runner_values)
adapt_delta <-  eval(sample_call_of_the_runner$adapt_delta)
max_treedepth <-  eval(sample_call_of_the_runner$max_treedepth)
## (the sim study's models are sampled as the paper sampled them, with no auto-chunking library)
NicoStan_args$auto_chunking <-  FALSE
## (Jones has no fused AVX file: NicoStan samples its own Stan file, R_fn_sample_model_with_NicoStan.R 40-47)
if (model_parameterisation == "Jones") {
      NicoStan_args$Stan_file_of_the_fused_version <-  system.file("stan_models", "stan_models_MA", "DTA_MA_Jones.stan",
                                                                    package = "MetaOrdDTA")
      stopifnot(nzchar(NicoStan_args$Stan_file_of_the_fused_version))
}
##
## ---- A saved 10-study data set of this test (the likelihood is off, so only its sizes matter):
##
files_of_the_test <-  list.files(folder_of_the_saved_results)
files_of_the_test <-  files_of_the_test[grepl(paste0("_test_", test_of_this_fit, "red_fact_"), files_of_the_test,
                                              fixed = TRUE) &
                                        grepl("n_studies_10param_", files_of_the_test, fixed = TRUE) &
                                        endsWith(files_of_the_test, "_results_list.RDS")]
stopifnot(length(files_of_the_test) > 0)
seeds_of_the_files <-  as.numeric(sub("^seed_([0-9]+)_.*$", "\\1", files_of_the_test))
saved_results_file <-  file.path(folder_of_the_saved_results, files_of_the_test[which.min(seeds_of_the_files)])
x <-  readRDS(saved_results_file)$x
n_studies_actual <-  nrow(x[[1]])
n_thr <-  ncol(x[[1]]) - 1
##
## ---- The model, prepared as the runner prepares it (run_MULTIPLE 160-213):
##
CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
AVX_FLAGS_of_this_machine <-  if (CPU_has_AVX512) {
      "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
} else {
      "-march=native -mtune=native -mavx2"
}
NicoSIMD_user_header <-  file.path(find.package("NicoSIMD"), "include", "NicoSIMD",
                                   "stan_external_functions_of_NicoSIMD.hpp")
if (is_a_pop_median_model) {
      ## (the same two includes as inst/stan_functions/RC_popmed_user_header.hpp, from the frozen NicoSIMD)
      NicoSIMD_user_header <-  file.path(frozen_library_of_this_sim, "RC_popmed_user_header.hpp")
      stopifnot(file.exists(NicoSIMD_user_header))
}
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
                                                 advanced_compile = use_custom_file,
                                                 force_recompile = FALSE,
                                                 quiet = TRUE,
                                                 set_custom_CXX_CPP_flags = use_custom_file,
                                                 custom_cpp_user_header_file_path = if (use_custom_file)
                                                                           NicoSIMD_user_header else NULL,
                                                 FMA_FLAGS = "-mfma",
                                                 AVX_FLAGS = AVX_FLAGS_of_this_machine)
##
## ---- The alpha = (1, ..., 1) priors, as the runner's alpha_pi == "one" branches set them (run_MULTIPLE 290-330,
##      433-452, 456-487, 598-600):
##
priors <-  model_prep_obj$priors
if (model_parameterisation %in% c("ord_bivariate", "R&G")) {
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
}
priors$K_fold_CV_indicator <-  rep(1, n_studies_actual)
if (model_parameterisation == "R&G") {
      priors$prior_raw_scale_mu_SD <-  1
      priors$prior_raw_scale_SD_SD <-  1
}
init_lists_per_chain <-  MetaOrdDTA::resize_init_list(init_lists_per_chain = model_prep_obj$init_lists_per_chain,
                                                      n_chains_new = n_chains)
if (random_thresholds == TRUE) {
      if (model_parameterisation == "ord_bivariate") {
            priors$prior_beta_SD_SD <- c(0.5, 0.5)
            if (kappa_prior_of_this_fit == "paper") {
                  priors$prior_kappa_mean <-  c(log(100), log(100))
                  priors$prior_kappa_SD <-  c(1.0, 1.0)
                  priors$prior_kappa_dof <-  c(5, 3)
            } else if (kappa_prior_of_this_fit == "wide") {
                  priors$prior_kappa_mean <-  c(log(150), log(150))
                  priors$prior_kappa_SD <-  c(1.25, 1.25)
                  priors$prior_kappa_dof <-  c(1000, 1000)
            } else {
            priors$prior_kappa_mean <- c(log(300), log(300))
            priors$prior_kappa_SD <- c(0.5, 0.5)
            priors$prior_kappa_dof <- c(1000, 1000) ## default = c(5, 5)#
            }
            for (kk in 1:n_chains) {
              init_lists_per_chain[[kk]]$kappa <- if (kappa_prior_of_this_fit == "wide") c(150, 150) else c(2000, 500)
            }
      } else {
            priors$prior_beta_SD_SD <- c(0.5, 0.5)
            if (kappa_prior_of_this_fit == "paper") {
                  priors$prior_kappa_mean <-  log(100)
                  priors$prior_kappa_SD <-  1.0
                  priors$prior_kappa_dof <-  5
            } else if (kappa_prior_of_this_fit == "wide") {
                  priors$prior_kappa_mean <-  log(150)
                  priors$prior_kappa_SD <-  1.25
                  priors$prior_kappa_dof <-  1000
            } else {
            priors$prior_kappa_mean <- log(300)
            priors$prior_kappa_SD <- 0.5
            priors$prior_kappa_dof <- 1000
            }
            for (kk in 1:n_chains) {
              init_lists_per_chain[[kk]]$kappa <- if (kappa_prior_of_this_fit == "wide") 150 else 1000
            }
      }
}
if (model_parameterisation == "ord_bivariate") {
      priors$prior_beta_SD_SD <- c(0.5, 0.5)
} else if (model_parameterisation == "R&G") {
      priors$prior_beta_SD_SD      <- 0.5
      priors$prior_raw_scale_SD_SD <- 0.5
} else if (model_parameterisation == "Jones") {
      priors$prior_beta_SD_SD      <- c(0.5, 0.5)
      priors$prior_raw_scale_SD_SD <- c(0.5, 0.5)
}
if (is_a_pop_median_model) {
      priors[names(prior_fields)] <- prior_fields
}
## (PRIOR_BETA_SD_SD, PRIOR_RAW_SCALE_SD_SD = a number: the SD of the half-normal prior of the between-study SDs)
if (nzchar(Sys.getenv("PRIOR_BETA_SD_SD"))) {
      priors$prior_beta_SD_SD <-  rep(as.numeric(Sys.getenv("PRIOR_BETA_SD_SD")),
                                      length(priors$prior_beta_SD_SD))
}
if (nzchar(Sys.getenv("PRIOR_RAW_SCALE_SD_SD"))) {
      priors$prior_raw_scale_SD_SD <-  rep(as.numeric(Sys.getenv("PRIOR_RAW_SCALE_SD_SD")),
                                           length(priors$prior_raw_scale_SD_SD))
}
##
## ---- The hyperparameters given (the runner's own values printed first):
##
if (model_parameterisation == "R&G") {
      message(paste0("  the runner's R&G priors: beta_mu ~ N(", priors$prior_beta_mu_mean[1], ", ",
                     priors$prior_beta_mu_SD[1], "), raw_scale_mu ~ N(", priors$prior_raw_scale_mu_mean[1], ", ",
                     priors$prior_raw_scale_mu_SD[1], ")"))
      if (SD_1_given != "runner") priors$prior_beta_mu_SD[1] <-  as.numeric(SD_1_given)
      if (SD_2_given != "runner") priors$prior_raw_scale_mu_SD[1] <-  as.numeric(SD_2_given)
      message(paste0("  the priors of this fit: beta_mu ~ N(", priors$prior_beta_mu_mean[1], ", ",
                     priors$prior_beta_mu_SD[1], "), raw_scale_mu ~ N(", priors$prior_raw_scale_mu_mean[1], ", ",
                     priors$prior_raw_scale_mu_SD[1], ")"))
}
if (model_parameterisation == "Jones") {
      message(paste0("  the runner's Jones priors: beta_mu SD ", paste(priors$prior_beta_mu_SD[, 1], collapse = " / "),
                     ", raw_scale_mu SD ", paste(priors$prior_raw_scale_mu_SD[, 1], collapse = " / "),
                     ", boxcox lambda SD ", paste(priors$prior_boxcox_lambda_SD, collapse = " / ")))
      if (SD_1_given != "runner") priors$prior_beta_mu_SD[, 1] <-  as.numeric(SD_1_given)
      if (SD_2_given != "runner") priors$prior_raw_scale_mu_SD[, 1] <-  as.numeric(SD_2_given)
      if (SD_3_given != "runner") priors$prior_boxcox_lambda_SD <-  rep(as.numeric(SD_3_given), 2)
      message(paste0("  the runner's Jones prior means: beta_mu ",
                     paste(priors$prior_beta_mu_mean[, 1], collapse = " / "), ", raw_scale_mu ",
                     paste(priors$prior_raw_scale_mu_mean[, 1], collapse = " / "), ", boxcox lambda ",
                     paste(priors$prior_boxcox_lambda_mean, collapse = " / ")))
      if (mean_1_given != "runner") priors$prior_beta_mu_mean[, 1] <-  as.numeric(mean_1_given)
      if (mean_2_given != "runner") priors$prior_raw_scale_mu_mean[, 1] <-  as.numeric(mean_2_given)
      if (mean_3_given != "runner") priors$prior_boxcox_lambda_mean <-  rep(as.numeric(mean_3_given), 2)
      message(paste0("  the priors of this fit: beta_mu N(",
                     paste(priors$prior_beta_mu_mean[, 1], collapse = " / "),
                     ", ", paste(priors$prior_beta_mu_SD[, 1], collapse = " / "),
                     "), raw_scale_mu N(", paste(priors$prior_raw_scale_mu_mean[, 1], collapse = " / "), ", ",
                     paste(priors$prior_raw_scale_mu_SD[, 1], collapse = " / "),
                     "), boxcox lambda N(", paste(priors$prior_boxcox_lambda_mean, collapse = " / "), ", ",
                     paste(priors$prior_boxcox_lambda_SD, collapse = " / "), ")"))
}
##
## ---- Every setting printed with its source:
##
message(paste0("Prior-only fit of ", model_of_this_fit, " (",
               if (is.null(custom_file_name)) "DTA_MA_Jones.stan" else custom_file_name, ") on ", test_of_this_fit,
               ": ", n_studies_actual, " studies, ", n_thr, " thresholds; data set ", basename(saved_results_file)))
message(paste0("  chains ", n_chains, ", burn-in ", n_burnin, ", sampling iterations ", n_iter, ", seed ",
               MCMC_seed, "; from the runner's $sample() call: adapt_delta ", adapt_delta, ", max_treedepth ",
               max_treedepth))
message(paste0("  nuisance blocks (the runner's fn_nuisance_params_of_a_MetaOrdDTA_model): ",
               paste(nuisance_params, collapse = ", ")))
if (model_parameterisation != "Jones") {
      message(paste0("  priors: alpha ", paste(unique(unlist(priors$prior_dirichlet_alpha)), collapse = " "),
                     "; kappa mean ", paste(round(priors$prior_kappa_mean, 3), collapse = ", "), ", SD ",
                     paste(priors$prior_kappa_SD, collapse = ", "), ", dof ",
                     paste(priors$prior_kappa_dof, collapse = ", "),
                     "; beta_SD_SD ", paste(priors$prior_beta_SD_SD, collapse = ", ")))
}
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
      algorithm = "NicoStan",
      NicoStan_args = NicoStan_args,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      probs = c(0.025, 0.975),
      show_what_the_alphas_alone_say = FALSE)
##
result <-  list(model = model_of_this_fit,
                test = test_of_this_fit,
                SDs_given = c(SD_1_given, SD_2_given, SD_3_given),
                means_given = c(mean_1_given, mean_2_given, mean_3_given),
                kappa_prior = kappa_prior_of_this_fit,
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
output_file <-  file.path(output_dir, paste0(model_of_this_fit, "_", test_of_this_fit, "_SDs_", SD_1_given, "_",
                                             SD_2_given, "_", SD_3_given,
                                             if (mean_1_given != "runner" || mean_2_given != "runner")
                                                   paste0("_means_", mean_1_given, "_", mean_2_given) else "",
                                             if (mean_3_given != "runner") paste0("_", mean_3_given) else "",
                                             if (random_thresholds) paste0("_kappa_", kappa_prior_of_this_fit) else "",
                                             "_seed_", MCMC_seed, ".rds"))
saveRDS(result, output_file)
utils::write.csv(result$table_of_the_priors, sub("\\.rds$", ".csv", output_file), row.names = FALSE)
message(paste0("DONE ", model_of_this_fit, " ", test_of_this_fit, " SDs ", SD_1_given, " ", SD_2_given, " ",
               SD_3_given, " means ", mean_1_given, " ", mean_2_given, " ", mean_3_given, " seed ", MCMC_seed, ": ",
               formatC(result$seconds_of_the_fit, format = "f", digits = 1), " s; min ESS bulk / tail / SD ",
               round(result$min_ESS_of_each_type$min_bulk_ESS), " / ", round(result$min_ESS_of_each_type$min_tail_ESS),
               " / ", round(result$min_ESS_of_each_type$min_SD_ESS), "; divergent ", result$n_divergent_transitions,
               " of ", result$n_draws))
