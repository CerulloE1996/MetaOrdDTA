##
## ---- Prior sets 1 and 2 on every sim-study model -----------------------------------------------------------------------------------
##
## Prior-only NicoStan fits of every sim-study model with the direct Se / Sp prior of prior set 1 and of prior
## set 2, for every test: a Normal on probit(Se) and on probit(Sp) at the main cutoffs of the test, nothing typed
## at the other cutoffs (or N(0, 1) at them); the runner's priors everywhere else. Prints the Se / Sp each model
## implies at every cutoff and checks the main cutoffs against the prior set within 0.05. Rscript arguments: Stan
## file names (otherwise every file of stan_model_files below), "tests=GAD_2,HADS", "prior_sets=1",
## "n_iter=1000", "n_chains=16", "cutoffs=N_0_1".
##
## ---- Settings ----------------------------------------------------------------------------------------------------------------------
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      ##
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
      checks_dir <-  file.path(sim_study_dir, "prior_derivation_tools_checks")
      outputs_dir <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study")
      ##
      tests_to_test <-  c("GAD_2", "HADS", "BAI")
      prior_sets_to_test <-  c(1, 2)
      stan_model_files <-  c("DTA_MA_Jones_Se_Sp_prior.stan",
                             "DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior_NicoStan_fused_AVX.stan",
                             "DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan",
                             "DTA_MA_ord_HSROC_RANDthr_kappa_Se_Sp_prior_NicoStan_fused_AVX.stan",
                             "DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan")
      prior_at_the_other_cutoffs <-  "blank"  ## "blank" (nothing typed) or "N_0_1"
      n_chains <-  2             ## a prior check: <= 5 minutes, 2 chains, <= 2 CPUs
      n_burnin <-  125           ## NicoStan
      n_burnin_of_Stan <-  250   ## cmdstanr (the stratified-bivariate model, whose z NicoStan does not hold)
      n_iter <-  250
      ##
      args <-  commandArgs(trailingOnly = TRUE)
      if (any(grepl("[.]stan$", args))) stan_model_files <-  args[grepl("[.]stan$", args)]
      ##
      ## ---- Rscript argument tests= -------------------------------------------------------------------------------------------------
      ##
      if (any(grepl("^tests=", args))) {
            tests_to_test <-  strsplit(sub("^tests=", "", args[grepl("^tests=", args)]), ",")[[1]]
      }
      ##
      ## ---- Rscript argument prior_sets= --------------------------------------------------------------------------------------------
      ##
      if (any(grepl("^prior_sets=", args))) {
            prior_sets_to_test <-  as.integer(strsplit(sub("^prior_sets=", "", args[grepl("^prior_sets=", args)]),
                                                       ",")[[1]])
      }
      if (any(grepl("^n_iter=", args))) n_iter <-  as.integer(sub("^n_iter=", "", args[grepl("^n_iter=", args)]))
      if (any(grepl("^n_burnin=", args))) {
            n_burnin <-  as.integer(sub("^n_burnin=", "", args[grepl("^n_burnin=", args)]))
            n_burnin_of_Stan <-  n_burnin
      }
      ##
      ## ---- Rscript argument cutoffs= -----------------------------------------------------------------------------------------------
      ##
      if (any(grepl("^cutoffs=", args))) {
            prior_at_the_other_cutoffs <-  sub("^cutoffs=", "", args[grepl("^cutoffs=", args)])
      }
      ##
      ## ---- Rscript argument n_chains= ----------------------------------------------------------------------------------------------
      ##
      if (any(grepl("^n_chains=", args))) {
            n_chains <-  as.integer(sub("^n_chains=", "", args[grepl("^n_chains=", args)]))
      }
      ##
      ## ---- Rscript argument algorithm= (NicoStan, or Stan = cmdstanr) ----------------------------------------------------------------
      ##
      sampling_algorithm <-  "NicoStan"
      if (any(grepl("^algorithm=", args))) {
            sampling_algorithm <-  sub("^algorithm=", "", args[grepl("^algorithm=", args)])
      }
      stopifnot(sampling_algorithm %in% c("NicoStan", "Stan"))
      stopifnot(prior_at_the_other_cutoffs %in% c("blank", "N_0_1"), all(prior_sets_to_test %in% c(1, 2)))
      ## (output names carry n_iter and n_chains when they are not the defaults)
      file_tag <-  paste0(if (sampling_algorithm == "Stan") "_cmdstanr" else "",
                          if (n_burnin != 125) paste0("_n_burnin_", n_burnin) else "",
                          if (n_iter != 250) paste0("_n_iter_", n_iter) else "",
                          if (n_chains != 2) paste0("_n_chains_", n_chains) else "")
      seed <-  1
      Se_Sp_prior_type <-  4
}
##
## ---- Prior sets 1 and 2 at the main cutoffs ----------------------------------------------------------------------------------------
##
## The median of Se and of Sp at each main cutoff (score >=): midway between the lowest and highest true value
## of the four DGMs (2 dp).
## Set 1: probit = N(a, b), a = qnorm(median), b = |qnorm(median) - qnorm(inner)| / 1.96, the inner 95% limit
## 0.15 from the median towards 0.5. Set 2: N(a + 1.96 b, b), its lower 95% limit the median of set 1.
##
{
      prior_set_medians <-  rbind(
            data.frame(test = "GAD_2", cutoff =  2, Se_median = 0.85, Sp_median = 0.68),
            data.frame(test = "GAD_2", cutoff =  3, Se_median = 0.65, Sp_median = 0.86),
            data.frame(test = "HADS",  cutoff =  8, Se_median = 0.82, Sp_median = 0.73),
            data.frame(test = "HADS",  cutoff = 11, Se_median = 0.57, Sp_median = 0.90),
            data.frame(test = "BAI",   cutoff =  8, Se_median = 0.92, Sp_median = 0.55),
            data.frame(test = "BAI",   cutoff = 16, Se_median = 0.75, Sp_median = 0.81),
            data.frame(test = "BAI",   cutoff = 26, Se_median = 0.47, Sp_median = 0.93))
}
##
## ---- The runner, its sample call and its parameters of interest --------------------------------------------------------------------
##
{
      ##
      ## ---- The sim runner of the local-HPC, the same file on both machines ---------------------------------------------------------
      ##
      file_of_the_runner <-  if (Sys.info()[["user"]] == "enzocerullo") {
            file.path(sim_study_dir, "run_MULTIPLE_MetaOrd_simulation.R")
      } else {
            "/tmp/run_MULTIPLE_MetaOrd_simulation_of_the_local_HPC.R"
      }
      source(file_of_the_runner)
}
##
## ---- Function: the call that has an argument ---------------------------------------------------------------------------------------
##
fn_find_the_call_with_the_argument <-  function(expression_to_search,
                                                name_of_the_argument) {

      ##
      ## ---- A call: itself if it has the argument, else each of its parts -----------------------------------------------------------
      ##
      if (is.call(expression_to_search)) {
            if (name_of_the_argument %in% names(expression_to_search)) return(expression_to_search)
            ##
            ## ---- Each part of the call ---------------------------------------------------------------------------------------------
            ##
            for (part in as.list(expression_to_search)) {
                  if (missing(part)) next
                  call_found <-  fn_find_the_call_with_the_argument(part, name_of_the_argument)
                  if (!is.null(call_found)) return(call_found)
            }
      }
      return(NULL)

}
##
## ---- Function: the value assigned to a variable ------------------------------------------------------------------------------------
##
fn_find_the_assignment_of <-  function(expression_to_search,
                                       name_of_the_variable) {

      ##
      ## ---- A call: itself if it assigns the variable, else each of its parts -------------------------------------------------------
      ##
      if (is.call(expression_to_search)) {
            ##
            ## ---- An assignment to the variable -------------------------------------------------------------------------------------
            ##
            if (identical(expression_to_search[[1]], as.name("<-")) &&
                identical(expression_to_search[[2]], as.name(name_of_the_variable))) {
                  return(expression_to_search[[3]])
            }
            ##
            ## ---- Each part of the call ---------------------------------------------------------------------------------------------
            ##
            for (part in as.list(expression_to_search)) {
                  if (missing(part)) next
                  found <-  fn_find_the_assignment_of(part, name_of_the_variable)
                  if (!is.null(found)) return(found)
            }
      }
      return(NULL)

}
##
## ---- The runner's sampler settings, compiler flags and user headers ----------------------------------------------------------------
##
{
      sample_call_of_the_runner <-  fn_find_the_call_with_the_argument(body(run_multiple_sims), "NicoStan_args")
      interest_only_of_the_runner <-  fn_find_the_assignment_of(body(run_multiple_sims), "interest_only_params")
      stopifnot(!is.null(sample_call_of_the_runner), !is.null(interest_only_of_the_runner))
      adapt_delta <-  eval(sample_call_of_the_runner$adapt_delta)
      ## (output names carry adapt_delta when ADAPT_DELTA is set)
      if (nzchar(Sys.getenv("ADAPT_DELTA"))) file_tag <-  paste0(file_tag, "_adapt_delta_", adapt_delta)
      max_treedepth <-  eval(sample_call_of_the_runner$max_treedepth)
      metric_shape <-  eval(sample_call_of_the_runner$metric_shape)
      ##
      CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
      ##
      ## ---- The SIMD flags of this machine ------------------------------------------------------------------------------------------
      ##
      AVX_FLAGS_of_this_machine <-  if (CPU_has_AVX512) {
            "-march=native -mtune=native -mavx2 -mavx512f -mavx512vl -mavx512dq"
      } else {
            "-march=native -mtune=native -mavx2"
      }
      NicoSIMD_user_header <-  normalizePath(file.path(local_pkg_dir, "..", "NicoStan", "inst", "NicoSIMD",
                                                       "inst", "include", "NicoSIMD",
                                                       "stan_external_functions_of_NicoSIMD.hpp"))
      ## (the pop-median models: NicoSIMD's header and the Beta median of each cutpoint)
      RC_popmed_user_header <-  normalizePath(file.path(local_pkg_dir, "inst", "stan_functions",
                                                        "RC_popmed_user_header.hpp"))
      machine <-  if (CPU_has_AVX512) "local_HPC" else "laptop"
      source(file.path(sim_study_dir, "sim_study_priors_helper_fns.R"))
}
##
## ---- Function: the prior set of a test ---------------------------------------------------------------------------------------------
##
## The probit Normals of Se and Sp at the main cutoffs of the test and the targets they imply (median, 2.5%,
## 97.5%), one column per cutoff:
##
fn_prior_set_of_a_test <-  function(test,
                                    prior_set) {

      rows <-  prior_set_medians[prior_set_medians$test == test, ]
      rows <-  rows[order(rows$cutoff), ]
      stopifnot(nrow(rows) >= 1)
      ##
      ## ---- The probit mean and SD of a median --------------------------------------------------------------------------------------
      ##
      fn_probit_mean_and_SD <-  function(median) {

            inner <-  ifelse(median > 0.5, median - 0.15, median + 0.15)
            SD <-  abs(stats::qnorm(median) - stats::qnorm(inner)) / 1.96
            mean <-  stats::qnorm(median) + if (prior_set == 2) 1.96 * SD else 0
            cbind(mean = mean, SD = SD)

      }
      ##
      ## ---- The median and 95% limits of the probit Normals -------------------------------------------------------------------------
      ##
      fn_targets <-  function(probit) {

            rbind(median = stats::pnorm(probit[, "mean"]),
                  lower = stats::pnorm(probit[, "mean"] - 1.96 * probit[, "SD"]),
                  upper = stats::pnorm(probit[, "mean"] + 1.96 * probit[, "SD"]))

      }
      probit_Se <-  fn_probit_mean_and_SD(rows$Se_median)
      probit_Sp <-  fn_probit_mean_and_SD(rows$Sp_median)
      list(cutoff = rows$cutoff,
           sources = paste0("prior set ", prior_set),
           probit_Se = probit_Se,
           probit_Sp = probit_Sp,
           target_Se = fn_targets(probit_Se),
           target_Sp = fn_targets(probit_Sp),
           n_sources = rep(1, nrow(rows)))

}
##
## ---- Function: the settings of a model from its Stan file name ---------------------------------------------------------------------
##
fn_settings_of_a_model <-  function(stan_model_file) {

      ord_HSROC <-  grepl("ord_HSROC", stan_model_file)
      strat_biv <-  grepl("strat_bivariate", stan_model_file)
      Jones <-  grepl("_Jones", stan_model_file)
      ## (a model file with no Se / Sp prior: its own priors, matched to the prior set by a prior model)
      priors_matched_by_a_prior_model <-  !grepl("_Se_Sp_prior", stan_model_file)
      random_thresholds <-  grepl("RANDthr", stan_model_file)
      model_parameterisation <-  if (strat_biv) "Bivariate" else if (Jones) "Jones" else
                                 if (ord_HSROC) "ord_HSROC" else "ord_biv"
      list(model = if (priors_matched_by_a_prior_model) {
                         paste0(sub("[.]stan$", "", sub("^DTA_MA_", "", stan_model_file)),
                                "_priors_matched_by_a_prior_model")
                   } else {
                         sub("_Se_Sp_prior.*$", "", sub("^DTA_MA_", "", stan_model_file))
                   },
           priors_matched_by_a_prior_model = priors_matched_by_a_prior_model,
           model_parameterisation = model_parameterisation,
           label_of_the_runner = if (strat_biv) "Bivariate" else if (Jones) "Jones" else if (ord_HSROC) "R&G" else
                                 "ord_bivariate",
           random_thresholds = random_thresholds,
           Dirichlet_random_effects_type = if (random_thresholds) "kappa" else "fixed",
           algorithm = if (strat_biv) "Stan" else "NicoStan",
           Se_first = strat_biv,
           cts = Jones,
           user_header = if (grepl("_popmed_", stan_model_file, fixed = TRUE)) RC_popmed_user_header else
                         NicoSIMD_user_header)

}
##
## ---- Function: the runner's priors of a model --------------------------------------------------------------------------------------
##
## (copies of run_MULTIPLE_MetaOrd_simulation.R lines 215-216, 278-305, 410-432, 434-464 and 520-522, alpha_pi =
##  "one"; the stratified-bivariate model: N(0, 1) on probit Se and Sp at every cutoff)
##
fn_runner_priors_of_a_model <-  function(priors,
                                         settings,
                                         n_thr,
                                         test) {

      model_parameterisation <-  settings$label_of_the_runner
      random_thresholds <-  settings$random_thresholds
      ##
      ## ---- The stratified-bivariate model ------------------------------------------------------------------------------------------
      ##
      if (model_parameterisation == "Bivariate") {
            priors$prior_mu_mean <-  list(rep(0, n_thr), rep(0, n_thr))
            priors$prior_mu_SD <-  list(rep(1, n_thr), rep(1, n_thr))
            return(priors)
      }
      ##
      priors$prior_dirichlet_alpha <- list( rep(1.0, n_thr + 1), rep(1.0, n_thr + 1))
      ##
      ##
      ## ---- The Dirichlet and kappa priors: ord_HSROC, else ord_biv and Jones -------------------------------------------------------
      ##
      if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
              priors$prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
              priors$prior_kappa_mean <- log(100)
              priors$prior_kappa_SD <- 1.0
              priors$kappa_lb <- 0.001
      } else {
              prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
              priors$prior_dirichlet_alpha <-  list(prior_dirichlet_alpha, prior_dirichlet_alpha)
              priors$prior_kappa_mean   <-  c(log(100), log(100))
              priors$prior_kappa_SD   <-  c(1.0, 1.0)
              priors$kappa_lb <-  0.001
      }
      ##
      priors$use_probit_link <- 1
      ##
      ##
      ## ---- The raw-scale priors of ord_HSROC ---------------------------------------------------------------------------------------
      ##
      if (model_parameterisation == "R&G") {
          priors$prior_raw_scale_mu_SD <- 1
          priors$prior_raw_scale_SD_SD <- 1
      }
      ##
      ##
      ## ---- The random-cutpoint priors ----------------------------------------------------------------------------------------------
      ##
      if (random_thresholds == TRUE) {
            ##
            ## ---- ord_biv -----------------------------------------------------------------------------------------------------------
            ##
            if (model_parameterisation == "ord_bivariate") {
                priors$prior_beta_SD_SD <- c(0.5, 0.5)
                priors$prior_kappa_mean <- c(log(100), log(100))
                priors$prior_kappa_SD <- c(1.0, 1.0)
                priors$prior_kappa_dof <- c(5, 3) ## default = c(5, 5)#
            } else {
                priors$prior_beta_SD_SD <- c(0.5, 0.5)
                priors$prior_kappa_mean <- log(100)
                priors$prior_kappa_SD <- 1.0
                priors$prior_kappa_dof <- 5
            }
      }
      ##
      ##
      ## ---- The SD priors of each model ---------------------------------------------------------------------------------------------
      ##
      if (model_parameterisation == "ord_bivariate") {
            priors$prior_beta_SD_SD <- c(0.5, 0.5)
      } else if (model_parameterisation == "R&G") {
            priors$prior_beta_SD_SD      <- 0.5
            priors$prior_raw_scale_SD_SD <- 0.5
      } else if (model_parameterisation == "Jones") {
            priors$prior_beta_SD_SD      <- c(0.5, 0.5)
            priors$prior_raw_scale_SD_SD <- c(0.5, 0.5)
            ## (as the runner: the location, log-scale and Box-Cox priors matched over the whole band of thresholds)
            mean_of_beta_mu_matched <-  c(HADS = 6.571, GAD_2 = 2.5, BAI = 9.094)[[test]]
            mean_of_raw_scale_mu_matched <-  c(HADS = 1.497, GAD_2 = 0.708, BAI = 1.432)[[test]]
            mean_of_boxcox_lambda_matched <-  c(HADS = 0.75, GAD_2 = 1.0, BAI = 0.5)[[test]]
            SD_of_beta_mu_matched <-  c(HADS = 4.062, GAD_2 = 2.256, BAI = 4.651)[[test]]
            SD_of_raw_scale_mu_matched <-  c(HADS = 0.25, GAD_2 = 0.5, BAI = 0.1875)[[test]]
            SD_of_boxcox_lambda_matched <-  c(HADS = 0.1, GAD_2 = 0.1, BAI = 0.05)[[test]]
            priors$prior_beta_mu_mean[, 1] <-  mean_of_beta_mu_matched
            priors$prior_raw_scale_mu_mean[, 1] <-  mean_of_raw_scale_mu_matched
            priors$prior_boxcox_lambda_mean <-  mean_of_boxcox_lambda_matched
            priors$prior_beta_mu_SD[, 1] <-  SD_of_beta_mu_matched
            priors$prior_raw_scale_mu_SD[, 1] <-  SD_of_raw_scale_mu_matched
            priors$prior_boxcox_lambda_SD <-  SD_of_boxcox_lambda_matched
      }
      ##
      ##
      ## ---- The fixed-cutpoint ordinal models: anchor at the location ---------------------------------------------------------------
      ##
      priors

}
##
## ---- Function: the runner's starting values of kappa -------------------------------------------------------------------------------
##
## (copies of run_MULTIPLE_MetaOrd_simulation.R lines 418-420 and 427-429)
##
fn_runner_starting_kappa <-  function(init_lists_per_chain,
                                      settings) {

      ##
      ## ---- The random-cutpoint models ----------------------------------------------------------------------------------------------
      ##
      if (settings$random_thresholds == TRUE) {
            ##
            ## ---- Each chain --------------------------------------------------------------------------------------------------------
            ##
            for (kk in seq_along(init_lists_per_chain)) {
                  ##
                  ## ---- ord_biv: one kappa per class, else one kappa ----------------------------------------------------------------
                  ##
                  if (settings$label_of_the_runner == "ord_bivariate") {
                        init_lists_per_chain[[kk]]$kappa <- c(2000, 500)
                  } else {
                        init_lists_per_chain[[kk]]$kappa <- 1000
                  }
            }
      }
      init_lists_per_chain

}
##
## ---- Function: the direct Se / Sp prior fields -------------------------------------------------------------------------------------
##
## The prior set at its cutoffs; the other cutoffs blank or N(0, 1). Class order: Sp then Se in the ordinal and
## Jones models, Se then Sp in the stratified-bivariate model. The Stan files take the mean of probit(Se) as it is
## (they flip its sign for 1 - Se themselves):
##
fn_Se_Sp_prior_fields <-  function(priors,
                                   settings,
                                   prior_set,
                                   n_thr) {

      keep <-  prior_set$cutoff <= n_thr
      ##
      ## ---- The cutoffs typed: every cutoff (N(0, 1) at the others), else the main cutoffs only -------------------------------------
      ##
      if (prior_at_the_other_cutoffs == "N_0_1") {
            chosen_Sp <-  chosen_Se <-  1:n_thr
            link_Sp <-  cbind(rep(0, n_thr), rep(1, n_thr))
            link_Se <-  cbind(rep(0, n_thr), rep(1, n_thr))
            link_Sp[prior_set$cutoff[keep], ] <-  prior_set$probit_Sp[keep, , drop = FALSE]
            link_Se[prior_set$cutoff[keep], ] <-  prior_set$probit_Se[keep, , drop = FALSE]
      } else {
            chosen_Sp <-  chosen_Se <-  prior_set$cutoff[keep]
            link_Sp <-  prior_set$probit_Sp[keep, , drop = FALSE]
            link_Se <-  prior_set$probit_Se[keep, , drop = FALSE]
      }
      n_Sp <-  length(chosen_Sp)
      n_Se <-  length(chosen_Se)
      n_chosen_cutoffs_max <-  max(n_Sp, n_Se)
      fn_pad_cutoffs <-  function(chosen) c(chosen, rep(1L, n_chosen_cutoffs_max - length(chosen)))
      fn_pad_link <-  function(link) rbind(link, cbind(rep(0, n_chosen_cutoffs_max - nrow(link)),
                                                       rep(1, n_chosen_cutoffs_max - nrow(link))))
      priors$use_Se_Sp_prior <-  1
      priors$Se_Sp_prior_type <-  Se_Sp_prior_type
      priors$n_chosen_cutoffs_max <-  n_chosen_cutoffs_max
      priors$prior_Se_Sp_beta <-  list(matrix(0, n_chosen_cutoffs_max, 2), matrix(0, n_chosen_cutoffs_max, 2))
      ##
      ## ---- Class order: Se then Sp, else Sp then Se --------------------------------------------------------------------------------
      ##
      if (settings$Se_first) {
            priors$n_chosen_cutoffs <-  c(n_Se, n_Sp)
            priors$chosen_cutoffs <-  rbind(fn_pad_cutoffs(chosen_Se), fn_pad_cutoffs(chosen_Sp))
            priors$prior_Se_Sp_link_mean_SD <-  list(fn_pad_link(link_Se), fn_pad_link(link_Sp))
      } else {
            priors$n_chosen_cutoffs <-  c(n_Sp, n_Se)
            priors$chosen_cutoffs <-  rbind(fn_pad_cutoffs(chosen_Sp), fn_pad_cutoffs(chosen_Se))
            ##
            ## ---- The ordinal models: the Dirichlet fields --------------------------------------------------------------------------
            ##
            if (!settings$cts) {
                  priors$prior_Se_Sp_dirichlet_A <-  list(rep(0, n_chosen_cutoffs_max + 1), rep(0, n_chosen_cutoffs_max + 1))
            }
            priors$prior_Se_Sp_link_mean_SD <-  list(fn_pad_link(link_Sp), fn_pad_link(link_Se))
      }
      message(paste0("Se/Sp prior typed at ", n_Sp, " Sp cutoffs (", paste(chosen_Sp, collapse = " "), ") and ",
                     n_Se, " Se cutoffs (", paste(chosen_Se, collapse = " "), ")"))
      priors

}
##
## ---- Function: the median and 95% limits of each column ----------------------------------------------------------------------------
##
fn_median_and_95_limits <-  function(draws) {

      rbind(median = apply(draws, 2, stats::median),
            lower = apply(draws, 2, stats::quantile, probs = 0.025),
            upper = apply(draws, 2, stats::quantile, probs = 0.975))

}
##
## ---- Function: one prior-only fit, a model on a test under a prior set -------------------------------------------------------------
##
fn_prior_only_fit <-  function(test,
                               stan_model_file,
                               number_of_the_prior_set) {

      settings <-  fn_settings_of_a_model(stan_model_file)
      lit <-  fn_prior_set_of_a_test(test, number_of_the_prior_set)
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
      x <-  fn_data_set_of_a_test(outputs_dir, test)
      n_studies <-  nrow(x[[1]])
      n_burnin_of_this_fit <-  if (settings$algorithm == "NicoStan") n_burnin else n_burnin_of_Stan
      fn_message_of_the_fit <-  function(n_thr,
                                         nuisance_params,
                                         ##
                                         ## ---- The message of the fit ---------------------------------------------------------------
                                         ##
                                         interest_only_params) {

            message(paste0("\033[36m", test, " x ", settings$model, ": Se/Sp prior type ", Se_Sp_prior_type,
                           ", ", lit$sources, " at cutoffs ", paste(lit$cutoff, collapse = " "),
                           ": probit(Se) ~ N(", paste(signif(lit$probit_Se[, "mean"], 3), collapse = " "),
                           "; SD ", paste(signif(lit$probit_Se[, "SD"], 3), collapse = " "), "), probit(Sp) ~ N(",
                           paste(signif(lit$probit_Sp[, "mean"], 3), collapse = " "), "; SD ",
                           paste(signif(lit$probit_Sp[, "SD"], 3), collapse = " "),
                           "); other cutoffs ", prior_at_the_other_cutoffs, "; ", n_studies, " studies, ", n_thr,
                           " cutoffs, all prior-only; ",
                           settings$algorithm, ", ", n_chains, " chains x (", n_burnin_of_this_fit, " + ", n_iter,
                           ") = ", n_chains * (n_burnin_of_this_fit + n_iter), " iterations; nuisance ",
                           paste(nuisance_params, collapse = ", "),
                           "; interest ", paste(interest_only_params, collapse = ", "), "\033[0m"))

      }
      ##
      ## ---- The stratified-bivariate model (cmdstanr), else NicoStan ----------------------------------------------------------------
      ##
      if (settings$algorithm == "Stan") {
            ##
            ## ---- The stratified-bivariate model: cmdstanr on its compiled executable -----------------------------------------------
            ##
            ## (its z block is array[n_thr, n_studies] vector[2], which NicoStan does not hold; every count -1,
            ##  prior-only)
            Stan_data_list <-  MetaOrdDTA:::R_fn_prep_MA_data(x = x,
                                             cov_data = NULL,
                                             intercept_only = TRUE,
                                             model_parameterisation = "Bivariate",
                                             n_index_tests_per_study = NULL,
                                             indicator_index_test_in_study = NULL)$stan_data_list
            n_thr <-  Stan_data_list$n_thr
            stopifnot(any(lit$cutoff <= n_thr))
            cutpoints_index_array <-  array(-1L, dim = c(2, n_studies, n_thr + 1))
            ##
            ## ---- Each class --------------------------------------------------------------------------------------------------------
            ##
            for (c in 1:2) {
                  cutpoints_index_array[c, , ] <-  Stan_data_list$cutpoint_index[[c]]
            }
            priors <-  list(n_studies = n_studies,
                            n_thr = n_thr,
                            n_obs_cutpoints = as.integer(Stan_data_list$n_obs_cutpoints),
                            x = array(-1L, dim = c(2, n_studies, n_thr + 1)),
                            cutpoint_index = cutpoints_index_array,
                            use_probit_link = 1L)
            priors <-  fn_runner_priors_of_a_model(priors, settings, n_thr, test)
            priors <-  fn_Se_Sp_prior_fields(priors, settings, lit, n_thr)
            initial_values <-  list(mu = replicate(n_thr, c(0, 0), simplify = FALSE),
                                    L_Omega = replicate(n_thr, diag(2), simplify = FALSE),
                                    sigma = replicate(n_thr, c(0.5, 0.5), simplify = FALSE),
                                    z = array(0, dim = c(n_thr, n_studies, 2)))
            nuisance_params <-  "none (cmdstanr)"
            interest_only_params <-  "mu"
            NicoStan_args <-  NULL
            fn_message_of_the_fit(n_thr, nuisance_params, interest_only_params)
            executable <-  system.file("stan_models", "stan_models_MA", sub("[.]stan$", "", stan_model_file),
                                       package = "MetaOrdDTA")
            stopifnot(file.exists(executable))
            model <-  cmdstanr::cmdstan_model(exe_file = executable)
            seconds_at_the_start <-  Sys.time()
            fit <-  model$sample(data = priors,
                                 init = replicate(n_chains, initial_values, simplify = FALSE),
                                 seed = seed,
                                 chains = n_chains,
                                 parallel_chains = n_chains,
                                 iter_warmup = n_burnin_of_Stan,
                                 iter_sampling = n_iter,
                                 adapt_delta = adapt_delta,
                                 max_treedepth = max_treedepth,
                                 refresh = 0,
                                 show_messages = TRUE)
            seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), seconds_at_the_start, units = "secs"))
            draws_array <-  fit$draws(c("Se_baseline", "Sp_baseline"), format = "draws_array")
            n_divergent <-  sum(fit$sampler_diagnostics()[, , "divergent__"])
            outs <-  list(divergences = list(n_divs = n_divergent))
      } else {
            ##
            ## ---- The ordinal and Jones models: NicoStan, the runner's settings -----------------------------------------------------
            ##
            fn_new_model <-  MetaOrdDTA::MetaOrd_model$new
            model_prep_obj <-  fn_new_model(debugging = FALSE,
                                            x = x,
                                            indicator_index_test_in_study = NULL,
                                            intercept_only = TRUE,
                                            cov_data = NULL,
                                            n_chains = n_chains,
                                            network = FALSE,
                                            compound_symmetry = NULL,
                                            hetero_sigma = NULL,
                                            prior_only = FALSE,
                                            model_parameterisation = settings$model_parameterisation,
                                            random_thresholds = settings$random_thresholds,
                                            Dirichlet_random_effects_type =
                                                  settings$Dirichlet_random_effects_type,
                                            custom_file_name = stan_model_file,
                                            advanced_compile = TRUE,
                                            force_recompile = FALSE,
                                            quiet = TRUE,
                                            set_custom_CXX_CPP_flags = TRUE,
                                            custom_cpp_user_header_file_path = settings$user_header,
                                            FMA_FLAGS = "-mfma",
                                            AVX_FLAGS = AVX_FLAGS_of_this_machine)
            priors <-  model_prep_obj$priors
            init_lists_per_chain <-  MetaOrdDTA::resize_init_list(
                                           init_lists_per_chain = model_prep_obj$init_lists_per_chain,
                                           n_chains_new = n_chains)
            init_lists_per_chain <-  fn_runner_starting_kappa(init_lists_per_chain, settings)
            n_thr <-  model_prep_obj$internal_obj$outs_data$n_thr
            stopifnot(any(lit$cutoff <= n_thr))
            priors <-  fn_runner_priors_of_a_model(priors, settings, n_thr, test)
            ##
            ## ---- The model's own priors matched by a prior model, else the direct Se / Sp prior ------------------------------------
            ##
            if (settings$priors_matched_by_a_prior_model) {
                  matched <-  readRDS(file.path(checks_dir, paste0(sub("_priors_matched_by_a_prior_model$", "",
                                                                       settings$model),
                                                                   "_priors_matched_to_prior_sets_1_and_2_by_a_",
                                                                   "prior_model.rds")))
                  matched <-  matched$matched_priors[[paste0("prior_set_", number_of_the_prior_set, "_", test)]]
                  priors$prior_beta_mu_mean[, 1] <-  matched$prior_beta_mu_mean
                  priors$prior_beta_mu_SD[, 1] <-  matched$prior_beta_mu_SD
                  priors$prior_raw_scale_mu_mean[, 1] <-  matched$prior_raw_scale_mu_mean
                  priors$prior_raw_scale_mu_SD[, 1] <-  matched$prior_raw_scale_mu_SD
                  priors$prior_boxcox_lambda_mean <-  matched$prior_boxcox_lambda_mean
                  priors$prior_boxcox_lambda_SD <-  matched$prior_boxcox_lambda_SD
            } else {
                  priors <-  fn_Se_Sp_prior_fields(priors, settings, lit, n_thr)
            }
            priors$K_fold_CV_indicator <-  rep(0, n_studies)
            ##
            ## (the runner's nuisance blocks and parameters of interest)
            model_parameterisation <-  settings$label_of_the_runner
            random_thresholds <-  settings$random_thresholds
            nuisance_params <-  MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                                     name_of_the_model = sub("[.]stan$", "", stan_model_file),
                                     model_parameterisation = settings$model_parameterisation,
                                     random_thresholds = random_thresholds)
            interest_only_params <-  eval(interest_only_of_the_runner)
            NicoStan_args <-  eval(sample_call_of_the_runner$NicoStan_args,
                                   envir = list2env(list(Stan_file_of_the_fused_version = NULL,
                                                         sampler_settings_to_test = list())))
            ## (the runner's own Stan file, e.g. DTA_MA_Jones.stan, replaced by the file checked here)
            if (!is.null(NicoStan_args$Stan_file_of_the_fused_version)) {
                  NicoStan_args$Stan_file_of_the_fused_version <-  system.file("stan_models", "stan_models_MA",
                                                                               stan_model_file,
                                                                               package = "MetaOrdDTA")
                  stopifnot(nzchar(NicoStan_args$Stan_file_of_the_fused_version))
            }
            fn_message_of_the_fit(n_thr, nuisance_params, interest_only_params)
            seconds_at_the_start <-  Sys.time()
            model_samples_obj <-  model_prep_obj$sample(seed = seed,
                                                        n_burnin = n_burnin,
                                                        n_iter = n_iter,
                                                        adapt_delta = adapt_delta,
                                                        max_treedepth = max_treedepth,
                                                        metric_shape = metric_shape,
                                                        priors = priors,
                                                        init_lists_per_chain = init_lists_per_chain,
                                                        n_chains = n_chains,
                                                        algorithm = sampling_algorithm,
                                                        NicoStan_args = NicoStan_args)
            seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), seconds_at_the_start, units = "secs"))
            outs <-  model_prep_obj$internal_obj$outs_stan_sampling
            ##
            ## ---- The draws and the divergences: NicoStan, else cmdstanr -------------------------------------------------------------
            ##
            if (sampling_algorithm == "NicoStan") {
                  draws_array <-  outs$draws_array
                  n_divergent <-  sum(unlist(outs$divergences$n_divs))
            } else {
                  draws_array <-  outs$stan_mod_samples$draws(variables = c("Se_baseline", "Sp_baseline"),
                                                              format = "draws_array")
                  n_divergent <-  sum(outs$stan_mod_samples$diagnostic_summary()$num_divergent)
                  outs$divergences <-  list(n_divs = n_divergent)
            }
      }
      ##
      ## ---- The Se / Sp the prior implies at every cutoff, against the prior set at its cutoffs -------------------------------------
      ##
      cutoffs_assessed <-  1:n_thr
      names_of_Se <-  paste0("Se_baseline[", cutoffs_assessed, "]")
      names_of_Sp <-  paste0("Sp_baseline[", cutoffs_assessed, "]")
      Se <-  matrix(draws_array[, , names_of_Se], ncol = length(cutoffs_assessed))
      Sp <-  matrix(draws_array[, , names_of_Sp], ncol = length(cutoffs_assessed))
      ##
      ## ---- ESS bulk, tail and SD of each Se / Sp (NicoStan's own summary) ----------------------------------------------------------
      ##
      summary_of_Se_and_Sp <-  NicoStan::generate_summary_tibble(n_threads = 1,
                                                                  trace = aperm(draws_array[, , c(names_of_Se, names_of_Sp),
                                                                                            drop = FALSE],
                                                                                c(3, 1, 2)),
                                                                  param_names = c(names_of_Se, names_of_Sp),
                                                                  n_to_compute = length(c(names_of_Se, names_of_Sp)),
                                                                  compute_nested_rhat = FALSE,
                                                                  n_chains = dim(draws_array)[2],
                                                                  n_superchains = 1)
      ESS <-  rbind(bulk = summary_of_Se_and_Sp$n_eff,
                    tail = summary_of_Se_and_Sp$n_eff_tail,
                    SD = summary_of_Se_and_Sp$n_eff_sd)
      fitted <-  list(Se = fn_median_and_95_limits(Se), Sp = fn_median_and_95_limits(Sp))
      colnames(fitted$Se) <-  colnames(fitted$Sp) <-  paste0("cutoff_", cutoffs_assessed)
      message(paste0("\033[36m", test, " x ", settings$model, ", ", lit$sources,
                     ": Se the prior implies at every cutoff:\033[0m"))
      print(round(fitted$Se, 3))
      message(paste0("\033[36m", test, " x ", settings$model, ", ", lit$sources,
                     ": Sp the prior implies at every cutoff:\033[0m"))
      print(round(fitted$Sp, 3))
      keep <-  lit$cutoff <= n_thr
      gaps <-  list(Se = abs(fitted$Se[, lit$cutoff[keep], drop = FALSE] - lit$target_Se[, keep, drop = FALSE]),
                    Sp = abs(fitted$Sp[, lit$cutoff[keep], drop = FALSE] - lit$target_Sp[, keep, drop = FALSE]))
      n_compared <-  length(unlist(gaps))
      largest_gap <-  max(unlist(gaps))
      mean_gap <-  mean(unlist(gaps))
      n_over_0_05 <-  sum(unlist(gaps) > 0.05)
      pass <-  largest_gap <= 0.05
      message(paste0(if (pass) "\033[32m" else "\033[31m", test, " x ", settings$model, ", ", lit$sources,
                     " at cutoffs ", paste(lit$cutoff[keep], collapse = " "), ": ", if (pass) "PASS" else "FAIL",
                     "; largest gap ", signif(largest_gap, 3), "; mean gap ", signif(mean_gap, 3), "; ",
                     n_over_0_05, " of ", n_compared, " over 0.05; largest Se gap per cutoff ",
                     paste(signif(apply(gaps$Se, 2, max), 3), collapse = " / "),
                     "; largest Sp gap per cutoff ",
                     paste(signif(apply(gaps$Sp, 2, max), 3), collapse = " / "), "; ", n_divergent,
                     " divergent; min ESS bulk ", round(min(ESS["bulk", ])), ", tail ", round(min(ESS["tail", ])),
                     "; ", round(seconds_of_the_fit, 1), " s\033[0m"))
      ## (the draws of every parameter other than Se / Sp are kept too, so that the geometry of a failed fit can
      ##  be read from the file without a new fit)
      names_of_the_other_params <-  setdiff(dimnames(draws_array)[[3]], c(names_of_Se, names_of_Sp))
      other_draws <-  draws_array[, , names_of_the_other_params, drop = FALSE]
      result <-  list(test = test,
                      model = settings$model,
                      prior_set = number_of_the_prior_set,
                      stan_model_file = stan_model_file,
                      literature = lit,
                      Se = Se,
                      Sp = Sp,
                      other_draws = other_draws,
                      fitted = fitted,
                      gaps = gaps,
                      pass = pass,
                      ESS = ESS,
                      divergences = outs$divergences,
                      seconds_of_the_fit = seconds_of_the_fit,
                      priors = priors,
                      NicoStan_args = NicoStan_args)
      saveRDS(result, file.path(checks_dir, paste0("PO_check_prior_set_", number_of_the_prior_set, "_", test, "_",
                                                   settings$model, "_", prior_at_the_other_cutoffs, "_", machine,
                                                   file_tag, ".rds")))
      fn_interval <-  function(limits,
                               ##
                               ## ---- A median and its 95% limits as text ------------------------------------------------------------
                               ##
                               j) {

            paste0(round(limits["median", j], 3), " [", round(limits["lower", j], 3), ", ",
                   round(limits["upper", j], 3), "]")

      }
      data.frame(prior_set = number_of_the_prior_set,
                 test = test,
                 model = settings$model,
                 cutoff = lit$cutoff[keep],
                 Se_target = sapply(which(keep), function(j) fn_interval(lit$target_Se, j)),
                 Se_fitted = sapply(lit$cutoff[keep], function(k) fn_interval(fitted$Se, k)),
                 Sp_target = sapply(which(keep), function(j) fn_interval(lit$target_Sp, j)),
                 Sp_fitted = sapply(lit$cutoff[keep], function(k) fn_interval(fitted$Sp, k)),
                 largest_gap = signif(largest_gap, 3),
                 mean_gap = signif(mean_gap, 3),
                 n_over_0_05 = n_over_0_05,
                 n_compared = n_compared,
                 pass = pass,
                 n_divergent = n_divergent,
                 min_ESS_bulk = round(min(ESS["bulk", ])),
                 min_ESS_tail = round(min(ESS["tail", ])),
                 n_iter = n_iter,
                 seconds = round(seconds_of_the_fit, 1))

}
##
## ---- Every model on every test under each prior set --------------------------------------------------------------------------------
##
{
      summary_rows <-  list()
      ##
      ## ---- Each prior set ----------------------------------------------------------------------------------------------------------
      ##
      for (number_of_the_prior_set in prior_sets_to_test) {
            ##
            ## ---- Each model --------------------------------------------------------------------------------------------------------
            ##
            for (stan_model_file in stan_model_files) {
                  ##
                  ## ---- Each test ---------------------------------------------------------------------------------------------------
                  ##
                  for (test in tests_to_test) {
                        summary_rows[[length(summary_rows) + 1]] <-  fn_prior_only_fit(test,
                                                                                      stan_model_file,
                                                                                      number_of_the_prior_set)
                  }
            }
      }
      summary_table <-  do.call(rbind, summary_rows)
      options(width = 220)
      message("\033[36mEvery model on every test, at the main cutoffs, against prior sets 1 and 2:\033[0m")
      print(summary_table, row.names = FALSE)
      write.csv(summary_table,
                file.path(checks_dir, paste0("PO_check_prior_sets_all_tests_", prior_at_the_other_cutoffs, "_",
                                             machine, file_tag, ".csv")),
                row.names = FALSE)
}
##
## ---- Every model against every other, per prior set and test -----------------------------------------------------------------------
##
## The implied Se and Sp (median, 2.5%, 97.5%) of the saved fits of this machine: largest gap, mean gap and the
## count over 0.05, at the main cutoffs (the prior set in force in every model) and at the other assessed cutoffs
## (MetaOrdDTA_run_sims.R lines 491-499: GAD-2 1:6, HADS 1:21, BAI 5:43), where each model runs its own prior:
##
{
      assessed_cutoffs <-  list(GAD_2 = 1:6, HADS = 1:21, BAI = 5:43)
      between_model_rows <-  list()
      ##
      ## ---- Each prior set ----------------------------------------------------------------------------------------------------------
      ##
      for (number_of_the_prior_set in prior_sets_to_test) {
            ##
            ## ---- Each test ---------------------------------------------------------------------------------------------------------
            ##
            for (test in tests_to_test) {
                  files_of_the_test <-  list.files(checks_dir,
                                                   pattern = paste0("^PO_check_prior_set_",
                                                                    number_of_the_prior_set, "_", test, "_.*_",
                                                                    prior_at_the_other_cutoffs, "_", machine,
                                                                    ".*[.]rds$"))
                  fits <-  lapply(file.path(checks_dir, files_of_the_test), readRDS)
                  if (length(fits) < 2) next
                  ## (one fit per model: the one with the most sampling iterations)
                  models_of_the_fits <-  sapply(fits, function(fit) fit$model)
                  n_iter_of_the_fits <-  sapply(fits, function(fit) nrow(fit$Se))
                  fits <-  fits[order(models_of_the_fits, -n_iter_of_the_fits)]
                  fits <-  fits[!duplicated(sapply(fits, function(fit) fit$model))]
                  if (length(fits) < 2) next
                  main_cutoffs <-  fits[[1]]$literature$cutoff
                  cutoff_sets <-  list(main_cutoffs = main_cutoffs,
                                       other_cutoffs = setdiff(assessed_cutoffs[[test]], main_cutoffs))
                  ##
                  ## ---- Each pair of models -----------------------------------------------------------------------------------------
                  ##
                  for (a in seq_len(length(fits) - 1)) {
                        ##
                        ## ---- The second model of the pair --------------------------------------------------------------------------
                        ##
                        for (b in (a + 1):length(fits)) {
                              all_gaps <-  rbind(abs(fn_median_and_95_limits(fits[[a]]$Se) -
                                                     fn_median_and_95_limits(fits[[b]]$Se)),
                                                 abs(fn_median_and_95_limits(fits[[a]]$Sp) -
                                                     fn_median_and_95_limits(fits[[b]]$Sp)))
                              ##
                              ## ---- The main cutoffs, then the other assessed cutoffs -----------------------------------------------
                              ##
                              for (set in names(cutoff_sets)) {
                                    gaps <-  all_gaps[, cutoff_sets[[set]], drop = FALSE]
                                    between_model_rows[[length(between_model_rows) + 1]] <-  data.frame(
                                          prior_set = number_of_the_prior_set,
                                          test = test,
                                          cutoffs = set,
                                          model_a = fits[[a]]$model,
                                          model_b = fits[[b]]$model,
                                          n_cutoffs = ncol(gaps),
                                          largest_gap = signif(max(gaps), 3),
                                          mean_gap = signif(mean(gaps), 3),
                                          n_over_0_05 = sum(gaps > 0.05),
                                          n_compared = length(gaps),
                                          largest_gap_of_medians = signif(max(gaps[c(1, 4), ]), 3),
                                          mean_gap_of_medians = signif(mean(gaps[c(1, 4), ]), 3))
                              }
                        }
                  }
            }
      }
      ##
      ## ---- The table of every pair -------------------------------------------------------------------------------------------------
      ##
      if (length(between_model_rows) > 0) {
            between_model_table <-  do.call(rbind, between_model_rows)
            message(paste0("\033[36mEvery model vs every other, the main cutoffs (prior set) and the other ",
                           "assessed cutoffs (each model's own prior), Se and Sp x median / 2.5% / 97.5%: ",
                           "largest, mean gap, count over 0.05:\033[0m"))
            print(between_model_table, row.names = FALSE)
            write.csv(between_model_table,
                      file.path(checks_dir, paste0("PO_check_prior_sets_models_vs_each_other_",
                                                   prior_at_the_other_cutoffs, "_", machine, ".csv")),
                      row.names = FALSE)
      }
}
