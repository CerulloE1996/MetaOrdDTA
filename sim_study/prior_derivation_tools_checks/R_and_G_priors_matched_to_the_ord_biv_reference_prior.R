##
## ---- Ordinal-HSROC (R&G) priors matched to the ordinal-bivariate reference prior ---------------------------------------------------
##
## The trace of Se and Sp at each threshold from the prior-only fit of the ordinal-bivariate (fixed-cutpoint) model under
## the flat prior (alpha = 1) or under prior set 1 or 2 is the data; each draw of it is given the ordinal-HSROC model's
## own parameters (the location, the log-scale and the cutpoints, the cutpoints under the model's own induced Dirichlet
## prior with alpha = 1), with the model's own priors on the location and the log-scale, and the priors that copy the
## trace are found by Stan sampling (prior_model_for_R_and_G_params.stan): beta_mu ~ normal(prior_beta_mu_mean,
## prior_beta_mu_SD) and raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD). The draws copied are
## n_draws_to_copy of the trace, evenly spaced (each draw carries its own cutpoints).
## Rscript arguments: "tests=GAD_2,HADS", "references=flat,1,2".
##
## ---- Settings ----------------------------------------------------------------------------------------------------------------------
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      ##
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      checks_dir <-  file.path(local_pkg_dir, "sim_study", "prior_derivation_tools_checks")
      tests_to_match <-  c("GAD_2", "HADS", "BAI")
      references_to_match <-  c("flat", "1", "2")
      n_iter_of_the_reference_fit <-  2500   ## (one chain of the prior-only ordinal-bivariate fit under prior set 1 or 2)
      n_draws_to_copy <-  1000
      n_chains <-  4
      n_burnin <-  500
      n_iter <-  500
      seed <-  1
      ##
      args <-  commandArgs(trailingOnly = TRUE)
      if (any(grepl("^tests=", args))) {
            tests_to_match <-  strsplit(sub("^tests=", "", args[grepl("^tests=", args)]), ",")[[1]]
      }
      if (any(grepl("^references=", args))) {
            references_to_match <-  strsplit(sub("^references=", "", args[grepl("^references=", args)]), ",")[[1]]
      }
}
##
## ---- The prior model ---------------------------------------------------------------------------------------------------------------
##
model <-  cmdstanr::cmdstan_model(file.path(checks_dir, "prior_model_for_R_and_G_params.stan"),
                                  include_paths = file.path(local_pkg_dir, "inst", "stan_functions"))
names_of_the_params <-  c("prior_beta_mu_mean", "prior_beta_mu_SD", "prior_raw_scale_mu_mean", "prior_raw_scale_mu_SD",
                          "sigma")
##
## ---- Function: the reference trace of a test (the prior-only ordinal-bivariate fit under the flat prior or a prior set) -----------
##
fn_reference_trace <-  function(reference,
                                test) {

      if (reference == "flat") {
            file <-  file.path(checks_dir, paste0("reference_prior_draws_of_the_fused_ord_biv_fixed_model_at_alpha_one_",
                                                  test, "_16_chains_lower_bound_minus_15.rds"))
            reference_fit <-  readRDS(file)
            return(list(Se = reference_fit$Se,
                        Sp = reference_fit$Sp,
                        source = basename(file)))
      }
      reference_fit <-  readRDS(file.path(checks_dir, paste0("PO_check_prior_set_", reference, "_", test,
                                                             "_ord_biv_FIXEDthr_blank_local_HPC_n_burnin_500_n_iter_2500_",
                                                             "n_chains_8.rds")))
      ## (rows of Se and Sp: iterations within chains, so chain 1 is the first n_iter_of_the_reference_fit rows)
      list(Se = reference_fit$Se[1:n_iter_of_the_reference_fit, , drop = FALSE],
           Sp = reference_fit$Sp[1:n_iter_of_the_reference_fit, , drop = FALSE],
           source = basename(reference_fit$stan_model_file))

}
##
## ---- Each reference and test -------------------------------------------------------------------------------------------------------
##
matched_priors <-  list()
rows <-  list()
for (reference in references_to_match) {
      for (test in tests_to_match) {

            trace <-  fn_reference_trace(reference, test)
            rows_to_copy <-  unique(round(seq(1, nrow(trace$Se), length.out = n_draws_to_copy)))
            trace$Se <-  trace$Se[rows_to_copy, , drop = FALSE]
            trace$Sp <-  trace$Sp[rows_to_copy, , drop = FALSE]
            n_thr <-  ncol(trace$Se)
            tictoc::tic()
            fit <-  model$sample(data = list(n_thr = n_thr,
                                             N = nrow(trace$Se),
                                             probit_Sp_draws_to_copy = stats::qnorm(trace$Sp),
                                             probit_Se_draws_to_copy = stats::qnorm(trace$Se),
                                             prior_dirichlet_alpha = rep(1.0, n_thr + 1),
                                             softplus = 0,
                                             use_probit_link = 1),
                                 init = function() list(C_raw = matrix(-2.0, nrow(trace$Se), n_thr)),
                                 seed = seed,
                                 chains = n_chains,
                                 parallel_chains = n_chains,
                                 iter_warmup = n_burnin,
                                 iter_sampling = n_iter,
                                 refresh = 0,
                                 show_messages = FALSE)
            elapsed <-  tictoc::toc(quiet = TRUE)
            seconds <-  unname(elapsed$toc - elapsed$tic)
            draws_array <-  fit$draws(names_of_the_params, format = "draws_array")
            summary_of_the_fit <-  NicoStan::generate_summary_tibble(n_threads = 1,
                                                                      trace = aperm(draws_array, c(3, 1, 2)),
                                                                      param_names = names_of_the_params,
                                                                      n_to_compute = length(names_of_the_params),
                                                                      compute_nested_rhat = FALSE,
                                                                      n_chains = n_chains,
                                                                      n_superchains = 1)
            fn_median_of <-  function(param) stats::median(draws_array[, , param])
            matched <-  list(prior_beta_mu_mean = fn_median_of("prior_beta_mu_mean"),
                             prior_beta_mu_SD = fn_median_of("prior_beta_mu_SD"),
                             prior_raw_scale_mu_mean = fn_median_of("prior_raw_scale_mu_mean"),
                             prior_raw_scale_mu_SD = fn_median_of("prior_raw_scale_mu_SD"),
                             sigma = fn_median_of("sigma"),
                             reference = trace$source,
                             n_draws_of_the_reference = nrow(trace$Se))
            matched_priors[[paste0(if (reference == "flat") "flat_" else paste0("prior_set_", reference, "_"), test)]] <-  matched
            rows[[length(rows) + 1]] <-  data.frame(reference = reference,
                                                    test = test,
                                                    beta_mu_mean = round(matched$prior_beta_mu_mean, 3),
                                                    beta_mu_SD = round(matched$prior_beta_mu_SD, 3),
                                                    raw_scale_mu_mean = round(matched$prior_raw_scale_mu_mean, 3),
                                                    raw_scale_mu_SD = round(matched$prior_raw_scale_mu_SD, 3),
                                                    sigma = round(matched$sigma, 3),
                                                    min_ESS = round(min(summary_of_the_fit$n_eff, summary_of_the_fit$n_eff_tail,
                                                                        summary_of_the_fit$n_eff_sd)),
                                                    largest_R_hat = round(max(summary_of_the_fit$Rhat), 4),
                                                    n_divergent = sum(fit$sampler_diagnostics()[, , "divergent__"]),
                                                    seconds = round(seconds))
            options(width = 250)
            print(rows[[length(rows)]], row.names = FALSE)
            saveRDS(list(matched_priors = matched_priors[length(matched_priors)],
                         table_of_the_matched_priors = rows[[length(rows)]],
                         seed = seed),
                    file.path(checks_dir, paste0("R_and_G_priors_matched_to_the_ord_biv_reference_prior_",
                                                 if (reference == "flat") "flat" else paste0("set_", reference), "_", test,
                                                 ".rds")))

      }
}
##
## ---- Every reference and test in one file -------------------------------------------------------------------------------------------
##
files_of_each <-  list.files(checks_dir, pattern = "^R_and_G_priors_matched_to_the_ord_biv_reference_prior_(flat|set_1|set_2)_.*[.]rds$",
                             full.names = TRUE)
each <-  lapply(files_of_each, readRDS)
table_of_the_matched_priors <-  do.call(rbind, lapply(each, `[[`, "table_of_the_matched_priors"))
print(table_of_the_matched_priors, row.names = FALSE)
saveRDS(list(matched_priors = do.call(base::c, lapply(each, `[[`, "matched_priors")),
             table_of_the_matched_priors = table_of_the_matched_priors,
             seed = seed),
        file.path(checks_dir, "R_and_G_priors_matched_to_the_ord_biv_reference_priors_by_a_prior_model.rds"))
