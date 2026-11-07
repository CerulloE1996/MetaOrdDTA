##
## ---- Direct Se/Sp-prior check: gradients, then a NicoStan prior-only fit:
##
root <- c("/home/enzocerullo/Documents/Work/PhD_work", "/home/enzo/Documents/Work/PhD_work")
root <- root[dir.exists(root)][1]
pkg_dir <- file.path(root, "R_packages", "MetaOrdDTA")
sim_dir <- file.path(pkg_dir, "sim_study")
source(file.path(sim_dir, "ord_biv_FC_Se_Sp_prior_fns.R"))
source(file.path(sim_dir, "sim_study_priors_helper_fns.R"))
##
## ---- GAD-2 typed-prior check:
##
test <- "GAD_2"
chosen_cutoffs <- 1:6
Se_Sp_prior_type <- 4
##
## ---- Existing prior-only check settings:
##
## PO_check_ord_biv_FC_lit_alpha.R:32-34.
n_chains <- 4
n_burnin <- 125
n_iter <- 250
seed <- 1
gradient_tolerance <- 1e-8
finite_difference_step <- 1e-4
probabilities <- c(0.025, 0.5, 0.975)
gap_limit <- 0.05
comparison_reference_file <- NULL
comparison_reference <- if (is.null(comparison_reference_file)) NULL else {

      readRDS(comparison_reference_file)

}
##
## ---- Original data preparation and seeded starts; no native model compilation:
##
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
x <- fn_data_set_of_a_test(file.path(root, "R_packages", "MetaOrdDTA_outputs", "sim_study"), test)
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(debugging = FALSE,
                                                 x = x,
                                                 indicator_index_test_in_study = NULL,
                                                 intercept_only = TRUE,
                                                 cov_data = NULL,
                                                 n_chains = n_chains,
                                                 network = FALSE,
                                                 compound_symmetry = NULL,
                                                 hetero_sigma = NULL,
                                                 prior_only = FALSE,
                                                 model_parameterisation = "ord_biv",
                                                 random_thresholds = FALSE,
                                                 Dirichlet_random_effects_type = "fixed",
                                                 softplus = FALSE,
                                                 advanced_compile = FALSE,
                                                 force_recompile = FALSE,
                                                 quiet = TRUE,
                                                 compile = FALSE,
                                                 custom_file_name = NULL)
priors <- model_prep_obj$priors
priors$prior_dirichlet_alpha <- list(rep(1, model_prep_obj$internal_obj$outs_data$n_cat),
                                    rep(1, model_prep_obj$internal_obj$outs_data$n_cat))
priors$prior_beta_mu_mean[, 1] <- 0
priors$prior_beta_mu_SD[, 1] <- 1
priors$prior_beta_SD_SD <- c(0.5, 0.5)
priors$use_probit_link <- 1
priors$K_fold_CV_indicator <- rep(0, nrow(x[[1]]))
outs_data <- model_prep_obj$internal_obj$outs_data
Stan_data <- outs_data$stan_data_list
Stan_data$softplus <- FALSE
Stan_data$n_studies <- outs_data$n_studies
Stan_data$n_thr <- outs_data$n_thr
Stan_data$n_cat <- outs_data$n_cat
Stan_data <- c(Stan_data, priors)
Stan_data <- MetaOrdDTA:::R_fn_remove_duplicates_from_list(Stan_data)
Stan_data$use_Se_Sp_prior <- 1
Stan_data$Se_Sp_prior_type <- Se_Sp_prior_type
Stan_data$n_chosen_cutoffs_max <- length(chosen_cutoffs)
Stan_data$n_chosen_cutoffs <- rep(length(chosen_cutoffs), 2)
Stan_data$chosen_cutoffs <- rbind(chosen_cutoffs, chosen_cutoffs)
Stan_data$prior_Se_Sp_dirichlet_A <- priors$prior_dirichlet_alpha
Stan_data$prior_Se_Sp_beta <- list(matrix(0, length(chosen_cutoffs), 2), matrix(0, length(chosen_cutoffs), 2))
Stan_data$prior_Se_Sp_link_mean_SD <- list(cbind(rep(0, length(chosen_cutoffs)),
                                                rep(1, length(chosen_cutoffs))),
                                           cbind(rep(0, length(chosen_cutoffs)),
                                                 rep(1, length(chosen_cutoffs))))
init_lists_per_chain <- MetaOrdDTA::resize_init_list(init_lists_per_chain = model_prep_obj$init_lists_per_chain,
                                                     n_chains_new = n_chains)
##
## ---- Sampler controls copied literally from the existing sim runner:
##
file_of_the_runner <- if (Sys.info()[["user"]] == "enzocerullo") {
      file.path(sim_dir, "run_MULTIPLE_MetaOrd_simulation.R")
} else {
      "/tmp/run_MULTIPLE_MetaOrd_simulation_of_the_local_HPC.R"
}
source(file_of_the_runner)
sample_call <- fn_sample_call_Se_Sp_prior(body(run_multiple_sims), "NicoStan_args")
stopifnot(!is.null(sample_call))
nuisance_params <- MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                                    name_of_the_model = "DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior",
                                    model_parameterisation = "ord_biv",
                                    random_thresholds = FALSE)
interest_only_params <- MetaOrdDTA::fn_interest_only_params_of_a_MetaOrdDTA_model(
                                       name_of_the_model = "DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior",
                                       model_parameterisation = "ord_biv",
                                       random_thresholds = FALSE)
sampler_args <- eval(sample_call$NicoStan_args)
Stan_file <- file.path(pkg_dir, "inst", "stan_models", "stan_models_MA",
                       "DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior.stan")
output_prefix <- file.path(sim_dir, "prior_derivation_tools_checks",
                           paste0("PO_ord_biv_FC_Se_Sp_prior_typed_type_", Se_Sp_prior_type, "_", test))
##
## ---- Gradients are a compulsory gate before sampling:
##
results <- fn_run_Se_Sp_prior_check(Stan_data = Stan_data,
                                     Stan_file = Stan_file,
                                     init_lists_per_chain = init_lists_per_chain,
                                     sampler_args = sampler_args,
                                     n_chains = n_chains,
                                     n_burnin = n_burnin,
                                     n_iter = n_iter,
                                     seed = seed,
                                     output_prefix = output_prefix,
                                     gradient_tolerance = gradient_tolerance,
                                     finite_difference_step = finite_difference_step,
                                     probabilities = probabilities,
                                     gap_limit = gap_limit,
                                     comparison_reference = comparison_reference)





















