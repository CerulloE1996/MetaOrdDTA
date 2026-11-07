##
## ==== check_lp_grad_NMA_ord_biv_FC_Se_Sp_fused_vs_Stan_math.R ================================================
##
## ---- DTA_NMA_Nyaga_ord_biv_FIXEDthr_Se_Sp_prior_NicoStan_fused_AVX.stan against
##      DTA_NMA_Nyaga_ord_biv_FIXEDthr_Se_Sp_prior.stan (Stan math): the log density and the gradient at
##      6 random unconstrained points, in 10 settings of the Se/Sp prior (off; off with the likelihood on; each
##      of the 4 prior types at each anchor), on NicoStan's worked example 4 (4 ordinal screening tests in a
##      network). The parameters blocks of the two files are identical, so the same points feed both libraries
##      and the gradients are compared column for column. Gate: log density within 1e-8 relative, gradient
##      within 1e-8 of its largest element.
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
      checks_dir <-  file.path(sim_study_dir, "prior_derivation_tools_checks")
      models_folder <-  file.path(local_pkg_dir, "inst", "stan_models", "stan_models_NMA")
      functions_folder <-  file.path(local_pkg_dir, "inst", "stan_functions")
      NicoSIMD_user_header <-  file.path(root_of_the_work_folder, "R_packages", "NicoStan", "inst", "NicoSIMD",
                                         "inst", "include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp")
      file_of_the_simulated_data <-  paste0(root_of_the_work_folder,
                                            "/R_packages/NicoStan/inst/examples/data/",
                                            "simulated_network_meta_analysis_of_four_ordinal_screening_",
                                            "tests_with_missing_thresholds.rds")
      model <-  "DTA_NMA_Nyaga_ord_biv_FIXEDthr_Se_Sp_prior"
      n_points <-  6
      SD_of_the_noise_of_the_points <-  0.3
      tolerance_of_the_rule <-  1e-8
      seed <-  1
      CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
      machine <-  if (CPU_has_AVX512) "local_HPC" else "laptop"
      AVX_label <-  if (CPU_has_AVX512) "AVX512" else "AVX2"
      compiler_flags <-  if (CPU_has_AVX512) {
            "-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq"
      } else {
            "-O3 -march=native -mtune=native -mfma -mavx2"
      }
      libraries_folder_of_the_Stan_math_checks <-  file.path(sim_study_dir,
                                                             paste0("libraries_of_the_NMA_ord_biv_FC_Se_Sp_",
                                                                    "checks_", machine))
      libraries_folder_of_the_fused_versions <-  file.path(sim_study_dir,
                                                           paste0("libraries_of_the_NicoStan_fused_AVX_",
                                                                  "versions_", machine))
}
##
options(warning.length = 8170)
loadNamespace("NicoStan")
fn_log_density_at_points <-  get("fn_log_density_and_gradient_of_compiled_Stan_model_at_points",
                                 envir = asNamespace("NicoStan"))
dir.create(libraries_folder_of_the_Stan_math_checks, showWarnings = FALSE)
dir.create(libraries_folder_of_the_fused_versions, showWarnings = FALSE)
##
## ---- The data of NicoStan's worked example 4 (as check_NicoStan_option_NMA_model_vs_cmdstanr.R):
##
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
simulated_NMA_data <-  readRDS(file_of_the_simulated_data)
x <-  lapply(simulated_NMA_data$x, unname)
indicator_index_test_in_study <-  simulated_NMA_data$indicator_index_test_in_study
model_prep_obj <-  MetaOrdDTA::MetaOrd_model$new(x = x,
                                                 indicator_index_test_in_study = indicator_index_test_in_study,
                                                 intercept_only = TRUE,
                                                 network = TRUE,
                                                 model_parameterisation = "ord_biv",
                                                 random_thresholds = FALSE,
                                                 Dirichlet_random_effects_type = "fixed",
                                                 n_chains = 4,
                                                 compile = FALSE,
                                                 quiet = TRUE)
priors <-  model_prep_obj$priors
outs_data <-  model_prep_obj$internal_obj$outs_data
Stan_data_of_the_cell <-  outs_data$stan_data_list
Stan_data_of_the_cell$softplus <-  0
Stan_data_of_the_cell <-  c(Stan_data_of_the_cell, priors)
Stan_data_of_the_cell <-  MetaOrdDTA:::R_fn_remove_duplicates_from_list(Stan_data_of_the_cell)
fields_the_model_needs <-  c("n_studies", "n_index_tests", "n_studies_per_test", "n_thr", "n_cat", "n_thr_max",
                             "n_cat_max", "n_obs_cutpoints", "indicator_index_test_in_study", "x",
                             "cutpoint_index", "n_covariates_nd", "n_covariates_d", "X_nd", "X_d",
                             "baseline_case_nd", "baseline_case_d", "prior_beta_mu_mean", "prior_beta_mu_SD",
                             "prior_beta_tau_SD", "prior_beta_sigma_SD", "beta_corr_lb", "beta_corr_ub",
                             "prior_beta_corr_LKJ", "prior_dirichlet_alpha", "softplus", "use_probit_link",
                             "hetero_sigma", "compound_symmetry", "K_fold_CV_indicator")
stopifnot(all(fields_the_model_needs %in% names(Stan_data_of_the_cell)))
##
## ---- One R shape for every field:
##
fn_coerce_the_shapes_of_the_NMA_data <- function(Stan_data_list) {

      for (name in c("baseline_case_nd", "baseline_case_d")) {
            if (is.list(Stan_data_list[[name]])) {
                  Stan_data_list[[name]] <-  do.call(rbind, Stan_data_list[[name]])
            }
      }
      if (!is.list(Stan_data_list$prior_dirichlet_alpha)) {
            alpha_arr <-  Stan_data_list$prior_dirichlet_alpha
            Stan_data_list$prior_dirichlet_alpha <-  list(alpha_arr[1, , , drop = FALSE][1, , ],
                                                          alpha_arr[2, , , drop = FALSE][1, , ])
      }
      Stan_data_list

}
Stan_data_of_the_cell <-  fn_coerce_the_shapes_of_the_NMA_data(Stan_data_of_the_cell)
n_tests <-  Stan_data_of_the_cell$n_index_tests
n_studies <-  Stan_data_of_the_cell$n_studies
n_thr <-  Stan_data_of_the_cell$n_thr
n_covariates_max <-  max(c(Stan_data_of_the_cell$n_covariates_nd, Stan_data_of_the_cell$n_covariates_d))
##
## ---- The settings of the Se/Sp prior (off; each type at each anchor; one cell with nothing chosen):
##
fn_segment_sums <- function(alpha,
                            chosen_cutoffs) {

      bounds <-  c(0, chosen_cutoffs, length(alpha))
      sapply(seq_len(length(bounds) - 1), function(j) sum(alpha[(bounds[j] + 1):bounds[j + 1]]))

}
fn_chosen_of_a_test <- function(n_thr_of_the_test) {

      if (n_thr_of_the_test >= 4) return(c(2, 4))
      if (n_thr_of_the_test >= 2) return(c(1, n_thr_of_the_test))
      return(1)

}
fn_Se_Sp_fields_NMA <- function(use_Se_Sp_prior,
                                Se_Sp_prior_type,
                                chosen_of_each_test_class,
                                A_of_each_test_class,
                                beta_ab_of_each_test_class,
                                link_mean_SD_of_each_test_class) {

      n_chosen_cutoffs_max <-  1
      for (t in seq_along(chosen_of_each_test_class)) {
            for (c in 1:2) {
                  n_chosen_cutoffs_max <-  max(n_chosen_cutoffs_max, length(chosen_of_each_test_class[[t]][[c]]))
            }
      }
      n_chosen_cutoffs <-  matrix(0L, n_tests, 2)
      for (t in seq_len(n_tests)) {
            for (c in 1:2) {
                  n_chosen_cutoffs[t, c] <-  length(chosen_of_each_test_class[[t]][[c]])
            }
      }
      chosen_cutoffs <-  array(1L, dim = c(n_tests, 2, n_chosen_cutoffs_max))
      A_arr <-  array(1, dim = c(n_tests, 2, n_chosen_cutoffs_max + 1))
      beta_ab_arr <-  array(1, dim = c(n_tests, 2, n_chosen_cutoffs_max, 2))
      link_arr <-  array(0, dim = c(n_tests, 2, n_chosen_cutoffs_max, 2))
      for (t in seq_len(n_tests)) {
            for (c in 1:2) {
                  m <-  length(chosen_of_each_test_class[[t]][[c]])
                  if (m == 0) next
                  chosen_cutoffs[t, c, 1:m] <-  as.integer(chosen_of_each_test_class[[t]][[c]])
                  A_arr[t, c, 1:(m + 1)] <-  A_of_each_test_class[[t]][[c]]
                  beta_ab_arr[t, c, 1:m, ] <-  beta_ab_of_each_test_class[[t]][[c]]
                  link_arr[t, c, 1:m, 2] <-  1
            }
      }
      list(use_Se_Sp_prior = use_Se_Sp_prior,
           Se_Sp_prior_type = Se_Sp_prior_type,
           n_chosen_cutoffs_max = n_chosen_cutoffs_max,
           n_chosen_cutoffs = n_chosen_cutoffs,
           chosen_cutoffs = chosen_cutoffs,
           prior_Se_Sp_dirichlet_A = A_arr,
           prior_Se_Sp_beta = beta_ab_arr,
           prior_Se_Sp_link_mean_SD = link_arr)

}
fn_field_lists_of_a_set_of_cutoffs <- function(chosen_of_each_test_class,
                                               alpha_of_each_test_class) {

      A_of_each_test_class <-  list()
      beta_ab_of_each_test_class <-  list()
      link_mean_SD_of_each_test_class <-  list()
      for (t in seq_len(n_tests)) {
            A_of_each_test_class[[t]] <-  list()
            beta_ab_of_each_test_class[[t]] <-  list()
            link_mean_SD_of_each_test_class[[t]] <-  list()
            for (c in 1:2) {
                  m <-  length(chosen_of_each_test_class[[t]][[c]])
                  if (m == 0) {
                        A_of_each_test_class[[t]][[c]] <-  c(1, 1)
                        beta_ab_of_each_test_class[[t]][[c]] <-  matrix(1, 1, 2)
                        link_mean_SD_of_each_test_class[[t]][[c]] <-  matrix(0, 1, 2)
                  } else {
                        A_of_each_test_class[[t]][[c]] <-
                              fn_segment_sums(alpha_of_each_test_class[[t]][[c]],
                                              chosen_of_each_test_class[[t]][[c]])
                        beta_ab_of_each_test_class[[t]][[c]] <-
                              if (m >= 2) matrix(c(4, 6, 6, 3), nrow = m, ncol = 2, byrow = TRUE) else
                                          matrix(c(4, 6), nrow = m)
                        link_mean_SD_of_each_test_class[[t]][[c]] <-  matrix(0, m, 2)
                  }
            }
      }
      list(A = A_of_each_test_class, beta_ab = beta_ab_of_each_test_class,
           link_mean_SD = link_mean_SD_of_each_test_class)

}
alpha_of_each_test_class <-  lapply(seq_len(n_tests), function(t) {

      list(Stan_data_of_the_cell$prior_dirichlet_alpha[[1]][t, 1:(n_thr[t] + 1)],
           Stan_data_of_the_cell$prior_dirichlet_alpha[[2]][t, 1:(n_thr[t] + 1)])

})
chosen_of_the_off_setting <-  lapply(seq_len(n_tests), function(t) list(integer(0), integer(0)))
chosen_of_the_on_settings <-  lapply(seq_len(n_tests), function(t) {

      list(fn_chosen_of_a_test(n_thr[t]), fn_chosen_of_a_test(n_thr[t]))

})
chosen_of_the_on_settings[[1]][[2]] <-  integer(0)   ## the else branch of the new file under use = 1
off_field_lists <-  fn_field_lists_of_a_set_of_cutoffs(chosen_of_the_off_setting, alpha_of_each_test_class)
on_field_lists <-  fn_field_lists_of_a_set_of_cutoffs(chosen_of_the_on_settings, alpha_of_each_test_class)
settings <-  list()
settings[["off"]] <-  fn_Se_Sp_fields_NMA(0, 1, chosen_of_the_off_setting, off_field_lists$A,
                                          off_field_lists$beta_ab, off_field_lists$link_mean_SD)
for (prior_type in 1:4) {

      {

            settings[[paste0("type_", prior_type)]] <-
                  fn_Se_Sp_fields_NMA(1, prior_type, chosen_of_the_on_settings, on_field_lists$A,
                                      on_field_lists$beta_ab, on_field_lists$link_mean_SD)

      }

}
K_fold_of_each_setting <-  c(off = 0, off_likelihood_on = 1, type_1 = 0, type_1 = 0,
                             type_2 = 0, type_2 = 0, type_3 = 0, type_3 = 0,
                             type_4 = 0, type_4 = 0)
settings[["off_likelihood_on"]] <-  settings[["off"]]
##
fn_Stan_data_of_a_setting <- function(name_of_the_setting,
                                      Stan_data_of_the_cell,
                                      settings,
                                      K_fold_of_each_setting) {

      Stan_data_list <-  c(Stan_data_of_the_cell, settings[[name_of_the_setting]])
      Stan_data_list <-  MetaOrdDTA:::R_fn_remove_duplicates_from_list(Stan_data_list)
      Stan_data_list$K_fold_CV_indicator <-  rep(K_fold_of_each_setting[[name_of_the_setting]],
                                                 Stan_data_list$n_studies)
      fn_coerce_the_shapes_of_the_NMA_data(Stan_data_list)

}
## (a field of one value declared as a container is written as an array of length 1)
lines_of_the_model <-  readLines(file.path(models_folder, paste0(model, ".stan")), warn = FALSE)
lines_of_the_data_block <-  sub("//.*$", "",
                                lines_of_the_model[grep("^data *\\{", lines_of_the_model)[1]:
                                                   grep("^transformed data", lines_of_the_model)[1]])
fn_containers_of_one_value_as_arrays <- function(Stan_data_list) {

      for (name in names(Stan_data_list)) {
            declared_as_a_container <-  any(grepl(paste0("^ *(vector|row_vector|array)[^;]*[] >]", name, " *;"),
                                                  lines_of_the_data_block))
            if (declared_as_a_container && (length(Stan_data_list[[name]]) == 1) &&
                is.null(dim(Stan_data_list[[name]]))) {
                  Stan_data_list[[name]] <-  array(Stan_data_list[[name]], dim = 1)
            }
      }
      Stan_data_list

}
fn_JSON_text_of <- function(Stan_data_list) {

      Stan_data_list <-  fn_containers_of_one_value_as_arrays(Stan_data_list)
      JSON_file <-  tempfile(fileext = ".json")
      cmdstanr::write_stan_json(data = Stan_data_list, file = JSON_file)
      paste(readLines(JSON_file, warn = FALSE), collapse = "\n")

}
Stan_data_of_the_off_setting <-  fn_Stan_data_of_a_setting("off", Stan_data_of_the_cell, settings,
                                                           K_fold_of_each_setting)
##
## ---- The two libraries (built unless built already):
##
fn_build_unless_built <- function(source_file,
                                  library_name,
                                  Stan_data_list) {

      Stan_model_file_path <-  file.path(libraries_folder_of_the_Stan_math_checks, paste0(library_name, ".stan"))
      so_file <-  sub("\\.stan$", "_model.so", Stan_model_file_path)
      if (file.exists(so_file) && file.exists(Stan_model_file_path) &&
          identical(readLines(Stan_model_file_path, warn = FALSE), readLines(source_file, warn = FALSE)) &&
          (file.info(so_file)$mtime > file.info(Stan_model_file_path)$mtime)) {
            message(paste0("\033[36mAlready built: ", so_file, "\033[0m"))
            return(so_file)
      }
      stopifnot(file.copy(source_file, Stan_model_file_path, overwrite = TRUE, copy.date = FALSE))
      Sys.setFileTime(Stan_model_file_path, Sys.time())
      make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                      paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
      stanc_args <-  c(paste0("--include-paths=", functions_folder))
      start <-  Sys.time()
      NicoStan:::initialise_model(Model_type = "Stan",
                                  sample_nuisance = FALSE,
                                  Stan_data_list = Stan_data_list,
                                  Stan_model_file_path = Stan_model_file_path,
                                  Stan_cpp_user_header = NULL,
                                  stanc_args = stanc_args,
                                  make_args = make_args,
                                  compile = TRUE)
      stopifnot(file.exists(so_file))
      message(paste0("\033[32mBuilt ", so_file, " in ",
                     formatC(as.numeric(difftime(Sys.time(), start, units = "secs")), format = "f", digits = 1),
                     " s\033[0m"))
      return(so_file)

}
Stan_data_list_of_the_build <-  fn_containers_of_one_value_as_arrays(Stan_data_of_the_off_setting)
fn_build_fused_unless_built <- function(source_file,
                                        library_name,
                                        user_header,
                                        AVX2_only) {

        Stan_model_file_path <-  file.path(libraries_folder_of_the_fused_versions, paste0(library_name, ".stan"))
        so_file <-  sub("\\.stan$", "_model.so", Stan_model_file_path)
        if (file.exists(so_file) && file.exists(Stan_model_file_path) &&
            identical(readLines(Stan_model_file_path, warn = FALSE), readLines(source_file, warn = FALSE)) &&
            (file.info(so_file)$mtime > file.info(Stan_model_file_path)$mtime)) {
              message(paste0("\033[36mAlready built: ", so_file, "\033[0m"))
              return(so_file)
        }
        stopifnot(file.copy(source_file, Stan_model_file_path, overwrite = TRUE, copy.date = FALSE))
        Sys.setFileTime(Stan_model_file_path, Sys.time())
        make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                        paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
        if (AVX2_only) make_args <-  c(make_args, "CPPFLAGS_OPTIM=-DBAYESMVP_FORCE_AVX2")
        stanc_args <-  c(paste0("--include-paths=", functions_folder))
        if (!is.null(user_header)) stanc_args <-  c(stanc_args, "--allow-undefined")
        start <-  Sys.time()
        NicoStan:::initialise_model( Model_type            = "Stan",
                                     sample_nuisance       = FALSE,
                                     Stan_data_list        = Stan_data_list_of_the_build,
                                     Stan_model_file_path  = Stan_model_file_path,
                                     Stan_cpp_user_header  = user_header,
                                     stanc_args            = stanc_args,
                                     make_args             = make_args,
                                     compile               = TRUE)
        stopifnot(file.exists(so_file))
        message(paste0("\033[32mBuilt ", so_file, " in ",
                       formatC(as.numeric(difftime(Sys.time(), start, units = "secs")), format = "f", digits = 1),
                       " s\033[0m"))
        return(so_file)

}
so_file_of_the_model <-  fn_build_unless_built(source_file = file.path(models_folder, paste0(model, ".stan")),
                                               library_name = paste0(model, "_Stan_math"),
                                               Stan_data_list = Stan_data_of_the_off_setting)
so_file_of_the_fused_version <-  fn_build_fused_unless_built(
                                       source_file   = file.path(models_folder,
                                                                 paste0(model, "_NicoStan_fused_AVX.stan")),
                                       library_name  = paste0(model, "_NicoStan_fused_", AVX_label),
                                       user_header   = NicoSIMD_user_header,
                                       AVX2_only     = !CPU_has_AVX512)
##
## ---- The points (centre 0; log(0.3) for the positive SD parameters):
##
fn_points <- function(parameter_names) {

      centre <-  rep(0, length(parameter_names))
      centre[grepl("beta_sigma|beta_tau_raw|_SD", parameter_names)] <-  log(0.3)
      lapply(seq_len(n_points), function(point_index) {
            centre + stats::rnorm(length(parameter_names), mean = 0, sd = SD_of_the_noise_of_the_points)
      })

}
fn_evaluate_in_a_fresh_R_process <- function(so_file,
                                             JSON_text,
                                             points) {

      callr::r(func = fn_log_density_at_points,
               args = list(model_so_file = so_file,
                           JSON_text_of_Stan_data = JSON_text,
                           points = points,
                           propto = FALSE,
                           jacobian = TRUE,
                           seed = 1))

}
JSON_text_of_the_build <-  fn_JSON_text_of(Stan_data_of_the_off_setting)
names_of_the_parameters <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_model, JSON_text_of_the_build,
                                                             list())$unconstrained_parameter_names
names_of_the_fused_version <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_fused_version,
                                                                JSON_text_of_the_build,
                                                                list())$unconstrained_parameter_names
stopifnot(identical(names_of_the_parameters, names_of_the_fused_version))
message(paste0("\033[36m", length(names_of_the_parameters), " unconstrained parameters; first: ",
               paste(head(names_of_the_parameters, 8), collapse = ", "), "; ",
               length(settings), " settings\033[0m"))
##
## ---- The log density and gradient of both at the same points, setting by setting:
##
rows <-  NULL
for (name_of_the_setting in names(settings)) {

      set.seed(seed + match(name_of_the_setting, names(settings)))
      Stan_data_of_the_setting <-  fn_Stan_data_of_a_setting(name_of_the_setting, Stan_data_of_the_cell,
                                                             settings, K_fold_of_each_setting)
      JSON_text <-  fn_JSON_text_of(Stan_data_of_the_setting)
      points <-  fn_points(names_of_the_parameters)
      start <-  Sys.time()
      of_the_model <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_model, JSON_text, points)
      of_the_fused_version <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_fused_version, JSON_text,
                                                                points)
      every_point_finite <-  all(is.finite(of_the_model$log_density)) &&
                             all(is.finite(of_the_fused_version$log_density)) &&
                             all(is.finite(of_the_model$gradient)) &&
                             all(is.finite(of_the_fused_version$gradient))
      log_density_relative <-  max(abs(of_the_fused_version$log_density - of_the_model$log_density) /
                                   abs(of_the_model$log_density))
      gradient_relative <-  max(apply(abs(of_the_fused_version$gradient - of_the_model$gradient), 1, max) /
                                apply(abs(of_the_model$gradient), 1, max))
      row <-  data.frame(setting = name_of_the_setting,
                         n_points = n_points,
                         every_point_finite = every_point_finite,
                         largest_abs_log_density = max(abs(of_the_model$log_density)),
                         log_density_relative = log_density_relative,
                         gradient_relative_of_largest_element = gradient_relative,
                         passed_1e_8 = every_point_finite && log_density_relative <= tolerance_of_the_rule &&
                                       gradient_relative <= tolerance_of_the_rule,
                         seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
      rows <-  rbind(rows, row)
      message(paste0(if (row$passed_1e_8) "\033[32m" else "\033[31m",
                     name_of_the_setting, ": log density ",
                     formatC(log_density_relative, format = "e", digits = 2),
                     " relative, gradient ", formatC(gradient_relative, format = "e", digits = 2),
                     " of its largest element; ", if (row$passed_1e_8) "PASS (1e-8)" else "FAIL",
                     "; ", formatC(row$seconds, format = "f", digits = 1), " s\033[0m"))

}
print(rows, digits = 4)
file_of_the_check <-  file.path(checks_dir, paste0("check_lp_grad_NMA_ord_biv_FC_Se_Sp_fused_vs_Stan_math_",
                                                   machine, ".csv"))
utils::write.csv(rows, file_of_the_check, row.names = FALSE)
message(paste0(if (all(rows$passed_1e_8)) "\033[32mPASS" else "\033[31mFAIL",
               ": ", sum(rows$passed_1e_8), " of ", nrow(rows), " settings within 1e-8; saved ",
               file_of_the_check, "\033[0m"))
stopifnot(all(rows$passed_1e_8))






















