##
## ==== check_lp_grad_ord_biv_FC_Se_Sp_fused_vs_Stan_math.R ===========================================================
##
## ---- DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior_NicoStan_fused_AVX.stan against DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior.stan
##      (Stan math): the log density (every constant) and the gradient at 6 random unconstrained points, in 12
##      settings of the Se/Sp prior (off; each of the 4 prior types at each anchor with the cutoffs 2 and 5 chosen;
##      prior type 1 at the location with every cutoff chosen, likelihood on and off; prior type 4 at the location
##      with every cutoff chosen, mean 0 and SD 1, the default); the GAD-2 data set of the sim study and
##      the priors of its literature cell. Gate: log density within 1e-8 (relative), gradient within 1e-8 of its
##      largest element; the fused AVX gate of 0.001 reported beside it.
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
      checks_dir <-  file.path(sim_study_dir, "prior_derivation_tools_checks")
      alpha_sets_dir <-  file.path(sim_study_dir, "prior_only_checks_of_the_induced_Dirichlet_alpha_sets")
      outputs_dir <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study")
      models_folder <-  file.path(local_pkg_dir, "inst", "stan_models", "stan_models_MA")
      functions_folder <-  file.path(local_pkg_dir, "inst", "stan_functions")
      NicoSIMD_user_header <-  file.path(root_of_the_work_folder, "R_packages", "NicoStan", "inst", "NicoSIMD",
                                         "inst", "include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp")
      model <-  "DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior"
      test <-  "GAD_2"
      n_points <-  6
      SD_of_the_noise_of_the_points <-  0.3
      tolerance_of_the_rule <-  1e-8
      tolerance_of_the_fused_gate <-  0.001
      seed <-  1
      CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
      machine <-  if (CPU_has_AVX512) "local_HPC" else "laptop"
      AVX_label <-  if (CPU_has_AVX512) "AVX512" else "AVX2"
      compiler_flags <-  if (CPU_has_AVX512) {
            "-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq"
      } else {
            "-O3 -march=native -mtune=native -mfma -mavx2"
      }
      libraries_folder <-  file.path(sim_study_dir, paste0("libraries_of_the_NicoStan_fused_AVX_versions_",
                                                           machine))
}
##
options(warning.length = 8170)
source(file.path(sim_study_dir, "sim_study_priors_helper_fns.R"))
loadNamespace("NicoStan")
fn_log_density_at_points <-  get("fn_log_density_and_gradient_of_compiled_Stan_model_at_points",
                                 envir = asNamespace("NicoStan"))
dir.create(libraries_folder, showWarnings = FALSE)
##
## ---- The data and priors of the literature cell of GAD-2 (as ord_biv_FC_Se_Sp_prior.R), likelihood on:
##
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
x <-  fn_data_set_of_a_test(outputs_dir, test)
model_prep_obj <-  MetaOrdDTA::MetaOrd_model$new(debugging = FALSE,
                                                 x = x,
                                                 indicator_index_test_in_study = NULL,
                                                 intercept_only = TRUE,
                                                 cov_data = NULL,
                                                 n_chains = 4,
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
derived <-  readRDS(file.path(alpha_sets_dir,
                              paste0("informative_from_the_literature_by_Stan_sampling_", test, ".rds")))
priors <-  model_prep_obj$priors
priors$prior_dirichlet_alpha <-  derived$alpha
priors$prior_beta_mu_mean[, 1] <-  derived$location_mean
priors$prior_beta_mu_SD[, 1] <-  derived$location_SD
priors$prior_beta_SD_SD <-  c(0.5, 0.5)
priors$use_probit_link <-  1
priors$K_fold_CV_indicator <-  rep(1, nrow(x[[1]]))
outs_data <-  model_prep_obj$internal_obj$outs_data
Stan_data_of_the_cell <-  outs_data$stan_data_list
Stan_data_of_the_cell$softplus <-  0
Stan_data_of_the_cell$n_studies <-  outs_data$n_studies
Stan_data_of_the_cell$n_thr <-  outs_data$n_thr
Stan_data_of_the_cell$n_cat <-  outs_data$n_cat
Stan_data_of_the_cell <-  c(Stan_data_of_the_cell, priors)
Stan_data_of_the_cell <-  MetaOrdDTA:::R_fn_remove_duplicates_from_list(Stan_data_of_the_cell)
n_thr <-  Stan_data_of_the_cell$n_thr
##
## ---- The 12 settings of the Se/Sp prior:
##
fn_Se_Sp_fields <-  function(use_Se_Sp_prior,
                             Se_Sp_prior_type,
                             chosen_cutoffs,
                             A_of_each_class,
                             beta_ab_of_each_class,
                             link_mean_SD_of_each_class) {

      m <-  length(chosen_cutoffs)
      list(use_Se_Sp_prior = use_Se_Sp_prior,
           Se_Sp_prior_type = Se_Sp_prior_type,
           n_chosen_max = m,
           n_chosen = c(m, m),
           chosen_cutoffs = rbind(chosen_cutoffs, chosen_cutoffs),
           prior_Se_Sp_dirichlet_A = A_of_each_class,
           prior_Se_Sp_beta = beta_ab_of_each_class,
           prior_Se_Sp_link_mean_SD = link_mean_SD_of_each_class)

}
fn_segment_sums <-  function(alpha,
                             chosen_cutoffs) {

      bounds <-  c(0, chosen_cutoffs, length(alpha))
      sapply(seq_len(length(bounds) - 1), function(j) sum(alpha[(bounds[j] + 1):bounds[j + 1]]))

}
chosen_two <-  c(2, 5)
A_of_two <-  list(fn_segment_sums(derived$alpha[[1]], chosen_two), fn_segment_sums(derived$alpha[[2]], chosen_two))
beta_ab_of_two <-  list(cbind(c(4, 6), c(6, 3)), cbind(c(4, 6), c(6, 3)))
link_mean_SD_of_two <-  list(cbind(c(0.3, -0.2), c(0.8, 1.2)), cbind(c(0.4, -0.5), c(0.8, 1.2)))
link_mean_SD_of_one <-  list(matrix(c(0, 1), 1, 2), matrix(c(0, 1), 1, 2))
link_mean_SD_of_all <-  list(cbind(rep(0, n_thr), rep(1, n_thr)), cbind(rep(0, n_thr), rep(1, n_thr)))
settings <-  list()
settings[["off"]] <-  fn_Se_Sp_fields(0, 1, 1, list(c(1, 1), c(1, 1)), list(matrix(1, 1, 2), matrix(1, 1, 2)),
                                      link_mean_SD_of_one)
for (prior_type in 1:4) {

      for (anchor in 0:1) {

            settings[[paste0("type_", prior_type, "_anchor_", anchor, "_cutoffs_2_5")]] <-
                  fn_Se_Sp_fields(1, anchor, prior_type, chosen_two, A_of_two, beta_ab_of_two,
                                  link_mean_SD_of_two)

      }

}
settings[["type_1_all_cutoffs"]] <-  fn_Se_Sp_fields(1, 1, 1:n_thr, derived$alpha,
                                                              list(matrix(1, n_thr, 2), matrix(1, n_thr, 2)),
                                                              link_mean_SD_of_all)
settings[["type_1_all_cutoffs_prior_only"]] <-  settings[["type_1_all_cutoffs"]]
settings[["type_4_all_cutoffs_default"]] <-  fn_Se_Sp_fields(1, 4, 1:n_thr, derived$alpha,
                                                                      list(matrix(1, n_thr, 2),
                                                                           matrix(1, n_thr, 2)),
                                                                      link_mean_SD_of_all)
##
fn_Stan_data_of_a_setting <-  function(name_of_the_setting) {

      Stan_data_list <-  c(Stan_data_of_the_cell, settings[[name_of_the_setting]])
      if (grepl("prior_only", name_of_the_setting)) {
            Stan_data_list$K_fold_CV_indicator <-  rep(0, Stan_data_list$n_studies)
      }
      Stan_data_list

}
## (a field of one value declared as a container is written as an array of length 1)
lines_of_the_model <-  readLines(file.path(models_folder, paste0(model, ".stan")), warn = FALSE)
lines_of_the_data_block <-  sub("//.*$", "",
                                lines_of_the_model[grep("^data *\\{", lines_of_the_model)[1]:
                                                   grep("^transformed data", lines_of_the_model)[1]])
fn_containers_of_one_value_as_arrays <-  function(Stan_data_list) {

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
fn_JSON_text_of <-  function(Stan_data_list) {

      Stan_data_list <-  fn_containers_of_one_value_as_arrays(Stan_data_list)
      JSON_file <-  tempfile(fileext = ".json")
      cmdstanr::write_stan_json(data = Stan_data_list, file = JSON_file)
      paste(readLines(JSON_file, warn = FALSE), collapse = "\n")

}
Stan_data_list <-  fn_containers_of_one_value_as_arrays(fn_Stan_data_of_a_setting("type_1_all_cutoffs"))
##
## ---- The two libraries (built unless built already):
##
fn_build_unless_built <-  function(source_file,
                                   library_name,
                                   user_header,
                                   AVX2_only) {

        Stan_model_file_path <-  file.path(libraries_folder, paste0(library_name, ".stan"))
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
                                     Stan_data_list        = Stan_data_list,
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
so_file_of_the_model <-  fn_build_unless_built(source_file   = file.path(models_folder, paste0(model, ".stan")),
                                               library_name  = paste0(model, "_Stan_math"),
                                               user_header   = NULL,
                                               AVX2_only     = FALSE)
so_file_of_the_fused_version <-  fn_build_unless_built(
                                       source_file   = file.path(models_folder,
                                                                 paste0(model, "_NicoStan_fused_AVX.stan")),
                                       library_name  = paste0(model, "_NicoStan_fused_", AVX_label),
                                       user_header   = NicoSIMD_user_header,
                                       AVX2_only     = !CPU_has_AVX512)
##
## ---- The points (the model's ordering of the parameters; the fused version's mapped by name):
##
fn_points <-  function(parameter_names) {

        centre <-  rep(0, length(parameter_names))
        centre[grepl("_SD", parameter_names)] <-  log(0.3)
        lapply(seq_len(n_points), function(point_index) {
              centre + stats::rnorm(length(parameter_names), mean = 0, sd = SD_of_the_noise_of_the_points)
        })

}
fn_evaluate_in_a_fresh_R_process <-  function(so_file,
                                              JSON_text,
                                              points) {

        callr::r(func  = fn_log_density_at_points,
                 args  = list(model_so_file           = so_file,
                              JSON_text_of_Stan_data  = JSON_text,
                              points                  = points,
                              propto                  = FALSE,
                              jacobian                = TRUE,
                              seed                    = 1))

}
fn_relabel_arrays_of_vectors <-  function(names_of_a_library,
                                          declarations_of_the_parameters) {

      block_names <-  unique(sub("\\..*$", "", names_of_a_library))
      for (block_name in block_names) {
            declaration <-  declarations_of_the_parameters[grepl(
                  paste0("(^|[^A-Za-z0-9_.])", block_name, " *;"),
                  declarations_of_the_parameters)]
            if (length(declaration) == 0) next
            is_an_array_of_row_vectors <-  grepl("array *\\[[^]]*\\] *row_vector", declaration[1])
            if (!is_an_array_of_row_vectors) next
            positions <-  grep(paste0("^", block_name, "\\."), names_of_a_library)
            if (length(positions) == 0) next
            first_index <-  unique(sub("\\..*$", "", sub(paste0("^", block_name, "\\."), "",
                                                          names_of_a_library[positions])))
            n_first <-  length(first_index)
            n_second <-  length(positions) / n_first
            if (n_first * n_second != length(positions)) next
            names_of_a_library[positions] <-  as.vector(t(outer(seq_len(n_first), seq_len(n_second),
                                                                function(i, j) paste0(block_name, ".", i, ".", j))))
      }
      names_of_a_library

}
JSON_text_of_the_build <-  fn_JSON_text_of(Stan_data_list)
names_of_the_parameters <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_model, JSON_text_of_the_build,
                                                             list())$unconstrained_parameter_names
names_of_the_fused_version <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_fused_version,
                                                                JSON_text_of_the_build,
                                                                list())$unconstrained_parameter_names
lines_of_the_fused_model <-  readLines(file.path(models_folder, paste0(model, "_NicoStan_fused_AVX.stan")),
                                       warn = FALSE)
declarations_of_the_fused_version <-  sub("//.*$", "",
                                          lines_of_the_fused_model[grep("^ *parameters *\\{",
                                                                        lines_of_the_fused_model)[1]:
                                                                   grep("^ *model *\\{",
                                                                        lines_of_the_fused_model)[1]])
names_of_the_fused_version <-  fn_relabel_arrays_of_vectors(names_of_the_fused_version,
                                                            declarations_of_the_fused_version)
stopifnot(setequal(names_of_the_parameters, names_of_the_fused_version))
mapping_into_the_fused_ordering <-  match(names_of_the_fused_version, names_of_the_parameters)
mapping_back <-  match(names_of_the_parameters, names_of_the_fused_version)
message(paste0("\033[36m", length(names_of_the_parameters), " unconstrained parameters; ", n_points,
               " points per setting; ", length(settings), " settings\033[0m"))
##
## ---- The log density and gradient of both at the same points, setting by setting:
##
rows <-  NULL
for (name_of_the_setting in names(settings)) {

      set.seed(seed + match(name_of_the_setting, names(settings)))
      JSON_text <-  fn_JSON_text_of(fn_Stan_data_of_a_setting(name_of_the_setting))
      points <-  fn_points(names_of_the_parameters)
      points_of_the_fused_version <-  lapply(points, function(point) point[mapping_into_the_fused_ordering])
      start <-  Sys.time()
      of_the_model <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_model, JSON_text, points)
      of_the_fused_version <-  fn_evaluate_in_a_fresh_R_process(so_file_of_the_fused_version, JSON_text,
                                                                points_of_the_fused_version)
      gradient_of_the_fused_version <-  of_the_fused_version$gradient[, mapping_back, drop = FALSE]
      every_point_finite <-  all(is.finite(of_the_model$log_density)) &&
                             all(is.finite(of_the_fused_version$log_density)) &&
                             all(is.finite(of_the_model$gradient)) &&
                             all(is.finite(of_the_fused_version$gradient))
      log_density_relative <-  max(abs(of_the_fused_version$log_density - of_the_model$log_density) /
                                   abs(of_the_model$log_density))
      gradient_relative <-  max(apply(abs(gradient_of_the_fused_version - of_the_model$gradient), 1, max) /
                                apply(abs(of_the_model$gradient), 1, max))
      row <-  data.frame(setting = name_of_the_setting,
                         n_points = n_points,
                         every_point_finite = every_point_finite,
                         largest_abs_log_density = max(abs(of_the_model$log_density)),
                         log_density_relative = log_density_relative,
                         gradient_relative_of_largest_element = gradient_relative,
                         passed_1e_8 = every_point_finite && log_density_relative <= tolerance_of_the_rule &&
                                       gradient_relative <= tolerance_of_the_rule,
                         passed_fused_gate_0_001 = every_point_finite &&
                                                   log_density_relative <= tolerance_of_the_fused_gate &&
                                                   gradient_relative <= tolerance_of_the_fused_gate,
                         seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
      rows <-  rbind(rows, row)
      message(paste0(if (row$passed_1e_8) "\033[32m" else if (row$passed_fused_gate_0_001) "\033[33m" else "\033[31m",
                     name_of_the_setting, ": log density ", formatC(log_density_relative, format = "e", digits = 2),
                     " relative, gradient ", formatC(gradient_relative, format = "e", digits = 2),
                     " of its largest element; ", if (row$passed_1e_8) "PASS (1e-8)" else
                     if (row$passed_fused_gate_0_001) "within 0.001 only" else "FAIL",
                     "; ", formatC(row$seconds, format = "f", digits = 1), " s\033[0m"))

}
print(rows, digits = 4)
file_of_the_check <-  file.path(checks_dir, paste0("check_lp_grad_ord_biv_FC_Se_Sp_fused_vs_Stan_math_", machine,
                                                   ".csv"))
utils::write.csv(rows, file_of_the_check, row.names = FALSE)
message(paste0(if (all(rows$passed_1e_8)) "\033[32mPASS" else if (all(rows$passed_fused_gate_0_001)) "\033[33mWITHIN 0.001 ONLY" else "\033[31mFAIL",
               ": ", sum(rows$passed_1e_8), " of ", nrow(rows), " settings within 1e-8, ",
               sum(rows$passed_fused_gate_0_001), " within 0.001; saved ", file_of_the_check, "\033[0m"))






















