##
## ==== check_lp_grad_ord_HSROC_FC_Se_Sp_fused_vs_Stan_math.R ====================================================
##
## ---- DTA_MA_ord_HSROC_FIXEDthr_NicoStan_fused_AVX_Se_Sp_prior.stan: (a) use_Se_Sp_prior = 0 against its parent
##      DTA_MA_ord_HSROC_FIXEDthr_NicoStan_fused_AVX.stan, likelihood on and prior-only; (b) every setting of
##      check_lp_grad_ord_HSROC_FC_Se_Sp.R against DTA_MA_ord_HSROC_FIXEDthr_Se_Sp_prior.stan (Stan math). Log
##      density (every constant) and gradient at 6 random unconstrained points; gate 1e-8 (relative; gradient of
##      its largest element), the fused AVX gate of 0.001 reported beside it.
##
{
      roots_of_the_work_folder <-  c("/home/enzocerullo/Documents/Work/PhD_work",
                                     "/home/enzo/Documents/Work/PhD_work")
      root_of_the_work_folder <-  roots_of_the_work_folder[dir.exists(roots_of_the_work_folder)][1]
      local_pkg_dir <-  file.path(root_of_the_work_folder, "R_packages", "MetaOrdDTA")
      sim_study_dir <-  file.path(local_pkg_dir, "sim_study")
      checks_dir <-  file.path(sim_study_dir, "prior_derivation_tools_checks")
      outputs_dir <-  file.path(dirname(local_pkg_dir), "MetaOrdDTA_outputs", "sim_study")
      models_folder <-  file.path(local_pkg_dir, "inst", "stan_models", "stan_models_MA")
      functions_folder <-  file.path(local_pkg_dir, "inst", "stan_functions")
      parent <-  "DTA_MA_ord_HSROC_FIXEDthr"
      model <-  "DTA_MA_ord_HSROC_FIXEDthr_Se_Sp_prior"
      test <-  "GAD_2"
      n_points <-  6
      SD_of_the_noise_of_the_points <-  0.3
      finite_difference_step <-  1e-3
      tolerance_of_the_rule <-  1e-8
      tolerance_of_the_fused_gate <-  0.001
      NicoSIMD_user_header <-  file.path(root_of_the_work_folder, "R_packages", "NicoStan", "inst", "NicoSIMD",
                                         "inst", "include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp")
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
##
## ---- The data of GAD-2 and the package's default priors of the ord_HSROC fixed-cutpoint model:
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
                                                 model_parameterisation = "ord_HSROC",
                                                 random_thresholds = FALSE,
                                                 Dirichlet_random_effects_type = "fixed",
                                                 softplus = FALSE,
                                                 advanced_compile = FALSE,
                                                 force_recompile = FALSE,
                                                 quiet = TRUE,
                                                 compile = FALSE,
                                                 custom_file_name = NULL)
outs_data <-  model_prep_obj$internal_obj$outs_data
Stan_data_of_the_cell <-  c(outs_data$stan_data_list, model_prep_obj$priors)
Stan_data_of_the_cell <-  MetaOrdDTA:::R_fn_remove_duplicates_from_list(Stan_data_of_the_cell)
Stan_data_of_the_cell$softplus <-  0
Stan_data_of_the_cell$use_probit_link <-  1
n_thr <-  Stan_data_of_the_cell$n_thr
n_studies <-  Stan_data_of_the_cell$n_studies
##
## ---- The settings of the Se/Sp prior (padded rows of the data are never read by the model):
##
fn_Se_Sp_fields <-  function(use_Se_Sp_prior,
                             Se_Sp_prior_type,
                             chosen_of_each_class,
                             A_of_each_class,
                             beta_ab_of_each_class,
                             mean_SD_of_each_class) {

      n_chosen <-  sapply(chosen_of_each_class, length)
      n_chosen_max <-  max(n_chosen)
      fn_pad_rows <-  function(values,
                               n_rows,
                               padding) {

            if (n_rows == 0) return(matrix(padding, n_chosen_max, length(padding), byrow = TRUE))
            rbind(values, matrix(padding, n_chosen_max - n_rows, length(padding), byrow = TRUE))

      }
      list(use_Se_Sp_prior = use_Se_Sp_prior,
           Se_Sp_prior_type = Se_Sp_prior_type,
           n_chosen_max = n_chosen_max,
           n_chosen = n_chosen,
           chosen_cutoffs = t(sapply(1:2, function(c) c(chosen_of_each_class[[c]],
                                                         rep(1, n_chosen_max - n_chosen[c])))),
           prior_Se_Sp_dirichlet_A = lapply(1:2, function(c) c(A_of_each_class[[c]],
                                                                rep(1, n_chosen_max + 1 -
                                                                       length(A_of_each_class[[c]])))),
           prior_Se_Sp_beta = lapply(1:2, function(c) fn_pad_rows(beta_ab_of_each_class[[c]], n_chosen[c],
                                                                   c(1, 1))),
           prior_Se_Sp_link_mean_SD = lapply(1:2, function(c) fn_pad_rows(mean_SD_of_each_class[[c]],
                                                                           n_chosen[c], c(0, 1))))

}
beta_ab_of_three <-  cbind(c(4, 6, 3), c(6, 3, 2))
mean_SD_of_three <-  cbind(c(-0.5, 0.2, 1.0), c(1.0, 0.8, 1.2))
settings <-  list()
settings[["off_likelihood_on"]] <-  fn_Se_Sp_fields(0, 4, list(c(2, 5), integer(0)), list(c(1, 1, 1), 1),
                                                    list(cbind(c(1, 1), c(1, 1)), NULL),
                                                    list(cbind(c(0, 0), c(1, 1)), NULL))
settings[["off_prior_only"]] <-  settings[["off_likelihood_on"]]
{

      settings[[paste0("Sp_all_type_4_N_0_1")]] <-
            fn_Se_Sp_fields(1, 4, list(1:n_thr, integer(0)), list(rep(1, n_thr + 1), 1),
                            list(matrix(1, n_thr, 2), NULL), list(cbind(rep(0, n_thr), rep(1, n_thr)), NULL))
      settings[[paste0("Sp_2_5_type_1")]] <-
            fn_Se_Sp_fields(1, 1, list(c(2, 5), integer(0)), list(c(2, 3, 4), 1),
                            list(beta_ab_of_three[1:2, ], NULL), list(mean_SD_of_three[1:2, ], NULL))

}
settings[["Sp_2_5_type_2"]] <-  fn_Se_Sp_fields(1, 2, list(c(2, 5), integer(0)), list(c(2, 3, 4), 1),
                                                         list(beta_ab_of_three[1:2, ], NULL),
                                                         list(mean_SD_of_three[1:2, ], NULL))
settings[["Sp_2_5_type_3"]] <-  fn_Se_Sp_fields(1, 3, list(c(2, 5), integer(0)), list(c(2, 3, 4), 1),
                                                         list(beta_ab_of_three[1:2, ], NULL),
                                                         list(mean_SD_of_three[1:2, ], NULL))
settings[["Se_1_3_6_type_4"]] <-  fn_Se_Sp_fields(1, 4, list(integer(0), c(1, 3, 6)),
                                                           list(1, c(2, 3, 4, 2)),
                                                           list(NULL, beta_ab_of_three),
                                                           list(NULL, mean_SD_of_three))
settings[["Sp_2_4_5_Se_4_type_4"]] <-  fn_Se_Sp_fields(1, 4, list(c(2, 4, 5), 4),
                                                                list(c(2, 3, 4, 2), c(3, 2)),
                                                                list(beta_ab_of_three,
                                                                     beta_ab_of_three[2, , drop = FALSE]),
                                                                list(mean_SD_of_three,
                                                                     mean_SD_of_three[2, , drop = FALSE]))
settings[["Sp_2_4_5_Se_2_5_type_4"]] <-  fn_Se_Sp_fields(1, 4, list(c(2, 4, 5), c(2, 5)),
                                                                  list(c(2, 3, 4, 2), c(3, 2, 2)),
                                                                  list(beta_ab_of_three, beta_ab_of_three[1:2, ]),
                                                                  list(mean_SD_of_three,
                                                                       mean_SD_of_three[c(3, 1), ]))
settings[["Sp_2_4_5_Se_2_5_type_2"]] <-  settings[["Sp_2_4_5_Se_2_5_type_4"]]
settings[["Sp_2_4_5_Se_2_5_type_2"]]$Se_Sp_prior_type <-  2
settings[["Se_1_3_6_Sp_1_6_type_1"]] <-  fn_Se_Sp_fields(1, 1, list(c(1, 6), c(1, 3, 6)),
                                                                  list(c(2, 3, 2), c(2, 3, 4, 2)),
                                                                  list(beta_ab_of_three[1:2, ], beta_ab_of_three),
                                                                  list(mean_SD_of_three[1:2, ], mean_SD_of_three))
settings[["Sp_2_4_5_Se_2_5_type_3_softplus"]] <-  settings[["Sp_2_4_5_Se_2_5_type_4"]]
settings[["Sp_2_4_5_Se_2_5_type_3_softplus"]]$Se_Sp_prior_type <-  3
settings[["Sp_2_4_5_Se_2_5_type_4_baseline_case_1_3"]] <-  settings[["Sp_2_4_5_Se_2_5_type_4"]]
settings[["Sp_2_4_5_Se_2_5_type_1_logit"]] <-  settings[["Sp_2_4_5_Se_2_5_type_4"]]
settings[["Sp_2_4_5_Se_2_5_type_1_logit"]]$Se_Sp_prior_type <-  1
##
fn_Stan_data_of_a_setting <-  function(name_of_the_setting) {

      Stan_data_list <-  c(Stan_data_of_the_cell, settings[[name_of_the_setting]])
      if (name_of_the_setting != "off_likelihood_on") {
            Stan_data_list$K_fold_CV_indicator <-  rep(0, n_studies)
      }
      if (grepl("softplus", name_of_the_setting)) Stan_data_list$softplus <-  1
      if (grepl("baseline_case_1_3", name_of_the_setting)) Stan_data_list$baseline_case <-  1.3
      if (grepl("logit", name_of_the_setting)) Stan_data_list$use_probit_link <-  0
      Stan_data_list

}
##
## ---- Only the fields the new model declares; a field of one value declared as a container is an array of 1:
##
lines_of_the_model <-  readLines(file.path(models_folder, paste0(model, ".stan")), warn = FALSE)
lines_of_the_data_block <-  sub("//.*$", "",
                                lines_of_the_model[grep("^data *\\{", lines_of_the_model)[1]:
                                                   grep("^transformed data", lines_of_the_model)[1]])
names_declared <-  sub("^.*[] >]([A-Za-z0-9_]+) *;.*$", "\\1",
                       grep(";", lines_of_the_data_block, value = TRUE))
fn_declared_fields_only <-  function(Stan_data_list) {

      Stan_data_list <-  Stan_data_list[intersect(names(Stan_data_list), names_declared)]
      for (name in names(Stan_data_list)) {
            declared_as_a_container <-  any(grepl(paste0("^ *(vector|row_vector|array)[^;]*[] >]", name, " *;"),
                                                  lines_of_the_data_block))
            if (declared_as_a_container && (length(Stan_data_list[[name]]) == 1) &&
                is.null(dim(Stan_data_list[[name]])) && !is.list(Stan_data_list[[name]])) {
                  Stan_data_list[[name]] <-  array(Stan_data_list[[name]], dim = 1)
            }
      }
      Stan_data_list

}
fn_JSON_text_of <-  function(Stan_data_list) {

      Stan_data_list <-  fn_declared_fields_only(Stan_data_list)
      JSON_file <-  tempfile(fileext = ".json")
      cmdstanr::write_stan_json(data = Stan_data_list, file = JSON_file)
      paste(readLines(JSON_file, warn = FALSE), collapse = "\n")

}
##
## ---- The three libraries: the Stan math copy (as built by check_lp_grad_ord_HSROC_FC_Se_Sp.R) and the two fused
##      AVX files (built unless built already, both at once):
##
fn_build_unless_built <-  function(source_file,
                                   library_name) {

      Stan_model_file_path <-  file.path(libraries_folder, paste0(library_name, ".stan"))
      so_file <-  sub("\\.stan$", "_model.so", Stan_model_file_path)
      if (file.exists(so_file) && file.exists(Stan_model_file_path) &&
          identical(readLines(Stan_model_file_path, warn = FALSE), readLines(source_file, warn = FALSE)) &&
          (file.info(so_file)$mtime > file.info(Stan_model_file_path)$mtime)) {
            message(paste0("\033[36mAlready built: ", so_file, "\033[0m"))
            return(NULL)
      }
      stopifnot(file.copy(source_file, Stan_model_file_path, overwrite = TRUE, copy.date = FALSE))
      Sys.setFileTime(Stan_model_file_path, Sys.time())
      make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                      paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
      if (!CPU_has_AVX512) make_args <-  c(make_args, "CPPFLAGS_OPTIM=-DBAYESMVP_FORCE_AVX2")
      callr::r_bg(func = function(Stan_model_file_path,
                                  Stan_data_list,
                                  user_header,
                                  stanc_args,
                                  make_args) {

                        NicoStan:::initialise_model( Model_type            = "Stan",
                                                     sample_nuisance       = FALSE,
                                                     Stan_data_list        = Stan_data_list,
                                                     Stan_model_file_path  = Stan_model_file_path,
                                                     Stan_cpp_user_header  = user_header,
                                                     stanc_args            = stanc_args,
                                                     make_args             = make_args,
                                                     compile               = TRUE)

                  },
                  args = list(Stan_model_file_path = Stan_model_file_path,
                              Stan_data_list = Stan_data_list_of_the_build,
                              user_header = NicoSIMD_user_header,
                              stanc_args = c(paste0("--include-paths=", functions_folder), "--allow-undefined"),
                              make_args = make_args),
                  stdout = paste0(Stan_model_file_path, "_build_log.txt"),
                  stderr = "2>&1")

}
Stan_data_list_of_the_build <-  fn_declared_fields_only(fn_Stan_data_of_a_setting("off_prior_only"))
library_of <-  c(fused_parent = paste0(parent, "_NicoStan_fused_", AVX_label),
                 fused_model = paste0(model, "_NicoStan_fused_", AVX_label))
so_file_of <-  c(Stan_math_model = file.path(libraries_folder, paste0(model, "_Stan_math_model.so")),
                 file.path(libraries_folder, paste0(library_of, "_model.so")))
names(so_file_of) <-  c("Stan_math_model", names(library_of))
stopifnot(file.exists(so_file_of[["Stan_math_model"]]),
          identical(readLines(file.path(libraries_folder, paste0(model, "_Stan_math.stan")), warn = FALSE),
                    readLines(file.path(models_folder, paste0(model, ".stan")), warn = FALSE)))
start <-  Sys.time()
builds <-  list(fused_parent = fn_build_unless_built(file.path(models_folder,
                                                               paste0(parent, "_NicoStan_fused_AVX.stan")),
                                                     library_of[["fused_parent"]]),
                fused_model = fn_build_unless_built(file.path(models_folder,
                                                              paste0(model, "_NicoStan_fused_AVX.stan")),
                                                    library_of[["fused_model"]]))
while (any(sapply(builds, function(build) !is.null(build) && build$is_alive()))) {
      Sys.sleep(5)
      message(paste0("\033[36mBuilding (", paste(names(builds)[sapply(builds, function(build) !is.null(build) &&
                                                                        build$is_alive())], collapse = ", "),
                     "): ", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), " s\033[0m"))
}
for (name in names(builds)) {
      if (!is.null(builds[[name]])) {
            cat(readLines(file.path(libraries_folder, paste0(library_of[[name]], ".stan_build_log.txt"))),
                sep = "\n")
      }
      stopifnot(file.exists(so_file_of[[name]]))
}
message(paste0("\033[32mLibraries ready in ", round(as.numeric(difftime(Sys.time(), start, units = "secs"))),
               " s\033[0m"))
##
## ---- The points (one ordering: the three parameters blocks are identical):
##
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
JSON_text_of_the_build <-  fn_JSON_text_of(fn_Stan_data_of_a_setting("off_prior_only"))
names_of_u <-  fn_evaluate_in_a_fresh_R_process(so_file_of[["Stan_math_model"]], JSON_text_of_the_build,
                                                list())$unconstrained_parameter_names
for (name in c("fused_parent", "fused_model")) {
      stopifnot(identical(fn_evaluate_in_a_fresh_R_process(so_file_of[[name]], JSON_text_of_the_build,
                                                           list())$unconstrained_parameter_names, names_of_u))
}
fn_points <-  function() {

      centre <-  rep(0, length(names_of_u))
      centre[grepl("_SD", names_of_u)] <-  log(0.3)
      C_raw_positions <-  grep("^C_raw_vec", names_of_u)
      centre[C_raw_positions] <-  stats::qlogis((c(-1.2, rep(-0.8, n_thr - 1)) + 15) / 17.5)
      lapply(seq_len(n_points), function(point_index) {
            centre + stats::rnorm(length(names_of_u), mean = 0, sd = SD_of_the_noise_of_the_points)
      })

}
message(paste0("\033[36m", length(names_of_u), " unconstrained parameters; ", n_points, " points per setting; ",
               length(settings), " settings\033[0m"))
##
## ---- Gate 4(a) (use_Se_Sp_prior = 0: fused copy against the fused parent) and gate 4(b) (fused copy against the
##      Stan math copy), per setting:
##
fn_compare <-  function(evaluation,
                        reference) {

      every_point_finite <-  all(is.finite(evaluation$log_density)) && all(is.finite(evaluation$gradient)) &&
                             all(is.finite(reference$log_density)) && all(is.finite(reference$gradient))
      c(every_point_finite = every_point_finite,
        log_density_relative = max(abs(evaluation$log_density - reference$log_density) /
                                   abs(reference$log_density)),
        gradient_relative = max(apply(abs(evaluation$gradient - reference$gradient), 1, max) /
                                apply(abs(reference$gradient), 1, max)))

}
rows <-  NULL
for (name_of_the_setting in names(settings)) {

      index_of_the_setting <-  match(name_of_the_setting, names(settings))
      set.seed(seed + index_of_the_setting)
      Stan_data <-  fn_Stan_data_of_a_setting(name_of_the_setting)
      JSON_text <-  fn_JSON_text_of(Stan_data)
      points <-  fn_points()
      start <-  Sys.time()
      of_the_fused_model <-  fn_evaluate_in_a_fresh_R_process(so_file_of[["fused_model"]], JSON_text, points)
      comparisons <-  list(fused_against_Stan_math = fn_compare(of_the_fused_model,
                                                                fn_evaluate_in_a_fresh_R_process(
                                                                      so_file_of[["Stan_math_model"]],
                                                                      JSON_text, points)))
      if (Stan_data$use_Se_Sp_prior == 0) {
            comparisons[["fused_against_the_fused_parent"]] <-  fn_compare(of_the_fused_model,
                                                                           fn_evaluate_in_a_fresh_R_process(
                                                                                 so_file_of[["fused_parent"]],
                                                                                 JSON_text, points))
      }
      for (name_of_the_comparison in names(comparisons)) {
            comparison <-  comparisons[[name_of_the_comparison]]
            row <-  data.frame(setting = name_of_the_setting,
                               comparison = name_of_the_comparison,
                               n_points = n_points,
                               every_point_finite = as.logical(comparison[["every_point_finite"]]),
                               log_density_relative = comparison[["log_density_relative"]],
                               gradient_relative_of_largest_element = comparison[["gradient_relative"]],
                               passed_1e_8 = as.logical(comparison[["every_point_finite"]]) &&
                                             comparison[["log_density_relative"]] <= tolerance_of_the_rule &&
                                             comparison[["gradient_relative"]] <= tolerance_of_the_rule,
                               passed_fused_gate_0_001 = as.logical(comparison[["every_point_finite"]]) &&
                                                         comparison[["log_density_relative"]] <=
                                                               tolerance_of_the_fused_gate &&
                                                         comparison[["gradient_relative"]] <=
                                                               tolerance_of_the_fused_gate,
                               seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
            rows <-  rbind(rows, row)
            message(paste0(if (row$passed_1e_8) "\033[32m" else if (row$passed_fused_gate_0_001) "\033[33m" else
                           "\033[31m", "[", index_of_the_setting, "/", length(settings), "] ",
                           name_of_the_setting,
                           " (", name_of_the_comparison, "): log density ",
                           formatC(row$log_density_relative, format = "e", digits = 2), " relative, gradient ",
                           formatC(row$gradient_relative_of_largest_element, format = "e", digits = 2),
                           " of its largest element; ", if (row$passed_1e_8) "PASS (1e-8)" else
                           if (row$passed_fused_gate_0_001) "within 0.001 only" else "FAIL", "; ",
                           formatC(row$seconds, format = "f", digits = 1), " s\033[0m"))
      }

}
print(rows, digits = 4)
file_of_the_check <-  file.path(checks_dir, paste0("check_lp_grad_ord_HSROC_FC_Se_Sp_fused_vs_Stan_math_",
                                                   machine, ".csv"))
utils::write.csv(rows, file_of_the_check, row.names = FALSE)
message(paste0(if (all(rows$passed_1e_8)) "\033[32mPASS" else if (all(rows$passed_fused_gate_0_001))
               "\033[33mWITHIN 0.001 ONLY" else "\033[31mFAIL", ": ", sum(rows$passed_1e_8), " of ", nrow(rows),
               " comparisons within 1e-8, ", sum(rows$passed_fused_gate_0_001), " within 0.001; saved ",
               file_of_the_check, "\033[0m"))





















