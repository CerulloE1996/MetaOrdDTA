##
## ==== check_lp_grad_ord_biv_RC_Se_Sp_prior_vs_parent.R =============================================
##
## ---- The five ord_biv random-cutpoint Se/Sp files against their own parents, at 6 saved points of
## ---- the GAD-2 literature cell: the log density (every constant) and the gradient of every one of
## ---- them, in both states of use_Se_Sp_prior.
## ----   G1  each child's own gradient against finite differences of its own log density;
## ----   G2  use_Se_Sp_prior = 0 against the parent: 1e-8 (relative), gradient 1e-8 of its largest
## ----       element, the gate that the extension switched off leaves the parent;
## ----   G3  use_Se_Sp_prior = 1 with every cutoff chosen and the anchor at the location, reported
## ----       against the parent: it must differ, so that G2 cannot pass by the new term never
## ----       being reached;
## ----   G4  each fused file against the first Stan::math file, extension on: 0.001.
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
      NicoSIMD_user_header <-  file.path(root_of_the_work_folder, "R_packages", "NicoStan", "inst", "NicoSIMD",
                                         "inst", "include", "NicoSIMD",
                                         "stan_external_functions_of_NicoSIMD.hpp")
      check_name <-  "check_lp_grad_ord_biv_RC_Se_Sp_prior_vs_parent"
      test <-  "GAD_2"
      n_points <-  6
      SD_of_the_noise_of_the_points <-  0.3
      tolerance_of_the_rule <-  1e-8
      tolerance_of_the_fused_gate <-  0.001
      seed <-  1
      CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
      machine <-  if (CPU_has_AVX512) "local_HPC" else "laptop"
      compiler_flags <-  if (CPU_has_AVX512) {
            "-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq"
      } else {
            "-O3 -march=native -mtune=native -mfma -mavx2"
      }
      builds_folder <-  file.path(tempdir(), paste0("libraries_of_", check_name))
}
##
options(warning.length = 8170)
source(file.path(sim_study_dir, "sim_study_priors_helper_fns.R"))
suppressMessages(library(bridgestan))
dir.create(builds_folder, showWarnings = FALSE, recursive = TRUE)
##
## ---- The five pairs, each child beside the parent it was written from:
##
parent_names <-  c("DTA_MA_ord_biv_RANDthr_kappa",
                   "DTA_MA_ord_biv_RANDthr_kappa_2",
                   "DTA_MA_ord_biv_RANDthr_kappa_NicoStan_fused_AVX",
                   "DTA_MA_ord_biv_RANDthr_kappa_with_C_raw_in_the_nuisance_block_NicoStan_fused_AVX",
                   "DTA_MA_ord_biv_RANDthr_kappa_NicoStan_fused_AVX_reduce_sum_chunked")
child_names <-  paste0(parent_names, "_Se_Sp_prior")
fused_child <-  c(FALSE, FALSE, TRUE, TRUE, TRUE)
##
## ---- The data and priors of the literature cell of GAD-2, as the ord_biv random-cutpoint model:
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
                                                 random_thresholds = TRUE,
                                                 Dirichlet_random_effects_type = "kappa",
                                                 softplus = FALSE,
                                                 advanced_compile = FALSE,
                                                 force_recompile = FALSE,
                                                 quiet = TRUE,
                                                 compile = FALSE,
                                                 custom_file_name = NULL)
derived <-  readRDS(file.path(sim_study_dir, "prior_only_checks_of_the_induced_Dirichlet_alpha_sets",
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
stopifnot(n_thr >= 5)
##
## ---- The nine fields of the extension, at chosen cutoffs 2 and 5 of both classes, with values that
## ---- differ per class and per cutoff so that a transposed table is caught:
##
cutoffs_chosen <-  c(2, 5)
fn_segment_sums <-  function(alpha,
                             chosen_cutoffs) {

      bounds <-  c(0, chosen_cutoffs, length(alpha))
      sapply(seq_len(length(bounds) - 1), function(j) sum(alpha[(bounds[j] + 1):bounds[j + 1]]))

}
prior_Se_Sp_dirichlet_A <-  list(fn_segment_sums(derived$alpha[[1]], cutoffs_chosen),
                                 fn_segment_sums(derived$alpha[[2]], cutoffs_chosen))
prior_Se_Sp_beta <-  list(matrix(c(4, 6, 5, 3), nrow = 2, ncol = 2, byrow = TRUE),
                          matrix(c(5, 2, 3, 7), nrow = 2, ncol = 2, byrow = TRUE))
prior_Se_Sp_link_mean_SD <-  list(matrix(c(0.9, 1.4, -0.7, 0.6), nrow = 2, ncol = 2, byrow = TRUE),
                                  matrix(c(1.6, 0.8, 0.4, 1.9), nrow = 2, ncol = 2, byrow = TRUE))
fn_Se_Sp_fields <-  function(use_Se_Sp_prior,
                             n_chosen_of_each_class,
                             chosen_cutoffs_of_each_class) {

      list(use_Se_Sp_prior =  as.integer(use_Se_Sp_prior),
           Se_Sp_prior_type =  1L,
           n_chosen_cutoffs_max =  as.integer(length(cutoffs_chosen)),
           n_chosen_cutoffs =  as.integer(n_chosen_of_each_class),
           chosen_cutoffs =  chosen_cutoffs_of_each_class,
           prior_Se_Sp_dirichlet_A =  prior_Se_Sp_dirichlet_A,
           prior_Se_Sp_beta =  prior_Se_Sp_beta,
           prior_Se_Sp_link_mean_SD =  prior_Se_Sp_link_mean_SD)

}
fn_chosen_cutoffs <-  function(cutoffs_of_a_class) {

      matrix(as.integer(rep(cutoffs_of_a_class, 2)), nrow = 2, byrow = TRUE)

}
data_off <-  c(Stan_data_of_the_cell,
               fn_Se_Sp_fields(0, c(0L, 0L), fn_chosen_cutoffs(cutoffs_chosen)))
data_on <-  c(Stan_data_of_the_cell,
              fn_Se_Sp_fields(1, c(2L, 2L), fn_chosen_cutoffs(cutoffs_chosen)))
##
## ---- Every library built once in its own folder, and rebuilt when its Stan file changes:
##
fn_build <-  function(name_of_the_model,
                      fused) {

      source_file <-  file.path(models_folder, paste0(name_of_the_model, ".stan"))
      stopifnot(file.exists(source_file))
      in_the_build_folder <-  file.path(builds_folder, paste0(name_of_the_model, ".stan"))
      so_file <-  sub("\\.stan$", "_model.so", in_the_build_folder)
      if (file.exists(so_file) && file.exists(in_the_build_folder) &&
          identical(readLines(in_the_build_folder, warn = FALSE), readLines(source_file, warn = FALSE)) &&
          (file.info(so_file)$mtime > file.info(in_the_build_folder)$mtime)) {
            message(paste0("\033[36mAlready built: ", so_file, "\033[0m"))
            return(so_file)
      }
      stopifnot(file.copy(source_file, in_the_build_folder, overwrite = TRUE, copy.date = FALSE))
      Sys.setFileTime(in_the_build_folder, Sys.time())
      stanc_args <-  paste0("--include-paths=", functions_folder)
      make_args <-  c("STAN_THREADS=true", paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
      if (fused) {
            stanc_args <-  c(stanc_args, "--allow-undefined")
            make_args <-  c(make_args, paste0("USER_HEADER=", NicoSIMD_user_header))
      }
      start <-  Sys.time()
      library_of_the_model <-  bridgestan::compile_model(stan_file = in_the_build_folder,
                                                         stanc_args = stanc_args,
                                                         make_args = make_args)
      stopifnot(file.exists(library_of_the_model))
      message(paste0("\033[32mBuilt ", name_of_the_model, " in ",
                     formatC(as.numeric(difftime(Sys.time(), start, units = "secs")), format = "f", digits = 1),
                     " s\033[0m"))
      as.character(library_of_the_model)

}
fn_model <-  function(so_file,
                      Stan_data_list,
                      tag) {

      JSON_file <-  file.path(builds_folder, paste0(basename(so_file), "_", tag, ".json"))
      cmdstanr::write_stan_json(Stan_data_list, JSON_file)
      bridgestan::StanModel$new(so_file, JSON_file, 1L)

}
libraries_of_the_parents <-  vapply(seq_along(parent_names), function(i) {
      fn_build(parent_names[i], FALSE)
}, character(1))
libraries_of_the_children <-  vapply(seq_along(child_names), function(i) {
      fn_build(child_names[i], fused_child[i])
}, character(1))
parents <-  lapply(libraries_of_the_parents, fn_model, Stan_data_list = Stan_data_of_the_cell,
                   tag = "parent")
children <-  lapply(libraries_of_the_children, fn_model, Stan_data_list = data_off, tag = "off")
children_on <-  lapply(libraries_of_the_children, fn_model, Stan_data_list = data_on, tag = "on")
message(paste0("\033[36m", length(parent_names), " parents, ", length(child_names),
               " children each in both states, ", n_points, " points, ",
               parents[[1]]$param_unc_num(), " unconstrained parameters\033[0m"))
##
## ---- The points, in each library's own ordering of the unconstrained parameters:
##
fn_points <-  function(names_of_the_parameters) {

      centre <-  rep(0, length(names_of_the_parameters))
      centre[grepl("_SD", names_of_the_parameters)] <-  log(0.3)
      lapply(seq_len(n_points), function(point_index) {
            centre + stats::rnorm(length(names_of_the_parameters), mean = 0,
                                  sd = SD_of_the_noise_of_the_points)
      })

}
fn_permutation <-  function(from_these_names,
                            to_these_names) {

      stopifnot(setequal(from_these_names, to_these_names))
      match(to_these_names, from_these_names)

}
set.seed(seed + 1)
points <-  fn_points(parents[[1]]$param_unc_names())
##
## ---- The rows of every gate:
##
rows <-  NULL
fn_row <-  function(gate,
                    model,
                    state,
                    point,
                    what,
                    value,
                    tolerance,
                    must_differ) {

      passed <-  if (must_differ) value > tolerance else value <= tolerance
      data.frame(gate = gate, model = model, state = state, point = point, target = what,
                 value = value, tolerance = tolerance, passed = passed, stringsAsFactors = FALSE)

}
fn_density <-  function(model,
                        theta) {

      model$log_density_gradient(theta, propto = FALSE, jacobian = TRUE)$val

}
fn_gradient_by_differences <-  function(model,
                                        theta,
                                        step) {

      vapply(seq_along(theta), function(i) {
            values <-  vapply(c(-2, -1, 1, 2) * step, function(h) {
                  moved <-  theta
                  moved[i] <-  moved[i] + h
                  fn_density(model, moved)
            }, numeric(1))
            (values[1] - 8 * values[2] + 8 * values[3] - values[4]) / (12 * step)
      }, numeric(1))

}
fn_relative <-  function(reference,
                         other,
                         what) {

      if (what == "density") return(abs(other$val - reference$val) / abs(reference$val))
      max(abs(other$gradient - reference$gradient)) / max(abs(reference$gradient))

}
fn_two_rows <-  function(gate,
                         name_of_the_model,
                         state,
                         point,
                         reference,
                         other,
                         tolerance,
                         must_differ) {

      rbind(fn_row(gate, name_of_the_model, state, point, "log density relative",
                   fn_relative(reference, other, "density"), tolerance, must_differ),
            fn_row(gate, name_of_the_model, state, point, "gradient of its largest element",
                   fn_relative(reference, other, "gradient"), tolerance, must_differ))

}
for (i in seq_along(parent_names)) {

      child <-  children[[i]]
      parent <-  parents[[i]]
      permutation_into_the_child <-  fn_permutation(parent$param_unc_names(),
                                                    child$param_unc_names())
      for (point_index in seq_along(points)) {

            point <-  paste0("point_", point_index)
            theta_of_the_parent <-  points[[point_index]]
            theta_of_the_child <-  theta_of_the_parent[permutation_into_the_child]
            of_the_parent <-  parent$log_density_gradient(theta_of_the_parent, propto = FALSE,
                                                          jacobian = TRUE)
            ##
            ## ---- G1: the child's own gradient against finite differences of its own log density:
            ##
            checked <-  child$log_density_gradient(theta_of_the_child, propto = FALSE, jacobian = TRUE)
            by_differences <-  fn_gradient_by_differences(child, theta_of_the_child, 1e-5)
            rows <-  rbind(rows, fn_row("G1 own gradient against own finite differences", child_names[i],
                                        "off", point, "largest error over the largest element",
                                        max(abs(checked$gradient - by_differences)) /
                                              max(abs(checked$gradient)), tolerance_of_the_rule, FALSE))
            ##
            ## ---- G2: the option off against the parent:
            ##
            rows <-  rbind(rows, fn_two_rows("G2 option off against the parent", child_names[i], "off",
                                             point, of_the_parent, checked, tolerance_of_the_rule, FALSE))
            ##
            ## ---- G3: the extension on, against the parent, which it must differ from:
            ##
            with_the_extension <-  children_on[[i]]$log_density_gradient(theta_of_the_child,
                                                                         propto = FALSE, jacobian = TRUE)
            rows <-  rbind(rows, fn_two_rows("G3 extension on against the parent, must differ",
                                             child_names[i], "on", point, of_the_parent,
                                             with_the_extension, tolerance_of_the_rule, TRUE))
      }
}
##
## ---- G4: each fused file against the first Stan::math file, extension on, 0.001:
##
for (i in which(fused_child)) {

      permutation_into_the_first <-  fn_permutation(children_on[[1]]$param_unc_names(),
                                                    children_on[[i]]$param_unc_names())
      permutation_back <-  fn_permutation(children_on[[i]]$param_unc_names(),
                                          children_on[[1]]$param_unc_names())
      for (point_index in seq_along(points)) {

            point <-  paste0("point_", point_index)
            theta_of_the_first <-  points[[point_index]]
            of_the_first <-  children_on[[1]]$log_density_gradient(theta_of_the_first, propto = FALSE,
                                                                   jacobian = TRUE)
            of_this_one <-  children_on[[i]]$log_density_gradient(theta_of_the_first[permutation_into_the_first],
                                                                  propto = FALSE, jacobian = TRUE)
            of_this_one$gradient <-  of_this_one$gradient[permutation_back]
            rows <-  rbind(rows,
                           fn_row("G4 fused against Stan::math, extension on", child_names[i], "on", point,
                                  "log density relative", fn_relative(of_the_first, of_this_one, "density"),
                                  tolerance_of_the_fused_gate, FALSE),
                           fn_row("G4 fused against Stan::math, extension on", child_names[i], "on", point,
                                  "gradient of its largest element",
                                  fn_relative(of_the_first, of_this_one, "gradient"),
                                  tolerance_of_the_fused_gate, FALSE))
      }
}
##
## ---- The rows to disk and to the console:
##
file_of_the_rows <-  file.path(checks_dir, paste0(check_name, "_", machine, ".csv"))
utils::write.csv(rows, file_of_the_rows, row.names = FALSE)
for (gate in unique(rows$gate)) {

      of_this_gate <-  rows[rows$gate == gate, ]
      message(paste0(if (all(of_this_gate$passed)) "\033[32m" else "\033[31m", gate, ": ",
                     sum(of_this_gate$passed), " of ", nrow(of_this_gate), " rows pass, largest ",
                     formatC(max(of_this_gate$value), format = "e", digits = 2), "\033[0m"))
}
message(paste0("\033[36mrows: ", file_of_the_rows, "\033[0m"))
