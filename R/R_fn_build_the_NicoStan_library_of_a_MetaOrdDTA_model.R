##
## ==== R_fn_build_the_NicoStan_library_of_a_MetaOrdDTA_model.R ================================================
##
## ---- Builds the NicoStan library of one MetaOrdDTA model once, without fitting it, so that no worker of a run
##      builds it; the preparation is the same as R_fn_sample_model_with_NicoStan.R, step for step.
##
#' Build the NicoStan library of one MetaOrdDTA model, once, without fitting
#'
#' @param stan_data_list The model's data list, exactly as it is given to the fit (the library's sizes come from
#'   it).
#' @param stan_model_file_name The model's file name, with or without ".stan" and with or without
#'   "_NicoStan_fused_AVX", e.g. "DTA_MA_ord_biv_FIXEDthr".
#' @param NicoStan_args The same NicoStan arguments the fits are given, e.g. list(sample_nuisance = TRUE);
#'   NicoStan_args$nuisance_params, if given, wins over the package's own choice.
#' @return The path of the library that was built, as a character string.
#' @export
fn_build_the_NicoStan_library_of_a_MetaOrdDTA_model <-  function( stan_data_list,
                                                                  stan_model_file_name,
                                                                  NicoStan_args ) {

      stopifnot(is.list(stan_data_list), length(stan_data_list) > 0)
      stopifnot(length(stan_model_file_name) == 1, nzchar(stan_model_file_name))
      stopifnot(is.list(NicoStan_args))
      ##
      ## ---- (1) The model and its NicoStan version with the fused AVX functions
      name_of_the_model <-  sub("(_NicoStan_fused_AVX)?\\.stan$", "", basename(stan_model_file_name))
      folder_of_the_models <-  if (grepl("^DTA_NMA_", name_of_the_model)) "stan_models_NMA" else "stan_models_MA"
      Stan_file_of_the_fused_version <-  system.file("stan_models", folder_of_the_models,
                                                     paste0(name_of_the_model, "_NicoStan_fused_AVX.stan"),
                                                     package = "MetaOrdDTA")
      if (!nzchar(Stan_file_of_the_fused_version)) {
            stop(paste0("fn_build_the_NicoStan_library_of_a_MetaOrdDTA_model: no NicoStan version (",
                        name_of_the_model, "_NicoStan_fused_AVX.stan) of the model ", name_of_the_model, " in ",
                        system.file("stan_models", folder_of_the_models, package = "MetaOrdDTA")))
      }
      if (!is.null(NicoStan_args$Stan_file_of_the_fused_version)) {
            Stan_file_of_the_fused_version <-  NicoStan_args$Stan_file_of_the_fused_version
            name_of_the_model <-  sub("_NicoStan_fused_AVX\\.stan$", "",
                                      basename(Stan_file_of_the_fused_version))
      }
      NicoStan_args$Stan_file_of_the_fused_version <-  NULL
      folder_of_the_Stan_functions <-  system.file("stan_functions", package = "MetaOrdDTA")
      user_header_of_NicoSIMD <-  system.file("include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp",
                                              package = "NicoSIMD")
      if (!nzchar(user_header_of_NicoSIMD)) {
            stop("fn_build_the_NicoStan_library_of_a_MetaOrdDTA_model: the NicoSIMD package is not installed")
      }
      ## (the same 0 / 1 switches and length-1 containers as the fit, so that the library is the same)
      Stan_data_list <-  stan_data_list
      for (name in names(Stan_data_list)) {
            if (is.logical(Stan_data_list[[name]])) Stan_data_list[[name]] <-  as.integer(Stan_data_list[[name]])
      }
      lines_of_the_model <-  readLines(Stan_file_of_the_fused_version, warn = FALSE)
      start_of_the_data_block <-  grep("^data *\\{", lines_of_the_model)[1]
      end_of_the_data_block <-  grep("^transformed data", lines_of_the_model)[1]
      lines_of_the_data_block <-  sub("//.*$", "",
                                      lines_of_the_model[start_of_the_data_block:end_of_the_data_block])
      for (name in names(Stan_data_list)) {
            declared_as_a_container <-  any(grepl(paste0("^ *(vector|row_vector|array)[^;]*[] >]", name, " *;"),
                                                  lines_of_the_data_block))
            if (declared_as_a_container && (length(Stan_data_list[[name]]) == 1) &&
                is.null(dim(Stan_data_list[[name]]))) {
                  Stan_data_list[[name]] <-  array(Stan_data_list[[name]], dim = 1)
            }
      }
      ##
      ## ---- (2) The library, in the same cache folder and under the same name as the fit
      CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
      SIMD_label <-  if (CPU_has_AVX512) "AVX512" else "AVX2"
      compiler_flags <-  if (CPU_has_AVX512) {
            paste("-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq",
                  "-DNICOSIMD_SIMD_LUB_CONSTRAIN")
      } else {
            "-O3 -march=native -mtune=native -mfma -mavx2 -DNICOSIMD_SIMD_LUB_CONSTRAIN"
      }
      library_of_MetaOrdDTA <-  basename(dirname(system.file(package = "MetaOrdDTA")))
      cache_folder <-  file.path(tools::R_user_dir("MetaOrdDTA", which = "cache"), "NicoStan_fused_AVX",
                                 library_of_MetaOrdDTA)
      lock_dir <-  file.path(cache_folder, paste0(name_of_the_model, "_", SIMD_label, "_build_lock"))
      ## (one lock per model and backend, as in the fit; a lock older than 15 minutes is left by a killed process)
      seconds_waited_for_the_lock <-  0
      while (!dir.create(lock_dir, showWarnings = FALSE)) {
            if (seconds_waited_for_the_lock %% 30 == 0) {
                  message(paste0("\033[36m", "waiting for the build lock of ", name_of_the_model, " (",
                                 SIMD_label, "): ", lock_dir, "\033[0m"))
            }
            Sys.sleep(1)
            seconds_waited_for_the_lock <-  seconds_waited_for_the_lock + 1
            age_of_the_lock_in_minutes <-  tryCatch(as.numeric(difftime(Sys.time(), file.info(lock_dir)$mtime,
                                                                       units = "mins")),
                                                    error = function(e) NA_real_)
            if (is.finite(age_of_the_lock_in_minutes) && (age_of_the_lock_in_minutes > 15)) {
                  message(paste0("\033[31m", "removing a stale build lock (",
                                 round(age_of_the_lock_in_minutes, 1), " minutes old) of ", name_of_the_model,
                                 "\033[0m"))
                  unlink(lock_dir, recursive = TRUE, force = TRUE)
            }
      }
      ## (released when the build is finished, or by on.exit on an error)
      this_fit_holds_the_lock <-  TRUE
      on.exit(if (this_fit_holds_the_lock) unlink(lock_dir, recursive = TRUE), add = TRUE)
      folder_of_the_copied_Stan_functions <-  file.path(cache_folder, "stan_functions")
      dir.create(folder_of_the_copied_Stan_functions, recursive = TRUE, showWarnings = FALSE)
      Stan_model_file_path <-  file.path(cache_folder, paste0(name_of_the_model, "_NicoStan_fused_",
                                                              SIMD_label, ".stan"))
      so_file <-  sub("\\.stan$", "_model.so", Stan_model_file_path)
      files_of_the_Stan_functions <-  list.files(folder_of_the_Stan_functions, pattern = "\\.stan$")
      something_has_changed <-  !file.exists(Stan_model_file_path) ||
            !identical(readLines(Stan_model_file_path, warn = FALSE), lines_of_the_model)
      for (file_name in files_of_the_Stan_functions) {
            copied_file <-  file.path(folder_of_the_copied_Stan_functions, file_name)
            lines_of_the_file <-  readLines(file.path(folder_of_the_Stan_functions, file_name), warn = FALSE)
            if (!file.exists(copied_file) ||
                !identical(readLines(copied_file, warn = FALSE), lines_of_the_file)) {
                  writeLines(lines_of_the_file, copied_file)
                  something_has_changed <-  TRUE
            }
      }
      if (something_has_changed) {
            writeLines(lines_of_the_model, Stan_model_file_path)
            if (file.exists(so_file)) {
                  file.rename(so_file, paste0(so_file, ".built_before_",
                                              format(Sys.time(), "%Y_%m_%d_%H%M%S")))
            }
      }
      make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                      paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
      if (!CPU_has_AVX512) make_args <-  c(make_args, "CPPFLAGS_OPTIM=-DBAYESMVP_FORCE_AVX2")
      stanc_args <-  c(paste0("--include-paths=", folder_of_the_copied_Stan_functions), "--allow-undefined")
      ##
      ## ---- (3) The nuisance blocks: the fit's own if given, the package's choice otherwise
      nuisance_params <-  NicoStan_args$nuisance_params
      NicoStan_args$nuisance_params <-  NULL
      if (is.null(nuisance_params)) {
            parameterisation_of_the_model <-  if (grepl("^DTA_NMA_", name_of_the_model)) "NMA" else
                  if (grepl("Jones", name_of_the_model)) "Jones" else
                  if (grepl("ord_HSROC", name_of_the_model)) "ordinal_HSROC" else
                  if (grepl("ord_biv", name_of_the_model)) "ord_biv" else "Bivariate"
            nuisance_params <-  fn_nuisance_params_of_a_MetaOrdDTA_model(
                  name_of_the_model        = name_of_the_model,
                  model_parameterisation   = parameterisation_of_the_model,
                  random_thresholds        = grepl("RANDthr", name_of_the_model))
            for (block in nuisance_params) {
                  message(paste0("\033[36m", "nuisance block ", block, ": ",
                                 fn_reason_of_a_nuisance_block_of_a_MetaOrdDTA_model(
                                       name_of_a_block   = block,
                                       name_of_the_model = name_of_the_model), "\033[0m"))
            }
      }
      sample_nuisance <-  isTRUE(NicoStan_args$sample_nuisance) || !is.null(nuisance_params)
      n_nuisance_override <-  NicoStan_args$n_nuisance_override
      ##
      ## ---- The named nuisance blocks: NicoStan's nico_reduce_sum version only for a model that calls reduce_sum --------------------
      ##
      ## (without n_nuisance_override the blocks go to NicoStan by name: a model that calls reduce_sum or map_rect
      ##  gets its nico_reduce_sum version; a model that calls neither keeps them in its own parameters block,
      ##  declared first (<model>_nuisance_first.stan), and NicoStan samples them as its nuisance block)
      nuisance_params_named_for_NicoStan <-  if (is.null(n_nuisance_override)) {
            nuisance_params
      } else {
            NULL
      }
      ##
      ## ---- (4) The build, with no fit
      ## (a library already built: the lock is released before NicoStan loads it; only a build holds it)
      library_to_load <-  if (!is.null(nuisance_params_named_for_NicoStan)) {
            fn_library_NicoStan_builds_for_a_MetaOrdDTA_model(so_file     = so_file,
                                                              stanc_args  = stanc_args)
      } else {
            so_file
      }
      if (!something_has_changed && file.exists(library_to_load)) {
            unlink(lock_dir, recursive = TRUE)
            this_fit_holds_the_lock <-  FALSE
      }
      start_of_the_build <-  Sys.time()
      model_NicoStan <-  NicoStan::Nico_model$new( Model_type            = "Stan",
                                                   Stan_data_list        = Stan_data_list,
                                                   Stan_model_file_path  = Stan_model_file_path,
                                                   Stan_cpp_user_header  = user_header_of_NicoSIMD,
                                                   stanc_args            = stanc_args,
                                                   make_args             = make_args,
                                                   nuisance_params       =
                                                         nuisance_params_named_for_NicoStan,
                                                   n_nuisance_override   = n_nuisance_override,
                                                   sample_nuisance       = sample_nuisance)
      library_built <-  if (!is.null(nuisance_params_named_for_NicoStan)) {
            fn_library_NicoStan_builds_for_a_MetaOrdDTA_model(so_file     = so_file,
                                                              stanc_args  = stanc_args)
      } else {
            so_file
      }
      message(paste0("\033[36m", "NicoStan library of ", name_of_the_model, " (fused AVX, ", SIMD_label, "): ",
                     library_built, " (",
                     formatC(as.numeric(difftime(Sys.time(), start_of_the_build, units = "secs")),
                             format = "f", digits = 1), " s)", "\033[0m"))
      ## (the build either made the library or it stops here)
      if (!file.exists(library_built)) {
            stop(paste0("fn_build_the_NicoStan_library_of_a_MetaOrdDTA_model: the build of ", name_of_the_model,
                        " finished but the library ", library_built, " is not there"))
      }
      if (this_fit_holds_the_lock) {
            unlink(lock_dir, recursive = TRUE)
            this_fit_holds_the_lock <-  FALSE
      }
      return(library_built)
}


#' The library of the nico_reduce_sum version of a MetaOrdDTA model
#'
#' The library built from <model>_nico_reduce_sum.stan is <model>_nico_reduce_sum_model.so, or
#' <model>_reduce_sum_chunked_model.so for a version written as <model>_reduce_sum_chunked.stan.
#'
#' @param so_file The library of the model itself, <model>_model.so.
#' @return The path of the library of its nico_reduce_sum version.
#' @noRd
fn_library_of_the_nico_reduce_sum_version_of_a_MetaOrdDTA_model <-  function(
                    so_file
) {

      file_under_the_former_name <-  sub("_model\\.so$", "_reduce_sum_chunked.stan", so_file)
      if (file.exists(file_under_the_former_name)) {
            return(sub("_model\\.so$", "_reduce_sum_chunked_model.so", so_file))
      }
      return(sub("_model\\.so$", "_nico_reduce_sum_model.so", so_file))

}



##
## ---- The library NicoStan builds for a MetaOrdDTA model whose nuisance blocks are named --------------------------------------------
##
#' The library NicoStan builds for a MetaOrdDTA model whose nuisance blocks are named
#'
#' A model that calls reduce_sum or map_rect: the library of its nico_reduce_sum version. A model that calls
#' neither: the library of <model>_nuisance_first.stan, the model with its nuisance blocks declared first.
#'
#' @param so_file The library of the model itself, <model>_model.so.
#' @param stanc_args The stanc arguments of the build (its include paths).
#' @return The path of the library NicoStan builds.
#' @noRd
fn_library_NicoStan_builds_for_a_MetaOrdDTA_model <-  function(
                    so_file,
                    stanc_args
) {

      Stan_model_file_path <-  sub("_model\\.so$", ".stan", so_file)
      ##
      ## ---- A model that calls neither reduce_sum nor map_rect: its copy with the nuisance blocks first -----------------------------
      ##
      if (!NicoStan:::R_fn_Stan_model_calls_reduce_sum_or_map_rect(Stan_model_file_path  = Stan_model_file_path,
                                                                    stanc_args            = stanc_args)) {

            return(sub("_model\\.so$", "_nuisance_first_model.so", so_file))

      }
      return(fn_library_of_the_nico_reduce_sum_version_of_a_MetaOrdDTA_model(so_file))

}
