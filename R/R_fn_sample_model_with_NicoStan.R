##
## ==== R_fn_sample_model_with_NicoStan.R ========================================================================
##
## ---- MetaOrd_model$sample(algorithm = "NicoStan"): the model sampled with NicoStan, with its version that uses
##      NicoSIMD's fused AVX functions (<model>_NicoStan_fused_AVX.stan): data, library, fit and draws.
##
R_fn_sample_model_with_NicoStan  <-  function(  stan_data_list,
                                                stan_model_file_name,
                                                n_chains,
                                                init_lists_per_chain,
                                                seed,
                                                n_burnin,
                                                n_iter,
                                                NicoStan_args) {
        ##
        ##
        ## ---- (1) The model, its NicoStan version with the fused AVX functions, and the data: ------------------
        ##
        start_of_the_preparation <-  Sys.time()
        prepare_only <-  isTRUE(NicoStan_args$prepare_only)
        NicoStan_args$prepare_only <-  NULL
        name_of_the_model <-  sub("(_NicoStan_fused_AVX)?\\.stan$", "", basename(stan_model_file_name))
        ## (the MA models in stan_models_MA, the NMA models in stan_models_NMA)
        folder_of_the_models <-  if (grepl("^DTA_NMA_", name_of_the_model)) {
              "stan_models_NMA"
        } else {
              "stan_models_MA"
        }
        Stan_file_of_the_fused_version <-  NicoStan_args$Stan_file_of_the_fused_version
        if (is.null(Stan_file_of_the_fused_version)) {
              Stan_file_of_the_fused_version <-  system.file(
                    "stan_models", folder_of_the_models,
                    paste0(name_of_the_model, "_NicoStan_fused_AVX.stan"), package = "MetaOrdDTA")
        }
        if (!nzchar(Stan_file_of_the_fused_version)) {
              stop(paste0("algorithm = \"NicoStan\": no NicoStan version (", name_of_the_model,
                          "_NicoStan_fused_AVX.stan) of the model ", name_of_the_model, " in ",
                          system.file("stan_models", folder_of_the_models, package = "MetaOrdDTA")))
        }
        ## (optional: NicoStan_args$Stan_file_of_the_fused_version, another file of the fused version)
        if (!is.null(NicoStan_args$Stan_file_of_the_fused_version)) {
              Stan_file_of_the_fused_version <-  NicoStan_args$Stan_file_of_the_fused_version
              name_of_the_model <-  sub("(_NicoStan_fused_AVX)?\\.stan$", "",
                                        basename(Stan_file_of_the_fused_version))
        }
        stopifnot(file.exists(Stan_file_of_the_fused_version))
        fused_AVX <-  grepl("_NicoStan_fused_AVX[.]stan$", basename(Stan_file_of_the_fused_version))
        math_label <-  if (fused_AVX) "fused AVX" else "native Stan::math"
        NicoStan_args$Stan_file_of_the_fused_version <-  NULL
        folder_of_the_Stan_functions <-  system.file("stan_functions", package = "MetaOrdDTA")
        user_header_of_NicoSIMD <-  system.file("include", "NicoSIMD", "stan_external_functions_of_NicoSIMD.hpp",
                                                package = "NicoSIMD")
        if (!nzchar(user_header_of_NicoSIMD)) {
              stop("algorithm = \"NicoStan\": the NicoSIMD package (NicoStan's AVX functions) is not installed")
        }
        ##
        ## (the 0 / 1 switches as numbers, as cmdstanr writes them)
        Stan_data_list <-  stan_data_list
        for (name in names(Stan_data_list)) {
              if (is.logical(Stan_data_list[[name]])) {
                    Stan_data_list[[name]] <-  as.integer(Stan_data_list[[name]])
              }
        }
        ## (NicoStan_args$N_chunks_of_auto_chunking: the number of chunks of AUTO-CHUNKING, read by the library
        ##  from the data; NicoStan puts 1 in the data when it is not given)
        if (!is.null(NicoStan_args$N_chunks_of_auto_chunking)) {
              Stan_data_list$N_chunks_of_auto_chunking <-  as.integer(NicoStan_args$N_chunks_of_auto_chunking)
        }
        NicoStan_args$N_chunks_of_auto_chunking <-  NULL
        ## (a field of one value declared as a vector, row_vector or array: an array of length 1)
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
        ## ---- (2) The library, in MetaOrdDTA's cache folder, built again only when its files change: ----------
        ##
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
        ## (the lock lives in the cache folder, so the folder must exist before the first lock is taken)
        dir.create(cache_folder, recursive = TRUE, showWarnings = FALSE)
        lock_dir <-  file.path(cache_folder, paste0(name_of_the_model, "_", SIMD_label, "_build_lock"))
        ## (one lock per model and backend: the first worker builds the library and the others wait; a lock older
        ##  than 15 minutes was left by a killed process and is removed)
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
        Stan_model_file_path <-  file.path(cache_folder,
                                           paste0(name_of_the_model,
                                                  if (fused_AVX) "_NicoStan_fused_" else "_NicoStan_",
                                                  SIMD_label, ".stan"))
        so_file <-  sub("\\.stan$", "_model.so", Stan_model_file_path)
        ##
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
              ## (the earlier libraries kept beside them: the model's own, and the one NicoStan built from its
              ##  copy with the nuisance blocks declared first, which NicoStan otherwise refuses as older than
              ##  that copy)
              for (library_built_before in c(so_file, sub("_model\\.so$", "_nuisance_first_model.so", so_file))) {
                    if (file.exists(library_built_before)) {
                          file.rename(library_built_before,
                                      paste0(library_built_before, ".built_before_",
                                             format(Sys.time(), "%Y_%m_%d_%H%M%S")))
                    }
              }
        }
        make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                        paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
        if (!CPU_has_AVX512) make_args <-  c(make_args, "CPPFLAGS_OPTIM=-DBAYESMVP_FORCE_AVX2")
        stanc_args <-  c(paste0("--include-paths=", folder_of_the_copied_Stan_functions), "--allow-undefined")
        ##
        ## (NicoStan_args$nuisance_params go to Nico_model$new, the other NicoStan_args to $sample)
        nuisance_params <-  NicoStan_args$nuisance_params
        NicoStan_args$nuisance_params <-  NULL
        ## (the package names the nuisance blocks of its own models)
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
                                         name_of_a_block     = block,
                                         name_of_the_model   = name_of_the_model), "\033[0m"))
              }
        }
        ## (the parameter families of interest, NicoStan_args$interest_only: the tau criterion and the min ESS use
        ##  only these; the package names them when they are not given)
        if (is.null(NicoStan_args$interest_only)) {
              NicoStan_args$interest_only <-  fn_interest_only_params_of_a_MetaOrdDTA_model(
                    name_of_the_model        = name_of_the_model,
                    model_parameterisation   = if (grepl("^DTA_NMA_", name_of_the_model)) "NMA" else
                                                     if (grepl("Jones", name_of_the_model)) "Jones" else
                                                     if (grepl("ord_HSROC", name_of_the_model)) "ordinal_HSROC" else
                                                     if (grepl("ord_biv", name_of_the_model)) "ord_biv" else "Bivariate",
                    random_thresholds        = grepl("RANDthr", name_of_the_model))
        }
        message(paste0("\033[36m", "parameter families of interest (tau criterion and min ESS): ",
                       paste(NicoStan_args$interest_only, collapse = ", "), "\033[0m"))
        ## (sample_nuisance = TRUE needs the nuisance blocks named in NicoStan_args$nuisance_params)
        # ## (sample_nuisance = TRUE without nuisance_params: the first parameter block of the model is the nuisance
        # ##  block, and the model file is used as it is)
        sample_nuisance <-  isTRUE(NicoStan_args$sample_nuisance) || !is.null(nuisance_params)
        if (isTRUE(NicoStan_args$sample_nuisance) && is.null(nuisance_params)) {
              stop(paste0("sample_model_with_NicoStan: sample_nuisance = TRUE needs the nuisance blocks of ",
                          name_of_the_model, " NAMED in NicoStan_args$nuisance_params (e.g. c(\"beta_z\")); ",
                          "NicoStan never takes the first parameter block as the nuisance block."))
        }
        NicoStan_args$sample_nuisance <-  NULL
        ## (optional: NicoStan_args$n_nuisance_override, how many first unconstrained values are nuisance)
        n_nuisance_override <-  NicoStan_args$n_nuisance_override
        NicoStan_args$n_nuisance_override <-  NULL
        ## (named nuisance blocks declared first, with n_nuisance_override: held in the model's own library,
        ##  without the nico_reduce_sum writer; this takes bounded blocks such as C_raw <lower=-15, upper=2.5>)
        ##
        ## ---- The named nuisance blocks: NicoStan's nico_reduce_sum version only for a model that calls reduce_sum ------------------
        ##
        ## (without n_nuisance_override the blocks go to NicoStan by name: a model that calls reduce_sum or map_rect
        ##  gets its nico_reduce_sum version; a model that calls neither keeps them in its own parameters block,
        ##  declared first (<model>_nuisance_first.stan), and NicoStan samples them as its nuisance block)
        nuisance_params_named_for_NicoStan <-  if (is.null(n_nuisance_override)) {
              nuisance_params
        } else {
              NULL
        }
        if (!is.null(n_nuisance_override)) {
              message(paste0("\033[36m", "The nuisance blocks ", paste(nuisance_params, collapse = ", "),
                             " (", n_nuisance_override, " values, declared first in ",
                             basename(Stan_file_of_the_fused_version), ") held in the model's own library",
                             "\033[0m"))
        }
        ## (NicoStan_args$auto_chunking: NicoStan's AUTO-CHUNKING of the model, TRUE unless it is FALSE; given to
        ##  a NicoStan whose Nico_model$new takes it)
        auto_chunking <-  !isFALSE(NicoStan_args$auto_chunking)
        NicoStan_args$auto_chunking <-  NULL
        NicoStan_takes_auto_chunking <-  "auto_chunking" %in%
              names(formals(NicoStan::Nico_model$public_methods$initialize))
        if (!NicoStan_takes_auto_chunking && !auto_chunking) {
              stop("NicoStan_args$auto_chunking = FALSE: this NicoStan has no AUTO-CHUNKING option")
        }
        ## (with AUTO-CHUNKING the data always hold its number of chunks: 1 unless it is given)
        if (NicoStan_takes_auto_chunking && auto_chunking &&
            is.null(Stan_data_list$N_chunks_of_auto_chunking)) {
              Stan_data_list$N_chunks_of_auto_chunking <-  1L
        }
        ## (a library already built: the lock is released before NicoStan loads it; only a build holds it; with
        ##  AUTO-CHUNKING, NicoStan checks its own receipt and may build, so the lock is kept)
        library_to_load <-  if (!is.null(nuisance_params_named_for_NicoStan)) {
              fn_library_NicoStan_builds_for_a_MetaOrdDTA_model(so_file     = so_file,
                                                                stanc_args  = stanc_args)
        } else {
              so_file
        }
        if (!something_has_changed && file.exists(library_to_load) &&
            !(NicoStan_takes_auto_chunking && auto_chunking)) {
              unlink(lock_dir, recursive = TRUE)
              this_fit_holds_the_lock <-  FALSE
        }
        start_of_the_build <-  Sys.time()
        arguments_of_the_model <-  list( Model_type            = "Stan",
                                         Stan_data_list        = Stan_data_list,
                                         Stan_model_file_path  = Stan_model_file_path,
                                         Stan_cpp_user_header  = user_header_of_NicoSIMD,
                                         stanc_args            = stanc_args,
                                         make_args             = make_args,
                                         nuisance_params       = nuisance_params_named_for_NicoStan,
                                         n_nuisance_override   = n_nuisance_override,
                                         sample_nuisance       = sample_nuisance)
        if (NicoStan_takes_auto_chunking) {
              arguments_of_the_model$auto_chunking <-  auto_chunking
        }
        model_NicoStan <-  do.call(NicoStan::Nico_model$new, arguments_of_the_model)
        seconds_of_the_build <-  as.numeric(difftime(Sys.time(), start_of_the_build, units = "secs"))
        # message(paste0("\033[36m", "NicoStan version of ", name_of_the_model, " (fused AVX, ", SIMD_label, "): ",
        #                so_file, " (", formatC(seconds_of_the_build, format = "f", digits = 1), " s)", "\033[0m"))
        ## (the library NicoStan loads when a nuisance block is named: its nico_reduce_sum version, or the library
        ##  of <model>_nuisance_first.stan for a model that calls neither reduce_sum nor map_rect; else so_file)
        library_built <-  if (!is.null(nuisance_params_named_for_NicoStan)) {
              fn_library_NicoStan_builds_for_a_MetaOrdDTA_model(so_file     = so_file,
                                                                stanc_args  = stanc_args)
        } else {
              so_file
        }
        ## (the library NicoStan itself loaded, e.g. the AUTO-CHUNKING version of the nuisance-first library)
        if (!is.null(model_NicoStan$init_object$model_so_file)) {
              library_built <-  model_NicoStan$init_object$model_so_file
        }
        message(paste0("\033[36m", "NicoStan version of ", name_of_the_model,
                       " (", math_label, ", ", SIMD_label, "): ",
                       library_built, " (", formatC(seconds_of_the_build, format = "f", digits = 1), " s)",
                       "\033[0m"))
        ## (a missing library stops here)
        if (!file.exists(library_built)) {
              stop(paste0("sample_model_with_NicoStan: the build of ", name_of_the_model,
                          " finished but the library ", library_built, " is not there"))
        }
        if (this_fit_holds_the_lock) {
              unlink(lock_dir, recursive = TRUE)
              this_fit_holds_the_lock <-  FALSE
        }
        if (prepare_only) {
              return(list(NicoStan_model_obj = model_NicoStan,
                          Stan_data_list = Stan_data_list,
                          so_file = library_built))
        }
        ##
        ## ---- (3) The fit: -----------------------------------------------------------------------------------
        ##
        start_of_the_fit <-  Sys.time()
        ## (n_chains burn-in and sampling chains, unless NicoStan_args gives n_chains_burnin / n_chains_sampling;
        ##  the initial values are those of the burn-in chains)
        n_chains_burnin <-  if (is.null(NicoStan_args$n_chains_burnin)) n_chains else
                                  NicoStan_args$n_chains_burnin
        n_chains_sampling <-  if (is.null(NicoStan_args$n_chains_sampling)) n_chains else
                                    NicoStan_args$n_chains_sampling
        NicoStan_args$n_chains_burnin <-  NULL
        NicoStan_args$n_chains_sampling <-  NULL
        ## (the nuisance values kept in memory for param_constrain: the random-cutpoint Se / Sp need C_raw)
        if (!is.null(NicoStan_args$n_nuisance_to_track)) {
              stop(paste0("sample_model_with_NicoStan: n_nuisance_to_track is now nuisance_to_track; give ",
                          "nuisance_to_track = c(\"", paste(nuisance_params, collapse = "\", \""),
                          "\") to keep the named nuisance blocks, or \"none\" to keep none."))
        }
        nuisance_to_track <-  if (!is.null(NicoStan_args$nuisance_to_track)) {
              NicoStan_args$nuisance_to_track
        } else if (!is.null(nuisance_params)) {
              nuisance_params
        } else {
              "none"
        }
        NicoStan_args$nuisance_to_track <-  NULL
        ## (NicoStan_args$compute_transformed_parameters = TRUE: the draws also hold the transformed parameters)
        compute_transformed_parameters <-  isTRUE(NicoStan_args$compute_transformed_parameters)
        NicoStan_args$compute_transformed_parameters <-  NULL
        arguments_of_the_fit <-  c(list(init_lists_per_chain  = init_lists_per_chain[seq_len(n_chains_burnin)],
                                        n_chains_burnin       = n_chains_burnin,
                                        n_chains_sampling     = n_chains_sampling,
                                        n_burnin              = n_burnin,
                                        n_iter                = n_iter,
                                        seed                  = seed,
                                        nuisance_to_track     = nuisance_to_track),
                                   NicoStan_args)
        do.call(model_NicoStan$sample, arguments_of_the_fit)
        seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), start_of_the_fit, units = "secs"))
        ##
        ## ---- (4) The draws, the divergent transitions and the summaries of Se_baseline and Sp_baseline: -----
        ##
        summary_NicoStan <-  model_NicoStan$summary( compute_main_params             = TRUE,
                                                     compute_transformed_parameters  =
                                                           compute_transformed_parameters,
                                                     compute_generated_quantities    = TRUE,
                                                     save_log_lik_trace              = FALSE)
        draws_array <-  summary_NicoStan$get_posterior_draws()
        divergences <-  summary_NicoStan$get_divergences()
        ##
        ## (Se_baseline and Sp_baseline from NicoStan's own summary of the fit)
        summary_of_every_output <-  rbind(summary_NicoStan$get_summary_main(),
                                          summary_NicoStan$get_summary_transformed(),
                                          summary_NicoStan$get_summary_generated_quantities())
        fn_rows_of_one_variable <-  function(name_of_the_variable) {

              return(summary_of_every_output[grepl(paste0("^", name_of_the_variable, "\\["),
                                                   summary_of_every_output$parameter), ])

        }
        summary_of_Se <-  fn_rows_of_one_variable("Se_baseline")
        summary_of_Sp <-  fn_rows_of_one_variable("Sp_baseline")
        message(paste0("\033[36m", "NicoStan fit of ", name_of_the_model, ": ", n_chains_burnin, " burn-in and ",
                       n_chains_sampling, " sampling chains, ", n_burnin,
                       " burn-in and ", n_iter, " sampling iterations, ",
                       formatC(seconds_of_the_fit, format = "f", digits = 1), " s; divergent transitions: ",
                       paste0(unlist(divergences), collapse = ", "), "\033[0m"))
        ##
        return(list( algorithm                               = "NicoStan",
                     NicoStan_model_obj                      = model_NicoStan,
                     NicoStan_summary_obj                    = summary_NicoStan,
                     draws_array                             = draws_array,
                     divergences                             = divergences,
                     summary_of_Se_baseline                  = summary_of_Se,
                     summary_of_Sp_baseline                  = summary_of_Sp,
                     ## (NicoStan's own summary of every main parameter, transformed parameter and generated
                     ##  quantity)
                     summary_of_every_output                 = summary_of_every_output,
                     so_file                                 = so_file,
                     ## (the library NicoStan loaded, as NicoStan records it: e.g. the nuisance-first,
                     ##  nico_reduce_sum or AUTO-CHUNKING version)
                     library_used                            =
                           if (!is.null(model_NicoStan$init_object$model_so_file))
                                 model_NicoStan$init_object$model_so_file else library_built,
                     seconds_of_the_fit                      = seconds_of_the_fit,
                     seconds_of_the_preparation_and_the_fit  =
                           as.numeric(difftime(Sys.time(), start_of_the_preparation, units = "secs"))))
##
}






















