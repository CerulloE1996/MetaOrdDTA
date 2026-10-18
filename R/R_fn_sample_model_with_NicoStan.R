##
## ==== R_fn_sample_model_with_NicoStan.R ========================================================================
##
## ---- MetaOrd_model$sample(algorithm = "NicoStan"): the model sampled with NicoStan rather than cmdstanr, with
##      the NicoStan version of the model whose likelihood and priors use NicoSIMD's fused AVX functions
##      (<model>_NicoStan_fused_AVX.stan, next to the model in the MetaOrdDTA that is loaded: the MA models
##      DTA_MA_Xu_FIXEDthr, DTA_MA_Xu_RANDthr_kappa, DTA_MA_Gat_FIXEDthr and DTA_MA_Gat_RANDthr_kappa, and the NMA
##      model DTA_NMA_Nyaga_Xu_FIXEDthr once its NicoStan version is in stan_models_NMA):
##        (1) the same data list and priors as the cmdstanr path (its 0 / 1 switches as numbers, and a field of
##            one value that the model's data block declares as a vector, row_vector or array written as an
##            array of length 1, as cmdstanr writes it; as in sim_study/check_the_NicoStan_fused_AVX_version_of_
##            a_MetaOrdDTA_MA_model_against_the_model_on_Stan_math.R);
##        (2) the library, built once on each machine (AVX-512 where the CPU has it, otherwise AVX2) in
##            MetaOrdDTA's cache folder, one folder for each MetaOrdDTA library (e.g. each copy of the simulation
##            study), with a copy of its Stan functions; built again only when the Stan file or one of its
##            functions has changed;
##        (3) the fit with NicoStan: the same initial values, n_chains burn-in and sampling chains, n_burnin and
##            n_iter, and NicoStan's own defaults otherwise (its arguments can be given with NicoStan_args);
##            these models have no reduce_sum(), so NicoStan runs one thread per chain, without chunks;
##        (4) the draws, the divergent transitions, and the summaries of Se_baseline and Sp_baseline with
##            BayesMVP's C++ functions (Rcpp_compute_chain_stats, Rcpp_compute_MCMC_diagnostics).
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
        name_of_the_model <-  sub("\\.stan$", "", basename(stan_model_file_name))
        ## (the MA models in stan_models_MA, the NMA models in stan_models_NMA)
        folder_of_the_models <-  if (grepl("^DTA_NMA_", name_of_the_model)) {
              "stan_models_NMA"
        } else {
              "stan_models_MA"
        }
        Stan_file_of_the_fused_version <-  system.file("stan_models", folder_of_the_models,
                                                       paste0(name_of_the_model, "_NicoStan_fused_AVX.stan"),
                                                       package = "MetaOrdDTA")
        if (!nzchar(Stan_file_of_the_fused_version)) {
              stop(paste0("algorithm = \"NicoStan\": no NicoStan version (", name_of_the_model,
                          "_NicoStan_fused_AVX.stan) of the model ", name_of_the_model, " in ",
                          system.file("stan_models", folder_of_the_models, package = "MetaOrdDTA")))
        }
        ## (optional: NicoStan_args$Stan_file_of_the_fused_version, another file of the fused version (e.g. with
        ##  its parameter blocks in another order), built and cached under its own name)
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
              stop("algorithm = \"NicoStan\": the NicoSIMD package (NicoStan's AVX functions) is not installed")
        }
        ##
        ## (the 0 / 1 switches as numbers: cmdstanr writes TRUE / FALSE as 1 / 0, the JSON for BridgeStan keeps
        ##  true / false, which Stan's int data refuse)
        Stan_data_list <-  stan_data_list
        for (name in names(Stan_data_list)) {
              if (is.logical(Stan_data_list[[name]])) {
                    Stan_data_list[[name]] <-  as.integer(Stan_data_list[[name]])
              }
        }
        ## (a field of one value declared as a vector, row_vector or array in the data block: an array of
        ##  length 1)
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
        ## ---- (2) The library, in MetaOrdDTA's cache folder (built again only when something it is built from
        ##          has changed): ------------------------------------------------------------------------------
        ##
        CPU_has_AVX512 <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
        SIMD_label <-  if (CPU_has_AVX512) "AVX512" else "AVX2"
        compiler_flags <-  if (CPU_has_AVX512) {
              "-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq"
        } else {
              "-O3 -march=native -mtune=native -mfma -mavx2"
        }
        library_of_MetaOrdDTA <-  basename(dirname(system.file(package = "MetaOrdDTA")))
        cache_folder <-  file.path(tools::R_user_dir("MetaOrdDTA", which = "cache"), "NicoStan_fused_AVX",
                                   library_of_MetaOrdDTA)
        folder_of_the_copied_Stan_functions <-  file.path(cache_folder, "stan_functions")
        dir.create(folder_of_the_copied_Stan_functions, recursive = TRUE, showWarnings = FALSE)
        Stan_model_file_path <-  file.path(cache_folder, paste0(name_of_the_model, "_NicoStan_fused_", SIMD_label,
                                                                ".stan"))
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
              if (file.exists(so_file)) {
                    ## (the library built from the earlier files kept beside it, never deleted)
                    file.rename(so_file, paste0(so_file, ".built_before_", format(Sys.time(), "%Y_%m_%d_%H%M%S")))
              }
        }
        make_args <-  c("STAN_THREADS=true", "PRECOMPILED_HEADERS=false",
                        paste0("CXXFLAGS_OPTIM=", shQuote(compiler_flags)))
        if (!CPU_has_AVX512) make_args <-  c(make_args, "CPPFLAGS_OPTIM=-DBAYESMVP_FORCE_AVX2")
        stanc_args <-  c(paste0("--include-paths=", folder_of_the_copied_Stan_functions), "--allow-undefined")
        ##
        ## (the nuisance blocks, given in NicoStan_args$nuisance_params, go to Nico_model$new for diffusion HMC;
        ##  the other NicoStan_args go to $sample)
        nuisance_params <-  NicoStan_args$nuisance_params
        NicoStan_args$nuisance_params <-  NULL
        ## (sample_nuisance = TRUE without nuisance_params: the first parameter block of the model is the nuisance
        ##  block, and the model file is used as it is)
        sample_nuisance <-  isTRUE(NicoStan_args$sample_nuisance) || !is.null(nuisance_params)
        NicoStan_args$sample_nuisance <-  NULL
        ## (optional: NicoStan_args$n_nuisance_override, the number of the model's first unconstrained values held
        ##  as the nuisance block, e.g. two blocks declared one after the other at the start)
        n_nuisance_override <-  NicoStan_args$n_nuisance_override
        NicoStan_args$n_nuisance_override <-  NULL
        start_of_the_build <-  Sys.time()
        model_NicoStan <-  NicoStan::Nico_model$new( Model_type            = "Stan",
                                                     Stan_data_list        = Stan_data_list,
                                                     Stan_model_file_path  = Stan_model_file_path,
                                                     Stan_cpp_user_header  = user_header_of_NicoSIMD,
                                                     stanc_args            = stanc_args,
                                                     make_args             = make_args,
                                                     nuisance_params       = nuisance_params,
                                                     n_nuisance_override   = n_nuisance_override,
                                                     sample_nuisance       = sample_nuisance)
        seconds_of_the_build <-  as.numeric(difftime(Sys.time(), start_of_the_build, units = "secs"))
        message(paste0("\033[36m", "NicoStan version of ", name_of_the_model, " (fused AVX, ", SIMD_label, "): ",
                       so_file, " (", formatC(seconds_of_the_build, format = "f", digits = 1), " s)", "\033[0m"))
        ##
        ## ---- (3) The fit: -----------------------------------------------------------------------------------
        ##
        start_of_the_fit <-  Sys.time()
        ## (n_chains burn-in and sampling chains, unless NicoStan_args gives n_chains_burnin and / or
        ##  n_chains_sampling, e.g. 4 burn-in chains and 180 sampling chains; the initial values are those of the
        ##  burn-in chains)
        n_chains_burnin <-  if (is.null(NicoStan_args$n_chains_burnin)) n_chains else
                                  NicoStan_args$n_chains_burnin
        n_chains_sampling <-  if (is.null(NicoStan_args$n_chains_sampling)) n_chains else
                                    NicoStan_args$n_chains_sampling
        NicoStan_args$n_chains_burnin <-  NULL
        NicoStan_args$n_chains_sampling <-  NULL
        arguments_of_the_fit <-  c(list(init_lists_per_chain  = init_lists_per_chain[seq_len(n_chains_burnin)],
                                        n_chains_burnin       = n_chains_burnin,
                                        n_chains_sampling     = n_chains_sampling,
                                        n_burnin              = n_burnin,
                                        n_iter                = n_iter,
                                        seed                  = seed),
                                   NicoStan_args)
        do.call(model_NicoStan$sample, arguments_of_the_fit)
        seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), start_of_the_fit, units = "secs"))
        ##
        ## ---- (4) The draws, the divergent transitions and the summaries of Se_baseline and Sp_baseline: -----
        ##
        summary_NicoStan <-  model_NicoStan$summary( compute_main_params             = TRUE,
                                                     compute_transformed_parameters  = FALSE,
                                                     compute_generated_quantities    = TRUE,
                                                     save_log_lik_trace              = FALSE)
        draws_array <-  summary_NicoStan$get_posterior_draws()
        divergences <-  summary_NicoStan$get_divergences()
        ##
        fn_summary_of_draws_with_BayesMVP <-  function(name_of_the_variable) {
              index_of_each_element <-  grep(paste0("^", name_of_the_variable, "\\["), dimnames(draws_array)[[3]])
              draws_as_list_of_matrices <-  list()
              for (p in seq_along(index_of_each_element)) {
                    draws_of_this_element <-  as.numeric(draws_array[, , index_of_each_element[p]])
                    draws_as_list_of_matrices[[p]] <-  matrix(draws_of_this_element, nrow = dim(draws_array)[1])
              }
              n_threads <-  min(length(index_of_each_element), 4)
              means <-  BayesMVP:::Rcpp_compute_chain_stats(draws_as_list_of_matrices, stat_type = "mean",
                                                            n_threads = n_threads)$statistics[, 1]
              SDs <-  BayesMVP:::Rcpp_compute_chain_stats(draws_as_list_of_matrices, stat_type = "sd",
                                                          n_threads = n_threads)$statistics[, 1]
              quantiles <-  BayesMVP:::Rcpp_compute_chain_stats(draws_as_list_of_matrices,
                                                                stat_type = "quantiles",
                                                                n_threads = n_threads)$statistics
              ESS <-  BayesMVP:::Rcpp_compute_MCMC_diagnostics(draws_as_list_of_matrices,
                                                               diagnostic = "split_ESS",
                                                               n_threads = n_threads)$diagnostics[, 1]
              return(data.frame(variable = dimnames(draws_array)[[3]][index_of_each_element], mean = means,
                                `2.5%` = quantiles[, 1], `50%` = quantiles[, 2], `97.5%` = quantiles[, 3],
                                sd = SDs, n_eff = ESS, check.names = FALSE))
        }
        summary_of_Se <-  fn_summary_of_draws_with_BayesMVP("Se_baseline")
        summary_of_Sp <-  fn_summary_of_draws_with_BayesMVP("Sp_baseline")
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
                     so_file                                 = so_file,
                     seconds_of_the_fit                      = seconds_of_the_fit,
                     seconds_of_the_preparation_and_the_fit  =
                           as.numeric(difftime(Sys.time(), start_of_the_preparation, units = "secs"))))
##
}






















