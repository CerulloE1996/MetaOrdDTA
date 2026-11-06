## ===============================================================================================================
## Build the four NicoStan fused AVX versions of MetaOrdDTA's MA models as cmdstanr executables (NicoSIMD's C++
## through user_header); an executable already built is skipped.
## Writes: sim_study/cmdstanr_executables_of_the_NicoStan_fused_AVX_MA_models_<machine>/<model>_NicoStan_fused_AVX
## ===============================================================================================================
##
##
{
    MetaOrdDTA_folder      <-  path.expand("~/Documents/Work/PhD_work/R_packages/MetaOrdDTA")
    models_folder          <-  file.path(MetaOrdDTA_folder, "inst", "stan_models", "stan_models_MA")
    functions_folder       <-  file.path(MetaOrdDTA_folder, "sim_study", "check_of_the_cutpoint_Jacobian_17_Oct",
                                        "library_with_the_corrected_cutpoint_Jacobian", "MetaOrdDTA",
                                        "stan_functions")
    NicoSIMD_user_header   <-  path.expand(paste0("~/Documents/Work/PhD_work/R_packages/NicoStan/inst/NicoSIMD/",
                                                 "inst/include/NicoSIMD/stan_external_functions_of_NicoSIMD.hpp"))
    models                 <-  c("DTA_MA_ord_HSROC_RANDthr_kappa",
                                "DTA_MA_ord_biv_RANDthr_kappa", "DTA_MA_ord_biv_FIXEDthr")
    CPU_has_AVX512         <-  any(grepl("\\bavx512f\\b", readLines("/proc/cpuinfo", warn = FALSE)))
    machine                <-  if (CPU_has_AVX512) "local_HPC" else "laptop"
    compiler_flags         <-  if (CPU_has_AVX512) {
                                  "-O3 -march=native -mtune=native -mfma -mavx2 -mavx512f -mavx512vl -mavx512dq"
                              } else {
                                  "-O3 -march=native -mtune=native -mfma -mavx2"
                              }
    executables_folder     <-  file.path(MetaOrdDTA_folder, "sim_study",
                                        paste0("cmdstanr_executables_of_the_NicoStan_fused_AVX_MA_models_",
                                               machine))
}
##
##
##
dir.create(executables_folder, showWarnings = FALSE, recursive = TRUE)
for (model in models) {
    stan_file  <-  file.path(executables_folder, paste0(model, "_NicoStan_fused_AVX.stan"))
    executable <-  file.path(executables_folder, paste0(model, "_NicoStan_fused_AVX"))
    if (file.exists(executable) && (file.mtime(executable) > file.mtime(file.path(models_folder,
                                                                         basename(stan_file))))) {
        message(paste0("Already built, skipped: ", executable))
        next
    }
    file.copy(file.path(models_folder, basename(stan_file)), stan_file, overwrite = TRUE)
    cpp_options <-  list(CXXFLAGS_OPTIM = compiler_flags)
    if (!CPU_has_AVX512) cpp_options$CPPFLAGS_OPTIM <-  "-DBAYESMVP_FORCE_AVX2"
    time_start <-  Sys.time()
    model_compiled <-  cmdstanr::cmdstan_model(  stan_file       = stan_file,
                                                include_paths   = functions_folder,
                                                user_header     = NicoSIMD_user_header,
                                                ## (cmdstanr adds --allow-undefined itself with a user_header)
                                                cpp_options     = cpp_options,
                                                force_recompile = TRUE)
    message(paste0("Built: ", model_compiled$exe_file(), " (",
                   round(as.numeric(difftime(Sys.time(), time_start, units = "secs")), 1), " s)"))
}






















