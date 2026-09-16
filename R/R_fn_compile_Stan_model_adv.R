

#### ------- R function to compile Stan model  -----------------------------------------------------------------------------------

R_fn_compile_stan_advanced_given_file_name <- function( stan_model_file_name, 
                                                        cts,
                                                        network,
                                                        prior_only,
                                                        ##
                                                        debugging ,
                                                        force_recompile,
                                                        quiet,
                                                        compile,
                                                        ##
                                                        set_custom_CXX_CPP_flags, 
                                                        CCACHE_PATH, 
                                                        custom_cpp_user_header_file_path, ## if wish to include custom C++ files
                                                        CXX_COMPILER_PATH,
                                                        CPP_COMPILER_PATH,
                                                        MATH_FLAGS,
                                                        FMA_FLAGS,
                                                        AVX_FLAGS,
                                                        THREAD_FLAGS
) {
  
                    ##
                    ## First, run the BASIC stan compile fn but WITHOUT compiling - just to get the model names!
                    ##
                    outs_stan_basic_wo_compiling <- R_fn_compile_stan_model_basic_given_file_name(
                                                                              stan_model_file_name = stan_model_file_name,
                                                                              cts = cts,
                                                                              network = network,
                                                                              prior_only = prior_only,
                                                                              ##
                                                                              debugging = debugging,
                                                                              force_recompile = force_recompile,
                                                                              quiet = quiet,
                                                                              compile = FALSE)
                    
                    stan_model_file_name <- outs_stan_basic_wo_compiling$stan_model_file_name
                    stan_model_file_path <- outs_stan_basic_wo_compiling$stan_model_file_path
                    ##
                    pkg_root_directory       <- outs_stan_basic_wo_compiling$pkg_root_directory
                    stan_models_directory    <- outs_stan_basic_wo_compiling$stan_models_directory
                    stan_functions_directory <- outs_stan_basic_wo_compiling$stan_functions_directory
                    ##
                    stan_MA_directory       <- outs_stan_basic_wo_compiling$stan_MA_directory
                    stan_MA_prior_directory <- outs_stan_basic_wo_compiling$stan_MA_prior_directory
                    ##
                    stan_NMA_directory       <- outs_stan_basic_wo_compiling$stan_NMA_directory
                    stan_NMA_prior_directory <- outs_stan_basic_wo_compiling$stan_NMA_prior_directory
  
  
                    ## Set C++ flags: 
                    ## Use C++ 17:
                    CXX_STD <- "CXX17"
                    ##
                    if (is.null(CXX_COMPILER_PATH)) {
                      ## CXX_COMPILER <- "/opt/AMD/aocc-compiler-5.0.0/bin/clang++"
                      CXX_COMPILER <- "g++"
                      ## CXX_COMPILER <- "clang++"
                    } else { 
                      CXX_COMPILER <- CXX_COMPILER_PATH
                    }
                    
                    if (is.null(CPP_COMPILER_PATH)) {
                      ## CPP_COMPILER <- "/opt/AMD/aocc-compiler-5.0.0/bin/clang"
                      CPP_COMPILER <- "gcc"
                      # CPP_COMPILER <- "clang"
                    } else { 
                      CPP_COMPILER <- CPP_COMPILER_PATH
                    }
                    
                    ##
                    CCACHE <- CCACHE_PATH
                    # if (CCACHE_PATH == " ") {
                    #     if (.Platform$OS.type == "unix") {
                    #       CCACHE <- "/usr/bin/ccache"
                    #     } else { 
                    #       CCACHE <- CCACHE_PATH
                    #     }
                    # }
                    ##
                    CPU_BASE_FLAGS <- "-O3" ##  -march=native  -mtune=native"
                    ## Check for CPU features (i.e., FMA and/or AVX and/or AVX2 and/or AVX-512)
                    if (is.null(FMA_FLAGS)) {
                      FMA_FLAGS <- " "
                    }
                    try({  
                        try({ 
                             check_CPU_features_using_BayesMVP <- BayesMVP:::checkCPUFeatures()
                          }, silent = TRUE)
                          ## FMA:
                          if (is.null(FMA_FLAGS)) {
                             has_FMA <- 0
                            try({ 
                              has_FMA <- check_CPU_features_using_BayesMVP$has_fma
                            })
                              if (has_FMA == 1) { 
                                FMA_FLAGS <- "-mfma"
                              }
                          }
                    })
                    
                    detect_SIMD_type <- NULL
                    
                    if (is.null(AVX_FLAGS)) {
                       AVX_FLAGS <- " "
                    }
                    try({  
                        if (is.null(AVX_FLAGS)) {
                            ## AVX Flags (using AVX2 on my laptop and AVX-512 on my local HPC):
                            AVX_512_FLAGS <- "-mavx -mavx2 -mavx512f -mavx512vl -mavx512dq" ## for AVX-512 - ONLY USE IF YOUR PC SUPPORTS AVX-512 (MANY DONT!) - ELSE COMMENT OUT!
                            AVX2_FLAGS <- "-mavx -mavx2" ## for AVX2 - ONLY USE IF YOUR PC SUPPORTS AVX2 - ELSE COMMENT OUT!
                            AVX1_FLAGS <- "-mavx"
                            #
                            detect_SIMD_type <- NULL
                            AVX_FLAGS <- NULL
                            try({ 
                              detect_SIMD_type <- BayesMVP:::detect_vectorization_support()
                            }, silent = TRUE)
                            ##
                            if (detect_SIMD_type == "AVX512") { 
                              AVX_FLAGS <- AVX_512_FLAGS
                            } else if  (detect_SIMD_type == "AVX2") { 
                              AVX_FLAGS <- AVX2_FLAGS
                            } else if  (detect_SIMD_type == "AVX") {  
                              AVX_FLAGS <- AVX1_FLAGS 
                            } else { 
                              AVX_FLAGS <- " "  ## No AVX fns at all otherwise
                            }
                        }
                    })
                    ##
                    CPU_FLAGS <- paste(CPU_BASE_FLAGS, FMA_FLAGS, AVX_FLAGS)
                    ## Math flags:
                    if (is.null(MATH_FLAGS)) {
                       MATH_FLAGS <- "-fno-math-errno  -fno-signed-zeros" ####  -fno-trapping-math"
                       ##  MATH_FLAGS <- "-fno-math-errno  -fno-signed-zeros  -fno-trapping-math"
                    }
                    ##
                    if (is.null(THREAD_FLAGS)) {
                        ## THREAD_FLAGS <- "-pthread -D_REENTRANT -DSTAN_THREADS -DSTAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP"
                        THREAD_FLAGS <- "-D_REENTRANT"
                        ## THREAD_FLAGS <- "-pthread -DSTAN_THREADS"
                        ##THREAD_FLAGS <- " "
                    }
                    ##
                    OTHER_FLAGS <- " "
                    # ## OTHER_FLAGS <- "-std=c++17"
                    # OTHER_FLAGS <- paste(OTHER_FLAGS, "-fPIC")  
                    # # OTHER_FLAGS <- paste(OTHER_FLAGS, "-DNDEBUG")  
                    # # OTHER_FLAGS <- paste(OTHER_FLAGS, "-fpermissive")  
                    # OTHER_FLAGS <- paste(OTHER_FLAGS, "-DBOOST_DISABLE_ASSERTS")
                    # OTHER_FLAGS <- paste(OTHER_FLAGS, "-Wno-deprecated-declarations") 
                    # OTHER_FLAGS <- paste(OTHER_FLAGS, "-Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess") 
                    # # OTHER_FLAGS <- paste0(OTHER_FLAGHS, " -ftls-model=global-dynamic") # doesnt fix clang issue
                    # # OTHER_FLAGS <- paste(OTHER_FLAGS, "-fno-gnu-unique") # doesnt fix clang issue
                    # ##OTHER_FLAGS <- paste(OTHER_FLAGS, "-ftls-model=initial-exec")  # doesnt fix clang issue
                    # OTHER_FLAGS <-  " "
                    # ##
                    # OMP_LIB_PATH = "/opt/AMD/aocc-compiler-5.0.0/lib"
                    # OMP_FLAGS = "-fopenmp"
                    # OMP_LIB_FLAGS = "-L'OMP_LIB_PATH' -lomp 'OMP_LIB_PATH/libomp.so'"
                    # SHLIB_OPENMP_CFLAGS <- OMP_FLAGS
                    # SHLIB_OPENMP_CXXFLAGS <- SHLIB_OPENMP_CFLAGS
                    # 
                    CMDSTAN_INCLUDE_PATHS <- " "
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/tbb_2020.3/include")
                    # # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/tbb")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I src")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/src")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/rapidjson_1.1.0/")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I lib/CLI11-1.9.1/")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/eigen_3.4.0")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/boost_1.84.0")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/sundials_6.1.1/include")
                    # CMDSTAN_INCLUDE_PATHS <- paste(CMDSTAN_INCLUDE_PATHS, "-I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials")
                
                    ##
                    ##
                    ## ----- Now set "BASE_FLAGS" --------------------
                    BASE_FLAGS <- paste(CPU_FLAGS, 
                                        ## `SHLIB_OPENMP_CFLAGS`,
                                        MATH_FLAGS, 
                                        OTHER_FLAGS, 
                                        THREAD_FLAGS, 
                                        CMDSTAN_INCLUDE_PATHS)
                    
                    ## Now set standard C++/C flags. 
                    CC <-  paste(CCACHE, CPP_COMPILER)
                    CXX <- paste(CCACHE, CXX_COMPILER)
                    PKG_CPPFLAGS <- BASE_FLAGS
                    PKG_CXXFLAGS <- BASE_FLAGS
                    CPPFLAGS <- BASE_FLAGS
                    CXXFLAGS <- BASE_FLAGS
                    CFLAGS <- CPPFLAGS
                    
                    ## Linking flags to match Stan's
                    LINKER_FLAGS <- paste(
                      "-Wl,-L,\"$(CMDSTAN_PATH)/stan/lib/stan_math/lib/tbb\"",
                      "-Wl,-rpath,\"$(CMDSTAN_PATH)/stan/lib/stan_math/lib/tbb\"",
                      ##  "-lpthread",
                      "-ltbb"
                    )
                    LDFLAGS <- LINKER_FLAGS
                    
                    
                    cmdstan_cpp_flags <- list(  paste0("CC = ", CC),
                                                paste0("CXX = ", CXX),
                                                paste0("PKG_CPPFLAGS = ", PKG_CPPFLAGS),
                                                paste0("PKG_CXXFLAGS = ", PKG_CXXFLAGS),
                                                paste0("CPPFLAGS = ", CPPFLAGS),
                                                paste0("CXXFLAGS = ", CXXFLAGS),
                                                paste0("CFLAGS = ", CFLAGS), 
                                                paste0("CXX_STD = ", CXX_STD), 
                                                paste0("LDFLAGS = ", LDFLAGS))
                    
                    ## Print the C++ flags being used to compile the Stan model:
                    message(paste("BASE_FLAGS set to:"))
                    print(BASE_FLAGS)
                    ##
                    message(paste("cmdstan_cpp_flags set to:"))
                    print(cmdstan_cpp_flags)
                    
                    
                    ## --------- Compile Stan model: -----------------------------------------------------------------------
                    print(custom_cpp_user_header_file_path)
                    ##
                    try({  
                        ##
                        if (.Platform$OS.type == "unix") {
                          system(paste("chmod +x", stan_model_file_path))
                          Sys.chmod(stan_model_file_path, mode = "0755")  # Read/write/execute for owner, read/execute for others
                        }
                        ##
                        if (set_custom_CXX_CPP_flags == TRUE) {
                          
                              message((paste("aaaaaaaaaaaaaaaaaaaaaaaaa")))
                          
                              if (is.null(custom_cpp_user_header_file_path)) {
                                    
                                    message((paste("bbbbbbbbbbbbbbbbbbbb")))
                                    stan_model_obj <- cmdstanr::cmdstan_model(  stan_model_file_path,
                                                                       force_recompile = force_recompile,
                                                                       quiet = quiet,
                                                                       cpp_options = cmdstan_cpp_flags,
                                                                       include_paths = stan_functions_directory)
                                
                              } else { 
                                
                                    stan_model_obj <- cmdstanr::cmdstan_model(  stan_model_file_path,
                                                                       force_recompile = force_recompile,
                                                                       quiet = quiet,
                                                                       user_header = custom_cpp_user_header_file_path,
                                                                       cpp_options = cmdstan_cpp_flags,
                                                                       include_paths = stan_functions_directory)
                                
                              }
                          
                        } else { 
                          
                              if (is.null(custom_cpp_user_header_file_path)) {
                                    
                                    stan_model_obj <- cmdstanr::cmdstan_model(  stan_model_file_path,
                                                                     force_recompile = force_recompile,
                                                                     quiet = quiet,
                                                                     include_paths = stan_functions_directory)
                                
                                                                
                              } else { 
                                    
                                    stan_model_obj <- cmdstanr::cmdstan_model(  stan_model_file_path,
                                                                     force_recompile = force_recompile,
                                                                     quiet = quiet,
                                                                     user_header = custom_cpp_user_header_file_path,
                                                                     include_paths = stan_functions_directory)
                                
                                                               
                              }
                          
                        }
                    })
                    
                    ## Store flags in lists for output:
                    FLAGS_STANDARD_MACROS <- list(  CXX_STD = CXX_STD,
                                                    CC = CC,
                                                    CXX = CXX,
                                                    PKG_CPPFLAGS = PKG_CPPFLAGS,
                                                    PKG_CXXFLAGS = PKG_CXXFLAGS,
                                                    CPPFLAGS = CPPFLAGS,
                                                    CXXFLAGS = CXXFLAGS,
                                                    CFLAGS = CFLAGS,
                                                    LDFLAGS = LDFLAGS)
                    ##
                    FLAGS_CUSTOM_MACROS <- list(BASE_FLAGS = BASE_FLAGS)
                    
                    
                    ##
                    ## output list:
                    ##
                    return(list( FLAGS_STANDARD_MACROS = FLAGS_STANDARD_MACROS,
                                 FLAGS_CUSTOM_MACROS = FLAGS_CUSTOM_MACROS,
                                 cmdstan_cpp_flags = cmdstan_cpp_flags,
                                 ##
                                 stan_model_obj = stan_model_obj, 
                                 ##
                                 stan_model_file_name = stan_model_file_name,
                                 stan_model_file_path = stan_model_file_path,
                                 ##
                                 pkg_root_directory = pkg_root_directory,
                                 stan_models_directory = stan_models_directory,
                                 stan_functions_directory = stan_functions_directory,
                                 ##
                                 stan_MA_directory = stan_MA_directory,
                                 stan_MA_prior_directory = stan_MA_prior_directory,
                                 ##
                                 stan_NMA_directory = stan_NMA_directory,
                                 stan_NMA_prior_directory = stan_NMA_prior_directory))
              
                    
}

    
    
   
    
    
    
    
    
    
    
    
    
    
    
    

    




 