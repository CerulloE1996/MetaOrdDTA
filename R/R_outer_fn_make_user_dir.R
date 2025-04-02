 


#' .make_user_dir
#' @keywords internal
#' @export
.make_user_dir <- function(libname, 
                           pkgname) {
  
  
              ## Sys.getenv("USERPROFILE")
             ##  Sys.getenv("R_USER") ## This is BAD on windows dont use this !! (often points e.g. to OneDrive even if OneDrive uninstalled!!)

              os <- .Platform$OS.type
              
              if (os == "unix") { 
                user_root_dir <- Sys.getenv("PWD")
              } else if (os == "windows") { 
                user_root_dir <- Sys.getenv("USERPROFILE")
              }
              
              cat("user_root_dir =", user_root_dir, "\n")
              ##
              user_MetaOrdinal_dir <- file.path(user_root_dir, "MetaOrdinal")
              ##
              if (!dir.exists(user_MetaOrdinal_dir)) {
                dir.create(user_MetaOrdinal_dir)
              }
              ##
              cat("user_MetaOrdinal_dir =", user_MetaOrdinal_dir, "\n")
              
              ## -------------- USER Examples dir:
              user_examples_dir <- file.path(user_MetaOrdinal_dir, "examples")
              cat("user_examples_dir =", user_examples_dir, "\n")
              
              system_examples_dir <- system.file("examples", package = "MetaOrdinal")
              cat("system_examples_dir =", system_examples_dir, "\n")
              
              if (!dir.exists(user_examples_dir)) {
                dir.create(user_examples_dir)
              }
              ## Copy the entire examples directory:
              file.copy(
                from = list.files(system_examples_dir, full.names = TRUE),
                to = user_examples_dir,
                recursive = TRUE
              )
              
              # ## -------------- USER src dir:
              # user_src_dir <- file.path(user_MetaOrdinal_dir, "src")
              # cat("user_src_dir =", user_src_dir, "\n")
              # 
              # system_src_dir <- system.file("MetaOrdinal/src", package = "MetaOrdinal")
              # cat("system_src_dir =", system_src_dir, "\n")
              # 
              # if (!dir.exists(user_src_dir)) {
              #   dir.create(user_src_dir)
              # }
              # {
              #   ## Copy the entire src directory:
              #   file.copy(
              #     from = list.files(system_src_dir, full.names = TRUE),
              #     to = user_src_dir,
              #     recursive = TRUE
              #   )
              # }
              
              ## -------------- USER stan_models dir:
              user_stan_models_dir <- file.path(user_MetaOrdinal_dir, "stan_models")
              cat("user_stan_models_dir =", user_stan_models_dir, "\n")
              
              system_stan_models_dir <- system.file("MetaOrdinal/inst/stan_models/", package = "MetaOrdinal")
              cat("system_stan_models_dir =", system_stan_models_dir, "\n")
              
              if (!dir.exists(user_stan_models_dir)) {
                dir.create(user_stan_models_dir)
              }
              {
                ## Copy the entire src directory:
                file.copy(
                  from = list.files(system_stan_models_dir, full.names = TRUE),
                  to = user_stan_models_dir,
                  recursive = TRUE
                )
              }
          
              # }
              
              
}
