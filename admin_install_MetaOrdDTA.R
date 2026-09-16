



rm(list = ls())
.rs.restartR()  # In RStudio 
##
{
  os <- .Platform$OS.type
  ##
  if (os == "unix") {
    user_root_dir <- Sys.getenv("PWD")
  } else if (os == "windows") {
    user_root_dir <- Sys.getenv("USERPROFILE")
  }
  ##
  local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
  ##
  {
        ## First remove any possible package fragments:
        ## Find user_pkg_ins  tall_dir:
        user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
        print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
        ##
        ## Find pkg_install_path + pkg_temp_install_path:
        pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
        pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
        ##
        ## Remove any (possible) MetaOrdDTA package fragments:
        remove.packages("MetaOrdDTA")
        unlink(pkg_install_path, recursive = TRUE, force = TRUE)
        unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
  }
  ##
  ## ---- ACTUAL (LOCAL) INSTALL:
  ##
  devtools::clean_dll(local_pkg_dir)
  devtools::document(local_pkg_dir)
  ## sroxygen2::roxygenize(local_pkg_dir)
  ##
  ## ---- Install:
  ##
  devtools::clean_dll(local_pkg_dir)
  devtools::install(local_pkg_dir,
                    ## build_vignettes = FALSE,
                    upgrade = "never",
                    quick = TRUE)
  ##
  ## ---- May need to restart R:
  ##
  .rs.restartR() # In RStudio
  ##
  # ?devtools::install
  ##
  require(MetaOrdDTA)
  ##
  # # Remove any lock files
  # unlink(file.path(local_pkg_dir, "00LOCK-MetaOrdDTA"), recursive = TRUE)
  #
  # # Clean compiled objects
  # pkgbuild::clean_dll(local_pkg_dir)
  #
  # # Then just load
  # devtools::load_all(local_pkg_dir)
}









