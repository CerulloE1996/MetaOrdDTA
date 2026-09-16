


#' set_pkg_eg_path_and_wd
#' @keywords internal
#' @export
set_pkg_eg_path_and_wd <- function() {
        
        ## Find user_root_dir:
        os <- .Platform$OS.type
        
        if (os == "unix") { 
              user_root_dir <- Sys.getenv("PWD")
        } else if (os == "windows") { 
              user_root_dir <- Sys.getenv("USERPROFILE")
        }
        
        print(paste("user_root_dir = ", user_root_dir))
        ##
        user_MetaOrdDTA_dir <- file.path(user_root_dir, "MetaOrdDTA")
        print(paste("user_MetaOrdDTA_dir = ", user_MetaOrdDTA_dir))
        ##
        pkg_example_path <- file.path(user_MetaOrdDTA_dir, "examples")
        print(paste("pkg_example_path = ", pkg_example_path))
        ##
        ## Set working directory:
        ## Create directory if it doesn't exist
        if (!dir.exists(pkg_example_path)) {
          dir.create(pkg_example_path, recursive = TRUE)
          print(paste("Created directory:", pkg_example_path))
        }
        #
        try ({ 
            setwd(pkg_example_path)
            message(paste("Working directory set to: ", pkg_example_path))
        })
        ##
        ## Find user_pkg_install_dir:
        ##
        user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
        print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
        ##
        ## Find pkg_install_path:
        ##
        pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
        print(paste("pkg_install_path = ", pkg_install_path))
        
        outs <- list(user_root_dir = user_root_dir,
                     user_MetaOrdDTA_dir = user_MetaOrdDTA_dir,
                     pkg_example_path = pkg_example_path)
        
        return(outs)
  
}

