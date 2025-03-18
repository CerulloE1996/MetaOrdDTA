
#' transform_stan_path
#' @keywords internal
#' @export
transform_stan_path <- function(stan_path) {
  
  # Remove any leading ~/ if present
  path_no_tilde <- sub("^~/", "", stan_path)
  
  # Extract directory and filename
  dir_path <- dirname(path_no_tilde)
  filename <- basename(path_no_tilde)
  
  # Remove .stan extension and add _model.so or _model.dll
  if (.Platform$OS.type == "windows") {
    
    model_name <- sub("\\.stan$", "_model.dll", filename)
    final_path <- file.path(dir_path, model_name)
    
  } else { 
    
    model_name <- sub("\\.stan$", "_model.so", filename)
    final_path <- file.path("/", dir_path, model_name)
    
  }
  
  # final_path <- file.path("~", dir_path, model_name)
  
  return(final_path)
  
}

