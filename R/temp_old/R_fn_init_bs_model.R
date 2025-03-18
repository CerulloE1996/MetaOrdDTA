
  
  
#' get_BayesMVP_Stan_paths
#' @keywords internal
#' @export
get_BayesMVP_Stan_paths <- function() {
  
      Sys.setenv(STAN_THREADS = "true")
      
      # Get package directory paths
      pkg_dir <- system.file(package = "BayesMVP")
      data_dir <- file.path(pkg_dir, "stan_data")  # directory to store data inc. JSON data files
      stan_dir <- file.path(pkg_dir, "stan_models")
      ##
      print(paste("pkg_dir = ", pkg_dir))
      print(paste("data_dir = ", data_dir))
      print(paste("stan_dir = ", stan_dir))
      
      return(list(pkg_dir = pkg_dir,
                  data_dir = data_dir,
                  stan_dir = stan_dir))
  
}





#' init_bs_model_external
#' @keywords internal
#' @export
init_bs_model_external <- function(Stan_data_list, 
                                   Stan_model_file_path) {
  
        # Get package directory paths
        outs <- get_BayesMVP_Stan_paths()
        pkg_dir <- outs$pkg_dir
        data_dir <- outs$data_dir
        stan_dir <- outs$stan_dir
        
        ##
        # # Create data directory if it doesn't exist
        # if (!dir.exists(data_dir)) {
        #   dir.create(data_dir, recursive = TRUE)
        # }
        # 
        # ## make persistent (non-temp) JSON data file path with unique identifier:
        # data_hash <- digest::digest(Stan_data_list)  # Hash the data to create unique identifier
        # json_filename <- paste0("data_", data_hash, ".json")
        # json_file_path <- file.path(data_dir, json_filename)
        # ##
        # ## write JSON data using cmdstanr:
        # cmdstanr::write_stan_json(Stan_data_list, json_file_path)
        # ##
        # json_file_path <- normalizePath(json_file_path)
        ##
        json_file_path <- convert_Stan_data_list_to_JSON(Stan_data_list, data_dir = data_dir)
        ##
        validated_json_string <- paste(readLines(json_file_path), collapse="")
        jsonlite::write_json(x = validated_json_string, path = json_file_path)
        ##
        ## Create bridgestan model:
        Stan_model <- file.path(Stan_model_file_path)
        ##
        bs_model <- bridgestan::StanModel$new(
          lib = Stan_model,
          data = validated_json_string,
          seed = 123)
        ##
        json_file_path <- normalizePath(json_file_path)
        Stan_model_file_path <- normalizePath(Stan_model_file_path)
        
        # Return both model and data path
        return(list(
          bs_model = bs_model,
          json_file_path = json_file_path,
          Stan_model_file_path = Stan_model_file_path))
  
}







#' init_bs_model_internal
#' @keywords internal
#' @export
init_bs_model_internal <- function( Stan_data_list, 
                                    Stan_model_name) {
  
        # Get package directory paths
        outs <- get_BayesMVP_Stan_paths()
        pkg_dir <- outs$pkg_dir
        data_dir <- outs$data_dir
        stan_dir <- outs$stan_dir
        
        ## Stan model path
        Stan_model_file_path <- file.path(stan_dir, Stan_model_name)
        
        outs_bs_model <- init_bs_model_external(Stan_data_list = Stan_data_list, 
                                                Stan_model_file_path = Stan_model_file_path)
        
        return(outs_bs_model)
  
}














