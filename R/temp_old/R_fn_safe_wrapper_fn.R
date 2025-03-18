

#' safe_test_wrapper_1
#' @keywords internal
#' @export
safe_test_wrapper_1 <- function(expr, 
                                ..., 
                                timeout = 5) {
          
          # Get the additional arguments
          args <- list(...)
          
          if (.Platform$OS.type == "unix") {
            library(parallel)
            result <- tryCatch({
              p <- mcparallel({
                # Create the environment with the passed arguments
                list2env(args, environment())
                expr
              })
              val <- mccollect(p, wait = FALSE, timeout = timeout)
              if (is.null(val)) {
                mckill(p, signal = 15)
                return(list(status = "timeout", error = "Function execution timed out"))
              }
              list(status = "success", result = val)
            }, error = function(e) {
              list(status = "error", error = conditionMessage(e))
            })
          } else {
            library(processx)
            
            # Fix 1: Use normalized paths
            temp_script <- normalizePath(tempfile(fileext = ".R"), winslash = "/", mustWork = FALSE)
            args_file <- normalizePath(tempfile(fileext = ".rds"), winslash = "/", mustWork = FALSE)
            
            # Save arguments to a temporary RDS
            saveRDS(args, args_file)
            
            # Fix 2: Escape the path properly for the R script
            args_file_escaped <- gsub("\\\\", "/", args_file)
            
            # Write the expression to a temporary file
            code_str <- deparse(substitute(expr))
            if (length(code_str) > 1) code_str <- paste(code_str, collapse = "\n")
            
            cat(sprintf('
                    tryCatch({
                        # Load the arguments
                        args <- readRDS("%s")
                        list2env(args, environment())
                        
                        result <- {
                            %s
                        }
                        saveRDS(list(status = "success", result = result), "temp_result.rds")
                    }, error = function(e) {
                        saveRDS(list(status = "error", error = conditionMessage(e)), "temp_result.rds")
                    })
                ', args_file_escaped, code_str), file = temp_script)
            
            # Run R in a separate process
            r_process <- processx::process$new(
              command = "R",
              args = c("--vanilla", "-f", temp_script),
              cleanup = TRUE,
              stdout = "|",
              stderr = "|"
            )
            
            # Wait with timeout
            timeout_occurred <- FALSE
            start_time <- Sys.time()
            while (r_process$is_alive()) {
              if (difftime(Sys.time(), start_time, units = "secs") > timeout) {
                r_process$kill()
                timeout_occurred <- TRUE
                break
              }
              Sys.sleep(0.1)
            }
            
            if (timeout_occurred) {
              result <- list(status = "timeout", error = "Function execution timed out")
            } else if (file.exists("temp_result.rds")) {
              result <- readRDS("temp_result.rds")
              file.remove("temp_result.rds")
            } else {
              result <- list(status = "error", 
                             error = paste("Process error.",
                                           "STDOUT:", r_process$read_all_output(),
                                           "STDERR:", r_process$read_all_error()))
            }
            
            # Cleanup
            unlink(temp_script)
            unlink(args_file)
          }
          
          return(result)
          
}







#' safe_test_wrapper_2
#' @keywords internal
#' @export
safe_test_wrapper_2 <- function(expr, 
                                ...,
                                timeout = 5) {
          
          args <- list(...)
          
          if (.Platform$OS.type == "unix") {
            
            library(parallel)
            result <- tryCatch({
              p <- mcparallel({
                list2env(args, environment())
                expr
              })
              val <- mccollect(p, wait = FALSE, timeout = timeout)
              if (is.null(val)) {
                mckill(p, signal = 15)
                return(list(status = "timeout", error = "Function execution timed out"))
              }
              list(status = "success", result = val)
            }, error = function(e) {
              list(status = "error", error = conditionMessage(e))
            })
            
          } else {
            
            library(processx)
            
            temp_script <- normalizePath(tempfile(fileext = ".R"), winslash = "/", mustWork = FALSE)
            args_file <- normalizePath(tempfile(fileext = ".rds"), winslash = "/", mustWork = FALSE)
            
            saveRDS(args, args_file)
            args_file_escaped <- gsub("\\\\", "/", args_file)
            
            code_str <- deparse(substitute(expr))
            if (length(code_str) > 1) code_str <- paste(code_str, collapse = "\n")
            
            cat(sprintf('
                    options(error = function() {
                        traceback(2)
                        saveRDS(
                            list(
                                status = "error",
                                error = geterrmessage(),
                                traceback = capture.output(traceback())
                            ),
                            "temp_result.rds"
                        )
                    })
                    
                    result <- tryCatch({
                       ## library(BayesMVP)
                        library(bridgestan)
                        library(cmdstanr)
                        
                        args <- readRDS("%s")
                        list2env(args, environment())
                        
                        output <- capture.output({
                            result <- {
                                %s
                            }
                        })
                        
                        saveRDS(
                            list(
                                status = "success",
                                result = result,
                                output = output
                            ),
                            "temp_result.rds"
                        )
                    }, error = function(e) {
                        saveRDS(
                            list(
                                status = "error",
                                error = conditionMessage(e),
                                call = deparse(e$call),
                                traceback = capture.output(traceback())
                            ),
                            "temp_result.rds"
                        )
                        stop(e)
                    })
                ', args_file_escaped, code_str), file = temp_script)
            
            r_process <- processx::process$new(
              command = "R",
              args = c("--vanilla", "-f", temp_script),
              cleanup = TRUE,
              stdout = "|",
              stderr = "|"
            )
            
            # Fixed: Initialize the timeout flag
            is_timeout <- FALSE
            start_time <- Sys.time()
            
            while (r_process$is_alive()) {
              if (difftime(Sys.time(), start_time, units = "secs") > timeout) {
                r_process$kill()
                is_timeout <- TRUE
                break
              }
              Sys.sleep(0.1)
            }
            
            if (is_timeout) {
              result <- list(status = "timeout", error = "Function execution timed out")
            } else if (file.exists("temp_result.rds")) {
              result <- readRDS("temp_result.rds")
              if (result$status == "error") {
                result$error <- paste0(
                  "Error: ", result$error, "\n",
                  "Traceback:\n", 
                  paste(result$traceback, collapse = "\n")
                )
              }
            } else {
              stdout <- r_process$read_all_output()
              stderr <- r_process$read_all_error()
              result <- list(
                status = "error",
                error = paste0(
                  "Process failed.\n",
                  "STDOUT:\n", stdout, "\n",
                  "STDERR:\n", stderr
                )
              )
            }
            
            # Cleanup
            unlink(temp_script)
            unlink(args_file)
          }
          
          return(result)
  
}

