
 







#' plot_param_group_batched
#' @keywords internal
#' @export
plot_param_group_batched <- function(  draws_array, 
                                       param_string,
                                       condition = "exact_match",
                                       plot_type = c("density", "trace"),
                                       batch_size = 9,
                                       debugging = FALSE
) {

    
         ## Make sure input array is in correct format:
         try({  
           draws_array <- MetaOrdDTA:::format_named_array_for_bayesplot(draws_array)
         })
  
         plot_type <- match.arg(plot_type)
        
         all_param_names_vec <- dimnames(draws_array)$variable
         if (debugging) print(paste(head(all_param_names_vec, 10)))
        
         str(draws_array)
   
         ## Get param names:
         valid_params_list <-  filter_param_names_string_batched(   all_param_names_vec = all_param_names_vec,
                                                                    param_string = param_string,
                                                                    condition = condition,
                                                                    debugging = debugging)

         # Choose plot function:
         plot_func <- switch(plot_type,
                             density = bayesplot::mcmc_dens,
                             trace =   bayesplot::mcmc_trace)

         plots <- list()
         for (l in 1:length(valid_params_list))  {
           
               valid_params_vec <- valid_params_list[[l]]
               ##
               n_batches <- ceiling(length(valid_params_vec) / batch_size) ## Plot in batches:
               
               plots[[l]] <- list()
               for (i in 1:n_batches) {
      
                         start_idx <- (i - 1) * batch_size + 1
                         end_idx <- min(i * batch_size, length(valid_params_vec))
                         batch_params <- valid_params_vec[start_idx:end_idx]
      
                         if (debugging) print(cat("batch_params = ", batch_params))
      
                         ### plot:
                         ##
                         if (debugging) str(draws_array)
                         ##
                         plots[[l]][[i]] <- plot_func(draws_array[,,batch_params, drop=FALSE]) +
                                       ggplot2::ggtitle(paste0(valid_params_vec,
                                                               " (batch ", i, " of ", n_batches, ")"))
      
               }
         }
         
         #str(plots)
         
         #plots[[2]][[1]]
         # 
         return(plots)
      
}


 



# named_draws_array = stan_draws_array
# param_strings_vec = grouped_params$names_main
# condition = "exact_match"
# exclude_NA = exclude_NA
# debugging = FALSE
# 
# 
# 
# draws_array_list = param_draws_filtered_list
# 
# 
# named_draws_array = stan_draws_array
# 
# 
# 
# 
# 
# 
# param_string <- "beta_mu"
# all_param_names_vec <-  dimnames(named_draws_array)[[3]]
# length(all_param_names_vec)
# 
# 
# filtered_param_names_string(all_param_names_vec = all_param_names_vec,
#                             param_string = "beta",
#                             condition = "exact_match")
# 
# 
# filter_param_names_string_batched(all_param_names_vec = all_param_names_vec,
#                             param_string = c("beta", "raw_scale"),
#                             condition = "exact_match")
# 
# 
# 
# find_all_stan_param_names(named_draws_array)
