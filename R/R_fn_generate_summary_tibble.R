

#' generate_summary_tibble
#' @export
generate_summary_tibble <- function(n_threads,
                                    input_array,
                                    iter_dim,
                                    chain_dim,
                                    param_names,
                                    compute_nested_rhat,
                                    n_superchains) {

              ##
              ## ---- NicoStan's summary of the trace (its ESS, R-hat and nR-hat) ---------------------------------------------------------------------
              ##
              params_dim <-  factorial(3) - (iter_dim + chain_dim)
              reformatted_trace <-  make_array_dims_equal_Bayesplot_order(  input_array = input_array,
                                                                            iter_dim = iter_dim,
                                                                            chain_dim = chain_dim,
                                                                            params_dim = params_dim)
              n_chains <- dim(reformatted_trace)[2]
              n_superchains <- if_null_then_set_to(n_superchains, n_chains)
              ## (NicoStan reads the trace as parameters x iterations x chains)
              summary_tibble <- NicoStan::generate_summary_tibble( n_threads = n_threads,
                                                                   trace = aperm(reformatted_trace, c(3, 1, 2)),
                                                                   param_names = param_names,
                                                                   n_to_compute = length(param_names),
                                                                   compute_nested_rhat = compute_nested_rhat,
                                                                   n_chains = n_chains,
                                                                   n_superchains = n_superchains)
              ##
              return(summary_tibble)

}



#' fn_draws_data_frame_of_an_array
#' @keywords internal
fn_draws_data_frame_of_an_array <- function(trace_array) {

              ##
              ## ---- One row per draw (iterations x chains x variables -> data frame) ----------------------------------------------------------------
              ##
              n_iter   <- dim(trace_array)[1]
              n_chains <- dim(trace_array)[2]
              draws_matrix <- matrix(trace_array, nrow = n_iter * n_chains)
              colnames(draws_matrix) <- dimnames(trace_array)[[3]]
              draws_data_frame <- data.frame(draws_matrix, check.names = FALSE)
              draws_data_frame$.chain <- rep(seq_len(n_chains), each = n_iter)
              draws_data_frame$.iteration <- rep(seq_len(n_iter), times = n_chains)
              draws_data_frame$.draw <- seq_len(n_iter * n_chains)
              ##
              return(draws_data_frame)

}

