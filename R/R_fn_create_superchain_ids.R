
## Create superchain IDs for chain grouping when many chains are used


#' #' create_superchain_ids
#' @keywords internal
#' @export
create_superchain_ids <- function(n_superchains,
                                  n_chains) {
  
      # rep(1:n_superchains, each = n_chains/n_superchains)
      # or more safely to handle non-perfect division:
      chains_per_superchain <- ceiling(n_chains/n_superchains)
      superchain_ids <- rep(1:n_superchains, each = chains_per_superchain)[1:n_chains]
      
      return(superchain_ids)
  
}












