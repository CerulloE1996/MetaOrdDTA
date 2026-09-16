



#' softplus
#' @export
softplus <- function(x) { 
  
    return(log(1.0 + exp(x)))
  
}


#' softplus_trans
#' @export
softplus_scaled <- function(x) { 
  
  log_2_recip <- 1.4426950408889634
  return(log_2_recip*softplus(x))
  
}




#' sp_scaled
#' @export
sp_scaled <- softplus_scaled





