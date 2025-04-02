



#' default_cts
#' @keywords internal
#' @export
default_cts <- function() { 
  
  cts <- FALSE
  ##
  # print(paste("default_cts = ", default_cts))
  return(list(cts))
  
}




#' default_network
#' @keywords internal
#' @export
default_network <- function(x) { 
  
  
      network <- NA
      ##
      try({ 
         is_x_list <-  MetaOrdDTA:::try_silent(is.list(x))
      })
      try({ 
         is_x_nested_list <- MetaOrdDTA:::try_silent(is.list(x[[1]]))
      })
      ##
      if (is_x_list) { 
        network <- FALSE
      } else if (is_x_nested_list) { 
        network <- TRUE
      }
      
      # print(paste("network = ", network))
      return(list(network))
  
  
}




#' default_parameterisation
#' @keywords internal
#' @export
default_model_parameterisation <- function(cts) { 
  
  if (cts) { 
    model_parameterisation <- "Jones"
  } else   {
    model_parameterisation <- "HSROC" ## ord "R&G"
  }
  
  # print(paste("model_parameterisation = ", model_parameterisation))
  return(list(model_parameterisation))
  
}





#' default_random_thresholds
#' @keywords internal
#' @export
default_random_thresholds <- function(cts) { 
  
  if (cts) { 
    random_thresholds <- FALSE
  } else   {
    random_thresholds <- FALSE ## default is fixed between-study threshold model)
  }
  
  # print(paste("random_thresholds = ", random_thresholds))
  return(list(random_thresholds))
  
}




#' default_Dirichlet_random_effects_type
#' @keywords internal
#' @export
default_Dirichlet_random_effects_type <- function(cts, 
                                                  random_thresholds) { 
  
  if (cts) { 
    Dirichlet_random_effects_type <- "none"
  } else   {
      if (random_thresholds == TRUE) { 
        Dirichlet_random_effects_type <- "alpha"
      } else { 
        Dirichlet_random_effects_type <- "none" ## only relvent if using random-effect (between-study) thresholds. 
      }
  }
  
  # print(paste("Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))
  return(list(Dirichlet_random_effects_type))
  
}




#' default_box_cox
#' @keywords internal
#' @export
default_box_cox <- function(cts) { 
  
    if (cts) { 
      box_cox <- TRUE
    } else   {
      box_cox <- FALSE ## only relvant for continuous models.
    }
    
    # print(paste("box_cox = ", box_cox))
    return(list(box_cox))
    
}




#' default_softplus
#' @keywords internal
#' @export
default_softplus <- function() { 

    softplus <- TRUE ## default is to use softplus on ALL models (instead of exp() / log-normal)
    ##
    # print(paste("softplus = ", softplus))
    return(list(softplus))
  
}














