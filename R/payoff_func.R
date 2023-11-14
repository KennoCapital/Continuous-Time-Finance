# Author: Nicklas Kenno Hansen

payoff_european <- function(S, K, type='C'){
  type <- toupper(type)
  
  if (type == 'C'){
    return(pmax(S-K, 0))
  }
  
  if (type == 'P'){
    return(pmax(K-S, 0))
  }
  else {
    stop("The passed option-type `", option_type, "` has not been implemented")
  }
}

