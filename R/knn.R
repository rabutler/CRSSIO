
#' Returns "k" and weights for the KNN algorithm
#' 
#' @param n The number of observations
#' 
#' @noRd
knn_params <- function(n)
{
  k <- sqrt(n) 
  
  # defines matrix for weights
  w <- 1 / matrix(seq_len(k), ncol = 1)
  w <- w / sum(w)
  
  list(k = k, weights = w)
}