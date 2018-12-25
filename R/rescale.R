#' Rescale a vector to [0, 1]
#' @description Rescales a vector to [0, 1]
#' @details Rescales a vector to [0, 1].
#' @usage rescale(x)
#' @param \code{x} A vector input
#' @return A vector rescaled to [0, 1]
#' @examples
#' rescale(rnorm(10, 0, 100))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
rescale <- function(x){

  minX <- min(x)
  rangeX <- max(x) - minX
  scaled <- (x - minX)/rangeX

  return(scaled)

}
