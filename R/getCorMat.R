#' Create a correlation matrix.
#'
#' \code{getCorMat} returns a correlation matrix.
#'
#' This creates a correlation matrix, \emph{R}, where \deqn{R_{ij} = \prod_{k =
#' 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}{R_ij = \prod_{k =
#' 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}
#'
#' @param x An \emph{n x d} matrix.
#' @param rho A vector of length d.
#' @return An \emph{n x n} correlation matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' getCorMat(x, rho)

getCorMat <- function(x, rho){
  n <- nrow(x)
  d <- ncol(x)

  R <- matrix(0, nrow = n, ncol = n)

  for(i in 1:(n-1)){
    R[i, i] <- 1.0
    for(j in (i+1):n){
      R[i, j] <- 1.0
      dist <- x[i, ] - x[j, ]
      R[i, j] <- prod(rho^(16*dist^2))
      R[j, i] <- R[i, j]
    }

  }
  R[n, n] <- 1.0
  return(R)
}
