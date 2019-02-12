#' Calculate the log of the determinant of a matrix.
#'
#' \code{logDet} returns the log of the determinant of a matrix..
#'
#' This calculates the log of the determinant of a matrix. In the \code{bcgp}
#' setting, \emph{C} is a covariance matrix resulting from
#' \code{\link{getCovMat}}.
#'
#' @param C An \emph{n x n} covariance matrix.
#' @return A scalar
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' R <- getCorMat(x, rho)
#' sig2 <- 0.01
#' V <- rlnorm(n, -0.1, 0.1)
#' C <- getCovMat(V, R, sig2)
#' logDet(C)
#' @export
logDet <- function(C){
  a <- try(chol(C), silent = TRUE)
  if(is.matrix(a)){
    logDeterminant <- 2 * sum(log(diag(a)))
  }else{
    a <- determinant(C, logarithm = TRUE)
    logDeterminant <- as.numeric(a$modulus) # since every C matrix will be pd
  }
  return(logDeterminant)
}

