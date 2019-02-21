#' Create a covariance matrix.
#'
#' \code{getCovMat} returns a covariance matrix (with nugget).
#'
#' This creates a covariance matrix, \emph{C}, where
#' \deqn{C =  V^{0.5}RV^{0.5} + \sigma^2 I}
#' where V is a matrix with variances on the diagonal. In the \code{bcgp}
#' setting, \emph{R} is a correlation matrix resulting from
#' \code{\link{combineCorMats}}, \emph{V} is a vector of process variances,
#' and \emph{sig2} is the variance of the noise (or nugget).
#'
#' @param V A positive vector of length \emph{n}.
#' @param R An \emph{n x n} correlation matrix.
#' @param sig2 A positive scalar representing the variance of the noise
#' (or nugget).
#' @return An \emph{n x n} covariance matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' R <- getCorMat(x, rho)
#' sig2 <- 0.01
#' V <- rlnorm(n, -0.1, 0.1)
#' getCovMat(V, R, sig2)
#' @export

getCovMat <- function(V, R, sig2){

  # If I choose not to export this function, then I'll skip the error-checking
  # since the only time this function would be called is if everything is correct.
  stopifnot(all(V > 0),
            checkValidCorMat(R),
            (sig2 >= 0))

  rootV <- sqrt(V)
  C <- t(rootV*R) * rootV + diag(sig2, length(V))
  return(C)

  # NOTE: Surprisingly, this method is substantially faster, even for relatively large
  # matrices (tested up to 5000 x 5000), than doing sparse matrix multiplication in the
  # Matrix package.

  # Sparse matrix multiplication was orders of magnitude slower for small matrices
  # than the method immplemented above or for diag(V)^0.5 %*% R %*% diag(V)^0.5, which
  # was roughly the same speed as above for small matrices, but much slower for large
  # matrices.

}
