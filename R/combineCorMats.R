#' Combine correlation matrices.
#'
#' \code{combineCorMats} returns a correlation matrix that is a weighted sum of two
#' other correlation matrices.
#'
#' This creates a correlation matrix, \emph{R}, where \deqn{R = wG + (1-w)L}In
#' the \code{bcgp} setting, \emph{G} is the global correlation matrix, \emph{L}
#' is the local correlation matrix, and \emph{w} is the weight.
#'
#' @param w A scalar
#' @param G An \emph{n x n} correlation matrix, often the result of \code{\link{getCorMat}}
#' @param L An \emph{n x n} correlation matrix, often the result of \code{\link{getCorMat}}
#' @return An \emph{n x n} correlation matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rhoG <- runif(d, 0, 1)
#' rhoL <- runif(d, 0, rhoG)
#' G <- getCorMat(x, rhoG)
#' L <- getCorMat(x, rhoL)
#' w <- runif(1, 0.5, 1)
#' combineCorMats(w, G, L)

combineCorMats <- function(w, G, L){
  R <- w*G + (1 - w)*L
  return(R)
}

# getCovMat <- function(V, R, sig2eps){
#
#   V <- diag(V)
#   C <- V^0.5 %*% R %*% V^0.5 + sig2eps*diag(1, nrow(V))
#   return(C)
#
# }
