initFunc <- function(initList, priors, x){

  d <- ncol(x)
  n <- nrow(x)
  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$w <- priors$w$lower + rbeta(1, priors$w$alpha, priors$w$beta)*
    (priors$w$upper - priors$w$lower)
  initReturn$rhoG <- rbeta(d, priors$rhoG$alpha, priors$rhoG$beta)
  initReturn$rhoL <- initReturn$rhoG * rbeta(d, priors$rhoL$alpha, priors$rhoL$beta)
  initReturn$sig2eps <- max(2* .Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$muV <- rnorm(1, priors$muV$betaV, sqrt(priors$muV$sig2))
  initReturn$rhoV <- rbeta(d, priors$rhoV$alpha, priors$rhoV$beta)
  initReturn$sig2V <- 1/rgamma(1, priors$sig2V$alpha, scale = priors$sig2V$beta)
  K <- initReturn$sig2V * getCorMat(x,initReturn$rhoV) + 1e-10*diag(n)
  initReturn$V <- exp(MASS::mvrnorm(1, initReturn$muV*rep(1, n), K))
  return(initReturn)

}

checkValidCorMat <- function(x){

  if(!is.matrix(x)) return(FALSE) # check that it's a matrix
  if(!is.numeric(x)) return(FALSE) # check that it's numeric
  if(!(nrow(x) == ncol(x))) return(FALSE) # check square
  if(isFALSE(all(diag(x) == 1))) return(FALSE) # check 1's on diagonal
  if(!(sum(x == t(x)) == (nrow(x)^2))) return(FALSE) # check symmetric

  # check positive semi-definite
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < 1e-8] <- 0
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)

}

#' Create a correlation matrix.
#'
#' \code{x1Ainvx2} returns a list where each element is the result of multiplying
#' a vector times a matrix inverse times another vector.
#'
#' This returns a list where each element is the result of \eqn{x_1^{\top}A^{-1}x_2}
#'
#' @param x1 A list of vectors of length \emph{n}. Generally, in the bcgp setting,
#' \emph{n} will be the number of training data locations.
#' @param A An \emph{n x n} matrix. In the \code{bcgp} setting, \emph{A} is a
#' covariance matrix.
#' @return x2 A list of vectors of length \emph{n}. Generally, in the bcgp setting,
#' \emph{n} will be the number of training data locations.
#' @examples
#' n <- 10
#' x <- matrix(sort(runif(n)), ncol = 1)
#' A <- getCorMat(x, rho = 0.6)
#' x1 <- list(rep(1, n), rep(1,n))
#' x2 <- list(rep(1, n), MASS::mvrnorm(1, rep(0, n), A))
#' x1Ainvx2(x1, A, x2)
#' @export
x1Ainvx2 <- function(x1, A, x2){

  cholAR <- try(chol(A), silent = TRUE)
  if(is.matrix(cholAR)){

    tmp1 <- lapply(x1, forwardsolve, l = t(cholAR), k = ncol(cholAR))
    tmp2 <- lapply(x2, forwardsolve, l = t(cholAR), k = ncol(cholAR))

    xAinvy <- colSums(mapply(`*`, tmp1, tmp2))

  }else{

    Ainvy <- try(lapply(x2, solve, a = A), silent = TRUE)
    if(all(sapply(Ainvy, is.numeric))){
      xAinvy <- colSums(mapply(`*`, x1, Ainvy))
    }else{
      svdA <- svd(A)
      middle <- svdA$v %*% diag(1/svdA$d) %*% t(svdA$u)
      end <- lapply(x2, `%*%`, middle)
      xAinvy <- colSums(mapply(`*`, x1, end))
    }
  }
  return(xAinvy)
}


