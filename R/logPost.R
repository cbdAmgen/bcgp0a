#' Calculate the log posterior density.
#'
#' \code{logPost} returns the log of the posterior density.
#'
#' This calculates the log of the posterior density in the \code{bcgp}
#' setting.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param y A vector containing the observed response values in the training
#' set.
#' @param params A named vector containing the current value of each of
#' the parameters.
#' @param priors A list containing the prior information.
#' @param C An \emph{n x n} covariance matrix for the training data.
#' @param K An \emph{n x n} covariance matrix for the variance process
#' at the training data locations.
#' @return A scalar
#' @examples

#' @export
logPost <- function(x, y, params, priors, C, K){

  return(logPost)
}

