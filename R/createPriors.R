#' Create a list with prior information.
#'
#' \code{createPriors} returns a list that contains default values for the priors.
#'
#' This creates a list that contains default values for the priors. The user
#' can change the values as they like prior to inputting the list into
#' \code{bcgp}.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param noise If the data is assumed to be noise-free, then
#' \code{noise} should be \code{FALSE}. Otherwise, it should be
#' \code{TRUE}.
#' @return A list containing the default values for all the prior parameters.
#' @family preprocessing functions
#' @seealso \code{\link{bcgp}}
#' @section TODO: Decide whether to add options for "heteroscedastic" and "composite"
#' @examples
#' x <- matrix(runif(40), ncol= 4, nrow = 10)
#' createPriors(x)
#' createPriors(x, noise = TRUE)
#' @export

createPriors  <- function(x, noise = FALSE){
  d <- ncol(x)
  priorList <- list(w = list(lower = 0.5,
                             upper = 1.0,
                             alpha = 1,
                             beta = 1),
                    rhoG = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    rhoL = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    muV = list(betaV = -0.1,
                               sig2 = 0.1),
                    rhoV = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    sig2V = list(alpha = 2 + sqrt(0.1),
                                 beta = 100/(1+sqrt(1/10))))

  if(!noise){
    priorList$sig2eps <- list(alpha = 1e-3,
                              beta = 1e-3)
  }else{
    priorList$sig2eps <- list(alpha = 1e-1,
                              beta = 1e-1)
  }

  return(priorList)
}
