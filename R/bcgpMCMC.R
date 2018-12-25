#' Draw samples from a bcgp model
#'
#' \code{bcgpMCMC} draws samples from the Bayesian Composite Gaussian Process model
#'
#' This draws samples from the posterior distribution for the Bayesian
#' Composite Gaussian Process (BCGP) model.
#'
#' @param X An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param yTrain Avector containing the observed response values in the training
#' set.
#' @param prior A list containing the values for the prior parameters.
#' @param numUpdates The number of updates in the proposal stepsize adaptation phase.
#' @param numAdapt The number of samples within each update in the proposal stepsize
#' adaptation phase.
#' @param burnin The number of burnin samples to discard after the stepsize
#' adaptation phase is finished
#' @return An object of S4 class \code{bcgp} representing the fitted results.
#' @family Major functions
#' @examples
#'
#' xTrain <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' yTrain <- xTrain[, 1] + sin(xTrain[, 2])
#' priors <- createPriors(noise = FALSE, d = 2)
#' bcgp(xTrain, yTrain, priors)

bcgpMCMC  <- function(xTrain, yTrain, priors, inits, numUpdates, numAdapt,
                      burnin, nmcmc){

  y <- scale(yTrain, center = TRUE, scale = TRUE)
  X <- apply(xTrain, 2, rescale)

  epsV <- 1e-10
  tau2 <- 0.08

  bfit <- new("bcgpfit",
  )
  return(bfit)
}
