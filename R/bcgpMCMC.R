#' Draw samples from a bcgp model
#'
#' \code{bcgpMCMC} draws samples from the Bayesian Composite Gaussian Process model
#'
#' This draws samples from the posterior distribution for the Bayesian
#' Composite Gaussian Process (BCGP) model.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param y A vector containing the observed response values in the training
#' set.
#' @param prior A list containing the values for the prior parameters.
#' @param numUpdates The number of updates in the proposal stepsize adaptation phase.
#' @param numAdapt The number of samples within each update in the proposal stepsize
#' adaptation phase.
#' @param burnin The number of burnin samples to discard after the stepsize
#' adaptation phase is finished
#' @param nmcmc The number of samples to be kept for each Markov chain.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#' @param cores The number of cores to use when executing the Markov chains in
#' parallel. The default is to use the value of the \code{mc.cores} option if it
#' has been set and otherwise to default to 1 core.
#' @return An object of S4 class \code{bcgp} representing the fitted results.
#' @family Major functions
#' @examples
#'
#' x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' bcgp(x, y, priors)

bcgpMCMC  <- function(x, y, priors, inits, numUpdates, numAdapt,
                      burnin, nmcmc, chains = 1, cores = 1){

  nTrain <- nrow(x)
  iterations <- numUpdates*numAdapt + burnin + nmcmc
  epsV <- 1e-10
  tau2 <- 0.08

  G <- getCorMat(x, inits[[1]]$rhoG)
  L <- getCorMat(x, inits[[1]]$rhoL)
  R <- combineCorMats(inits[[1]]$w, G, L)
  C <- getCovMat(inits[[1]]$V, R, inits[[1]]$sig2eps)
  K <- inits[[1]]$sig2V * getCorMat(x, inits[[1]]$rhoV) + diag(epsV, nTrain)

  bfit <- list(G = G, L = L, R = R, C = C, K = K)
  # GC = getGPredC(train.xt,rhoGC);
  # LC = getGPredC(train.xt,rhoLC);
  # RC = getRPred(wC,GC,LC);
  # CC = getCPred(VC,RC,sig2epsC);
  # KC = sig2KC * getGPredC(train.xt,rhoVC) + epsV*eye(size(CC,1));

  # bfit <- new("bcgpfit",
  # )
  return(bfit)
}
