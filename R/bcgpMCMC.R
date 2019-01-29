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
                      burnin, nmcmc, chains = 4, cores = 1){

  nTrain <- nrow(x)
  d <- ncol(x)
  iterations <- numUpdates*numAdapt + burnin + nmcmc
  epsV <- 1e-10
  tau2 <- 0.08

  rhoNames <- rep("rho", d)
  rhoGNames <- paste0(rhoNames, paste0("G", 1:d))
  rhoLNames <- paste0(rhoNames, paste0("L", 1:d))
  rhoVNames <- paste0(rhoNames, paste0("V", 1:d))
  rm(rhoNames)

  samples <- vector("list", chains)
  warmup <- vector("list", chains)
  acceptances <- vector("list", chains)

  ## TODO: Currently not parallelized. Need to make that happen
  for(i in 1:chains){

    ## TODO: right now I'm keeping all draws in a single matrix with column names
    ## Should I have each parameter group in its own separate matrix? Either way,
    ## at the end, they'll be put into a named list

    ## Side note, working with matrices will be faster than working with
    ## data frames or lists

    allDraws <- matrix(NA, nrow = iterations, ncol = 5 + 3*d + nTrain)
    row1 <- unlist(inits[[i]])
    colnames(allDraws) <- names(row1)
    allDraws[1, ] <- row1

    allAcceptances <- matrix(0, nrow = iterations, ncol = 5 + 3*d + nTrain)
    colnames(allAcceptances) <- names(row1)
    allAcceptances[1, ] <- 1

    rm(row1)

    G <- getCorMat(x, inits[[i]]$rhoG)
    L <- getCorMat(x, inits[[i]]$rhoL)
    R <- combineCorMats(inits[[i]]$w, G, L)
    C <- getCovMat(inits[[i]]$V, R, inits[[i]]$sig2eps)
    K <- inits[[i]]$sig2V * getCorMat(x, inits[[i]]$rhoV) + diag(epsV, nTrain)

    propWidths <- c(0.4, rep(0.07, d), rep(0.03, d), 1e-3, rep(0.25, d))
    names(propWidths) <- c("w", rhoGNames, rhoLNames, "sig2eps", rhoVNames)

    ## TODO: This is where the work goes
    for(j in 2:iterations){
      allDraws[j, ] <- runif(ncol(allDraws), 0, 1)
      allAcceptances[j, ] <- sample.int(2, size = ncol(allDraws), replace = TRUE) - 1
    }


    warmup[[i]] <- as.list(data.frame(allDraws[1:(numUpdates*numAdapt + burnin), ]))
    samples[[i]] <- as.list(data.frame(allDraws[(iterations - nmcmc + 1):iterations, ]))
    acceptances[[i]] <- as.list(
      data.frame(allAcceptances[(iterations - nmcmc + 1):iterations, ]))

  }


  sim <- list(samples = samples,
              warmup = warmup,
              acceptances = acceptances,
              chains = chains,
              numUpdates = numUpdates,
              numAdapt = numAdapt,
              burnin = burnin,
              nmcmc = nmcmc)

  bfit <- new("bcgpfit",
              model_pars = c("beta0", "w", "rhoG", "rhoL", "sig2eps", "sig2Y",
                             "muV", "rhoV", "sig2V"),
              par_dims = list(beta0 = 1,
                              w = 1,
                              rhoG = d,
                              rhoL = d,
                              sig2eps = 1,
                              sig2Y = nTrain,
                              muV = 1,
                              rhoV = d,
                              sig2V = 1),
              sim = sim,
              priors = priors,
              inits = inits,
              args = list(chains = chains,
                          numUpdates = numUpdates,
                          numAdapt = numAdapt,
                          burnin = burnin,
                          nmcmc = nmcmc),
              algorithm = "M-H and Gibbs")

  return(bfit)
}
