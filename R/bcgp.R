#' Draw samples from a bcgp model
#'
#' \code{bcgp} draws samples from the Bayesian Composite Gaussian Process model
#'
#' This draws samples from the posterior distribution for the Bayesian
#' Composite Gaussian Process (BCGP) model.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param y A vector containing the observed response values in the training
#' set.
#' @param priors Can be either the string "default" or a list containing the values
#' for the prior parameters.
#'
#' \describe{
#' \code{priors = "default"} (default):
#' The priors are given default values.
#'
#' \code{priors} via list:
#' Set prior values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createPriors()} will assist in the correct creation of this list.
#' }
#'
#' @param inits Can be either the string "random" or a list of length \code{chains}.
#' The elements of this list should be named lists, where each of these has the
#' name of a parameter and is used to specify the initial values for that parameter
#' for the corresponding chain.
#'
#' \describe{
#' \code{inits = "random"} (default):
#' The initial values will be generated randomly from their respective prior
#' distributions.
#'
#' \code{inits} via list:
#' Set initial values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createInits()} will assist in the correct creation of this list.
#' }
#'
#' @param noise If the data is assumed to be noise-free, then
#' \code{noise} should be \code{FALSE}. Otherwise, it should be
#' \code{TRUE}.
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
#' @return An object of S4 class \code{bcgpfit} representing the fitted results.
#' @family Major functions
#' @seealso \code{\link{createPriors}} \code{\link{createInits}}
#' @examples
#'
#' x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' bcgp(x, y, priors)
#' @export

bcgp  <- function(x, y, priors = "default",
                  inits = "random", noise = FALSE,
                  numUpdates = 5,
                  numAdapt = 1000,
                  burnin = 1000,
                  nmcmc = 10000,
                  chains = 4,
                  cores = getOption("mc.cores", 1L)){

  if(is.character(priors) && priors == "default"){
    priorList <- createPriors(x, noise = noise)
  }else if(is.list(priors)){
    ## TODO: Check to make sure the prior list is in the correct form
    priorList <- priors
  }else{
    stop("Incorrect specification of prior parameter values. Either use
         'default' or try calling createPriors() for correct specification
         before inputting your list.")
  }

  if(is.character(inits) && inits == "random"){
    initList <- createInits(x, priors = priorList, chains = chains)
  }else if(is.list(inits)){
    ## TODO: Check to make sure the inits list is in the correct form
    initList <- inits
  }else{
    stop("Incorrect specification of initial parameter values. Either use
         'random' or try calling createInits() for correct specification
         before inputting your list.")
  }

  yScaled <- scale(y, center = TRUE, scale = TRUE)
  xScaled <- scaleX(x)


  bfit <- bcgpMCMC(x = xScaled, y = yScaled, priors = priorList, inits = initList,
                   numUpdates, numAdapt,
                   burnin, nmcmc, chains = chains, cores = cores)

  slot(bfit, "scale") <- xScaled
  return(bfit)

}
