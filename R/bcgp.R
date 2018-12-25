#' Draw samples from a bcgp model
#'
#' \code{bcgp} draws samples from the Bayesian Composite Gaussian Process model
#'
#' This draws samples from the posterior distribution for the Bayesian
#' Composite Gaussian Process (BCGP) model.
#'
#' @param xTrain An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param yTrain A vector containing the observed response values in the training
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
#' @param numUpdates The number of updates in the proposal stepsize adaptation phase.
#' @param numAdapt The number of samples within each update in the proposal stepsize
#' adaptation phase.
#' @param burnin The number of burnin samples to discard after the stepsize
#' adaptation phase is finished
#' @param nmcmc The number of samples to be kept for each Markov chain.
#' @param chains The number of Markov chains. The default is 4.
#' @return An object of S4 class \code{bcgpfit} representing the fitted results.
#' @family Major functions
#' @seealso \code{\link{createPriors}} \code{\link{createInits}}
#' @examples
#'
#' xTrain <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' yTrain <- xTrain[, 1] + sin(xTrain[, 2])
#' priors <- createPriors(noise = FALSE, d = 2)
#' bcgp(xTrain, yTrain, priors)
#' @export

bcgp  <- function(xTrain, yTrain, priors = "default",
                  inits = "random",
                  numUpdates = 5,
                  numAdapt = 1000,
                  burnin = 1000,
                  nmcmc = 10000,
                  chains = 4,
                  cores = getOption("mc.cores", 1L),
                  noise = FALSE){

  if(priors == "default"){
    priorList <- createPriors(xTrain, noise = noise)
  }else if(is.list(priors)){
    ## FIX: Check to make sure the prior list is in the correct form
    priorList <- priors
  }else{
    stop("Incorrect specification of prior parameter values. Either use
         'default' or try calling createPriors() for correct specification.")
  }

  if(init == "random"){
    initList <- createInits(xTrain, priors = priorList, chains = chains)
  }else if(is.list(priors)){
    ## FIX: Check to make sure the init list is in the correct form
    initList <- init
  }else{
    stop("Incorrect specification of initial parameter values. Either use
         'random' or try calling createInits() for correct specification.")
  }

  bcgpMCMC(xTrain = xTrain, yTrain = yTrain, priors = priorList, inits = initList)

  # bfit <- new("bcgpfit",
  #             )
  return(bfit)
}
