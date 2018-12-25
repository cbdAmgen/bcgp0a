#' Create a list with initial values.
#'
#' \code{createInits} returns a list that contains randomly generated initial
#' values.
#'
#' This creates a list of length \code{chains} that contains randomly generated
#' initial values. The intention is to specify initial values for the parameters
#' for the Markov chains. The user can change the values as they like prior to
#' inputting the list into \code{bcgp}.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param priors A list that contains the parameter values for the priors.
#' @param chains The number of Markov chains. The default is 4.
#' @return A list of length \code{chains} The elements of this list will be named
#' lists, where each of these has the name of a parameter.
#' @family preprocessing functions
#' @seealso \code{\link{bcgp}}
#' @section TODO: Decide whether to add options for "heteroscedastic" and "composite"
#' @examples
#' x <- matrix(runif(40), ncol= 4, nrow = 10)
#' createInits(x)
#' createInits(x, priors = createPriors(x, noise = TRUE), chains = 2)
#' @export

createInits  <- function(x, priors = createPriors(x), chains = 4){
  initList <- vector("list", length = chains)
  initList <- lapply(initList, initFunc, priors = priors, x = x)
  return(initList)
}
