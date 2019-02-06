#' Determine if a Metropolis–Hastings step should be accepted
#'
#' \code{acceptProposal} is a utility function to determine if a proposal should
#' be accepted in a Metropolis or Metropolis-Hastings step. This is shamelessly
#' stolen from \link[overture]{AcceptProposal} in the \code{overture} package.
#'
#' The function uses the Metropolis choice for a Metropolis/Metropolis-Hastings
#' sampler, which accepts a proposed value \eqn{x'} with probability \deqn{
#' A(x', x) = min(1, P(x')/P(x) g(x|x')/g(x'|x)) } where \eqn{P(x)} is the
#' target distribution and \eqn{g(x'|x)} is the proposal distribution.
#'
#' @param logCurr log density of the target at the current value,
#'   \eqn{log(P(x))}
#' @param logProp log density of the target at the proposed value,
#'   \eqn{log(P(x'))}
#' @param logCurrToProp log of transition distribution from current value to
#'   proposed value, \eqn{log(g(x'|x))}
#' @param logPropToCurr log of transition distribution from proposed value to
#'   current value, \eqn{log(g(x|x'))}
#' @return \code{TRUE/FALSE} for whether the proposal should be accepted or
#'   rejected, respectively
#' @examples
#' curr <- rnorm(1, 0, 1)
#' prop <- curr + runif(1, -0.5, 0.5)
#' acceptProposal(logCurr = dnorm(curr, 0, 1, log = TRUE),
#'                logProp = dnorm(prop, 0, 1, log = TRUE))
#' @export
acceptProposal <- function(logCurr, logProp, logCurrToProp=0,
                           logPropToCurr=0) {
  u <- stats::runif(1)
  log(u) <= (logProp - logCurr + logPropToCurr - logCurrToProp)
}
