#' Calculate the log posterior density.
#'
#' \code{logPost} returns the log of the posterior density up to a constant.
#'
#' This calculates the log of the posterior density in the \code{bcgp}
#' setting up to a constant.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param y A vector containing the observed response values in the training
#' set.
#' @param params A named vector containing the current value of each of
#' the parameters.
#' @param priors A named vector containing the prior information. In the
#' \code{bcgp} setting this will generally result from \code{unlist(priors)}.
#' @param C An \emph{n x n} covariance matrix for the training data.
#' @param K An \emph{n x n} covariance matrix for the variance process.
#' at the training data locations.
#' @return A scalar
#' @examples
#' x <- matrix(runif(20, 0, 1), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' inits <- createInits(x, priors, chains = 1)
#' G <- getCorMat(x, inits[[1]]$rhoG)
#' L <- getCorMat(x, inits[[1]]$rhoL)
#' R <- combineCorMats(inits[[1]]$w, G, L)
#' C <- getCovMat(inits[[1]]$V, R, inits[[1]]$sig2eps)
#' K <- inits[[1]]$sig2V * getCorMat(x, inits[[1]]$rhoV) + diag(1e-10, length(y))
#' params <- unlist(inits[[1]])
#' priors <- unlist(priors)
#' logPost((x, y, params, priors, C, K))

#' @export
logPost <- function(x, y, params, priors, C, K){

  yMinusMu <- y - params["beta0"]
  RC <- chol(C)
  tmpC <- forwardsolve(t(RC), yMinusMu )

  paramNames <- names(params)
  priorNames <- names(priors)

  rhoG <- params[startsWith(paramNames, "rhoG")]
  rhoL <- params[startsWith(paramNames, "rhoL")]
  rhoV <- params[startsWith(paramNames, "rhoV")]
  V <- params[startsWith(paramNames, "V")]

  logVMinusMuV <- log(V) - params["muV"]
  RK <- try(chol(K), silent = TRUE)
  if(is.matrix(RK)){
    tmpK <- forwardsolve(t(RK), logVMinusMuV)
    tmpK2 <- sum(tmpK^2)
  }else{
    # tmpK2 <- try(as.numeric(t(logVMinusMuV) %*% try(solve(K, logVMinusMuV))))
    KinvlogVMinusMuV <- try(solve(K, logVMinusMuV), silent = TRUE)
    if(is.numeric(KinvlogVMinusMuV)){
      tmpK2 <- sum(logVMinusMuV * KinvlogVMinusMuV)
    }else{
      R <- svd(K)
      tmpK2 <- t(logVMinusMuV) %*% R$v %*% diag(1/R$d) %*% t(R$u) %*% logVMinusMuV
      # warning("The covariance matrix for the variance process is ill-conditioned.
      #         A possible solution (not necessarily a good solution) is to use less
      #         data.")
    }
  }

  rhoGAlpha <- priors[startsWith(priorNames, "rhoG.alpha")]
  rhoGBeta <- priors[startsWith(priorNames, "rhoG.beta")]
  rhoLAlpha <- priors[startsWith(priorNames, "rhoL.alpha")]
  rhoLBeta <- priors[startsWith(priorNames, "rhoL.beta")]
  rhoVAlpha <- priors[startsWith(priorNames, "rhoV.alpha")]
  rhoVBeta <- priors[startsWith(priorNames, "rhoV.beta")]

  like <- -0.5*logDet(C) - 0.5 * sum(tmpC^2) # = -.5*logDet(C) -
  #   0.5 * t(yMinusMu) %*% solve(C) %*% yMinusMu
  # This is the fastest way I've found

  ## TODO: Check this and make sure the part for V is right.
  ## Have I gotten the Jacobian right? Check line 73
  prior <- (priors["w.alpha"] - 1) * log(params["w"] - priors["w.lower"]) +
    (priors["w.beta"] - 1) * log(priors["w.upper"] - params["w"]) +
    sum((rhoLAlpha - 1) * log(rhoL)) + sum((rhoLBeta - 1) * log(rhoG - rhoL)) +
    sum((rhoGAlpha - rhoLAlpha - rhoLBeta) * log(rhoG)) +
    sum((rhoGBeta - 1) * log(1 - rhoG)) +
    (priors["sig2eps.alpha"] - 1) * log(params["sig2eps"]) -
    params["sig2eps"]/priors["sig2eps.beta"] -
    0.5 * logDet(K) - 0.5 * tmpK2 - sum(log(V)) -
    1/(2*priors["muV.sig2"]) * (params["muV"] - priors["muV.betaV"])^2 +
    sum((rhoVAlpha - 1) * log(rhoV)) + sum((rhoVBeta - 1) * log(1 - rhoV)) -
    (priors["sig2V.alpha"] + 1) * log(params["sig2V"]) -
    1/(params["sig2V"]*priors["sig2V.beta"])

  logPost <- like + prior

  return(as.numeric(logPost))
}

