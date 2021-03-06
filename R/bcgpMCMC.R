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
#' x <- matrix(runif(20, 0, 1), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' inits <- createInits(x, priors, chains = 4)
#' numUpdates <- 3
#' numAdapt <- 500
#' burnin <- 500
#' nmcmc <- 5000
#' chains <- 4
#' cores <- 1
#' bcgpMCMC(x, y, priors, inits, numUpdates, numAdapt,
#'          burnin, nmcmc, chains, cores)

bcgpMCMC  <- function(x, y, priors, inits, numUpdates, numAdapt,
                      burnin, nmcmc, chains = 4, cores = 1){

  nTrain <- nrow(x)
  d <- ncol(x)
  iterations <- numUpdates*numAdapt + burnin + nmcmc
  epsV <- 1e-10
  tau2 <- 0.08
  priorVec <- unlist(priors)

  rhoNames <- rep("rho", d)
  rhoGNames <- paste0(rhoNames, paste0("G", 1:d))
  rhoLNames <- paste0(rhoNames, paste0("L", 1:d))
  rhoVNames <- paste0(rhoNames, paste0("V", 1:d))
  rm(rhoNames)

  samples <- vector("list", chains)
  warmup <- vector("list", chains)
  acceptances <- vector("list", chains)

  if(nTrain >= 20){
    nProp <- 10 # number of training locations to sample near the focal point for V proposals
    m <- ceiling(nTrain/nProp) + 1 # Number of times to cycle through for V proposals
  }

  onesNTrain <- rep(1, nTrain)

  ## TODO: Currently not parallelized. Need to make that happen
  for(i in seq_len(chains)){

    ## Side note: working with named matrices and vectors is faster
    ## than working with data frames or lists

    allDraws <- matrix(NA, nrow = iterations, ncol = 5 + 3*d + nTrain)
    row1 <- unlist(inits[[i]])
    colnames(allDraws) <- names(row1)
    allDraws[1, ] <- row1

    allAcceptances <- matrix(0, nrow = iterations, ncol = 5 + 3*d + nTrain)
    colnames(allAcceptances) <- names(row1)
    allAcceptances[1, ] <- 1

    logProb <- vector(mode = "numeric", length = iterations)

    if(d == 1){
      colnames(allDraws)[startsWith(colnames(allDraws), "rho")] <-
        paste0(colnames(allDraws)[startsWith(colnames(allDraws), "rho")],"1")
      colnames(allAcceptances)[startsWith(colnames(allAcceptances), "rho")] <-
        paste0(colnames(allAcceptances)[startsWith(colnames(allAcceptances), "rho")],"1")
    }

    rm(row1)

    propWidths <- c(0.4, rep(0.07, d), rep(0.03, d), 1e-3, rep(0.25, d))
    names(propWidths) <- c("w", rhoGNames, rhoLNames, "sig2eps", rhoVNames)
    calAccept <- rep(0, length(propWidths)) # container for acceptances
                                            # during prop width calibration
    names(calAccept) <- names(propWidths)


    for(j in 2:iterations){

      if(j %% 100 == 0) print(j)
      G <- getCorMat(x, allDraws[j - 1, rhoGNames])
      L <- getCorMat(x, allDraws[j - 1, rhoLNames])
      R <- combineCorMats(allDraws[j - 1, "w"], G, L)
      C <- getCovMat(allDraws[j - 1, startsWith(colnames(allDraws), "V")], R,
                     allDraws[j - 1, "sig2eps"])
      K <- allDraws[j - 1, "sig2V"] * getCorMat(x, allDraws[j - 1, rhoVNames]) +
        diag(epsV, nTrain)

      allDraws[j, ] = allDraws[j - 1, ]

      ## Sample for beta0. Gibbs step
      beta0Calcs <- x1Ainvx2(x1 = list(onesNTrain, onesNTrain), A = C,
                             x2 = list(onesNTrain, y))
      oneCinvone <- beta0Calcs[1]
      oneCinvY <- beta0Calcs[2]

      allDraws[j, "beta0"] <- rnorm(1, oneCinvY/oneCinvone, sqrt(1/oneCinvone))
      allAcceptances[j, "beta0"] <- 1

      ## Propose for w, accept or reject in Metropolis step
      wP <- allDraws[j, "w"] + runif(1, -propWidths["w"], propWidths["w"])
      if(wP < priors$w$lower || wP > priors$w$upper){
        accept <- FALSE
      }else{
        RP <- combineCorMats(wP, G, L)
        CP <- getCovMat(allDraws[j, startsWith(colnames(allDraws), "V")],
                        RP, allDraws[j, "sig2eps"])
        paramsP <- allDraws[j, ]
        paramsP["w"] <- wP
        accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                   priorVec, C, K),
                                 logProp = logPost(x, y, params = paramsP,
                                                   priorVec, CP, K))


      }

      if(accept){
        allDraws[j, "w"] <- wP
        R <- RP
        C <- CP
        allAcceptances[j, "w"] <- 1
        calAccept["w"] <- calAccept["w"] + 1
      }

      ## Propose for sig2eps, accept or reject in Metropolis step
      sig2epsP <- allDraws[j, "sig2eps"] + runif(1, -propWidths["sig2eps"],
                                                 propWidths["sig2eps"])
      if(isTRUE(sig2epsP <= 0)){
        accept <- FALSE
      }else{
        CP <- getCovMat(allDraws[j, startsWith(colnames(allDraws), "V")],
                        R, sig2epsP)
        paramsP <- allDraws[j, ]
        paramsP["sig2eps"] <- sig2epsP
        accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                   priorVec, C, K),
                                 logProp = logPost(x, y, params = paramsP,
                                                   priorVec, CP, K))
      }

      if(isTRUE(accept)){
        allDraws[j, "sig2eps"] <- sig2epsP
        C <- CP
        allAcceptances[j, "sig2eps"] <- 1
        calAccept["sig2eps"] <- calAccept["sig2eps"] + 1
      }

      ## Propose for rhoGk, accept or reject in Metropolis step
      ## Then propose for rhoLk given rhoGk, accept or reject in Metropolis step
      for(k in 1:d){

        rhoGx <- paste0("rhoG", k)
        rhoLx <- paste0("rhoL", k)

        rhoGkP <- allDraws[j, rhoGx] + runif(1, -propWidths[rhoGx],
                                             propWidths[rhoGx])

        if(rhoGkP < allDraws[j, rhoLx] || rhoGkP > 1){
          accept <- FALSE
        }else{
          paramsP <- allDraws[j, ]
          paramsP[rhoGx] <- rhoGkP

          GP <- getCorMat(x, paramsP[rhoGNames])
          RP <- combineCorMats(allDraws[j, "w"], GP, L)
          CP <- getCovMat(allDraws[j, startsWith(colnames(allDraws), "V")],
                          RP, allDraws[j, "sig2eps"])

          accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                     priorVec, C, K),
                                   logProp = logPost(x, y, params = paramsP,
                                                     priorVec, CP, K))

        }

        if(accept){
          allDraws[j, rhoGx] <- rhoGkP
          G <- GP
          R <- RP
          C <- CP
          allAcceptances[j, rhoGx] <- 1
          calAccept[rhoGx] <- calAccept[rhoGx] + 1
        }



        rhoLkP <- allDraws[j, rhoLx] + runif(1, -propWidths[rhoLx],
                                             propWidths[rhoLx])

        if(rhoLkP > allDraws[j, rhoGx] || rhoLkP < 0){
          accept <- FALSE
        }else{
          paramsP <- allDraws[j, ]
          paramsP[rhoLx] <- rhoLkP

          LP <- getCorMat(x, paramsP[rhoLNames])
          RP <- combineCorMats(allDraws[j, "w"], G, LP)
          CP <- getCovMat(allDraws[j, startsWith(colnames(allDraws), "V")],
                          RP, allDraws[j, "sig2eps"])

          accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                     priorVec, C, K),
                                   logProp = logPost(x, y, params = paramsP,
                                                     priorVec, CP, K))

        }

        if(accept){
          allDraws[j, rhoLx] <- rhoLkP
          L <- LP
          R <- RP
          C <- CP
          allAcceptances[j, rhoLx] <- 1
          calAccept[rhoLx] <- calAccept[rhoLx] + 1
        }


      }

      ## Sample for muV and sig2V (one-at-a-time). Gibbs steps
      Rt <- getCorMat(x, allDraws[j, rhoVNames]) + diag(epsV, nTrain)
      W <- log(allDraws[j, startsWith(colnames(allDraws), "V")])
      WMinusMuV <- W - allDraws[j, "muV"]
      muVSig2VCalcs <- x1Ainvx2(x1 = list(onesNTrain, onesNTrain, WMinusMuV),
                                A = Rt,
                                x2 = list(onesNTrain, W, WMinusMuV))
      oneRtinvone <- muVSig2VCalcs[1]
      oneRtinvW <- muVSig2VCalcs[2]
      WMinusMuVRtinvWMinusMuV <- muVSig2VCalcs[3]


      condMeanNum <- priorVec["muV.betaV"]/priorVec["muV.sig2"] +
        oneRtinvW/allDraws[j, "sig2V"]
      condMeanDenom <- 1/priorVec["muV.sig2"] + oneRtinvone/allDraws[j, "sig2V"]
      condMean <- condMeanNum/condMeanDenom # The conditional mean for muV
      condVar <- 1/condMeanDenom            # the conditional variance for muV
      allDraws[j, "muV"] <- rnorm(1, condMean, sqrt(condVar))
      allAcceptances[j, "muV"] <- 1

      newAlpha <- priorVec["sig2V.alpha"] + nTrain/2
      newBeta <- 1/(0.5 * WMinusMuVRtinvWMinusMuV + 1/priorVec["sig2V.beta"])
      allDraws[j, "sig2V"] <- 1/rgamma(1, shape = newAlpha, scale = newBeta)
      allAcceptances[j, "sig2V"] <- 1
      K <- allDraws[j, "sig2V"] * getCorMat(x, allDraws[j, rhoVNames]) +
        diag(epsV, nTrain)


      ## Propose for rhoV, accept or reject in Metropolis step
      for(k in 1:d){

        rhoVx <- paste0("rhoV", k)

        rhoVkP <- allDraws[j, rhoVx] + runif(1, -propWidths[rhoVx],
                                             propWidths[rhoVx])

        if(rhoVkP < 0 || rhoVkP > 1){
          accept <- FALSE
        }else{
          paramsP <- allDraws[j, ]
          paramsP[rhoVx] <- rhoVkP

          KP <- allDraws[j, "sig2V"]*getCorMat(x, paramsP[rhoVNames]) +
            diag(epsV, nTrain)

          accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                     priorVec, C, K),
                                   logProp = logPost(x, y, params = paramsP,
                                                     priorVec, CP, K))

        }

        if(accept){
          allDraws[j, rhoVx] <- rhoVkP
          K <- KP
          allAcceptances[j, rhoVx] <- 1
          calAccept[rhoVx] <- calAccept[rhoVx] + 1
        }
      }

      ## An attempt at adapting all the proposal widths at the same time
      if(j %% numAdapt == 0 && j < (numUpdates*numAdapt + 1)){
        acceptRate <- calAccept/numAdapt
        propWidths <- propWidths*acceptRate/0.33
        calAccept <- rep(0, length(propWidths))
        names(calAccept) <- names(propWidths)
      }


      ## Sample for sig2(x) (V)
      if(nTrain >= 20){ # GET VALUES FOR SIG2X FOR "nProp" LOCATIONS AT A TIME
        VC <- allDraws[j, startsWith(colnames(allDraws), "V")]
        for(k in seq_len(m)){

          focalPoint <- runif(d, 0, 1)
          distMat <- as.matrix(dist(rbind(focalPoint, x),
                                    method = "euclidean",
                                    diag = FALSE, upper = FALSE))[, 1]
          idxIn <- sort(order(distMat[-1])[1:nProp])
          trainIn <- as.matrix(x[idxIn, ])
          trainOut <- as.matrix(x[-idxIn, ])

          allTrain <- rbind(trainIn, trainOut)

          VCIn <- VC[idxIn]
          propMeanWIn <- log(VCIn)

          KW <- tau2 * getCorMat(allTrain, allDraws[j, rhoVNames]) + diag(epsV, nTrain)

          KWIn <- KW[seq_len(nProp), seq_len(nProp)]
          KWOut <- KW[-seq_len(nProp), -seq_len(nProp)]
          KWBetween <- KW[-seq_len(nProp), seq_len(nProp)]

          RKWOut <- try(chol(KWOut), silent = TRUE)
          if(is.matrix(RKWOut)){
            halfVar <- forwardsolve(t(RKWOut), KWBetween)
            propVar <- KWIn - t(halfVar) %*% halfVar
          }else{

            KWOinvKWB <- try(solve(KWOut, KWBetween), silent = TRUE)
            if(is.numeric(KWOinvKWB)){
              propVar <- KWIn - t(KWBetween) %*% KWOinvKWB
            }else{
              svdKWO <- svd(KWOut)
              propVar <- KWIn - t(KWBetween) %*% svdKWO$v %*% diag(1/svdKWO$d) %*%
                t(svdKWO$u) %*% KWBetween
            }
          }

          VP <- allDraws[j, startsWith(colnames(allDraws), "V")]
          VP[idxIn] <- exp(MASS::mvrnorm(1, propMeanWIn, propVar))

          CP <- getCovMat(VP, R, allDraws[j, "sig2eps"])
          paramsP <- allDraws[j, ]
          paramsP[startsWith(colnames(allDraws), "V")] <- VP

          curr <- logPost(x, y, params = allDraws[j, ],
                         priorVec, C, K)
          prop <- logPost(x, y, params = paramsP,
                         priorVec, CP, K)

          accept <- acceptProposal(logCurr = curr,
                                   logProp = prop,
                                   logCurrToProp = -sum(log(VP)),
                                   logPropToCurr = -sum(log(allDraws[j, startsWith(colnames(allDraws), "V")])))

          if(accept){
            allDraws[j, startsWith(colnames(allDraws), "V")] <- VP
            C <- CP
            logProb[j] <- prop
          }else{
            logProb[j] <- curr
          }

        }

        allAcceptances[j, startsWith(colnames(allDraws), "V")] <-
          (allDraws[j, startsWith(colnames(allDraws), "V")] != VC)



      }else{ # GET VALUES FOR SIG2X FOR THE ENTIRE VECTOR AT THE SAME TIME
        KW <- tau2 * getCorMat(x, allDraws[j, rhoVNames]) + diag(epsV, nTrain)
        VP <- exp(MASS::mvrnorm(1, log(allDraws[j, startsWith(colnames(allDraws), "V")]),
                                KW))
        CP <- getCovMat(VP, R, allDraws[j, "sig2eps"])
        paramsP <- allDraws[j, ]
        paramsP[startsWith(colnames(allDraws), "V")] <- VP

        # cholKWR <- try(chol(C), silent = TRUE)
        # if(is.matrix(cholKWR)){
        #   logVMinusLogVP <- log(VP) - log(allDraws[j, startsWith(colnames(allDraws), "V")])
        #   logVPMinusLogV <- -logVMinusLogVP
        #   tmpKWN <- forwardsolve(t(cholKWR), logVMinusLogVP)
        #   tmpKWD <- forwardsolve(t(cholKWR), logVPMinusLogV)
        #
        #   propToCurr <- -0.5*sum(tmpKWN^2) -
        #     sum(log(allDraws[j, startsWith(colnames(allDraws), "V")]))
        #   currToProp <- -0.5*sum(tmpKWD^2) - sum(log(VP))
        # }else{
        #   stop("The covariance matrix is ill-conditioned. Possible solutions (not
        #        necessarily good solutions) are to use less data or to set the priors
        #        for 'sig2eps' in such a way that 'sig2eps' will be larger.")
        # }

        curr <- logPost(x, y, params = allDraws[j, ],
                        priorVec, C, K)
        prop <- logPost(x, y, params = paramsP,
                        priorVec, CP, K)

        accept <- acceptProposal(logCurr = curr,
                                 logProp = prop,
                                 logCurrToProp = -sum(log(VP)),
                                 logPropToCurr = -sum(log(allDraws[j, startsWith(colnames(allDraws), "V")])))

        if(accept){
          allDraws[j, startsWith(colnames(allDraws), "V")] <- VP
          C <- CP
          logProb[j] <- prop
          allAcceptances[j, startsWith(colnames(allDraws), "V")] <- 1
        }else{
          logProb[j] <- curr
        }
      }

    }


    warmup[[i]] <- as.list(data.frame(allDraws[1:(numUpdates*numAdapt + burnin),
                                    !startsWith(colnames(allDraws), "V")]))
    warmup[[i]]$V <- allDraws[1:(numUpdates*numAdapt + burnin),
                                    startsWith(colnames(allDraws), "V")]
    warmup[[i]]$lp__ <- logProb[1:(numUpdates*numAdapt + burnin)]

    samples[[i]] <- as.list(data.frame(allDraws[(iterations - nmcmc + 1):iterations,
                                               !startsWith(colnames(allDraws), "V")]))
    samples[[i]]$V <- allDraws[(iterations - nmcmc + 1):iterations,
                              startsWith(colnames(allDraws), "V")]
    samples[[i]]$lp__ <- logProb[(iterations - nmcmc + 1):iterations]

    # warmup[[i]] <- as.list(data.frame(allDraws[1:(numUpdates*numAdapt + burnin), ]))
    # samples[[i]] <- as.list(data.frame(allDraws[(iterations - nmcmc + 1):iterations, ]))
    acceptances[[i]] <- as.list(
      data.frame(allAcceptances[(iterations - nmcmc + 1):iterations, ]))

  }

  params <- c("beta0", "w", "rhoG", "rhoL", "sig2eps", "V",
              "muV", "rhoV", "sig2V", "lp__")
  sim <- list(samples = samples,
              warmup = warmup,
              acceptances = acceptances,
              chains = chains,
              numUpdates = numUpdates,
              numAdapt = numAdapt,
              burnin = burnin,
              nmcmc = nmcmc,
              pars_oi = params)

  bfit <- new("bcgpfit",
              model_pars = params,
              par_dims = list(beta0 = 1,
                              w = 1,
                              rhoG = d,
                              rhoL = d,
                              sig2eps = 1,
                              V = nTrain,
                              muV = 1,
                              rhoV = d,
                              sig2V = 1,
                              lp__ = 1),
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
