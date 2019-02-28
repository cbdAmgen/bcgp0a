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

    if(d == 1){
      colnames(allDraws)[startsWith(colnames(allDraws), "rho")] <-
        paste0(colnames(allDraws)[startsWith(colnames(allDraws), "rho")],"1")
      colnames(allAcceptances)[startsWith(colnames(allAcceptances), "rho")] <-
        paste0(colnames(allAcceptances)[startsWith(colnames(allAcceptances), "rho")],"1")
    }

    rm(row1)

    # G <- getCorMat(x, inits[[i]]$rhoG)
    # L <- getCorMat(x, inits[[i]]$rhoL)
    # R <- combineCorMats(inits[[i]]$w, G, L)
    # C <- getCovMat(inits[[i]]$V, R, inits[[i]]$sig2eps)
    # K <- inits[[i]]$sig2V * getCorMat(x, inits[[i]]$rhoV) + diag(epsV, nTrain)

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
      cholCR <- try(chol(C), silent = TRUE)
      if(is.matrix(cholCR)){
        tmpC <- forwardsolve(t(cholCR), rep(1, nTrain))
        tmpC2 <- forwardsolve(t(cholCR), y)
        oneCinvone <- sum(tmpC^2)
        oneCinvY <- sum(tmpC * tmpC2)
        allDraws[j, "beta0"] <- rnorm(1, oneCinvY/oneCinvone, sqrt(1/oneCinvone))
        allAcceptances[j, "beta0"] <- 1
        rm(tmpC, tmpC2)
      }else{
        stop("The covariance matrix is ill-conditioned. Possible solutions (not
             necessarily good solutions) are to use less data or to set the priors
             for 'sig2eps' in such a way that 'sig2eps' will be larger.")
      }

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
      cholRtR <- try(chol(Rt), silent = TRUE)
      if(is.matrix(cholRtR)){
        tmpRt <- forwardsolve(t(cholRtR), rep(1, nTrain))
        tmpRt2 <- forwardsolve(t(cholRtR),
                               log(allDraws[j, startsWith(colnames(allDraws), "V")]))
        oneRtinvone <- sum(tmpRt^2)
        oneRtinvW <- sum(tmpRt * tmpRt2)

        condmNum <- priorVec["muV.betaV"]/priorVec["muV.sig2"] +
          oneRtinvW/allDraws[j, "sig2V"]
        condmDenom <- 1/priorVec["muV.sig2"] + oneRtinvone/allDraws[j, "sig2V"]
        condm <- condmNum/condmDenom # The conditional mean for muV
        condv <- 1/condmDenom        # the conditional variamce for muV
        allDraws[j, "muV"] <- rnorm(1, condm, sqrt(condv))
        allAcceptances[j, "muV"] <- 1


        WMinusMuV <- log(allDraws[j, startsWith(colnames(allDraws), "V")]) -
          allDraws[j, "muV"]
        tmpWMinusMuV <- forwardsolve(t(cholRtR), WMinusMuV)
        WMinusMuVRtinvWMinusMuV <- sum(tmpWMinusMuV^2)
        newAlpha <- priorVec["sig2V.alpha"] + nTrain/2
        newBeta <- 1/(0.5 * WMinusMuVRtinvWMinusMuV + 1/priorVec["sig2V.beta"])
        allDraws[j, "sig2V"] <- 1/rgamma(1, shape = newAlpha, scale = newBeta)
        allAcceptances[j, "sig2V"] <- 1
        K <- allDraws[j, "sig2V"] * getCorMat(x, allDraws[j, rhoVNames]) +
          diag(epsV, nTrain)

        rm(tmpRt, tmpRt2)
      }else{
        stop("The covariance matrix for the variance process is ill-conditioned.
               A possible solution (not necessarily a good solution) is to use less
               data.")
      }

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
        # print(paste0('The most recent acceptance rates are ', round(acceptRate,4)))
        # print(paste0('The old proposal widths are ', round(propWidths,4)))
        propWidths <- propWidths*acceptRate/0.33
        # print(paste0('The new proposal widths are ', round(propWidths,4)))
        calAccept <- rep(0, length(propWidths))
        names(calAccept) <- names(propWidths)
      }


      ## Sample for sig2(x) (V)
      if(nTrain >= 20){ # GET VALUES FOR SIG2X FOR "howManyClose" LOCATIONS AT A TIME
        for(k in seq_len(m)){

          focalPoint <- runif(d, 0, 1)
          distMat <- as.matrix(dist(rbind(focalPoint, x),
                                    method = "euclidean",
                                    diag = FALSE, upper = FALSE))[, 1]
          idxIn <- sort(order(distMat[-1])[1:nProp])
          trainIn <- as.matrix(x[idxIn, ])
          trainOut <- as.matrix(x[-idxIn, ])

          allTrain <- rbind(trainIn, trainOut)

          VC <- allDraws[j, startsWith(colnames(allDraws), "V")][idxIn]
          propMeanWIn <- log(VC[idxIn])

          KW <- tau2 * getCorMat(allTrain, allDraws[j, rhoVNames]) + diag(epsV, nTrain)

          KWIn <- KW[seq_len(nProp), seq_len(nProp)]
          KWOut <- KW[-seq_len(nProp), -seq_len(nProp)]
          KWBetween <- KW[-seq_len(nProp), seq_len(nProp)]

          RKWOut <- try(chol(KWOut), silent = TRUE)
          if(is.matrix(RKWOut)){
            halfVar <- forwardsolve(t(RKWOut), KWBetween)
            propVar <- KWIn - t(halfVar) %*% halfVar
          }else{
            propVar <- KWIn - t(KWBetween) %*% solve(KWOut, KWBetween)
          }

          # VP <- VC
          # VP[idxIn] <- exp(MASS::mvrnorm(1, propMeanWIn, propVar))

          # VP <- exp(MASS::mvrnorm(1, log(allDraws[j, startsWith(colnames(allDraws), "V")]),
          #                         KW))
          # CP <- getCovMat(VP, R, allDraws[j, "sig2eps"])
          # paramsP <- allDraws[j, ]
          # paramsP[startsWith(colnames(allDraws), "V")] <- VP

          # KIdxW = tau2*(RIdxW - RBetweenW'/(RNoIdxW + epsV*eye(size(RNoIdxW,1)))*RBetweenW) + epsV*eye(sum(idx~=0));
          # wIdx = randnorm(1,muIdxW',[],KIdxW);
          # vIdx = exp(wIdx);
          # VP = VC;
          # VP(idx~=0) = vIdx;
          #
          # RC = getRPred(wC,GC,LC);
          # CP = getCPred(VP,RC,sig2epsC);

        }
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

        accept <- acceptProposal(logCurr = logPost(x, y, params = allDraws[j, ],
                                                   priorVec, C, K),
                                 logProp = logPost(x, y, params = paramsP,
                                                   priorVec, CP, K),
                                 logCurrToProp = -sum(log(VP)),
                                 logPropToCurr = -sum(log(allDraws[j, startsWith(colnames(allDraws), "V")])))


        if(accept){
          allDraws[j, startsWith(colnames(allDraws), "V")] <- VP
          C <- CP
          allAcceptances[j, startsWith(colnames(allDraws), "V")] <- 1
        }
      }

      # %%% sample for sig2(x) %%%
      #   if size(train.xt,1) >= 20 %%%%%% THIS BLOCK GETS VALUES FOR SIG2X AT
      # %%%%%% TRAINING LOCATIONS "howManyClose" AT A
      # %%%%%% TIME
      #
      # for k = 1:numPropose
      #
      # midPoint = unifrnd(0,1,1,d);
      # distance = pdist2(train.xt,midPoint);
      # [~,index] = getNElements(distance,howManyClose);
      # idx = zeros(size(train.xt,1),1);
      # idx(index) = 1;
      #
      # trainIn = train.xt(idx~=0,:);
      # trainOut = train.xt(idx==0,:);
      #
      # muIdxW = log(VC(idx~=0));
      #
      # RIdxW = getGPredC(train.xt(idx~=0,:),rhoVC);
      # RNoIdxW = getGPredC(train.xt(idx==0,:),rhoVC);
      #
      # RBetweenW = ones(sum(idx==0),sum(idx~=0));
      # for j = 1:sum(idx==0),
      # for h = 1:sum(idx~=0),
      # for s = 1:size(train.xt,2)
      # RBetweenW(j,h) = RBetweenW(j,h) * rhoVC(s).^(16*((trainOut(j,s)-trainIn(h,s))^2));
      # end
      # end
      # end
      #
      # KIdxW = tau2*(RIdxW - RBetweenW'/(RNoIdxW + epsV*eye(size(RNoIdxW,1)))*RBetweenW) + epsV*eye(sum(idx~=0));
      # wIdx = randnorm(1,muIdxW',[],KIdxW);
      # vIdx = exp(wIdx);
      # VP = VC;
      # VP(idx~=0) = vIdx;
      #
      # RC = getRPred(wC,GC,LC);
      # CP = getCPred(VP,RC,sig2epsC);
      #
      # t_f = log(unifrnd(0,1)) < forVPred(VP,VC,CP,CC,train,muC,muVC,KC);
      # if t_f == 1
      # draws.Vt(i,:) = VP;
      # VC = VP;
      # accept.V(i,idx~=0) = 1;
      # else
      #   draws.Vt(i,:) = VC;
      # end
      #
      # end
      #
      # else %%% GET VALUES FOR SIG2X FOR THE ENTIRE VECTOR AT THE SAME TIME %%%
      #
      #   KW = tau2 * getGPredC(train.xt,rhoVC) + epsV*eye(size(CC,1));
      # WP = mvnrnd(log(VC)',KW);
      #   VP = exp(WP);
      #
      #   proposals.Vt(i,:) = VP;
      #   RC = getRPred(wC,GC,LC);
      #   CP = getCPred(VP,RC,sig2epsC);
      #
      #   t_f = log(unifrnd(0,1)) < forVPred(VP,VC,CP,CC,train,muC,muVC,KC);
      #   if t_f == 1
      #       draws.Vt(i,:) = VP';
      #             else
      #               draws.Vt(i,:) = VC;
      #             accept.V(i) = 0;
      #             end
      #             end
      #
      #             end


    }


    warmup[[i]] <- as.list(data.frame(allDraws[1:(numUpdates*numAdapt + burnin),
                                    !startsWith(colnames(allDraws), "V")]))
    warmup[[i]]$V <- allDraws[1:(numUpdates*numAdapt + burnin),
                                    startsWith(colnames(allDraws), "V")]
    samples[[i]] <- as.list(data.frame(allDraws[(iterations - nmcmc + 1):iterations,
                                               !startsWith(colnames(allDraws), "V")]))
    samples[[i]]$V <- allDraws[(iterations - nmcmc + 1):iterations,
                              startsWith(colnames(allDraws), "V")]

    # warmup[[i]] <- as.list(data.frame(allDraws[1:(numUpdates*numAdapt + burnin), ]))
    # samples[[i]] <- as.list(data.frame(allDraws[(iterations - nmcmc + 1):iterations, ]))
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
