#' Obtain the parameters for mutation signatures and memberships
#' 
#' @param G the mutation data processed in the function (readMutFile or readRawMutfeatFile)
#' @param K the number of mutation signatures
#' @param isBG the logical value showing whether a background mutaiton features is included or not
#' @param BG0 a background mutaiton features
#' @param numInit the number of performing calculations with different initial values
#' @useDynLib pmsignature
#' @importFrom Rcpp sourceCpp
#' @importFrom turboEM turboem
#' @export
getPMSignature <- function(G, K, isBG = FALSE, BG0 = 0, numInit = 10) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }

  sampleNum <- G[[1]];
  fdim <- G[[2]];
  
  tempL <- -Inf;
  tempPar <- c();
  for (kkk in 1:numInit) {
    
    F <- array(0, c(varK, length(fdim), max(fdim)));
    
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        F[k,kk,1:fdim[kk]] <- rgamma(fdim[kk], rep(1, fdim[kk]));
        F[k,kk,1:fdim[kk]] <- F[k,kk,1:fdim[kk]] / sum(F[k,kk,1:fdim[kk]]);
      }
    }
    

    Q <- matrix(rgamma(sampleNum * K, 1, 1), K, sampleNum);
    Q <- sweep(Q, 2, apply(Q, 2, sum), `/`)
    
    p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(t(Q)), K, sampleNum));
    Y <- list(G, K, isBG, BG0);  
    
    res1 <- turboem(par=p0, y=Y, fixptfn=updatePMSParam, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y), control.run = list(convtype = "objfn", tol = 1e-4));
    
    print(c(kkk, res1$itr, res1$runtime[3], res1$value.objfn));
    
    if (res1$value.objfn > tempL) {
      tempL <- res1$value.objfn;
      tempPar <- res1$par;
    }
    
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- sampleNum * (K - 1);
  F <- convertFromTurbo_F(tempPar[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(tempPar[(lenF + 1):(lenF + lenQ)], K, sampleNum);
  dim(F) <- c(varK, length(fdim), max(fdim));
  dim(Q) <- c(sampleNum, K);
  
  return(list(F, Q, tempL))
  
}


#' Obtain the standard error estimates for parameters for mutation signatures and memberships
#' 
#' @param G the matrix of mutation feature data
#' @param K the number of mutation signatures
#' @param isBG the logical value showing whether a background mutaiton features is included or not
#' @param BG0 a background mutaiton features
#' @param F0 the initial value for the parameter of mutation signatures used for bootstraped parameter estimations
#' @param Q0 the initial value for the parameter of memberships used for bootstraped parameter estimations
#' @param bootNum the number of performing bootstrap calculations
#' @export
bootPMSignature <- function(G, K, isBG = FALSE, BG0 = 0, F0, Q0, bootNum = 0) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  sampleNum <- G[[1]];
  fdim <- G[[2]];
  
  sqF <- array(0, c(bootNum, varK, length(fdim), max(fdim)));
  sqQ <- array(0, c(bootNum, N, K));
  
  for (bbb in 1:bootNum) {
    
    ##########
    # This part is under construction!!!!
    bootG <- matrix(0, sampleNum, M);
    for (n in 1:sampleNum) {
      tempG <- table(sample(1:M, sum(G[n,]), replace=TRUE, prob=G[n,] / sum(G[n,])));
      bootG[n,as.integer(names(tempG))] <- as.integer(tempG);
    }
    ##########
    
    p0 <- c(convertToTurbo_F(as.vector(F0), fdim, K, isBG), convertToTurbo_Q(as.vector(Q0), K, sampleNum));
    Y <- list(bootG, fdim, K, N, M, isBG, BG0);  
    
    res1 <- turboEM::turboem(par=p0, y=Y, fixptfn=updatePMSParam, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y), control.run = list(convtype = "objfn", tol = 1e-4));
    
    lenF <- varK * (sum(fdim) - length(fdim));
    lenQ <- sampleNum * (K - 1);
    F <- convertFromTurbo_F(res1$par[1:lenF], fdim, K, isBG);
    Q <- convertFromTurbo_Q(res1$par[(lenF + 1):(lenF + lenQ)], K, sampleNum);
    dim(F) <- c(varK, length(fdim), max(fdim));
    dim(Q) <- c(N, K);
    
    for (k in 1:varK) {
      sqF[bbb,k,,] <- (F[k,,] - F0[k,,])^2;
    }
    
    for (n in 1:sampleNum) {
      sqQ[bbb,,] <- (Q[n,] - Q0[n,])^2;
    }
    
    print(c(bbb, res1$itr, res1$runtime[3], res1$value.objfn));
    
  }
  
  return(list(sqF, sqQ))
  
}


#' A function for updating parameters using EM-algorithm
#' 
#' @param p this variable includes the parameters for mutation signatures and membership parameters
#' @param y this variable includes the information on the mutation features, 
#' the number of mutation signatures specified and so on
#' @export
updatePMSParam <- function(p, y) {
  
  sampleNum <- y[[1]][[1]];
  fdim <- y[[1]][[2]];
  patternList <- y[[1]][[3]];
  sparseCount <- y[[1]][[4]];
  K <- y[[2]];
  isBG <- y[[3]];
  F0 <- y[[4]];  
  
  patternNum <- ncol(patternList);
  samplePatternNum <- ncol(sparseCount);

  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- (K - 1) * sampleNum;
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, sampleNum);
  
  dim(Q) <- c(sampleNum, K);
  Q <- t(Q);
  ####################
  # E-step
  Theta <- updateTheta_NormalizedC(as.vector(patternList), as.vector(sparseCount), as.vector(F), as.vector(Q), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, F0);
  dim(Theta) <- c(K, samplePatternNum);
  
  ####################
  # M-step 
  F_Q <- updateMstepFQC(as.vector(patternList), as.vector(sparseCount), as.vector(Theta), fdim, K, sampleNum, patternNum, samplePatternNum, isBG);
  #########################################
  F <- F_Q[1:(varK * length(fdim) * max(fdim))];
  Q <- F_Q[(varK * length(fdim) * max(fdim) + 1):(varK * length(fdim) * max(fdim) + K * sampleNum)];
  dim(F) <- c(varK, length(fdim), max(fdim));
  dim(Q) <- c(K, sampleNum);
  Q <- t(Q);
  
  return(c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, sampleNum)));
  
}


#' A function for calculating the log-likelihood from the data and parameters
#' 
#' @param p this variable includes the parameters for mutation signatures and membership parameters
#' @param y this variable includes the information on the mutation features, 
#' the number of mutation signatures specified and so on
#' @export
calcPMSLikelihood <- function(p, y) {
  
  sampleNum <- y[[1]][[1]];
  fdim <- y[[1]][[2]];
  patternList <- y[[1]][[3]];
  sparseCount <- y[[1]][[4]];
  K <- y[[2]];
  isBG <- y[[3]];
  F0 <- y[[4]];  
  
  patternNum <- ncol(patternList);
  samplePatternNum <- ncol(sparseCount);
  
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- (K - 1) * sampleNum;
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, sampleNum);
  
  dim(Q) <- c(sampleNum, K);
  Q <- t(Q);
  ####################
  
  return(getLogLikelihoodC(as.vector(patternList), as.vector(sparseCount), as.vector(F), as.vector(Q), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, F0));

}


#' A functional for generating the function checking the parameter (p) is within the restricted conditions or not
#' 
#' @param y this variable includes the information on the mutation features, 
#' the number of mutation signatures specified and so on
#' @export
PMSboundary <- function(y) {
  
  sampleNum <- y[[1]][[1]];
  fdim <- y[[1]][[2]];
  patternList <- y[[1]][[3]];
  sparseCount <- y[[1]][[4]];
  K <- y[[2]];
  isBG <- y[[3]];
  F0 <- y[[4]];  
  
  patternNum <- ncol(patternList);
  samplePatternNum <- ncol(sparseCount);
  
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- (K - 1) * sampleNum;

  
  function(p) {
    return(all(boundaryTurbo_F(p[1:lenF], fdim, varK), boundaryTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, sampleNum)));
  }
  
}

