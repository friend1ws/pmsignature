
getPMSignature_sparse <- function(G, K, isBG, BG0 = 0, numInit = 10) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }

  N <- G[[1]];
  fdim <- G[[2]];
  # M <- prod(fdim);
  # N <- nrow(G);
  
  tempL <- -Inf;
  tempPar <- c();
  for (kkk in 1:numInit) {
    
    F <- array(0, c(varK, length(fdim), max(fdim)));
    Q <- matrix(0, K, N);
    
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        F[k,kk,1:fdim[kk]] <- rgamma(fdim[kk], rep(1, fdim[kk]));
        F[k,kk,1:fdim[kk]] <- F[k,kk,1:fdim[kk]] / sum(F[k,kk,1:fdim[kk]]);
      }
    }
    

    Q <- matrix(rgamma(N * K, 1, 1), K, N);
    Q <- sweep(Q, 2, apply(Q, 2, sum), `/`)
    
    p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(t(Q)), K, N));
    Y <- list(G, K, isBG, BG0);  
    
    res1 <- turboEM::turboem(par=p0, y=Y, fixptfn=updatePMSParam_sparse, objfn=calcPMSLikelihood_sparse, method=c("squarem"), pconstr=PMSboundary_sparse(Y), control.run = list(convtype = "objfn", tol = 1e-4));
    
    print(c(kkk, res1$itr, res1$runtime[3], res1$value.objfn));
    
    if (res1$value.objfn > tempL) {
      tempL <- res1$value.objfn;
      tempPar <- res1$par;
    }
    
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- N * (K - 1);
  F <- convertFromTurbo_F(tempPar[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(tempPar[(lenF + 1):(lenF + lenQ)], K, N);
  dim(F) <- c(varK, length(fdim), max(fdim));
  dim(Q) <- c(N, K);
  
  return(list(F, Q, tempL))
  
}


updatePMSParam_sparse <- function(p, y) {
  
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
  lenQ <- (K - 1) * N;
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N);
  
  dim(Q) <- c(N, K);
  Q <- t(Q);
  ####################
  # E-step
  Theta <- updateTheta_NormalizedSparseC(as.vector(patternList), as.vector(sparseCount), as.vector(F), as.vector(Q), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, F0);
  # L <- updateLikelihoodC(vTheta, as.vector(G), K, N, M);
  # Theta <- normalizeTheta(vTheta, K, N, M);
  dim(Theta) <- c(K, samplePatternNum);
  
  ####################
  # M-step 
  F <- updateMstepFSparseC(as.vector(patternList), as.vector(sparseCount), as.vector(Theta), fdim, K, sampleNum, patternNum, samplePatternNum, isBG);
  Q <- updateMstepQSparseC(as.vector(patternList), as.vector(sparseCount), as.vector(Theta), K, sampleNum, patternNum, samplePatternNum);
  #########################################
  dim(F) <- c(K, length(fdim), max(fdim));
  dim(Q) <- c(K, N);
  Q <- t(Q);
  
  return(c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, N)));
  
}



calcPMSLikelihood_sparse <- function(p, y) {
  
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
  lenQ <- (K - 1) * N;
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N);
  
  dim(Q) <- c(N, K);
  Q <- t(Q);
  ####################
  
  return(getLogLikelihoodSparseC(as.vector(patternList), as.vector(sparseCount), as.vector(F), as.vector(Q), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, F0));

}


PMSboundary_sparse <- function(y) {
  
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
  lenQ <- (K - 1) * N;

  
  function(p) {
    return(all(boundaryTurbo_F(p[1:lenF], fdim, varK), boundaryTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N)));
  }
  
}

