
getPMSignature <- function(G, K, fdim, isBG, BG0, numInit) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  M <- prod(fdim);
  N <- nrow(G);
  
  tempL <- -Inf;
  tempPar <- c();
  for (kkk in 1:numInit) {
    
    F <- array(0, c(varK, length(fdim), max(fdim)));
    Q <- matrix(0, N, K);
    
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        F[k,kk,1:fdim[kk]] <- rgamma(fdim[kk], rep(1, fdim[kk]));
        F[k,kk,1:fdim[kk]] <- F[k,kk,1:fdim[kk]] / sum(F[k,kk,1:fdim[kk]]);
      }
    }
    
    for (i in 1:N) {
      Q[i,] <- rgamma(K, 1, 1);
      Q[i,] <- Q[i,] / sum(Q[i,]);
    }
    
    p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, N));
    Y <- list(G, fdim, K, N, M, isBG, BG0);  
    
    res1 <- turboem(par=p0, y=Y, fixptfn=updatePMSParam, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y), control.run = list(convtype = "objfn", tol = 1e-4));
    
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


updatePMSParam <- function(p, y) {
  
  G <- y[[1]];
  fdim <- y[[2]];
  K <- y[[3]];
  N <- y[[4]];
  M <- y[[5]];
  isBG <- y[[6]];
  F0 <- y[[7]];  
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- N * (K - 1);
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N);
  
  ####################
  # E-step
  Theta <- updateTheta_NormalizedC(as.vector(F), as.vector(Q), fdim, K, N, M, isBG, F0);
  # L <- updateLikelihoodC(vTheta, as.vector(G), K, N, M);
  # Theta <- normalizeTheta(vTheta, K, N, M);
  dim(Theta) <- c(N, K, M);
  
  ####################
  # M-step 
  F <- updateMstepFC(as.vector(Theta), as.vector(G), fdim, K, N, M, isBG);
  Q <- updateMstepQC(as.vector(Theta), as.vector(G), K, N, M);
  #########################################
  
  return(c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, N)));
  
}

bootPMSignature <- function(G, K, fdim, isBG, BG0, F0, Q0, bootNum) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  M <- prod(fdim);
  N <- nrow(G);
  
  sqF <- array(0, c(bootNum, varK, length(fdim), max(fdim)));
  sqQ <- array(0, c(bootNum, N, K));
  
  for (bbb in 1:bootNum) {
    
    bootG <- matrix(0, N, M);
    for (n in 1:N) {
      tempG <- table(sample(1:M, sum(G[n,]), replace=TRUE, prob=G[n,] / sum(G[n,])));
      bootG[n,as.integer(names(tempG))] <- as.integer(tempG);
    }
    
    p0 <- c(convertToTurbo_F(as.vector(F0), fdim, K, isBG), convertToTurbo_Q(as.vector(Q0), K, N));
    Y <- list(bootG, fdim, K, N, M, isBG, BG0);  
    
    res1 <- turboem(par=p0, y=Y, fixptfn=updatePMSParam, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y), control.run = list(convtype = "objfn", tol = 1e-4));
    
    lenF <- varK * (sum(fdim) - length(fdim));
    lenQ <- N * (K - 1);
    F <- convertFromTurbo_F(res1$par[1:lenF], fdim, K, isBG);
    Q <- convertFromTurbo_Q(res1$par[(lenF + 1):(lenF + lenQ)], K, N);
    dim(F) <- c(varK, length(fdim), max(fdim));
    dim(Q) <- c(N, K);
    
    for (k in 1:varK) {
      sqF[bbb,k,,] <- (F[k,,] - F0[k,,])^2;
    }
    
    for (n in 1:N) {
      sqQ[bbb,,] <- (Q[n,] - Q0[n,])^2;
    }
    
    print(c(bbb, res1$itr, res1$runtime[3], res1$value.objfn));
    
  }
  
  return(list(sqF, sqQ))
  
}



calcPMSLikelihood <- function(p, y) {
  
  G <- y[[1]];
  fdim <- y[[2]];
  K <- y[[3]];
  N <- y[[4]];
  M <- y[[5]];
  isBG <- y[[6]];
  BG0 <- y[[7]];  
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- N * (K - 1);
  F <- convertFromTurbo_F(p[1:lenF], fdim, K, isBG);
  Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N);
  
  
  ####################
  
  return(getLogLikelihoodC(as.vector(G), as.vector(F), as.vector(Q), fdim, K, N, M, isBG, BG0));
  
}


PMSboundary <- function(y) {
  
  G <- y[[1]];
  fdim <- y[[2]];
  K <- y[[3]];
  N <- y[[4]];
  M <- y[[5]];
  isBG <- y[[6]];
  F0 <- y[[7]];  
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  lenF <- varK * (sum(fdim) - length(fdim));
  lenQ <- N * (K - 1);
  
  function(p) {
    return(all(boundaryTurbo_F(p[1:lenF], fdim, varK), boundaryTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, N)));
  }
  
}
