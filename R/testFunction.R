
crossPMSignature <- function(G, K, fdim, isBG, BG0, F0, Q0, crossNum) {
  
  if (isBG) {
    varK <- K - 1;
  } else {
    varK <- K;
  }
  
  M <- prod(fdim);
  N <- nrow(G);
  
  Ls <- rep(0, crossNum);
  for (bbb in 1:crossNum) {
    
    N_train <- N - floor(N / 5);
    N_test <- floor(N / 5);
    testInd <- sort(order(runif(N, 0, 1))[1:N_test]);
    testG <- G[testInd, ];
    trainG <- G[-testInd,];
    
    
    p_train <- c(convertToTurbo_F(as.vector(F0), fdim, K, isBG), convertToTurbo_Q(as.vector(Q0[-testInd,]), K, N_train));
    Y_train <- list(trainG, fdim, K, N_train, M, isBG, BG0);   
    res_train <- turboem(par=p_train, y=Y_train, fixptfn=updatePMSParam, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y_train), control.run = list(convtype = "objfn", tol = 1e-4));
    
    
    lenF <- varK * (sum(fdim) - length(fdim));
    lenQ <- N_train * (K - 1);
    F <- convertFromTurbo_F(res_train$par[1:lenF], fdim, K, isBG);
    Q <- convertFromTurbo_Q(res_train$par[(lenF + 1):(lenF + lenQ)], K, N_train);
    # Q <- rep(1 / K, N_test * K);
    dim(F) <- c(varK, length(fdim), max(fdim));
    dim(Q) <- c(N_train, K);
    # dim(Q) <- c(N_test, K);
    
    # p_test <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, N_test));
    # Y_test <- list(testG, fdim, K, N_test, M, isBG, BG0);   
    # res_test <- turboem(par=p_test, y=Y_test, fixptfn=updatePMSParam_nonF, objfn=calcPMSLikelihood, method=c("squarem"), pconstr=PMSboundary(Y_test), control.run = list(convtype = "objfn", tol = 1e-4));
    
    # lenF <- varK * (sum(fdim) - length(fdim));
    # F_test <- convertFromTurbo_F(res_test$par[1:lenF], fdim, K, isBG);
    # dim(F_test) <- c(varK, length(fdim), max(fdim));
    tempL <- 0;
    for (nnn in 1:N_train) {
      Q_emp <- matrix(Q[nnn,], N_test, K, byrow=T);
      tempL <- tempL + getLogLikelihoodC(as.vector(testG), as.vector(F), as.vector(Q_emp), fdim, K, N_test, M, isBG, BG0) / N_train;
    }
    
    Ls[bbb] <- tempL / sum(testG);
    
    print(c(bbb, res_train$itr, res_train$runtime[3], res_test$itr, res_test$runtime[3], Ls[bbb]));
    
    
  }
  
  return(Ls);
  
}


updatePMSParam_nonF <- function(p, y) {
  
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
  # F <- updateMstepFC(as.vector(Theta), as.vector(G), fdim, K, N, M, isBG);
  Q <- updateMstepQC(as.vector(Theta), as.vector(G), K, N, M);
  #########################################
  
  return(c(convertToTurbo_F(as.vector(F), fdim, K, isBG), convertToTurbo_Q(as.vector(Q), K, N)));
  
}


visPMS_ind5_mod <- function(V1 = matrix(1 / 4, 2, 4), V2 = matrix(1 / 4, 4, 4)) {
  
  A <- matrix(0, 5, 4);
  B <- matrix(0, 4, 4);
  
  A[1,] <- V2[1,];
  A[2,] <- V2[2,];
  A[4,] <- V2[3,];
  A[5,] <- V2[4,];
  A[3,2] <- V1[1,1];
  A[3,4] <- V1[1,2];
  
  B[2,] <- V1[2,];
  B[4,] <- V1[2,];
  
  num2base <- c("A", "C", "G", "T");
  
  plot.window(xlim=c(-0.25, 6.5), ylim=c(-0.25, 3.25));
  
  tcols <- c(4, 7, 3, 2);
  
  startx <- 0;
  for(l in 1:5) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col=tcols[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=1.2)
      }
      startx <- endx;
    }
    startx <- startx + 0.25;
  }
  
  startx <- 2 * 1.25;
  for (w in 1:4) {
    starty <- 2;
    endx <- startx + A[3,w];
    for(ww in 1:4) {
      endy <- starty + B[w,ww];
      polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=tcols[ww], border=F);
      if ((endy - starty > 1 / 4) & (endx - startx > 1 / 4)) {
        text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col="white", cex=1.2)
      }
      starty <- endy;
    }
    startx <- endx
    starty <- endy;
  }
  
  ##########
  # draw arrow
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + 2.5;
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1;
  
  polygon(xs, ys, col=8, border=F);
  ##########
  
}

