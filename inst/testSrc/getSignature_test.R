#' A function for estimating parameters using Squared EM algorithm
#' 
#' @param p this variable includes the parameters for mutation signatures and membership parameters
#' @param y this variable includes the information on the mutation features, 
#' the number of mutation signatures specified and so on
#' @param tol tolerance for the estimation
#' (when the difference of log-likelihoods become below this value, stop the estimation)
#' @param maxIter the maximum number of iteration of estimation
#' @export
mySquareEM_test <- function(p, y, tol = 1e-4, maxIter = 10000) {
  
  prevL <- -Inf;
  step.min <- 1;
  step.max <- 1;
  step.max0 <- 1;
  mstep <- 4;
  objfn.inc <- 1;
  updEvalNum <- 0;
  LEvalNum <- 0;
  useSquareEM <- 0;
  iterNum <- 0;
  convFlag <- FALSE;
  startTime <- proc.time();
  
  upd <- updatePMSParam_test(p, y);
  p <- upd[[1]];
  newL <- upd[[2]];
  LEvalNum <- LEvalNum + 1;
  
  for (iterNum in 1:maxIter) {
    
    upd1 <- updatePMSParam_test(p, y);
    p1 <- upd1[[1]];
    L1 <- upd1[[2]];
    updEvalNum <- updEvalNum + 1;
    if ( any(is.nan(unlist(p1))) ) {
      stop("Error in function evaluation");
    }
    
    q1 <- p1 - p;
    sr2 <- crossprod(q1);
    
    upd2 <- updatePMSParam_test(p1, y);
    p2 <- upd2[[1]];
    L2 <- upd2[[2]];
    updEvalNum <- updEvalNum + 1;
    if ( any(is.nan(unlist(p2))) ) {
      stop("Error in function evaluation");
    }
    
    q2 <- p2 - p1;
    sq2 <- sqrt(crossprod(q2));
    sv2 <- crossprod(q2 - q1);
    srv <- crossprod(q1, q2 - q1);
    
    # alpha <- switch(ctrl$version, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
    alpha <- -srv/sv2;
    alpha <- max(step.min, min(step.max, alpha));
    p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1);
    
    newFlag <- 0;
    # This step is done in the original turboEM code...
    # but I cannot understand why this step is necessary....
    if (isTRUE(abs(alpha - 1) > 0.01) ) {
      upd3 <- updatePMSParam_test(p.new, y);
      p.new <- upd3[[1]];
      newL <- upd3[[2]];  
      updEvalNum <- updEvalNum + 1;
      newFlag <- 1;
    } else {
      p.new <- p2;
      newL <- L2;
    }
    
    # when p.new has some problems...
    if (any(is.nan(p.new)) | !PMSboundary(y)(p.new) ) {
      
      p.new <- p2;
      newL <- L2;
      
      # since there was a problem, consider to reduce the amount of step max
      if (isTRUE(all.equal(alpha, step.max))) {
        step.max <- max(step.max0, step.max / mstep);
      }
      alpha <- 1
      
      # when p.new is O.K....
    } else if (newFlag == 1) {
      
      # upd4 <- updatePMSParam_test(p.new, y);
      # p.new <- upd4[[1]];
      # newL <- upd4[[2]];     
      
      # when the calculated log-likelihood has some problems
      # or the difference betwen the calculated log-likelihood is large...
      if (is.nan(newL) | (newL > prevL + objfn.inc)) {
        
        p.new <- p2
        newL <- L2;
        
        # since there was a problem, consider to reduce the amount of step max
        if (alpha == step.max) {
          step.max <- max(step.max0, step.max / mstep);
        }
        alpha <- 1
        
      } else {
        useSquareEM <- useSquareEM + 1;
      }
      
    }
    
    if (isTRUE(all.equal(alpha, step.max))) {
      step.max <- mstep * step.max;
    }
    
    if (step.min < 0 & isTRUE(all.equal(alpha, step.min))) {
      step.min <- mstep * step.min;
    }
    
    p <- p.new;
    
    # for debugging
    # print(c(updEvalNum, LEvalNum, useSquareEM, step.min, step.max, newL));
    
    if (abs(prevL - newL) < tol) {
      convFlag <- TRUE;
      break;
    }
    
    if (!is.nan(newL)) {
      prevL <- newL;
    }
    
  }
  
  calcTime <- proc.time() - startTime;
  
  return(list(par = p,
              value.objfn = newL,
              itr = iterNum,
              fpeval = updEvalNum,
              convergence = convFlag,
              elapsed.time = calcTime[3]));
  
}


#' A function for updating parameters using EM-algorithm
#' 
#' @param p this variable includes the parameters for mutation signatures and membership parameters
#' @param y this variable includes the information on the mutation features, 
#' the number of mutation signatures specified and so on
#' @export
updatePMSParam_test <- function(p, y) {
  
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
  Theta_L <- updateTheta_NormalizedC_new(as.vector(patternList), as.vector(sparseCount), as.vector(F), as.vector(Q), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, F0);
  Theta <- Theta_L[-length(Theta_L)];
  L <- Theta_L[length(Theta_L)];
  
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
  
  return(list(c(convertToTurbo_F(as.vector(F), fdim, K, isBG), 
                convertToTurbo_Q(as.vector(Q), K, sampleNum)
  ), 
  L));
  
}





