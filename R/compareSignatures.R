#' From the list of Parameters from the getPMSignature function, 
#' return the order of signatures so that signatures from the several lists are well aligned.
#' 
#' @param Params the vector of the result of getPMSignature (EstimatedParameters class)

#' @export
alignmentSignature <- function(Params) {
  
  paramNum <- length(Params);
  if (paramNum <= 1) {
    stop("the length of the input should be at least 2!");
  } 
  
  # get the order of the number of signatures
  sigNums <- c();
  for (l in 1:paramNum) {
    sigNums <- c(sigNums, slot(Params[[l]], "signatureNum") - as.integer(slot(Params[[l]], "isBackGround")));
  }
  calcOrder <- order(sigNums);
  
  
  alignmentList <- list();
  alignmentList[[calcOrder[1]]] <- 1:sigNums[calcOrder[1]];
  for (l in 1:(length(calcOrder) - 1)) {
    
    # get the signature feature dist. for two parameters in sequence
    F_prev <- slot(Params[[calcOrder[l]]], "signatureFeatureDistribution");
    F_next <- slot(Params[[calcOrder[l + 1]]], "signatureFeatureDistribution");
    
    # calculate the distance for eatch signature
    distMat <- matrix(0, sigNums[calcOrder[l]], sigNums[calcOrder[l + 1]]);
    for (k1 in 1:sigNums[calcOrder[l]]) {
      for (k2 in 1:sigNums[calcOrder[l + 1]]) {
        distMat[k1, k2] <- sum((F_prev[k1,,] - F_next[k2,,])^2);
      }
    }
    
    # get the alignment order
    tempOrder <- rep(0, sigNums[calcOrder[l + 1]]);
    for (k1 in 1:sigNums[calcOrder[l]]) {
      sigInd <- which(distMat[k1, ] == min(distMat[k1, ]));
      tempOrder[sigInd] <- alignmentList[[calcOrder[l]]][k1]
      # tempOrder[k] <- c(tempOrder, alignmentList[[calcOrder[l]]][sigInd]);
      distMat[, sigInd] <- Inf;
    }
    tempOrder[tempOrder == 0] <- (sigNums[calcOrder[l]] + 1):sigNums[calcOrder[l + 1]];
    alignmentList[[calcOrder[l + 1]]] <- tempOrder;
    
  }
  
  return(alignmentList);
  
}


getSignaturesForMultipleK <- function(inputFile, numBases, trDir) {
 
  G <- readMPFile(inputFile, numBases = numBases, trDir = trDir);

  BG_prob <- readBGFile(G);
  Param2 <- getPMSignature(G, K = 2, BG = BG_prob);
  Param3 <- getPMSignature(G, K = 3, BG = BG_prob);
  Param4 <- getPMSignature(G, K = 4, BG = BG_prob);
  Param5 <- getPMSignature(G, K = 5, BG = BG_prob);
  Param6 <- getPMSignature(G, K = 6, BG = BG_prob);
  Params <- c(Param2, Param3, Param4, Param5, Param6);


  alignOrder <- alignmentSignature(Params);
  layMat <- matrix(0, 5, 5);
  num <- 1;
  for (l1 in 1:5) {
    for (l2 in 1:l1) {
      layMat[l1, l2] <- num;
      num <- num + 1;
    }
  }
  layout(layMat);

  tempMar <- par("mar");
  par(mar = 0.4 * tempMar);

  for (l1 in 1:5) {
    for (l2 in 1:l1) {
      visPMSignature(Params[[l1]], which(alignOrder[[l1]] == l2));
      # visPMSignature(Params[[l1]], l2);
    }
  }

}


