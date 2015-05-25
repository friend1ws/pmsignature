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


#' Generate the figure of mutation signature for multiple K aligning the signatures.

#' @param inputFile the path for the input file (MP format)
#' @param numBases the number of flanking bases
#' @param trDir the index for the usage of transcription strand bias
#' @param startK the minimum number of signatures
#' @param endK the maximum number of signatures
#' 
#' @export
getSignaturesForMultipleK <- function(inputFile, numBases = 5, trDir = TRUE, startK = NULL, endK = 6) {
 
  G <- readMPFile(inputFile, numBases = numBases, trDir = trDir);

  if (is.null(startK)) {
    if (is.BackGround == TRUE) {
      startK <- 2;
    } else {
      startK <- 1;
    }
  }
  
  Params <- c();
  if (isBackGround == TRUE) {
    BG_prob <- readBGFile(G);
    Params <- c();
    for (k in startK:endK) {
      Params <- c(Params, getPMSignature(G, K = k, BG = BG_prob)); 
    }
  } else {
    
    Params <- c();
    for (k in startK:endK) {
      Params <- c(Params, getPMSignature(G, K = k)); 
    }  
    
  }
  

  
  alignOrder <- alignmentSignature(Params);
  layMat <- matrix(0, length(Params), length(Params));
  num <- 1;
  for (l1 in 1:length(Params)) {
    for (l2 in 1:l1) {
      layMat[l1, l2] <- num;
      num <- num + 1;
    }
  }
  layout(layMat);

  tempMar <- par("mar");
  par(mar = 0.4 * tempMar);

  for (l1 in 1:length(Params)) {
    # aInds <- c()
    for (l2 in 1:l1) {
      aInd <- which(alignOrder[[l1]] == l2);
      # aInds <- c(aInds, aInd);
      visPMSignature(Params[[l1]], aInd);
      # visPMSignature(Params[[l1]], l2);
    }
    # visMembership(G, Params[[l1]], toSample = 100, multiplySampleNum = FALSE, colourBrewer = "Set2", reorder = 1:6)
    
  }

}


