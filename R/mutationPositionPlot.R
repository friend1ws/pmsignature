getMembershipParam <- function(mutationFeatureData, Param, BG = NULL) {
  
  K <- slot(Param, "signatureNum");
  isBG <- slot(Param, "isBackGround");
  
  if (isBG == TRUE) {
    if (!is.null(BG)) {
      varK <- K - 1;
    } else {
      stop(paste("The input parameter is estimated using a background signature.\n",
                 "Please specify the same background signature."));
    }
  } else {
    if (!is.null(BG)) {      
      warning(paste("The input parameter is estimated without using a background signature.\n",
                    "Specified background signature is ignored."));
    }
    varK <- K;
    BG <- 0;
  }
  
  sampleNum <- length(slot(mutationFeatureData, "sampleList"));
  fdim <- slot(mutationFeatureData, "possibleFeatures");
  
  F0 <- slot(Param, "signatureFeatureDistribution");
  Q0 <- slot(Param, "sampleSignatureDistribution");
  
  patternList <- slot(mutationFeatureData, "featureVectorList");
  sparseCount <- slot(mutationFeatureData, "countData");
  mutationPosition <- slot(mutationFeatureData, "mutationPosition");

  patternNum <- ncol(patternList);
  samplePatternNum <- ncol(sparseCount);

  dim(Q0) <- c(sampleNum, K);
  Q0 <- t(Q0);
  ####################
  # E-step
  Theta <- updateTheta_NormalizedC(as.vector(patternList), as.vector(sparseCount), as.vector(F0), as.vector(Q0), fdim, K, sampleNum, patternNum, samplePatternNum, isBG, BG);
  dim(Theta) <- c(K, samplePatternNum);
  
  tTheta <- t(Theta);
  rownames(tTheta) <- paste(sparseCount[2,], sparseCount[1,], sep=",");
  tTheta_maxSig <- apply(tTheta, 1, function(x) {which(x == max(x))});  
  tTheta_entropy <- apply(tTheta, 1, function(x) {return( - sum(x * log2(x)))});

  
  membershipParam <- tTheta[paste(mutationPosition[,3], mutationPosition[,4], sep=","), ];
  colnames(membershipParam) <- paste("sig", 1:K, sep="_");
  rownames(membershipParam) <- NULL;

  maxSig <- tTheta_maxSig[paste(mutationPosition[,3], mutationPosition[,4], sep=",")];
  entropy <- tTheta_entropy[paste(mutationPosition[,3], mutationPosition[,4], sep=",")];
  names(maxSig) <- NULL;
  names(entropy) <- NULL;
  
  # obtaining the coordinate for the plot
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
  startCoordinate <- cumsum(seqlengths(genome) / 100) - seqlengths(genome) / 100;
  plotCoordinate <- startCoordinate[mutationPosition[,1]] + mutationPosition[,2] / 100;
  names(plotCoordinate) <- NULL;
  

  return(data.frame(sampleID = mutationPosition[,3],
                    chr = mutationPosition[,1],
                    pos = mutationPosition[,2],
                    plotCoordinate = plotCoordinate,
                    membershipParam,
                    maxSig = maxSig,
                    entropy = entropy,
                    stringsAsFactors = FALSE));
        
}


# ggplot(a, aes(x = plotCoordinate, y = sampleID, color = factor(maxSig), alpha = entropy)) + geom_point(position = position_jitter(width=0, height=0.3), size = rel(1)) + scale_alpha(range = c(0, 1), limits = c(0, 1.6), trans = "reverse");


