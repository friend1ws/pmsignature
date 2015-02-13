
makeSimData <- function(type = "independent", numBases = 3, trDir = FALSE, K = 3, sampleNum = 10, mutationNum = 100, param_alpha, param_gamma, isBG = FALSE) {

  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)));
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)));
  } else {
    stop('for reading mutation position format, the type argument has to be "independent" or "full"');
  }
  
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers");
  }
  
  if (isBG == TRUE) {
    if (numBases > 5 | numBases < 3) {
      stop('Background data whose number of flanking bases is other than 3 or 5 is not available');
    }
    
    if (type == "independent") {
      tempType <- "ind";
    } else if (type == "full") {
      tempType <- "full";
    } else {
      stop('Background data for types other than "independent" or "full" is not available');
    }
  
    if (trDir == TRUE) {
      bgfile <- paste("bgdata/bg.", tempType, numBases, "_dir.txt", sep="");
    } else {
      bgfile <- paste("bgdata/bg.", tempType, numBases, ".txt", sep="");
    }
  
    bdata <- read.table(system.file(bgfile, package = "pmsignature"), header = FALSE, sep="\t");
    bprob <- bdata[,2];
    names(bprob) <- bdata[,1];
    
    varK <- K - 1;
    
  } else {
    varK <- K;
    BG <- 0;
  }
  

  
  F <- array(0, c(varK, length(fdim), max(fdim)));
  for (k in 1:varK) {
    for (kk in 1:length(fdim)) {
      F[k,kk,1:fdim[kk]] <- rgamma(fdim[kk], param_alpha);
      F[k,kk,1:fdim[kk]] <- F[k,kk,1:fdim[kk]] / sum(F[k,kk,1:fdim[kk]]);
    }
  }
  
  Q <- matrix(0, N, K); 
  for (i in 1:N) {
    Q[i,] <- rgamma(K, param_gamma);
    Q[i,] <- Q[i,] / sum(Q[i,]);
  }
  
  currentInd <- 0;
  sampleName_str <- rep(0, sampleNum * mutationNum);
  mutFeatures <- matrix(0, sampleNum * mutationNum, length(fdim));
  for (n in 1:sampleNum) {
    
    sampleName_str[currentInd + 1:mutationNum] <- paste("sample_", n, sep="");
    Z <- sample(1:K, mutationNum, replace = TRUE, prob = Q[n,]);
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        mutFeatures[currentInd + which(Z == k), kk] <- sample(max(fdim), sum(Z == k), replace = TRUE, prob = F[k, kk, ]);
      }
    }
    
    if (isBG == TRUE) {
      tempBG_Str <- names(bprob)[sample(1:length(bprob), sum(Z == K), replace = TRUE, prob = bprob)];
      
      if (length(fdim) == 1) {
        mutFeatures[currentInd + which(Z == K), 1] <- as.integer(tempBG_Str);
      } else {
        mutFeatures[currentInd + which(Z == K), ] <- t(vapply(tempBG_Str, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))));
      } 
      
    }
    
    currentInd <- currentInd + mutationNum;
  }
  
  
  suSampleStr <- sort(unique(sampleName_str));
  lookupSampleInd <- 1:length(suSampleStr);
  names(lookupSampleInd) <- suSampleStr;
  sampleIDs <- lookupSampleInd[sampleName_str];
  
  
  featStr <- apply(mutFeatures, 1, paste0, collapse=",");
  
  suFeatStr <- sort(unique(featStr));
  lookupFeatInd <- 1:length(suFeatStr);
  names(lookupFeatInd) <- suFeatStr;
  
  rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[featStr]);
  
  tableCount <- table(rawCount);
  w <- which(tableCount > 0, arr.ind=TRUE);
  procCount <- cbind(w[,2], w[,1], tableCount[w]);
  
  
  if (length(fdim) == 1) {
    # mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(1)));
    mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1);
  } else {
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(numBases + as.integer(trDir))));
  } 
  
  rownames(mutFeatList) <- NULL;
  rownames(procCount) <- NULL;
  
  return(list(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = data.frame()
  ), F, Q)
  )
  
  
}