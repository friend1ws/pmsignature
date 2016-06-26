#' Generate simulation data
#' 
#' @param type this argument can take either independent, full
#' @param numBases the number of upstream and downstream flanking bases
#' @param trDir the index representing whether transcription direction is considered or not
#' @param K the number of mutation signature
#' @param sampleNum the number of cancer genomes in the simulation data
#' @param mutationNum the number of mutations in each cancer genome
#' @param param_alpha the parameter of the Dirichlet distribution for the sample signature distribution
#' @param param_gamma the parameter of the Diriculet distribution for signature feature distribution
#' @param isBG the logical value showing whether the background signature is used or not
#'
#' @export
makeSimData <- function(type = "independent", numBases = 3, trDir = FALSE, fdim = NULL, K = 3, sampleNum = 10, mutationNum = 100, param_alpha = 1, param_gamma = 1, isBG = FALSE) {

  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
  } else if (type == "custom") {
    if (is.null(fdim)) stop('for type "custom"", argument fdim should be specified')
  } else {
    stop('the type argument has to be "independent", "full" or "custom"')
  }
  
  if (type != "custom" & numBases %% 2 != 1) {
    stop("numBases should be odd numbers")
  }
  
  ##########
  # obtaining background signature
  if (type != "custom" & isBG == TRUE) {
    if (numBases > 5 | numBases < 3) {
      stop('Background data whose number of flanking bases is other than 3 or 5 is not available')
    }
    
    if (type == "independent") {
      tempType <- "ind"
    } else if (type == "full") {
      tempType <- "full"
    } else {
      stop('Background data for types other than "independent" or "full" is not available')
    }
  
    if (trDir == TRUE) {
      bgfile <- paste("bgdata/bg.", tempType, numBases, "_dir.txt", sep="")
    } else {
      bgfile <- paste("bgdata/bg.", tempType, numBases, ".txt", sep="")
    }
  
    bdata <- read.table(system.file(bgfile, package = "pmsignature"), header = FALSE, sep="\t")
    bprob <- bdata[,2]
    names(bprob) <- bdata[,1]
    
    varK <- K - 1
    
  } else {
    varK <- K
    BG <- 0
  }
  ##########

  ##########
  # generating the 'true' parameter
  F <- array(0, c(varK, length(fdim), max(fdim)))
  for (k in 1:varK) {
    for (kk in 1:length(fdim)) {
      F[k,kk,1:fdim[kk]] <- rgamma(fdim[kk], param_alpha)
      F[k,kk,1:fdim[kk]] <- F[k,kk,1:fdim[kk]] / sum(F[k,kk,1:fdim[kk]])
    }
  }
  
  Q <- matrix(0, sampleNum, K)
  for (i in 1:sampleNum) {
    Q[i,] <- rgamma(K, param_gamma)
    Q[i,] <- Q[i,] / sum(Q[i,])
  }
  ##########
  
  
  ##########
  # based on the true parameters above, generate mutation feature for each sample and mutation
  currentInd <- 0
  sampleName_str <- rep(0, sampleNum * mutationNum)
  mutFeatures <- matrix(0, sampleNum * mutationNum, length(fdim))
  for (n in 1:sampleNum) {
    
    sampleName_str[currentInd + 1:mutationNum] <- paste("sample_", n, sep="")
    Z <- sample(1:K, mutationNum, replace = TRUE, prob = Q[n,])
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        mutFeatures[currentInd + which(Z == k), kk] <- sample(max(fdim), sum(Z == k), replace = TRUE, prob = F[k, kk, ])
      }
    }
    
    if (isBG == TRUE) {
      tempBG_Str <- names(bprob)[sample(1:length(bprob), sum(Z == K), replace = TRUE, prob = bprob)]
      
      if (length(fdim) == 1) {
        mutFeatures[currentInd + which(Z == K), 1] <- as.integer(tempBG_Str)
      } else {
        mutFeatures[currentInd + which(Z == K), ] <- t(vapply(tempBG_Str, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))))
      } 
      
    }
    
    currentInd <- currentInd + mutationNum
  }
  ##########
  
  
  suSampleStr <- sort(unique(sampleName_str))
  lookupSampleInd <- 1:length(suSampleStr)
  names(lookupSampleInd) <- suSampleStr
  sampleIDs <- lookupSampleInd[sampleName_str]
  
  
  featStr <- apply(mutFeatures, 1, paste0, collapse=",")
  
  suFeatStr <- sort(unique(featStr))
  lookupFeatInd <- 1:length(suFeatStr)
  names(lookupFeatInd) <- suFeatStr
  
  rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[featStr])
  
  tableCount <- table(rawCount)
  w <- which(tableCount > 0, arr.ind=TRUE)
  procCount <- cbind(w[,2], w[,1], tableCount[w])
  
  
  if (type == "full") {
    # mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(1)))
    mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1)
  } else {
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))))
  } 
  
  rownames(mutFeatList) <- NULL
  rownames(procCount) <- NULL
  
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

#' the function for converting the mutation signature matrix to a vector
#' 
#' @param Fmat a matrix for mutation signature
#' @param fdim a vector specifying the number of possible values for each mutation signature
convertSignatureMatrixToVector <- function(Fmat, fdim) {
  
  M <- prod(fdim)
  Fvec <- rep(1, M)
  
  temp1 <- 1
  temp2 <- 1
  for (i in 1:length(fdim)) {
    temp1 <- temp1 * fdim[i]
    divInd <- (1:M - 1) %% temp1 + 1
    for (j in 1:fdim[i]) {
      targetInd <- divInd > temp2 * (j - 1) & divInd <= temp2 * j
      Fvec[targetInd] <- Fvec[targetInd] * Fmat[i,j]
    }
    temp2 <- temp2 * fdim[i]
  }
  
  return(Fvec)
}


#' the function for calculating cosine distances between two mutation signature matrices.
#' 
#' @param F_1 the first matrix for mutation signature
#' @param F_2 the second matrix for mutation signature
#' @param fdim a vector specifying the number of possible values for each mutation signature
#' @export
getCosDistance <- function(F_1, F_2, fdim) {
   
  if (any(dim(F_1) != dim(F_2))) {
    stop("possible features for the two input parameters are different")
  }
  
  K_1 <- dim(F_1)[1]
  K_2 <- dim(F_2)[1]
  M <- prod(fdim)
  
  # I don't like this way of writing... but I have no idea currently... (Y.S. 20150215)
  for (k in 1:K_1) {
    for (i in 1:length(fdim)) {
      if (fdim[i] < dim(F_1)[3]) {
        if (any(F_1[k,i,(fdim[i] + 1):dim(F_1)[3]] != 0)) {
          stop("the first input matrix is not consistent to the possible feature vector")
        }
      }
    }
  }
  for (k in 1:K_2) {
    for (i in 1:length(fdim)) {
      if (fdim[i] < dim(F_2)[3]) {
        if (any(F_2[k,i,(fdim[i] + 1):dim(F_2)[3]] != 0)) {
          stop("the first input matrix is not consistent to the possible feature vector")
        }
      }
    }
  }
  
  F_1_mat <- matrix(0, K_1, M)
  F_2_mat <- matrix(0, K_2, M)
  
  for (k in 1:K_1) {
    inputF_1 <- F_1[k,,]
    dim(inputF_1) <- c(length(fdim), max(fdim))
    F_1_mat[k,] <- convertSignatureMatrixToVector(inputF_1, fdim)
  }
  for (k in 1:K_2) {
    inputF_2 <- F_2[k,,]
    dim(inputF_2) <- c(length(fdim), max(fdim))
    F_2_mat[k,] <- convertSignatureMatrixToVector(inputF_2, fdim)
  }
  
  F_1_mat_nor <- diag(sqrt(rowSums(F_1_mat^2))^(-1)) %*% F_1_mat
  F_2_mat_nor <- diag(sqrt(rowSums(F_2_mat^2))^(-1)) %*% F_2_mat
  
  cos_F1_F2 <- F_1_mat_nor %*% t(F_2_mat_nor)
  
  tempSims <- rep(0, K_1)
  tcos_F1_F2 <- cos_F1_F2
  for (k in 1:K_1) {
    maxV <- max(tcos_F1_F2[k,])
    maxInd <- which(tcos_F1_F2[k,] == maxV)
    tcos_F1_F2[,-maxInd]
    tempSims[k] <- maxV
  }
  
  return(tempSims)
  
}
