#' Read and format the raw mutation data
#' 
#' @param infile the path for the input file for the raw mutation data
#' @export
readRawMutFile <- function(infile) {
  
  adata <- read.table(infile, sep="\t", header=FALSE, fill=TRUE);
  N <- adata[1,1];
  fdim <- as.integer(adata[2,1:adata[1,2]]);
  M <- prod(fdim);
  L <- length(fdim);
  
  sampleIDs <- as.integer(adata[3:nrow(adata),1]);
  mutFeatures <- adata[3:nrow(adata),2:(length(fdim) + 1)];
  
  # I don't like this implementation..
  # Ideally I would like to remove the step calculating the value for each mutation feature,
  # and directly order the mutation features like
  # order(uniqueMutFeatures[,6], uniqueMutFeatures[,5], uniqueMutFeatures[,4], uniqueMutFeatures[,3], uniqueMutFeatures[,2], uniqueMutFeatures[,1])
  # But, currently, I don't know how to deal the unspecified dimension of mutation feature..
  G <- matrix(0, N, M);
  mutIDs <- rowSums(matrix(c(1, cumprod(fdim)[1:(L-1)]), length(sampleIDs), L, byrow=TRUE) * (mutFeatures - 1)) + 1;
  
  for (i in 1:length(mutIDs)) {
    G[sampleIDs[i], mutIDs[i]] = G[sampleIDs[i], mutIDs[i]] + 1;
  }
  
  return(G);

}


#' Read and format the raw mutation data
#' 
#' @param infile the path for the input file for the raw mutation data
#' @export
readRawMutFile_sparse <- function(infile) {
  
  adata <- read.table(infile, sep="\t", header=FALSE, fill=TRUE);
  N <- adata[1,1];
  fdim <- as.integer(adata[2,1:adata[1,2]]);
  M <- prod(fdim);
  L <- length(fdim);

  sampleIDs <- as.integer(adata[3:nrow(adata),1]);
  mutFeatures <- adata[3:nrow(adata),2:(length(fdim) + 1)];

  featStr <- apply(mutFeatures, 1, paste, collapse=",");
  suFeatStr <- sort(unique(featStr));
  lookupFeatInd <- 1:length(suFeatStr);
  names(lookupFeatInd) <- suFeatStr;

  rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[apply(mutFeatures, 1, paste, collapse=",")]);

  tableCount <- table(rawCount);
  w <- which(tableCount > 0, arr.ind=TRUE);
  procCount <- cbind(w[,2], w[,1], tableCount[w]);

  mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(6)));
  rownames(mutFeatList) <- NULL;

  return(list(N, fdim, t(mutFeatList), t(procCount)));

}



#' Read and format the background vector data
#' 
#' @param bgfile the path for the background mutation signature file
#' @export
readBGFile <- function(bgfile) {
  
  adata <- read.table(bgfile, sep="\t");
  return(adata[,2]);

}
