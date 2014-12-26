#' Read and format the raw mutation data (deprecated and not used anymore)
#' 
#' @param infile the path for the input file for the raw mutation data
#' @export
readRawMutFile_old <- function(infile) {
  
  adata <- read.table(infile, sep="\t", header=FALSE, fill=TRUE);
  N <- adata[1,1];
  fdim <- as.integer(adata[2,1:adata[1,2]]);
  M <- prod(fdim);
  L <- length(fdim);
  
  sampleIDs <- as.integer(adata[3:nrow(adata),1]);
  mutFeatures <- adata[3:nrow(adata),2:(length(fdim) + 1)];
  
  G <- matrix(0, N, M);
  mutIDs <- rowSums(matrix(c(1, cumprod(fdim)[1:(L-1)]), length(sampleIDs), L, byrow=TRUE) * (mutFeatures - 1)) + 1;
  
  for (i in 1:length(mutIDs)) {
    G[sampleIDs[i], mutIDs[i]] = G[sampleIDs[i], mutIDs[i]] + 1;
  }
  
  return(G);

}
