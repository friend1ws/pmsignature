

#' Read the raw mutation data with the mutation feature format.
#' 
#' @param infile the path for the input file for the mutation feature data.
#' @export
readRawMutfeatFile <- function(infile) {
  
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


#' Read the raw mutation data.
#' 
#' @param infile the path for the input file for the mutation data.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' 
#' VariantAnnotation::VRanges
#' @export
readMutFile <- function(infile, numBases = 3, trDir = FALSE) {

  # library(VariantAnnotation);
  # library(BSgenome.Hsapiens.UCSC.hg19);

  fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)));
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers");
  }
  centerInd <- (numBases + 1) / 2;

  mutFile <- read.table(infile, sep="\t", header=FALSE);

  vr = VariantAnnotation::VRanges(mutFile[,2], IRanges::IRanges(mutFile[,3], mutFile[,3]),
               ref = mutFile[,4],
               alt = mutFile[,5],
               sampleNames = mutFile[,1]
                )

  gr = GenomicRanges::granges(vr) ## drop mcols

  ranges = GenomicRanges::resize(gr, numBases, fix = "center")
  context = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, ranges);

  ref_base = Biostrings::DNAStringSet(VariantAnnotation::ref(vr));
  alt_base = Biostrings::DNAStringSet(VariantAnnotation::alt(vr));


  removeInd <- which(XVector::subseq(context, start = centerInd, end = centerInd) != ref_base);
  if (sum(removeInd) > 0) {
    warning(paste("The central bases are inconsistent in", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(alt_base);
  removeInd <- which(rowSums(alphabetFreq[,1:4]) != 1);
  if (sum(removeInd) > 0) {
    warning(paste("The characters other than (A, C, G, T) are included in alternate bases of", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(context);
  removedInd <- which(rowSums(alphabetFreq[,1:4]) != numBases);
  if (sum(removeInd) > 0) {
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    warning(paste("The characters other than (A, C, G, T) are included in flanking bases of", sum(removeInd), "mutations. We have removed them."));
  }

  revCompInd <- which(as.character(XVector::subseq(context, start = centerInd, end = centerInd)) %in% c("A", "G"));
  context[revCompInd] <- Biostrings::reverseComplement(context[revCompInd]);
  ref_base[revCompInd] <- Biostrings::reverseComplement(ref_base[revCompInd]);
  alt_base[revCompInd] <- Biostrings::reverseComplement(alt_base[revCompInd]);

  # Obtaining transcription strand information using VariantAnnotation packages
  # I feel it's a bit slow to get the strand information, 
  # and would like to improve the implementation in near future (Y.S., 20141223)
  if (trDir == TRUE) {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene;
    vr_txdb <- VariantAnnotation::locateVariants(vr, txdb, VariantAnnotation::AllVariants(), ignore.strand=TRUE);
    vr_strand <- cbind(vr_txdb@elementMetadata@listData$QUERYID, as.vector(vr_txdb@strand));
    uvr_strand <- unique(vr_strand[vr_strand[, 2] != "*" ,], MARGIN=1);
    rmdup_uvr_strand <- uvr_strand[!duplicated(uvr_strand[, 1]), ];
    txdb_plus_vr_ind <- as.integer(rmdup_uvr_strand[rmdup_uvr_strand[, 2] == "+", 1]);
    txdb_minus_vr_ind <- as.integer(rmdup_uvr_strand[rmdup_uvr_strand[, 2] == "-", 1]);
    
    strandInfo <- rep("*", length(vr));
    strandInfo[setdiff(txdb_plus_vr_ind, revCompInd)] <- "+";
    strandInfo[intersect(txdb_plus_vr_ind, revCompInd)] <- "-";    
    strandInfo[setdiff(txdb_minus_vr_ind, revCompInd)] <- "-";
    strandInfo[intersect(txdb_minus_vr_ind, revCompInd)] <- "+";   
    
    warning(paste("Out of", length(context), "mutations, we could obtain transcription direction information for", 
                  length(txdb_plus_vr_ind) + length(txdb_minus_vr_ind), "mutation. Other mutations are removed."));
    context <- context[strandInfo != "*"];
    ref_base <- ref_base[strandInfo != "*"];
    alt_base <- alt_base[strandInfo != "*"];
    strandInfo <- strandInfo[strandInfo != "*"];
    
  }
  
  
  mutFeatures <- matrix(0, length(ref_base), length(fdim));

  mutFeatures[which(ref_base == "C" & alt_base == "A"), 1] <- 1;
  mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] <- 2;
  mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] <- 3;
  mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] <- 4;
  mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] <- 5;
  mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] <- 6;

  columnInd <- 2;
  for (baseInd in 1:numBases) {
    if (baseInd == centerInd) {
      next;
    }
    mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "A"), columnInd] <- 1;
    mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), columnInd] <- 2;
    mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), columnInd] <- 3;
    mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), columnInd] <- 4;
    columnInd <- columnInd + 1;
  }

  if (trDir == TRUE) {
    mutFeatures[which(strandInfo == "+"), length(fdim)] <- 1;
    mutFeatures[which(strandInfo == "-"), length(fdim)] <- 2;    
  }

  suSampleStr <- sort(unique(VariantAnnotation::sampleNames(vr)));
  lookupSampleInd <- 1:length(suSampleStr);
  sampleIDs = lookupSampleInd[VariantAnnotation::sampleNames(vr)];


  featStr <- apply(mutFeatures, 1, paste0, collapse=",");

  suFeatStr <- sort(unique(featStr));
  lookupFeatInd <- 1:length(suFeatStr);
  names(lookupFeatInd) <- suFeatStr;

  rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[featStr]);

  tableCount <- table(rawCount);
  w <- which(tableCount > 0, arr.ind=TRUE);
  procCount <- cbind(w[,2], w[,1], tableCount[w]);

  mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(5)));
  rownames(mutFeatList) <- NULL;

  return(list(length(suSampleStr), fdim, t(mutFeatList), t(procCount)));
  
}



#' Read and format the background vector data
#' 
#' @param bgfile the path for the background mutation signature file
#' @export
readBGFile <- function(bgfile) {
  
  adata <- read.table(bgfile, sep="\t");
  return(adata[,2]);

}
