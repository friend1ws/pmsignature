

#' Read the raw mutation data with the mutation feature vector format.
#' 
#' @param infile the path for the input file for the mutation feature data.
#' @param infile the path for the input file for the mutation data.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is considered or not
#' @param type this argument can take either independent, full, or custom.
#' @export
readMFVFile <- function(infile, numBases = 3, trDir = FALSE, type = "custom") {
  
  if (!(type %in% c("independent", "full", "custom"))) {
    stop('The argument type should be eigher "independent", "full" or "custom"');
  }
  
  if (type == "custom") {
    if (!missing(numBases)) {
      warning('the argument numBaeses is applicable only when type = "custom"');
    }
    if (!missing(trDir)) {
      warning('the argument trDir is applicable only when type = "custom"');
    }    
  }

  mutFile <- read.table(infile, sep="\t", header=FALSE);
  sampleName_str <- as.character(mutFile[,1]);
  mutFeatures <- mutFile[,2:ncol(mutFile)];
  fdim <- unname(apply(mutFeatures, 2, max));
  
  if (type == "independent") {

    if (numBases %% 2 != 1) {
      stop("numBases should be odd numbers");
    }
    
    # check the dimension of the mutation feature vector
    if (ncol(mutFeatures) < numBases + as.integer(trDir)) {
      stop('When type = "independent", the dimension of feature vectors should be at least (numBases + as.integer(trDir))');
    }
    
    # check the 1-st column of the input feature vector (substitution type, should be 1:6)
    if (any(!(mutFeatures[,1] %in% 1:6))) {
      stop(paste('When type = "independent", the 1st feature values should be 1 to 6\n',
                 'In the 1st column, the following rows have values other than 1 to 6;\n',
                 paste0(which(!(mutFeatures[,1] %in% 1:6)), collapse = ","), sep=""));
    }
    
    # check the 2nd to (numBases)-th columns of the input feature vector (flanking bases, should be 1:4)
    for (i in 2:numBases) {
      if (any(!(mutFeatures[,i] %in% 1:4))) {
        stop(paste('When type = "independent", the 2 to (numBases)-th feature values should be 1 to 4\n',
                   'In the ', i, '-th column, the following rows have values other than 1 to 4;\n', 
                   paste0(which(!(mutFeatures[,i] %in% 1:4)), collapse= ","), sep=""));
      }
    }
    
    # check the (numBases + 1)-th column of the input feature vector (the flag for transcription direction, should be 1:2)
    if (trDir == TRUE) {
      if (any(!(mutFeatures[,numBases + 1] %in% 1:2))) {
        stop(paste('When type = "independent" and trDir = TRUE, ', 
                   'the (numBases + 1)-th feature values should be 1 to 2\n',
                   'In the ', (numBases + 1), '-th column, following rows have values other than 1 to 2;\n', 
                   paste0(which(!(mutFeatures[,numBases + 1] %in% 1:2)), collapse = ","), sep=""));
      }
    }

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
  
  mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))));
  rownames(mutFeatList) <- NULL;
  rownames(procCount) <- NULL;

  
  return(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount)
  ))
  
}


#' Read the raw mutation data of mutation position format.
#' 
#' @param infile the path for the input file for the mutation data.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is considered or not
#' @param type this argument can take either independent, full, or custom.
#' 
#' @export
readMPFile <- function(infile, numBases = 3, trDir = FALSE, type = "independent") {

  fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)));
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers");
  }
  centerInd <- (numBases + 1) / 2;

  mutFile <- read.table(infile, sep="\t", header=FALSE);

  vr <- VariantAnnotation::VRanges(mutFile[,2], IRanges::IRanges(mutFile[,3], mutFile[,3]),
               ref = mutFile[,4],
               alt = mutFile[,5],
               sampleNames = mutFile[,1]
                )

  gr <- GenomicRanges::granges(vr) ## drop mcols

  ranges <- GenomicRanges::resize(gr, numBases, fix = "center")
  context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, ranges);

  ref_base <- Biostrings::DNAStringSet(VariantAnnotation::ref(vr));
  alt_base <- Biostrings::DNAStringSet(VariantAnnotation::alt(vr));
  sampleName_str <- as.character(VariantAnnotation::sampleNames(vr));

  removeInd <- which(XVector::subseq(context, start = centerInd, end = centerInd) != ref_base);
  if (sum(removeInd) > 0) {
    warning(paste("The central bases are inconsistent in", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(alt_base);
  removeInd <- which(rowSums(alphabetFreq[,1:4]) != 1);
  if (sum(removeInd) > 0) {
    warning(paste("The characters other than (A, C, G, T) are included in alternate bases of", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(context);
  removeInd <- which(alphabetFreq[,"A"] + alphabetFreq[,"C"] + alphabetFreq[,"G"] + alphabetFreq[,"T"] != numBases);
  if (sum(removeInd) > 0) {
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
    warning(paste("The characters other than (A, C, G, T) are included in flanking bases of", length(removeInd), "mutations. We have removed them."));
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
    vr_strand <- cbind(S4Vectors::mcols(vr_txdb)@listData$QUERYID, as.character(S4Vectors::as.data.frame(strand(vr_txdb))$value));
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
    sampleName_str <- sampleName_str[strandInfo != "*"];
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

  mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(numBases + as.integer(trDir))));
  rownames(mutFeatList) <- NULL;
  rownames(procCount) <- NULL;

  return(new(Class = "MutationFeatureData", 
             type = "independent",
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount)
             )
  )

}



#' Read and format the background vector data
#' 
#' @param bgfile the path for the background mutation signature file
#' @export
readBGFile <- function(bgfile) {
  
  adata <- read.table(bgfile, sep="\t");
  return(adata[,2]);

}
