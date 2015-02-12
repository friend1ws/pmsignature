

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
      warning('the argument numBases is applicable only when type = "independent", or "full"');
    }
    if (!missing(trDir)) {
      warning('the argument trDir is applicable only when type = "independent", or "full"');
    }    
  }

  mutFile <- read.table(infile, sep="\t", header=FALSE);
  sampleName_str <- as.character(mutFile[,1]);
  mutFeatures <- mutFile[,2:ncol(mutFile), drop = FALSE];
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
  
  if (length(fdim) == 1) {
    mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1);
  } else {
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))));
  }
  
  rownames(mutFeatList) <- NULL;
  rownames(procCount) <- NULL;

  
  return(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = NULL
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
  centerInd <- (numBases + 1) / 2;

  mutFile <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors = FALSE);

  chrInfo <- mutFile[,2];
  posInfo <- mutFile[,3];
  ref_base <- Biostrings::DNAStringSet(mutFile[,4]);
  alt_base <- Biostrings::DNAStringSet(mutFile[,5]);
  sampleName_str <- as.character(mutFile[,1]);
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo, 
                                            start = posInfo, 
                                            end = posInfo), ignore.strand = TRUE);

  ranges <- GenomicRanges::resize(gr, numBases, fix = "center")
  context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, ranges);


  removeInd <- which(XVector::subseq(context, start = centerInd, end = centerInd) != ref_base);
  if (length(removeInd) > 0) {
    warning(paste("The central bases are inconsistent in", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
    chrInfo <- chrInfo[-removeInd];
    posInfo <- posInfo[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(alt_base);
  removeInd <- which(rowSums(alphabetFreq[,1:4]) != 1);
  if (length(removeInd) > 0) {
    warning(paste("The characters other than (A, C, G, T) are included in alternate bases of", length(removeInd), "mutations. We have removed them."));
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
    chrInfo <- chrInfo[-removeInd];
    posInfo <- posInfo[-removeInd];
  }

  alphabetFreq <- Biostrings::alphabetFrequency(context);
  removeInd <- which(alphabetFreq[,"A"] + alphabetFreq[,"C"] + alphabetFreq[,"G"] + alphabetFreq[,"T"] != numBases);
  if (length(removeInd) > 0) {
    context <- context[-removeInd];
    ref_base <- ref_base[-removeInd];
    alt_base <- alt_base[-removeInd];
    sampleName_str <- sampleName_str[-removeInd];
    chrInfo <- chrInfo[-removeInd];
    posInfo <- posInfo[-removeInd];
    warning(paste("The characters other than (A, C, G, T) are included in flanking bases of", length(removeInd), "mutations. We have removed them."));
  }

  revCompInd <- which(as.character(XVector::subseq(context, start = centerInd, end = centerInd)) %in% c("A", "G"));
  context[revCompInd] <- Biostrings::reverseComplement(context[revCompInd]);
  ref_base[revCompInd] <- Biostrings::reverseComplement(ref_base[revCompInd]);
  alt_base[revCompInd] <- Biostrings::reverseComplement(alt_base[revCompInd]);

  # Obtaining transcription strand information using VariantAnnotation packages
  if (trDir == TRUE) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo, 
                                                             start = posInfo, 
                                                             end = posInfo), ignore.strand = TRUE);
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene;
    # exons_txdb <- GenomicFeatures::exons(txdb);
    # gr_txdb <- GenomicRanges::findOverlaps(gr, exons_txdb, ignore.strand = FALSE);
    # gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(exons_txdb[gr_txdb@subjectHits]))));
    txdb_bed <- GenomicFeatures::asBED(txdb);
    gr_txdb <- GenomicRanges::findOverlaps(gr, txdb_bed, ignore.strand = FALSE)
    gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(txdb_bed[gr_txdb@subjectHits]))));
    ugr_strand <- unique(gr_strand[gr_strand[, 2] != "*" ,], MARGIN=1);
    
    rmdup_ugr_strand <- ugr_strand[!duplicated(ugr_strand[, 1]), ];
    txdb_plus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "+", 1]);
    txdb_minus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "-", 1]);
    
    strandInfo <- rep("*", length(gr));
    strandInfo[setdiff(txdb_plus_gr_ind, revCompInd)] <- "+";
    strandInfo[intersect(txdb_plus_gr_ind, revCompInd)] <- "-";    
    strandInfo[setdiff(txdb_minus_gr_ind, revCompInd)] <- "-";
    strandInfo[intersect(txdb_minus_gr_ind, revCompInd)] <- "+";   
    
    warning(paste("Out of", length(context), "mutations, we could obtain transcription direction information for", 
                  length(txdb_plus_gr_ind) + length(txdb_minus_gr_ind), "mutation. Other mutations are removed."));
    
    # this may be unstable... need to review thoroughly later...
    context <- context[strandInfo != "*"];
    ref_base <- ref_base[strandInfo != "*"];
    alt_base <- alt_base[strandInfo != "*"];
    sampleName_str <- sampleName_str[strandInfo != "*"];
    chrInfo <- chrInfo[strandInfo != "*"];
    posInfo <- posInfo[strandInfo != "*"];
    strandInfo <- strandInfo[strandInfo != "*"];
  }
  
  
  if (trDir == FALSE) {
    strandInfo <- NULL;
  }
  
  
  mutFeatures <- getMutationFeatureVector(context, ref_base, alt_base, strandInfo, numBases, type);

    
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

  return(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = data.frame(chr = chrInfo, pos = posInfo, sampleID = unname(sampleIDs), mutID = unname(lookupFeatInd[featStr]), stringsAsFactors = FALSE)
             )
  )

}



#' get mutation feature vector from context sequence data and reference and alternate allele information
#'
#' @param context the context sequence data around the mutated position. This shoud be Biostrings::DNAStringSet class
#' @param ref_base the reference bases at the mutated position.
#' @param alt_base the alternate bases at the mutated position.
#' @param strandInfo transcribed strand information at the mutated position. (this is optional)
#' @param numBases the number of flanking bases around the mutated position.
#' @param type the type of mutation feature vecotr (should be "independent" or "full").
getMutationFeatureVector <- function(context, ref_base, alt_base, strandInfo = NULL, numBases, type) {
  
  trDir <- !is.null(strandInfo);
  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)));
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)));
  } else {
    stop('the type argument has to be "independent" or "full"');
  }
  
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers");
  }
  centerInd <- (numBases + 1) / 2;
  
  
  ##########
  # whether we should add additional check function for the bases of context sequences, reference and alternate allele...?
  ##########
  
  if (type == "independent") {
    
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
    
    if (trDir ==TRUE) {
      mutFeatures[which(strandInfo == "+"), length(fdim)] <- 1;
      mutFeatures[which(strandInfo == "-"), length(fdim)] <- 2;    
    }
    
  } else {
    
    mutFeatures <- matrix(1, length(ref_base), length(fdim));
    
    tempDigits <- 1;
    for (i in 1:((numBases - 1) / 2)) {
      
      baseInd <- numBases + 1 - i;
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] + tempDigits * 1;
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] + tempDigits * 2;
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] + tempDigits * 3;
      
      tempDigits <- tempDigits * 4;
      baseInd <- i;      
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] + tempDigits * 1;
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] + tempDigits * 2;
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] + tempDigits * 3;
      
      tempDigits <- tempDigits * 4;
    }
    
    mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] <- mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] + tempDigits * 1;
    mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] <- mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] + tempDigits * 2;
    mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] + tempDigits * 3;
    mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] + tempDigits * 4;
    mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] + tempDigits * 5;
    
    if (trDir ==TRUE) {
      tempDigits <- tempDigits * 6;
      mutFeatures[which(strandInfo == "-"), 1] <- mutFeatures[which(strandInfo == "-"), 1] + tempDigits * 1;
    }
    
  }
  
  return(mutFeatures);
  
}

#' Read and format the background vector data
#' 
#' @param mutationFeatureData the mutation data processed in the function (readMutFile or readRawMutfeatFile)
#' @export
readBGFile <- function(mutationFeatureData) {
  
  if (slot(mutationFeatureData, "type") == "independent") {
    tempType <- "ind";
  } else if (slot(mutationFeatureData, "type") == "full") {
    tempType <- "full";
  } else {
    stop('Background data for types other than "independent" or "full" is not available');
  }
  
  if (slot(mutationFeatureData, "flankingBasesNum") %in% c(3, 5)) {
    tempNumBase <- slot(mutationFeatureData, "flankingBasesNum");
  } else {
    stop('Background data whose number of flanking bases is other than 3 or 5 is not available');
  }
  
  if (slot(mutationFeatureData, "transcriptionDirection") == TRUE) {
    bgfile <- paste("bgdata/bg.", tempType, tempNumBase, "_dir.txt", sep="");
  } else {
    bgfile <- paste("bgdata/bg.", tempType, tempNumBase, ".txt", sep="");
  }
  
  bdata <- read.table(system.file(bgfile, package = "pmsignature"), sep="\t");
  
  tempFeatureVectorList <- apply(slot(mutationFeatureData, "featureVectorList"), 2, paste0, collapse=",");
  bprob <- bdata[,2];
  names(bprob) <- bdata[,1];
  
  if (!all(tempFeatureVectorList %in% names(bprob))) {
    noNameInd <- which(! (tempFeatureVectorList %in% names(bprob)));
    stop(paste('The information of following mutation features are not included in the specified background file:\n', 
               paste0(tempFeatureVectorList[noNameInd], collapse= ","), sep=" "));
  }
  
  return(bprob[tempFeatureVectorList]);

}
