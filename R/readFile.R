#' Read the raw mutation data with the mutation feature vector format.
#' 
#' @description
#' The mutation feature vector format is tab-delimited format, 
#' where the 1st column shows the name of samples, and the 2nd-last columns show the values of mutation features.
#' 
#' When type = "independent", the mutation feature should be equal to \code{(numBases + as.integer(trDir))}-dimensional vector,
#' representing substitution patterns (1 to 6, C>A, C>G, C>T, T>A, T>C and T>G), 
#' 5' and 3' flanking bases (1 to 4, A, C, G and T),
#' and transcription direction (1 to 2, + and -), in this order.
#' 
#' When type = "full", the mutation feature should be equal to \code{1}-dimensional vector,
#' taking the value of 1 to \code{6 * 4^{(numBases - 1)} * 2^{as.integer(trDir)}}.
#' For the coding of the number to the actual mutation pattern (substitution patterns, flanking bases, transcription direction),
#' seeing the source code of \code{getMutationFeatureVector} may be helpful.
#' 
#' Also, this function usually can accept compressed files (e.g., by gzip, bzip2 and so on) when using recent version of R.
#' 
#' @param infile the path for the input file for the mutation feature data.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is considered or not.
#' @param type this argument can take either "independent", "full", or "custom".
#' The values "independent" or "full" can be set only when the input file is decently formatted.
#' 
#' 
#' @examples 
#' We have example data in the pmsignature package;
#' inputFile <- system.file("extdata/Hoang_MFVF.ind.txt", package="pmsignature")
#' G <- readMFVFile(inputFile, numBases = 5, type="independent", trDir=TRUE)
#' 
#' @export
readMFVFile <- function(infile, numBases = 3, trDir = FALSE, type = "custom") {
  
  if (!(type %in% c("independent", "full", "custom"))) {
    stop('The argument type should be eigher "independent", "full" or "custom"')
  }
  
  if (type == "custom") {
    if (!missing(numBases)) {
      warning('the argument numBases is applicable only when type = "independent", or "full"')
    }
    if (!missing(trDir)) {
      warning('the argument trDir is applicable only when type = "independent", or "full"')
    }    
  }

  mutFile <- read.table(infile, sep="\t", header=FALSE)
  sampleName_str <- as.character(mutFile[,1])
  mutFeatures <- mutFile[,2:ncol(mutFile), drop = FALSE]
  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
  } else {
    mutFeatures <- mutFile[,2:ncol(mutFile), drop = FALSE]
    fdim <- unname(apply(mutFeatures, 2, max))
  }
  
  # check the consistency of the dimension of mutation features
  if (length(fdim) != ncol(mutFile) - 1) {
    stop(paste("The dimension of mutation features should be", length(fdim), "for numBases:", numBases, ",trDir:", trDir, ",type:", type))
  }

  
  if (type == "independent") {

    if (numBases %% 2 != 1) {
      stop("numBases should be odd numbers")
    }
    
    # check the dimension of the mutation feature vector
    if (ncol(mutFeatures) < numBases + as.integer(trDir)) {
      stop('When type = "independent", the dimension of feature vectors should be at least (numBases + as.integer(trDir))')
    }
    
    # check the 1-st column of the input feature vector (substitution type, should be 1:6)
    if (any(!(mutFeatures[,1] %in% 1:6))) {
      stop(paste('When type = "independent", the 1st feature values should be 1 to 6\n',
                 'In the 1st column, the following rows have values other than 1 to 6;\n',
                 paste0(which(!(mutFeatures[,1] %in% 1:6)), collapse = ","), sep=""))
    }
    
    # check the 2nd to (numBases)-th columns of the input feature vector (flanking bases, should be 1:4)
    for (i in 2:numBases) {
      if (any(!(mutFeatures[,i] %in% 1:4))) {
        stop(paste('When type = "independent", the 2 to (numBases)-th feature values should be 1 to 4\n',
                   'In the ', i, '-th column, the following rows have values other than 1 to 4;\n', 
                   paste0(which(!(mutFeatures[,i] %in% 1:4)), collapse= ","), sep=""))
      }
    }
    
    # check the (numBases + 1)-th column of the input feature vector (the flag for transcription direction, should be 1:2)
    if (trDir == TRUE) {
      if (any(!(mutFeatures[,numBases + 1] %in% 1:2))) {
        stop(paste('When type = "independent" and trDir = TRUE, ', 
                   'the (numBases + 1)-th feature values should be 1 to 2\n',
                   'In the ', (numBases + 1), '-th column, following rows have values other than 1 to 2;\n', 
                   paste0(which(!(mutFeatures[,numBases + 1] %in% 1:2)), collapse = ","), sep=""))
      }
    }

  }
    
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
  
  if (length(fdim) == 1) {
    mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1)
  } else {
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(length(fdim))))
  }
  
  rownames(mutFeatList) <- NULL
  rownames(procCount) <- NULL

  
  return(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = data.frame()
  ))
  
}


#' Read the raw mutation data of Mutation Position Format. 
#' 
#' @description
#' The mutation position format is tab-delimited text file, where 
#' the 1st-5th columns shows sample names, chromosome names, 
#' coordinates, reference bases (A, C, G, or T) and
#' the alternate bases (A, C, G, or T), respectively. An example is as follows;
#' 
#' ---
#' 
#' sample1 chr1 100 A C
#' 
#' sample1 chr1 200 A T
#' 
#' sample1 chr2 100 G T
#' 
#' sample2 chr1 300 T C
#' 
#' sample3 chr3 400 T C
#' 
#' ---
#' 
#' Also, this function usually can accept compressed files (e.g., by gzip, bzip2 and so on) when using recent version of R.
#' Currently, only UCSC hg19 (BSgenome.Hsapiens.UCSC.hg19) is supported. 
#' 
#' @param infile the path for the input file for the mutation data of Mutation Position Format.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is considered or not.
#' The gene annotation information is given by UCSC knownGene (TxDb.Hsapiens.UCSC.hg19.knownGene object).
#' When trDir is TRUE, the mutations located in intergenic region are excluded from the analysis.
#' @param type this argument can take either "independent", "full", or "custom".
#' @param bs_genome this argument specifies the reference genome (e.g., BSgenome.Mmusculus.UCSC.mm10 can be used for the mouse genome).
#' See https://bioconductor.org/packages/release/bioc/html/BSgenome.html for the available genome list
#' @param txdb_transcript this argument specified the transcript database (e.g., TxDb.Mmusculus.UCSC.mm10.knownGene can be used for the mouse genome).
#' See https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html for details.
#' @return The output is an instance of MutationFeatureData S4 class (which stores summarized information on mutation data).
#' This will be typically the input of \code{getPMSignature} function for estimating the parameters of mutation signatures and memberships.
#' 
#' @examples 
#' We have example data in the pmsignature package;
#' inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
#' G <- readMPFile(inputFile, numBases = 5)
#' 
#' When assuming non-independent model;
#' G <- readMPFile(inputFile, numBases = 5, type = "full")
#' 
#' When adding transcription direction information;
#' G <- readMPFile(inputFile, numBases = 5, trDir = TRUE)
#' 
#' @export
readMPFile <- function(infile, numBases = 3, trDir = FALSE, type = "independent", 
                       bs_genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                       txdb_transcript = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene) {

  
  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
  } else {
    stop('for reading mutation position format, the type argument has to be "independent" or "full"')
  }

  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers")
  }
  centerInd <- (numBases + 1) / 2

  mutFile <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors = FALSE)

  chrInfo <- mutFile[,2]
  posInfo <- mutFile[,3]
  ref_base <- Biostrings::DNAStringSet(mutFile[,4])
  alt_base <- Biostrings::DNAStringSet(mutFile[,5])
  sampleName_str <- as.character(mutFile[,1])
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo, 
                                            start = posInfo, 
                                            end = posInfo), ignore.strand = TRUE)

  ranges <- GenomicRanges::resize(gr, numBases, fix = "center")
  context <- Biostrings::getSeq(bs_genome, ranges)


  # check the consistency between the input reference base and the obtained base using hg19 reference genome.
  removeInd <- which(XVector::subseq(context, start = centerInd, end = centerInd) != ref_base)
  if (length(removeInd) > 0) {
    warning(paste("The central bases are inconsistent in", length(removeInd), "mutations. We have removed them."))
    context <- context[-removeInd]
    ref_base <- ref_base[-removeInd]
    alt_base <- alt_base[-removeInd]
    sampleName_str <- sampleName_str[-removeInd]
    chrInfo <- chrInfo[-removeInd]
    posInfo <- posInfo[-removeInd]
  }

  # check the characters on alternative base
  alphabetFreq <- Biostrings::alphabetFrequency(alt_base)
  removeInd <- which(rowSums(alphabetFreq[,1:4]) != 1)
  if (length(removeInd) > 0) {
    warning(paste("The characters other than (A, C, G, T) are included in alternate bases of", length(removeInd), "mutations. We have removed them."))
    context <- context[-removeInd]
    ref_base <- ref_base[-removeInd]
    alt_base <- alt_base[-removeInd]
    sampleName_str <- sampleName_str[-removeInd]
    chrInfo <- chrInfo[-removeInd]
    posInfo <- posInfo[-removeInd]
  }

  # check the characters on flanking bases
  alphabetFreq <- Biostrings::alphabetFrequency(context)
  removeInd <- which(alphabetFreq[,"A"] + alphabetFreq[,"C"] + alphabetFreq[,"G"] + alphabetFreq[,"T"] != numBases)
  if (length(removeInd) > 0) {
    context <- context[-removeInd]
    ref_base <- ref_base[-removeInd]
    alt_base <- alt_base[-removeInd]
    sampleName_str <- sampleName_str[-removeInd]
    chrInfo <- chrInfo[-removeInd]
    posInfo <- posInfo[-removeInd]
    warning(paste("The characters other than (A, C, G, T) are included in flanking bases of", length(removeInd), "mutations. We have removed them."))
  }
  
  # check the characters on alternative base
  alphabetFreq <- Biostrings::alphabetFrequency(alt_base)
  removeInd <- which(ref_base == alt_base)
  if (length(removeInd) > 0) {
    warning(paste("The reference base and alternative bases are equal for", length(removeInd), "mutations. We have removed them."))
    context <- context[-removeInd]
    ref_base <- ref_base[-removeInd]
    alt_base <- alt_base[-removeInd]
    sampleName_str <- sampleName_str[-removeInd]
    chrInfo <- chrInfo[-removeInd]
    posInfo <- posInfo[-removeInd]
  }

  
  revCompInd <- which(as.character(XVector::subseq(context, start = centerInd, end = centerInd)) %in% c("A", "G"))
  context[revCompInd] <- Biostrings::reverseComplement(context[revCompInd])
  ref_base[revCompInd] <- Biostrings::reverseComplement(ref_base[revCompInd])
  alt_base[revCompInd] <- Biostrings::reverseComplement(alt_base[revCompInd])

  # Obtaining transcription strand information using GenomicRanges packages
  if (trDir == TRUE) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo, 
                                                             start = posInfo, 
                                                             end = posInfo), ignore.strand = TRUE)
    
    # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    txdb <- txdb_transcript
    # exons_txdb <- GenomicFeatures::exons(txdb)
    # gr_txdb <- GenomicRanges::findOverlaps(gr, exons_txdb, ignore.strand = FALSE)
    # gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(exons_txdb[gr_txdb@subjectHits]))))
    
    # txdb_bed <- GenomicFeatures::asBED(txdb)
    txdb_bed <- GenomicFeatures::transcripts(txdb)
    gr_txdb <- GenomicRanges::findOverlaps(gr, txdb_bed, ignore.strand = FALSE)
    # gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(txdb_bed[gr_txdb@subjectHits]))))
    gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(txdb_bed[S4Vectors::subjectHits(gr_txdb)]))))
    ugr_strand <- unique(gr_strand[gr_strand[, 2] != "*" ,], MARGIN=1)
    
    rmdup_ugr_strand <- ugr_strand[!duplicated(ugr_strand[, 1]), ]
    txdb_plus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "+", 1])
    txdb_minus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "-", 1])
    
    strandInfo <- rep("*", length(gr))
    strandInfo[setdiff(txdb_plus_gr_ind, revCompInd)] <- "+"
    strandInfo[intersect(txdb_plus_gr_ind, revCompInd)] <- "-"
    strandInfo[setdiff(txdb_minus_gr_ind, revCompInd)] <- "-"
    strandInfo[intersect(txdb_minus_gr_ind, revCompInd)] <- "+"
    
    warning(paste("Out of", length(context), "mutations, we could obtain transcription direction information for", 
                  length(txdb_plus_gr_ind) + length(txdb_minus_gr_ind), "mutation. Other mutations are removed."))
    
    # this may be unstable... need to review thoroughly later...
    context <- context[strandInfo != "*"]
    ref_base <- ref_base[strandInfo != "*"]
    alt_base <- alt_base[strandInfo != "*"]
    sampleName_str <- sampleName_str[strandInfo != "*"]
    chrInfo <- chrInfo[strandInfo != "*"]
    posInfo <- posInfo[strandInfo != "*"]
    strandInfo <- strandInfo[strandInfo != "*"]
  }
  
  
  if (trDir == FALSE) {
    strandInfo <- NULL
  }
  
  
  mutFeatures <- getMutationFeatureVector(context, ref_base, alt_base, strandInfo, numBases, type)

    
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
  
  
  if (length(fdim) == 1) {
    # mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(1)))
    mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1)
  } else {
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(numBases + as.integer(trDir))))
  } 

  rownames(mutFeatList) <- NULL
  rownames(procCount) <- NULL

  if (trDir == FALSE) {
    strandInfo_for_class <- rep(NA, length(chrInfo))
  } else {
    strandInfo_for_class <- strandInfo
  }
  
  return(new(Class = "MutationFeatureData", 
             type = type,
             flankingBasesNum = as.integer(numBases),
             transcriptionDirection = trDir,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = data.frame(chr = chrInfo, pos = posInfo, ref = ref_base, alt = alt_base, strand = strandInfo_for_class, context = context, sampleID = unname(sampleIDs), mutID = unname(lookupFeatInd[featStr]), stringsAsFactors = FALSE)
             )
  )

}



#' Get mutation feature vector from context sequence data and reference and alternate allele information
#'
#' @param context the context sequence data around the mutated position. This shoud be Biostrings::DNAStringSet class
#' @param ref_base the reference bases at the mutated position.
#' @param alt_base the alternate bases at the mutated position.
#' @param strandInfo transcribed strand information at the mutated position. (this is optional)
#' @param numBases the number of flanking bases around the mutated position.
#' @param type the type of mutation feature vecotr (should be "independent" or "full").
#' 
#' @export
getMutationFeatureVector <- function(context, ref_base, alt_base, strandInfo = NULL, numBases, type) {
  
  trDir <- !is.null(strandInfo)
  if (type == "independent") {
    fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
  } else if (type == "full") {
    fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
  } else {
    stop('the type argument has to be "independent" or "full"')
  }
  
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers")
  }
  centerInd <- (numBases + 1) / 2
  
  
  ##########
  # whether we should add additional check function for the bases of context sequences, reference and alternate allele...?
  ##########
  
  if (type == "independent") {
    
    mutFeatures <- matrix(0, length(ref_base), length(fdim))
    
    mutFeatures[which(ref_base == "C" & alt_base == "A"), 1] <- 1
    mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] <- 2
    mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] <- 3
    mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] <- 4
    mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] <- 5
    mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] <- 6
    
    columnInd <- 2
    for (baseInd in 1:numBases) {
      if (baseInd == centerInd) {
        next
      }
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "A"), columnInd] <- 1
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), columnInd] <- 2
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), columnInd] <- 3
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), columnInd] <- 4
      columnInd <- columnInd + 1
    }
    
    if (trDir ==TRUE) {
      mutFeatures[which(strandInfo == "+"), length(fdim)] <- 1
      mutFeatures[which(strandInfo == "-"), length(fdim)] <- 2
    }
    
  } else {
    
    mutFeatures <- matrix(1, length(ref_base), length(fdim))
    
    tempDigits <- 1
    for (i in 1:((numBases - 1) / 2)) {
      
      baseInd <- numBases + 1 - i
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] + tempDigits * 1
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] + tempDigits * 2
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] + tempDigits * 3
      
      tempDigits <- tempDigits * 4
      baseInd <- i;      
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "C"), 1] + tempDigits * 1
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "G"), 1] + tempDigits * 2
      mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] <- mutFeatures[which(XVector::subseq(context, start = baseInd, end = baseInd) == "T"), 1] + tempDigits * 3
      
      tempDigits <- tempDigits * 4
    }
    
    mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] <- mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] + tempDigits * 1
    mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] <- mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] + tempDigits * 2
    mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] + tempDigits * 3
    mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] + tempDigits * 4
    mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] <- mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] + tempDigits * 5
    
    if (trDir ==TRUE) {
      tempDigits <- tempDigits * 6
      mutFeatures[which(strandInfo == "-"), 1] <- mutFeatures[which(strandInfo == "-"), 1] + tempDigits * 1
    }
    
  }
  
  return(mutFeatures)
  
}

#' Read and format the background vector data
#' 
#' @description
#' This function provides a background probability for each mutation feature,
#' that is available only when the model type is "independent" or "full",
#' and the numBases is either of 3, 5, 7 or 9.
#'
#' The background probability vectors are calculated by the function \code{getBackgroudSignature},
#' checking the frequencies of consecutive nucleotides on \strong{exonic regions}. 
#' Therefore, when you are using whole genome sequencing data, 
#' using the background signatures provided by this function may not be appropriate.
#' 
#' @param mutationFeatureData the mutation data processed in the function (\code{readMPFile} or \code{readMFVFile})
#' 
#' @return
#' The output is a background probability vector corresponding to the input mutation data.
#' 
#' @examples 
#' inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
#' G <- readMPFile(inputFile, numBases = 7)
#' BG_prob <- readBGFile(G)
#' 
#' @export
readBGFile <- function(mutationFeatureData) {
  
  if (slot(mutationFeatureData, "type") == "independent") {
    tempType <- "ind"
  } else if (slot(mutationFeatureData, "type") == "full") {
    tempType <- "full"
  } else {
    stop('Background data for types other than "independent" or "full" is not available')
  }
  
  if (slot(mutationFeatureData, "flankingBasesNum") %in% c(3, 5, 7, 9)) {
    tempNumBase <- slot(mutationFeatureData, "flankingBasesNum")
  } else {
    stop('Background data whose number of flanking bases is other than 3 or 5 is not available')
  }
  
  if (slot(mutationFeatureData, "transcriptionDirection") == TRUE) {
    bgfile <- paste("bgdata/bg.", tempType, tempNumBase, "_dir.txt", sep="")
  } else {
    bgfile <- paste("bgdata/bg.", tempType, tempNumBase, ".txt", sep="")
  }
  
  bdata <- read.table(system.file(bgfile, package = "pmsignature"), header = FALSE, sep="\t")
  
  tempFeatureVectorList <- apply(slot(mutationFeatureData, "featureVectorList"), 2, paste0, collapse=",")
  bprob <- bdata[,2] / sum(bdata[,2]) # normalize sum to one
  names(bprob) <- bdata[,1]
  
  if (!all(tempFeatureVectorList %in% names(bprob))) {
    noNameInd <- which(! (tempFeatureVectorList %in% names(bprob)))
    stop(paste('The information of following mutation features are not included in the specified background file:\n', 
               paste0(tempFeatureVectorList[noNameInd], collapse= ","), sep=" "))
  }
  
  return(bprob[tempFeatureVectorList])

}
