#' get mutation feature vector from context sequence data and reference and alternate allele information
#'
#' @param type the type of mutation feature vecotr (should be "independent" or "full").
#' @param numBases the number of flanking bases around the mutated position.
#' @param trDir the index representing whether transcription direction is considered or not
#' @param trial the number of randome site generations
#' @export
getBackgroudSignature <- function(type = "independent", numBases = 3, trDir = FALSE, trial = 1000000) {

  if (!(type %in% c("independent", "full"))) {
    stop('the parameter type should be "independent" or "full"')
  }
  
  if (numBases %% 2 != 1) {
    stop("numBases should be odd numbers")
  }
  centerInd <- (numBases + 1) / 2
  
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  targetChr <- paste("chr", c(1:22, "X", "Y"), sep="")
  seqlen_chr <- GenomeInfoDb::seqlengths(genome)[targetChr]
  
  prob_chr <- seqlen_chr / 100 # for avoiding over flow..
  prob_chr <- prob_chr / sum(prob_chr)
  
  
  chr_trial <- sample(targetChr, trial, replace = TRUE, prob = prob_chr)
  start_trial <- rep(0, trial)
  for (seq in targetChr) {
    start_trial[chr_trial == seq] <- sample.int(seqlen_chr[seq] - numBases + 1, sum(chr_trial == seq), replace = TRUE)
  }
  end_trial <- start_trial + numBases - 1
  context_trial <- Biostrings::getSeq(genome, name = chr_trial, start = start_trial, end = end_trial)
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr_trial, 
                                            start = start_trial, 
                                            end = end_trial), ignore.strand = TRUE)

  
  alphabetFreq <- Biostrings::alphabetFrequency(context_trial)
  removeInd <- which(alphabetFreq[,"A"] + alphabetFreq[,"C"] + alphabetFreq[,"G"] + alphabetFreq[,"T"] != numBases)
  if (length(removeInd) > 0) {
    context_trial <- context_trial[-removeInd]
    chr_trial <- chr_trial[-removeInd]
    start_trial <- start_trial[-removeInd]
    end_trial <- end_trial[-removeInd]
    warning(paste("The characters other than (A, C, G, T) are included in flanking bases of", length(removeInd), "mutations. We have removed them."))
  }
  
  
  revCompInd <- which(as.character(XVector::subseq(context_trial, start = centerInd, end = centerInd)) %in% c("A", "G"))
  context_trial[revCompInd] <- Biostrings::reverseComplement(context_trial[revCompInd])
  
  # if (trDir == TRUE) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr_trial, 
                                            start = start_trial, 
                                            end = end_trial), ignore.strand = TRUE)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    exons_txdb <- GenomicFeatures::exons(txdb)
    gr_txdb <- GenomicRanges::findOverlaps(gr, exons_txdb, ignore.strand = FALSE)
    gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(exons_txdb[gr_txdb@subjectHits]))))
    # txdb_bed <- GenomicFeatures::asBED(txdb)
    # gr_txdb <- GenomicRanges::findOverlaps(gr, txdb_bed, ignore.strand = FALSE)
    # gr_strand <- cbind(S4Vectors::queryHits(gr_txdb), as.character(S4Vectors::as.factor(BiocGenerics::strand(txdb_bed[gr_txdb@subjectHits]))))
    
    ugr_strand <- unique(gr_strand[gr_strand[, 2] != "*" ,], MARGIN=1)
  
    rmdup_ugr_strand <- ugr_strand[!duplicated(ugr_strand[, 1]), ]
    txdb_plus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "+", 1])
    txdb_minus_gr_ind <- as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "-", 1])
  
    strandInfo_trial <- rep("*", length(gr))
    strandInfo_trial[setdiff(txdb_plus_gr_ind, revCompInd)] <- "+"
    strandInfo_trial[intersect(txdb_plus_gr_ind, revCompInd)] <- "-"  
    strandInfo_trial[setdiff(txdb_minus_gr_ind, revCompInd)] <- "-"
    strandInfo_trial[intersect(txdb_minus_gr_ind, revCompInd)] <- "+" 
  
    warning(paste("Out of", length(context_trial), "mutations, we could obtain transcription direction information for", 
                length(txdb_plus_gr_ind) + length(txdb_minus_gr_ind), "mutation. Other mutations are removed."))
    context_trial <- context_trial[strandInfo_trial != "*"]
    chr_trial <- chr_trial[strandInfo_trial != "*"]
    start_trial <- start_trial[strandInfo_trial != "*"]
    end_trial <- end_trial[strandInfo_trial != "*"]
    strandInfo_trial <- strandInfo_trial[strandInfo_trial != "*"]
  # }
  
  if (trDir == FALSE) {
    strandInfo_trial <- NULL
  }
  
  ref_base_trial <- XVector::subseq(context_trial, start = centerInd, end = centerInd)
  
  alt_base_trial <- sample.int(3, length(ref_base_trial), replace = TRUE)
  alt_base_trial[ref_base_trial == "C" & alt_base_trial == 1] <- "A"
  alt_base_trial[ref_base_trial == "C" & alt_base_trial == 2] <- "G"
  alt_base_trial[ref_base_trial == "C" & alt_base_trial == 3] <- "T"
  alt_base_trial[ref_base_trial == "T" & alt_base_trial == 1] <- "A"
  alt_base_trial[ref_base_trial == "T" & alt_base_trial == 2] <- "C"
  alt_base_trial[ref_base_trial == "T" & alt_base_trial == 3] <- "G"
  alt_base_trial <- Biostrings::DNAStringSet(alt_base_trial)
  
  mutFeatures <- getMutationFeatureVector(context_trial, ref_base_trial, alt_base_trial, strandInfo_trial, numBases, type)
 

  featStr <- apply(mutFeatures, 1, paste0, collapse=",")
  
  suFeatStr <- sort(unique(featStr))
  lookupFeatInd <- 1:length(suFeatStr)
  names(lookupFeatInd) <- suFeatStr
  rawCount <- table(lookupFeatInd[featStr])
  
  tempBgProb <- rawCount / sum(rawCount)
  names(tempBgProb) <- suFeatStr
  
  
  # first all the possible mutation features are listed
  if (type == "independent") {
    allFeatStr <- 1:6
    for(i in 1:(numBases - 1)) {
      allFeatStr <- as.vector(outer(allFeatStr, 1:4, FUN = function(x, y) {paste(x, y, sep=",")}))
    }
    if (trDir == TRUE) {
      allFeatStr <- as.vector(outer(allFeatStr, 1:2, FUN = function(x, y) {paste(x, y, sep=",")}))     
    }
  } else {
    allFeatStr <- 1:(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
  }
  
  # then, make the vector whose dimension is the number of all the possible features and put the probabilities obtained above
  bgProb <- rep(0, length(allFeatStr))
  names(bgProb) <- sort(allFeatStr)
  bgProb[names(tempBgProb)] <- tempBgProb
  
  return(bgProb)
  
  
}


#' update the background signature file in the package system.file.
#'
#' @param type the type of mutation feature vecotr (should be "independent" or "full").
#' @param numBases the number of flanking bases around the mutated position.
#' @param trDir the index representing whether transcription direction is considered or not
#' @param trial the number of randome site generations
updateBackgroudSignatureFile <- function(type = "independent", numBases = 3, trDir = FALSE, trial = 1000000) {
  
  if (type == "independent") {
    tempType <- "ind"
  } else if (type == "full") {
    tempType <- "full"
  } else {
    stop('the parameter type should be "independent" or "full"')
  }
  
  if (!(numBases %in% c(3, 5, 7, 9))) {
    stop('Background data whose number of flanking bases is other than 3, 5, 7 or 9 is not available')
  }
  
  if (trDir == TRUE) {
    bgfile <- paste("inst/bgdata/bg.", tempType, numBases, "_dir.txt", sep="")
  } else {
    bgfile <- paste("inst/bgdata/bg.", tempType, numBases, ".txt", sep="")
  }
  
  bgProb <- getBackgroudSignature(type = type, numBases = numBases, trDir = trDir, trial = trial)
  
  write.table(format(round(bgProb, 8), digit = 8, nsmall = 8), file = bgfile, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

}
  
  