#' read the mutation position format file and extract the genomic position as a mutation feature
#' @param infile the path for the input mutation position format file
#' 
#' @export
readMPFile_localMutRate <- function(infile) {
  
  type = "custom";
  
  # read the genome length information and calculate the start position for each chromosom
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
  targetChr <- paste("chr", c(1:22, "X", "Y"), sep="");
  seqlen_chr <- GenomeInfoDb::seqlengths(genome)[targetChr];
  cseqlen_chr <- ceiling(seqlen_chr / 1000000);
  genomeStartInd <- cumsum(cseqlen_chr) - cseqlen_chr;
  fdim <- c(cumsum(cseqlen_chr)[length(cseqlen_chr)]); # possible number of mutation features

  # read the input file
  mutFile <- read.table(infile, sep="\t", header=FALSE);
  chrInfo <- mutFile[,2];
  posInfo <- mutFile[,3];
  sampleName_str <- as.character(mutFile[,1]);
  
  # check the 2nd column of the input feature vector (the name of chromosome)
  removeInd <- which(!(chrInfo %in% targetChr));
  if (length(removeInd) > 0) {
    warning(paste("The chromosome name is not any of chr1 ~ chr22, chrX, chrY for ", length(removeInd), " records. We have removed them."));
    sampleName_str <- sampleName_str[-removeInd];
    chrInfo <- chrInfo[-removeInd];
    posInfo <- posInfo[-removeInd];
  }
  
  suSampleStr <- sort(unique(sampleName_str));
  lookupSampleInd <- 1:length(suSampleStr);
  names(lookupSampleInd) <- suSampleStr;
  sampleIDs <- lookupSampleInd[sampleName_str];


  mutFeatures <- as.integer(genomeStartInd[as.character(chrInfo)] + ceiling(posInfo / 1000000));
  featStr <- as.character(mutFeatures);
  suFeatStr <- sort(unique(featStr));
  lookupFeatInd <- 1:length(suFeatStr);
  names(lookupFeatInd) <- suFeatStr;


  rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[featStr]);

  tableCount <- table(rawCount);
  w <- which(tableCount > 0, arr.ind=TRUE);
  procCount <- cbind(w[,2], w[,1], tableCount[w]);


  mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1);


  rownames(mutFeatList) <- NULL;
  rownames(procCount) <- NULL;
  
  return(new(Class = "MutationFeatureData", 
             type = "custom",
             flankingBasesNum = as.integer(1),
             transcriptionDirection = FALSE,
             possibleFeatures = as.integer(fdim),
             featureVectorList = t(mutFeatList),
             sampleList = suSampleStr,
             countData = t(procCount),
             mutationPosition = data.frame()
  ))
  
}




