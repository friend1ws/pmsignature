
#' Read the raw mutation data of mutation position format.
#' 
#' @param infile the path for the input file for the mutation data.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is considered or not
#' @param type this argument can take either independent, full, or custom.
#' 
#' @export
readMPFile_DNP <- function(infile, numBases = 3, trDir = FALSE, type = "independent") {
  
  
  if (type == "independent") {
    fdim <- c(90, rep(4, numBases - 2), rep(2, as.integer(trDir)));
  } else if (type == "full") {
    fdim <- c(90 * 4^(numBases - 2) * 2^(as.integer(trDir)));
  } else {
    stop('for reading mutation position format, the type argument has to be "independent" or "full"');
  }
  
  if (numBases %% 2 != 0) {
    stop("numBases should be even numbers");
  }
  centerInd <- c(numBases / 2, numBases / 2 + 1);
  
  mutFile <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors = FALSE);
  
  chrInfo <- mutFile[,2];
  posInfo <- mutFile[,3];
  ref_base <- Biostrings::DNAStringSet(mutFile[,4]);
  alt_base <- Biostrings::DNAStringSet(mutFile[,5]);
  sampleName_str <- as.character(mutFile[,1]);
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo, 
                                                           start = posInfo - numBases / 2 + 1 , 
                                                           end = posInfo + numBases / 2), ignore.strand = TRUE);
  
  context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr);
  
  
  removeInd <- which(XVector::subseq(context, start = centerInd[1], end = centerInd[2]) != ref_base);
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
  removeInd <- which(rowSums(alphabetFreq[,1:4]) != 2);
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
  
  DNPTable <- generateDNPInd();
  revCompInd <- which(!(as.character(XVector::subseq(context, start = centerInd[1], end = centerInd[2])) %in% DNPTable[,1]));
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
  
  
  mutFeatures <- getMutationFeatureVector_DNP(context, ref_base, alt_base, strandInfo, numBases, type);
  
  
  
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
    mutFeatList <- t(vapply(suFeatStr, function(x) as.numeric(unlist(strsplit(x, ","))), numeric(numBases + as.integer(trDir) - 1)));
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



generateDNPInd <- function() {
  
  ref <- c("CC", "CT", "CA", "CG", "TC", "TT", "TA", "AC", "AT", "GC");
  
  DNPTable <- matrix(0, 10 * 9, 2);
  tempInd <- 1;
  for (i in 1:length(ref)) {
    
    for (c1 in c("C", "T", "A", "G")) {
      for (c2 in c("C", "T", "A", "G")) {
        
        if (substr(ref[i], 1, 1) != c1 && substr(ref[i], 2, 2) != c2) {
          DNPTable[tempInd, 1] <- ref[i];
          DNPTable[tempInd, 2] <- paste(c1, c2, sep="");
          tempInd <- tempInd + 1;
        }
        
      }
    }
  }
  
  colnames(DNPTable) <- c("ref", "alt");
  return(DNPTable);
}


getMutationFeatureVector_DNP <- function(context, ref_base, alt_base, strandInfo = NULL, numBases, type) {
  
  trDir <- !is.null(strandInfo);
  if (type == "independent") {
    fdim <- c(90, rep(4, numBases - 2), rep(2, as.integer(trDir)));
  } else if (type == "full") {
    fdim <- c(90 * 4^(numBases - 2) * 2^(as.integer(trDir)));
  } else {
    stop('for reading mutation position format, the type argument has to be "independent" or "full"');
  }
  
  if (numBases %% 2 != 0) {
    stop("numBases should be even numbers");
  }
  centerInd <- c(numBases / 2, numBases / 2 + 1);
  
  DNPTable <- generateDNPInd();
  changePattern <- paste(DNPTable[,1], DNPTable[,2], sep=">");
  lookupTableChangePattern <- 1:length(changePattern);
  names(lookupTableChangePattern) <- changePattern;
  ##########
  # whether we should add additional check function for the bases of context sequences, reference and alternate allele...?
  ##########
  
  if (type == "independent") {
    
    mutFeatures <- matrix(0, length(ref_base), length(fdim));
    mutFeatures[, 1] <- lookupTableChangePattern[paste(as.character(ref_base), as.character(alt_base), sep=">")];

    
    columnInd <- 2;
    for (baseInd in 1:numBases) {
      if (baseInd %in% centerInd) {
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
    for (i in 1:((numBases - 2) / 2)) {
      
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
    
    mutFeatures[, 1] <- mutFeatures[, 1] + tempDigits * lookupTableChangePattern[paste(as.character(ref_base), as.character(alt_base), sep=">")];

    if (trDir ==TRUE) {
      tempDigits <- tempDigits * 90;
      mutFeatures[which(strandInfo == "-"), 1] <- mutFeatures[which(strandInfo == "-"), 1] + tempDigits * 1;
    }
    
  }
  
  return(mutFeatures);
  
}


visPMS_ind_DNP <- function(vF, numBases, baseCol = NA, trDir, charSize = 1.2) {
  
  if (FALSE) {
    
  if (is.na(baseCol)) {
    gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
    baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1], gg_color_hue6[4], gg_color_hue6[6]);
  }
  centerBase <- numBases / 2;
  
  DNPTable <- generateDNPInd();
  
  A_flank <- vF[2:(numBases - 1),1:4];
  
  A_center <- rep(0, 10);
  lookup_ref <- 1:length(unique(DNPTable[,1]));
  names(lookup_ref) <- unique(DNPTable[,1]);
  unique_ref <- unique(DNPTable[,1]);
  
  for (i in 1:length(DNPTable[,1])) {
    A_center[lookup_ref[DNPTable[i,1]]] <- A_center[lookup_ref[DNPTable[i,1]]] + vF[1, i];
  }
  
  B <- matrix(0, 10, 16) 
  
  if (trDir == TRUE) {
    v3 <- vF[(numBases + 1),1:2];
  }
  
  
  DNPTable <- generateDNPInd();
  num2base <- c("A", "C", "G", "T");
  base2num <- 1:4;
  names(base2num) <- num2base;
  
  frame();
  plot.window(xlim=c(-0.25, 1.25 * numBases), ylim=c(-0.25, 3.25));
  
  startx <- 0;
  for(l in 1:(numBases - 2)) {
    
    if (l == centerBase) {
      
      for (w in 1:10) {
        endx <- startx + 2 * A_center[w];
        tcol1 <- baseCol[base2num[substr(names(lookup_ref)[w], 1, 1)]];
        tcol2 <- baseCol[base2num[substr(names(lookup_ref)[w], 2, 2)]];
        polygon_zebra(c(startx, endx, endx, startx), c(0, 0, 1, 1), col1 = tcol1, col2 = tcol2, margin = 1, zwidth = 0.1);
        startx <- endx;
      }
      
      startx <- startx + 0.25;
    } 
    
    for(w in 1:4) {
       
      endx <- startx + A_flank[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col = baseCol[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=charSize)
      }
      startx <- endx;
        
    }
    startx <- startx + 0.25;
  }
  
  startx <- (centerBase - 1) * 1.25;
  for (w in 1:4) {
    starty <- 2;
    endx <- startx + A[centerBase,w];
    for(ww in 1:4) {
      endy <- starty + B[w,ww];
      polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[ww], border=F);
      if ((endy - starty > 1 / 4) & (endx - startx > 1 / 4)) {
        text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col="white", cex=charSize)
      }
      starty <- endy;
    }
    startx <- endx
    starty <- endy;
  }
  
  ##########
  # draw arrow
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + (centerBase - 1) * 1.25;
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1;
  
  polygon(xs, ys, col=8, border=F);
  ##########
  
  ##########
  if (trDir == TRUE) {
    # draw direction bias
    # startx <- (numBases - 1) * 1.25 + 0.5;
    # endx <- (numBases - 1) * 1.25 + 0.75;
    # starty <- 1.9;
    # endy <- starty + v3[1];
    # polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[5], border=F);
    
    # if (endy - starty > 1 / 8) {
    #   text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col="white", cex=1.2)
    # }
    # starty <- endy;
    # endy <- 2.9;
    # polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[6], border=F);
    # if (endy - starty > 1 / 8) {
    #   text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col="white", cex=1.2)
    # } 
    
    # draw direction bias
    startx <- (numBases - 1) * 1.25 + 0.24;
    endx <- (numBases - 1) * 1.25 + 0.49;
    starty <- 2;
    endy <- starty + v3[1];
    polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[5], border=F);
    if (endy - starty > 1 / 8) {
      text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col="white", cex=charSize)
    }
    
    startx <- (numBases - 1) * 1.25 + 0.51;
    endx <- (numBases - 1) * 1.25 + 0.76;
    starty <- 2;
    endy <- starty + v3[2];
    polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[6], border=F);
    if (endy - starty > 1 / 8) {
      text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col="white", cex=charSize)
    }
    
  }
  ##########
  
}
  
}



polygon_zebra <- function(posx, posy, col1, col2, margin, zwidth) {
  
  if (margin == 1) {
    
    starty <- 0;
    endy <- 0.5 * zwidth;
    tempCol = col1
    polygon(posx, c(starty, starty, endy, endy), col = tempCol, border = FALSE);
    starty <- endy;
    endy <- endy + zwidth;
    if (tempCol == col1) {
      tempCol <- col2;
    } else {
      tempCol <- col1;
    }
    
    while(endy < posy[3]) {
      polygon(posx, c(starty, starty, endy, endy), col = tempCol, border = FALSE);
      starty <- endy;
      endy <- endy + zwidth;
      if (tempCol == col1) {
        tempCol <- col2;
      } else {
        tempCol <- col1;
      }
    }
    
    endy <- posy[3];
    polygon(posx, c(starty, starty, endy, endy), col = tempCol, border = FALSE);
    
  } else {
    
    startx <- 0;
    endx <- 0.5 * zwidth;
    tempCol = col1
    polygon(c(startx, endx, endx, startx), posy, col = tempCol, border = FALSE);
    startx <- endx;
    endx <- endx + zwidth;
    if (tempCol == col1) {
      tempCol <- col2;
    } else {
      tempCol <- col1;
    }
    
    while(endx < posx[3]) {
      polygon(c(startx, endx, endx, startx), posy, col = tempCol, border = FALSE);
      startx <- endx;
      endx <- endx + zwidth;
      if (tempCol == col1) {
        tempCol <- col2;
      } else {
        tempCol <- col1;
      }
    }
    endx <- posx[3];
    polygon(c(startx, endx, endx, startx), posy, col = tempCol, border = FALSE); 
    
    
  }
  
}
