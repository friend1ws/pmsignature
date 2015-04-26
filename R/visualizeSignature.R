#' @export
setGeneric("visPMSignature", function(object, sinInd = 1, baseCol) {
  standardGeneric("visPMSignature")
})

#' @export
setGeneric("visPMSignature", function(object, sinInd = 1, ...) {
  standardGeneric("visPMSignature")
})

#' @export
setGeneric("getMutNum", function(object) {
  standardGeneric("getMutNum")
})

#' @export
setGeneric("visMembership", function(object1, object2, ylog = FALSE, sortSampleNum = TRUE, multiplySampleNum = TRUE, fromSample = NULL, toSample = NULL, reorderSig = NULL, colourBrewer = NULL) {
  standardGeneric("visMembership")
})

setMethod("visPMSignature", 
          signature = c(object = "EstimatedParameters", sinInd = "numeric"), 
          function(object, sinInd = 1, ...) {
            
            vF <- object@signatureFeatureDistribution[sinInd,,];
            
            if (object@type == "independent") {
              visPMS_ind(vF, numBases = object@flankingBasesNum, trDir = object@transcriptionDirection, ...);
            } else if (object@type == "full") {
              visPMS_full(vF, numBases = object@flankingBasesNum, object@transcriptionDirection);
            }
            
          }
)



#' @title visualize probabisitic mutaiton signature for the independent model
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and flanking bases represented by the indepenent
#'   representation.
#'   
#' @param vF a matrix for mutation signature
#' @param numBases the number of flanking bases 
#' @param baseCol the colour of the bases (A, C, G, T, plus strand, minus strand)
#' @param trDir the index whether the strand direction is plotted or not
#' @param charSize the size of the character
#' @param isScale the index whether the height of the flanking base is changed or not
#' @param alpha the parameter for the Renyi entropy (applicable only if the isScale is TRUE)
visPMS_ind <- function(vF, numBases, baseCol = NA, trDir = FALSE, charSize = 1.2, isScale = FALSE, alpha = 2, charLimit = 0.25) {
  
  if (is.na(baseCol)) {
    gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
    baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1], gg_color_hue6[4], gg_color_hue6[6]);
  }
  
  centerBase <- (1 + numBases) / 2;
  
  v1 <- vF[1,1:6];
  V2 <- vF[2:(numBases),1:4];
  A <- matrix(0, numBases, 4);
  B <- matrix(0, 4, 4);
  
  if (trDir == TRUE) {
    v3 <- vF[(numBases + 1),1:2];
  }
  
  for (l in 1:numBases) {
    if (l < centerBase) {
      A[l, ] <- V2[l, ];
    } else if (l > centerBase) {
      A[l, ] <- V2[l - 1, ];
    }
  }
  A[centerBase,2] <- sum(v1[1:3]);
  A[centerBase,4] <- sum(v1[4:6]);
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3]);
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6]);
  
  num2base <- c("A", "C", "G", "T");
  
  frame();
  plot.window(xlim=c(-0.25, 1.25 * numBases + 0.25), ylim=c(-0.25, 3.25));
  
  # function for calculating renyi entropy
  renyi = function(p, tAlpha = alpha) {
    if (tAlpha == 1) {
      return(- sum(p * log2(p), na.rm = TRUE));
    } else {
      return( log(sum(p^tAlpha)) / (1 - tAlpha));
    }
  }
  
  if (isScale == FALSE) {
    sizes <- rep(1, numBases)
  } else {
    sizes <- 0.5 * (2 - apply(A, MARGIN = 1, FUN = renyi));
  }
  
  startx <- 0;
  for(l in 1:numBases) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      # polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col = baseCol[w], border=F);
      polygon(c(startx, endx, endx, startx), c(0, 0, sizes[l], sizes[l]), col = baseCol[w], border=F);
      # if (endx - startx > 1 / 4 & charSize > 0) {
      #   text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=charSize)
      # }
      if (endx - startx > charLimit & sizes[l] > 0.5 & charSize > 0) {
        text(0.5 * (endx + startx), 0.5 * sizes[l], num2base[w], col="white", cex=charSize)
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
      if ((endy - starty > charLimit) & (endx - startx > charLimit) & charSize > 0) {
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
  startx <- (numBases - 1) * 1.25 + 0.24;
  endx <- (numBases - 1) * 1.25 + 0.49;
  starty <- 2;
  endy <- starty + v3[1];
  polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[5], border=F);
  if (endy - starty > 1 / 8 & charSize > 0) {
    text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col="white", cex=charSize)
  }
  
  startx <- (numBases - 1) * 1.25 + 0.51;
  endx <- (numBases - 1) * 1.25 + 0.76;
  starty <- 2;
  endy <- starty + v3[2];
  polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=baseCol[6], border=F);
  if (endy - starty > 1 / 8 & charSize > 0) {
    text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col="white", cex=charSize)
  }
  
  }
  ##########
  
}


#' @title visualize probabisitic mutaiton signature for the full model
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and flanking bases represented by the full
#'   representation.
#'   
#' @param vF a vector for mutation signature
#' @param numBases the number of flanking bases 
#' @param trDir the index showing whether the transcription direction is used or not
visPMS_full <- function(Fvec, numBases, trDir) {
  
  flankingPatternNum <- 4^(numBases - 1);
  subPattern <- c(rep("C>A", flankingPatternNum), 
                  rep("C>G", flankingPatternNum), 
                  rep("C>T", flankingPatternNum), 
                  rep("T>A", flankingPatternNum), 
                  rep("T>C", flankingPatternNum), 
                  rep("T>G", flankingPatternNum)
                  );
  
  X <- data.frame(probability = Fvec);
  # X$strand <- factor(rep(c("plus", "minus"), 1536), levels=c("plus", "minus"));
  
  if (trDir == TRUE) {
    X$strand <- factor(c(rep("plus", flankingPatternNum * 6), rep("minus", flankingPatternNum * 6)), levels=c("plus", "minus"));
    X$subtype <- factor(rep(subPattern, 2), levels=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"));
    X$flank <- rep(1:flankingPatternNum, 2);
  } else {
    X$subtype <- factor(subPattern, levels=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"));
    X$flank <- 1:flankingPatternNum
  }
  
  
  gp <- ggplot(X, aes(x=flank, y=probability, fill=subtype)) +
    geom_bar(stat="identity", position="identity", width = 0.8) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.y = element_text(size=rel(1.2)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text= element_text(face="bold", size=rel(1.2))) +
    guides(fill=FALSE);
  
  if (trDir == TRUE) {
    gp <- gp + facet_grid(strand ~ subtype);
  } else {
    gp <- gp + facet_grid(. ~ subtype);
  }
  
  gp
  
}


setMethod("getMutNum", 
          signature = c(object = "MutationFeatureData"), 
          function(object) {         
            mutData <- data.frame(type = object@countData[1,], sampleName = object@sampleList[object@countData[2,]], count = object@countData[3,]);
            sample2mutNum <- summarize(group_by(mutData, sampleName), mutationNum = sum(count))
            return(as.data.frame(sample2mutNum))
          }
)


setMethod("visMembership",
          signature = c(object1 = "MutationFeatureData", object2 = "EstimatedParameters"),
          function(object1, object2, ylog = FALSE, sortSampleNum = TRUE, multiplySampleNum = TRUE, fromSample = NULL, toSample = NULL, reorderSig = NULL, colourBrewer = NULL) {
          
            snum <- getMutNum(object1);
            if (ylog == TRUE) {
              snum[,2] <- log10(snum[,2]);
            }  
            
            sampleList <- object2@sampleList;
            signatureNum <- object2@signatureNum;
            Q <- as.data.frame(object2@sampleSignatureDistribution)

            if (is.null(fromSample)) {
              fromSample <- 1;
            }
            if (is.null(toSample)) {
              toSample <- length(sampleList);
            }
            
            if (is.null(reorderSig)) {
              reorderSig <- 1:signatureNum;
            }
            sigOrder <- reorderSig;
            
            vMutationNum <- c()
            vSample <- c()
            vSignature <- c()

            if (sortSampleNum == TRUE) {
              mutNumOrder <- order(snum$mutationNum, decreasing = TRUE)[fromSample:toSample];
            } else {
              mutNumOrder <- fromSample:toSample;
            }

            for (k in 1:signatureNum) {
              vSample <- c(vSample, sampleList[mutNumOrder]);
              vSignature <- c(vSignature ,rep(sigOrder[k], length(mutNumOrder)));
              if (multiplySampleNum == TRUE) {
                vMutationNum <- c(vMutationNum, snum$mutationNum[mutNumOrder] * Q[mutNumOrder,k]);
              } else {
                vMutationNum <- c(vMutationNum, Q[mutNumOrder,k]);                
              }
            }
            vSample <- factor(vSample, levels = sampleList[mutNumOrder]);
            
            membership <- data.frame(sample = vSample, signature = as.factor(vSignature), mutationNum = vMutationNum);
            # membership <- data.frame(sample = reorder(vSample, -vIntensity), signature = as.factor(vSignature), intensity = vIntensity);
            
            
            gg <- ggplot(membership, aes(x = sample, y = mutationNum, fill = signature)) +
              geom_bar(width = 0.8, stat = "identity") +
              theme_bw() +
              theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank());
            if (multiplySampleNum == TRUE) {
              if (ylog == TRUE) {
                gg <- gg + ylab("log10(#mutation)");
              } else {
                gg <- gg + ylab("#mutation");
              }
            } else {
                gg <- gg + ylab("#membershipRatio");             
            }
            if (!is.null(colourBrewer)) {
              gg <- gg + scale_fill_brewer(palette = colourBrewer);
            }
            
            gg
          }

)


