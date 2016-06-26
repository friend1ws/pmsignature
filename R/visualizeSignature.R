#' Visualize estimated probabilistic mutation signatures
#' 
#' @description By checking the meta-information on the mutation signature 
#' (model type, #flanking-bases, transcription strand usages), 
#' appropriate visualization will be automatically produced
#' (but currently, we are supporting just for "independent" and "full" type models).
#' 
#' @param object EstimatedParameters class
#' @param sigInd the index of the mutation signature plotted
#' @param baseCol the colors of the bases (A, C, G, T, plus strand, minus strand), (default: ggplot2 default colors, active only if type = "independent")
#' @param charSize the size of the character passed to geom_text function (default: 5 (active only if type = "independent")
#' @param isScale the index whether the height of the flanking base is changed or not (default: FALSE, active only if type = "independent")
#' @param alpha the parameter for the Renyi entropy (default: 2, applicable only if the isScale == TRUE)
#' 
#' @return a figure of estimated probabilistic mutation signature via ggplot2 is generated
#' (therefore, can be saved using \code{ggsave} function).
#' @examples 
#' After obtaining EstimatedParameters (typically by \code{getPMSignature}) as Param,
#' visPMSignature(Param, 1)
#' 
#' You can change the heights of flanking bases according to information contents,
#' visPMSignature(Param, 1, isScale = TRUE)
#' 
#' @export
setGeneric("visPMSignature", function(object, sinInd = 1, ...) {
  standardGeneric("visPMSignature")
})


#' Obtain somatic mutation count for each sample
#' 
#' @param object MutationFeatureData class 
#' 
#' @return the number of somatic mutations for each cancer sample 
#' stored into the instance of mutationFeatureData class. 
#' 
#' @examples 
#' After obtaining mutationFeatureData (see e.g., readMPFile function) as G,
#' mutNum <- getMutNum(G)
#' print(mutNum)
#' 
#' @export
setGeneric("getMutNum", function(object) {
  standardGeneric("getMutNum")
})



#' Visualize estimated membership parameters
#' 
#' @param object1 MutationFeatureData class 
#' @param object2 EstimatedParameters class
#' @param multiplySampleNum barplot height for each sample is multiplied by #mutations (default: TRUE).
#' @param ylog barplot height is an logarithms of #mutations or not 
#' (default: FALSE, active only if multiplySampleNum = TRUE).
#' @param sortSampleNum samples are sorted according to #mutations (default: TRUE).
#' @param fromSample only samples from the specified index will be plotted (default: NULL).
#' @param fromSample only samples until the specified index will be plotted (default: NULL).
#' @param reorderSig the order of signatures are reordered according to the specified order (default: NULL).
#' @param colourBrewer colourBrewer palette set passed to scale_fill_brewer function.
#' See, e.g., \url{http://docs.ggplot2.org/current/scale_brewer.html} for detail.
#' 
#' @return a figure of estimated membership parameter via ggplot2 is generated
#' (therefore, can be saved using \code{ggsave} function).
#' 
#' @examples 
#' After obtaining EstimatedParameters (typically by \code{getPMSignature}) as Param,
#' visPMSignature(G, Param)
#' 
#' You can equate the heights of barplot
#' visPMSignature(G,Param, multiplySampleNum = TRUE)
#' 
#' Use colourBrewer palette,
#' visPMSignature(G,Param, colourBrewer = "Set2")
#' 
#' @export
setGeneric("visMembership", function(object1, object2, ylog = FALSE, sortSampleNum = TRUE, multiplySampleNum = TRUE, fromSample = NULL, toSample = NULL, reorderSig = NULL, colourBrewer = NULL) {
  standardGeneric("visMembership")
})

setMethod("visPMSignature", 
          signature = c(object = "EstimatedParameters", sinInd = "numeric"), 
          function(object, sinInd = 1, ...) {
            
            vF <- object@signatureFeatureDistribution[sinInd,,]
            
            if (object@type == "independent") {
              visPMS_ind(vF, numBases = object@flankingBasesNum, trDir = object@transcriptionDirection, ...)
            } else if (object@type == "full") {
              visPMS_full(vF, numBases = object@flankingBasesNum, object@transcriptionDirection)
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
visPMS_ind <- function(vF, numBases, baseCol = NA, trDir = FALSE, charSize = 5, isScale = FALSE, alpha = 2, charLimit = 0.25) {
  
  if (is.na(baseCol)) {
    gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
    baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1], gg_color_hue6[4], gg_color_hue6[6])
  }
  
  centerBase <- (1 + numBases) / 2
  
  v1 <- vF[1,1:6]
  V2 <- vF[2:(numBases),1:4]
  A <- matrix(0, numBases, 4)
  B <- matrix(0, 4, 4)
  
  if (trDir == TRUE) {
    v3 <- vF[(numBases + 1),1:2]
  }
  
  for (l in 1:numBases) {
    if (l < centerBase) {
      A[l, ] <- V2[l, ]
    } else if (l > centerBase) {
      A[l, ] <- V2[l - 1, ]
    }
  }
  A[centerBase,2] <- sum(v1[1:3])
  A[centerBase,4] <- sum(v1[4:6])
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3])
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6])
 
  renyi <- function(p, tAlpha = alpha) {
    if (tAlpha == 1) {
      return(- sum(p * log2(p), na.rm = TRUE))
    } else {
      return( log(sum(p^tAlpha)) / (1 - tAlpha))
    }
  }
  
  if (isScale == FALSE) {
    fheight <- rep(1, numBases)
  } else {
    fheight <- 0.5 * (2 - apply(A, MARGIN = 1, FUN = renyi))
  }

  ##########
  # collecting data for ggplot
  x_start <- c()
  x_end <- c()
  y_start <- c()
  y_end <- c()
  text_x <- c()
  text_y <- c()
  text_lab <- c()
  text_col <- c()
  rectType <- c()
  num2base <- c("A", "C", "G", "T")
  
  # flanking bases
  tempStartX <- 0
  for (i in 1:numBases) {
    x_start <- c(x_start, tempStartX + c(0, cumsum(A[i,1:3])))
    x_end <- c(x_end, tempStartX + cumsum(A[i,1:4]))
    y_start <- c(y_start, rep(0, 4))
    y_end <- c(y_end, rep(fheight[i], 4))
    rectType <- c(rectType, c("A", "C", "G", "T")) 
    for (j in 1:4) {
      tempPos <- c(0, cumsum(A[i,1:4]))
      if (A[i,j] > charLimit && fheight[i] > charLimit) {
        text_x <- c(text_x, tempStartX + 0.5 * (tempPos[j] + tempPos[j + 1]))
        text_y <- c(text_y, 0.5 * (0 + fheight[i]))
        text_lab <- c(text_lab, num2base[j])
        text_col <- c(text_col, "w")
      }
    }
    tempStartX <- tempStartX + 1.25
  }  
  
  # alternative bases from C
  tempStartX <- (centerBase - 1) * 1.25
  x_start <- c(x_start, rep(tempStartX, 4))
  x_end <- c(x_end, rep(tempStartX + A[centerBase, 2], 4))
  y_start <- c(y_start, 2 + c(0, cumsum(B[2,1:3])))
  y_end <- c(y_end, 2 + cumsum(B[2,1:4]))
  rectType <- c(rectType, c("A", "C", "G", "T"))
  
  tempPos <- c(0, cumsum(B[2,1:4]))
  for (j in 1:4) {
    if (A[centerBase, 2] > charLimit && B[2,j] > charLimit) {
      text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 2])
      text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
      text_lab <- c(text_lab, num2base[j])
      text_col <- c(text_col, "w")
    }
  }
  
  # alternative bases from T
  tempStartX <- tempStartX + A[centerBase, 2]
  x_start <- c(x_start, rep(tempStartX, 4))
  x_end <- c(x_end, rep(tempStartX + A[centerBase, 4], 4))
  y_start <- c(y_start, 2 + c(0, cumsum(B[4,1:3])))
  y_end <- c(y_end, 2 + cumsum(B[4,1:4]))
  rectType <- c(rectType, c("A", "C", "G", "T"))
  
  tempPos <- c(0, cumsum(B[4,1:4]))
  for (j in 1:4) {
    if (A[centerBase, 4] > charLimit && B[4,j] > charLimit) {
      text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 4])
      text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
      text_lab <- c(text_lab, num2base[j])
      text_col <- c(text_col, "w")
    }
  }
  
  if (trDir == TRUE) {
    
     # draw direction bias
     x_start <- c(x_start, (numBases - 1) * 1.25 + 0.24)
     x_end <- c(x_end, (numBases - 1) * 1.25 + 0.49)
     y_start <- c(y_start, 2)
     y_end <- c(y_end, 2 + v3[1])
     rectType <- c(rectType, c("+"))
     
     if (v3[1] > 0.125) {
       text_x <- c(text_x, (numBases - 1) * 1.25 + 0.5 * (0.24 + 0.49))
       text_y <- c(text_y, 2 + 0.5 * v3[1])
       text_lab <- c(text_lab, "+")
       text_col <- c(text_col, "w")
     }
     

     x_start <- c(x_start, (numBases - 1) * 1.25 + 0.51)
     x_end <- c(x_end, (numBases - 1) * 1.25 + 0.76)
     y_start <- c(y_start, 2)
     y_end <- c(y_end, 2 + v3[2])
     rectType <- c(rectType, c("-"))
     
     if (v3[2] > 0.125) {
       text_x <- c(text_x, (numBases - 1) * 1.25 + 0.5 * (0.51 + 0.76))
       text_y <- c(text_y, 2 + 0.5 * v3[2])
       text_lab <- c(text_lab, "-")
       text_col <- c(text_col, "w")
     }
     
  }
  
  # arrow
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + (centerBase - 1) * 1.25
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1
  vs <- rep("arrow", length(xs))
  
  arrow_poly <- data.frame(x = xs, y = ys, v = vs)
  rect_data <- data.frame(x_start = x_start, x_end = x_end, y_start = y_start, y_end = y_end, rectType = rectType)
  text_data <- data.frame(x = text_x, y = text_y, label = text_lab, text_col = text_col)
  
  ggplot() + 
    geom_rect(data = rect_data, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end, fill = rectType)) + 
    geom_polygon(data = arrow_poly, aes(x = x, y = y, fill = v))  +
    geom_text(data = text_data, aes(label = label, x = x, y = y, colour = text_col), size = charSize) + 
    scale_colour_manual(values = c("#FFFFFF")) + 
    scale_fill_manual(
      values = c("A" = baseCol[1], "C" = baseCol[2], "G" = baseCol[3], 
                 "T" = baseCol[4], "+" = baseCol[5], "-" = baseCol[6], arrow = "#A8A8A8")) + 
    guides(fill=FALSE) +
    guides(colour=FALSE) +
    guides(size=FALSE) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          axis.title = element_blank())
  
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
  
  flankingPatternNum <- 4^(numBases - 1)
  subPattern <- c(rep("C>A", flankingPatternNum), 
                  rep("C>G", flankingPatternNum), 
                  rep("C>T", flankingPatternNum), 
                  rep("T>A", flankingPatternNum), 
                  rep("T>C", flankingPatternNum), 
                  rep("T>G", flankingPatternNum)
                  )
  
  X <- data.frame(probability = Fvec)
  # X$strand <- factor(rep(c("plus", "minus"), 1536), levels=c("plus", "minus"))
  
  if (trDir == TRUE) {
    X$strand <- factor(c(rep("plus", flankingPatternNum * 6), rep("minus", flankingPatternNum * 6)), levels=c("plus", "minus"))
    X$subtype <- factor(rep(subPattern, 2), levels=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
    X$flank <- rep(1:flankingPatternNum, 2)
  } else {
    X$subtype <- factor(subPattern, levels=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
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
    guides(fill=FALSE)
  
  if (trDir == TRUE) {
    gp <- gp + facet_grid(strand ~ subtype)
  } else {
    gp <- gp + facet_grid(. ~ subtype)
  }
  
  gp
  
}

  
setMethod("getMutNum", 
          signature = c(object = "MutationFeatureData"), 
          function(object) {         
            mutData <- data.frame(type = object@countData[1,], sampleName = object@sampleList[object@countData[2,]], count = object@countData[3,])
            sample2mutNum <- summarize(group_by(mutData, sampleName), mutationNum = sum(count))
            return(as.data.frame(sample2mutNum))
          }
)


setMethod("visMembership",
          signature = c(object1 = "MutationFeatureData", object2 = "EstimatedParameters"),
          function(object1, object2, ylog = FALSE, sortSampleNum = TRUE, multiplySampleNum = TRUE, fromSample = NULL, toSample = NULL, reorderSig = NULL, colourBrewer = NULL) {
          
            snum <- getMutNum(object1)
            if (ylog == TRUE) {
              snum[,2] <- log10(snum[,2])
            }  
            
            sampleList <- object2@sampleList
            signatureNum <- object2@signatureNum
            Q <- as.data.frame(object2@sampleSignatureDistribution)

            if (is.null(fromSample)) {
              fromSample <- 1
            }
            if (is.null(toSample)) {
              toSample <- length(sampleList)
            }
            
            if (is.null(reorderSig)) {
              reorderSig <- as.character(1:signatureNum)
              if (object2@isBackGround) {
                reorderSig[signatureNum] <- "BG"
              }
            }
            sigOrder <- reorderSig
            
            
            vMutationNum <- c()
            vSample <- c()
            vSignature <- c()

            if (sortSampleNum == TRUE) {
              mutNumOrder <- order(snum$mutationNum, decreasing = TRUE)[fromSample:toSample]
            } else {
              mutNumOrder <- fromSample:toSample
            }

            for (k in 1:signatureNum) {
              vSample <- c(vSample, sampleList[mutNumOrder])
              vSignature <- c(vSignature ,rep(sigOrder[k], length(mutNumOrder)))
              if (multiplySampleNum == TRUE) {
                vMutationNum <- c(vMutationNum, snum$mutationNum[mutNumOrder] * Q[mutNumOrder,k])
              } else {
                vMutationNum <- c(vMutationNum, Q[mutNumOrder,k]);                
              }
            }
            vSample <- factor(vSample, levels = sampleList[mutNumOrder])
            
            membership <- data.frame(sample = vSample, signature = as.factor(vSignature), mutationNum = vMutationNum)
            # membership <- data.frame(sample = reorder(vSample, -vIntensity), signature = as.factor(vSignature), intensity = vIntensity)
            
            
            gg <- ggplot(membership, aes(x = sample, y = mutationNum, fill = signature)) +
              geom_bar(width = 0.8, stat = "identity") +
              theme_bw() +
              theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
            if (multiplySampleNum == TRUE) {
              if (ylog == TRUE) {
                gg <- gg + ylab("log10(#mutation)")
              } else {
                gg <- gg + ylab("#mutation")
              }
            } else {
                gg <- gg + ylab("#membershipRatio");             
            }
            if (!is.null(colourBrewer)) {
              gg <- gg + scale_fill_brewer(palette = colourBrewer)
            }
            
            gg
          }

)


