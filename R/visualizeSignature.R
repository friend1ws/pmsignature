#' @export
setGeneric("visPMSignature", function(object, K = 1, baseCol) {
  standardGeneric("visPMSignature")
})

setMethod("visPMSignature", 
          signature = c(object = "EstimatedParameters", K = "numeric", baseCol = "numeric"), 
          function(object, K = 1, baseCol) {

            vF <- object@signatureFeatureDistribution[K,,];
            numBases <- object@flankingBasesNum;
            centerBase <- (1 + numBases) / 2;
            
            visPMS_ind(vF, numBases, baseCol);
          }
)

setMethod("visPMSignature", 
          signature = c(object = "EstimatedParameters", K = "numeric"), 
          function(object, K = 1) {
            
            vF <- object@signatureFeatureDistribution[K,,];
            numBases <- object@flankingBasesNum;
            centerBase <- (1 + numBases) / 2;
            
            # soft colors
            visPMS_ind(vF, numBases, c(rgb(57, 168, 105, max = 256), rgb(71, 132, 191, max = 256), rgb(242, 229, 92, max = 256), rgb(221, 102, 115, max = 256)));

            # humble colors
            # visPMS_ind(vF, numBases, c(rgb(0, 148, 83, max = 256), rgb(19, 110, 171, max = 256), rgb(223, 210, 56, max = 256), rgb(202, 71, 92, max = 256)));
            
          }
)

#' @title visualize probabisitic mutaiton signature for the independent model
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and flanking bases represented by the indepenent
#'   representation.
#'   
#' @param vF a matrix for mutation signature
#' @param numBases the number of flanking bases 
#' @param centerBases the number of flanking bases
visPMS_ind <- function(vF, numBases, baseCol) {
  
  centerBase <- (1 + numBases) / 2;
  
  v1 <- vF[1,1:6];
  V2 <- vF[2:(numBases),1:4];
  A <- matrix(0, numBases, 4);
  B <- matrix(0, 4, 4);
  
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
  
  startx <- 0;
  for(l in 1:numBases) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col = baseCol[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=1.2)
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
        text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col="white", cex=1.2)
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
  
}


vispMS_full <- function(Fvec, numBases, trDir) {
  
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
  
  # X$flank <- rep(1:flankingPatternNum, 2);              
  # X$flank <- rep(as.vector(rbind(1:flankingPatternNum, 1:flankingPatternNum)), 6);
  
  gp <- ggplot(X, aes(x=flank, y=probability, fill=subtype)) +
    geom_bar(stat="identity", position="identity") + 
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


