setGeneric("visPMSignature", function(object, K, baseCol) {
  standardGeneric("visPMSignature")
})

setMethod("visPMSignature", signature(object = "EstimatedParameters"), 
          function(object, K = 1, baseCol = c(4, 7, 3, 1)) {

            vF <- object@signatureFeatureDistribution[K,,];
            numBases <- object@flankingBasesNum;
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
)

  
#' @title visualize probabisitic mutaiton signature for the independent model
#' with two 5' and 3' bases 
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and two 5' and 3' bases represented by the indepenent
#'   representation.
#'   
#' @param vF a matrix for mutation signature
#' @export
visPMS_ind5 <- function(vF = matrix(0, 5, 6)) {
  
  v1 <- vF[1,1:6];
  V2 <- vF[2:5,1:4];
  
  A <- matrix(0, 5, 4);
  B <- matrix(0, 4, 4);
  
  A[1,] <- V2[1,];
  A[2,] <- V2[2,];
  A[4,] <- V2[3,];
  A[5,] <- V2[4,];
  A[3,2] <- sum(v1[1:3]);
  A[3,4] <- sum(v1[4:6]);
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3]);
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6]);
  
  num2base <- c("A", "C", "G", "T");
  
  plot.window(xlim=c(-0.25, 6.5), ylim=c(-0.25, 3.25));
  # plot.window(xlim=c(0, 6.25), ylim=c(0, 3));
  
  # tcols <- c(4, 7, 3, 2);
  # tcols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99");
  # tcols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99"); # accent
  # tcols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"); # darks
  # tcols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"); # set 1
  tcols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"); # set2
  # tcols <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072"); # set3
  # tcols <- c("#d7191c", "#fdae61", "#abdda4", "#2b83ba"); # spectral
  
  startx <- 0;
  for(l in 1:5) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col=tcols[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=1.2)
      }
      startx <- endx;
    }
    startx <- startx + 0.25;
  }
  
  startx <- 2 * 1.25;
  for (w in 1:4) {
    starty <- 2;
    endx <- startx + A[3,w];
    for(ww in 1:4) {
      endy <- starty + B[w,ww];
      polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=tcols[ww], border=F);
      if ((endy - starty > 1 / 2) & (endx - startx > 1 / 2)) {
        text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col="white", cex=1.2)
      }
      starty <- endy;
    }
    startx <- endx
    starty <- endy;
  }
  
  ##########
  # draw arrow
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + 2.5;
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1;
  
  polygon(xs, ys, col=8, border=F);
  ##########
  
}

#' @title visualize probabisitic mutaiton signature for the independent model
#' with two 5' and 3' bases and transcription direction
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns, two 5' and 3' bases and transcription directions represented by the indepenent
#'   representation.
#'   
#' @param vF a matrix for mutation signature
#' @export
visPMS_ind5_dir <- function(vF = matrix(0, 6, 6)) {
  
  v1 <- vF[1,1:6];
  V2 <- vF[2:5,1:4];
  v3 <- vF[6,1:2];
  
  A <- matrix(0, 5, 4);
  B <- matrix(0, 4, 4);
  
  A[1,] <- V2[1,];
  A[2,] <- V2[2,];
  A[4,] <- V2[3,];
  A[5,] <- V2[4,];
  A[3,2] <- sum(v1[1:3]);
  A[3,4] <- sum(v1[4:6]);
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3]);
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6]);
  
  num2base <- c("A", "C", "G", "T");
  
  plot.window(xlim=c(-0.25, 6.5), ylim=c(-0.25, 3.25));
  
  tcols <- c(4, 7, 3, 2);
  
  startx <- 0;
  for(l in 1:5) {
    
    for(w in 1:4) {
      endx <- startx + A[l,w]
      polygon(c(startx, endx, endx, startx), c(0, 0, 1, 1), col=tcols[w], border=F);
      if (endx - startx > 1 / 4) {
        text(0.5 * (endx + startx), 0.5, num2base[w], col="white", cex=1.2)
      }
      startx <- endx;
    }
    startx <- startx + 0.25;
  }
  
  startx <- 2 * 1.25;
  for (w in 1:4) {
    starty <- 2;
    endx <- startx + A[3,w];
    for(ww in 1:4) {
      endy <- starty + B[w,ww];
      polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col=tcols[ww], border=F);
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
  xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + 2.5;
  ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1;
  
  polygon(xs, ys, col=8, border=F);
  ##########
  
  ##########
  # draw direction bias
  startx <- 4 * 1.25 + 0.5;
  endx <- 4 * 1.25 + 0.75;
  starty <- 1.9;
  endy <- starty + v3[1];
  polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col="magenta", border=F);
  
  if (endy - starty > 1 / 8) {
    text(5.625, 0.5 * (starty + endy), "+", col="white", cex=1.2)
  }
  starty <- endy;
  endy <- 2.9;
  polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col="cyan", border=F);
  if (endy - starty > 1 / 8) {
    text(5.625, 0.5 * (starty + endy), "-", col="white", cex=1.2)
  }  
  
  
}


#' @title visualize probabisitic mutaiton signature for the full model
#' with one 5' and 3' bases
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and one 5' and 3' bases represented by the full
#'   representation.
#'   
#' @param v1 a vector for mutation signature
#' @export
visPMS_full3 <- function(v1 = rep(1 / 96, 96)) {
  barplot(v1, col=c(rep(1, 16), rep(2, 16), rep(3, 16), rep(4, 16), rep(5, 16), rep(6, 16)));
}

#' @title visualize probabisitic mutaiton signature for the full model
#' with two 5' and 3' bases
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and two 5' and 3' bases represented by the full
#'   representation.
#'   
#' @param v1 a vector for mutation signature
#' @export
visPMS_full5 <- function(v1 = rep(1 / 1536, 1536)) {
  barplot(v1, col=c(rep(1, 256), rep(2, 256), rep(3, 256), rep(4, 256), rep(5, 256), rep(6, 256)));
}
