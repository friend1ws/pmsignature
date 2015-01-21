#' An S4 class to represent a mutation meta information common to many data types
#'
#'  @slot type type of data format (independent, full, custom)
#'  @slot flankingBasesNum the number of flanking bases to consider (only apprecable for independent and full types)
#'  @slot transcriptionDirection the flag representing whether transcription direction is considered or not
#'  @slot possibleFeatures a vector representing the numbers of possible values for each mutation feature
setClass("MetaInformation",
         representation = representation(
           "VIRTUAL",
           type = "character",
           flankingBasesNum = "integer",
           transcriptionDirection = "logical",
           possibleFeatures = "integer"
         )
)
           
#' An S4 class representing the mutation data
#' 
#' @slot featureVectorList a list of feature vectors actually observed in the input mutation data
#' @slot sampleList a list of sample names observed in the input mutation data
#' @slot countData a matrix representing the number of mutations and samples.
#'  The (first, second, third) colums are for (mutation pattern indice, sample indice, frequencies).
setClass(
  Class = "MutationFeatureData",
  contains = "MetaInformation",
  representation = representation(
    featureVectorList = "matrix",
    sampleList = "character",
    countData = "matrix"
    ),
  validity = function(object) {
    errors <- character();
    # check for the consistency about the possible feature vector.
    for (i in 1:length(object@possibleFeatures)) {
      if (any(!(object@featureVectorList[i] %in% seq_len(object@possibleFeatures[i])))) {
        errors <- c(errors, paste("Inconsistency in the ", i, "-th feature", sep = ""));
      }
    }
    # check for the number of mutation patterns
    if (any(!sort(unique(object@countData[1,])) %in% seq_len(ncol(object@featureVectorList)))) {
      errors <- c(errors, "Inconsistency about the number of mutation patterns");
    } 
    # check for the number of samples
    if (any(!sort(unique(object@countData[2,])) %in% seq_len(length(object@sampleList)))) {
      errors <- c(errors, "Inconsistency about the sample indice");
    }
    # check for the zero count
    if (any(!object@countData[3,] != 0)) {
      errors <- c(errors, "The count data should not include 0");      
    }
    if (length(errors) == 0) TRUE else errors
  }
)

#' An S4 class representing the estimated parameters
#' 
#' @slot sampleList a list of sample names observed in the input mutation data
#' @slot SignatureNum the number of mutation signatures specified at the time of estimation
#' @slot isBackGround the flag showing whether the background signature data is used or not.
#' @slot signatureFeatureDistribution estimated parameters for mutation signatures
#' @slot sampleSignatureDistribution estimated parameters for memberships of mutation signatures for each sample
setClass(
  Class = "EstimatedParameters",
  contains = "MetaInformation",
  representation = representation(
    sampleList = "character",
    SignatureNum = "integer",
    isBackGround = "logical",
    signatureFeatureDistribution = "array",
    sampleSignatureDistribution = "matrix"
  ),
  validity = function(object) {
    errors <- character();
    variantSignatureNum = object@SignatureNum - as.integer(object@isBackGround);
    # check for the estimated signature feature distribution
    if (any(object@signatureFeatureDistribution < 0) || any(object@signatureFeatureDistribution > 1)) {
      errors <- c(errors, "The estimated signature feature distribution value should be between 0 to 1");    
    }
    for (k in 1:variantSignatureNum) {
      for (l in 1:length(object@possibleFeatures)) {
        if (abs(sum(object@signatureFeatureDistribution[k, l, 1:object@possibleFeatures[l]]) - 1) > 1e-10) {
          errors <- c(errors, paste("The estimated values should sum to 1 at the ", k, "-th signature and the ", l, "-th feature", sep=""));
        }
      }
    }
    # check for estimated sample signature distribution
    if (any(object@sampleSignatureDistribution < 0) || any(object@sampleSignatureDistribution > 1)) {
      errors <- c(errors, "The estimated signature feature distribution value should be between 0 to 1");    
    }    
    if (nrow(object@sampleSignatureDistribution) != length(object@sampleList)) {
      errors <- c(errors, "Inconsistency in the number of samples and the estimated sample signature distibution");
    }
    if (ncol(object@sampleSignatureDistribution) != object@SignatureNum) {
      errors <- c(errors, "Inconsistency in the number of signatures and the estimated sample signature distibution");
    }
    for (n in 1:nrow(object@sampleSignatureDistribution)) {
      if (abs(sum(object@sampleSignatureDistibutrion[n, ]) - 1) > 1e-10) {
        errors <- c(errors, paste("The estimated values should sum to 1 at the ", n, "-th sample", sep=""));
      }
    }      
        
    if (length(errors) == 0) TRUE else errors
  }
)


