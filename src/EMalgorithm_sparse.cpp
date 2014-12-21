#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector updateTheta_NormalizedSparseC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground, NumericVector vF0) {

  NumericVector vTheta(signatureNum * samplePatternNum);
  NumericVector vF_full(signatureNum * patternNum);

  
  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
    for (int m = 0; m < patternNum; m++) {
      vF_full[signatureNum - 1 + m * signatureNum] = vF0[m];
    }
  } else {
    variableSigNum = signatureNum;
  }

  
  for (int m = 0; m < patternNum; m++) {
    for (int k = 0; k < variableSigNum; k++) {
      vF_full[k + m * signatureNum] = 1;
    }
  }
    
    
  int featureInd;
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      featureInd = vPatternList[l + m * fdim.size()] - 1;
      for (int k = 0; k < variableSigNum; k++) {
        vF_full[k + m * signatureNum] = vF_full[k + m * signatureNum] * vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()];
      }
    }
    
  }
  

  double tempSum = 0;
  double invTempSum = 0;
  int sampleInd;
  int patternInd;
  

  NumericVector nTheta(signatureNum * samplePatternNum);
  for (int nm = 0; nm < samplePatternNum; nm++) {
 
    patternInd = vSparseCount[0 + nm * 3] - 1;
    sampleInd = vSparseCount[1 + nm * 3] - 1;
      
    tempSum = 0;
    for(int k = 0; k < signatureNum; k++) {
      vTheta[k + nm * signatureNum] = vQ[k + sampleInd * signatureNum] * vF_full[k + patternInd * signatureNum];
      tempSum = tempSum + vTheta[k + nm * signatureNum];
    }
  
    if (tempSum > 1e-10) {
      invTempSum = 1 / tempSum;
      for(int k = 0; k < signatureNum; k++) {
        nTheta[k + nm * signatureNum] = vTheta[k + nm * signatureNum] * invTempSum;
      }
    } else {
      for(int k = 0; k < signatureNum; k++) {
        nTheta[k + nm * signatureNum] = 1 / signatureNum;
      }
    }
       
  }

  return(nTheta);


}


// [[Rcpp::export]]
double getLogLikelihoodSparseC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground, NumericVector vF0) {

  NumericVector vTheta(signatureNum * samplePatternNum);
  NumericVector vF_full(signatureNum * patternNum);

  
  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
    for (int m = 0; m < patternNum; m++) {
      vF_full[signatureNum - 1 + m * signatureNum] = vF0[m];
    }
  } else {
    variableSigNum = signatureNum;
  }

  
  for (int m = 0; m < patternNum; m++) {
    for (int k = 0; k < variableSigNum; k++) {
      vF_full[k + m * signatureNum] = 1;
    }
  }
    
    
  int featureInd;
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      featureInd = vPatternList[l + m * fdim.size()] - 1;
      for (int k = 0; k < variableSigNum; k++) {
        vF_full[k + m * signatureNum] = vF_full[k + m * signatureNum] * vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()];
      }
    }
    
  }
  

  double tempSum = 0;
  double invTempSum = 0;
  int sampleInd;
  int patternInd;
  
  double logLikelihood = 0;
  for (int nm = 0; nm < samplePatternNum; nm++) {
 
    patternInd = vSparseCount[0 + nm * 3] - 1;
    sampleInd = vSparseCount[1 + nm * 3] - 1;
      
    tempSum = 0;
    for(int k = 0; k < signatureNum; k++) {
      tempSum = tempSum + vQ[k + sampleInd * signatureNum] * vF_full[k + patternInd * signatureNum];
 
    }
  
    /* 
    if (nm == 1824) {
      std::cout << nm << "\t" << tempSum << "\t" << vSparseCount[2 + nm * 3] << "\t" << logLikelihood << "\n";
    }
    */ 
    
    logLikelihood = logLikelihood + log(tempSum) * vSparseCount[2 + nm * 3];
       
  }

  return(logLikelihood);
  
}



// [[Rcpp::export]]
NumericVector updateMstepFSparseC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector nTheta, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground) {

  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
  } else {
    variableSigNum = signatureNum;
  }
  
  NumericVector vF(variableSigNum * fdim.size() * max(fdim)); 
  IntegerVector currentDigits(fdim.size());
  
  NumericVector vG_Theta_sum(signatureNum * patternNum);
  int patternInd;
  for (int nm = 0; nm < samplePatternNum; nm++) {
    
      patternInd = vSparseCount[0 + nm * 3] - 1;
      
      for (int k = 0; k < signatureNum; k++) {
        vG_Theta_sum[k + patternInd * signatureNum] = vG_Theta_sum[k + patternInd * signatureNum] + vSparseCount[2 + nm * 3] * nTheta[k + nm * signatureNum];
      }
      
  }
  

  int featureInd;
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      featureInd = vPatternList[l + m * fdim.size()] - 1;
      for (int k = 0; k < variableSigNum; k++) {
        vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] = vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] + vG_Theta_sum[k + m * signatureNum];
      }
    } 
    
  }
    
  double tempSum;
  double invTempSum;
  for (int k = 0; k < variableSigNum; k++) {
    
    for (int l = 0; l < fdim.size(); l++) {
    
      tempSum = 0;
      for (int ll = 0; ll < max(fdim); ll++) {
        tempSum = tempSum + vF[k + l * variableSigNum + ll * variableSigNum * fdim.size()];
      }
  
      invTempSum = 1 / tempSum;
      for (int ll = 0; ll < max(fdim); ll++) {
        vF[k + l * variableSigNum + ll * variableSigNum * fdim.size()] = vF[k + l * variableSigNum + ll * variableSigNum * fdim.size()] * invTempSum;
      }
    
    } 
      
  }
  
  return(vF);
  
}



// [[Rcpp::export]]
NumericVector updateMstepQSparseC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector nTheta, int signatureNum, int sampleNum, int patternNum, int samplePatternNum) {

  NumericVector vQ(signatureNum * sampleNum); 
  
  int sampleInd;
  for (int nm = 0; nm < samplePatternNum; nm++) {
    
    sampleInd = vSparseCount[1 + nm * 3] - 1;
      
    for (int k = 0; k < signatureNum; k++) {
      vQ[k + sampleInd * signatureNum] = vQ[sampleInd + k * sampleNum] + vSparseCount[2 + nm * 3] * nTheta[k + nm * signatureNum];
    }
      
  }

  
  double tempSum;
  double invTempSum;
  for (int n = 0; n < sampleNum; n++) {
    
    tempSum = 0;
    for (int k = 0; k < signatureNum; k++) {
      tempSum = tempSum + vQ[k + n * signatureNum];
    }
  
    invTempSum = 1 / tempSum;
    for (int k = 0; k < signatureNum; k++) {
      vQ[k + n * signatureNum] = vQ[k + n * signatureNum] * invTempSum;
    }
    
  } 
  return(vQ);
  
}

