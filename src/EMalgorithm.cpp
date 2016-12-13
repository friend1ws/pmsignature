#include <Rcpp.h>
using namespace Rcpp;


//' Update the auxiliary parameters theta and normalize them so that the summation of each group sums to 1 (E-step),
//' also calculate the current log-likelihood value
//' 
//' @param vPatternList The list of possible mutation features (converted to a vector)
//' @param vSparseCount The table showing (mutation feature, sample, the number of mutation) (converted to a vector)
//' @param vF F (converted to a vector)
//' @param vQ Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param samplePatternNum the number of possible combination of samples and mutation patternns
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @param vF0 a background mutaiton features
// [[Rcpp::export]]
NumericVector updateTheta_NormalizedC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground, NumericVector vF0) {

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
        /*
        if (vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] < 0) {
          std::cout << l << "\t" << k << "\t" << featureInd << "\t" << vF[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] << "\n";
        }
        */
      }
    }
    
  }
  

  double tempSum = 0;
  double invTempSum = 0;
  int sampleInd;
  int patternInd;
  // double tempQsum = 0;

  NumericVector nTheta(signatureNum * samplePatternNum);
  for (int nm = 0; nm < samplePatternNum; nm++) {
 
    patternInd = vSparseCount[0 + nm * 3] - 1;
    sampleInd = vSparseCount[1 + nm * 3] - 1;
      
    tempSum = 0;
    // tempQsum = 0;
    for(int k = 0; k < signatureNum; k++) {
      vTheta[k + nm * signatureNum] = vQ[k + sampleInd * signatureNum] * vF_full[k + patternInd * signatureNum];
      tempSum = tempSum + vTheta[k + nm * signatureNum];
      // tempQsum = tempQsum + vQ[k + sampleInd * signatureNum];
      
    }

    /*
    if (tempQsum < 0.99 or tempQsum > 1.01) {
      std::cout << nm << "\t" << tempQsum << "\t" << "\n";
    }
    */
    
    
    
    if (tempSum > 1e-8) {
      invTempSum = 1 / tempSum;
      for(int k = 0; k < signatureNum; k++) {
        nTheta[k + nm * signatureNum] = vTheta[k + nm * signatureNum] * invTempSum;
      }
    } else {
      for(int k = 0; k < signatureNum; k++) {
        nTheta[k + nm * signatureNum] = 1.0 / signatureNum;
        // std::cout << nm << "\t" << k << "\t" << signatureNum << "\t" << nTheta[k + nm * signatureNum] << "\t" << tempSum << "\t" << "\n";
      }
    }
       
  }

  return(nTheta);


}


//' Calculate the value of the log-likelihood for given parameters
//' 
//' @param vPatternList The list of possible mutation features (converted to a vector)
//' @param vSparseCount The table showing (mutation feature, sample, the number of mutation) (converted to a vector)
//' @param vF F (converted to a vector)
//' @param vQ Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param samplePatternNum the number of possible combination of samples and mutation patternns
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @param vF0 a background mutaiton features
// [[Rcpp::export]]
double getLogLikelihoodC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground, NumericVector vF0) {

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
    if (tempSum > 1e-10 and vSparseCount[2 + nm * 3] > 0) {
      logLikelihood = logLikelihood + log(tempSum) * vSparseCount[2 + nm * 3];
    }
       
  }

  return(logLikelihood);
  
}


//' Update the parameter F and Q (M-step in the EM-algorithm)
//' 
//' @param vPatternList The list of possible mutation features (converted to a vector)
//' @param vSparseCount The table showing (mutation feature, sample, the number of mutation) (converted to a vector)
//' @param vF F (converted to a vector)
//' @param vQ Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param samplePatternNum the number of possible combination of samples and mutation patternns
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @param vF0 a background mutaiton features
// [[Rcpp::export]]
NumericVector updateMstepFQC(NumericVector vPatternList, NumericVector vSparseCount, NumericVector nTheta, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, int samplePatternNum, bool isBackground) {

  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
  } else {
    variableSigNum = signatureNum;
  }
  
  int Fsize = variableSigNum * fdim.size() * max(fdim);
  int Qsize = signatureNum * sampleNum;
  NumericVector vF_Q(Fsize + Qsize); 
  IntegerVector currentDigits(fdim.size());
  
  NumericVector vG_Theta_sum_F(signatureNum * patternNum);
  NumericVector vG_Theta_sum_Q(signatureNum * sampleNum);
  NumericVector vQ(signatureNum * sampleNum); 
  int patternInd;
  int sampleInd;
  double tempCountTheta;
  for (int nm = 0; nm < samplePatternNum; nm++) {
    
      patternInd = vSparseCount[0 + nm * 3] - 1;
      sampleInd = vSparseCount[1 + nm * 3] - 1;
      
      for (int k = 0; k < signatureNum; k++) {
        tempCountTheta = vSparseCount[2 + nm * 3] * nTheta[k + nm * signatureNum];
        vG_Theta_sum_F[k + patternInd * signatureNum] = vG_Theta_sum_F[k + patternInd * signatureNum] + tempCountTheta;
        vG_Theta_sum_Q[k + sampleInd * signatureNum] = vG_Theta_sum_Q[k + sampleInd * signatureNum] + tempCountTheta;
      }
      
  }
  

  int featureInd;
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      featureInd = vPatternList[l + m * fdim.size()] - 1;
      for (int k = 0; k < variableSigNum; k++) {
        vF_Q[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] = vF_Q[k + l * variableSigNum + featureInd * variableSigNum * fdim.size()] + vG_Theta_sum_F[k + m * signatureNum];
      }
    } 
    
  }
    
  double tempSum;
  double invTempSum;
  for (int k = 0; k < variableSigNum; k++) {
    
    for (int l = 0; l < fdim.size(); l++) {
    
      tempSum = 0;
      for (int ll = 0; ll < max(fdim); ll++) {
        tempSum = tempSum + vF_Q[k + l * variableSigNum + ll * variableSigNum * fdim.size()];
      }
  
      invTempSum = 1 / tempSum;
      for (int ll = 0; ll < max(fdim); ll++) {
        vF_Q[k + l * variableSigNum + ll * variableSigNum * fdim.size()] = vF_Q[k + l * variableSigNum + ll * variableSigNum * fdim.size()] * invTempSum;
      }
    
    } 
      
  }
  
  

  for (int n = 0; n < sampleNum; n++) {
    
    tempSum = 0;
    for (int k = 0; k < signatureNum; k++) {
      tempSum = tempSum + vG_Theta_sum_Q[k + n * signatureNum];
    }
  
    if (tempSum > 1e-8) {
      invTempSum = 1 / tempSum;
        for (int k = 0; k < signatureNum; k++) {
          vF_Q[Fsize + k + n * signatureNum] = vG_Theta_sum_Q[k + n * signatureNum] * invTempSum;
        }
    } else {
        invTempSum = 1 / signatureNum;
        for (int k = 0; k < signatureNum; k++) {
          vF_Q[Fsize + k + n * signatureNum] = invTempSum;
        }        
    }
  }
  return(vF_Q);
  
  
}

