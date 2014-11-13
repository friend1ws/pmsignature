#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

//' Update the auxiliary parameters theta and normalize them so that the summation of each group sums to 1 (E-step)
//' 
//' @param vF F (converted to a vector)
//' @param vQ Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @param vF0 a background mutaiton features
//' @export
// [[Rcpp::export]]
NumericVector updateTheta_NormalizedC(NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, bool isBackground, NumericVector vF0) {

  NumericVector vTheta(sampleNum * signatureNum * patternNum);
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
    
  IntegerVector currentDigits(fdim.size());
  for (int l = 0; l < fdim.size(); l++) {
    currentDigits[l] = 0;
  }
  
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      for (int k = 0; k < variableSigNum; k++) {
        vF_full[k + m * signatureNum] = vF_full[k + m * signatureNum] * vF[k + l * variableSigNum + currentDigits[l] * variableSigNum * fdim.size()];
      }
    }
          
    int tl = 0;
    while(tl < fdim.size() and currentDigits[tl] + 1 >= fdim[tl]) {
      currentDigits[tl] = 0;
      tl = tl + 1;
    }
      
    if (tl < fdim.size()) {
      currentDigits[tl] = currentDigits[tl] + 1;
    }
    
  }

  double tempSum = 0;
  double invTempSum = 0;
  NumericVector nTheta(sampleNum * signatureNum * patternNum);

  for (int m = 0; m < patternNum; m++) {
     for (int n = 0; n < sampleNum; n++) {
       
      tempSum = 0;
      for(int k = 0; k < signatureNum; k++) {
        vTheta[n + k * sampleNum + m * sampleNum * signatureNum] = vQ[n + k * sampleNum] * vF_full[k + m * signatureNum];
        tempSum = tempSum + vTheta[n + k * sampleNum + m * sampleNum * signatureNum];
      }
  
      if (tempSum > 1e-10) {
        invTempSum = 1 / tempSum;
        for(int k = 0; k < signatureNum; k++) {
          nTheta[n + k * sampleNum + m * sampleNum * signatureNum] = vTheta[n + k * sampleNum + m * sampleNum * signatureNum] * invTempSum;
        }
      } else {
        for(int k = 0; k < signatureNum; k++) {
          nTheta[n + k * sampleNum + m * sampleNum * signatureNum] = 1 / signatureNum;
        }
      }
        
    }
  }

  return(nTheta);

}


//' Calculate the value of the log-likelihood for given parameters
//' 
//' @param vG the matrix of mutation feature data
//' @param vF the parameter F (converted to a vector)
//' @param vQ the parameter Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @param vF0 a background mutaiton features
//' @export
// [[Rcpp::export]]
double getLogLikelihoodC(NumericVector vG, NumericVector vF, NumericVector vQ, NumericVector fdim, int signatureNum, int sampleNum, int patternNum, bool isBackground, NumericVector vF0) {

  NumericVector Theta(sampleNum * signatureNum * patternNum);
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
    
  IntegerVector currentDigits(fdim.size());
  for (int l = 0; l < fdim.size(); l++) {
    currentDigits[l] = 0;
  }
  
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      for (int k = 0; k < variableSigNum; k++) {
        vF_full[k + m * signatureNum] = vF_full[k + m * signatureNum] * vF[k + l * variableSigNum + currentDigits[l] * variableSigNum * fdim.size()];
      }
    }
          
    int tl = 0;
    while(tl < fdim.size() and currentDigits[tl] + 1 >= fdim[tl]) {
      currentDigits[tl] = 0;
      tl = tl + 1;
    }
      
    if (tl < fdim.size()) {
      currentDigits[tl] = currentDigits[tl] + 1;
    }
    
  }


  for(int k = 0; k < signatureNum; k++) {
    for (int m = 0; m < patternNum; m++) {
      for (int n = 0; n < sampleNum; n++) {
        Theta[n + k * sampleNum + m * sampleNum * signatureNum] = vQ[n + k * sampleNum] * vF_full[k + m * signatureNum];
      }
    }
  }


  double logLikelihood = 0;
  double tempSum = 0;

  for (int m = 0; m < patternNum; m++) {
    for (int n = 0; n < sampleNum; n++) {
      
      if (vG[n + m * sampleNum] > 0) {
        tempSum = 0;
        for(int k = 0; k < signatureNum; k++) {
          tempSum = tempSum + Theta[n + k * sampleNum + m * sampleNum * signatureNum];
        }
        // if (tempSum < 1e-10) {
        //   std::cout << m << "\t" << n << "\t" << vG[n + m * sampleNum] << "\t" << "\t" << tempSum << "\n";
        //   for (int k = 0; k < signatureNum; k++) {
        //     std::cout << vQ[n + k * sampleNum] << "\t" << vF_full[k + m * signatureNum] << "\n";
        //   }
        // }
        
        logLikelihood = logLikelihood + log(tempSum) * vG[n + m * sampleNum];
      }
      
    }
  }

  return(logLikelihood);
}


//' Update the parameter F (M-step in the EM-algorithm)
//' 
//' @param vTheta theparameter Theta (converted to a vector)
//' @param vG the matrix of mutation feature data
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
//' @export
// [[Rcpp::export]]
NumericVector updateMstepFC(NumericVector vTheta, NumericVector vG,NumericVector fdim, int signatureNum, int sampleNum, int patternNum, bool isBackground) {

  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
  } else {
    variableSigNum = signatureNum;
  }
  
  NumericVector vF(variableSigNum * fdim.size() * max(fdim)); 
  IntegerVector currentDigits(fdim.size());
  
  NumericVector vG_Theta_sum(signatureNum * patternNum);
  for (int m = 0; m < patternNum; m++) {
    for (int n = 0; n < sampleNum; n++) {
      for (int k = 0; k < signatureNum; k++) {
        vG_Theta_sum[k + m * signatureNum] = vG_Theta_sum[k + m * signatureNum] + vG[n + m * sampleNum] * vTheta[n + k * sampleNum + m * sampleNum * signatureNum];
      }
    }
  }
  
  for (int l = 0; l < fdim.size(); l++) {
    currentDigits[l] = 0;
  }

  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      for (int k = 0; k < variableSigNum; k++) {
        vF[k + l * variableSigNum + currentDigits[l] * variableSigNum * fdim.size()] = vF[k + l * variableSigNum + currentDigits[l] * variableSigNum * fdim.size()] + vG_Theta_sum[k + m * signatureNum];
      }
    }
      
    int tl = 0;
    while(tl < fdim.size() and currentDigits[tl] + 1 >= fdim[tl]) {
      currentDigits[tl] = 0;
      tl = tl + 1;
    }
      
    if (tl < fdim.size()) {
      currentDigits[tl] = currentDigits[tl] + 1;
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


//' Update the parameter Q (M-step in the EM-algorithm)
//' 
//' @param vTheta theparameter Theta (converted to a vector)
//' @param vG the matrix of mutation feature data
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation features
//' @export
// [[Rcpp::export]]
NumericVector updateMstepQC(NumericVector vTheta, NumericVector vG, int signatureNum, int sampleNum, int patternNum) {

  NumericVector vQ(sampleNum * signatureNum); 
  
  for (int k = 0; k < signatureNum; k++) {
    
    for (int m = 0; m < patternNum; m++) {
      
      for (int n = 0; n < sampleNum; n++) {
        vQ[n + k * sampleNum] = vQ[n + k * sampleNum] + vG[n + m * sampleNum] * vTheta[n + k * sampleNum + m * sampleNum * signatureNum];
      } 
      
    }

  }
  
  double tempSum;
  double invTempSum;
  for (int n = 0; n < sampleNum; n++) {
    
    tempSum = 0;
    for (int k = 0; k < signatureNum; k++) {
      tempSum = tempSum + vQ[n + k * sampleNum];
    }
  
    invTempSum = 1 / tempSum;
    for (int k = 0; k < signatureNum; k++) {
      vQ[n + k * sampleNum] = vQ[n + k * sampleNum] * invTempSum;
    }
    
  } 
  return(vQ);
  
}

