#include <Rcpp.h>
using namespace Rcpp;


//' Convert the parameter Q so that turboEM can treat
//' 
//' @param vQ Q (converted to a vector)
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
// [[Rcpp::export]]
NumericVector convertToTurbo_Q(NumericVector vQ, int signatureNum, int sampleNum) {
  
  NumericVector turboQ(sampleNum * (signatureNum - 1)); 
    
  for (int n = 0; n < sampleNum; n++) {

    for (int k = 0; k < (signatureNum - 1); k++) {
      turboQ[n + k * sampleNum] = vQ[n + (k + 1) * sampleNum];
    }
    
  } 
  
  return(turboQ);

}


//' Convert the parameter F so that turboEM can treat
//' 
//' @param vF F (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
// [[Rcpp::export]]
NumericVector convertToTurbo_F(NumericVector vF, NumericVector fdim, int signatureNum, bool isBackground) {
  
  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
  } else {
    variableSigNum = signatureNum;
  }
  
  NumericVector turboF(variableSigNum * (sum(fdim) - fdim.size())); 
  int cumFdim;
  for (int k = 0; k < variableSigNum; k++) {

    cumFdim = 0;
    for (int l = 0; l < fdim.size(); l++) {
      for (int ll = 0; ll < (fdim[l] - 1); ll++) {
        turboF[k + (ll + cumFdim) * variableSigNum] = vF[k + l * variableSigNum + (ll + 1) * variableSigNum * fdim.size()];
      }
      cumFdim = cumFdim + fdim[l] - 1;
    }

  } 
  
  return(turboF);

}


//' Restore the converted parameter Q for turboEM
//' 
//' @param turboQ Q (converted for turboEM)
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
// [[Rcpp::export]]
NumericVector convertFromTurbo_Q(NumericVector turboQ, int signatureNum, int sampleNum) {
  
  NumericVector rQ(sampleNum * signatureNum); 
  double tempSum;
  for (int n = 0; n < sampleNum; n++) {

    tempSum = 0;
    for (int k = 1; k < signatureNum; k++) {
      rQ[n + k * sampleNum] = turboQ[n + (k - 1) * sampleNum];
      tempSum = tempSum + rQ[n + k * sampleNum];
    }
    
    if (1 - tempSum < 0) {
      rQ[n] = 0;
    } else {
      rQ[n] = 1 - tempSum;
    }
 
  } 
  
  return(rQ);

}


//' Restore the converted parameter F for turboEM
//' 
//' @param turboF F (converted for turboEM)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
//' @param isBackground the logical value showing whether a background mutaiton features is included or not
// [[Rcpp::export]]
NumericVector convertFromTurbo_F(NumericVector turboF, NumericVector fdim, int signatureNum, bool isBackground) {
  
  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
  } else {
    variableSigNum = signatureNum;
  }
  
  NumericVector rF(variableSigNum * fdim.size() * max(fdim));
  double tempSum;
  int cumFdim;
  
  for (int k = 0; k < variableSigNum; k++) {
    
    cumFdim = 0;
    for (int l = 0; l < fdim.size(); l++) {
      
      tempSum = 0;
      for (int ll = 1; ll < fdim[l]; ll++) {
        rF[k + l * variableSigNum + ll * variableSigNum * fdim.size()] = turboF[k + (cumFdim + ll - 1) * variableSigNum];
        tempSum = tempSum + rF[k + l * variableSigNum + ll * variableSigNum * fdim.size()];
      }
      
      if (1 - tempSum < 0) {
        rF[k + l * variableSigNum + 0 * variableSigNum * fdim.size()] = 0;
      } else {
        rF[k + l * variableSigNum + 0 * variableSigNum * fdim.size()] = 1 - tempSum;
      }
      
      cumFdim = cumFdim + fdim[l] - 1;
      
    }
    
  }
  
  return(rF);
  
}





