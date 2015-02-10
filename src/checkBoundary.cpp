#include <Rcpp.h>
using namespace Rcpp;


//' Check whether the parameter F is within the appropriate range
//' 
//' @param turboF F (converted for turboEM)
//' @param fdim a vector specifying the number of possible values for each mutation signature
//' @param signatureNum the number of mutation signatures
// [[Rcpp::export]]
bool boundaryTurbo_F(NumericVector turboF, NumericVector fdim, int signatureNum) {
  
  bool isBoundary =TRUE;
  double tempSum;
  int cumFdim;
  
  for (int k = 0; k < signatureNum; k++) {
    
    cumFdim = 0;
    for (int l = 0; l < fdim.size(); l++) {
      
      tempSum = 0;
      for (int ll = 1; ll < fdim[l]; ll++) {
        if (turboF[k + (cumFdim + ll - 1) * signatureNum] < 0) {
          isBoundary = FALSE;
          // std::cout << k << "\t" << l << "\t" << ll << "\n";
        }
        tempSum = tempSum + turboF[k + (cumFdim + ll - 1) * signatureNum];
      }
      
      if (1 - tempSum < 0) {
        isBoundary = FALSE;
      }
      
      cumFdim = cumFdim + fdim[l] - 1;
      
    }
    
  }
  
  return(isBoundary);
  
}


//' Check whether the parameter Q is within the appropriate range
//' 
//' @param turboQ Q (converted for turboEM)
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
// [[Rcpp::export]]
bool boundaryTurbo_Q(NumericVector turboQ, int signatureNum, int sampleNum) {
  
  bool isBoundary = TRUE; 
  double tempSum;
  for (int n = 0; n < sampleNum; n++) {

    tempSum = 0;
    for (int k = 1; k < signatureNum; k++) {
      if (turboQ[n + (k - 1) * sampleNum] < 0) {
        isBoundary = FALSE;
      }
      tempSum = tempSum + turboQ[n + (k - 1) * sampleNum];
    }
    
    if (1 - tempSum < 0) {
      isBoundary = FALSE;
    }
 
  } 
  
  return(isBoundary);

}

