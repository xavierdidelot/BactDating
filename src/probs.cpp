#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double coalpriorC(NumericVector leaves, NumericVector intnodes, double neg) {
  int n=leaves.length();
  NumericVector nodes = clone(intnodes);
  std::sort(nodes.begin(), nodes.end(), std::greater<double>());
  double p = -log(neg) * (n - 1);
  int i1=1,i2=0,k=1;
  double prev=leaves[0];
  for (int i=1;i<(n+n-1);i++) {
    if (i1<n&&leaves[i1]>nodes[i2]) {
      p = p - (k * (k - 1.0) / (2.0 * neg) * (prev-leaves[i1]));
      prev=leaves[i1++];
      k++;
    } else {
      p = p - (k * (k - 1.0) / (2.0 * neg) * (prev-nodes[i2]));
      prev=nodes[i2++];
      k--;
    }
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodGammaC(NumericMatrix tab, double rate) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  for (int i=0;i<(n+n-1);i++) {
    if (i!=n) p+=R::dgamma(tab(i,1),rate*(tab(i,2)-tab(tab(i,3)-1,2)),1,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodPoissonC(NumericMatrix tab, double rate) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  for (int i=0;i<(n+n-1);i++) {
    if (i!=n) p+=R::dpois(tab(i,1),rate*(tab(i,2)-tab(tab(i,3)-1,2)),1);
  }
  return(p);
}

