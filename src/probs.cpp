#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double coalpriorC(NumericVector leaves, NumericVector intnodes, double alpha) {
  int n=leaves.length();
  NumericVector nodes = clone(intnodes);
  std::sort(nodes.begin(), nodes.end(), std::greater<double>());
  double p = -log(alpha) * (n - 1);
  int i1=1,i2=0,k=1;
  double prev=leaves[0];
  for (int i=1;i<(n+n-1);i++) {
    if (i1<n&&leaves[i1]>nodes[i2]) {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-leaves[i1]));
      prev=leaves[i1++];
      k++;
    } else {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-nodes[i2]));
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
  double unrec=1;
  for (int i=0;i<(n+n-1);i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    p+=R::dgamma(unrec*tab(i,1),unrec*rate*(tab(i,2)-tab(tab(i,3)-1,2)),1,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodRelaxedgammaC(NumericMatrix tab, double rate, double ratevar) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double l=0;
  double unrec=1;
  for (int i=0;i<(n+n-1);i++) {
    if (i==n) continue;
    l=tab(i,2)-tab(tab(i,3)-1,2);
    if (tab.ncol()==5) {unrec=tab(i,4);l=l*unrec;}
    p+=R::dgamma(unrec*tab(i,1),l*rate*rate/(rate+l*ratevar),1.0+l*ratevar/rate,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodPoissonC(NumericMatrix tab, double rate) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  for (int i=0;i<(n+n-1);i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    p+=R::dpois(round(unrec*tab(i,1)),unrec*rate*(tab(i,2)-tab(tab(i,3)-1,2)),1);
  }
  return(p);
}

