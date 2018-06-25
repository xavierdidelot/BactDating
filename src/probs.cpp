#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double coalpriorC(NumericVector leaves, NumericVector nodes, double alpha) {
  int n=leaves.length();
  //NumericVector nodes = clone(intnodes);
  //std::sort(nodes.begin(), nodes.end(), std::greater<double>());
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
double likelihoodGammaC(NumericMatrix tab, double mu) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    p+=R::dgamma(unrec*tab(i,1),unrec*mu*(tab(i,2)-tab(tab(i,3)-1,2)),1,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodRelaxedgammaC(NumericMatrix tab, double mu, double sigma) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double l=0;
  double unrec=1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    l=tab(i,2)-tab(tab(i,3)-1,2);
    if (tab.ncol()==5) {unrec=tab(i,4);l=l*unrec;}
    double ratevar=sigma*sigma;
    p+=R::dgamma(unrec*tab(i,1),l*mu*mu/(mu+l*ratevar),1.0+l*ratevar/mu,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodPoissonC(NumericMatrix tab, double mu) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    p+=R::dpois(round(unrec*tab(i,1)),unrec*mu*(tab(i,2)-tab(tab(i,3)-1,2)),1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodNegbinC(NumericMatrix tab, double mu, double sigma) {
  double k=mu*mu/sigma/sigma;
  double theta=sigma*sigma/mu;
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  double lengths;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    lengths=unrec*(tab(i,2)-tab(tab(i,3)-1,2));
    p+=R::dnbinom(round(unrec*tab(i,1)),k,1/(1.0+theta*lengths),1);
  }
  return(p);
}

// [[Rcpp::export]]
void changeinorderedvec(NumericVector vec,double old,double n) {
  //NumericVector res = clone(vec);
  int i=0;
  while (vec(i)!=old) i++;
  vec(i)=n;
  while (1) {
    if (i>0               &&vec(i-1)<vec(i)) {vec(i)=vec(i-1);vec(i-1)=n;i--;continue;}
    if (i<(vec.length()-1)&&vec(i+1)>vec(i)) {vec(i)=vec(i+1);vec(i+1)=n;i++;continue;}
    break;
  }
  //  return(res);
}
