#' Poisson likelihood function
#' @param tab Table of nodes
#' @param rate Substitution rate
#' @return log-likelihood
likelihoodPoisson = function(tab, rate) {
  n = ceiling(nrow(tab)/2)
  if (rate < 0) stop('error1')
  t2 = tab[-(n + 1), ]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  return(sum(-lengths * rate + muts * log(lengths * rate)))
}

#' Negative-binomial likelihood function
#' @param tab Table of nodes
#' @param r r parameter
#' @param phi phi parameter
#' @return log-likelihood
likelihoodNegbin = function(tab, r, phi) {
  n = ceiling(nrow(tab)/2)
  if (r < 0 || phi < 0) stop('error1')
  t2 = tab[-(n + 1), ]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  return(sum(dnbinom(muts,r,1-phi*lengths/(1+phi*lengths),log=T)))
}

#' Gamma likelihood function
#' @param tab Table of nodes
#' @param rate rate parameter
#' @return log-likelihood
likelihoodGamma = function(tab, rate) {
  n = ceiling(nrow(tab)/2)
  if (rate < 0) stop('error1')
  t2 = tab[-(n + 1), ]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  return(sum(dgamma(muts,shape=rate*lengths,scale=1,log=T)))
}

#' Coalescent prior function
#' @param tab Table of nodes
#' @param neg Coalescent rate
#' @return The log-prior in Eq (1) of Drummond et al (2002) Genetics
coalpriorR = function(tab, neg) {
  n = ceiling(nrow(tab)/2)
  p = -log(neg) * (n - 1)
  l=nrow(tab)
  s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
  k=cumsum(2*(s$ix<=n)-1)
  difs=s$x[1:(l-1)]-s$x[2:l]#faster than -diff(s$x)
  p=p-sum(k[1:(l-1)]*(k[1:(l-1)]-1)*difs)/(2*neg)
  if (k[length(k)] != 1)
    stop('error')
  return(p)
}

#Provide equivalent R function to the C prior function
  coalprior = function(leaves, intnodes, neg) {
  n = length(leaves)
  tab=matrix(0,n*2-1,4)
  tab[,3]=c(leaves,intnodes)
  coalpriorR(tab,neg)
  }

