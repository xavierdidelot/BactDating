#' Poisson likelihood function
#' @param tab Table of nodes
#' @param mu Substitution clock rate
#' @return log-likelihood
likelihoodPoisson = function(tab, mu) {
  n = ceiling(nrow(tab)/2)
  if (mu < 0) stop('error1')
  t2 = tab[-(n + 1), ,drop=F]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  #return(sum(-lengths * mu + round(muts) * log(lengths * mu)))
  return(sum(dpois(round(muts),lengths*mu,log=T)))
}

#' Negative-binomial likelihood function
#' @param tab Table of nodes
#' @param r r parameter
#' @param phi phi parameter
#' @return log-likelihood
likelihoodNegbin = function(tab, r, phi) {
  n = ceiling(nrow(tab)/2)
  if (r < 0 || phi < 0) stop('error1')
  t2 = tab[-(n + 1), ,drop=F]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  return(sum(dnbinom(round(muts),r,1-phi*lengths/(1+phi*lengths),log=T)))
}

#' Gamma likelihood function
#' @param tab Table of nodes
#' @param mu Clock rate parameter
#' @return log-likelihood
likelihoodGamma = function(tab, mu) {
  n = ceiling(nrow(tab)/2)
  if (mu < 0) stop('error1')
  t2 = tab[-(n + 1), ,drop=F]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  return(sum(dgamma(muts,shape=mu*lengths,scale=1,log=T)))
}

#' Relaxed gamma likelihood function
#' @param tab Table of nodes
#' @param mu Clock rate parameter
#' @param sigma Std of per branch clock rate
#' @return log-likelihood
likelihoodRelaxedgamma = function(tab, mu, sigma) {
  n = ceiling(nrow(tab)/2)
  if (mu < 0) stop('error1')
  t2 = tab[-(n + 1), ,drop=F]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0) stop('error2')
  muts = t2[, 2]
  if (min(muts)<0) stop('error3')
  if (ncol(tab)==5) {lengths=lengths*t2[,5];muts=muts*t2[,5]}
  ratevar=sigma^2
  return(sum(dgamma(muts,shape=mu*mu*lengths/(mu+ratevar*lengths),scale=1+lengths*ratevar/mu,log=T)))
}

#' Coalescent prior function
#' @param tab Table of nodes
#' @param alpha Coalescent time unit
#' @return The log-prior in Eq (1) of Drummond et al (2002) Genetics
coalpriorR = function(tab, alpha) {
  n = ceiling(nrow(tab)/2)
  p = -log(alpha) * (n - 1)
  l=nrow(tab)
  s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
  k=cumsum(2*(s$ix<=n)-1)
  difs=s$x[1:(l-1)]-s$x[2:l]#faster than -diff(s$x)
  p=p-sum(k[1:(l-1)]*(k[1:(l-1)]-1)*difs)/(2*alpha)
  if (k[length(k)] != 1)
    stop('error')
  return(p)
}

#Provide equivalent R function to the C prior function
  coalprior = function(leaves, intnodes, alpha) {
  n = length(leaves)
  tab=matrix(0,n*2-1,4)
  tab[,3]=c(leaves,intnodes)
  coalpriorR(tab,alpha)
  }

