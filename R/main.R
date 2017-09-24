#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param rate Rate of substitutions per genome (not per site)
#' @param nbIts Number of MCMC iterations to perform
#' @param usePrior Whether or not to use a coalescent prior on the tree
#' @param updateRate Whether or not to update the substitution rate
#' @return Dating results
#' @export
credating = function(tree, date, rate = 1, nbIts = 1000, useCoalPrior = F, updateRate = F)
{
  prior=function(tab,neg,n) return(0)
  if (useCoalPrior) prior=coalprior
  n = length(tree$tip.label)

  #Create table of nodes (col1=name,col2=expected mutation on branch above,col3=date,col4=father)
  tab = matrix(NA, n + tree$Nnode, 4)
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    tab[i, 1] = i
    tab[i, 4] = tree$edge[r, 1]
    tab[i, 2] = round(tree$edge.length[r])
    if (i <= n)
      tab[i, 3] = date[i]
  }

  #Starting point
  tree = reorder.phylo(tree, 'postorder')
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    if (i <= n)
      next
    children = which(tab[, 4] == i)
    if (length(children) != 0)
      tab[i, 3] = min(tab[children, 3]) - 1
  }
  i = n + 1
  children = which(tab[, 4] == i)
  tab[i, 1] = i
  tab[i, 3] = min(tab[children, 3]) - 1

  #MCMC
  l = likelihood(tab, rate, n)
  p = prior(tab, 1, n)
  record = matrix(NA, nbIts / 10, 2 + nrow(tab))
  for (i in 1:nbIts) {
    #Record
    if (i %% 10 == 0) {
      record[i / 10, nrow(tab) + 1] = l
      record[i / 10, nrow(tab) + 2] = rate
      for (j in 1:nrow(tab))
        record[i / 10, j] = tab[j, 3]
    }

    if (updateRate == 1) {
      #MH move assuming flat prior
      rate2=abs(rate+runif(1)-0.5)
      l2=likelihood(tab,rate2,n)
      if (log(runif(1))<l2-l) {l=l2;rate=rate2}
    }

    if (updateRate == 2) {
      #Gibbs move assuming Exp(1) prior on rate
      lengths=tab[-(n+1),3]-tab[tab[-(n+1),4],3]
      muts=tab[-(n+1),2]
      rate=rgamma(1,1+sum(muts),sum(lengths))
      l=likelihood(tab,rate,n)
    }

    #MH to update dates
    for (j in (n + 1):nrow(tab)) {
      old = tab[j, 3]
      tab[j, 3] = old + runif(1) - 0.5
      l2 = likelihood(tab, rate, n)
      p2 = prior(tab, 1, n)
      if (log(runif(1)) < l2 - l + p2 - p)
        {l = l2; p = p2}
      else
        tab[j, 3] = old
    }
  }

  #Output
  meanRec = colMeans(record[(nrow(record) / 2):nrow(record), ])
  for (i in 1:nrow(tree$edge))
    tree$edge.length[i] = meanRec[tree$edge[i, 2]] - meanRec[tree$edge[i, 1]]
  return(list(
    tree = tree,
    record = record,
    rootdate = meanRec[n + 1]
  ))
}

#' Likelihood function
#' @return log-likelihood
likelihood = function(tab, rate, n) {
  if (rate < 0)
    return(-Inf)
  #n=ceiling(nrow(tab)/2)
  t2 = tab[-(n + 1), ]
  lengths = t2[, 3] - tab[t2[, 4], 3]
  if (min(lengths) < 0)
    return(-Inf)
  muts = t2[, 2]
  return(sum(-lengths * rate + muts * log(lengths * rate)))
}

#' Coalescent prior function
#'@return The log-prior in Eq (1) of Drummond et al (2002) Genetics
coalprior = function(tab, neg, n) {
  p = -log(neg) * (n - 1)
  MySort <- sort(tab[, 3], decreasing = T, index.return = TRUE)
  ind <- MySort$ix
  s <- MySort$x
  k = 1
  for (i in 2:length(s)) {
    p = p - (k * (k - 1) / (2 * neg) * (s[i-1] - s[i]))
    if (ind[i] <= n)
      k = k + 1
    else
      k = k - 1
  }
  if (k != 1)
    print('error')
  return(p)
}

#' Compute dates of leaves for a given tree and date of root
#' @param phy Tree
#' @param rootdate Date of root
#' @return Dates of leaves
#' @export
leafDates = function (phy,rootdate=0) {
nsam=length(phy$tip.label)
dates=rep(rootdate,nsam)
for (i in 1:nsam) {
  w=i
  while (1) {
    r=which(phy$edge[,2]==w)
    if (length(r)==0) break
    dates[i]=dates[i]+phy$edge.length[r]
    w=phy$edge[r,1]
  }
}
}
