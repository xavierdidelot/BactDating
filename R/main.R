#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param initRate Initial rate of substitutions per genome (not per site)
#' @param nbIts Number of MCMC iterations to perform
#' @param useCoalPrior Whether or not to use a coalescent prior on the tree
#' @param updateRate Whether or not to update the substitution rate
#' @param initNeg Initial rate of coalescence, equal to Ne*g
#' @param updateNeg Whether or not to update the neg parameter
#' @param model Which model to use (1=Poisson, 2=NegBin, 3=Gamma)
#' @param findRoot Root finding algorithm (0=none, 1=on preset branch, 2=anywhere)
#' @return Dating results
#' @export
credate = function(tree, date, initRate = 1, nbIts = 1000, useCoalPrior = T, updateRate = 2, initNeg = 1, updateNeg = 2, model = 1, findRoot = 0)
{
  n = Ntip(tree)
  rate = initRate
  neg = initNeg

  prior=function(tab,neg) return(0)
  if (useCoalPrior) prior=coalprior

  if (model == 1) likelihood=likelihoodPoisson
  if (model>1 && updateRate==2) updateRate=1
  if (model == 2) likelihood=function(tab,rate) return(likelihoodNegbin(tab,r=rate,phi=1))
  if (model == 3) likelihood=function(tab,rate) return(likelihoodGamma(tab,rate))

  #Deal with missing dates
  misDates=which(is.na(date))
  date[misDates]=mean(date,na.rm = T)

  #Create table of nodes (col1=name,col2=expected mutation on branch above,col3=date,col4=father)
  tab = matrix(NA, n + tree$Nnode, 4)
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    tab[i, 1] = i
    tab[i, 4] = tree$edge[r, 1]
    if (model == 1 || model == 2) tab[i, 2] = round(tree$edge.length[r])
    else tab[i, 2] = pmax(1e-7,tree$edge.length[r])#NOTE THIS NEEDS TO BE CONFIRMED
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
  l = likelihood(tab, rate)
  p = prior(tab, neg)
  thin = nbIts/1000
  record = matrix(NA, nbIts / thin, 4 + nrow(tab))
  for (i in 1:nbIts) {
    #Record
    if (i %% thin == 0) {
      record[i / thin, nrow(tab) + 1] = l
      record[i / thin, nrow(tab) + 2] = rate
      record[i / thin, nrow(tab) + 3] = neg
      record[i / thin, nrow(tab) + 4] = p
      record[i / thin, 1:nrow(tab)] = tab[1:nrow(tab), 3]
    }

    if (updateRate == 1) {
      #MH move assuming flat prior
      rate2=abs(rate+runif(1)-0.5)
      l2=likelihood(tab,rate2)
      if (log(runif(1))<l2-l) {l=l2;rate=rate2}
    }

    if (updateRate == 2) {
      #Gibbs move assuming Poisson model and Exp(1) prior on rate
      lengths=tab[-(n+1),3]-tab[tab[-(n+1),4],3]
      muts=tab[-(n+1),2]
      rate=rgamma(1,1+sum(muts),sum(lengths))
      l=likelihood(tab,rate)
    }

    if (updateNeg == 1) {
      #MH move assuming flat prior
      neg2=abs(neg+runif(1)-0.5)
      p2=prior(tab,neg2)
      if (log(runif(1))<p2-p) {p=p2;neg=neg2}
    }

    if (updateNeg == 2) {
      #Gibbs move assuming inverse-gamma prior
      s <- sort(tab[, 3], decreasing = T, index.return = TRUE)
      k=cumsum(2*(s$ix<=n)-1)
      difs=-diff(s$x)#s$x[1:(nrow(tab)-1)]-s$x[2:nrow(tab)]
      b=sum(k[2:nrow(tab)]*(k[2:nrow(tab)]-1)*difs)/2
      neg=1/rgamma(1,shape=n-1,scale=1/b)
      p=prior(tab,neg)
    }

    #MH to update dates
    for (j in c(misDates,(n + 1):nrow(tab))) {
      old = tab[j, 3]
      ru=runif(1)
      tab[j, 3] = old + ru - 0.5
      if (ru<0.5&&(!is.na(tab[j,4])&&tab[j,3]<tab[tab[j,4],3])) {tab[j,3]=old;next}
      if (ru>0.5&&j>n&&tab[j,3]>min(tab[which(tab[,4]==j),3])) {tab[j,3]=old;next}
      if (j<=n&&(tab[j,3]==max(tab[,3])||tab[j,3]==min(tab[,3]))) {tab[j,3]=old;next}
      l2 = likelihood(tab, rate)
      p2 = prior(tab, neg)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2}
      else
        tab[j, 3] = old
    }

    if (findRoot>0) {
      #Move root on its branch
      root=which(tab[,3]==max(tab[,3]))
      sides=which(tab[,4]==root)
      old=tab[sides,2]
      r=runif(1)
      tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
      if (model==1||model==2) tab[sides,2]=round(tab[sides,2])
      l2=likelihood(tab,rate)
      if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
    }

    if (findRoot==2 && i<(nbIts/4)) {
      #Move root branch
      root=which(tab[,3]==min(tab[,3]))
      sides=which(tab[,4]==root)
      if (tab[sides[1],3]<tab[sides[2],3]) {left=sides[1];right=sides[2]} else {left=sides[2];right=sides[1]}
      if (left>n) {
        oldtab=tab
        ab=which(tab[,4]==left)
        if (runif(1)<0.5) {a=ab[1];b=ab[2]} else {a=ab[2];b=ab[1]}
        tab[a,4]=root
        tab[right,4]=left
        r=runif(1)
        tab[a,2]=oldtab[a,2]*r
        tab[left,2]=oldtab[a,2]*(1-r)
        if (model==1||model==2) {tab[a,2]=round(tab[a,2]);tab[left,2]=round(tab[left,2])}
        tab[right,2]=oldtab[right,2]+oldtab[left,2]
        l2=likelihood(tab,rate)
        if (log(runif(1))<l2-l) l=l2 else tab=oldtab
      }
    }
  }

  #Output
  meanRec = colMeans(record[(nrow(record) / 2):nrow(record), ])
  for (i in 1:nrow(tree$edge)) {
    tree$edge[i,1]=tab[tree$edge[i,2],4]
    tree$edge.length[i] = meanRec[tree$edge[i, 2]] - meanRec[tree$edge[i, 1]]
  }
  tree$root.time = max(date)-max(leafDates(tree))
  CI = matrix(NA, nrow(tab), 2)
  for (i in 1:nrow(tab)) {
    s=sort(record[(nrow(record)/2):nrow(record),i])
    CI[i,1]=s[floor(length(s)*0.025)]
    CI[i,2]=s[ceiling(length(s)*0.975)]
  }
  out = list(
    tree = tree,
    record = record,
    rootdate = meanRec[n + 1],
    CI = CI
  )
  class(out) <- 'resCreDating'
  return(out)
}

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
  return(sum(dgamma(muts,shape=rate*lengths,scale=1,log=T)))
}

#' Coalescent prior function
#' @param tab Table of nodes
#' @param neg Coalescent rate
#'@return The log-prior in Eq (1) of Drummond et al (2002) Genetics
coalprior = function(tab, neg) {
  n = ceiling(nrow(tab)/2)
  p = -log(neg) * (n - 1)
  l=nrow(tab)
  s <- sort(tab[, 3], decreasing = T, index.return = TRUE)
  k=cumsum(2*(s$ix<=n)-1)
  difs=-diff(s$x)#s$x[1:(l-1)]-s$x[2:l]
  p=p-sum(k[2:l]*(k[2:l]-1)*difs)/(2*neg)
  if (k[length(k)] != 1)
    stop('error')
  return(p)
}

#' Plotting methods
#' @param x Output from running function credating
#' @param type Type of plot to do. Currently either 'tree' or 'treeCI' or 'trace'
#' @param ... Additional parameters are passed on
#' @return Nothing
#' @importFrom grDevices rgb
#' @export
plot.resCreDating = function(x, type='tree', ...) {

  if (type=='tree') {
    plot.phylo(x$tree, ...)
    axisPhylo(backward = F)
  }

  if (type=='treeCI') {
    plot.phylo(x$tree, x.lim=c(min(x$CI),max(x$CI))-x$tree$root.time,...)
    axisPhylo(backward = F)
    obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    transblue=rgb(0,0,1,0.4)
    for(i in 1:(Nnode(x$tree)+Ntip(x$tree)))
      if (x$CI[i,1]!=x$CI[i,2])
      lines(x=c(x$CI[i,1],x$CI[i,2])-x$tree$root.time,
            y=rep(obj$yy[i],2),lwd=11,lend=0,
            col=transblue)
    points(obj$xx[1:x$tree$Nnode+Ntip(x$tree)],
           obj$yy[1:x$tree$Nnode+Ntip(x$tree)],pch=19,col="blue",
           cex=1.8)
  }

  if (type=='trace') {
    nc=ncol(x$record)
    par(mfrow=c(2,3))
    plot(x$record[,nc-3],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,nc  ],main='Prior',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,Ntip(x$tree)+1],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,nc-2],main='Substitution rate',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,nc-1],main='Coalescent rate',type='l',xlab='Sampled iterations',ylab='')
  }
}
