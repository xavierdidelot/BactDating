#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param initRate Initial rate of substitutions per genome (not per site)
#' @param nbIts Number of MCMC iterations to perform
#' @param thin Thining interval between recorded MCMC samples
#' @param useCoalPrior Whether or not to use a coalescent prior on the tree
#' @param updateRate Whether or not to update the substitution rate
#' @param initNeg Initial rate of coalescence, equal to Ne*g
#' @param updateNeg Whether or not to update the neg parameter
#' @param initRatevar Initial variance on per-branch substituion rate (only used in relaxedgamma model)
#' @param updateRatevar Whether or not to per-branch substituion rate (only used in relaxedgamma model)
#' @param model Which model to use (poisson or negbin or gamma or relaxedgamma)
#' @param updateRoot Root finding algorithm (0=none, 1=on preset branch, 2=anywhere)
#' @param showProgress Whether or not to show a progress bar
#' @return Dating results
#' @export
credate = function(tree, date, initRate = 1, nbIts = 10000, thin=ceiling(nbIts/1000), useCoalPrior = T, updateRate = 1, initNeg = 1, updateNeg = 2, initRatevar = 1, updateRatevar = F,  model = 'gamma', updateRoot = 0, showProgress = T)
{
  n = Ntip(tree)
  rate = initRate
  ratevar = initRatevar
  neg = initNeg
  if (is.rooted(tree)==F) {
    first=which(date==min(date,na.rm = T))[1]
    tree=root(tree,outgroup=first,resolve.root=T)
    w=which(tree$edge[,1]==Ntip(tree)+1)
    tree$edge.length[w]=rep(sum(tree$edge.length[w])/2,2)
  }

  prior=function(...) return(0)
  if (useCoalPrior) prior=coalpriorC else updateNeg=0

  if ((model == 'gamma'||model == 'negbin'||model=='relaxedgamma') && updateRate==2) updateRate=1#Gibbs move on rate is only available for poisson model
  if (model == 'gamma' ||model == 'relaxedgamma') tree$edge.length=pmax(tree$edge.length,1e-7)
  if (model == 'poisson') likelihood=function(tab,rate,ratevar) return(likelihoodPoissonC(tab,rate))
  if (model == 'negbin') likelihood=function(tab,rate,ratevar) return(likelihoodNegbin(tab,r=rate,phi=1))
  if (model == 'gamma') likelihood=function(tab,rate,ratevar) return(likelihoodGammaC(tab,rate))
  if (model == 'relaxedgamma') likelihood=function(tab,rate,ratevar) return(likelihoodRelaxedgammaC(tab,rate,ratevar))

  #Deal with missing dates
  misDates=which(is.na(date))
  date[misDates]=mean(date,na.rm = T)
  ordereddate=sort(date,decreasing = T)

  #Create table of nodes (col1=name,col2=expected mutation on branch above,col3=date,col4=father)
  tab = matrix(NA, n + tree$Nnode, 4)
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    tab[i, 1] = i
    tab[i, 4] = tree$edge[r, 1]
    if (model == 'poisson' || model == 'negbin') tab[i, 2] = round(tree$edge.length[r])
    else tab[i, 2] = tree$edge.length[r]
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
  l = likelihood(tab, rate, ratevar)
  p = prior(ordereddate, tab[(n+1):nrow(tab),3], neg)
  record = matrix(NA, floor(nbIts / thin), nrow(tab)*2 + 6)
  colnames(record)<-c(rep(NA,nrow(tab)*2),'likelihood','rate','ratevar','neg','prior','root')
  if (showProgress) pb <- txtProgressBar(min=0,max=nbIts,style = 3)
  for (i in 1:nbIts) {
    #Record
    if (i %% thin == 0) {
      if (showProgress) setTxtProgressBar(pb, i)
      rootchildren=which(tab[,4]==(n+1))
      curroot=NA
      for (j in 1:nrow(tree$edge)) if (setequal(rootchildren,tree$edge[j,])) {curroot=j;break}
      if (is.na(curroot)) curroot=which(tree$edge[,1]==(n+1))[1]
      record[i / thin, 1:nrow(tab)] = tab[, 3]
      record[i / thin, (1:nrow(tab))+nrow(tab)]=tab[, 4]
      record[i / thin, 'likelihood'] = l
      record[i / thin, 'rate'] = rate
      record[i / thin, 'ratevar'] = ratevar
      record[i / thin, 'neg'] = neg
      record[i / thin, 'prior'] = p
      record[i / thin, 'root'] = curroot
    }

    if (updateRate == 1) {
      #MH move assuming flat prior
      rate2=abs(rnorm(1,rate,0.1*initRate))
      l2=likelihood(tab,rate2,ratevar)
      if (log(runif(1))<l2-l) {l=l2;rate=rate2}
    }

    if (updateRate == 2) {
      #Gibbs move assuming Poisson model and Exp(1) prior on rate
      lengths=tab[-(n+1),3]-tab[tab[-(n+1),4],3]
      muts=tab[-(n+1),2]
      rate=rgamma(1,1+sum(muts),sum(lengths))
      l=likelihood(tab,rate,ratevar)
    }

    if (updateRatevar == 1) {
      #MH move assuming flat prior
      ratevar2=abs(rnorm(1,ratevar,0.1*initRatevar))
      l2=likelihood(tab,rate,ratevar2)
      if (log(runif(1))<l2-l) {l=l2;ratevar=ratevar2}
    }

    if (updateNeg == 1) {
      #MH move assuming flat prior
      neg2=abs(rnorm(1,neg,0.1*initNeg))
      p2=prior(ordereddate, tab[(n+1):nrow(tab),3],neg2)
      if (log(runif(1))<p2-p) {p=p2;neg=neg2}
    }

    if (updateNeg == 2) {
      #Gibbs move assuming inverse-gamma prior
      s <- sort(tab[, 3], decreasing = T, index.return = TRUE)
      k=cumsum(2*(s$ix<=n)-1)
      difs=-diff(s$x)#s$x[1:(nrow(tab)-1)]-s$x[2:nrow(tab)]
      b=sum(k[2:nrow(tab)]*(k[2:nrow(tab)]-1)*difs)/2
      neg=1/rgamma(1,shape=n-1,scale=1/b)
      p=prior(ordereddate, tab[(n+1):nrow(tab),3],neg)
    }

    #MH to update internal dates
    for (j in c((n + 1):nrow(tab))) {
      old = tab[j, 3]
      tab[j, 3] = old + rnorm(1)
      if (tab[j,3]-old<0&&(!is.na(tab[j,4])&&tab[j,3]<tab[tab[j,4],3])) {tab[j,3]=old;next}#can't be older than father
      if (tab[j,3]-old>0&&j>n&&tab[j,3]>min(tab[which(tab[,4]==j),3])) {tab[j,3]=old;next}#can't be younger than sons
      l2 = likelihood(tab, rate, ratevar)
      p2 = prior(ordereddate, tab[(n+1):nrow(tab),3], neg)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2}
      else
        tab[j, 3] = old
    }

    #MH to update missing leaf dates
    for (j in misDates) {
      old = tab[j, 3]
      tab[j, 3] = old + rnorm(1)
      if (tab[j,3]-old<0&&(!is.na(tab[j,4])&&tab[j,3]<tab[tab[j,4],3])) {tab[j,3]=old;next}#can't be older than father
      if (tab[j,3]==max(tab[,3])||tab[j,3]==min(tab[,3])) {tab[j,3]=old;next}#stay within prior range
      l2 = likelihood(tab, rate, ratevar)
      ordereddate2=sort(tab[1:n,3],decreasing = T)
      p2 = prior(ordereddate2, tab[(n+1):nrow(tab),3], neg)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2; ordereddate = ordereddate2}
      else
        tab[j, 3] = old
    }

    if (updateRoot>0) {
      #Move root on current branch
      root=which(tab[,3]==max(tab[,3]))
      sides=which(tab[,4]==root)
      old=tab[sides,2]
      r=runif(1)
      tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
      if (model=='poisson'||model=='negbin') tab[sides,2]=round(tab[sides,2])
      l2=likelihood(tab,rate, ratevar)
      if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
    }

    if (updateRoot==2) {
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
        if (model=='poisson'||model=='negbin') {tab[a,2]=round(tab[a,2]);tab[left,2]=round(tab[left,2])}
        tab[right,2]=oldtab[right,2]+oldtab[left,2]
        l2=likelihood(tab,rate, ratevar)
        if (log(runif(1))<l2-l+log(oldtab[a,2]/tab[right,2])) l=l2 else tab=oldtab
      }
    }
  }

  #Output
  inputtree = tree
  bestroot = as.numeric(names(sort(table(record[floor(nrow(record) / 2):nrow(record),'root']),decreasing=T)[1]))
  bestrows = intersect(floor(nrow(record) / 2):nrow(record),which(record[,'root']==bestroot))
  meanRec = colMeans(record[bestrows, ])
  for (i in 1:nrow(tree$edge)) {
    tree$edge[i,1]=record[bestrows[1],nrow(tab)+tree$edge[i,2]]#tab[tree$edge[i,2],4]
    tree$edge.length[i] = meanRec[tree$edge[i, 2]] - meanRec[tree$edge[i, 1]]
  }
  tree$root.time = max(date)-max(leafDates(tree))
  CI = matrix(NA, nrow(tab), 2)
  for (i in 1:nrow(tab)) {
    s=sort(record[bestrows,i])
    CI[i,1]=s[ceiling(length(s)*0.025)]
    CI[i,2]=s[floor(length(s)*0.975)]
  }
  out = list(
    inputtree = inputtree,
    tree = tree,
    record = record,
    rootdate = unname(meanRec[n + 1]),
    rootprob = length(bestrows)*2/nrow(record),
    CI = CI
  )
  class(out) <- 'resCreDating'
  return(out)
}

