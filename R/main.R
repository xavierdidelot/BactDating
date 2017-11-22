#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param initRate Initial rate of substitutions per genome (not per site), or zero to use root-to-tip estimate
#' @param initNeg Initial coalescent time unit
#' @param initRatevar Initial variance on per-branch substituion rate (only used in relaxedgamma model)
#' @param updateRate Whether or not to update the substitution rate
#' @param updateNeg Whether or not to update the coalescent time unit
#' @param updateRatevar Whether or not to per-branch substituion rate (only used in relaxedgamma model)
#' @param updateRoot Root finding algorithm (0=none, 1=on preset branch, 2=anywhere)
#' @param nbIts Number of MCMC iterations to perform
#' @param thin Thining interval between recorded MCMC samples
#' @param useCoalPrior Whether or not to use a coalescent prior on the tree
#' @param model Which model to use (poisson or negbin or gamma or relaxedgamma)
#' @param useRec Whether or not to use results from previous recombination analysis
#' @param showProgress Whether or not to show a progress bar
#' @return Dating results
#' @export
credate = function(tree, date, initRate = NA, initNeg = NA, initRatevar = NA, updateRate = T, updateNeg = T, updateRatevar = T, updateRoot = T, nbIts = 10000, thin=ceiling(nbIts/1000), useCoalPrior = T,  model = 'gamma', useRec = F, showProgress = F)
{
  #Rooting of tree without recombination
  if (is.rooted(tree)==F && useRec==F) {
    tree=initRoot(tree,date)#TODO USE THIS FUNCTION ALSO WITH RECOMBINATION
    #first=which(date==min(date,na.rm = T))[1]
    #tree=root(tree,outgroup=first,resolve.root=T)
    #w=which(tree$edge[,1]==Ntip(tree)+1)
    #tree$edge.length[w]=rep(sum(tree$edge.length[w])/2,2)
  }

  #Rooting of tree with recombination
  if (useRec==T && is.null(tree$unrec)) stop("To use recombination, the proportion of unrecombined needs to be input.\n")
  if (is.rooted(tree)==F && useRec==T) {
    first=which(date==min(date,na.rm = T))[1]
    tree$node.label=sprintf('n%d',1:Nnode(tree))
    edgenames=cbind(c(tree$tip.label,tree$node.label)[tree$edge[,1]],c(tree$tip.label,tree$node.label)[tree$edge[,2]])
    unrec=tree$unrec
    unrecfirst=unrec[which(tree$edge[,2]==first)]
    tree=root(tree,outgroup=first,resolve.root=T,edgelabel=F)
    w=which(tree$edge[,1]==Ntip(tree)+1)
    tree$edge.length[w]=rep(sum(tree$edge.length[w])/2,2)
    tree$unrec=rep(NA,nrow(tree$edge))
    tree$unrec[w]=unrecfirst
    for (i in 1:nrow(tree$edge)) {
      nams=c(tree$tip.label,tree$node.label)[tree$edge[i,]]
      for (j in 1:nrow(edgenames)) if (setequal(nams,edgenames[j,])) {tree$unrec[i]=unrec[j];break}
    }
  }

  #If the initial rate was not specified, start with the rate implied by root-to-tip analysis
  if (is.na(initRate)) {
    r=unname(roottotip(tree,date,showFig=F)$rate)
    if (is.na(r) || r<0) r=1
    initRate=r
  }
  rate=initRate
  if (is.na(initRatevar)) initRatevar=rate*rate
  ratevar=initRatevar

  #Select prior function
  prior=function(...) return(0)
  if (useCoalPrior) prior=coalpriorC else updateNeg=F

  #Selection likelihood function
  if (is.element(model,c('gamma','relaxedgamma','gammaR','relaxedgammaR'))) tree$edge.length=pmax(tree$edge.length,1e-7)
  if (model == 'poisson') likelihood=function(tab,rate,ratevar) return(likelihoodPoissonC(tab,rate))
  if (model == 'poissonR') likelihood=function(tab,rate,ratevar) return(likelihoodPoisson(tab,rate))
  if (model == 'negbin') likelihood=function(tab,rate,ratevar) return(likelihoodNegbin(tab,r=rate,phi=1))
  if (model == 'gamma') likelihood=function(tab,rate,ratevar) return(likelihoodGammaC(tab,rate))
  if (model == 'gammaR') likelihood=function(tab,rate,ratevar) return(likelihoodGamma(tab,rate))
  if (model == 'relaxedgamma') likelihood=function(tab,rate,ratevar) return(likelihoodRelaxedgammaC(tab,rate,ratevar))
  if (model == 'relaxedgammaR') likelihood=function(tab,rate,ratevar) return(likelihoodRelaxedgamma(tab,rate,ratevar))
  if (model != 'relaxedgamma'&&model!='relaxedgammaR') updateRatevar=F
  if (model == 'null') {updateRate=0;likelihood=function(tab,rate,ratevar) return(0)}
  if (!exists('likelihood')) stop('Unknown model.')

  #Deal with missing dates
  misDates=which(is.na(date))
  date[misDates]=mean(date,na.rm = T)
  ordereddate=sort(date,decreasing = T)

  #Create table of nodes (col1=name,col2=observed substitutions on branch above,col3=date,col4=father,col5=unrecombined proportion, only if useRec=T)
  n = Ntip(tree)
  tab = matrix(NA, n + tree$Nnode, ifelse(useRec,5,4))
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    tab[i, 1] = i
    tab[i, 4] = tree$edge[r, 1]
    tab[i, 2] = tree$edge.length[r]
    if (useRec) tab[i,5]=tree$unrec[r]
    if (i <= n)
      tab[i, 3] = date[i]
  }

  #Create initial tree
  tree = reorder.phylo(tree, 'postorder')
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    if (i <= n)
      next
    children = which(tab[, 4] == i)
    if (length(children) != 0)
      tab[i, 3] = min(tab[children, 3]-tab[children,2]/rate)
  }
  i = n + 1
  children = which(tab[, 4] == i)
  tab[i, 1] = i
  tab[i, 3] = min(tab[children, 3]-tab[children,2]/rate)

  #Sample Neg if no initial point provided
  if (is.na(initNeg)) {
    s <- sort(tab[, 3], decreasing = T, index.return = TRUE)
    k=cumsum(2*(s$ix<=n)-1)
    difs=-diff(s$x)
    b=sum(k[2:nrow(tab)]*(k[2:nrow(tab)]-1)*difs)/2
    initNeg=1/rgamma(1,shape=n-1,scale=1/b)
  }
  neg=initNeg

  #Start of MCMC algorithm
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
      if (is.na(curroot)) curroot=min(which(tree$edge[,1]==(n+1)))
      record[i / thin, 1:nrow(tab)] = tab[, 3]
      record[i / thin, (1:nrow(tab))+nrow(tab)]=tab[, 4]
      record[i / thin, 'likelihood'] = l
      record[i / thin, 'rate'] = rate
      record[i / thin, 'ratevar'] = ratevar
      record[i / thin, 'neg'] = neg
      record[i / thin, 'prior'] = p
      record[i / thin, 'root'] = curroot
    }

    if (updateRate == T) {
      #MH move using flat prior
      rate2=abs(rnorm(1,rate,0.1*initRate))
      l2=likelihood(tab,rate2,ratevar)
      if (log(runif(1))<l2-l) {l=l2;rate=rate2}
    }

    if (updateRatevar == T) {
      #MH move using flat prior
      ratevar2=abs(rnorm(1,ratevar,0.1*initRatevar))
      l2=likelihood(tab,rate,ratevar2)
      if (log(runif(1))<l2-l) {l=l2;ratevar=ratevar2}
    }

    if (updateNeg == T) {
      #Gibbs move using inverse-gamma prior
      s <- sort(tab[, 3], decreasing = T, index.return = TRUE)
      k=cumsum(2*(s$ix<=n)-1)
      difs=-diff(s$x)
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

    if (updateRoot) {
      #Move root on current branch
      root=which(is.na(tab[,4]))
      sides=which(tab[,4]==root)
      old=tab[sides,2]
      r=runif(1)
      tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
      l2=likelihood(tab,rate, ratevar)
      if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
    }

    if (updateRoot) {
      #Move root branch
      root=which(is.na(tab[,4]))
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
        tab[right,2]=oldtab[right,2]+oldtab[left,2]
        if (useRec) tab[left,5]=tab[a,5]
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

  #Calculate DIC
  meantab=tab
  meantab[,3]=meanRec[1:nrow(tab)]
  meantab[,4]=meanRec[(1:nrow(tab))+nrow(tab)]
  dic=-2*likelihood(meantab,meanRec['rate'],meanRec['ratevar'])+var(-2*record[bestrows,'likelihood'])

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
    CI = CI,
    dic = dic
  )
  class(out) <- 'resCreDating'
  return(out)
}

#' Perform Bayesian comparison between two models based on DIC values
#' @param res1 An output from the credate function
#' @param res2 Another output from the credate function
#' @return Prints out results of model comparison
#' @export
modelcompare = function(res1,res2) {
  dic1=res1$dic
  dic2=res2$dic
  dif=dic2-dic1
  cat(sprintf('The first model has DIC=%.2f and the second model has DIC=%.2f.\n',dic1,dic2))
  if (dif>10) cat('Model 1 is definitely better.\n')
  if (dif>5 && dif<10) cat('Model 1 is slightly better.\n')
  if (abs(dif)<5) cat('The difference is not significant.\n')
  if (dif< -5 && dif>-10) cat('Model 2 is slightly better.\n')
  if (dif< -10) cat('Model 2 is definitely better.\n')
}
