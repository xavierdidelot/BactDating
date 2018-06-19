#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param initMu Initial rate of substitutions per genome (not per site), or zero to use root-to-tip estimate
#' @param initAlpha Initial coalescent time unit
#' @param initSigma Initial std on per-branch substitution rate (only used in relaxedgamma and mixedgamma models)
#' @param updateMu Whether or not to update the substitution rate
#' @param updateAlpha Whether or not to update the coalescent time unit
#' @param updateSigma Whether or not to per-branch substitution rate (only used in relaxedgamma and mixedgamma models)
#' @param updateRoot Whether or not to use the root finding algorithm
#' @param nbIts Number of MCMC iterations to perform
#' @param thin Thining interval between recorded MCMC samples
#' @param useCoalPrior Whether or not to use a coalescent prior on the tree
#' @param model Which model to use (poisson or negbin or strictgamma or relaxedgamma or mixedgamma)
#' @param useRec Whether or not to use results from previous recombination analysis
#' @param showProgress Whether or not to show a progress bar
#' @return Dating results
#' @export
bactdate = function(tree, date, initMu = NA, initAlpha = NA, initSigma = NA, updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T, nbIts = 10000, thin=ceiling(nbIts/1000), useCoalPrior = T,  model = 'mixedgamma', useRec = F, showProgress = F)
{
  #Rerranging of dates, if needed
  if (is.matrix(date)) {
    rangedate=date
    if (!is.null(rownames(rangedate))) rangedate=cbind(findDates(tree,rangedate[,1]),findDates(tree,rangedate[,2]))
    date=(rangedate[,1]+rangedate[,2])/2
  } else rangedate=NULL
  if (is.vector(date) && !is.null(names(date))) date=findDates(tree,date)
  if (is.null(rangedate)) rangedate=cbind(date,date)

  #Initial rooting of tree, if needed
  if (useRec==T && is.null(tree$unrec)) stop("To use recombination, the proportion of unrecombined needs to be input.")
  if (is.rooted(tree)==F) tree=initRoot(tree,date,useRec=useRec)
  testSignal=F

  #If the initial rate was not specified, start with the rate implied by root-to-tip analysis
  if (is.na(initMu)) {
    r=suppressWarnings(unname(roottotip(tree,date,showFig=F)$rate))
    if (is.na(r) || r<0) r=1
    initMu=r
  }
  mu=initMu
  if (is.na(initSigma)) initSigma=mu
  sigma=initSigma

  #Select prior function
  prior=function(...) return(0)
  if (useCoalPrior) prior=coalpriorC else updateAlpha=F

  #Selection likelihood function
  if (!is.element(model,c('poisson','poissonR','negbin','negbinR'))) tree$edge.length=pmax(tree$edge.length,1e-7)
  if (!is.element(model,c('relaxedgamma','relaxedgammaR','mixedgamma','negbin','negbinR'))) {updateSigma=F;sigma=0}
  if (model == 'mixedgamma') sigma=0
  if (model == 'poisson') likelihood=function(tab,mu,sigma) return(likelihoodPoissonC(tab,mu))
  if (model == 'poissonR') likelihood=function(tab,mu,sigma) return(likelihoodPoisson(tab,mu))
  if (model == 'negbinR') likelihood=function(tab,mu,sigma) return(likelihoodNegbin(tab,mu,sigma))
  if (model == 'strictgamma') likelihood=function(tab,mu,sigma) return(likelihoodGammaC(tab,mu))
  if (model == 'strictgammaR') likelihood=function(tab,mu,sigma) return(likelihoodGamma(tab,mu))
  if (model == 'relaxedgamma'||model == 'mixedgamma') likelihood=function(tab,mu,sigma) return(likelihoodRelaxedgammaC(tab,mu,sigma))
  if (model == 'relaxedgammaR') likelihood=function(tab,mu,sigma) return(likelihoodRelaxedgamma(tab,mu,sigma))
  if (model == 'null') {updateMu=0;likelihood=function(tab,mu,sigma) return(0)}
  if (!exists('likelihood')) stop('Unknown model.')

  #Deal with missing dates
  rangedate[which(is.na(rangedate[,1])),1]=min(rangedate[,1],na.rm = T)
  rangedate[which(is.na(rangedate[,2])),2]=max(rangedate[,2],na.rm = T)
  misDates=which(rangedate[,1]!=rangedate[,2])
  date[misDates]=(rangedate[misDates,1]+rangedate[misDates,2])/2
  orderedleafdates=sort(date,decreasing = T)

  #Create table of nodes (col1=name,col2=observed substitutions on branch above,col3=date,col4=father,col5=unrecombined proportion, only if useRec=T)
  n = Ntip(tree)
  tab = matrix(NA, n + tree$Nnode, ifelse(useRec,5,4))
  nrowtab = nrow(tab)
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
      tab[i, 3] = min(tab[children,3]-0.01,mean(tab[children, 3]-tab[children,2]/mu))
  }
  i = n + 1
  children = which(tab[, 4] == i)
  tab[i, 1] = i
  tab[i, 3] = min(tab[children, 3]-tab[children,2]/mu)

  #Sample Alpha if no initial point provided
  if (is.na(initAlpha)) {
    s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
    k=cumsum(2*(s$ix<=n)-1)
    difs=s$x[1:(nrowtab-1)]-s$x[2:nrowtab]
    su=sum(k[1:(nrowtab-1)]*(k[1:(nrowtab-1)]-1)*difs)
    initAlpha=1/rgamma(1,shape=n+0.001-1,scale=2000/(su*1000+2))
  }
  alpha=initAlpha
  initHeight=max(tab[,3])-min(tab[,3])

  #Start of MCMC algorithm
  l = likelihood(tab, mu, sigma)
  orderednodedates=sort(tab[(n+1):nrow(tab),3],method='quick',decreasing = T)
  p = prior(orderedleafdates,orderednodedates,alpha)
  record = matrix(NA, floor(nbIts / thin), nrow(tab)*3 + 6)
  colnames(record)<-c(rep(NA,nrow(tab)*3),'likelihood','mu','sigma','alpha','prior','root')
  if (showProgress) pb <- txtProgressBar(min=0,max=nbIts,style = 3)
  children=vector("list", max(tab[,4],na.rm = T))
  for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
  curroot=NA
  rootchildren=which(tab[,4]==(n+1))
  for (j in 1:nrow(tree$edge)) if (setequal(rootchildren,tree$edge[j,])) {curroot=j;break}
  if (is.na(curroot)) curroot=min(which(tree$edge[,1]==(n+1)))

  for (i in 1:nbIts) {
    #Record
    if (i %% thin == 0) {
      if (showProgress) setTxtProgressBar(pb, i)
      record[i / thin, 1:nrow(tab)] = tab[, 3]
      record[i / thin, (1:nrow(tab))+nrow(tab)]=tab[, 4]
      record[i / thin, (1:nrow(tab))+2*nrow(tab)]=tab[, 2]
      record[i / thin, 'likelihood'] = l
      record[i / thin, 'mu'] = mu
      record[i / thin, 'sigma'] = sigma
      record[i / thin, 'alpha'] = alpha
      record[i / thin, 'prior'] = p
      record[i / thin, 'root'] = curroot
      if (testSignal) record[i / thin, 'signal'] = all(tab[1:n,3]==date)
    }

    if (updateMu == T) {
      #MH move using Gamma(1e-3,1e3) prior
      mu2=abs(rnorm(1,mu,0.1*ifelse(initMu==0,1e-3,initMu)))
      l2=likelihood(tab,mu2,sigma)
      if (log(runif(1))<l2-l+dgamma(mu2,shape=1e-3,scale=1e3,log=T)-dgamma(mu,shape=1e-3,scale=1e3,log=T)) {l=l2;mu=mu2}
    }

    if (updateSigma == T && sigma>0) {
      #MH move using Gamma(1e-3,1e3) prior
      sigma2=abs(rnorm(1,sigma,0.1*ifelse(initSigma==0,1e-3,initSigma)))
      l2=likelihood(tab,mu,sigma2)
      if (log(runif(1))<l2-l+dgamma(sigma2,shape=1e-3,scale=1e3,log=T)-dgamma(sigma,shape=1e-3,scale=1e3,log=T)) {l=l2;sigma=sigma2}
    }

    if (model == 'mixedgamma') {
      #Reversible-jump move
      if (sigma==0) {
        sigma2=rexp(1)
        qratio=sigma2#-dexp(sigma2,1,log=T)
        pratio=dgamma(sigma2,shape=1e-3,scale=1e3,log=T)
      } else {
        sigma2=0
        qratio=-sigma#dexp(sigma,1,log=T)
        pratio=-dgamma(sigma,shape=1e-3,scale=1e3,log=T)
      }
      l2=likelihood(tab,mu,sigma2)
      if (log(runif(1))<l2-l+pratio+qratio) {l=l2;sigma=sigma2}
    }

    if (updateAlpha == T) {
      #Gibbs move using inverse-gamma prior
      s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
      k=cumsum(2*(s$ix<=n)-1)
      difs=s$x[1:(nrowtab-1)]-s$x[2:nrowtab]
      su=sum(k[1:(nrowtab-1)]*(k[1:(nrowtab-1)]-1)*difs)
      alpha=1/rgamma(1,shape=n+0.001-1,scale=2000/(su*1000+2))
      p=prior(orderedleafdates,orderednodedates,alpha)
    }

    #MH to update internal dates
    rn=rnorm(nrow(tab)-n,0,ifelse(initHeight==0,1,initHeight)*0.05)
    for (j in c((n + 1):nrow(tab))) {
      r=rn[j-n]
      old=tab[j,3]
      new=old+r
      if (r<0&&(!is.na(tab[j,4])&&new<tab[tab[j,4],3])) next#can't be older than father
      if (r>0&&new>min(tab[children[[j]],3])) next#can't be younger than sons
      if (j>(n+1)) {mintab=rbind(tab[children[[j]],],tab[tab[j,4],],tab[j,]);mintab[,4]=c(4,4,NA,3)}
      else {mintab=rbind(tab[children[[j]],],tab[j,]);mintab[,4]=c(3,3,NA)}
      l2=l-likelihood(mintab,mu,sigma)
      mintab[nrow(mintab),3]=new
      l2=l2+likelihood(mintab,mu,sigma)
      #tab2=tab;tab2[j,3]=new;l2full=likelihood(tab2,mu,sigma);if (abs(l2-l2full)>1e-10) print(sprintf('error %f %f',l2,l2full))
      changeinorderedvec(orderednodedates,old,new)
      p2 = prior(orderedleafdates, orderednodedates, alpha)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2;tab[j,3]=new}
      else
        changeinorderedvec(orderednodedates,new,old)
    }

    #MH to update missing leaf dates
    for (j in misDates) {
      old = tab[j, 3]
      new = rnorm(1,old,(rangedate[j,2]-rangedate[j,1])*0.05)
      if (new-old<0&&(!is.na(tab[j,4])&&new<tab[tab[j,4],3])) next#can't be older than father
      if (new>rangedate[j,2]||new<rangedate[j,1]) next#stay within prior range
      mintab=rbind(tab[j,],tab[tab[j,4],])
      mintab[,4]=c(2,NA)
      l2=l-likelihood(mintab,mu,sigma)
      mintab[1,3]=new
      l2=l2+likelihood(mintab,mu,sigma)
      #tab2=tab;tab2[j,3]=new;l2full=likelihood(tab2,mu,sigma);if (abs(l2-l2full)>1e-10) print(sprintf('error %f %f',l2,l2full))
      changeinorderedvec(orderedleafdates,old,new)
      p2 = prior(orderedleafdates, orderednodedates, alpha)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2;tab[j,3]=new}
      else
        changeinorderedvec(orderedleafdates,new,old)
    }

    if (updateRoot == T || updateRoot == 'branch') {
      #Move root on current branch
      root=which(is.na(tab[,4]))
      sides=which(tab[,4]==root)
      old=tab[sides,2]
      r=runif(1)
      tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
      l2=likelihood(tab,mu, sigma)
      if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
  }

    if (updateRoot == T) {
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
        l2=likelihood(tab,mu,sigma)
        if (log(runif(1))<l2-l+log(oldtab[a,2]/tab[right,2]))
          {l=l2
          children=vector("list", max(tab[,4],na.rm = T))
          for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
          curroot=NA
          rootchildren=which(tab[,4]==(n+1))
          for (j in 1:nrow(tree$edge)) if (setequal(rootchildren,tree$edge[j,])) {curroot=j;break}
          if (is.na(curroot)) curroot=min(which(tree$edge[,1]==(n+1)))
        } else tab=oldtab
      }
    }

    if (testSignal) {
      #Model jump to test strength of temporal signal
      maxdates=max(date)
      if (all(tab[1:n,3]==maxdates)) {
        tab[1:n,3]=date
        l2=likelihood(tab,mu,sigma)
        if (log(runif(1))<l2-l) l=l2 else tab[1:n,3]=maxdates
      } else {
        tab[1:n,3]=maxdates
        l2=likelihood(tab,mu,sigma)
        if (log(runif(1))<l2-l) l=l2 else tab[1:n,3]=date
      }
    }

  }#End of MCMC loop

  #Output
  if (testSignal) tab[1:n,3]=date
  inputtree = tree
  bestroot = as.numeric(names(sort(table(record[floor(nrow(record) / 2):nrow(record),'root']),decreasing=T)[1]))
  bestrows = intersect(floor(nrow(record) / 2):nrow(record),which(record[,'root']==bestroot))
  meanRec = colMeans(record[bestrows, ,drop=F])
  for (i in 1:nrow(tree$edge)) {
    tree$edge[i,1]=record[bestrows[1],nrow(tab)+tree$edge[i,2]]
    tree$edge.length[i] = meanRec[tree$edge[i, 2]] - meanRec[tree$edge[i, 1]]
    tree$subs[i] = meanRec[nrow(tab)*2+tree$edge[i,2]]
  }

  #Calculate DIC
  meantab=tab
  meantab[,3]=meanRec[1:nrow(tab)]
  meantab[,4]=meanRec[(1:nrow(tab))+nrow(tab)]
  dic=-2*likelihood(meantab,meanRec['mu'],meanRec['sigma'])+var(-2*record[bestrows,'likelihood'])

  tree$root.time = max(date)-max(leafDates(tree))
  CI = matrix(NA, nrow(tab), 2)
  for (i in 1:nrow(tab)) {
    s=sort(record[bestrows,i])
    CI[i,1]=s[ceiling(length(s)*0.025)]
    CI[i,2]=s[max(1,floor(length(s)*0.975))]
  }
  out = list(
    inputtree = inputtree,
    tree = tree,
    model = model,
    record = record,
    rootdate = unname(meanRec[n + 1]),
    rootprob = length(bestrows)*2/nrow(record),
    CI = CI,
    dic = dic
  )
  class(out) <- 'resBactDating'
  if (model=='mixedgamma') out$pstrict=length(which(record[,'sigma']==0))/nrow(record)
  if (testSignal) out$psignal=length(which(record[,'signal']==1))/nrow(record)
  return(out)
}

#' Perform Bayesian comparison between two models based on DIC values
#' @param res1 An output from the bactdate function
#' @param res2 Another output from the bactdate function
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
