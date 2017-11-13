#' Root to tip correlation
#' @param tree Phylogenetic tree
#' @param date Dates of sampling
#' @param permTest Number of permutations to perform to compute the p-value using a permutation test
#' @param showFig Whether or not to show the root-to-tip regression figure
#' @param showPredInt To show 95percent confidence intervals, can be 'poisson' or 'gamma'
#' @param showTree Whether to show the tree or not
#' @return List containing estimated clock rate, date of origin and p-value
#' @importFrom graphics abline
#' @export
roottotip = function(tree,date,permTest=10000,showFig=T,showPredInt='gamma',showTree=T)
{
  if (var(date,na.rm=T)==0) {warning('All dates are identical.');return(list(rate=NA,ori=NA,pvalue=NA))}
  n=length(date)
  ys=leafDates(tree)
  res=lm(ys~date)
  ori=-coef(res)[1]/coef(res)[2]
  rate=coef(res)[2]
  r2=summary(res)$r.squared
  correl=cor(date,ys,use='complete.obs')
  #pvalue=summary(res)$coefficients[,4][2]
  #print(c(r2,correl^2))#Equal

  pvalue=0
  for (i in 1:permTest) {
    date2=sample(date,n,replace=F)
    correl2=cor(date2,ys,use='complete.obs')
    if (correl2>correl) pvalue=pvalue+1/permTest
  }

  if (rate<0) {warning('The linear regression suggests a negative rate.');return(list(rate=rate,ori=ori,pvalue=pvalue))}
  if (showFig==F) return(list(rate=rate,ori=ori,pvalue=pvalue))
  par(xpd=NA,oma = c(0, 0, 2, 0))
  if (showTree) {
    par(mfrow=c(1,2))
    plot(tree)
    axisPhylo(1,backward = F)
  }
  plot(date,ys,xlab='Sampling date',ylab='Root-to-tip distance',xaxs='i',yaxs='i',pch=19,ylim=c(0,max(ys)),xlim=c(ori,max(date,na.rm = T)))
  par(xpd=F)
  abline(res,lwd=2)
  xs=seq(ori,max(date,na.rm = T),0.1)
  plim=0.05
  if (showPredInt=='poisson') {
  lines(xs,qpois(  plim/2,(xs-ori)*rate),lty='dashed')
  lines(xs,qpois(1-plim/2,(xs-ori)*rate),lty='dashed')
  }
  if (showPredInt=='gamma') {
    lines(xs,qgamma(  plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
    lines(xs,qgamma(1-plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
  }
  if (pvalue==0) mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p<%.2e',rate,ori,r2,1/permTest), outer = TRUE, cex = 1.5)
  else           mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p=%.2e',rate,ori,r2,pvalue), outer = TRUE, cex = 1.5)
  return(list(rate=rate,ori=ori,pvalue=pvalue))
}

#' Compute dates of leaves for a given tree and date of root
#' @param phy Tree
#' @return Dates of leaves
#' @export
leafDates = function (phy) {
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
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
  return(dates)
}
