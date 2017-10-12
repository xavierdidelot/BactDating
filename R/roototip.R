#' Root to tip correlation
#' @param tree Phylogenetic tree
#' @param date Dates of sampling
#' @importFrom graphics abline
#' @export
roottotip = function(tree,date,showTree=F)
{
  n=length(date)
  ys=leafDates(tree)
  res=lm(ys~date)
  ori=-coef(res)[1]/coef(res)[2]
  rate=coef(res)[2]
  par(xpd=NA)
  if (showTree) {
    par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
    plot(obsphy)
    axisPhylo(1,backward = T)
  }
  plot(date,ys,xlab='Sampling date',ylab='Root-to-tip distance',xaxs='i',yaxs='i',pch=19,ylim=c(0,max(ys)),xlim=c(ori,max(date,na.rm = T)))
  par(xpd=F)
  abline(res,lwd=2)
  mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p=%.2e',rate,ori,summary(res)$r.squared,summary(res)$coefficients[,4][2]), outer = TRUE, cex = 1.5)
  return(list(rate=rate,ori=ori))
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
