#' Root to tip correlation
#' @param tree Phylogenetic tree
#' @param date Dates of sampling
#' @export
roottotip = function(tree,date)
{
  n=length(date)
  ys=leafDates(tree)
  res=lm(ys~date)
  ori=-coef(res)[1]/coef(res)[2]
  rate=coef(res)[2]
  par(xpd=NA)
  plot(date,ys,xlab='Sampling date',ylab='Root-to-tip distance',xaxs='i',yaxs='i',pch=19,ylim=c(0,max(ys)),xlim=c(ori,max(date)))
  par(xpd=F)
  abline(res,lwd=2)
  text(ori+(max(date)-ori)*0.5,max(ys)*0.9,sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p=%.2e',rate,ori,summary(res)$r.squared,summary(res)$coefficients[,4][2]))
  return(list(rate=rate,ori=ori))
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
  return(dates)
}
