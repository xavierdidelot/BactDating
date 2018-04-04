#' Plotting methods
#' @param x Output from running function credating
#' @param type Type of plot to do. Currently either 'tree' or 'treeCI' or 'trace' or 'treeRoot' or 'scatter'
#' @param ... Additional parameters are passed on
#' @return Plot of BactDating results
#' @importFrom grDevices rgb
#' @export
plot.resBactDating = function(x, type='tree', ...) {

  if (type=='tree') {
    plot.phylo(x$tree, ...)
    axisPhylo(backward = F)
  }

  if (type=='treeCI') {
    plot.phylo(x$tree, x.lim=c(min(x$CI),max(x$CI))-x$tree$root.time,...)
    axisPhylo(backward = F)
    obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    transblue=rgb(0,0,1,0.4)
    transred =rgb(1,0,0,0.4)
    for(i in 1:(Nnode(x$tree)+Ntip(x$tree)))
      if (x$CI[i,1]!=x$CI[i,2])
        lines(x=c(x$CI[i,1],x$CI[i,2])-x$tree$root.time,
              y=rep(obj$yy[i],2),lwd=5,lend=0,
              col=ifelse(i<=Ntip(x$tree),transred,transblue))
    #points(obj$xx[1:x$tree$Nnode+Ntip(x$tree)],
    #       obj$yy[1:x$tree$Nnode+Ntip(x$tree)],pch=19,col="blue",
    #       cex=1)
  }

  if (type=='trace') {
    par(mfrow=c(2,4))
    plot(x$record[,'likelihood'],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'prior'],main='Prior',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,Ntip(x$tree)+1],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'rate'],main='Clock rate',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'alpha'],main='Coalescent time unit',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'ratestd'],main='Clock rate std',type='l',xlab='Sampled iterations',ylab='')
    v=x$record[,'root']
    u=sort(unique(v))
    for (i in 1:length(u)) v[which(x$record[,'root']==u[i])]=i
    plot(v,main='Root branch',type='p',pch=20,cex=0.2,xlab='Sampled iterations',ylab='',yaxt='n')
    axis(2,1:length(u),u)
  }

  if (type=='treeRoot') {
    roots = table(x$record[floor(nrow(x$record) / 2):nrow(x$record),'root'])
    plot.phylo(x$inputtree, ...)
    probs=as.vector(roots)
    probs=probs/sum(probs)
    probs=round(probs*1000)/10
    edgelabels(probs,as.numeric(labels(roots)[[1]]),frame='none',adj=c(0.5,-0.1))
    axisPhylo(backward = F)
  }

  if (type=='scatter') {
#    normed=(date-min(date,na.rm=T))/(max(date,na.rm=T)-min(date,na.rm=T))
#    cols=rgb(ifelse(is.na(normed),0.5,normed),0,1-ifelse(is.na(normed),0.5,normed),ifelse(is.na(normed),1,0.5))
    xs=x$tree$edge.length
    ys=x$tree$subs
    ma=max(xs)*1.05
    rate=mean(x$record[(nrow(x$record)/2):nrow(x$record),'rate'])
    ratestd=mean(x$record[(nrow(x$record)/2):nrow(x$record),'ratestd'])
    par(mfrow=c(1,2))
    plot(c(0,ma),c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,max(xs)*1.05),ylim=c(0,max(ys)*1.05))
    par(xpd=F)
    ori=0
    xss=seq(ori,ma,0.1)
    plim=0.01
#    lines(xs,xs*rate)
    if (x$model=='poisson') {
      lines(xss,qpois(  plim/2,(xss-ori)*rate),lty='dashed')
      lines(xss,qpois(1-plim/2,(xss-ori)*rate),lty='dashed')
      ll=dpois(round(ys),xs*rate,log=T)
    }
    else if (x$model=='gamma') {
      lines(xss,qgamma(  plim/2,shape=(xss-ori)*rate,scale=1),lty='dashed')
      lines(xss,qgamma(1-plim/2,shape=(xss-ori)*rate,scale=1),lty='dashed')
      ll=dgamma(ys,shape=xs*rate,scale=1,log=T)
    } else {
      lines(xss,qgamma(  plim/2,shape=(xss-ori)*rate^2/(rate+(xss-ori)*ratestd^2),scale=1+(xss-ori)*ratestd^2/rate),lty='dashed')
      lines(xss,qgamma(1-plim/2,shape=(xss-ori)*rate^2/(rate+(xss-ori)*ratestd^2),scale=1+(xss-ori)*ratestd^2/rate),lty='dashed')
      ll=dgamma(ys,shape=xs*rate^2/(rate+xs*ratestd^2),scale=1+xs*ratestd^2/rate,log=T)
    }
    normed=(ll-min(ll))/(max(ll)-min(ll))
    cols=rgb(normed,0,1-normed)
    points(xs,ys,pch=19,col=cols)
    base=seq(0,1,0.25)
    legend("topleft",cex=0.5,legend=sprintf('%.2e',base*(max(ll)-min(ll))+min(ll)),pch=19,col=rgb(base,0,1-base))

    plot(x$tree,show.tip.label = F,edge.color=cols)
    axisPhylo(1,backward = F)

  }
}

#' Print function for resBactDating objects
#' @param x output from bactdate
#' @param ... Passed on to print.phylo
#' @return Print out details of BactDating results
#' @export
print.resBactDating <- function(x, ...)
{
  stopifnot(inherits(x, "resBactDating"))
  cat( 'Phylogenetic tree dated using BactDating\n')
  print(x$tree,...)
  cat(sprintf('Probability of root branch=%.2f\n', x$rootprob))
  for (nam in c('likelihood','prior','rate','ratestd','alpha')) {
    v=x$record[,nam]
    v=v[(1+length(v)/2):length(v)]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    cat(sprintf('%s=%.2e [%.2e;%.2e]\n',nam,vals[1],vals[2],vals[3]))
  }
  v=x$record[,Ntip(x$tree)+1]
  v=sort(v[(1+length(v)/2):length(v)])
  vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
  cat(sprintf('Root date=%.2f [%.2f;%.2f]\n',vals[1],vals[2],vals[3]))
  invisible(x)
}

#' Summary function for resBactDating objects
#' @param object output from bactdate
#' @param ... Passed on to print.phylo
#' @return Print out details of BactDating results
#' @export
summary.resBactDating <- function(object, ...){
  stopifnot(inherits(object, "resBactDating"))
  cat( 'Phylogenetic tree dated using BactDating\n')
  print(object$tree,...)
  invisible(object)
}
