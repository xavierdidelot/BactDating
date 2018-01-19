#' Plotting methods
#' @param x Output from running function credating
#' @param type Type of plot to do. Currently either 'tree' or 'treeCI' or 'trace' or 'treeRoot'
#' @param ... Additional parameters are passed on
#' @return Plot of CreDating results
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
    plot(x$record[,'ratevar'],main='Clock rate variance',type='l',xlab='Sampled iterations',ylab='')
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
}

#' Print function for resCreDating objects
#' @param x output from credate
#' @param ... Passed on to print.phylo
#' @return Print out details of CreDating results
#' @export
print.resCreDating <- function(x, ...)
{
  stopifnot(inherits(x, "resCreDating"))
  cat( 'Phylogenetic tree dated using CreDating\n')
  print(x$tree,...)
  cat(sprintf('Probability of root branch=%.2f\n', x$rootprob))
  for (nam in c('likelihood','prior','rate','ratevar','alpha')) {
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

#' Summary function for resCreDating objects
#' @param object output from credate
#' @param ... Passed on to print.phylo
#' @return Print out details of CreDating results
#' @export
summary.resCreDating <- function(object, ...){
  stopifnot(inherits(object, "resCreDating"))
  cat( 'Phylogenetic tree dated using CreDating\n')
  print(object$tree,...)
  invisible(object)
}
