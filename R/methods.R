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
    par(mfrow=c(2,4))
    plot(x$record[,'likelihood'],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'prior'],main='Prior',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,Ntip(x$tree)+1],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'rate'],main='Clock rate',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'neg'],main='Coalescent rate',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'ratevar'],main='Clock rate variance',type='l',xlab='Sampled iterations',ylab='')
    v=x$record[,'root']
    u=sort(unique(v))
    for (i in 1:length(u)) v[which(x$record[,'root']==u[i])]=i
    plot(v,main='Root branch',type='l',xlab='Sampled iterations',ylab='',yaxt='n')
    axis(2,1:length(u),u)
  }

  if (type=='treeRoot') {
    roots = table(x$record[floor(nrow(x$record) / 2):nrow(x$record),'root'])
    plot.phylo(x$inputtree, ...)
    probs=as.vector(roots)
    probs=probs/sum(probs)
    probs=round(probs*1000)/10
    edgelabels(probs,as.numeric(labels(roots)[[1]]))
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
  cat(paste( 'Mean date of root:', x$rootdate, '\n'))
  cat(paste('Probability of root branch:', x$rootprob, '\n'))
  cat(paste('Mean clock rate:',mean(x$record[round(nrow(x$record)/2):nrow(x$record),ncol(x$record)-2])))
  invisible(x)
}

#' Summary function for resCreDating objects
#' @param object output from credate
#' @param ... Passed on to print.phylo
#' @return Print out details of CreDating results
#' @export
summary.resCreDating <- function(object, ...){
  print.resCreDating(object,...)
}
