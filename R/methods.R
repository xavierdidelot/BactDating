#' Plotting methods
#' @param x Output from bactdate
#' @param type Type of plot to do. Currently either 'tree' or 'treeCI' or 'trace' or 'treeRoot' or 'scatter'
#' @param show.axis Whether or not to show the tree axis in mode tree, treeCI or treeRoot
#' @param ... Additional parameters are passed on
#' @return Plot of BactDating results
#' @importFrom grDevices rgb
#' @export
plot.resBactDating = function(x, type='tree', show.axis=T, ...) {

  old.par=par(no.readonly = T)
  if (type=='tree') {
    plot.phylo(x$tree, ...)
    if (show.axis) axisPhylo(backward = F)
  }

  if (type=='treeCI') {
    xl=c(min(x$CI),max(x$CI))-x$tree$root.time
    plot.phylo(x$tree, x.lim=xl,...)
    pre=pretty(xl+x$tree$root.time)
    if (show.axis) axis(1,pre-x$tree$root.time,pre)
#    axisPhylo(backward = F)
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
    plot(x$record[,'likelihood']+x$record[,'prior'],main='Posterior',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'likelihood'],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'prior'],main='Prior',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,Ntip(x$tree)+1],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'mu'],main='Clock rate',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'alpha'],main='Coalescent time unit',type='l',xlab='Sampled iterations',ylab='')
    plot(x$record[,'sigma'],main='Relaxation parameter',type='l',xlab='Sampled iterations',ylab='')
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
    probs=round(probs*100)
    labs=as.numeric(labels(roots)[[1]])
    labs=labs[which(probs>0)]
    probs=probs[which(probs>0)]
    edgelabels(probs,labs,frame='none',adj=c(0.5,-0.1))
    if (show.axis) axisPhylo(backward = F)
  }

  if (type=='scatter') {
#    normed=(date-min(date,na.rm=T))/(max(date,na.rm=T)-min(date,na.rm=T))
#    cols=rgb(ifelse(is.na(normed),0.5,normed),0,1-ifelse(is.na(normed),0.5,normed),ifelse(is.na(normed),1,0.5))
    xs=x$tree$edge.length
    ys=x$tree$subs
    ma=max(xs)*1.05
    rate=mean(x$record[(nrow(x$record)/2):nrow(x$record),'mu'])
    sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
    par(mfrow=c(1,2))
    plot(c(0,ma),c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,max(xs)*1.05),ylim=c(0,max(ys)*1.05))
    par(xpd=F)
    xss=seq(0,ma,0.1)
    plim=0.05
#    lines(xs,xs*rate)
    if (x$model=='poisson') {
      lines(xss,qpois(  plim/2,xss*rate),lty='dashed')
      lines(xss,qpois(1-plim/2,xss*rate),lty='dashed')
      ll=dpois(round(ys),xs*rate,log=T)
    }
    else if (x$model=='negbin') {
      k=rate*rate/sigma/sigma
      theta=sigma*sigma/rate
      lines(xss,qnbinom(  plim/2,size=k,prob=1/(1+theta*xss)),lty='dashed')
      lines(xss,qnbinom(1-plim/2,size=k,prob=1/(1+theta*xss)),lty='dashed')
      ll=dnbinom(round(ys),size=k,prob=1/(1+theta*xs),log=T)
    }
    else if (x$model=='strictgamma') {
      lines(xss,qgamma(  plim/2,shape=xss*rate,scale=1),lty='dashed')
      lines(xss,qgamma(1-plim/2,shape=xss*rate,scale=1),lty='dashed')
      ll=dgamma(ys,shape=xs*rate,scale=1,log=T)
    } else {
      lines(xss,qgamma(  plim/2,shape=xss*rate^2/(rate+xss*sigma^2),scale=1+xss*sigma^2/rate),lty='dashed')
      lines(xss,qgamma(1-plim/2,shape=xss*rate^2/(rate+xss*sigma^2),scale=1+xss*sigma^2/rate),lty='dashed')
      ll=dgamma(ys,shape=xs*rate^2/(rate+xs*sigma^2),scale=1+xs*sigma^2/rate,log=T)
    }
    normed=(ll-min(ll))/(max(ll)-min(ll))
    cols=rgb(normed,0,1-normed)
    base=seq(0,1,0.25)
    legend("topleft",cex=0.5,legend=sprintf('%.2e',base*(max(ll)-min(ll))+min(ll)),pch=19,col=rgb(base,0,1-base))
    points(xs,ys,pch=19,col=cols)
    #w=which(ll<sort(ll)[round(0.01*length(ll))])
    #print(cbind(x$tree$edge[w,2],ll[w]))
    #text(xs[w],ys[w]+0.5,x$tree$edge[w,2],cex=0.5)

    plot(x$tree,show.tip.label = F,edge.color=cols)
    axisPhylo(1,backward = F)

  }
  par(old.par)
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
  for (nam in c('likelihood','prior','mu','sigma','alpha')) {
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

#' Convert to coda mcmc format
#' @param x Output from bactdate
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return mcmc object from coda package
#' @export
as.mcmc.resBactDating <- function(x,burnin=0.5) {
  record=x$record
  record=record[max(1,round(nrow(record)*burnin)):nrow(record),]
  mat=cbind(
    record[,'mu'],
    record[,'sigma'],
    record[,'alpha'])
  colnames(mat)<-c('mu','sigma','alpha')
  return(coda::as.mcmc(mat))
}

#' Convert to treedata format from ggtree package
#' @param x Output from bactdate
#' @return treedata object from ggtree package
#' @importFrom methods new
#' @importFrom dplyr tbl_df
#' @export
as.treedata.resBactDating <- function(x) {
  t=x$tree
  h=node.depth.edgelength(t)
  meta=data.frame(node=sprintf('%5d',1:(Ntip(t)+Nnode(t))))
  meta$length_0.95_HPD=list(c(-1,1))
  for (i in 1:nrow(meta)) {
    w=which(t$edge[,2]==i)
    if (length(w)==0) len=0 else len=t$edge.length[w]
    interval=x$CI[i,]-t$root.time-h[i]+len
    meta$length_0.95_HPD[i]=list(interval)
  }
  obj=new("treedata",treetext='',phylo=t,data=tbl_df(as.data.frame(meta)),file='')
  return(obj)
  #can then plot with:
  #ggtree(obj) + geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2)
}

#' Modified drop.tip() from 'ape' package that also updates the recombination
#' proportions added when loaded from CFML or Gubbins
#' @param phy Tree from loadCFML() or loadGubbins()
#' @param tip a vector of mode numeric or character specifying the tips to delete.
#' @param trim.internal a logical specifying whether to delete the corresponding internal branches.
#' @param subtree a logical specifying whether to output in the tree how many tips have been deleted and where.
#' @param root.edge an integer giving the number of internal branches to be used to build the new root edge. This has no effect if \code{trim.internal FALSE}.
#' @param rooted a logical indicating whether the tree must be treated as rooted or not. This allows to force the tree to be considered as unrooted.
#' @param collapse.singles a logical specifying whether to delete the internal nodes of degree 2.
#' @param interactive if \code{TRUE} the user is asked to select the tips or the node by clicking on the tree which must be plotted.
#' @return tree with rec data
#' @importFrom methods new
#' @export
drop.tip.useRec = function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0,
          rooted = is.rooted(phy), collapse.singles = TRUE, interactive = FALSE)
{
  if (!inherits(phy, "phylo"))
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)
  if (interactive) {
    cat("Left-click close to the tips you want to drop; right-click when finished...\n")
    xy <- locator()
    nToDrop <- length(xy$x)
    tip <- integer(nToDrop)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    for (i in 1:nToDrop) {
      d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
      tip[i] <- which.min(d)
    }
  }
  else {
    if (is.character(tip))
      tip <- which(phy$tip.label %in% tip)
  }
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip))
    return(phy)
  if (length(tip) == Ntip) {
    if (Nnode(phy) < 3 || trim.internal) {
      warning("drop all tips of the tree: returning NULL")
      return(NULL)
    }
  }
  wbl <- !is.null(phy$edge.length)
  wbl.unrec <- !is.null(phy$unrec)
  if (length(tip) == Ntip - 1 && trim.internal) {
    i <- which(phy$edge[, 2] == (1:Ntip)[-tip])
    res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i,
                                                                             2]], Nnode = 1L)
    class(res) <- "phylo"
    if (wbl)
      res$edge.length <- phy$edge.length[i]
    if (wbl.unrec)
      res$unrec <- phy$unrec[i]
    if (!is.null(phy$node.label))
      res$node.label <- phy$node.label[phy$edge[i, 1] -
                                         Ntip]
    return(res)
  }
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  if (subtree) {
    trim.internal <- TRUE
    tr <- reorder(phy, "postorder")
    N <- ape::node_depth(as.integer(Ntip), as.integer(tr$edge[,1]), as.integer(tr$edge[, 2]), as.integer(Nedge),double(Ntip + Nnode), 1L)[[5]]
  }
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  if (trim.internal) {
    ints <- edge2 > Ntip
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep
      if (!sum(sel))
        break
      keep[sel] <- FALSE
    }
    if (subtree) {
      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
      keep[subt] <- TRUE
    }
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          degree <- tabulate(edge1[keep])
          if (degree[NEWROOT] > 1)
            break
        }
        keep[j] <- FALSE
        if (length(j) > root.edge)
          j <- 1:root.edge
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge))
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  if (!root.edge)
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl)
    phy$edge.length <- phy$edge.length[keep]
  if (wbl.unrec)
    phy$unrec <- phy$unrec[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  if (subtree) {
    i <- which(tip %in% oldNo.ofNewTips)
    if (length(i)) {
      phy$tip.label[tip[i]] <- "[1_tip]"
      tip <- tip[-i]
    }
  }
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  if (length(tip))
    phy$tip.label <- phy$tip.label[-tip]
  if (subtree || !trim.internal) {
    node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
    new.tip.label <- if (!length(node2tip)) {
      character(0)
    }
    else if (subtree) {
      paste("[", N[node2tip], "_tips]", sep = "")
    }
    else {
      if (is.null(phy$node.label))
        rep("NA", length(node2tip))
      else phy$node.label[node2tip - Ntip]
    }
    phy$tip.label <- c(phy$tip.label, new.tip.label)
  }
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label))
    phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
  if (collapse.singles)
    phy <- collapse.singles.useRec(phy)
  phy
}

#' Modified collapse.singles() from package ape which deletes the single nodes (i.e., with a single
#' descendant) in a tree, respecting the rec function. Internal, used by drop.tip.rec
#' @param tree an object of class \code{"phylo"}.
#' @param root.edge whether to get the singleton edges from the root until the first bifurcating node and put them as \code{root.edge} of the returned tree. By default, this is ignored or if the tree has no edge lengths.
#' @return tree with rec data
#' @importFrom methods new
collapse.singles.useRec <- function (tree, root.edge = FALSE)
{
  n <- length(tree$tip.label)
  tree <- reorder(tree)
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  tab <- tabulate(e1)
  if (all(tab[-c(1:n)] > 1))
    return(tree)
  if (is.null(tree$edge.length)) {
    root.edge <- FALSE
    wbl <- FALSE
  }
  else {
    wbl <- TRUE
    el <- tree$edge.length
  }
  if (is.null(tree$unrec)) {
    wbl.unrec <- FALSE
  }
  else {
    wbl.unrec <- TRUE
    unrec <- tree$unrec
  }
  if (root.edge)
    ROOTEDGE <- 0
  ROOT <- n + 1L
  while (tab[ROOT] == 1) {
    i <- which(e1 == ROOT)
    ROOT <- e2[i]
    if (wbl) {
      if (root.edge)
        ROOTEDGE <- ROOTEDGE + el[i]
      el <- el[-i]
      unrec <- unrec[-i]
    }
    e1 <- e1[-i]
    e2 <- e2[-i]
  }
  singles <- which(tabulate(e1) == 1)
  if (length(singles) > 0) {
    ii <- sort(match(singles, e1), decreasing = TRUE)
    jj <- match(e1[ii], e2)
    for (i in 1:length(singles)) {
      e2[jj[i]] <- e2[ii[i]]
      if (wbl.unrec)
        unrec[jj[i]] <- (el[jj[i]]*unrec[jj[i]] + el[ii[i]]*unrec[ii[i]])/(el[jj[i]] + el[ii[i]]) # correct?
      if (wbl)
        el[jj[i]] <- el[jj[i]] + el[ii[i]]
    }
    e1 <- e1[-ii]
    e2 <- e2[-ii]
    if (wbl)
      el <- el[-ii]
    if (wbl.unrec)
      unrec <- unrec[-ii]
  }
  Nnode <- length(e1) - n + 1L
  oldnodes <- unique(e1)
  if (!is.null(tree$node.label))
    tree$node.label <- tree$node.label[oldnodes - n]
  newNb <- integer(max(oldnodes))
  newNb[ROOT] <- n + 1L
  sndcol <- e2 > n
  e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
  e1 <- newNb[e1]
  tree$edge <- cbind(e1, e2, deparse.level = 0)
  tree$Nnode <- Nnode
  if (wbl) {
    if (root.edge)
      tree$root.edge <- ROOTEDGE
    tree$edge.length <- el
    tree$unrec <- unrec
  }
  tree
}

