#' Simulation of dated tree using rtree and rexp branch lengths
#' @param nsam Number of leaves in the tree
#' @param dateroot Date of the root
#' @return A simulated dated tree
#' @export
simdatedtree = function(nsam = 20, dateroot = 2000) {
  phy <- rtree(nsam, br = rexp)
  phy$root.time <- 2000
  return(phy)
}

#' Simulation of coalescent dated tree
#' @param dates Vector of dates at which samples are taken
#' @param alpha Coalescent time unit
#' @return A simulated dated tree
#' @export
simcoaltree = function(dates=NA,alpha=10) {
  prob <- 0
  if (is.na(dates[1])) dates=1990:2010
    MySort <- sort(dates,decreasing=TRUE,index.return = TRUE); tim <- MySort$x; ind <- MySort$ix
    n <- length(tim)
    nodes <- cbind(0,ind[1],-Inf)#Start with one node at time -Inf and with the first isolate connected to it
    i <- 2
    while (i <= n) {#Graft branches one by one
      r <- -log(runif(1)) * alpha
      curt <- tim[i];#Current time:start with date of isolate and go back in time until coalescence happens
      fi <- which( nodes[ ,1] < curt ) ;fi<-fi[1]
      if (fi<=nrow(nodes)) for (j in (fi:nrow(nodes)))  {
        if (r > (curt-nodes[j,1]) * (i-j))  {
          prob <- prob + log(1-pexp((curt-nodes[j,1]) * (i-j),alpha^(-1)))
          r <- r-(curt-nodes[j,1]) * (i-j)
          curt <- nodes[j,1]
        } else {
          curt <- curt-r/(i-j);#Found the time for grafting
          prob <- prob + log(dexp(r,alpha^(-1)))
          r <- 0
          break
        }
      }
      if (r>0) stop('This should not happen.\n')
      #Create new node
      a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line
      nodes2=c(curt,ind[i],0)
      if (1<=j-1) nodes2=rbind(nodes[1:(j-1), ],nodes2)
      if (j<=nrow(nodes)) nodes2=rbind(nodes2,nodes[j:nrow(nodes),])
      nodes=unname(nodes2)
      #Now choose on which branch to regraft amongst the branches alive at time curt
      no <- j
      side <- 2
      #prob <- prob + log(1/(nrow(nodes)-j))
      w <- 1 + floor(runif(1) * (nrow(nodes)-j))
      while (w > 0)  {
        no <- no + side-1
        side <- 3-side
        if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  {
          w <- w-1
        }
      }
      nodes[j,3] <- nodes[no,side + 1]
      nodes[no,side + 1] <- n + j
      i <- i + 1
    }
    v=nrow(nodes)-1
    if (nrow(nodes)>2) v=c(v,1:(nrow(nodes)-2))
    nodes=nodes[v,,drop=F]
    m=nodes[,2:3]
    m[which(m>n)]=m[which(m>n)]+1
    nodes[,2:3]=m
    nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes)
    nodes[1:n,1] <- dates

    #Convert into tree
    t=list()
    t$Nnode=n-1
    t$tip.label=as.character(1:n)
    t$edge=matrix(NA,2*n-2,2)
    t$edge.length=rep(NA,n*2-2)
    t$root.time=nodes[n+1,1]
    t$prob=prob
    c=1
    for (i in (n+1):nrow(nodes)) for (j in 2:3) {
      t$edge[c,1]=i
      t$edge[c,2]=nodes[i,j]
      t$edge.length[c]=nodes[nodes[i,j],1]-nodes[i,1]
      c=c+1
    }
    class(t)='phylo'
    return(t)
}

#' Simulation of observed phylogeny given a dated tree
#' @param tree Dated tree
#' @param rate Substitution clock rate
#' @param ratestd Per-branch std on the clock rate (used only by relaxed gamma model)
#' @param model Which model to use (poisson or gamma or relaxedgamma)
#' @return An observed phylogenetic tree
#' @export
simobsphy = function(tree, rate = 10, ratestd = 10, model = 'gamma') {
  obsphy=tree
  obsphy$prob=0
  obsphy$root.time=NULL
  if (model=='poisson') {
    e=unlist(lapply(obsphy$edge.length*rate,function (x) rpois(1,x)))
    obsphy$prob=sum(mapply(dpois,e,rate*obsphy$edge.length,rep(T,length(obsphy$edge.length))))
    obsphy$edge.length=e
  }
  if (model=='gamma') {
    e=unlist(lapply(obsphy$edge.length*rate,function (x) rgamma(1,shape=x,scale=1)))
    obsphy$prob=sum(log(mapply(dgamma,e,rate*obsphy$edge.length,rep(1,length(obsphy$edge.length)))))
    obsphy$edge.length=e
  }
  if (model=='relaxedgamma')  for (i in 1:length(obsphy$edge.length)) {
    l=obsphy$edge.length[i]
    ratevar=ratestd^2
    obsphy$edge.length[i]=rgamma(1,shape=l*rate*rate/(rate+l*ratevar),scale=1+l*ratevar/rate)
    obsphy$prob=obsphy$prob+dgamma(obsphy$edge.length[i],shape=l*rate*rate/(rate+l*ratevar),scale=1+l*ratevar/rate,log=T)
  }
  return(obsphy)
}
