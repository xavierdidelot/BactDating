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
#' @param neg Reciprocal of coalescent rate
#' @return A simulated dated tree
#' @export
simcoaltree = function(dates=NA,neg=10) {
  if (is.na(dates[1])) dates=1990:2010
    MySort <- sort(dates,decreasing=TRUE,index.return = TRUE); tim <- MySort$x; ind <- MySort$ix
    n <- length(tim)
    nodes <- cbind(0,ind[1],0)#Start with one node at time 0 and with the first isolate connected to it
    i <- 2
    while (i <= n) {#Graft branches one by one
      r <- -log(runif(1)) * neg
      curt <- tim[i];#Current time:start with date of isolate and go back in time until coalescence happens
      fi <- which( nodes[ ,1] < curt ) ;fi<-fi[1]
      if (fi<=nrow(nodes)) for (j in (fi:nrow(nodes)))  {
        if (r > (curt-nodes[j,1]) * (i-j))  {
          r <- r-(curt-nodes[j,1]) * (i-j)
          curt <- nodes[j,1]
        } else {
          curt <- curt-r/(i-j);#Found the time for grafting
          r <- 0
          break
        }
      }
      if (r>0) {next}
      #Create new node
      a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line
      nodes2=c(curt,ind[i],0)
      if (1<=j-1) nodes2=rbind(nodes[1:(j-1), ],nodes2)
      if (j<=nrow(nodes)) nodes2=rbind(nodes2,nodes[j:nrow(nodes),])
      nodes=unname(nodes2)
      #Now choose on which branch to regraft amongst the branches alive at time curt
      no <- j
      side <- 2
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
    nodes=nodes[c(nrow(nodes)-1,1:(nrow(nodes)-2)),]
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
#' @param ratevar Per-branch variance on the clock rate (used only by relaxed gamma model)
#' @param model Which model to use (poisson or gamma or gamma2)
#' @return The observed phylogenetic tree
#' @export
simobsphy = function(tree, rate = 10, ratevar = 0, model = 'gamma') {
  obsphy=tree
  obsphy$root.time=NULL
  if (model=='poisson') obsphy$edge.length=unlist(lapply(obsphy$edge.length*rate,function (x) rpois(1,x)))
  if (model=='gamma')   obsphy$edge.length=unlist(lapply(obsphy$edge.length*rate,function (x) rgamma(1,shape=x,scale=1)))
  if (model=='relaxedgamma')  for (i in 1:length(obsphy$edge.length)) {
    l=obsphy$edge.length[i]
    obsphy$edge.length[i]=rgamma(1,shape=l*rate*rate/(rate+l*ratevar),scale=1+l*ratevar/rate)
  }
  return(obsphy)
}
