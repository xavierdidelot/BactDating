#' Main function
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param rate Rate per genome (not per site)
#' @return Dating results
#' @export
credating = function(tree,date,rate)
{
  n=length(tree$tip.label)

  #Create table of nodes (col1=name,col2=expected mutation on branch above,col3=date,col4=father)
  tab=matrix(NA,n+tree$Nnode,4)
  for (r in 1:(nrow(tree$edge))) {
    i=tree$edge[r,2]
    tab[i,1]=i
    tab[i,4]=tree$edge[r,1]
    tab[i,2]=round(tree$edge.length[r])
    if (i<=n) tab[i,3]=date[i]
  }

  #Starting point
  tree=reorder.phylo(tree,'postorder')
  for (r in 1:(nrow(tree$edge))) {
    i=tree$edge[r,2]
    if (i<=n) next
    children=which(tab[,4]==i)
    if (length(children)!=0) tab[i,3]=min(tab[children,3])-1
  }
  i=n+1;children=which(tab[,4]==i);tab[i,1]=i;tab[i,3]=min(tab[children,3])-1

  #MCMC
  l=likelihood(tab,rate,n)
  nbIts=1000#TODO
  record=matrix(NA,nbIts/10,2+nrow(tab))
  for (i in 1:nbIts) {
    #Record
    if (i%%10==0) {
      record[i/10,nrow(tab)+1]=l
      record[i/10,nrow(tab)+2]=rate
      for (j in 1:nrow(tab)) record[i/10,j]=tab[j,3]
    }

    #Gibbs move to update rate which assumes Exp(1) prior on rate
    #lengths=tab[-(n+1),3]-tab[tab[-(n+1),4],3]
    #muts=tab[-(n+1),2]
    #rate=rgamma(1,1+sum(muts),sum(lengths))
    #l=likelihood(tab,rate,n)

    #MH to update dates
    for (j in (n+1):nrow(tab)) {
      old=tab[j,3]
      tab[j,3]=old+runif(1)-0.5
      l2=likelihood(tab,rate,n)
      if (log(runif(1))<l2-l) l=l2
      else tab[j,3]=old
    }
  }

  #Plot
#  par(mfrow=c(2,2))
#  plot(record[,nrow(tab)+1],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
#  plot(record[,n+1],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
#  plot(record[,nrow(tab)+2],main='Rate',type='l',xlab='Sampled iterations',ylab='')
#  par(mfrow=c(1,1))

  #Output
  meanRec=colMeans(record[(nrow(record)/2):nrow(record),])
  for (i in 1:nrow(tree$edge)) tree$edge.length[i]=meanRec[tree$edge[i,2]]-meanRec[tree$edge[i,1]]
  return(list(tree=tree,record=record,rootdate=meanRec[n+1]))
}

#' Likelihood function
#' @return log-likelihood
likelihood = function(tab,rate,n) {
if (rate<0) return(-Inf)
#n=ceiling(nrow(tab)/2)
t2=tab[-(n+1),]
lengths=t2[,3]-tab[t2[,4],3]
if (min(lengths)<0) return(-Inf)
muts=t2[,2]
return(sum(-lengths*rate+muts*log(lengths*rate)))
}
