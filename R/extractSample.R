#' Extract a sample of trees
#' @param x Output from bactdate
#' @param sampleSize Number of trees to sample
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return list of trees
#' @export
extractSample <- function(x,sampleSize=500,burnin=0.5) {
  record=x$record
  record=record[max(1,round(nrow(record)*burnin)):nrow(record),]
  out=vector("list", sampleSize)
  rows=round(seq(1,nrow(record),length.out=sampleSize))
  for (j in 1:sampleSize) {
    row=rows[j]
    curRec=record[row,]
    tree=x$inputtree
    for (i in 1:nrow(tree$edge)) {
      tree$edge[i,1]=curRec[Ntip(tree)+Nnode(tree)+tree$edge[i,2]]
      tree$edge.length[i] = curRec[tree$edge[i, 2]] - curRec[tree$edge[i, 1]]
    }
    tree$root.time=curRec[Ntip(tree)+1]
    out[[j]]=tree
  }
  return(out)
}
