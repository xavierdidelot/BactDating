#' Load output from Gubbins
#' @param prefix Prefix of the ClonalFrameML output files
#' @return Tree representing the output from Gubbins
#' @export
loadGubbins = function(prefix)
{
  tree=read.tree(sprintf('%s.final_tree.tre',prefix))
  tree2=read.tree(sprintf('%s.node_labelled.tre',prefix))
  tree$node.label=tree2$node.label
  tree=unroot(tree)
  n=length(tree$tip.label)
  t=read.table(sprintf('%s.per_branch_statistics.csv',prefix),sep='\t',as.is=T,header=T)
  tree$unrec=rep(NA,length(tree$edge.length))
  for (i in 1:length(tree$edge.length)) {
    j=tree$edge[i,2]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]
    w=which(t[,1]==nam)
    tree$unrec[i]=t[w,10]/t[w,9]
  }
  return(tree)
}
