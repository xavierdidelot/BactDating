#' Load output from Gubbins
#' @param prefix Prefix of the ClonalFrameML output files
#' @return Tree representing the output from Gubbins
#' @export
loadGubbins = function(prefix)
{
  tree=read.tree(sprintf('%s.final_tree.tre',prefix))
  nametree2=sprintf('%s.node_labelled.tre',prefix)
  if (!file.exists(nametree2)) nametree2=sprintf('%s.node_labelled.final_tree.tre',prefix)
  tree2=read.tree(nametree2)
  tree$node.label=tree2$node.label
  tree=unroot(tree)
  n=length(tree$tip.label)
  t=utils::read.table(sprintf('%s.per_branch_statistics.csv',prefix),sep='\t',as.is=T,header=T,comment.char="")
  tree$unrec=rep(NA,length(tree$edge.length))
  for (i in 1:length(tree$edge.length)) {
    j=tree$edge[i,2]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]
    w=which(t[,1]==nam)
    j=tree$edge[i,1]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]
    w2=which(t[,1]==nam)
    tree$unrec[i]=1-(t[w,6]-t[w2,6])/max(t[,9])
  }
  return(tree)
}
