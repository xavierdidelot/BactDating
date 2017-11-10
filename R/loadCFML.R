#' Load output from ClonalFrameML
#' @param prefix Prefix of the ClonalFrameML output files
#' @importFrom utils read.table
#' @export
loadCFML = function(prefix)
{
  tree=read.tree(sprintf('%s.labelled_tree.newick',prefix))
  tree=unroot(tree)
  n=length(tree$tip.label)
  imports=read.table(sprintf('%s.importation_status.txt',prefix),header=T,as.is=T,sep="\t")
  params=read.table(sprintf('%s.em.txt',prefix),header=T,as.is=T,sep="\t")
  L=length(scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet = T))
  tree$unrec=rep(NA,length(tree$edge.length))
  for (i in 1:length(tree$edge.length)) {
    j=tree$edge[i,2]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]
    w=which(params[,1]==nam)
    except=F
    if (length(w)==0 && except==F) {w=nrow(param);except=T}
    if (length(w)==0 && except==T) {warning('Error in ClonalFrameML output file.')}
    tree$unrec[i]=params[w,5]/L
  }
  if (max(tree$unrec)>0) tree$unrec=tree$unrec/max(tree$unrec)
  return(tree)
}
