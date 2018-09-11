#' Load output from ClonalFrameML
#' @param prefix Prefix of the ClonalFrameML output files
#' @param priorMeanM Mean of prior on M used when running ClonalFrameML
#' @param priorSdM Sd of prior on M used when running ClonalFrameML
#' @importFrom utils read.table
#' @return Tree representing the output from ClonalFrameML
#' @export
loadCFML = function(prefix,priorMeanM=0.0001,priorSdM=0.0001)
{
  priorA=priorMeanM^2/priorSdM^2
  priorB=priorMeanM/priorSdM^2
  tree=read.tree(sprintf('%s.labelled_tree.newick',prefix))
  tree=unroot(tree)
  n=length(tree$tip.label)
  imports=read.table(sprintf('%s.importation_status.txt',prefix),header=T,as.is=T,sep="\t",comment.char="")
  params=read.table(sprintf('%s.em.txt',prefix),header=T,as.is=T,sep="\t",comment.char="")
  L=length(scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet = T))
  tree$edge.length=tree$edge.length*L
  tree$unrec=rep(NA,length(tree$edge.length))
  for (i in 1:length(tree$edge.length)) {
    j=tree$edge[i,2]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]
    w=which(params[,1]==nam)
    except=F
    if (length(w)==0 && except==F) {w=nrow(params);except=T}
    if (length(w)==0 && except==T) {warning('Error in ClonalFrameML output file.')}
    tree$unrec[i]=(params[w,5]-priorB)/L
    if (tree$unrec[i]>1) {warning('Error in ClonalFrameML output file, are you sure the prior is correct?')}
  }
  return(tree)
}
