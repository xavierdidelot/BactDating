#' Load output from ClonalFrameML
#' @param prefix Prefix of the ClonalFrameML output files
#' @param priorMeanM Mean of prior on M used when running ClonalFrameML
#' @param priorSdM Sd of prior on M used when running ClonalFrameML
#' @return Tree representing the output from ClonalFrameML
#' @export
loadCFML = function(prefix,priorMeanM=0.0001,priorSdM=0.0001)
{
  priorA=priorMeanM^2/priorSdM^2
  priorB=priorMeanM/priorSdM^2
  tree=read.tree(sprintf('%s.labelled_tree.newick',prefix))
  tree=unroot(tree)
  n=length(tree$tip.label)
  imports=utils::read.table(sprintf('%s.importation_status.txt',prefix),header=T,as.is=T,sep="\t",comment.char="")
  params=utils::read.table(sprintf('%s.em.txt',prefix),header=T,as.is=T,sep="\t",comment.char="")
  L=length(scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet = T))

  tree$edge.length=tree$edge.length*L
  tree$unrec=rep(NA,length(tree$edge.length))
  for (i in 1:length(tree$edge.length)) {
    j=tree$edge[i,2]
    if (j<=n) nam=tree$tip.label[j] else nam=tree$node.label[j-n]

    if (ncol(params==5)) {
      #default em option was used
      w=which(params[,1]==nam)
      if (length(w)==0) {
        #Branch of minimal length
        tree$unrec[i]=1
      } else {
        tree$unrec[i]=(params[w,5]-priorB)/L
        if (tree$unrec[i]>1 || tree$unrec[i]<0) warning('Error in ClonalFrameML output file, are you sure the prior is correct?')
      }
    }

    if (ncol(params==6)) {
      #embranch was used
      clonal=L
      w=which(imports[,1]==nam)
      clonal=clonal-sum(imports[w,3]-imports[w,2]+1)
      tree$unrec[i]=clonal/L
    }
  }
  return(tree)
}
