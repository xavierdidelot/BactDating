#' Clustered permutation test
#' @param tree Phylogenetic tree
#' @param dates Dates of sampling
#' @return Results of root-to-tip analysis after clustering
#' @export
clusteredTest = function(tree,dates)
{
  #Rerranging of dates, if needed
  if (!is.null(names(dates))) dates=findDates(tree,dates)
  #Remove missing dates
  mis=which(is.na(dates))
  if (length(mis)>0) {dates=dates[-mis];tree=drop.tip(tree,mis)}

  curdates=dates
  curtree=tree
  window=0
  while (T) {
    tim=dist(curdates)
    gen=cophenetic.phylo(curtree)
    res=suppressMessages(vegan::mantel(gen,tim))
    if (res$signif>0.01) break
    window=window+1
    n=Ntip(curtree)
    drop=rep(0,n)
    order=sample(1:n,n,replace = F)
    for (i in 2:n)
      if (min(abs(curdates[which(drop[order[1:(i-1)]]==0)]-curdates[order[i]]))<window) drop[order[i]]=1
    curtree=drop.tip(curtree,which(drop==1))
    curdates=curdates[which(drop==0)]
  }

  roottotip(curtree,curdates)
}
