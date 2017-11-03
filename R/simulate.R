#' Simulation of dated tree
#' @param nsam Number of leaves in the tree
#' @param dateroot Date of the root
#' @return A simulated dated tree
#' @export
simdatedtree = function(nsam = 20, dateroot = 2000) {
  phy <- rtree(nsam, br = rexp)
  phy$root.time <- 2000
  return(phy)
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
