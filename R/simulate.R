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
#' @param model Which model to use (poisson or gamma)
#' @return The observed phylogenetic tree
#' @export
simobsphy = function(tree, rate = 10, model = 'gamma') {
  obsphy=tree
  obsphy$root.time=NULL
  if (model=='poisson') obsphy$edge.length=unlist(lapply(obsphy$edge.length*rate,function (x) rpois(1,x)))
  if (model=='gamma')   obsphy$edge.length=unlist(lapply(obsphy$edge.length*rate,function (x) rgamma(1,x,1)))
  return(obsphy)
}

