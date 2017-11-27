context("Test likelihood functions")

makeTab=function(timedTree,phyTree) {
  n=length(timedTree$tip.label)
  d1=leafDates(timedTree)
  d2=nodeDates(timedTree)
  tab=matrix(NA,n*2-1,4)
  tab[,3]=c(d1,d2)
  for (i in c(1:n,(n+2):(2*n-1))) {
    e=which(timedTree$edge[,2]==i)
    tab[i,4]=timedTree$edge[e,1]
    tab[i,2]=phyTree$edge.length[e]}
  return(tab)
}

test_that("Likelihood is equal to probability of simulation.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree,model='poisson',rate=5.5)
  expect_equal(likelihoodPoissonC(makeTab(tree,phy),5.5),phy$prob)
  phy=simobsphy(tree,model='gamma',rate=5.5)
  expect_equal(likelihoodGammaC(makeTab(tree,phy),5.5),phy$prob)
  phy=simobsphy(tree,model='relaxedgamma',rate=5.5,ratevar=4.1)
  expect_equal(likelihoodRelaxedgammaC(makeTab(tree,phy),5.5,4.1),phy$prob)
})

test_that("Likelihood in C++ and R give identical results.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  tab=makeTab(tree,phy)
  expect_equal(likelihoodPoisson(tab,5.5),likelihoodPoissonC(tab,5.5))
  expect_equal(likelihoodGamma(tab,5.5),likelihoodGammaC(tab,5.5))
  expect_equal(likelihoodRelaxedgamma(tab,5.5,4.1),likelihoodRelaxedgammaC(tab,5.5,4.1))
})

test_that("Likelihood of strict gamma is equal to relaxed gamma with no variance.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  tab=makeTab(tree,phy)
  expect_equal(likelihoodRelaxedgammaC(tab,5.5,0),likelihoodGammaC(tab,5.5))
})

test_that("Short MCMC runs give same results for C++ and R likelihoods.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  set.seed(0)
  res=credate(phy,leaves,nbIts=10,model='poisson')
  set.seed(0)
  res2=credate(phy,leaves,nbIts=10,model='poissonR')
  expect_equal(res$record[10,'likelihood'],res2$record[10,'likelihood'])
  set.seed(0)
  res=credate(phy,leaves,nbIts=10,model='gamma')
  set.seed(0)
  res2=credate(phy,leaves,nbIts=10,model='gammaR')
  expect_equal(res$record[10,'likelihood'],res2$record[10,'likelihood'])
  set.seed(0)
  res=credate(phy,leaves,nbIts=10,model='relaxedgamma')
  set.seed(0)
  res2=credate(phy,leaves,nbIts=10,model='relaxedgammaR')
  expect_equal(res$record[10,'likelihood'],res2$record[10,'likelihood'])
})

test_that("Likelihood in C++ and R give identical results on recombinant tree.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  phy$unrec=runif(length(phy$edge.length))
  res=credate(phy,leaves,nbIts=2,model='poisson',useRec=T)
  res2=credate(phy,leaves,nbIts=2,model='poissonR',useRec=T)
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
  res=credate(phy,leaves,nbIts=2,model='gamma',useRec=T)
  res2=credate(phy,leaves,nbIts=2,model='gammaR',useRec=T)
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
  res=credate(phy,leaves,nbIts=2,model='relaxedgamma',useRec=T)
  res2=credate(phy,leaves,nbIts=2,model='relaxedgammaR',useRec=T)
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
})
