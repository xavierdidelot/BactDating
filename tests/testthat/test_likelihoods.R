context("Test likelihood functions")

test_that("Poisson likelihood in C++ and R give identical results.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  res=credate(phy,leaves,nbIts=2,model='poisson')
  res2=credate(phy,leaves,nbIts=2,model='poissonR')
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
  res=credate(phy,leaves,nbIts=2,model='gamma')
  res2=credate(phy,leaves,nbIts=2,model='gammaR')
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
  res=credate(phy,leaves,nbIts=2,model='relaxedgamma')
  res2=credate(phy,leaves,nbIts=2,model='relaxedgammaR')
  expect_equal(res$record[1,'likelihood'],res2$record[1,'likelihood'])
})

test_that("Poisson likelihood in C++ and R give identical results on recombinant tree.", {
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
