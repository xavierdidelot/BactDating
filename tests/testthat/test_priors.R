context("Test prior functions")

test_that("Coalescent prior in C++ gives expected result on small example.", {
  a=coalpriorC(c(2000,2000,2000),c(1999,1998),2)
  b=dexp(1,3/2,T)+dexp(1,1/2,T)+log(1/3)
  expect_equal(a,b)
})

test_that("Coalescent prior in C++ and R give identical results.", {
  set.seed(0)
  leaves=2010:2000
  intnodes=2008.5:1999.5
  neg=1.1
  pCpp=coalpriorC(leaves,intnodes,neg)
  pR=coalprior(leaves,intnodes,neg)
  expect_equal(pCpp,pR)

  leaves=2010:1910
  tree=simcoaltree(leaves,5.1)
  intnodes=nodeDates(tree)
  pCpp=coalpriorC(leaves,intnodes,neg)
  pR=coalprior(leaves,intnodes,neg)
  expect_equal(pCpp,pR)
})

test_that("Probabilities of a coalescent tree are the same when simulating and evaluating.",{
  set.seed(0)
  leaves=2010:1910
  tree=simcoaltree(leaves,1.1)
  intnodes=nodeDates(tree)
  pCpp=coalpriorC(leaves,intnodes,1.1)
  expect_equal(pCpp,tree$prob)
})
