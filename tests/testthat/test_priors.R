context("Test prior functions")

test_that("Coalescent prior in C++ and R give identical results.", {
  skip_on_cran()
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
