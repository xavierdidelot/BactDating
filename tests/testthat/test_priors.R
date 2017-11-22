context("Test prior functions")

test_that("Coalescent prior in C++ and R give identical results.", {
  skip_on_cran()
  leaves=2010:2000
  intnodes=2008.5:1999.5
  neg=1.1
  pCpp=coalpriorC(leaves,intnodes,neg)
  pR=coalprior(leaves,intnodes,neg)
  expect_equal(pCpp,pR)
})
