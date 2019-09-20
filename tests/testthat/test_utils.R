#Test utility functions
context("Test utility functions")

test_that("The function changeinorderedvec behaves as expected.", {
  v=100:1
  v=v+0.1
  for (i in 1:100) BactDating:::changeinorderedvec(v,v[sample.int(100,1)],runif(1)*100)
  expect_false(is.unsorted(rev(v)))
})

test_that("The function as.treedata.resBactDating works.", {
  set.seed(0)
  res=bactdate(ape::rtree(10),rep(2000,10),nbIts=100)
  l=as.treedata.resBactDating(res)
  expect_is(l,'list')
})

test_that("The function drop.tip.useRec works.", {
  set.seed(0)
  leaves=2010:2000
  tree=simcoaltree(leaves,5.1)
  phy=simobsphy(tree)
  phy$unrec=runif(length(phy$edge.length))
  expect_silent(initRoot(phy,2010:2020,useRec = T))
  expect_is(drop.tip.useRec(phy,'1'),'phylo')
})

test_that("loadCFML and loadGubbins return errors if missing files.", {
  expect_error(expect_warning(loadCFML('nothing')))
  expect_error(expect_warning(loadGubbins('nothing')))
})
