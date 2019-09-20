#Test all is running without error
context("Test running without error")

test_that("Basic functions are running without error.", {
  set.seed(0)
  expect_silent(tree<-simdatedtree(10,2005))
  expect_silent(tree<-simcoaltree(2010:2020,5.1))
  expect_equal(allDates(tree)[1:11],2010:2020)
  d=2010:2020;names(d)=1:11
  expect_equal(findDates(tree,d),2010:2020)
  d[5]=NA
  expect_silent(phy<-simobsphy(tree))
  expect_silent(phy<-initRoot(phy,d))
  expect_silent(res<-roottotip(phy,d))
  expect_silent(res<-clusteredTest(phy,d))
  set.seed(0)
  res<-bactdate(phy,d,nbIts=10)
  expect_silent(res<-bactdate(phy,d,nbIts=10))
  expect_is(res,'resBactDating')
  expect_is(capture_output(print(res)),'character')
  expect_is(capture_output(summary(res)),'character')
  expect_silent(plot(res,type='tree'))
  expect_silent(plot(res,type='treeCI'))
  expect_silent(plot(res,type='trace'))
  expect_silent(plot(res,type='treeRoot'))
  expect_silent(plot(res,type='scatter'))
  expect_silent(plotDualScale(res$tree))
  expect_is(as.mcmc.resBactDating(res),'mcmc')
  expect_silent(capture_output(modelcompare(res,res)))
})

