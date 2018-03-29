context("Test utility functions")

test_that("The function changeinorderedvec behaves as expected.", {
  v=100:1
  v=v+0.1
  for (i in 1:100) BactDating:::changeinorderedvec(v,v[sample.int(100,1)],runif(1)*100)
  expect_false(is.unsorted(rev(v)))
})

