test_that("sim.one.trial() works", {
  set.seed(1)
  expect_equal(sim.one.trial(accrual = "poisson")$n,c(19,12,3))
})

test_that("get.oc.bf() works", {
  expect_equal(get.oc.bf(ntrial = 10, seed = 9, accrual = "poisson")$npatients,c(17.6,16.8,4.2))
})
