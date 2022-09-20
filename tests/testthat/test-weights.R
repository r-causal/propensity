test_that("ATE works", {
  expect_equal(
    wt_ate(c(.1, .3, .4, .3), .exposure = c(0, 0, 1, 0)),
    c(1.11, 1.43, 2.50, 1.43),
    tolerance = .01
  )
})
