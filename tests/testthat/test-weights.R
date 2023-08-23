test_that("ATE works", {
  expect_message(
    weights <- wt_ate(c(.1, .3, .4, .3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )
  expect_equal(
    weights,
    c(1.11, 1.43, 2.50, 1.43),
    tolerance = .01
  )
})
