# The tidier generics belong to the generics package, which is also where broom
# takes them from. Re-exporting them, rather than defining copies, is what lets a
# user who has only loaded propensity call `tidy()` on an `ipw` object and reach
# the same generic every other tidier method in their session dispatches from.

test_that("propensity exports the tidier generics", {
  expect_contains(
    getNamespaceExports("propensity"),
    c("tidy", "glance", "augment")
  )
})

test_that("the exported tidier generics are the ones generics owns", {
  expect_identical(propensity::tidy, generics::tidy)
  expect_identical(propensity::glance, generics::glance)
  expect_identical(propensity::augment, generics::augment)
})
