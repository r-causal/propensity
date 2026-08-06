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

# The two generics that move a result between its readings belong to
# causalgenerics, which owns the result class and the `effects` field they set.
# Re-exporting them is what lets a user who has loaded only propensity write
# `as_conditional(res)` unqualified on a result `ipw()` built.

test_that("propensity exports the presentation mode generics", {
  expect_contains(
    getNamespaceExports("propensity"),
    c("as_marginal", "as_conditional")
  )
})

test_that("the exported mode generics are the ones causalgenerics owns", {
  expect_identical(propensity::as_marginal, causalgenerics::as_marginal)
  expect_identical(propensity::as_conditional, causalgenerics::as_conditional)
})

# `ipw()` is the generic an estimating package supplies a method for, and
# `estimand()`, `estimand<-`, and `is_causal_wt()` read and set what the weights
# record about themselves. All four belong to causalgenerics, which owns the
# result class and the weight contract they describe. Re-exporting them is what
# lets a user who has loaded only propensity call `ipw()` on their models and
# reach the same generic another package building on those classes dispatches
# from.

test_that("propensity exports the shared causal generics", {
  expect_contains(
    getNamespaceExports("propensity"),
    c("ipw", "estimand", "estimand<-", "is_causal_wt")
  )
})

test_that("the exported causal generics are the ones causalgenerics owns", {
  expect_identical(propensity::ipw, causalgenerics::ipw)
  expect_identical(propensity::estimand, causalgenerics::estimand)
  expect_identical(propensity::`estimand<-`, causalgenerics::`estimand<-`)
  expect_identical(propensity::is_causal_wt, causalgenerics::is_causal_wt)
})
