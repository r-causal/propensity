# Summary test to understand which base R methods work for each class
library(testthat)

test_that("Summary: Base R methods status for psw class", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "ate")

  # Methods that work (through vctrs defaults)
  expect_s3_class(x[1], "psw") # Subsetting works
  expect_s3_class(sort(x), "psw") # Sort works
  expect_s3_class(unique(x), "psw") # Unique works
  expect_s3_class(rev(x), "psw") # Rev works
  expect_s3_class(head(x, 2), "psw") # Head works
  expect_s3_class(tail(x, 2), "psw") # Tail works
  expect_s3_class(rep(x, 2), "psw") # Rep works

  # Methods that return plain types (as expected)
  expect_type(duplicated(x), "logical")
  expect_type(is.na(x), "logical")
  expect_type(which(as.numeric(x) > 0.15), "integer")
  expect_type(order(x), "integer")
  expect_type(rank(x), "double")

  # na.omit works and preserves class
  x_na <- psw(c(0.1, NA, 0.3), estimand = "ate")
  expect_s3_class(na.omit(x_na), "psw")
})

test_that("Summary: Base R methods status for ps_trim class", {
  # The draw decides whether anything falls outside the trimming bounds, so it
  # is seeded: an unseeded sample that lands entirely inside them leaves nothing
  # trimmed and no missing value for `anyNA()` to find.
  withr::local_seed(1)
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Methods that preserve class (implemented)
  expect_s3_class(x[1:5], "ps_trim") # [ method implemented

  # Methods that work through vctrs but may have issues
  sorted <- sort(x) # Drops NAs by default
  expect_true(length(sorted) <= length(x))

  # These preserve class through vctrs
  expect_s3_class(unique(x), "ps_trim")
  expect_s3_class(rev(x), "ps_trim")
  expect_s3_class(head(x, 5), "ps_trim")
  expect_s3_class(tail(x, 5), "ps_trim")
  expect_s3_class(rep(x, 2), "ps_trim")

  # Methods that return plain types
  expect_type(duplicated(x), "logical")
  expect_type(is.na(x), "logical")
  expect_true(anyNA(x)) # ps_trim has NAs for trimmed values

  # na.omit works
  expect_s3_class(na.omit(x), "ps_trim")
})

test_that("Summary: Base R methods status for ps_trunc class", {
  # Seeded for the reason the trimming test above is: the draw decides how many
  # values the bounds have anything to do with.
  withr::local_seed(1)
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Methods that preserve class (implemented)
  expect_s3_class(x[1:5], "ps_trunc") # [ method implemented

  # These preserve class through vctrs
  expect_s3_class(sort(x), "ps_trunc")
  expect_s3_class(unique(x), "ps_trunc")
  expect_s3_class(rev(x), "ps_trunc")
  expect_s3_class(head(x, 5), "ps_trunc")
  expect_s3_class(tail(x, 5), "ps_trunc")
  expect_s3_class(rep(x, 2), "ps_trunc")

  # Methods that return plain types
  expect_type(duplicated(x), "logical")
  expect_type(is.na(x), "logical")
  expect_false(anyNA(x)) # ps_trunc doesn't have NAs

  # na.omit works
  expect_s3_class(na.omit(x), "ps_trunc")
})

test_that("Summary: Arithmetic behavior differences", {
  psw_obj <- psw(c(0.1, 0.2), estimand = "ate")
  ps_trim_obj <- ps_trim(
    c(0.1, 0.2, 0.9),
    method = "ps",
    lower = 0.15,
    upper = 0.8
  )
  ps_trunc_obj <- ps_trunc(
    c(0.1, 0.2, 0.9),
    method = "ps",
    lower = 0.15,
    upper = 0.8
  )

  # psw preserves class for unary operations
  expect_s3_class(-psw_obj, "psw")
  expect_s3_class(+psw_obj, "psw")

  # ps_trim and ps_trunc return numeric for unary operations
  expect_type(-ps_trim_obj, "double")
  expect_type(+ps_trim_obj, "double")
  expect_type(-ps_trunc_obj, "double")
  expect_type(+ps_trunc_obj, "double")

  # Binary operations with numeric
  expect_s3_class(psw_obj * 2, "psw") # psw preserves class
  expect_type(ps_trim_obj * 2, "double") # ps_trim returns numeric
  expect_type(ps_trunc_obj * 2, "double") # ps_trunc returns numeric
})

test_that("Summary: Key differences and missing methods", {
  defines_method <- function(generic, class) {
    !is.null(utils::getS3method(generic, class, optional = TRUE))
  }
  defined_for <- function(generics, class) {
    generics[vapply(generics, defines_method, logical(1), class = class)]
  }

  # psw takes `[`, summary(), min(), max(), range(), median(), quantile(),
  # anyDuplicated(), and diff() from causalgenerics' causal_wts methods;
  # ps_trim and ps_trunc each define their own.
  inherited <- c(
    "[",
    "summary",
    "min",
    "max",
    "range",
    "median",
    "quantile",
    "anyDuplicated",
    "diff"
  )
  expect_equal(defined_for(inherited, "psw"), character(0))
  expect_equal(defined_for(inherited, "causal_wts"), inherited)
  expect_equal(defined_for(inherited, "ps_trim"), inherited)
  expect_equal(defined_for(inherited, "ps_trunc"), inherited)

  # ps_trim and ps_trunc implement sort(), unique(), and rep() so that the
  # record of which units were modified follows the result. psw implements
  # none of the three and takes the vctrs defaults.
  record_keeping <- c("sort", "unique", "rep")
  expect_equal(defined_for(record_keeping, "psw"), character(0))
  expect_equal(defined_for(record_keeping, "causal_wts"), character(0))
  expect_equal(defined_for(record_keeping, "ps_trim"), record_keeping)
  expect_equal(defined_for(record_keeping, "ps_trunc"), record_keeping)

  # na.omit() has a method for ps_trim alone. psw and ps_trunc reach the
  # default, which keeps their class, as the tests above assert.
  expect_true(defines_method("na.omit", "ps_trim"))
  expect_false(defines_method("na.omit", "ps_trunc"))
  expect_false(defines_method("na.omit", "psw"))

  # Each class defines vec_restore() exactly once.
  expect_equal(
    sort(grep("^vec_restore\\.", ls(asNamespace("propensity")), value = TRUE)),
    c("vec_restore.ps_trim", "vec_restore.ps_trunc", "vec_restore.psw")
  )

  # All three classes return plain numeric from min(), max(), range(),
  # median(), and quantile(), an integer from anyDuplicated(), and plain
  # numeric from diff(). ps_trim marks trimmed units NA, so its summaries
  # need na.rm.
  psw_obj <- psw(c(0.1, 0.2, 0.3), estimand = "ate")
  ps_trim_obj <- ps_trim(
    c(0.1, 0.2, 0.9),
    method = "ps",
    lower = 0.15,
    upper = 0.8
  )
  ps_trunc_obj <- ps_trunc(
    c(0.1, 0.2, 0.9),
    method = "ps",
    lower = 0.15,
    upper = 0.8
  )

  for (x in list(psw_obj, ps_trim_obj, ps_trunc_obj)) {
    expect_identical(class(min(x, na.rm = TRUE)), "numeric")
    expect_identical(class(max(x, na.rm = TRUE)), "numeric")
    expect_identical(class(range(x, na.rm = TRUE)), "numeric")
    expect_length(range(x, na.rm = TRUE), 2)
    expect_identical(class(median(x, na.rm = TRUE)), "numeric")
    expect_identical(class(quantile(x, na.rm = TRUE)), "numeric")
    expect_identical(class(anyDuplicated(x)), "integer")
    expect_identical(class(diff(x)), "numeric")
  }
})
