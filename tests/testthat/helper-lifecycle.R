# lifecycle signals a given deprecation once per session and then goes quiet,
# unless it judges the call "direct". Setting `lifecycle_verbosity` to "warning"
# only forces the signal for calls it already judged direct, and for a
# deprecation raised inside a package it makes that judgment from the
# `TESTTHAT_PKG` environment variable, which testthat sets to the package under
# test when it knows it and leaves empty otherwise. Under a bare
# `testthat::test_file()` the variable is empty, so the second test to reach a
# deprecation id sees nothing at all. Forcing the option and the variable
# together makes every signal unconditional, which is what lets an assertion
# hold whatever reached the deprecation first and whichever runner drove the
# file. Assertions that a deprecation is *not* raised need the same treatment:
# without it they pass on a silence that says nothing.
with_always_deprecated <- function(code) {
  withr::local_options(lifecycle_verbosity = "warning")
  withr::local_envvar(TESTTHAT_PKG = "propensity")
  code
}
