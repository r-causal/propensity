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

# lifecycle decides whether a deprecation belongs to the caller or to the
# package that raised it from the environment it is handed, and says so in the
# warning: one it judges indirect carries a bullet naming the package and asking
# the reader to report an issue against it. A helper that signals a deprecation
# without passing the calling environment along is judged to have raised it
# itself, so that bullet arrives on calls the user wrote.
#
# `expr` is a quoted call, evaluated in an environment that inherits from the
# global environment, which is where a user's call lives, and `data` supplies
# the names it reads. `TESTTHAT_PKG` is unset so that lifecycle judges the
# attribution from the environment alone, which is what a user's session offers
# it. Every deprecation warning the call raises comes back, muffled so that none
# reaches the test run.
#
# A deprecation lifecycle judges indirect is signaled once per session, so a
# call that raises none has been attributed to this package just as surely as
# one that raises the bullet. Assert the count as well as the text.
deprecation_warnings_from_user <- function(expr, data = list()) {
  withr::local_options(lifecycle_verbosity = "warning")
  withr::local_envvar(TESTTHAT_PKG = NA)

  messages <- character()
  withCallingHandlers(
    eval(expr, rlang::new_environment(data, parent = globalenv())),
    lifecycle_warning_deprecated = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )

  messages
}

# The sentence lifecycle adds to a deprecation it judges indirect.
deprecation_misattributed <- function(messages) {
  any(grepl(
    "The deprecated feature was likely used in the propensity package",
    messages,
    fixed = TRUE
  ))
}
