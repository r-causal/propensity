expect_propensity_error <- function(expr) {
  testthat::expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    expr
  )
}

# Snapshot the conditions raised by `expr`. The expression is spliced into the
# `expect_snapshot()` call rather than passed as a promise, so the recorded Code
# section reads as the call the test made. Unless `print` is set, the value is
# assigned inside the helper, which keeps it out of the recorded Output while
# still handing it back, so a caller that needs the result can assign the
# expectation itself. Evaluation happens in a child of the caller's frame, so
# nothing the helper binds reaches the test.
expect_propensity_warning <- function(expr, print = FALSE) {
  quoted <- substitute(expr)
  if (!print) {
    quoted <- call("<-", quote(out), quoted)
  }

  env <- new.env(parent = parent.frame())
  eval(bquote(testthat::expect_snapshot(cnd_class = TRUE, .(quoted))), env)

  invisible(env$out)
}

# Replace the saturated-observation count in a separation error with a
# placeholder. The count comes from IRLS on a fit with no finite maximum
# likelihood estimate, so it is a coefficient in the tens of thousands read
# against a fixed threshold, and the observations that fall the near side of it
# are few enough that a change in the last bits of that coefficient moves the
# count. It is pinned where it can be recomputed from the fixture, in the tests
# either standard error method reports it to; what a snapshot is for here is the
# sentence around it.
mask_saturated_count <- function(lines) {
  sub(
    "exactly 0 or 1 for [0-9]+ observations",
    "exactly 0 or 1 for <n> observations",
    lines
  )
}
