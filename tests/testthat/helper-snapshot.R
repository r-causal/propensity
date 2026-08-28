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

# Replace the printed magnitudes of a degenerate fit with a placeholder. A fit
# the estimating equations have no unique root at reports standard errors, and
# the test statistics built from them, that are numerical noise rather than
# quantities: their digits are whatever the platform's arithmetic left behind,
# so they differ between machines, in magnitude, in sign, and in the number of
# digits they print to, and the column widths differ with them.
#
# A magnitude below `degenerate_floor`, or above `degenerate_ceiling`, is
# masked wherever it appears among the estimate, the standard error, the test
# statistic, and the interval bounds. The first is a value that the arithmetic
# that produced it cannot tell from zero, and which side of zero it landed on is
# not a property of the fit, so the sign is masked with it. The second is a
# statistic divided by such a value, which is as much noise as its denominator.
#
# An estimate can carry either shape too, and does: a counterfactual mean the
# fit pins at zero prints a few multiples of machine precision away from it and
# is masked alongside the standard error beside it. That is the right reading of
# it rather than a mask reaching too far, since such an estimate is noise on the
# same terms. An estimate of an ordinary magnitude is untouched, and so is the
# p-value column, which the mask never reaches.
#
# The test statistic of a row whose estimate and standard error are both noise
# is masked with them however it prints: it is a ratio of two numerical zeros
# and takes whatever value their last bits imply, which need be neither small
# nor large.
#
# What survives is normalized for width. A magnitude keeps its digits but loses
# the trailing zeros that a neighboring masked value would otherwise decide the
# count of, and each run of spaces in the table collapses to one, so the column
# padding the masked tokens would have set does not reach the snapshot either.
#
# What these snapshots pin is the sentence the fit is described in and the shape
# of the report, not the noise.
degenerate_floor <- 1e-12

degenerate_ceiling <- 1e16

# The columns the mask reads, in the order the table prints them: the estimate,
# its standard error, the test statistic, and the two interval bounds.
degenerate_n_columns <- 5L

degenerate_magnitude <- "-?[0-9]+\\.[0-9]+e[+-][0-9]+"

mask_degenerate_magnitudes <- function(lines) {
  in_table <- grepl(degenerate_magnitude, lines) |
    grepl("^\\s+estimate\\s+std\\.err\\b", lines)

  lines[in_table] <- vapply(
    lines[in_table],
    mask_degenerate_row,
    character(1),
    USE.NAMES = FALSE
  )

  lines
}

mask_degenerate_row <- function(line) {
  fields <- strsplit(trimws(line), " +")[[1]]

  at <- grep(paste0("^", degenerate_magnitude, "$"), fields)
  at <- at[seq_len(min(length(at), degenerate_n_columns))]
  magnitudes <- abs(as.numeric(fields[at]))
  noise <- magnitudes <= degenerate_floor | magnitudes >= degenerate_ceiling

  if (length(at) >= 3L && all(noise[1:2])) {
    noise[[3]] <- TRUE
  }

  fields[at[noise]] <- "<degenerate>"
  fields[at[!noise]] <- sub("\\.?0+e", "e", fields[at[!noise]])

  paste(fields, collapse = " ")
}
