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
# masked. The first is a value that the arithmetic that produced it cannot tell
# from zero, and which side of zero it landed on is not a property of the fit,
# so the sign is masked with it. The second is a statistic divided by such a
# value, which is as much noise as its denominator.
#
# An estimate can carry either shape too, and does: a counterfactual mean the
# fit pins at zero prints a few multiples of machine precision away from it and
# is masked alongside the standard error beside it. That is the right reading of
# it rather than a mask reaching too far, since such an estimate is noise on the
# same terms. An estimate of an ordinary magnitude is untouched.
#
# The test statistic of a row whose standard error is noise is masked with it
# however the statistic itself prints: it is a quantity divided by a numerical
# zero and takes whatever value the denominator's last bits imply, which need be
# neither small nor large.
#
# The mask reads the value columns up to the confidence level, which the table
# prints second from last, so it never reaches the p-value beside it. Reading
# them from the right rather than counting them from the left is what keeps that
# true of a row whose estimate prints in fixed rather than scientific notation.
#
# What survives is normalized for width. A magnitude keeps its digits but loses
# the trailing zeros that a neighboring masked value would otherwise decide the
# count of, and each run of spaces in the table collapses to one, so the column
# padding the masked tokens would have set does not reach the snapshot either.
# Only the table is normalized, from its column header to the rule that closes
# it, or to the first blank line if no rule was printed, so a line printed above
# or below it keeps the indentation it was printed with.
#
# What these snapshots pin is the sentence the fit is described in and the shape
# of the report, not the noise.
degenerate_floor <- 1e-12

degenerate_ceiling <- 1e16

# A printed magnitude, in fixed or scientific notation. The decimal point is
# required, so the integers a row label is written with are not read as values.
degenerate_magnitude <- "-?[0-9]+\\.[0-9]+(e[+-][0-9]+)?"

# The columns the mask leaves alone at the right of every row: the confidence
# level and the p-value.
degenerate_n_trailing <- 2L

mask_degenerate_magnitudes <- function(lines) {
  in_table <- degenerate_table_lines(lines)

  lines[in_table] <- vapply(
    lines[in_table],
    mask_degenerate_row,
    character(1),
    USE.NAMES = FALSE
  )

  lines
}

# The estimates table runs from the line naming its columns to the rule that
# closes it, and nothing outside those bounds is a row of it. A table no row of
# which earns a significance star is printed without that rule, so the blank
# line after the last row closes the table too.
degenerate_table_lines <- function(lines) {
  opens <- grepl("^\\s+estimate\\s+std\\.err\\b", lines)
  closes <- grepl("^\\s*(---)?\\s*$", lines)

  in_table <- logical(length(lines))
  open <- FALSE

  for (i in seq_along(lines)) {
    if (opens[[i]]) {
      open <- TRUE
    } else if (closes[[i]]) {
      open <- FALSE
    }

    in_table[[i]] <- open
  }

  in_table
}

mask_degenerate_row <- function(line) {
  fields <- strsplit(trimws(line), " +")[[1]]

  at <- grep(paste0("^", degenerate_magnitude, "$"), fields)
  at <- utils::head(at, -degenerate_n_trailing)

  magnitudes <- abs(as.numeric(fields[at]))
  noise <- magnitudes <= degenerate_floor | magnitudes >= degenerate_ceiling

  if (length(at) >= 3L && noise[[2]]) {
    noise[[3]] <- TRUE
  }

  fields[at[noise]] <- "<degenerate>"
  fields[at[!noise]] <- drop_trailing_zeros(fields[at[!noise]])

  paste(fields, collapse = " ")
}

# Trailing zeros a column's own format decided the count of, in either notation.
drop_trailing_zeros <- function(magnitudes) {
  scientific <- grepl("e", magnitudes, fixed = TRUE)
  magnitudes[scientific] <- sub("\\.?0+e", "e", magnitudes[scientific])
  magnitudes[!scientific] <- sub("\\.?0+$", "", magnitudes[!scientific])

  magnitudes
}
