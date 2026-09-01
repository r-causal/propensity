# Every row of a printed result is printed exactly once, under the label the
# accessors give it.
#
# The printed label occupies a fixed-width field, so reading that field back and
# counting exact matches is what pins the count. A prefix test would accept the
# row for a longer label as the row for the shorter label it extends, which for
# a `.by` fit is the difference between the row for one subgroup and the row for
# a comparison of two. The width is set wide enough that the estimate table is
# printed in one block rather than wrapped into several.
expect_printed_labels_once <- function(res, labels) {
  withr::local_options(width = 200)
  printed <- utils::capture.output(print(res))
  fields <- trimws(
    substr(printed, 1L, max(nchar(labels))),
    which = "right"
  )
  counts <- vapply(labels, function(label) sum(fields == label), integer(1))

  testthat::expect_identical(unname(counts), rep(1L, length(labels)))
}
