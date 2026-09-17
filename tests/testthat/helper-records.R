# The fields of a trimming or truncation record that name units by position.
record_position_fields <- c("keep_idx", "trimmed_idx", "truncated_idx", "n_obs")

# A record whose positions were dropped keeps everything it says about the
# modification itself, its method, bounds, and refit flag, and nothing about
# the units.
without_positions <- function(meta) {
  for (field in record_position_fields) {
    meta[[field]] <- NULL
  }
  meta
}

expect_positions_dropped <- function(meta, original, info = NULL) {
  expect_false(is.null(meta), info = info)
  expect_identical(meta, without_positions(original), info = info)
}
