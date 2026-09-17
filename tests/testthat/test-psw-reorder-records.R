# Reordering weights that carry position records -----------------------------

# A trimming, truncation, or weight truncation record names units by position.
# A subscript moves units to new positions, so `[` and everything built on it
# (`rev()`, `sort()`, `x[order(x)]`) carry each record through the subscript.
# A slice that reaches the restore without its subscript, such as `vec_slice()`
# or `dplyr::arrange()`, cannot place the record and drops it silently.

# Units 1 and 4 fall outside (0.1, 0.9); their weights are missing.
reorder_trimmed_psw <- function() {
  w <- psw(
    c(NA, 2, 1.6, NA, 3.3),
    estimand = "ate; trimmed",
    trimmed = TRUE
  )
  attr(w, "ps_trim_meta") <- list(
    method = "ps",
    lower = 0.1,
    upper = 0.9,
    focal_inverted = FALSE,
    keep_idx = c(2L, 3L, 5L),
    trimmed_idx = c(1L, 4L),
    n_obs = 5L
  )
  w
}

# Units 1 and 4 were pinned to the bounds.
reorder_truncated_psw <- function() {
  w <- psw(
    c(10, 2, 1.6, 10, 3.3),
    estimand = "ate; truncated",
    truncated = TRUE
  )
  attr(w, "ps_trunc_meta") <- list(
    method = "ps",
    lower_bound = 0.1,
    upper_bound = 0.9,
    truncated_idx = c(1L, 4L),
    n_obs = 5L
  )
  w
}

# Units 3 and 10 were brought down to the second largest weight, which unit 2
# already held, so the bound value appears at a unit the truncation never moved.
reorder_wt_truncated_psw <- function() {
  wt_trunc(psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100)), "count", upper = 2)
}

reorder_fixtures <- function() {
  list(
    trimmed = list(
      weights = reorder_trimmed_psw(),
      attr = "ps_trim_meta",
      query = is_unit_trimmed
    ),
    truncated = list(
      weights = reorder_truncated_psw(),
      attr = "ps_trunc_meta",
      query = is_unit_truncated
    ),
    wt_truncated = list(
      weights = reorder_wt_truncated_psw(),
      attr = "psw_trunc_meta",
      query = is_unit_wt_truncated
    )
  )
}

# The units a record names after the subscript `i`, read off the units it named
# before. A missing position names no unit.
expected_units <- function(units, i) {
  out <- units[i]
  out[is.na(out)] <- FALSE
  out
}

expect_reindexed <- function(out, fixture, i, label) {
  units <- fixture$query(fixture$weights)
  expect_s3_class(out, "psw")
  expect_identical(
    fixture$query(out),
    expected_units(units, i),
    info = label
  )
  expect_identical(
    attr(out, fixture$attr)$n_obs,
    length(i),
    info = label
  )
}

expect_record_dropped <- function(out, fixture, label) {
  expect_s3_class(out, "psw")
  expect_null(attr(out, fixture$attr), info = label)
  expect_error(
    fixture$query(out),
    class = "propensity_missing_meta_error",
    info = label
  )
}

test_that("the fixtures name the units they were built to name", {
  fixtures <- reorder_fixtures()

  expect_identical(
    which(fixtures$trimmed$query(fixtures$trimmed$weights)),
    c(1L, 4L)
  )
  expect_identical(
    which(fixtures$truncated$query(fixtures$truncated$weights)),
    c(1L, 4L)
  )
  expect_identical(
    which(fixtures$wt_truncated$query(fixtures$wt_truncated$weights)),
    c(3L, 10L)
  )
})

test_that("rev() carries each position record onto the reversed weights", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(rev(w))
    expect_reindexed(out, fixture, rev(seq_along(w)), label)
    expect_identical(vec_data(out), rev(vec_data(w)), info = label)
  }
})

test_that("rev() on the reported examples names the units that moved", {
  t <- reorder_wt_truncated_psw()
  expect_identical(which(is_unit_wt_truncated(rev(t))), c(1L, 8L))

  tr <- reorder_trimmed_psw()
  expect_identical(which(is_unit_trimmed(rev(tr))), c(2L, 5L))
  expect_identical(which(is.na(rev(tr))), c(2L, 5L))
})

test_that("sort() carries each position record through the order it sorts by", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    data <- vec_data(w)

    # `sort()` drops missing weights by default, so the record follows a
    # subscript shorter than the weights.
    out <- expect_silent(sort(w))
    expect_reindexed(out, fixture, order(data, na.last = NA), label)

    out <- expect_silent(sort(w, decreasing = TRUE, na.last = TRUE))
    expect_reindexed(
      out,
      fixture,
      order(data, decreasing = TRUE, na.last = TRUE),
      label
    )
  }
})

test_that("subsetting by an order carries each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    ord <- order(-vec_data(w))

    out <- expect_silent(w[ord])
    expect_reindexed(out, fixture, ord, label)
  }
})

test_that("a positive permutation carries each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    i <- c(2L, seq_along(w)[-2])

    out <- expect_silent(w[i])
    expect_reindexed(out, fixture, i, label)
    expect_identical(vec_data(out), vec_data(w)[i], info = label)
  }
})

test_that("a negative subscript carries each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(w[-1])
    expect_reindexed(out, fixture, seq_along(w)[-1], label)

    # A zero in a negative subscript removes nothing.
    out <- expect_silent(w[-0])
    expect_length(out, 0)
  }
})

test_that("a logical subscript carries each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    keep <- rep_len(c(TRUE, FALSE), length(w))
    keep[4] <- TRUE

    out <- expect_silent(w[keep])
    expect_reindexed(out, fixture, which(keep), label)

    out <- expect_silent(w[TRUE])
    expect_identical(attr(out, fixture$attr), attr(w, fixture$attr))
  }
})

test_that("a character subscript carries each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    names(w) <- letters[seq_along(w)]
    i <- rev(names(w))

    out <- expect_silent(w[i])
    expect_reindexed(out, fixture, rev(seq_along(w)), label)
    expect_identical(names(out), i, info = label)
  }
})

test_that("a subscript naming a unit twice or no unit reindexes element by element", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    # Unit 2 is untouched in every fixture.
    flagged <- which(fixture$query(w))
    i <- c(flagged[2], flagged[2], NA, 2L, flagged[1])

    out <- expect_silent(w[i])
    expect_reindexed(out, fixture, i, label)
    expect_identical(fixture$query(out), c(TRUE, TRUE, FALSE, FALSE, TRUE))
  }
})

test_that("an identity subscript keeps each position record unchanged", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(w[seq_along(w)])
    expect_identical(attr(out, fixture$attr), attr(w, fixture$attr))

    out <- expect_silent(w[])
    expect_identical(attr(out, fixture$attr), attr(w, fixture$attr))
  }
})

test_that("an empty subscript leaves the prototype as the restore builds it", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights

    expect_identical(
      attributes(w[integer()]),
      attributes(vctrs::vec_slice(w, integer())),
      info = label
    )
    expect_identical(
      attributes(w[0]),
      attributes(vctrs::vec_ptype(w)),
      info = label
    )
  }
})

test_that("a record that does not cover the weights is dropped by a subscript", {
  w <- reorder_wt_truncated_psw()
  w[12] <- 2
  expect_length(w, 12)

  out <- expect_silent(w[rev(seq_along(w))])
  expect_null(attr(out, "psw_trunc_meta"))
  expect_true(is_wt_truncated(out))
  expect_error(
    is_unit_wt_truncated(out),
    class = "propensity_missing_meta_error"
  )
})

test_that("vec_slice() at the same length drops each position record silently", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(vctrs::vec_slice(w, rev(seq_along(w))))
    expect_record_dropped(out, fixture, label)
    expect_identical(estimand(out), estimand(w), info = label)
    expect_identical(vec_data(out), rev(vec_data(w)), info = label)
  }

  w <- reorder_wt_truncated_psw()
  out <- vctrs::vec_slice(w, rev(seq_along(w)))
  expect_true(is_wt_truncated(out))
})

test_that("dplyr::arrange() on a psw column drops each position record silently", {
  skip_if_not_installed("dplyr")

  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    df <- data.frame(o = rev(seq_along(w)))
    df$w <- w

    out <- expect_silent(dplyr::arrange(df, o))
    expect_record_dropped(out$w, fixture, label)
    expect_identical(vec_data(out$w), rev(vec_data(w)), info = label)
  }
})

test_that("reordering a data frame with base subsetting carries each record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    df <- data.frame(o = rev(seq_along(w)))
    df$w <- w

    out <- expect_silent(df[order(df$o), , drop = FALSE])
    expect_reindexed(out$w, fixture, rev(seq_along(w)), label)
  }
})

test_that("same-length arithmetic keeps each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    meta <- attr(w, fixture$attr)

    results <- list(
      `w * 2` = expect_silent(w * 2),
      `2 * w` = expect_silent(2 * w),
      `w * 1L` = expect_silent(w * 1L),
      `-w` = expect_silent(-w),
      `w * w` = expect_silent(w * w),
      `w / sum(w, na.rm = TRUE)` = expect_silent(w / sum(w, na.rm = TRUE))
    )

    for (op in names(results)) {
      out <- results[[op]]
      expect_s3_class(out, "psw")
      expect_identical(attr(out, fixture$attr), meta, info = c(label, op))
    }
  }
})

test_that("same-length subassignment keeps each position record", {
  for (label in names(reorder_fixtures())) {
    fixture <- reorder_fixtures()[[label]]
    w <- fixture$weights
    meta <- attr(w, fixture$attr)
    units <- fixture$query(w)

    expect_silent({
      w[2] <- 5
    })
    expect_identical(attr(w, fixture$attr), meta, info = label)
    expect_identical(fixture$query(w), units, info = label)

    expect_silent({
      w[] <- rev(vec_data(w))
    })
    expect_identical(attr(w, fixture$attr), meta, info = label)
  }
})

test_that("truncating reordered weights records the units at their new positions", {
  w <- psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100))
  out <- wt_trunc(rev(w), "count", upper = 2)

  expect_identical(which(is_unit_wt_truncated(out)), c(1L, 8L))
  expect_identical(
    which(is_unit_wt_truncated(rev(out))),
    which(is_unit_wt_truncated(wt_trunc(w, "count", upper = 2)))
  )
})

test_that("a combine of reordered weights still compares their bounds", {
  a <- reorder_wt_truncated_psw()
  b <- wt_trunc(psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100)), "count", upper = 1)

  out <- expect_silent(vctrs::vec_c(rev(a), a))
  expect_true(is_wt_truncated(out))
  expect_null(attr(out, "psw_trunc_meta"))
  expect_null(attr(out, "psw_trunc_bound"))
  expect_null(attr(out, "psw_conflicted_attrs"))

  expect_warning(
    out <- vctrs::vec_c(rev(a), a, rev(b)),
    class = "propensity_coercion_warning"
  )
  expect_false(inherits(out, "psw"))
})
