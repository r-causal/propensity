# The record a weight truncation leaves on a psw ------------------------------

# Truncating the weights themselves is recorded apart from truncating the
# propensity scores they were built from. The `wt_truncated` flag says the
# weights were bounded, and the `psw_trunc_meta` record says at what bound and
# which units it moved. The record is positional, so it follows the rules the
# trimming record follows: kept through elementwise arithmetic and
# subassignment, carried through the subscript of `[`, dropped without comment
# by any other slice or length change, and left off the prototype a combine
# builds. The flag describes the vector as a whole and survives all of those.

# Units 2 and 5 sit at the bound, so the positional query has something to
# report and a result that lost the record reports something visibly different.
wt_trunc_record <- function(
  method = "wt",
  upper_value = 3,
  truncated_idx = c(2L, 5L),
  n_obs = 5L
) {
  new_psw_trunc_meta(
    method = method,
    lower = NULL,
    upper = upper_value,
    lower_value = NULL,
    upper_value = upper_value,
    truncated_idx = truncated_idx,
    n_obs = n_obs
  )
}

wt_truncated_psw <- function(
  x = c(1.2, 3, 1.5, 2.1, 3),
  record = wt_trunc_record()
) {
  w <- psw(x, estimand = "ate", wt_truncated = TRUE)
  attr(w, "psw_trunc_meta") <- record
  w
}

# The same five observations bounded lower, so the two records disagree on the
# bound as well as on the units it moved.
alt_wt_truncated_psw <- function() {
  wt_truncated_psw(
    x = c(1.2, 2.5, 1.5, 2.1, 2.5),
    record = wt_trunc_record(upper_value = 2.5, truncated_idx = c(2L, 4L, 5L))
  )
}

wt_truncated_units <- c(FALSE, TRUE, FALSE, FALSE, TRUE)

test_that("the truncation record holds the bound and the units it moved", {
  meta <- new_psw_trunc_meta(
    method = "pctl",
    lower = 0.01,
    upper = 0.99,
    lower_value = 0.2,
    upper_value = 7.5,
    truncated_idx = c(1L, 4L),
    n_obs = 10L
  )

  expect_s3_class(meta, "propensity_psw_trunc_meta")
  expect_named(
    meta,
    c(
      "method",
      "lower",
      "upper",
      "lower_value",
      "upper_value",
      "truncated_idx",
      "n_obs"
    )
  )
  expect_identical(meta$method, "pctl")
  expect_identical(meta$lower, 0.01)
  expect_identical(meta$upper, 0.99)
  expect_identical(meta$lower_value, 0.2)
  expect_identical(meta$upper_value, 7.5)
  expect_identical(meta$truncated_idx, c(1L, 4L))
  expect_identical(meta$n_obs, 10L)
})

test_that("psw() records whether the weights were truncated", {
  expect_true(is_wt_truncated(psw(c(1, 2), wt_truncated = TRUE)))
  expect_false(is_wt_truncated(psw(c(1, 2))))
  expect_false(is_wt_truncated(new_psw(c(1, 2))))
  expect_true(is_wt_truncated(new_psw(c(1, 2), wt_truncated = TRUE)))
})

test_that("weight truncation is recorded apart from score truncation", {
  # A score truncation and a weight truncation are different operations, so a
  # psw can carry either one or both without one answering for the other.
  weights_only <- psw(c(1, 2), wt_truncated = TRUE)
  expect_true(is_wt_truncated(weights_only))
  expect_false(is_ps_truncated(weights_only))

  scores_only <- psw(c(1, 2), truncated = TRUE)
  expect_true(is_ps_truncated(scores_only))
  expect_false(is_wt_truncated(scores_only))

  both <- psw(c(1, 2), truncated = TRUE, wt_truncated = TRUE)
  expect_true(is_ps_truncated(both))
  expect_true(is_wt_truncated(both))
})

test_that("is_unit_wt_truncated() reads the units the record names", {
  w <- wt_truncated_psw()

  expect_true(is_wt_truncated(w))
  expect_identical(attr(w, "psw_trunc_meta"), wt_trunc_record())
  expect_identical(expect_silent(is_unit_wt_truncated(w)), wt_truncated_units)

  # The score-scale query knows nothing of a weight truncation.
  expect_identical(is_unit_truncated(w), rep(FALSE, 5))
})

test_that("is_unit_wt_truncated() answers once per unit of untruncated weights", {
  expect_identical(
    expect_silent(is_unit_wt_truncated(psw(c(1, 2, 3)))),
    c(FALSE, FALSE, FALSE)
  )

  # No observations, no answers.
  expect_identical(is_unit_wt_truncated(psw(double())), logical(0))
  expect_identical(
    is_unit_wt_truncated(wt_truncated_psw()[0]),
    logical(0)
  )
})

test_that("length-preserving psw arithmetic keeps the weight truncation record", {
  w <- wt_truncated_psw()
  meta <- attr(w, "psw_trunc_meta")

  results <- list(
    `w * 1` = expect_silent(w * 1),
    `w / sum(w)` = expect_silent(w / sum(w)),
    `-w` = expect_silent(-w),
    `2 * w` = expect_silent(2 * w),
    `w * 1L` = expect_silent(w * 1L)
  )

  for (label in names(results)) {
    out <- results[[label]]
    expect_s3_class(out, "psw")
    expect_identical(attr(out, "psw_trunc_meta"), meta, info = label)
    expect_true(is_wt_truncated(out), info = label)
    expect_identical(
      is_unit_wt_truncated(out),
      wt_truncated_units,
      info = label
    )
  }
})

test_that("a full-length psw subset keeps the weight truncation record", {
  w <- wt_truncated_psw()
  meta <- attr(w, "psw_trunc_meta")

  whole <- expect_silent(w[seq_along(w)])
  expect_identical(attr(whole, "psw_trunc_meta"), meta)
  expect_identical(is_unit_wt_truncated(whole), wt_truncated_units)

  # A slice is not handed to anything that knows its subscript, so even one
  # that leaves every unit in place cannot vouch for the record's positions.
  sliced <- expect_silent(vec_slice(w, seq_along(w)))
  expect_positions_dropped(attr(sliced, "psw_trunc_meta"), meta)
  expect_true(is_wt_truncated(sliced))
})

test_that("shortening a psw re-indexes the weight truncation record through `[`", {
  w <- wt_truncated_psw()

  sub <- expect_silent(w[1:2])
  expect_s3_class(sub, "psw")
  expect_length(sub, 2)
  expect_identical(
    attr(sub, "psw_trunc_meta"),
    wt_trunc_record(truncated_idx = 2L, n_obs = 2L)
  )
  expect_true(is_wt_truncated(sub))
  expect_identical(estimand(sub), "ate")
  expect_identical(is_unit_wt_truncated(sub), c(FALSE, TRUE))
})

test_that("shortening a psw by a slice drops the record's positions and keeps the flag", {
  w <- wt_truncated_psw()

  sliced <- expect_silent(vec_slice(w, 1:2))
  expect_positions_dropped(
    attr(sliced, "psw_trunc_meta"),
    attr(w, "psw_trunc_meta")
  )
  expect_true(is_wt_truncated(sliced))
  expect_identical(estimand(sliced), "ate")

  # The positions are gone, so the positional query has nothing to answer from
  # and refuses rather than reporting every unit as untouched.
  expect_error(
    is_unit_wt_truncated(sliced),
    class = "propensity_missing_meta_error"
  )
})

test_that("is_unit_wt_truncated() refuses flagged weights with no record", {
  w <- psw(c(1, 2, 3), estimand = "ate", wt_truncated = TRUE)

  cnd <- expect_error(
    is_unit_wt_truncated(w),
    class = "propensity_missing_meta_error"
  )
  expect_s3_class(cnd, "propensity_error")
  expect_propensity_error(is_unit_wt_truncated(w))
})

test_that("growing a weight-truncated psw leaves a record that no longer covers it", {
  # `[<-` casts the replacement and then leaves base R to preserve the target's
  # attributes, so the record crosses a length change no restore ever sees. It
  # describes five observations and the weights now hold seven.
  w <- wt_truncated_psw()
  expect_silent({
    w[7] <- 2
  })

  expect_length(w, 7)
  expect_identical(attr(w, "psw_trunc_meta"), wt_trunc_record())
  expect_true(is_wt_truncated(w))
  expect_error(
    is_unit_wt_truncated(w),
    class = "propensity_missing_meta_error"
  )
  expect_propensity_error(is_unit_wt_truncated(w))
})

test_that("weights read out of an outcome model's frame refuse the weight query", {
  # `model.frame()` drops the `NA`-weighted rows in C and re-attaches the
  # original attributes to the shortened column, so the weights arrive carrying
  # a record written for rows that are gone.
  w <- wt_truncated_psw(x = c(1.2, 3, NA, 2.1, 3))
  dat <- data.frame(y = c(1, 2, 3, 4, 5), x = c(1, 3, 2, 5, 4))

  fit <- lm(y ~ x, data = dat, weights = w)
  model_wts <- model.frame(fit)[["(weights)"]]

  expect_s3_class(model_wts, "psw")
  expect_length(model_wts, 4)
  expect_identical(attr(model_wts, "psw_trunc_meta"), wt_trunc_record())
  expect_true(is_wt_truncated(model_wts))
  expect_error(
    is_unit_wt_truncated(model_wts),
    class = "propensity_missing_meta_error"
  )
})

test_that("a psw product carries a weight truncation record only one operand has", {
  w <- wt_truncated_psw()
  meta <- attr(w, "psw_trunc_meta")
  cens <- wt_cens(
    c(0.8, 0.7, 0.9, 0.6, 0.75),
    c(1, 1, 0, 1, 1),
    exposure_type = "binary"
  )

  out <- expect_silent(w * cens)
  expect_s3_class(out, "psw")
  expect_identical(attr(out, "psw_trunc_meta"), meta)
  expect_true(is_wt_truncated(out))
  expect_identical(is_unit_wt_truncated(out), wt_truncated_units)

  # The record and the flag travel from whichever operand holds them.
  reversed <- expect_silent(cens * w)
  expect_identical(attr(reversed, "psw_trunc_meta"), meta)
  expect_true(is_wt_truncated(reversed))
  expect_identical(is_unit_wt_truncated(reversed), wt_truncated_units)

  ones <- expect_silent(w * psw(rep(1, 5), estimand = "ate"))
  expect_identical(attr(ones, "psw_trunc_meta"), meta)
  expect_true(is_wt_truncated(ones))
})

test_that("a psw product carries a weight truncation record both operands share", {
  w <- wt_truncated_psw()

  out <- expect_silent(w * wt_truncated_psw())
  expect_identical(attr(out, "psw_trunc_meta"), wt_trunc_record())
  expect_true(is_wt_truncated(out))
  expect_identical(is_unit_wt_truncated(out), wt_truncated_units)
})

test_that("a psw product drops conflicting weight truncation records with one warning", {
  w <- wt_truncated_psw()
  alt <- alt_wt_truncated_psw()

  classes <- character()
  out <- withCallingHandlers(
    w * alt,
    warning = function(cnd) {
      classes <<- c(classes, class(cnd)[[1]])
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(classes, "propensity_metadata_conflict_warning")

  expect_s3_class(out, "psw")
  expect_null(attr(out, "psw_trunc_meta"))
  expect_true(is_wt_truncated(out))
  expect_error(
    is_unit_wt_truncated(out),
    class = "propensity_missing_meta_error"
  )

  cnd <- expect_warning(
    w * alt,
    class = "propensity_metadata_conflict_warning"
  )
  expect_true(grepl("psw_trunc_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("combining a weight-truncated and an untruncated psw downgrades to numeric", {
  w <- wt_truncated_psw()
  untruncated <- psw(c(1, 2), estimand = "ate")

  cnd <- expect_warning(
    out <- c(w, untruncated),
    class = "propensity_coercion_warning"
  )
  expect_true(grepl(
    "different weight truncation status",
    conditionMessage(cnd),
    fixed = TRUE
  ))
  expect_false(is_psw(out))
  expect_type(out, "double")
  expect_equal(out, c(1.2, 3, 1.5, 2.1, 3, 1, 2))

  expect_warning(
    reversed <- c(untruncated, w),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(reversed))
})

test_that("weights that differ in weight truncation status do not cast to each other", {
  w <- wt_truncated_psw()
  untruncated <- psw(c(1, 2), estimand = "ate")

  cnd <- expect_error(
    vec_cast(untruncated, w),
    class = "vctrs_error_incompatible_type"
  )
  expect_true(grepl(
    "different weight truncation status",
    conditionMessage(cnd),
    fixed = TRUE
  ))
  expect_error(
    vec_cast(w, untruncated),
    class = "vctrs_error_incompatible_type"
  )
})

test_that("combining weight-truncated psw objects keeps the flag and the bound", {
  w <- wt_truncated_psw()

  out <- expect_silent(c(w, w))
  expect_s3_class(out, "psw")
  expect_length(out, 10)
  expect_true(is_wt_truncated(out))
  expect_positions_dropped(attr(out, "psw_trunc_meta"), wt_trunc_record())

  # The record's positions describe one input, so the combined weights have no
  # positional answer to give.
  expect_error(
    is_unit_wt_truncated(out),
    class = "propensity_missing_meta_error"
  )

  proto <- expect_silent(vec_ptype2(w, w))
  expect_s3_class(proto, "psw")
  expect_true(is_wt_truncated(proto))
  expect_positions_dropped(attr(proto, "psw_trunc_meta"), wt_trunc_record())
})

test_that("subassigning into weight-truncated weights keeps the record", {
  # Subassignment casts the replacement to the target's own type, which has
  # already dropped nothing it compares, so it must not read as a disagreement.
  w <- wt_truncated_psw()
  expect_silent({
    w[1] <- w[3]
  })

  expect_s3_class(w, "psw")
  expect_identical(attr(w, "psw_trunc_meta"), wt_trunc_record())
  expect_identical(is_unit_wt_truncated(w), wt_truncated_units)
})

# The reason a bound mismatch gives, which tells it apart from weights that
# differ in whether they were truncated at all. Returns the combined value.
expect_bound_mismatch <- function(expr) {
  value <- NULL
  cnd <- expect_warning(
    value <- expr,
    class = "propensity_coercion_warning"
  )
  expect_true(grepl(
    "different weight truncation bounds",
    conditionMessage(cnd),
    fixed = TRUE
  ))
  invisible(value)
}

test_that("combining weights truncated at different bounds downgrades to numeric", {
  w <- wt_truncated_psw()
  alt <- alt_wt_truncated_psw()

  out <- expect_bound_mismatch(c(w, alt))
  expect_false(is_psw(out))
  expect_type(out, "double")

  reversed <- expect_bound_mismatch(c(alt, w))
  expect_false(is_psw(reversed))

  proto <- expect_bound_mismatch(vec_ptype2(w, alt))
  expect_identical(proto, double())
})

test_that("each bound parameter is compared when weight-truncated psw objects combine", {
  w <- wt_truncated_psw()

  other_method <- wt_truncated_psw(
    record = wt_trunc_record(method = "count")
  )
  expect_bound_mismatch(c(w, other_method))

  lower_bounded <- wt_truncated_psw(
    record = new_psw_trunc_meta(
      method = "wt",
      lower = 0.5,
      upper = 3,
      lower_value = 0.5,
      upper_value = 3,
      truncated_idx = c(2L, 5L),
      n_obs = 5L
    )
  )
  expect_bound_mismatch(c(w, lower_bounded))
})

test_that("a bound mismatch is caught whichever pair of the fold meets it", {
  # vctrs folds the prototype over the inputs two at a time, and the prototype
  # carries no positional record, so the bound the earlier inputs agreed on has
  # to reach the pair that meets the disagreeing one.
  w <- wt_truncated_psw()
  alt <- alt_wt_truncated_psw()

  out <- expect_bound_mismatch(c(w, w, alt))
  expect_false(is_psw(out))

  # A slice has no record, so its pair agrees and passes the bound on.
  out <- expect_bound_mismatch(c(w[1:2], w, alt))
  expect_false(is_psw(out))

  agreed <- expect_silent(c(w, w[1:2], w))
  expect_s3_class(agreed, "psw")
  expect_true(is_wt_truncated(agreed))
})

test_that("weights bounded differently do not cast to each other", {
  # Subassignment casts the value to the target and keeps the target's record,
  # so a value bounded differently would be written under a bound that does not
  # describe it. A prototype carries its inputs' bound without positions, so
  # weights still cast to their own.
  w <- wt_truncated_psw()
  alt <- alt_wt_truncated_psw()

  expect_error(vec_cast(alt, w), class = "vctrs_error_cast")
  expect_error(
    {
      w[seq_along(w)] <- alt
    },
    class = "vctrs_error_cast"
  )
  expect_identical(attr(w, "psw_trunc_meta"), wt_trunc_record())
  expect_true(is_wt_truncated(w))

  expect_silent(vec_cast(w, vec_ptype(w)))
})

test_that("weights bounded alike combine whichever units the bound moved", {
  # The comparison is of the bound, not of the units it moved: two samples
  # truncated at the same value have moved different units and still describe
  # the same operation.
  w <- wt_truncated_psw()
  other_units <- wt_truncated_psw(
    x = c(3, 1.1, 2.2, 3, 0.9, 1.4),
    record = wt_trunc_record(truncated_idx = c(1L, 4L), n_obs = 6L)
  )

  out <- expect_silent(c(w, other_units))
  expect_s3_class(out, "psw")
  expect_length(out, 11)
  expect_true(is_wt_truncated(out))
  expect_positions_dropped(attr(out, "psw_trunc_meta"), wt_trunc_record())
})

test_that("a weight-truncated psw with no record agrees with any bound", {
  # Weights flagged as truncated with no record, as an earlier combine leaves
  # them, have no bound to disagree with, so only the flag is compared.
  w <- wt_truncated_psw()
  alt_slice <- psw(c(1.2, 2.5), estimand = "ate", wt_truncated = TRUE)
  expect_null(attr(alt_slice, "psw_trunc_meta"))
  expect_true(is_wt_truncated(alt_slice))

  out <- expect_silent(c(w, alt_slice))
  expect_s3_class(out, "psw")
  expect_length(out, 7)
  expect_true(is_wt_truncated(out))
  expect_positions_dropped(attr(out, "psw_trunc_meta"), wt_trunc_record())

  reversed <- expect_silent(c(alt_slice, w))
  expect_s3_class(reversed, "psw")
  expect_true(is_wt_truncated(reversed))
  expect_positions_dropped(attr(reversed, "psw_trunc_meta"), wt_trunc_record())

  # The prototype takes the one record there is.
  bare <- psw(c(1, 2), estimand = "ate", wt_truncated = TRUE)
  proto <- expect_silent(vec_ptype2(bare, w))
  expect_s3_class(proto, "psw")
  expect_true(is_wt_truncated(proto))
  expect_positions_dropped(attr(proto, "psw_trunc_meta"), wt_trunc_record())
})

test_that("is_wt_truncated() answers FALSE for anything that is not a psw", {
  expect_false(is_wt_truncated(c(1, 2)))
  expect_false(is_wt_truncated(NULL))

  scores <- ps_trunc(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  expect_false(is_wt_truncated(scores))
})

test_that("is_unit_wt_truncated() refuses anything that is not a psw", {
  expect_error(
    is_unit_wt_truncated(c(1, 2)),
    class = "propensity_method_error"
  )
  expect_propensity_error(is_unit_wt_truncated(c(1, 2)))
})
