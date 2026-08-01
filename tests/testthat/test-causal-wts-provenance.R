# The abstract `causal_wts` weight layer is owned by causalgenerics. `psw`
# inherits from `causal_wts`, so every method that layer supplies must resolve
# to causalgenerics, and the behavior a `psw` user sees must be identical
# whichever package supplies the method.

# The namespace whose S3 methods table records registrations for each generic.
# Reading these tables reports which package owns a method independently of what
# is on the search path. It does not, on its own, prove a method is gone:
# `UseMethod()` searches the environment the generic was called from as well as
# the registration table, so a method left in `asNamespace("propensity")` but no
# longer registered still serves every call made from inside propensity. The
# absence checks below therefore pair a table read with a namespace read.
generic_namespaces <- c(
  estimand = "causalgenerics",
  `estimand<-` = "causalgenerics",
  is_causal_wt = "causalgenerics",
  vec_math = "vctrs",
  median = "stats",
  quantile = "stats",
  Summary = "base",
  min = "base",
  max = "base",
  range = "base",
  summary = "base",
  anyDuplicated = "base",
  diff = "base",
  `[` = "base",
  `==` = "base",
  `!=` = "base",
  `<` = "base",
  `>` = "base",
  `<=` = "base",
  `>=` = "base",
  # Generics behind the methods propensity must keep. `vec_arith.psw` is
  # propensity's own second-level dispatcher, so its methods register in
  # propensity's table rather than in vctrs'.
  vec_restore = "vctrs",
  vec_arith = "vctrs",
  `vec_arith.numeric` = "vctrs",
  `vec_arith.psw` = "propensity",
  vec_cast = "vctrs",
  vec_ptype2 = "vctrs",
  vec_ptype_abbr = "vctrs",
  vec_ptype_full = "vctrs",
  # The trimming, truncation, and calibration accessors are propensity's own
  # generics, so their methods register in propensity's table too.
  is_ps_trimmed = "propensity",
  is_unit_trimmed = "propensity",
  is_ps_truncated = "propensity",
  is_unit_truncated = "propensity",
  is_refit = "propensity",
  is_ps_calibrated = "propensity"
)

# Package that owns the method registered for one generic and class, or NA when
# no such method is registered.
method_owner <- function(generic, cls) {
  methods_table <- get(
    ".__S3MethodsTable__.",
    envir = asNamespace(generic_namespaces[[generic]]),
    inherits = FALSE
  )
  method <- get0(
    paste0(generic, ".", cls),
    envir = methods_table,
    inherits = FALSE
  )
  if (is.null(method)) {
    return(NA_character_)
  }
  environmentName(environment(method))
}

# Package that owns the method S3 dispatch would select for `object`.
dispatch_owner <- function(generic, object) {
  for (cls in class(object)) {
    owner <- method_owner(generic, cls)
    if (!is.na(owner)) {
      return(owner)
    }
  }
  NA_character_
}

# Owners for a set of methods given as classes named by their generic, keyed by
# full method name so a failure names the method that moved.
method_owners <- function(methods) {
  owners <- vapply(
    seq_along(methods),
    function(i) method_owner(names(methods)[[i]], methods[[i]]),
    character(1)
  )
  stats::setNames(owners, paste0(names(methods), ".", methods))
}

# Accessors the abstract layer registers directly on `causal_wts`.
causal_wts_accessors <- c("estimand", "estimand<-", "is_causal_wt")

# Methods the abstract layer supplies on `causal_wts` and that a `psw` object
# reaches by inheritance rather than through a method of its own.
psw_layer_generics <- c(
  "vec_math",
  "Summary",
  "min",
  "max",
  "range",
  "median",
  "quantile",
  "summary",
  "anyDuplicated",
  "diff",
  "[",
  "==",
  "!=",
  "<",
  ">",
  "<=",
  ">="
)

# Methods propensity must keep on `psw`. `vec_ptype2`, `vec_cast`,
# `vec_ptype_abbr`, and `vec_ptype_full` cannot be inherited at all: vctrs
# resolves them on `class(x)[[1]]` and never walks the class vector.
# `vec_restore` and `vec_arith` do dispatch through `causal_wts`, but no
# abstract implementation of either preserves propensity's behavior. The
# trimming, truncation, calibration, and refit accessors are propensity's own
# and have no upstream counterpart. All of these sit interleaved with the
# inherited methods in `R/psw.R`, so an over-wide deletion is the failure mode
# this pins. `is_unit_truncated.psw` in particular is deletable in silence: no
# other test passes a `psw` to it, so its loss only shows up as
# `is_unit_truncated.default` aborting in a user's session.
psw_retained_methods <- c(
  vec_restore = "psw",
  vec_arith = "psw",
  `vec_arith.numeric` = "psw",
  `vec_arith.psw` = "MISSING",
  `vec_arith.psw` = "default",
  `vec_arith.psw` = "integer",
  `vec_arith.psw` = "numeric",
  `vec_arith.psw` = "psw",
  vec_cast = "character.psw",
  vec_cast = "double.psw",
  vec_cast = "integer.psw",
  vec_cast = "ps_trim.psw",
  vec_cast = "ps_trunc.psw",
  vec_cast = "psw.double",
  vec_cast = "psw.integer",
  vec_cast = "psw.ps_trim",
  vec_cast = "psw.ps_trunc",
  vec_cast = "psw.psw",
  vec_ptype2 = "character.psw",
  vec_ptype2 = "double.psw",
  vec_ptype2 = "integer.psw",
  vec_ptype2 = "ps_calib.psw",
  vec_ptype2 = "ps_trim.psw",
  vec_ptype2 = "ps_trunc.psw",
  vec_ptype2 = "psw.character",
  vec_ptype2 = "psw.double",
  vec_ptype2 = "psw.integer",
  vec_ptype2 = "psw.ps_calib",
  vec_ptype2 = "psw.ps_trim",
  vec_ptype2 = "psw.ps_trunc",
  vec_ptype2 = "psw.psw",
  vec_ptype_abbr = "psw",
  vec_ptype_full = "psw",
  is_ps_trimmed = "psw",
  is_unit_trimmed = "psw",
  is_ps_truncated = "psw",
  is_unit_truncated = "psw",
  is_refit = "psw",
  is_ps_calibrated = "psw"
)

# Provenance ----

test_that("causalgenerics supplies the abstract layer on causal_wts", {
  owners <- vapply(
    psw_layer_generics,
    method_owner,
    character(1),
    cls = "causal_wts"
  )

  expect_identical(
    owners,
    stats::setNames(
      rep("causalgenerics", length(psw_layer_generics)),
      psw_layer_generics
    )
  )
})

test_that("the causal_wts accessors resolve to causalgenerics", {
  owners <- vapply(
    causal_wts_accessors,
    method_owner,
    character(1),
    cls = "causal_wts"
  )

  expect_identical(
    owners,
    stats::setNames(
      rep("causalgenerics", length(causal_wts_accessors)),
      causal_wts_accessors
    )
  )
})

test_that("propensity registers no psw-level copy of the abstract layer", {
  owners <- vapply(psw_layer_generics, method_owner, character(1), cls = "psw")

  expect_identical(
    owners,
    stats::setNames(
      rep(NA_character_, length(psw_layer_generics)),
      psw_layer_generics
    )
  )
})

test_that("no definition of the moved layer survives in propensity", {
  # An unregistered method still defined in the namespace is invisible to the
  # registration tables and to `R CMD check`, which only errors the other way
  # round, on a declared-but-missing method. Every internal call site would
  # keep reaching propensity's copy while users reach causalgenerics'.
  candidates <- c(
    paste0(psw_layer_generics, ".psw"),
    paste0(causal_wts_accessors, ".causal_wts"),
    "psw_compare"
  )
  survivors <- vapply(
    candidates,
    exists,
    logical(1),
    envir = asNamespace("propensity"),
    inherits = FALSE
  )

  expect_identical(candidates[survivors], character())
})

test_that("propensity still owns the psw methods that cannot move upstream", {
  owners <- method_owners(psw_retained_methods)

  expect_identical(
    owners,
    stats::setNames(rep("propensity", length(owners)), names(owners))
  )
})

test_that("a psw object dispatches to causalgenerics for the abstract layer", {
  w <- psw(c(1.2, 0.8, 1.5), estimand = "ate")
  owners <- vapply(
    psw_layer_generics,
    dispatch_owner,
    character(1),
    object = w
  )

  expect_identical(
    owners,
    stats::setNames(
      rep("causalgenerics", length(psw_layer_generics)),
      psw_layer_generics
    )
  )
})

# Behavior that must survive the change of owner ----

test_that("subsetting a psw preserves trimming, truncation, and calibration", {
  w <- psw(
    c(1.2, 0.8, 1.5, 0.4, 2),
    estimand = "ate",
    trimmed = TRUE,
    truncated = TRUE,
    calibrated = TRUE,
    stabilization_score = 0.37
  )

  metadata <- function(x) {
    list(
      trimmed = is_ps_trimmed(x),
      truncated = is_ps_truncated(x),
      calibrated = is_ps_calibrated(x),
      score = stabilization_score(x)
    )
  }
  expected <- list(
    trimmed = TRUE,
    truncated = TRUE,
    calibrated = TRUE,
    score = 0.37
  )

  expect_equal(metadata(w[c(2, 4)]), expected)
  expect_equal(metadata(w[c(TRUE, FALSE, TRUE, FALSE, TRUE)]), expected)
  expect_equal(metadata(w[-1]), expected)
  expect_equal(metadata(w[integer(0)]), expected)
  expect_equal(metadata(w[]), expected)
})

test_that("subsetting a trimmed psw drops the unit-level trimming index", {
  set.seed(42)
  n <- 60
  z <- rnorm(n)
  x <- rbinom(n, 1, plogis(2 * z))
  fit <- glm(x ~ z, family = binomial())
  trimmed <- ps_trim(
    predict(fit, type = "response"),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  wt <- wt_ate(
    ps_refit(trimmed, model = fit),
    .exposure = x,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_true(is_ps_trimmed(wt))
  expect_length(is_unit_trimmed(wt), n)

  # The trimming index is written against the full vector and nothing rebuilding
  # a psw is handed the subscript, so a subset cannot re-index it and the index
  # is dropped. The drop is silent because `model.frame()` shortens these weights
  # inside every outcome model fit on them, with no subscript in the user's code.
  # What the assertion guards on the other side is the blind-attribute-copy
  # regression, which would re-attach an index built against the original vector
  # and return a length-`n` logical for a length-2 subset. With no index left,
  # `is_unit_trimmed()` reports that it cannot answer rather than reporting every
  # unit as retained.
  first_two <- expect_silent(wt[1:2])
  expect_true(is_ps_trimmed(first_two))
  expect_null(ps_trim_meta(first_two))
  expect_error(
    is_unit_trimmed(first_two),
    class = "propensity_missing_meta_error"
  )
})

test_that("the causal_wts accessors read and write the estimand attribute", {
  w <- psw(c(1.2, 0.8), estimand = "ate")

  expect_identical(estimand(w), "ate")
  expect_true(is_causal_wt(w))
  expect_false(is_causal_wt(c(1.2, 0.8)))
  expect_null(estimand(psw(1)))

  estimand(w) <- "att"
  expect_identical(estimand(w), "att")
  expect_identical(attr(w, "estimand"), "att")
  expect_s3_class(w, "psw")
})

test_that("psw sum(), min(), max(), and range() honor na.rm", {
  w <- psw(c(1.2, NA, 0.5), estimand = "ate")

  expect_identical(sum(w), NA_real_)
  expect_identical(min(w), NA_real_)
  expect_identical(max(w), NA_real_)
  expect_identical(range(w), c(NA_real_, NA_real_))

  expect_equal(sum(w, na.rm = TRUE), 1.7)
  expect_equal(min(w, na.rm = TRUE), 0.5)
  expect_equal(max(w, na.rm = TRUE), 1.2)
  expect_equal(range(w, na.rm = TRUE), c(0.5, 1.2))

  expect_false(is_psw(sum(w, na.rm = TRUE)))
  expect_false(is_psw(range(w, na.rm = TRUE)))
})

test_that("psw median() and quantile() return plain numerics", {
  w <- psw(c(1.2, 0.8, 1.5, 0.4, 2), estimand = "ate", stabilized = TRUE)

  med <- median(w)
  expect_false(is_psw(med))
  expect_type(med, "double")
  expect_identical(med, 1.2)

  quartiles <- quantile(w)
  expect_false(is_psw(quartiles))
  expect_type(quartiles, "double")
  expect_identical(names(quartiles), c("0%", "25%", "50%", "75%", "100%"))
  expect_equal(unname(quartiles), c(0.4, 0.8, 1.2, 1.5, 2))

  w_na <- psw(c(1.2, NA, 0.5), estimand = "ate")
  expect_equal(median(w_na, na.rm = TRUE), 0.85)
  expect_equal(
    unname(quantile(w_na, probs = c(0.25, 0.75), na.rm = TRUE)),
    c(0.675, 1.025)
  )
})

test_that("psw anyDuplicated() returns an integer index, not a logical", {
  repeated <- psw(c(1, 2, 2, 3), estimand = "ate")
  index <- anyDuplicated(repeated)
  expect_type(index, "integer")
  expect_identical(index, 3L)

  distinct <- psw(c(1, 2, 3), estimand = "ate")
  expect_identical(anyDuplicated(distinct), 0L)
})

test_that("psw diff() and summary() return plain numerics", {
  w <- psw(c(1, 3, 6, 10), estimand = "ate")

  differences <- diff(w)
  expect_false(is_psw(differences))
  expect_equal(differences, c(2, 3, 4))
  expect_equal(diff(w, lag = 2), c(5, 7))

  overview <- summary(w)
  expect_false(is_psw(overview))
  expect_s3_class(overview, "summaryDefault")
  expect_equal(unname(overview[["Median"]]), 4.5)
})

test_that("every psw comparison operator returns logical under vctrs sizing", {
  # `test-bracket-psw.R` pins the returned values and both dispatch sides for
  # these operators. What is added here is uniformity across all six, including
  # the three that file does not exercise for size errors.
  w <- psw(c(1.2, 0.8, 1.5, 0.4), estimand = "ate")
  shorter <- psw(c(1, 2), estimand = "ate")

  for (op_name in c("==", "!=", "<", ">", "<=", ">=")) {
    op <- match.fun(op_name)
    out <- op(w, 1)
    expect_type(out, "logical")
    expect_length(out, length(w))
    expect_error(op(w, shorter), class = "vctrs_error_incompatible_size")
  }
})

test_that("psw cumulative math keeps the class and metadata but sqrt() does not", {
  w <- psw(
    c(1, 4, 9),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.37
  )

  expected <- list(
    cumsum = c(1, 5, 14),
    cumprod = c(1, 4, 36),
    cummin = c(1, 1, 1),
    cummax = c(1, 4, 9)
  )
  for (fn_name in names(expected)) {
    cumulative <- match.fun(fn_name)(w)
    expect_s3_class(cumulative, "psw")
    expect_identical(estimand(cumulative), "ate")
    expect_true(is_stabilized(cumulative))
    expect_equal(stabilization_score(cumulative), 0.37)
    expect_equal(as.numeric(cumulative), expected[[fn_name]])
  }

  roots <- sqrt(w)
  expect_false(is_psw(roots))
  expect_type(roots, "double")
  expect_equal(roots, c(1, 2, 3))
})
