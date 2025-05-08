#' Trim Propensity Scores
#'
#' `ps_trim()` applies trimming methods to a propensity-score vector, returning
#' a new vector of the *same length*, with trimmed entries replaced by `NA.` You
#' can inspect further metadata in `ps_trim_meta(x)`. After running `ps_trim()`,
#' you should refit the model with `ps_refit()`.
#'
#' @param ps The propensity score, a numeric vector between 0 and 1.
#' @param .exposure For methods like `"pref"` or `"cr"`, a vector for a binary
#'   exposure.
#' @param method One of `c("ps", "adaptive", "pctl", "pref", "cr")`.
#' @param lower,upper Numeric cutoffs or quantiles. If `NULL`, defaults vary by
#'   method.
#' @inheritParams wt_ate
#'
#' @details The returned object is a **`ps_trim`** vector of the same length as
#' `ps`, but with trimmed entries replaced by `NA`. An attribute `ps_trim_meta`
#' contains:
#'
#' - `method`: Which trimming method was used
#' - `keep_idx`: Indices retained
#' - `trimmed_idx`: Indices replaced by `NA`
#' - Possibly other fields such as final cutoffs, etc.
#'
#' @return A `ps_trim` object (numeric vector). The attribute `ps_trim_meta`
#' stores metadata.
#'
#' @seealso [ps_trunc()] for bounding/winsorizing instead of discarding,
#'   [is_refit()], [is_ps_trimmed()]
#' @examples
#'
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.3 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' ps_trim(ps, method = "adaptive")
#'
#' @export
ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  method <- rlang::arg_match(method)
  check_ps_range(ps)

  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    check_lower_upper(lower, upper)
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warn("For {.code method = 'adaptive'}, {.code lower} and {.code upper} are ignored.")
    }
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
  } else if (method == "pref") {
    if (is.null(lower)) lower <- 0.3
    if (is.null(upper)) upper <- 0.7
  } else {
    if (!is.null(lower) || !is.null(upper)) {
      warn("For {.code method = 'cr'}, {.code lower} and {.code upper} are ignored.")
    }
  }

  n <- length(ps)
  keep_idx <- integer(0)
  trimmed_idx <- integer(0)
  meta_list <- list(method = method, lower = lower, upper = upper)

  # Decide which indices are kept
  if (method == "ps") {
    keep_idx <- which(ps >= lower & ps <= upper)
  } else if (method == "adaptive") {
    sum_wt <- 1 / (ps * (1 - ps))
    k <- 2 * mean(sum_wt) - max(sum_wt)

    if (k >= 0) {
      cutoff <- 0
    } else {
      trim_fun <- function(x) {
        sum_wt_trim <- sum_wt[sum_wt <= x]
        x - 2 * mean(sum_wt_trim)
      }
      rng <- range(sum_wt)
      lambda <- uniroot(trim_fun, lower = rng[1], upper = rng[2])$root
      cutoff <- 0.5 - sqrt(0.25 - 1 / lambda)
    }
    meta_list$cutoff <- cutoff
    keep_idx <- which(pmin(ps, 1 - ps) > cutoff)
  } else if (method == "pctl") {
    q_lower <- quantile(ps, probs = lower)
    q_upper <- quantile(ps, probs = upper)
    meta_list$q_lower <- q_lower
    meta_list$q_upper <- q_upper
    keep_idx <- which(ps >= q_lower & ps <= q_upper)
  } else if (method == "pref") {
    if (is.null(.exposure)) {
      abort("For {.code method = 'pref'}, must supply {.arg exposure}.")
    }
    .exposure <- transform_exposure_binary(
      .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
    prop_exposure <- mean(.exposure)
    pref_score <- plogis(qlogis(ps) - qlogis(prop_exposure))
    meta_list$P <- prop_exposure
    keep_idx <- which(pref_score >= lower & pref_score <= upper)
  } else {
    if (is.null(.exposure)) {
      abort("For {.code method = 'cr'}, must supply {.arg exposure}.")
    }
    .exposure <- transform_exposure_binary(
      .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
    ps_treat <- ps[.exposure == 1]
    ps_untrt <- ps[.exposure == 0]
    cr_lower <- min(ps_treat)
    cr_upper <- max(ps_untrt)
    meta_list$cr_lower <- cr_lower
    meta_list$cr_upper <- cr_upper

    keep_idx <- which(ps >= cr_lower & ps <= cr_upper)
  }

  trimmed_idx <- setdiff(seq_len(n), keep_idx)

  # Replace trimmed entries with NA
  ps_na <- ps
  ps_na[trimmed_idx] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx    = keep_idx,
        trimmed_idx = trimmed_idx
      )
    )
  )
}

new_trimmed_ps <- function(x, ps_trim_meta = list()) {
  vec_assert(x, ptype = double())
  new_vctr(
    x,
    ps_trim_meta = ps_trim_meta,
    class = "ps_trim",
    inherit_base_type = TRUE
  )
}

#' Check if object is trimmed
#'
#' @description `is_ps_trimmed()` is an S3 generic that returns `TRUE` if its
#'   argument represents a `ps_trim` object or `psw` object created from trimmed
#'   propensity scores. `is_ps_trimmed()` is a question about whether or not the
#'   propensity scores *have* been trimmed, as opposed to [is_unit_trimmed()],
#'   which is a question about which *units* have been trimmed.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_ps_trimmed <- function(x) {
  UseMethod("is_ps_trimmed")
}

#' @export
is_ps_trimmed.default <- function(x) FALSE

is_ps_trimmed <- function(x) {
  UseMethod("is_ps_trimmed")
}

#' @export
is_ps_trimmed.default <- function(x) {
  FALSE
}

#' @export
is_ps_trimmed.ps_trim <- function(x) {
  TRUE
}

#' Check if units have been trimmed
#'
#' @description `is_unit_trimmed()` is an that vector of `TRUE` or `FALSE`
#'   values, representing if the unit was trimmed. `is_unit_trimmed()` is a
#'   question about which *units* have been trimmed, as opposed to
#'   [is_ps_trimmed()], which is a question about whether or not the propensity
#'   scores *have* been trimmed.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_unit_trimmed <- function(x) {
  UseMethod("is_unit_trimmed")
}

#' @export
is_unit_trimmed.default <- function(x) {
  abort(
    "{.code is_unit_trimmed()} not supported for class {.val {class(x)}}"
  )
}

#' @export
is_unit_trimmed.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  out <- vector("logical", length = length(x))
  out[meta$trimmed_idx] <- TRUE

  out
}


# vctrs machinery for ps_trim

#' @export
vec_ptype_abbr.ps_trim <- function(x, ...) "ps_trim"

#' @export
vec_ptype_full.ps_trim <- function(x, ...) {
  paste(
    "ps_trim;",
    "trimmed",
    length(ps_trim_meta(x)$trimmed_idx),
    "of "
  )
}

#' @export
#' @method vec_arith ps_trim
vec_arith.ps_trim <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trim", y)
}

#' @export
#' @method vec_arith.ps_trim default
vec_arith.ps_trim.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim ps_trim
vec_arith.ps_trim.ps_trim <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim numeric
vec_arith.ps_trim.numeric <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trim
vec_arith.numeric.ps_trim <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim integer
vec_arith.ps_trim.integer <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
vec_ptype2.ps_trim.ps_trim <- function(x, y, ...) {
  stop_incompatible_type(
    x,
    y,
    x_arg = "x",
    y_arg = "y",
    message = "Can't combine two ps_trim objects"
  )
}
#' @export
vec_ptype2.ps_trim.double <- function(x, y, ...) double()
#' @export
vec_ptype2.double.ps_trim <- function(x, y, ...) double()

#' @export
vec_cast.double.ps_trim <- function(x, to, ...) {
  # degrade to numeric with NAs
  vec_data(x)
}

#' @export
vec_cast.ps_trim.double <- function(x, to, ...) {
  # create a default ps_trim with no trimming
  new_trimmed_ps(
    x,
    ps_trim_meta = list(
      method      = "unknown",
      keep_idx    = seq_along(x),
      trimmed_idx = integer(0)
    )
  )
}

#' @export
vec_ptype2.psw.ps_trim <- function(x, y, ...) character()

#' @export
vec_ptype2.ps_trim.psw <- function(x, y, ...) character()

#' @export
vec_cast.character.ps_trim <- function(x, to, ...) as.character(vec_data(x))

#' @export
vec_ptype2.ps_trim.integer <- function(x, y, ...) integer()
#' @export
vec_ptype2.integer.ps_trim <- function(x, y, ...) integer()

#' @export
vec_cast.integer.ps_trim <- function(x, to, ...) {
  as.integer(vec_data(x))
}

#' @export
vec_cast.ps_trim.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_trimmed_ps(
    xx,
    ps_trim_meta = list(
      method      = "unknown",
      keep_idx    = seq_along(xx),
      trimmed_idx = integer(0)
    )
  )
}

#' Refit the Propensity Score Model on Retained Observations
#'
#' Takes a `ps_trim` object and the original model
#' used to calculate the propensity score, then:
#' 1. Retrieves data from the model (or from `.df` argument if provided)
#' 2. Subsets rows to the non‐trimmed indices
#' 3. Refits the model
#' 4. Predicts new propensity scores for all rows (trimmed rows -> `NA`)
#' 5. Returns a new `ps_trim` object with `refit = TRUE`.
#'
#' @param trimmed_ps A `ps_trim` object (same length as data, NAs for trimmed).
#' @param model The fitted model used to get the original PS (e.g. a glm).
#' @param .df Optional. A data frame. If `NULL`, we try to retrieve from `model`.
#' @param ... Additional arguments passed to `update()`.
#'
#' @return
#' A new `ps_trim` object with updated propensity scores and
#' `ps_trim_meta(x)$refit` set to `TRUE`.
#'
#' @seealso [ps_trim()], [is_refit()], [is_ps_trimmed()]
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # trim and refit
#' refit <- ps_trim(ps, lower = .2, upper = .8) |>
#'   ps_refit(fit)
#'
#' is_refit(refit)
#'
#' @export
ps_refit <- function(trimmed_ps, model, .df = NULL, ...) {
  assert_class(trimmed_ps, "ps_trim")
  meta <- ps_trim_meta(trimmed_ps)

  if (length(meta$keep_idx) == 0) {
    abort("No retained rows to refit on (all were trimmed).")
  }

  if (is.null(.df)) {
    .df <- model.frame(model)
  }

  if (nrow(.df) != length(trimmed_ps)) {
    abort(c(
      "{.arg .df} must have the same number of rows as \\
      {.code length(trimmed_ps)}.",
      x = "{.arg .df} has {nrow(.df)} row{?s}.",
      x = "{.arg trimmed_ps} has length {length(trimmed_ps)}."
    ))
  }

  # refit on untrimmed rows
  data_sub <- .df[meta$keep_idx, , drop = FALSE]
  refit_model <- stats::update(model, data = data_sub, ...)

  # predict new PS for all rows
  new_ps <- rep(NA_real_, length(trimmed_ps))
  new_ps[meta$keep_idx] <- stats::predict(
    refit_model,
    newdata = data_sub,
    type = "response"
  )

  meta$refit <- TRUE

  new_trimmed_ps(
    x = new_ps,
    ps_trim_meta = meta
  )
}

#' Check if an object has been refit
#'
#' @description
#' **`is_refit()`** is an S3 generic that returns `TRUE` if its argument
#' represents a [ps_trim] object (or a weighting object) that has had the
#' propensity model refit on the retained subset.
#'
#' @param x An R object (e.g. a [ps_trim] or [psw]).
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_refit <- function(x) {
  UseMethod("is_refit")
}

#' @export
is_refit.default <- function(x) {
  FALSE
}

#' @export
is_refit.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  isTRUE(meta$refit)
}

#' @title Extract `ps_trim` metadata
#' @description Returns the internal metadata list for a `ps_trim` object.
#' @param x A **`ps_trim`** object.
#' @return A named list of metadata.
#' @export
ps_trim_meta <- function(x) {
  attr(x, "ps_trim_meta")
}
