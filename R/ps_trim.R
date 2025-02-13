#' Trim Propensity Scores With NAs (Preserving Original Length)
#'
#' `ps_trim()` applies various trimming methods to a propensity-score vector,
#' returning a new vector of the *same length*, with trimmed entries replaced by NA.
#' The object is given class `"ps_trim"`, and you can inspect further metadata in
#' `ps_trim_meta(x)`. After running `ps_trim()`, you should refit the model with
#' `ps_refit()`.
#'
#' @param ps A numeric vector in (0,1).
#' @param exposure For methods like `"pref"` or `"cr"`, a 0/1 vector.
#' @param method One of `c("ps", "adaptive", "pctl", "pref", "cr")`.
#' @param lower,upper Numeric cutoffs or quantiles. If `NULL`, defaults vary by method.
#'
#' @details
#' The returned object is a **`ps_trim`** vector of the same length as `ps`, but
#' with trimmed entries replaced by `NA`.
#' An attribute `ps_trim_meta` contains:
#'
#' - `method`: Which trimming method was used
#' - `keep_idx`: Indices retained
#' - `trimmed_idx`: Indices replaced by `NA`
#' - Possibly other fields such as final cutoffs, etc.
#'
#' @return
#' A `ps_trim` object (numeric vector). The attribute `ps_trim_meta` stores metadata.
#'
#' @seealso [ps_trunc()] for bounding/winsorizing instead of discarding,
#'   [is_refit()], [is_trimmed()]
#'
#' @export
ps_trim <- function(
  ps,
  exposure = NULL,
  method = c("ps", "adaptive", "pctl", "pref", "cr"),
  lower = NULL,
  upper = NULL
) {
  method <- match.arg(method)

  if (any(ps <= 0 | ps >= 1)) {
    stop("All propensity scores must be strictly between 0 and 1.")
  }

  # Possibly set defaults or warn
  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    if (lower >= upper) {
      stop("For method='ps', need lower < upper. Got lower=", lower, ", upper=", upper)
    }
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warning("For method='adaptive', `lower`/`upper` are ignored.")
    }
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
  } else if (method == "pref") {
    if (is.null(lower)) lower <- 0.3
    if (is.null(upper)) upper <- 0.7
  } else {
    # "cr"
    if (!is.null(lower) || !is.null(upper)) {
      warning("For method='cr', `lower`/`upper` are ignored.")
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
    if (is.null(exposure)) {
      stop("For method='pref', must supply a binary `exposure`.")
    }
    if (!all(exposure %in% c(0, 1))) {
      stop("'exposure' must be 0/1 for preference scores.")
    }
    P <- mean(exposure)
    if (P <= 0 || P >= 1) {
      stop("Proportion of exposure is 0 or 1; cannot compute preference score.")
    }
    pref_score <- plogis(qlogis(ps) - qlogis(P))
    meta_list$P <- P
    meta_list$pref_formula <- "pref = expit(logit(ps) - logit(P))"
    keep_idx <- which(pref_score >= lower & pref_score <= upper)
  } else {
    # cr
    if (is.null(exposure)) {
      stop("For method='cr', must supply a binary `exposure`.")
    }
    if (!all(exposure %in% c(0, 1))) {
      stop("'exposure' must be 0/1 for common range method.")
    }
    if (all(exposure == 0) || all(exposure == 1)) {
      stop("Proportion of exposure is 0 or 1; cannot compute common range.")
    }
    ps_treat <- ps[exposure == 1]
    ps_untrt <- ps[exposure == 0]
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
  vctrs::vec_assert(x, ptype = double())
  vctrs::new_vctr(
    x,
    ps_trim_meta = ps_trim_meta,
    class = "ps_trim"
  )
}

#' Check if object is trimmed
#'
#' @description `is_trimmed()` is an S3 generic that returns `TRUE` if its
#' argument represents a ps_trim object or psw object flagged as trimmed.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_trimmed <- function(x) {
  UseMethod("is_trimmed")
}

#' @export
is_trimmed.default <- function(x) FALSE

is_trimmed <- function(x) {
  UseMethod("is_trimmed")
}

#' @export
is_trimmed.default <- function(x) {
  FALSE
}

#' @export
is_trimmed.ps_trim <- function(x) {
  TRUE
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
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim ps_trim
vec_arith.ps_trim.ps_trim <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim numeric
vec_arith.ps_trim.numeric <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trim
vec_arith.numeric.ps_trim <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim integer
vec_arith.ps_trim.integer <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
vec_ptype2.ps_trim.ps_trim <- function(x, y, ...) {
  vctrs::stop_incompatible_type(
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
  vctrs::vec_data(x)
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
vec_ptype2.ps_trim.integer <- function(x, y, ...) integer()
#' @export
vec_ptype2.integer.ps_trim <- function(x, y, ...) integer()

#' @export
vec_cast.integer.ps_trim <- function(x, to, ...) {
  as.integer(vctrs::vec_data(x))
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

#' Refit the PS Model on Retained Observations
#'
#' Takes a `ps_trim` object (with `NA` where trimmed) and the original model
#' used to get the PS, then:
#' 1. Retrieves data from the model (or from `data` argument if provided)
#' 2. Subsets rows to the non‐trimmed indices
#' 3. Refits the model
#' 4. Predicts new propensity scores for all rows (trimmed rows -> NA)
#' 5. Returns a new `ps_trim` object with `refit=TRUE`.
#'
#' @param ps_trim_obj A `ps_trim` object (same length as data, NAs for trimmed).
#' @param model The fitted model used to get the original PS (e.g. a glm).
#' @param data Optional. A data frame. If `NULL`, we try to retrieve from `model`.
#' @param ... Additional arguments passed to `update()`.
#'
#' @return
#' A new `ps_trim` object with updated propensity scores and
#' `ps_trim_meta(x)$refit` set to `TRUE`.
#'
#' @seealso [ps_trim()], [is_refit()], [is_trimmed()]
#'
#' @export
ps_refit <- function(ps_trim_obj, model, data = NULL, ...) {
  if (!inherits(ps_trim_obj, "ps_trim")) {
    stop("`ps_refit()` expects a `ps_trim` object.")
  }
  meta <- attr(ps_trim_obj, "ps_trim_meta")
  keep_idx <- meta$keep_idx

  if (length(keep_idx) == 0) {
    stop("No retained rows to refit on (all were trimmed).")
  }

  # 1) If data is not provided, attempt to retrieve from the model
  if (is.null(data)) {
    data <- tryCatch(
      stats::model.frame(model),
      error = function(e) {
        stop(
          "Couldn't retrieve data from `model`. ",
          "Please supply `data` explicitly."
        )
      }
    )
  }

  # 2) Check row counts
  if (nrow(data) != length(ps_trim_obj)) {
    stop("`data` must have the same number of rows as `length(ps_trim_obj)`.")
  }

  data_sub <- data[keep_idx, , drop = FALSE]

  # 3) Refit
  # For a typical glm, `update(model, data=data_sub, ...)` often works
  refit_model <- stats::update(model, data = data_sub, ...)

  # 4) Predict new PS for all rows
  new_ps <- rep(NA_real_, length(ps_trim_obj))
  new_ps[keep_idx] <- stats::predict(refit_model, newdata = data_sub, type = "response")

  # 5) Create a brand-new ps_trim object with refit=TRUE
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
