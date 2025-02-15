#' Calculate propensity score weights
#'
#' @param .propensity Either a vector of the predicted value of `.exposure` or a
#'   `data.frame` where each column is the predicted probability of a level of
#'   `.exposure`.
#' @param .exposure The exposure for which `.propensity` is calculated.
#' @param exposure_type The type of exposure. By default, automatically detected
#'   based on `.exposure`.
#' @param .sigma If `exposure_type` is continuous, a vector of observation-level
#'   standard errors passed to `dnorm()`. For an `lm` model this is
#'   `influence(model)$sigma`. For data frames produced by broom's `augment()`,
#'   this is the `.sigma` column.
#' @param .treated The treatment level of the exposure. Automatically detected
#'   by default.
#' @param .untreated The control level of the exposure. Automatically detected
#'   by default.
#' @param ... Passed to other functions Not currently used.
#' @param stabilize Logical. Stabilize the weights? By default, stabilizes with
#'   the mean of `.exposure`.
#' @param stabilization_score if `stabilize` is `TRUE`, optionally include a
#'   score by which to stabilize the score, e.g. the predicted values from a
#'   regression model with no predictors.
#'
#' @return A vector of propensity score weights
#' @export
#'
#' @examples
#' propensity_scores <- c(.1, .3, .4, .3)
#' x <- c(0, 0, 1, 0)
#' wt_ate(propensity_scores, .exposure = x)
wt_ate <- function(.propensity, ...) {
  UseMethod("wt_ate")
}

#' @export
#' @rdname wt_ate
wt_ate.numeric <- function(.propensity, .exposure, .sigma = NULL, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    wts <- ate_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, stabilize = stabilize, stabilization_score = stabilization_score, ...)
  } else if (exposure_type == "continuous") {
    wts <- ate_continuous(.propensity = .propensity, .exposure = .exposure, .sigma = .sigma, .treated = .treated, .untreated = .untreated, stabilize = stabilize, stabilization_score = stabilization_score, ...)
  } else {
    abort_unsupported(exposure_type, "ATE")
  }

  psw(wts, "ate")
}

ate_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  wt <- (.exposure / .propensity) + ((1 - .exposure) / (1 - .propensity))


  if (isTRUE(stabilize) && is.null(stabilization_score)) {
    wt <- wt * mean(.exposure, na.rm = TRUE)
  } else if (isTRUE(stabilize) && !is.null(stabilization_score)) {
    wt <- wt * stabilization_score
  }

  wt
}

ate_continuous <- function(.propensity, .exposure, .sigma, stabilize = FALSE, stabilization_score = NULL, ...) {
  .propensity <- stats::dnorm(.exposure, mean = .propensity, sd = mean(.sigma))
  wt <- 1 / .propensity

  if (isTRUE(stabilize) && is.null(stabilization_score)) {
    wt <- wt * stats::dnorm(
      .exposure,
      mean(.exposure, na.rm = TRUE),
      stats::sd(.exposure, na.rm = TRUE)
    )
  } else if (isTRUE(stabilize) && !is.null(stabilization_score)) {
    wt <- wt * stabilization_score
  }

  wt
}

#' @export
#' @rdname wt_ate
wt_att <- function(.propensity, ...) {
  UseMethod("wt_att")
}

#' @export
#' @rdname wt_ate
wt_att.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    wts <- att_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATT")
  }

  psw(wts, "att")
}

att_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  ((.propensity * .exposure) / .propensity) +
    ((.propensity * (1 - .exposure)) / (1 - .propensity))
}

#' @export
#' @rdname wt_ate
wt_atu <- function(.propensity, ...) {
  UseMethod("wt_atu")
}

#' @export
#' @rdname wt_ate
wt_atu.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    wts <- atu_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATU")
  }

  psw(wts, "atu")
}

atu_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  wt <- (((1 - .propensity) * .exposure) / .propensity) +
    (((1 - .propensity) * (1 - .exposure)) / (1 - .propensity))

  wt
}

#' @export
#' @rdname wt_ate
wt_atm <- function(.propensity, ...) {
  UseMethod("wt_atm")
}

#' @export
#' @rdname wt_ate
wt_atm.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    wts <- atm_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATM")
  }

  psw(wts, "atm")
}

atm_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  pmin(.propensity, 1 - .propensity) /
    (.exposure * .propensity + (1 - .exposure) * (1 - .propensity))
}


#' @export
#' @rdname wt_ate
wt_ato <- function(.propensity, ...) {
  UseMethod("wt_ato")
}

#' @export
#' @rdname wt_ate
wt_ato.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    wts <- ato_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATO")
  }

  psw(wts, "ato")
}


ato_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

# --------------------------------------------------------------------
#  Add methods for ps_trim
#  (All 5 estimands: ate, att, atu, atm, ato)
# --------------------------------------------------------------------

#' @export
wt_ate.ps_trim <- function(.propensity, .exposure, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity) # includes NA for trimmed rows
  base_wt <- wt_ate.numeric(numeric_ps, .exposure, ...)

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_att.ps_trim <- function(.propensity, .exposure, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_att.numeric(numeric_ps, .exposure, ...)

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_atu.ps_trim <- function(.propensity, .exposure, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atu.numeric(numeric_ps, .exposure, ...)

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_atm.ps_trim <- function(.propensity, .exposure, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atm.numeric(numeric_ps, .exposure, ...)

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_ato.ps_trim <- function(.propensity, .exposure, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ato.numeric(numeric_ps, .exposure, ...)

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_ate.ps_trunc <- function(.propensity, .exposure, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ate.numeric(numeric_ps, .exposure, ...)

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_att.ps_trunc <- function(.propensity, .exposure, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_att.numeric(numeric_ps, .exposure, ...)

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_atu.ps_trunc <- function(.propensity, .exposure, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atu.numeric(numeric_ps, .exposure, ...)

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_atm.ps_trunc <- function(.propensity, .exposure, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atm.numeric(numeric_ps, .exposure, ...)

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_ato.ps_trunc <- function(.propensity, .exposure, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ato.numeric(numeric_ps, .exposure, ...)

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}
