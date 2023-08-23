# propensity <- function(.propensity, .exposure, estimand = c("ate", "att", "atu", "atm", "ato"), exposure_type = c("auto", "binary", "categorical", "continuous")) {
#
# }
#
# add_weights <- function(.model, .exposure, estimand = c("ate", "att", "atu", "atm", "ato"), exposure_type = c("auto", "binary", "categorical", "continuous")) {
#   # add via broom
#   # check installed
# }

#' Calculate propensity score weights
#'
#' @param .propensity Either a vector of the predicted value of `.exposure` or a
#'   `data.frame` where each column is the predicted probability of a level of
#'   `.exposure`.
#' @param .exposure The exposure for which `.propensity` is calculated.
#' @param exposure_type The type of exposure. By default, automatically detected
#'   based on `.exposure`.
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
wt_ate.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    ate_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, stabilize = stabilize, stabilization_score = stabilization_score, ...)
  } else {
    abort_unsupported(exposure_type, "ATE")
  }
}

ate_binary <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
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
    att_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATT")
  }
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
    atu_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATU")
  }
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
    atm_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATM")
  }
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
    ato_binary(.propensity = .propensity, .exposure = .exposure, .treated = .treated, .untreated = .untreated, ...)
  } else {
    abort_unsupported(exposure_type, "ATO")
  }
}


ato_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, ...) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

