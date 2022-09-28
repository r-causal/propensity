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
#' @param .propensity Either a vector of the predicted value of `.exposure` or a `data.frame` where each column is the predicted probability of a level of `.exposure`.
#' @param .exposure The exposure for which `.propensity` is calculated.
#' @param exposure_type The type of exposure. By default, automatically detected based on `.exposure`.
#' @param .treated The treatment level of the exposure. Automatically detected by default.
#' @param .untreated The control level of the exposure. Automatically detected by default.
#' @param ... Passed to other functions Not currently used.
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
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

transform_exposure_binary <- function(.exposure, .treated = NULL, .untreated = NULL) {

   if (is_binary(.exposure)) {
    return(.exposure)
   }

  if (is.logical(.exposure)) {
    return(as.numeric(.exposure))
  }

  if (!is.null(.treated)) {
    return(ifelse(.exposure == .treated, 1, 0))
  }

  if (!is.null(.untreated)) {
    return(ifelse(.exposure != .untreated, 1, 0))
  }

  if (is.null(.treated) && is.null(.untreated) && has_two_levels(.exposure)) {
    levels <- if (is.factor(.exposure)) levels(.exposure) else sort(unique(.exposure))
    alert_info("Setting treatment to {.var {levels[[2]]}}")
    return(ifelse(.exposure == levels[[2]], 1, 0))
  } else {
    rlang::abort(c(
      "Don't know how to transform `.exposure` to 0/1 binary variable.",
      i = "Specify `.treated` and `.untreated.`"
    ))
  }
}

is_binary <- function(.exposure) {
  identical(sort(unique(.exposure)), c(0, 1))
}

has_two_levels <- function(.x) {
  length(unique(.x)) == 2
}
