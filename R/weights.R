#' Calculate Propensity Score Weights for Causal Inference
#'
#' @description This family of functions computes propensity score weights for
#'   various causal estimands:
#'
#' - **ATE** (Average Treatment Effect)
#' - **ATT** (Average Treatment Effect on the Treated)
#' - **ATU** (Average Treatment Effect on the Untreated, sometimes called
#'   the **ATC**, where the "C" stands for "control")
#' - **ATM** (Average Treatment Effect for the Evenly Matchable)
#' - **ATO** (Average Treatment Effect for the Overlap population)
#'
#'   The propensity score can be provided as a numeric vector of predicted
#'   probabilities or as a `data.frame` where each column represents the
#'   predicted probability for a level of the exposure. They can also be
#'   propensity score objects created by [ps_trim()], [ps_refit()], or
#'   [ps_trunc()]
#'
#'   The returned weights are encapsulated in a `psw` object, which is a numeric
#'   vector with additional attributes that record the estimand, and whether the
#'   weights have been stabilized, trimmed, or truncated.
#'
#' @details
#' ## Exposure Types
#'
#'   The functions support different types of exposures:
#'
#' - **`binary`**: For dichotomous treatments (e.g. 0/1).
#' - **`continuous`**: For numeric exposures. Here, weights are calculated via the normal density using
#'   `dnorm()`.
#' - **`categorical`**: Currently not supported (an error will be raised).
#' - **`auto`**: Automatically detects the exposure type based on `.exposure`.
#'
#'   ## Stabilization
#'
#'   For ATE weights, stabilization can improve the performance of the estimator
#'   by reducing variance. When `stabilize` is `TRUE` and no
#'   `stabilization_score` is provided, the weights are multiplied by the mean
#'   of `.exposure`. Alternatively, if a `stabilization_score` is provided, it
#'   is used as the multiplier.
#'
#'   ## Trimmed and Truncated Weights
#'
#'   In addition to the standard weight functions, versions exist for trimmed
#'   and truncated propensity score weights created by [ps_trim()],
#'   [ps_trunc()], and [ps_refit()]. These variants calculate the weights using
#'   modified propensity scores (trimmed or truncated) and update the estimand
#'   attribute accordingly.
#'
#'   The main functions (`wt_ate`, `wt_att`, `wt_atu`, `wt_atm`, and `wt_ato`)
#'   dispatch on the class of `.propensity`. For binary exposures, the weights
#'   are computed using inverse probability formulas. For continuous exposures
#'   (supported only for ATE), weights are computed as the inverse of the
#'   density function evaluated at the observed exposure.
#'
#' @param .propensity Either a numeric vector of predicted probabilities or a
#'   `data.frame` where each column corresponds to a level of the exposure.
#' @param .exposure The exposure variable. For binary exposures, a vector of 0s
#'   and 1s; for continuous exposures, a numeric vector.
#' @param exposure_type Character string specifying the type of exposure.
#'   Options are `"auto"`, `"binary"`, `"categorical"`, and `"continuous"`.
#'   Defaults to `"auto"`, which detects the type automatically.
#' @param .sigma For continuous exposures, a numeric vector of standard errors
#'   used with `dnorm()`. For example, this can be derived from the influence
#'   measures of a model (e.g., `influence(model)$sigma`).
#' @param .treated The value representing the treatment group. If not provided,
#'   it is automatically detected.
#' @param .untreated The value representing the control group. If not provided,
#'   it is automatically detected.
#' @param ... Reserved for future expansion. Not currently used.
#' @param stabilize Logical indicating whether to stabilize the weights. For ATE
#'   weights, stabilization multiplies the weight by either the mean of
#'   `.exposure` or the supplied `stabilization_score`.
#' @param stabilization_score Optional numeric value for stabilizing the weights
#'   (e.g., a predicted value from a regression model without predictors). Only
#'   used when `stabilize` is `TRUE`.
#'
#' @return A `psw` object (a numeric vector) with additional attributes:
#'   - **estimand**: A description of the estimand (e.g., "ate", "att").
#'   - **stabilized**: A logical flag indicating if stabilization was applied.
#'   - **trimmed**: A logical flag indicating if the weights are based on trimmed propensity scores.
#'   - **truncated**: A logical flag indicating if the weights are based on truncated propensity scores.
#'
#' @examples
#' ## ATE Weights with a Binary Exposure
#'
#' # Simulate a binary treatment and corresponding propensity scores
#' propensity_scores <- c(0.2, 0.7, 0.5, 0.8)
#' treatment <- c(0, 1, 0, 1)
#'
#' # Compute ATE weights (unstabilized)
#' weights_ate <- wt_ate(propensity_scores, .exposure = treatment)
#' weights_ate
#'
#' # Compute ATE weights with stabilization using the mean of the exposure
#' weights_ate_stab <- wt_ate(propensity_scores, .exposure = treatment, stabilize = TRUE)
#' weights_ate_stab
#'
#' ## ATT Weights for a Binary Exposure
#'
#' propensity_scores <- c(0.3, 0.6, 0.4, 0.7)
#' treatment <- c(1, 1, 0, 0)
#'
#' # Compute ATT weights
#' weights_att <- wt_att(propensity_scores, .exposure = treatment)
#' weights_att
#'
#' @seealso
#' - [psw()] for details on the structure of the returned weight objects.
#'
#' @export
wt_ate <- function(.propensity, .exposure, .sigma = NULL, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  UseMethod("wt_ate")
}

#' @export
wt_ate.numeric <- function(.propensity, .exposure, .sigma = NULL, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- ate_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated,
      stabilize = stabilize,
      stabilization_score = stabilization_score
    )
  } else if (exposure_type == "continuous") {
    wts <- ate_continuous(
      .propensity = .propensity,
      .exposure = .exposure,
      .sigma = .sigma,
      stabilize = stabilize,
      stabilization_score = stabilization_score
    )
  } else {
    abort_unsupported(exposure_type, "ATE")
  }

  psw(wts, "ate", stabilized = isTRUE(stabilize))
}

ate_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  if (isTRUE(stabilize) && is.null(stabilization_score)) {
    p1 <- mean(.exposure, na.rm = TRUE)
    p0 <- 1 - p1

    wt <- (.exposure * p1 / .propensity) +
      ((1 - .exposure) * p0 / (1 - .propensity))
  } else if (isTRUE(stabilize) && !is.null(stabilization_score)) {
    wt <- (.exposure / .propensity) + ((1 - .exposure) / (1 - .propensity))
    wt <- wt * stabilization_score
  } else {
    wt <- (.exposure / .propensity) + ((1 - .exposure) / (1 - .propensity))
  }

  wt
}

ate_continuous <- function(.propensity, .exposure, .sigma, stabilize = FALSE, stabilization_score = NULL) {
  # Compute population mean & variance of A
  un_mean  <- mean(.exposure, na.rm = TRUE)
  # sum(A−μ)^2 / n
  un_var   <- mean((.exposure - un_mean)^2, na.rm = TRUE)

  # compute population residual variance E[(A − E[A|X])^2]
  # sum(residual^2) / n
  cond_var <- mean((.exposure - .propensity)^2, na.rm = TRUE)

  # standardize into Z-scores
  z_num <- (.exposure - un_mean) / sqrt(un_var)
  z_den <- (.exposure - .propensity) / sqrt(cond_var)

  # evaluate standard normal densities on those Z's
  # f_A(A_i)
  f_num <- stats::dnorm(z_num)
  # f_{A|X}(A_i | X_i)
  f_den <- stats::dnorm(z_den)

  # build base weight = 1 / f_{A|X}
  wt <- 1 / f_den

  if (isTRUE(stabilize) && is.null(stabilization_score)) {
    wt <- wt * f_num
  } else if (isTRUE(stabilize) && !is.null(stabilization_score)) {
    wt <- wt * stabilization_score
  } else {
    alert_info(
      "Using unstabilized weights for continuous exposures is not recommended."
    )
  }

  wt
}


#' @export
#' @rdname wt_ate
wt_att <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  UseMethod("wt_att")
}

#' @export
wt_att.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(exposure_type, .exposure)

  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- att_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "ATT")
  }

  psw(wts, "att")
}

att_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL) {
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
wt_atu <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  UseMethod("wt_atu")
}

#' @export
wt_atu.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- atu_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "ATU")
  }

  psw(wts, "atu")
}

atu_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL) {
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
wt_atm <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  UseMethod("wt_atm")
}

#' @export
wt_atm.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- atm_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "ATM")
  }

  psw(wts, "atm")
}

atm_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL) {
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
wt_ato <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  UseMethod("wt_ato")
}

#' @export
wt_ato.numeric <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- ato_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "ATO")
  }

  psw(wts, "ato")
}


ato_binary <- function(.propensity, .exposure, .treated = NULL, .untreated = NULL) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

# --------------------------------------------------------------------
#  methods for `ps_trim()`
# --------------------------------------------------------------------

#' @export
wt_ate.ps_trim <- function(.propensity, .exposure, .sigma = NULL, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ate.numeric(
    numeric_ps,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_att.ps_trim <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_att.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_atu.ps_trim <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atu.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_atm.ps_trim <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atm.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_ato.ps_trim <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ato.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_ate.ps_trunc <- function(.propensity, .exposure, .sigma = NULL, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, stabilize = FALSE, stabilization_score = NULL, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ate.numeric(
    numeric_ps,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_att.ps_trunc <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_att.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_atu.ps_trunc <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atu.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_atm.ps_trunc <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_atm.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}

#' @export
wt_ato.ps_trunc <- function(.propensity, .exposure, exposure_type = c("auto", "binary", "categorical", "continuous"), .treated = NULL, .untreated = NULL, ...) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_ato.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}
