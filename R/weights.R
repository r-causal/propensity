#' Calculate propensity score weights
#'
#' @description
#' Compute inverse probability weights for causal inference under different
#' estimands. Each function targets a different population:
#'
#' - `wt_ate()`: **Average Treatment Effect** -- the full population.
#' - `wt_att()`: **Average Treatment Effect on the Treated** -- the treated
#'   (focal) group.
#' - `wt_atu()`: **Average Treatment Effect on the Untreated** -- the
#'   untreated (reference) group. `wt_atc()` is an alias.
#' - `wt_atm()`: **Average Treatment Effect for the Evenly Matchable** --
#'   units with the most overlap.
#' - `wt_ato()`: **Average Treatment Effect for the Overlap Population** --
#'   weights proportional to overlap.
#' - `wt_entropy()`: **Entropy-weighted Average Treatment Effect** --
#'   an entropy-balanced population.
#' - `wt_cens()`: **Inverse probability of censoring weights** -- uses the
#'   same formula as `wt_ate()` but labels the estimand `"uncensored"`. Use
#'   these to adjust for censoring in survival analysis, not for treatment
#'   weighting.
#'
#' `.propensity` accepts a numeric vector of predicted probabilities, a
#' `data.frame` of per-level probabilities, a fitted `glm` object, or a
#' modified propensity score created by [ps_trim()], [ps_trunc()],
#' [ps_refit()], or [ps_calibrate()].
#'
#' All functions return a [`psw`] object -- a numeric vector that tracks the
#' estimand, stabilization status, and any trimming or truncation applied.
#'
#' @details
#' ## Exposure types
#'
#' All weight functions support binary exposures. `wt_ate()` and `wt_cens()`
#' also support continuous exposures. All except `wt_cens()` support
#' categorical exposures.
#'
#' - **Binary**: `.exposure` is a two-level vector (e.g., 0/1, logical, or a
#'   two-level factor). `.propensity` is a numeric vector of P(treatment | X).
#' - **Categorical**: `.exposure` is a factor or character vector with 3+
#'   levels. `.propensity` must be a matrix or data frame with one column per
#'   level, where rows sum to 1.
#' - **Continuous**: `.exposure` is a numeric vector. `.propensity` is a
#'   vector of conditional means (fitted values). Weights use a normal density
#'   ratio; stabilization is strongly recommended.
#' - **Auto** (default): Detects the exposure type from `.exposure`.
#'
#' ## Stabilization
#'
#' Setting `stabilize = TRUE` multiplies the base weight by an estimate of
#' P(A) (binary) or f_A(A) (continuous), reducing variance. When no
#' `stabilization_score` is supplied, the marginal mean of `.exposure` is
#' used. Stabilization is supported for ATE and censoring weights
#' (`wt_ate()` and `wt_cens()`) and is strongly recommended for continuous
#' exposures.
#'
#' ## Handling extreme weights
#'
#' Extreme weights signal positivity violations, poor model fit, or limited
#' overlap. You can address them by:
#'
#' - Choosing an overlap-focused estimand (`wt_ato()`, `wt_atm()`,
#'   `wt_entropy()`), which down-weight units in regions of poor overlap.
#' - Trimming ([ps_trim()]) or truncating ([ps_trunc()]) propensity scores
#'   before computing weights.
#' - Calibrating weights with [ps_calibrate()].
#' - Stabilizing ATE weights (`stabilize = TRUE`).
#'
#' See the [halfmoon](https://CRAN.R-project.org/package=halfmoon) package
#' for weight diagnostics and visualization.
#'
#' ## Weight formulas
#'
#' ### Binary exposures
#'
#' For binary treatments (\eqn{A \in \{0, 1\}}), with propensity score
#' \eqn{e(X) = P(A=1 \mid X)}:
#'
#' - **ATE**: \eqn{w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}}
#' - **ATT**: \eqn{w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}}
#' - **ATU**: \eqn{w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)}
#' - **ATM**: \eqn{w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}
#' - **ATO**: \eqn{w = A \cdot (1-e(X)) + (1-A) \cdot e(X)}
#' - **Entropy**: \eqn{w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}, where \eqn{h(e) = -[e \log(e) + (1-e) \log(1-e)]}
#'
#' ### Continuous exposures
#'
#' Weights use the density ratio
#' \eqn{w = f_A(A) / f_{A|X}(A \mid X)}, where \eqn{f_A} is the marginal
#' density and \eqn{f_{A|X}} is the conditional density (both assumed
#' normal). Only `wt_ate()` and `wt_cens()` support continuous exposures.
#'
#' ### Categorical exposures
#'
#' For \eqn{K}-level treatments, weights take the tilting-function form
#' \eqn{w_i = h(\mathbf{e}_i) / e_{i,Z_i}}, where \eqn{e_{i,Z_i}} is the
#' propensity for unit \eqn{i}'s observed level and \eqn{h(\cdot)} depends
#' on the estimand:
#'
#' - **ATE**: \eqn{h(\mathbf{e}) = 1}
#' - **ATT**: \eqn{h(\mathbf{e}) = e_{\text{focal}}}
#' - **ATU**: \eqn{h(\mathbf{e}) = 1 - e_{\text{focal}}}
#' - **ATM**: \eqn{h(\mathbf{e}) = \min(e_1, \ldots, e_K)}
#' - **ATO**: \eqn{h(\mathbf{e}) = \bigl(\sum_k 1/e_k\bigr)^{-1}}
#' - **Entropy**: \eqn{h(\mathbf{e}) = -\sum_k e_k \log(e_k)}
#'
#' @param .propensity Propensity scores in one of several forms:
#'   * A **numeric vector** of predicted probabilities (binary/continuous).
#'   * A **data frame** or matrix with one column per exposure level
#'     (categorical), or two columns for binary (see `.propensity_col`).
#'   * A fitted **`glm`** object -- fitted values are extracted automatically.
#'   * A modified propensity score created by [ps_trim()], [ps_trunc()],
#'     [ps_refit()], or [ps_calibrate()].
#' @param .exposure The exposure (treatment) variable. For binary exposures, a
#'   numeric 0/1 vector, logical, or two-level factor. For categorical
#'   exposures, a factor or character vector. For continuous exposures, a
#'   numeric vector. Optional when `.propensity` is a `glm` object (extracted
#'   from the model).
#' @param exposure_type Type of exposure: `"auto"` (default), `"binary"`,
#'   `"categorical"`, or `"continuous"`. `"auto"` detects the type from
#'   `.exposure`.
#' @param .sigma Numeric vector of observation-level standard deviations for
#'   continuous exposures (e.g., `influence(model)$sigma`). Extracted
#'   automatically when `.propensity` is a `glm` object.
#' @param .treated `r lifecycle::badge("deprecated")` Use `.focal_level` instead.
#' @param .untreated `r lifecycle::badge("deprecated")` Use `.reference_level` instead.
#' @param .focal_level The value of `.exposure` representing the focal
#'   (treated) group. For binary exposures, defaults to the higher value.
#'   Required for `wt_att()` and `wt_atu()` with categorical exposures.
#' @param .reference_level The value of `.exposure` representing the reference
#'   (control) group. Automatically detected if not supplied.
#' @param ... These dots are for future extensions and must be empty.
#' @param stabilize If `TRUE`, multiply weights by an estimate of the marginal
#'   treatment probability (binary) or density (continuous). Only supported by
#'   `wt_ate()` and `wt_cens()`. See **Stabilization** in Details.
#' @param stabilization_score Optional numeric value to use as the
#'   stabilization multiplier instead of the default (the marginal mean of
#'   `.exposure`). Ignored when `stabilize = FALSE`.
#' @param .propensity_col Column to use when `.propensity` is a data frame
#'   with a binary exposure. Accepts a column name (quoted or unquoted) or
#'   numeric index. Defaults to the second column. Ignored for categorical
#'   exposures, where all columns are used.
#'
#' @return A [`psw`] vector (a double vector with class `psw`) carrying
#'   these attributes:
#'   - `estimand`: character, e.g. `"ate"`, `"att"`, `"uncensored"`.
#'   - `stabilized`: logical, whether stabilization was applied.
#'   - `trimmed`: logical, whether the propensity scores were trimmed.
#'   - `truncated`: logical, whether the propensity scores were truncated.
#'   - `calibrated`: logical, whether the propensity scores were calibrated.
#'
#' @examples
#' # -- Binary exposure, numeric propensity scores ----------------------
#' set.seed(123)
#' ps <- runif(100, 0.1, 0.9)
#' trt <- rbinom(100, 1, ps)
#'
#' wt_ate(ps, trt)
#' wt_att(ps, trt)
#' wt_atu(ps, trt)
#' wt_atm(ps, trt)
#' wt_ato(ps, trt)
#' wt_entropy(ps, trt)
#'
#' # Stabilized ATE weights (reduces variance)
#' wt_ate(ps, trt, stabilize = TRUE)
#'
#' # Inspect the result
#' w <- wt_ate(ps, trt)
#' estimand(w)
#' summary(w)
#'
#' # -- Overlap-focused estimands handle extreme PS better --------------
#' ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
#' trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)
#'
#' max(wt_ate(ps_extreme, trt_extreme))
#' max(wt_ato(ps_extreme, trt_extreme))
#'
#' # -- From a fitted GLM -----------------------------------------------
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' trt2 <- rbinom(100, 1, plogis(0.5 * x1 + 0.3 * x2))
#' ps_model <- glm(trt2 ~ x1 + x2, family = binomial)
#'
#' # Exposure is extracted from the model automatically
#' wt_ate(ps_model)
#'
#' # -- Data frame input ------------------------------------------------
#' ps_df <- data.frame(
#'   control = c(0.9, 0.7, 0.3, 0.1),
#'   treated = c(0.1, 0.3, 0.7, 0.9)
#' )
#' exposure <- c(0, 0, 1, 1)
#' wt_ate(ps_df, exposure)
#' wt_ate(ps_df, exposure, .propensity_col = "treated")
#'
#' # -- Censoring weights -----------------------------------------------
#' cens_ps <- runif(50, 0.6, 0.95)
#' cens_ind <- rbinom(50, 1, cens_ps)
#' wt_cens(cens_ps, cens_ind)
#' estimand(wt_cens(cens_ps, cens_ind))  # "uncensored"
#'
#' @references
#' Barrett, M., D'Agostino McGowan, L., & Gerke, T. *Causal Inference in R*.
#' \url{https://www.r-causal.org/}
#'
#' Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity
#' score in observational studies for causal effects. *Biometrika*, 70(1),
#' 41--55.
#'
#' Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
#' propensity score analysis. *The International Journal of Biostatistics*,
#' 9(2), 215--234. (ATM weights)
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. *Journal of the American Statistical
#' Association*, 113(521), 390--400. (ATO weights)
#'
#' Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score weighting
#' under limited overlap and model misspecification. *Statistical Methods in
#' Medical Research*, 29(12), 3721--3756. (Entropy weights)
#'
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. In *Applied Bayesian Modeling and Causal Inference from
#' Incomplete-Data Perspectives* (pp. 73--84).
#'
#' Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
#' using inverse probability of treatment weighting (IPTW). *Statistics in
#' Medicine*, 34(28), 3661--3679.
#'
#' @seealso
#' * [psw()] for the returned weight vector class.
#' * [ps_trim()], [ps_trunc()], [ps_refit()], and [ps_calibrate()] for
#'   modifying propensity scores before weighting.
#' * [ipw()] for inverse-probability-weighted estimation of causal effects.
#'
#' @export
wt_ate <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_ate")
}

#' @export
wt_ate.default <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_ate.numeric <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_ate"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(exposure_type, .exposure)

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- ate_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      stabilize = stabilize,
      stabilization_score = stabilization_score
    )
  } else if (exposure_type == "continuous") {
    check_lengths_match(.propensity, .exposure)
    wts <- ate_continuous(
      .propensity = .propensity,
      .exposure = .exposure,
      .sigma = .sigma,
      stabilize = stabilize,
      stabilization_score = stabilization_score
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    # including the more specific matrix checks
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "ate",
      .focal_level = NULL,
      stabilize = stabilize,
      stabilization_score = stabilization_score,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "ATE")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "ate", stabilized = isTRUE(stabilize))

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_ate.data.frame <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_ate.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
    .propensity_col_quo = col_quo,
    .sigma = .sigma,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ate.glm <- function(
  .propensity,
  .exposure = NULL,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)
  exposure_type <- match_exposure_type(exposure_type, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # For continuous exposures, extract sigma if not provided
  if (is.null(.sigma) && exposure_type == "continuous") {
    .sigma <- stats::influence(.propensity)$sigma
  }

  # Call the numeric method
  wt_ate.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}


ate_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
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

ate_continuous <- function(
  .propensity,
  .exposure,
  .sigma,
  stabilize = FALSE,
  stabilization_score = NULL
) {
  # Compute population mean & variance of A
  un_mean <- mean(.exposure, na.rm = TRUE)
  # sum(A−μ)^2 / n
  un_var <- mean((.exposure - un_mean)^2, na.rm = TRUE)

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
wt_att <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_att")
}

#' @export
wt_att.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_att.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_att"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- att_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "att",
      .focal_level = .focal_level,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "ATT")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "att")

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_att.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_att.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_att.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_att.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

att_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  ((.propensity * .exposure) / .propensity) +
    ((.propensity * (1 - .exposure)) / (1 - .propensity))
}

#' @export
#' @rdname wt_ate
wt_atu <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_atu")
}

#' @export
wt_atu.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_atu.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_atu"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- atu_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atu",
      .focal_level = .focal_level,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "ATU")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "atu")

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_atu.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atu.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_atu.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

atu_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  wt <- (((1 - .propensity) * .exposure) / .propensity) +
    (((1 - .propensity) * (1 - .exposure)) / (1 - .propensity))

  wt
}

#' @export
#' @rdname wt_ate
wt_atm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_atm")
}

#' @export
wt_atm.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_atm.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_atm"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- atm_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atm",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "ATM")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "atm")

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_atm.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atm.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_atm.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

atm_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  pmin(.propensity, 1 - .propensity) /
    (.exposure * .propensity + (1 - .exposure) * (1 - .propensity))
}


#' @export
#' @rdname wt_ate
wt_ato <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_ato")
}

#' @export
wt_ato.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_ato.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_ato"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- ato_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "ato",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "ATO")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "ato")

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_ato.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_ato.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_ato.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

ato_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

#' @export
#' @rdname wt_ate
wt_entropy <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_entropy")
}

#' @export
wt_entropy.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_entropy.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  rlang::check_dots_empty()

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_entropy"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- entropy_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "entropy",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = rlang::caller_env()
    )
  } else {
    abort_unsupported(exposure_type, "entropy")
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "entropy")

  # Preserve categorical attributes if they exist
  preserve_categorical_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_entropy.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_entropy.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_entropy.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

entropy_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  # Entropy tilting function: h(e) = -[e*log(e) + (1-e)*log(1-e)]
  h_e <- -.propensity *
    log(.propensity) -
    (1 - .propensity) * log(1 - .propensity)

  # Calculate weights using equation approach
  # w = h(e)/e for treated (exposure=1), w = h(e)/(1-e) for control (exposure=0)
  weights <- h_e /
    (.exposure * .propensity + (1 - .exposure) * (1 - .propensity))

  weights
}

# --------------------------------------------------------------------
#  wt_atc() - alias for wt_atu()
# --------------------------------------------------------------------

#' @export
#' @rdname wt_ate
wt_atc <- wt_atu

# --------------------------------------------------------------------
#  wt_cens() - censoring weights
# --------------------------------------------------------------------

#' @export
#' @rdname wt_ate
wt_cens <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_cens")
}

#' @export
wt_cens.default <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_cens.numeric <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Get weights using ATE formula
  wts <- wt_ate.numeric(
    .propensity = .propensity,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )

  # Update estimand to "uncensored"
  estimand(wts) <- "uncensored"
  wts
}

#' @export
#' @rdname wt_ate
wt_cens.data.frame <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_cens.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
    .propensity_col_quo = col_quo,
    .sigma = .sigma,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_cens.glm <- function(
  .propensity,
  .exposure = NULL,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)
  exposure_type <- match_exposure_type(exposure_type, .exposure)

  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # For continuous exposures, extract sigma if not provided
  if (is.null(.sigma) && exposure_type == "continuous") {
    .sigma <- stats::influence(.propensity)$sigma
  }

  # Call the numeric method
  wt_cens.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

# --------------------------------------------------------------------
#  methods for `ps_trim()`
# --------------------------------------------------------------------

#' @export
wt_ate.ps_trim <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trim",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_att.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ate.ps_trunc <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trunc",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_att.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atc.ps_trim <- wt_atu.ps_trim

#' @export
wt_atc.ps_trunc <- wt_atu.ps_trunc

#' @export
wt_cens.ps_trim <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trim",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_cens.ps_trunc <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trunc",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

# --------------------------------------------------------------------
#  Categorical exposure weight calculations
# --------------------------------------------------------------------

calculate_categorical_weights <- function(
  ps_matrix,
  .exposure,
  estimand,
  .focal_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  call = rlang::caller_env()
) {
  # Ensure exposure is a factor
  .exposure <- transform_exposure_categorical(
    .exposure,
    .focal_level,
    call = call
  )

  # Validate propensity score matrix
  ps_matrix <- check_ps_matrix(ps_matrix, .exposure, call = call)

  # Get dimensions
  n <- length(.exposure)
  k <- nlevels(.exposure)
  levels_exp <- levels(.exposure)

  # Create indicator matrix for exposures
  # Each row i has a 1 in column j if observation i has exposure level j
  Z <- matrix(0, nrow = n, ncol = k)
  for (j in 1:k) {
    Z[.exposure == levels_exp[j], j] <- 1
  }

  # Extract the propensity score for each unit's actual exposure
  # e_{i,Z_i} = sum over j of Z_{ij} * e_{ij}
  e_actual <- rowSums(Z * ps_matrix)

  # Calculate tilting function h(e_i) based on estimand
  h_e <- switch(
    estimand,
    "ate" = rep(1, n), # h(e) = 1 for ATE
    "att" = {
      if (is.null(.focal_level)) {
        abort(
          "Focal category must be specified for ATT with categorical exposures.",
          error_class = "propensity_focal_required_error"
        )
      }
      # h(e) = e_focal
      focal_idx <- which(levels_exp == .focal_level)
      ps_matrix[, focal_idx]
    },
    "atu" = {
      if (is.null(.focal_level)) {
        abort(
          "Focal category must be specified for ATU with categorical exposures.",
          error_class = "propensity_focal_required_error"
        )
      }
      # For ATU, we want weights for all non-focal categories
      # h(e) = 1 - e_focal
      focal_idx <- which(levels_exp == .focal_level)
      1 - ps_matrix[, focal_idx]
    },
    "ato" = {
      # h(e) = 1 / sum(1/e_k) - harmonic mean denominator
      1 / rowSums(1 / ps_matrix)
    },
    "atm" = {
      # h(e) = min(e_1, ..., e_K)
      apply(ps_matrix, 1, min)
    },
    "entropy" = {
      # h(e) = -sum(e_k * log(e_k)) - entropy
      # Need to handle e_k = 0 case to avoid -Inf
      ps_safe <- ps_matrix
      ps_safe[ps_matrix == 0] <- .Machine$double.eps
      -rowSums(ps_safe * log(ps_safe))
    },
    abort(
      "Unknown estimand: {estimand}",
      error_class = "propensity_unknown_estimand_error"
    )
  )

  # Calculate weights: w_i = h(e_i) / e_{i,Z_i}
  weights <- h_e / e_actual

  # Apply stabilization only for ATE
  if (isTRUE(stabilize)) {
    if (estimand != "ate") {
      abort(
        "Stabilization is only supported for ATE with categorical exposures.",
        error_class = "propensity_stabilize_categorical_error"
      )
    }

    if (!is.null(stabilization_score)) {
      weights <- weights * stabilization_score
    } else {
      # For categorical, use marginal probabilities
      p_marginal <- table(.exposure) / n

      # Create stabilization weights based on marginal probabilities
      stab_wts <- numeric(n)
      for (j in 1:k) {
        stab_wts[.exposure == levels_exp[j]] <- p_marginal[j]
      }

      weights <- weights * stab_wts
    }
  }

  # Add attributes for categorical weights
  attr(weights, "n_categories") <- k
  attr(weights, "category_names") <- levels_exp
  if (!is.null(.focal_level)) {
    attr(weights, "focal_category") <- .focal_level
  }

  weights
}

# ps_calib methods ----

#' @export
wt_ate.ps_calib <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "calib",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_att.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atc.ps_calib <- wt_atu.ps_calib

#' @export
wt_cens.ps_calib <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "calib",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}
