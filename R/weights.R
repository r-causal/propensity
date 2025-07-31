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
#' - **Entropy** (Average Treatment Effect for the Entropy-weighted population)
#'
#'   The propensity score can be provided as a numeric vector of predicted
#'   probabilities, as a `data.frame` where each column represents the
#'   predicted probability for a level of the exposure, or as a fitted
#'   GLM object. They can also be propensity score objects created by
#'   [ps_trim()], [ps_refit()], or [ps_trunc()]
#'
#'   The returned weights are encapsulated in a `psw` object, which is a numeric
#'   vector with additional attributes that record the estimand, and whether the
#'   weights have been stabilized, trimmed, or truncated.
#'
#' @details
#' ## Theoretical Background
#'
#' Propensity score weighting is a method for estimating causal effects by
#' creating a pseudo-population where the exposure is independent of measured
#' confounders. The propensity score, \eqn{e(X)}, is the probability of receiving
#' treatment given observed covariates \eqn{X}. By weighting observations inversely
#' proportional to their propensity scores, we can balance the distribution of
#' covariates between treatment groups. Other weights allow for different target populations.
#'
#' ## Mathematical Formulas
#'
#' ### Binary Exposures
#'
#' For binary treatments (\eqn{A = 0} or \eqn{1}), the weights are:
#'
#' - **ATE**: \eqn{w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}}
#' - **ATT**: \eqn{w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}}
#' - **ATU**: \eqn{w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)}
#' - **ATM**: \eqn{w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}
#' - **ATO**: \eqn{w = A \cdot (1-e(X)) + (1-A) \cdot e(X)}
#' - **Entropy**: \eqn{w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}, where \eqn{h(e) = -[e \cdot \log(e) + (1-e) \cdot \log(1-e)]}
#'
#' ### Continuous Exposures
#'
#' For continuous treatments, weights use the density ratio:
#' \eqn{w = \frac{f_A(A)}{f_{A|X}(A|X)}}, where \eqn{f_A} is the marginal density of \eqn{A}
#' and \eqn{f_{A|X}} is the conditional density given \eqn{X}.
#'
#' ## Exposure Types
#'
#' The functions support different types of exposures:
#'
#' - **`binary`**: For dichotomous treatments (e.g. 0/1).
#' - **`continuous`**: For numeric exposures. Here, weights are calculated via the normal density using
#'   `dnorm()`.
#' - **`categorical`**: Currently not supported (an error will be raised).
#' - **`auto`**: Automatically detects the exposure type based on `.exposure`.
#'
#' ## Stabilization
#'
#' For ATE weights, stabilization can improve the performance of the estimator
#' by reducing variance. When `stabilize` is `TRUE` and no
#' `stabilization_score` is provided, the weights are multiplied by the mean
#' of `.exposure`. Alternatively, if a `stabilization_score` is provided, it
#' is used as the multiplier. Stabilized weights have the form:
#' \eqn{w_s = f_A(A) \times w}, where \eqn{f_A(A)} is the marginal probability or density.
#'
#' ## Weight Properties and Diagnostics
#'
#' Extreme weights can indicate:
#' - Positivity violations (near 0 or 1 propensity scores)
#' - Poor model specification
#' - Lack of overlap between treatment groups
#'
#' See the halfmoon package for tools to diagnose and visualize weights.
#'
#' You can address extreme weights in several ways. The first is to modify the target population:
#' use trimming, truncation, or alternative estimands (ATM, ATO, entropy).
#' Another technique that can help is stabilization, which reduces variance of the weights.
#'
#' ## Trimmed and Truncated Weights
#'
#' In addition to the standard weight functions, versions exist for trimmed
#' and truncated propensity score weights created by [ps_trim()],
#' [ps_trunc()], and [ps_refit()]. These variants calculate the weights using
#' modified propensity scores (trimmed or truncated) and update the estimand
#' attribute accordingly.
#'
#' @param .propensity Either a numeric vector of predicted probabilities, a
#'   `data.frame` where each column corresponds to a level of the exposure,
#'   or a fitted GLM object. For data frames, the second column is used by
#'   default for binary exposures unless specified otherwise with
#'   `.propensity_col`. For GLM objects, fitted values are extracted
#'   automatically.
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
#' @param .propensity_col When `.propensity` is a data frame, specifies which
#'   column to use for propensity scores. Can be a column name (quoted or
#'   unquoted) or a numeric index. Defaults to the second column if available,
#'   otherwise the first.
#'
#' @return A `psw` object (a numeric vector) with additional attributes:
#'   - **estimand**: A description of the estimand (e.g., "ate", "att").
#'   - **stabilized**: A logical flag indicating if stabilization was applied.
#'   - **trimmed**: A logical flag indicating if the weights are based on trimmed propensity scores.
#'   - **truncated**: A logical flag indicating if the weights are based on truncated propensity scores.
#'
#' @examples
#' ## Basic Usage with Binary Exposures
#'
#' # Simulate a simple dataset
#' set.seed(123)
#' n <- 100
#' propensity_scores <- runif(n, 0.1, 0.9)
#' treatment <- rbinom(n, 1, propensity_scores)
#'
#' # Calculate different weight types
#' weights_ate <- wt_ate(propensity_scores, treatment)
#' weights_att <- wt_att(propensity_scores, treatment)
#' weights_atu <- wt_atu(propensity_scores, treatment)
#' weights_atm <- wt_atm(propensity_scores, treatment)
#' weights_ato <- wt_ato(propensity_scores, treatment)
#' weights_entropy <- wt_entropy(propensity_scores, treatment)
#'
#' # Compare weight distributions
#' summary(weights_ate)
#' summary(weights_ato)  # Often more stable than ATE
#'
#' ## Stabilized Weights
#'
#' # Stabilization reduces variance
#' weights_ate_stab <- wt_ate(propensity_scores, treatment, stabilize = TRUE)
#'
#' # Compare coefficient of variation
#' sd(weights_ate) / mean(weights_ate)      # Unstabilized
#' sd(weights_ate_stab) / mean(weights_ate_stab)  # Stabilized (lower is better)
#'
#' ## Handling Extreme Propensity Scores
#'
#' # Create data with positivity violations
#' ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
#' trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)
#'
#' # Standard ATE weights can be extreme
#' wt_extreme <- wt_ate(ps_extreme, trt_extreme)
#' # Very large!
#' max(wt_extreme)
#'
#' # ATO weights are bounded
#' wt_extreme_atm <- wt_ato(ps_extreme, trt_extreme)
#' # Much more reasonable
#' max(wt_extreme_atm)
#' # but they target a different population
#' estimand(wt_extreme_atm) # "ato"
#'
#' ## Working with Data Frames
#'
#' # Example with custom data frame
#' ps_df <- data.frame(
#'   control = c(0.9, 0.7, 0.3, 0.1),
#'   treated = c(0.1, 0.3, 0.7, 0.9)
#' )
#' exposure <- c(0, 0, 1, 1)
#'
#' # Uses second column by default (treated probabilities)
#' wt_ate(ps_df, exposure)
#'
#' # Explicitly specify column by name
#' wt_ate(ps_df, exposure, .propensity_col = "treated")
#'
#' # Or by position
#' wt_ate(ps_df, exposure, .propensity_col = 2)
#'
#' ## Working with GLM Objects
#'
#' # Fit a propensity score model
#' set.seed(123)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' treatment <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))
#'
#' ps_model <- glm(treatment ~ x1 + x2, family = binomial)
#'
#' # Use GLM directly for weight calculation
#' weights_from_glm <- wt_ate(ps_model, treatment)
#'
#' @references
#'
#' For detailed guidance on causal inference in R, see [*Causal Inference in R*](https://www.r-causal.org/)
#' by Malcolm Barrett, Lucy D'Agostino McGowan, and Travis Gerke.
#'
#' ## Foundational Papers
#'
#' Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity
#' score in observational studies for causal effects. *Biometrika*, 70(1), 41-55.
#'
#' ## Estimand-Specific Methods
#'
#' Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
#' propensity score analysis. *The International Journal of Biostatistics*, 9(2),
#' 215-234. (ATM weights)
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. *Journal of the American Statistical Association*,
#' 113(521), 390-400. (ATO weights)
#'
#' Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score weighting
#' under limited overlap and model misspecification. *Statistical Methods in
#' Medical Research*, 29(12), 3721-3756. (Entropy weights)
#'
#' ## Continuous Exposures
#'
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. *Applied Bayesian Modeling and Causal Inference from
#' Incomplete-Data Perspectives*, 226164, 73-84.
#'
#' ## Practical Guidance
#'
#' Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
#' using inverse probability of treatment weighting (IPTW) using the propensity
#' score to estimate causal treatment effects in observational studies.
#' *Statistics in Medicine*, 34(28), 3661-3679.
#'
#' @seealso
#' - [psw()] for details on the structure of the returned weight objects.
#' - [ps_trim()], [ps_trunc()], and [ps_refit()] for handling extreme weights.
#' - [ps_calibrate()] for calibrating weights.
#'
#' @export
wt_ate <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  UseMethod("wt_ate")
}

#' @export
wt_ate.default <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_ate.numeric <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
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

#' @export
#' @rdname wt_ate
wt_ate.data.frame <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_ate.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

#' @export
wt_ate.glm <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # For continuous exposures, extract sigma if not provided
  if (
    is.null(.sigma) &&
      match_exposure_type(exposure_type, .exposure) == "continuous"
  ) {
    .sigma <- stats::influence(.propensity)$sigma
  }

  # Call the numeric method
  wt_ate.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
}

# Helper function to extract propensity scores from data frames
extract_propensity_from_df <- function(
  .propensity,
  .propensity_col_quo = NULL
) {
  if (!rlang::quo_is_null(.propensity_col_quo)) {
    col_pos <- tryCatch(
      tidyselect::eval_select(
        .propensity_col_quo,
        data = .propensity
      ),
      error = function(e) {
        abort(
          paste0("Column selection failed: ", e$message),
          error_class = "propensity_df_column_error"
        )
      }
    )

    if (length(col_pos) != 1) {
      abort(
        "`.propensity_col` must select exactly one column.",
        error_class = "propensity_df_column_error"
      )
    }

    ps_vec <- .propensity[[col_pos]]
  } else {
    # Default behavior: use second column if available, otherwise first
    if (ncol(.propensity) >= 2) {
      ps_vec <- .propensity[[2]]
    } else if (ncol(.propensity) == 1) {
      ps_vec <- .propensity[[1]]
    } else {
      abort(
        "`.propensity` data frame must have at least one column.",
        error_class = "propensity_df_ncol_error"
      )
    }
  }

  ps_vec
}

# Helper function to extract propensity scores from GLM objects
extract_propensity_from_glm <- function(.propensity) {
  # Check if it's a valid GLM object
  if (!inherits(.propensity, "glm")) {
    abort(
      "`.propensity` must be a GLM object.",
      error_class = "propensity_glm_type_error"
    )
  }

  # Check if it's a binomial GLM for binary propensity scores
  if (
    !is.null(.propensity$family) &&
      .propensity$family$family == "binomial"
  ) {
    # Get predicted probabilities
    ps_vec <- stats::predict(.propensity, type = "response")
  } else {
    # For non-binomial GLMs, get linear predictor
    ps_vec <- stats::fitted(.propensity)
  }

  ps_vec
}

ate_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL
) {
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_att")
}

#' @export
wt_att.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_att.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )

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

#' @export
#' @rdname wt_ate
wt_att.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_att.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_att.glm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_att.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

att_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
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
wt_atu <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_atu")
}

#' @export
wt_atu.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_atu.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )
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

#' @export
#' @rdname wt_ate
wt_atu.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_atu.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.glm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_atu.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

atu_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
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
wt_atm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_atm")
}

#' @export
wt_atm.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_atm.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )
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

#' @export
#' @rdname wt_ate
wt_atm.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_atm.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.glm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_atm.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

atm_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
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
wt_ato <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_ato")
}

#' @export
wt_ato.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_ato.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )
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

#' @export
#' @rdname wt_ate
wt_ato.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_ato.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.glm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_ato.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

ato_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  (1 - .propensity) * .exposure + .propensity * (1 - .exposure)
}

#' @export
#' @rdname wt_ate
wt_entropy <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_entropy")
}

#' @export
wt_entropy.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
wt_entropy.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical")
  )
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- entropy_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "entropy")
  }

  psw(wts, "entropy")
}

#' @export
#' @rdname wt_ate
wt_entropy.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  # Capture the column selection expression
  col_quo <- rlang::enquo(.propensity_col)

  # Extract propensity scores from data frame
  ps_vec <- extract_propensity_from_df(.propensity, col_quo)

  # Call the numeric method
  wt_entropy.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.glm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  # Extract fitted values (propensity scores) from GLM
  ps_vec <- extract_propensity_from_glm(.propensity)

  # Call the numeric method
  wt_entropy.numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

entropy_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
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
#  methods for `ps_trim()`
# --------------------------------------------------------------------

#' @export
wt_ate.ps_trim <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
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
wt_att.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_atu.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_atm.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_ato.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_ate.ps_trunc <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
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
wt_att.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_atu.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_atm.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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
wt_ato.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
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

#' @export
wt_entropy.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_entropy.numeric(
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
wt_entropy.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_entropy.numeric(
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
