#' Calculate Propensity Score Weights for Causal Inference
#'
#' @description This family of functions computes propensity score weights for
#'   various causal estimands:
#'
#' - **ATE** (Average Treatment Effect)
#' - **ATT** (Average Treatment Effect on the Treated)
#' - **ATU** (Average Treatment Effect on the Untreated, sometimes called
#'   the **ATC**, where the "C" stands for "control"). `wt_atc()` is provided
#'   as an alias for `wt_atu()`
#' - **ATM** (Average Treatment Effect for the Evenly Matchable)
#' - **ATO** (Average Treatment Effect for the Overlap population)
#' - **Entropy** (Average Treatment Effect for the Entropy-weighted population)
#' - **Censoring weights** can be calculated using `wt_cens()`, which uses
#'   the same formula as ATE weights but with estimand "uncensored". These
#'   are useful for handling censoring in survival analysis
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
#' ### Categorical Exposures
#'
#' For categorical treatments with \eqn{K} levels, weights use a tilting function approach:
#' \eqn{w_i = \frac{h(e_i)}{e_{i,Z_i}}}, where \eqn{e_{i,Z_i}} is the propensity score for unit \eqn{i}'s
#' observed treatment level, and \eqn{h(e_i)} is a tilting function that depends on the estimand:
#'
#' - **ATE**: \eqn{h(e) = 1}
#' - **ATT**: \eqn{h(e) = e_{focal}} (propensity score for the focal category)
#' - **ATU**: \eqn{h(e) = 1 - e_{focal}} (complement of focal category propensity)
#' - **ATM**: \eqn{h(e) = \min(e_1, ..., e_K)}
#' - **ATO**: \eqn{h(e) = 1 / \sum_k(1/e_k)} (reciprocal of harmonic mean denominator)
#' - **Entropy**: \eqn{h(e) = -\sum_k[e_k \cdot \log(e_k)]} (entropy of propensity scores)
#'
#' ## Exposure Types
#'
#' The functions support different types of exposures:
#'
#' - **`binary`**: For dichotomous treatments (e.g. 0/1).
#' - **`continuous`**: For numeric exposures. Here, weights are calculated via the normal density using
#'   `dnorm()`.
#' - **`categorical`**: For exposures with more than 2 categories. Requires `.propensity` to be a
#'   matrix or data frame with columns representing propensity scores for each category.
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
#'   and 1s; for continuous exposures, a numeric vector. When `.propensity` is
#'   a GLM object, this argument is optional and will be extracted from the
#'   model if not provided.
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
#' @param focal For categorical exposures with ATT or ATU estimands, specifies
#'   the focal category. Must be one of the levels of the exposure variable.
#'   Required for `wt_att()` and `wt_atu()` with categorical exposures.
#' @param ... Reserved for future expansion. Not currently used.
#' @param stabilize Logical indicating whether to stabilize the weights. For ATE
#'   weights, stabilization multiplies the weight by either the mean of
#'   `.exposure` or the supplied `stabilization_score`. Note: stabilization is only
#'   supported for ATE and continuous exposures.
#' @param stabilization_score Optional numeric value for stabilizing the weights
#'   (e.g., a predicted value from a regression model without predictors). Only
#'   used when `stabilize` is `TRUE`.
#' @param .propensity_col With a binary exposure, when `.propensity` is a data frame, specifies which
#'   column to use for propensity scores. Can be a column name (quoted or
#'   unquoted) or a numeric index. Defaults to the second column if available,
#'   otherwise the first. For categorical exposures, the entire data frame is
#'   used as a matrix of propensity scores.
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
#' wt_extreme_ato <- wt_ato(ps_extreme, trt_extreme)
#' # Much more reasonable
#' max(wt_extreme_ato)
#' # but they target a different population
#' estimand(wt_extreme_ato) # "ato"
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
#' # Or omit the exposure argument (it will be extracted from the GLM)
#' weights_from_glm_auto <- wt_ate(ps_model)
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
  abort_no_method(.propensity)
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
    check_lengths_match(.propensity, .exposure)
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
      focal = NULL,
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL
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
  .exposure = NULL,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)

  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

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
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    ...
  )
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
  ...,
  focal = NULL
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
  ...,
  focal = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_att.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  focal = NULL
) {
  rlang::check_dots_empty()
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
      .treated = .treated,
      .untreated = .untreated
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "att",
      focal = focal,
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
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL,
  focal = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_att.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .treated = .treated,
    .untreated = .untreated,
    focal = focal,
    ...
  )
}

#' @export
wt_att.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  focal = NULL
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
    .treated = .treated,
    .untreated = .untreated,
    focal = focal,
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
  ...,
  focal = NULL
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
  ...,
  focal = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_atu.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  focal = NULL
) {
  rlang::check_dots_empty()
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
      .treated = .treated,
      .untreated = .untreated
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atu",
      focal = focal,
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
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL,
  focal = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atu.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .treated = .treated,
    .untreated = .untreated,
    focal = focal,
    ...
  )
}

#' @export
wt_atu.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...,
  focal = NULL
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
    .treated = .treated,
    .untreated = .untreated,
    focal = focal,
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
  abort_no_method(.propensity)
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
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- atm_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atm",
      focal = NULL,
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
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atm.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
  abort_no_method(.propensity)
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
  rlang::check_dots_empty()
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
      .treated = .treated,
      .untreated = .untreated
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "ato",
      focal = NULL,
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
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_ato.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
  abort_no_method(.propensity)
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
    check_lengths_match(.propensity, .exposure)
    check_ps_range(.propensity)
    wts <- entropy_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else if (exposure_type == "categorical") {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "entropy",
      focal = NULL,
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
  .treated = NULL,
  .untreated = NULL,
  ...,
  .propensity_col = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_entropy.numeric,
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  UseMethod("wt_cens")
}

#' @export
wt_cens.default <- function(
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
  abort_no_method(.propensity)
}

#' @export
wt_cens.numeric <- function(
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
  # Get weights using ATE formula
  wts <- wt_ate.numeric(
    .propensity = .propensity,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  exposure_type <- match_exposure_type(exposure_type, .exposure)

  # Handle optional exposure argument
  .exposure <- extract_exposure_from_glm(.propensity, .exposure)

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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trim",
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
wt_att.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trunc",
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
wt_att.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
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
  .treated = NULL,
  .untreated = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trim",
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
wt_cens.ps_trunc <- function(
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
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trunc",
    .sigma = .sigma,
    exposure_type = exposure_type,
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
  focal = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  call = rlang::caller_env()
) {
  # Ensure exposure is a factor
  .exposure <- transform_exposure_categorical(.exposure, focal, call = call)

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
      if (is.null(focal)) {
        abort(
          "Focal category must be specified for ATT with categorical exposures.",
          error_class = "propensity_focal_required_error"
        )
      }
      # h(e) = e_focal
      focal_idx <- which(levels_exp == focal)
      ps_matrix[, focal_idx]
    },
    "atu" = {
      if (is.null(focal)) {
        abort(
          "Focal category must be specified for ATU with categorical exposures.",
          error_class = "propensity_focal_required_error"
        )
      }
      # For ATU, we want weights for all non-focal categories
      # h(e) = 1 - e_focal
      focal_idx <- which(levels_exp == focal)
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
  if (!is.null(focal)) {
    attr(weights, "focal_category") <- focal
  }

  weights
}
