#' Inverse Probability Weighted Estimation
#'
#' @description
#' `ipw()` is a bring-your-own-model (BYOM) inverse probability weighted
#' estimator for causal inference. You supply a fitted propensity score model
#' and a fitted weighted outcome model; `ipw()` computes causal effect estimates
#' with standard errors that correctly account for the two-step estimation
#' process.
#'
#' `ipw()` supports binary, categorical, and continuous exposures. For a binary
#' or categorical exposure, a binary outcome returns the risk difference, log
#' risk ratio, and log odds ratio, and a continuous outcome returns the
#' difference in means. A continuous exposure is supplied through an [stats::lm()]
#' or identity-link gaussian [stats::glm()] propensity score model with a weighted
#' marginal structural outcome model whose link is identity, logit, or log, and
#' `ipw()` reports the single exposure coefficient of that model.
#'
#' @param ps_mod The weighting object: a fitted propensity score model that
#'   produced the weights, typically a logistic regression of class
#'   [stats::glm()] with the exposure as the left-hand side of the formula.
#'   `ipw()` is an S3 generic that dispatches on this object.
#' @param outcome_mod A fitted weighted outcome model of class [stats::glm()]
#'   or [stats::lm()], with the outcome as the dependent variable and
#'   propensity score weights supplied via the `weights` argument. The weights
#'   should be created with a propensity weight function such as [wt_ate()].
#'   Supported outcome models are an [stats::lm()], a gaussian
#'   [stats::glm()] with an identity link, and a binomial or quasibinomial
#'   [stats::glm()]; any other family (such as poisson or Gamma) or a
#'   non-identity gaussian link errors. A factor or logical outcome response is
#'   converted to `0`/`1` following glm's coding (the first factor level is the
#'   failure, every other level is the success).
#' @param .data A data frame containing the exposure, outcome, and covariates.
#'   If `NULL` (the default), `ipw()` extracts data from the model objects.
#'   Supply `.data` explicitly if the outcome model formula contains
#'   transformations that prevent extraction of the exposure variable from
#'   [stats::model.frame()], or if the propensity score model cannot reconstruct
#'   its design (for example, fit with `model = FALSE` with the fitting data no
#'   longer available); `ipw()` then rebuilds the propensity design from `.data`.
#' @param estimand A character string specifying the causal estimand: one of
#'   `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`. The available
#'   estimands depend on the exposure type: a binary or categorical exposure
#'   supports all six, while a continuous exposure supports only `"ate"`. For a
#'   categorical exposure, `"att"` and `"atu"` require a focal level (see
#'   `.focal_level`). See the **Estimands** section for the full support matrix.
#'   If `NULL`, the estimand is inferred from the weights in `outcome_mod`, which
#'   requires weights created with a propensity weight function such as
#'   [wt_ate()].
#' @param ps_link A character string specifying the link function used in the
#'   propensity score model: `"logit"`, `"probit"`, or `"cloglog"`. Defaults to
#'   the link used by `ps_mod`. This applies only to a binomial [stats::glm()]
#'   propensity score model on the binary path; leave it `NULL` for a
#'   multinomial or continuous propensity score model.
#' @param conf_level Confidence level for intervals. Default is `0.95`.
#' @param se_method Method for standard error estimation. `"mestimation"` (the
#'   default) stacks the propensity score, outcome, and estimand estimating
#'   equations and returns the empirical sandwich variance. `"linearization"`
#'   uses the influence-function linearization of Kostouraki et al. (2024). Both
#'   account for the uncertainty of estimating the propensity scores.
#'   `"linearization"` supports only an outcome model of the exposure alone; a
#'   covariate-adjusted outcome model requires `"mestimation"`.
#' @param ... Additional arguments. The `as.data.frame()` method passes these
#'   to [base::as.data.frame()]; the estimation methods do not currently use
#'   them and accept `...` for consistency with the `ipw()` generic.
#'
#' @details
#' # Workflow
#'
#' `ipw()` is designed around a three-step workflow:
#'
#' 1. Fit a propensity score model (e.g., logistic regression of exposure on
#'    confounders).
#' 2. Calculate propensity score weights for your estimand (e.g., [wt_ate()])
#'    and fit a weighted outcome model.
#' 3. Pass both models to `ipw()` to obtain causal effect estimates with
#'    correct standard errors.
#'
#' You are responsible for specifying and fitting both models. `ipw()` handles
#' the variance estimation.
#'
#' # Estimands
#'
#' The available estimands depend on the exposure type:
#'
#' | Estimand    | Binary | Categorical             | Continuous |
#' |-------------|:------:|:-----------------------:|:----------:|
#' | `"ate"`     | yes    | yes                     | yes        |
#' | `"att"`     | yes    | yes (needs focal level) | no         |
#' | `"atu"`     | yes    | yes (needs focal level) | no         |
#' | `"atm"`     | yes    | yes                     | no         |
#' | `"ato"`     | yes    | yes                     | no         |
#' | `"entropy"` | yes    | yes                     | no         |
#'
#' A continuous exposure supports only `"ate"`; any other request errors. For a
#' categorical exposure, `"att"` and `"atu"` require a focal level, supplied
#' through `.focal_level` or taken from the weights (see the `multinom` method).
#'
#' # Effect measures
#'
#' For a binary exposure, the reported measures depend on the outcome model. A
#' binary outcome ([stats::glm()] with `family = binomial()`) returns three
#' measures:
#' - `rd`: Risk difference (marginal risk in exposed minus unexposed)
#' - `log(rr)`: Log risk ratio
#' - `log(or)`: Log odds ratio
#'
#' A continuous outcome ([stats::lm()] or [stats::glm()] with
#' `family = gaussian()`) returns only the difference in means (`diff`).
#'
#' For a categorical exposure, the same measures are reported for each
#' non-reference level against the reference (first) factor level. The estimates
#' table gains a `comparison` column identifying each contrast (for example
#' `"b vs a"`), so a K-level exposure produces one block of measures per
#' non-reference level.
#'
#' For a continuous exposure, `ipw()` reports the single exposure coefficient of
#' the weighted marginal structural outcome model. Its label follows the outcome
#' link: `slope` for an identity link, `log(or)` for a logit link, and
#' `log(rr)` for a log link. The reported coefficient is the one attached to the
#' exposure term of the marginal structural model. When that term is a
#' transformation of the exposure (for example an `I(A^2)`-only model), the
#' coefficient of the transformed term is reported under the same link-based
#' label.
#'
#' Use [as.data.frame()] with `exponentiate = TRUE` to obtain risk ratios and
#' odds ratios on their natural scale.
#'
#' # Variance estimation
#'
#' By default (`se_method = "mestimation"`), standard errors are computed by
#' M-estimation. The propensity score model, the weighted outcome model, the
#' marginal means, and the effect contrasts are stacked as a single system of
#' estimating equations, and the empirical sandwich variance of that joint
#' system is used. The stacked marginal means are standardized to the estimand's
#' tilted target population, so a non-`ate` estimand with a covariate-adjusted
#' outcome model reports the contrast for that target population rather than the
#' full-sample average. Stacking the models accounts for the uncertainty
#' introduced by estimating the propensity scores, avoiding the underestimated
#' standard errors that arise from treating estimated weights as fixed. See
#' Stefanski
#' and Boos (2002) for the M-estimation framework.
#'
#' M-estimation standard errors are available for all exposure types: binary
#' (from a [stats::glm()] propensity score model), categorical (from a
#' [nnet::multinom()] model), and continuous (from an [stats::lm()] or
#' gaussian-family [stats::glm()] model). The `atm` weight `pmin(e, 1 - e)` is
#' not differentiable at a propensity score of `0.5`; deli's central finite
#' difference straddles the kink there and averages the one-sided slopes. Its
#' effect on the variance is negligible unless many observations sit at exactly
#' `0.5`.
#'
#' Setting `se_method = "linearization"` instead uses the influence-function
#' linearization of Kostouraki et al. (2024). It is available only for a binary
#' exposure and only for the `ate`, `att`, `ato`, and `atm` estimands; a
#' categorical or continuous exposure, or the `atu` or `entropy` estimand,
#' requires `se_method = "mestimation"`. The linearization also supports only an
#' outcome model whose formula contains the exposure alone. Its point estimates
#' come from g-computation on the fitted outcome model, while its influence
#' functions are those of the Hajek weighted-mean estimator; the two agree only
#' for an exposure-only outcome model. A covariate-adjusted outcome model (any
#' term beyond the exposure, including covariates, interactions, or transformed
#' terms) errors on the linearization path and requires
#' `se_method = "mestimation"`, which stacks adjusted outcome models correctly.
#'
#' # Model requirements
#'
#' `ipw()` cannot yet account for propensity scores that were trimmed
#' ([ps_trim()]), truncated ([ps_trunc()]), or calibrated ([ps_calibrate()])
#' before weighting. The stacked estimating equations rebuild the weights from
#' the propensity score model on every evaluation, so a weight that is no longer
#' a deterministic function of that model breaks the sandwich variance. Supplying
#' weights built from a modified score errors on either standard error method;
#' refit the weights from the unmodified propensity score model. An outcome model
#' fit without weights also errors on either method. The outcome model must not
#' carry an offset term on either method, since neither the stacked outcome score
#' nor the linearization influence functions thread an offset; supplying one
#' errors.
#'
#' With the default `se_method = "mestimation"`, two further requirements
#' apply. The outcome model weights must match the values implied by the
#' propensity score model; a mismatch errors. The propensity score model must be
#' fit without case weights, since the stacked propensity score equations are
#' unweighted and a weighted fit would not sit at the score root. Neither of
#' these two checks fires on the linearization path, which treats the weights as
#' fixed.
#'
#' @references
#' Stefanski LA, Boos DD. The calculus of M-estimation. *The American
#' Statistician*. 2002;56(1):29--38. \doi{10.1198/000313002753631330}
#'
#' Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
#' inverse probability-of-treatment weighting estimator: A tutorial for
#' different types of propensity score weights. *Statistics in Medicine*.
#' 2024;43(13):2672--2694. \doi{10.1002/sim.10078}
#'
#' @return Methods of `ipw()` return an S3 object of class `ipw`. This return
#'   contract is shared by every method and has the following components:
#' \describe{
#'   \item{`estimand`}{The causal estimand: one of `"ate"`, `"att"`, `"atu"`,
#'     `"atm"`, `"ato"`, or `"entropy"`.}
#'   \item{`ps_mod`}{The weighting object: the fitted model that produced the
#'     weights.}
#'   \item{`outcome_mod`}{The fitted outcome model.}
#'   \item{`estimates`}{A data frame with one row per effect measure and the
#'     following columns: `effect` (the measure name), `estimate` (point
#'     estimate), `std.err` (standard error), `z` (z-statistic),
#'     `ci.lower` and `ci.upper` (confidence interval bounds),
#'     `conf.level`, and `p.value`. For a categorical exposure the data frame
#'     also has a `comparison` column, placed after `effect`, naming the
#'     non-reference level and reference level of each contrast.}
#'   \item{`se_method`}{The standard error method used, `"mestimation"` or
#'     `"linearization"`.}
#'   \item{`fit`}{The fitted M-estimator object when `se_method` is
#'     `"mestimation"`, otherwise `NULL`. Calling [stats::vcov()] or
#'     [generics::tidy()] on this object exposes the full stacked system of
#'     estimating equations, including the propensity score, outcome, and
#'     estimand parameters.}
#' }
#'
#' @examples
#' # Simulate data with a confounder, binary exposure, and binary outcome
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.5 * x1))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
#' dat <- data.frame(x1, z, y)
#'
#' # Step 1: Fit a propensity score model
#' ps_mod <- glm(z ~ x1, data = dat, family = binomial())
#'
#' # Step 2: Calculate ATE weights and fit a weighted outcome model
#' wts <- wt_ate(ps_mod)
#' outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)
#'
#' # Step 3: Estimate causal effects with correct standard errors
#' result <- ipw(ps_mod, outcome_mod)
#' result
#'
#' # Exponentiate log-RR and log-OR to get RR and OR
#' as.data.frame(result, exponentiate = TRUE)
#'
#' # Continuous outcome example
#' y_cont <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
#' dat$y_cont <- y_cont
#' outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
#' ipw(ps_mod, outcome_cont)
#'
#' # Continuous exposure: an lm propensity model of the dose on covariates,
#' # stabilized weights, and a weighted marginal structural outcome model
#' a <- 0.5 + 0.8 * x1 + rnorm(n)
#' y_dose <- 1 + 0.6 * a + 0.3 * x1 + rnorm(n)
#' dat$a <- a
#' dat$y_dose <- y_dose
#' ps_cont <- lm(a ~ x1, data = dat)
#' wts_cont <- wt_ate(
#'   fitted(ps_cont),
#'   a,
#'   exposure_type = "continuous",
#'   stabilize = TRUE
#' )
#' msm <- lm(y_dose ~ a, data = dat, weights = wts_cont)
#' ipw(ps_cont, msm)
#'
#' @examplesIf requireNamespace("nnet", quietly = TRUE)
#' # Categorical exposure: a multinomial propensity model and per-level contrasts
#' set.seed(2)
#' n <- 300
#' x1 <- rnorm(n)
#' a <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.3 * x1))
#' dat_cat <- data.frame(x1, a, y)
#'
#' ps_cat <- nnet::multinom(a ~ x1, data = dat_cat, trace = FALSE)
#' ps_mat <- predict(ps_cat, type = "probs")
#' wts_cat <- wt_ate(ps_mat, dat_cat$a, exposure_type = "categorical")
#' outcome_cat <- glm(y ~ a, data = dat_cat, family = binomial(), weights = wts_cat)
#' ipw(ps_cat, outcome_cat)
#'
#' @seealso
#' [wt_ate()], [wt_att()], [wt_atu()], [wt_atm()], [wt_ato()], [wt_entropy()]
#'   for calculating propensity score weights.
#'
#' [ps_trim()], [ps_trunc()] for handling extreme propensity scores before
#'   weighting.
#'
#' @name ipw-methods
#' @export
#' @importFrom stats dnorm family formula model.frame model.matrix model.weights
#' @importFrom stats pnorm predict printCoefmat qnorm var
ipw.glm <- function(
  ps_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization")
) {
  se_method <- rlang::arg_match(se_method)
  assert_class(ps_mod, "glm")
  assert_class(outcome_mod, c("glm", "lm"))

  # A gaussian-family propensity model indicates a continuous exposure. Route it
  # to the shared continuous path, which an lm propensity model also uses; the two
  # share fitted values and so produce identical estimates and standard errors.
  if (identical(ps_mod$family$family, "gaussian")) {
    return(ipw_continuous_estimate(
      ps_mod,
      outcome_mod,
      .data = .data,
      estimand = estimand,
      conf_level = conf_level,
      se_method = se_method
    ))
  }

  # Guards first, on the weights that fit the outcome model. These carry the
  # psw attributes, so a modified propensity score is detected here before any
  # length reconciliation or estimand parsing that would otherwise fail
  # obliquely. Both standard error methods share these guards.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts)

  if (identical(se_method, "mestimation")) {
    spec <- ipw_spec_binary(
      ps_mod,
      outcome_mod,
      .data = .data,
      estimand = estimand,
      ps_link = ps_link
    )
    fit <- ipw_mestimation(spec, conf_level = conf_level)
    return(new_ipw(
      estimand = spec$estimand,
      ps_mod = ps_mod,
      outcome_mod = outcome_mod,
      estimates = fit$estimates,
      se_method = "mestimation",
      fit = fit$fit
    ))
  }

  weight_matrix <- model.matrix(ps_mod)
  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  check_ipw_offset(outcome_mod)
  check_ipw_linearization_outcome(outcome_mod, exposure_name)
  check_ipw_outcome_family(outcome_mod)

  if (is.null(.data)) {
    exposure <- fmla_extract_left_vctr(ps_mod)
    outcome <- fmla_extract_left_vctr(outcome_mod)
  } else {
    assert_class(exposure_name, "character", .length = 1)
    assert_class(outcome_name, "character", .length = 1)
    assert_columns_exist(.data, c(exposure_name, outcome_name))

    exposure <- .data[[exposure_name]]
    outcome <- .data[[outcome_name]]
  }

  # Convert a factor or logical response to its 0/1 coding so the influence
  # values compute Y - mu on numeric values rather than factor level codes.
  # Covers both extraction routes above.
  outcome <- ipw_outcome_numeric(outcome)

  ps <- predict(ps_mod, type = "response", newdata = .data)

  if (is.null(ps_link)) {
    ps_link <- ps_mod$family$link
  }

  if (!identical(length(exposure), length(outcome))) {
    abort(c(
      "{.arg exposure} and {.arg outcome} must be the same length.",
      x = "{.arg exposure} is length {length(exposure)}",
      x = "{.arg outcome} is length {length(outcome)}"
    ))
  }

  estimand <- check_estimand(wts, estimand)

  # The linearization influence functions are derived only for ate, att, ato, and
  # atm; atu and entropy require the mestimation path. Reject them here with the
  # documented classed error rather than letting the request reach derive_weights,
  # whose internal arg_match raises a bare, misleading error.
  if (!estimand %in% c("ate", "att", "ato", "atm")) {
    abort(
      c(
        "{.fun ipw} does not support {.val linearization} standard errors for \\
        the {.val {estimand}} estimand.",
        i = "Use {.code se_method = \"mestimation\"} for the {.val {estimand}} \\
        estimand."
      ),
      error_class = "propensity_method_error"
    )
  }

  marginal_means <- estimate_marginal_means(
    outcome_mod = outcome_mod,
    wts = wts,
    exposure = exposure,
    exposure_name = exposure_name,
    .data = .data
  )

  # Recode the exposure to 0/1 for the influence functions, in parity with the
  # M-estimation recode. The marginal means above keep the original exposure so
  # they set the counterfactual levels the outcome model expects; the influence
  # values below need the numeric indicator (wts * Z, exposure == 1, Z - ps).
  exposure_binary <- ipw_recode_binary_exposure(exposure)

  uncorrected_lin_vars <- linearize_variables_for_wts(
    exposure_binary,
    outcome,
    wts,
    marginal_means
  )

  lin_vars <- linearize_variables_for_ps(
    exposure = exposure_binary,
    outcome = outcome,
    wts = wts,
    ps = ps,
    estimand = estimand,
    weight_matrix = weight_matrix,
    marginal_means = marginal_means,
    uncorrected_lin_vars = uncorrected_lin_vars,
    ps_link = ps_link
  )

  estimates <- calculate_estimates(
    lin_vars = lin_vars,
    marginal_means = marginal_means,
    n = length(outcome),
    linear_regression = is_linear_regression(outcome_mod),
    conf_level = conf_level
  )

  new_ipw(
    estimand = estimand,
    ps_mod = ps_mod,
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = "linearization",
    fit = NULL
  )
}

# Assemble the ipw() return object. The four core fields are shared by every
# method; se_method records which variance estimator ran and fit holds the
# M-estimator object (NULL on the linearization path).
new_ipw <- function(estimand, ps_mod, outcome_mod, estimates, se_method, fit) {
  structure(
    list(
      estimand = estimand,
      ps_mod = ps_mod,
      outcome_mod = outcome_mod,
      estimates = estimates,
      se_method = se_method,
      fit = fit
    ),
    class = "ipw"
  )
}

# Shared continuous-exposure route for ipw.lm and the gaussian-family branch of
# ipw.glm. An lm and a gaussian glm propensity model share fitted values, so both
# reach the same M-estimation fit. Standard errors come only from M-estimation
# for a continuous exposure; the linearization path is not available.
ipw_continuous_estimate <- function(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  conf_level = 0.95,
  se_method = "mestimation",
  call = rlang::caller_env()
) {
  if (identical(se_method, "linearization")) {
    abort(
      c(
        "{.fun ipw} does not support {.val linearization} standard errors for \\
        continuous exposures.",
        i = "Use {.code se_method = \"mestimation\"} for a continuous exposure."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  # Guard the weights that fit the outcome model before building the spec, so a
  # modified propensity score is detected before any estimand parsing.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts, call = call)

  spec <- ipw_spec_continuous(
    ps_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    call = call
  )
  fit <- ipw_mestimation(spec, conf_level = conf_level, call = call)

  new_ipw(
    estimand = spec$estimand,
    ps_mod = ps_mod,
    outcome_mod = outcome_mod,
    estimates = fit$estimates,
    se_method = "mestimation",
    fit = fit$fit
  )
}

# Guard the weights that fit the outcome model. ipw() cannot yet account for
# propensity scores that were trimmed, truncated, or calibrated before
# weighting, so detect them here and direct the user to refit from the
# unmodified model. An outcome model fit without weights cannot yield an IPW
# estimate at all.
check_ipw_weights <- function(wts, call = rlang::caller_env()) {
  if (is.null(wts)) {
    abort(
      c(
        "{.arg outcome_mod} must be fit with propensity score weights.",
        x = "No {.arg weights} were found in {.arg outcome_mod}.",
        i = "Fit {.arg outcome_mod} with weights from a propensity weight \\
        function such as {.fun wt_ate}."
      ),
      error_class = "propensity_ipw_weights_missing_error",
      call = call
    )
  }

  if (is_ps_trimmed(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with trimmed propensity score weights.",
        x = "{.fun ipw} cannot yet account for trimmed propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_trimmed_error",
      call = call
    )
  }

  if (is_ps_truncated(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with truncated propensity score weights.",
        x = "{.fun ipw} cannot yet account for truncated propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_truncated_error",
      call = call
    )
  }

  if (is_ps_calibrated(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with calibrated propensity score weights.",
        x = "{.fun ipw} cannot yet account for calibrated propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_calibrated_error",
      call = call
    )
  }

  invisible(wts)
}

# Reject an outcome model that carries an offset. The stacked estimating
# equations rebuild the weighted outcome score without threading an offset, so
# an offset in the outcome model would converge to a root that disagrees with
# the fitted model. Detect an offset supplied through the model formula (the
# terms offset attribute) or through the `offset` argument (the model frame) and
# direct the user to refit without it. Shared by every mestimation spec path.
check_ipw_offset <- function(outcome_mod, call = rlang::caller_env()) {
  has_offset <- !is.null(attr(stats::terms(outcome_mod), "offset")) ||
    !is.null(stats::model.offset(stats::model.frame(outcome_mod)))

  if (has_offset) {
    abort(
      c(
        "{.arg outcome_mod} was fit with an offset term.",
        x = "{.fun ipw} does not support offset terms in the outcome model; \\
        the stacked estimating equations do not thread an offset.",
        i = "Refit {.arg outcome_mod} without the offset."
      ),
      error_class = "propensity_ipw_offset_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject an outcome model whose family or link the estimating equations cannot
# stack. Both standard error paths classify the outcome model as either a
# gaussian identity linear score or a binomial score with the model's link, so
# only three shapes are supported: an lm, a gaussian glm with an identity link,
# and a binomial or quasibinomial glm with a link the outcome score reconstructs
# (the links ipw_inv_link handles). Any other family (poisson, quasipoisson,
# Gamma, inverse gaussian, ...) would otherwise be silently stacked as a binomial
# score, solving to a different root than the fitted model; a non-identity
# gaussian link would be silently treated as identity. Shared by every mestimation
# spec path and the linearization entry.
check_ipw_outcome_family <- function(outcome_mod, call = rlang::caller_env()) {
  # An lm (not a glm) is a gaussian identity linear model.
  if (!inherits(outcome_mod, "glm")) {
    return(invisible(TRUE))
  }

  family <- outcome_mod$family$family
  link <- outcome_mod$family$link
  supported_links <- c("logit", "probit", "cloglog", "log", "identity")

  if (identical(family, "gaussian")) {
    if (!identical(link, "identity")) {
      abort(
        c(
          "{.fun ipw} supports only an identity link for a gaussian outcome \\
          model.",
          x = "{.arg outcome_mod} is a gaussian model with a {.val {link}} link.",
          i = "Refit {.arg outcome_mod} with a gaussian identity link, or use a \\
          binomial or quasibinomial outcome model."
        ),
        error_class = "propensity_ipw_family_error",
        call = call
      )
    }
    return(invisible(TRUE))
  }

  if (family %in% c("binomial", "quasibinomial")) {
    if (!link %in% supported_links) {
      abort(
        c(
          "{.fun ipw} does not support a {.val {link}} link for a \\
          {.val {family}} outcome model.",
          x = "{.arg outcome_mod} was fit with a {.val {link}} link.",
          i = "Refit {.arg outcome_mod} with one of {.val {supported_links}}."
        ),
        error_class = "propensity_ipw_family_error",
        call = call
      )
    }
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.val {family}} outcome models.",
      x = "{.arg outcome_mod} was fit with the {.val {family}} family.",
      i = "Fit {.arg outcome_mod} with a binomial or quasibinomial family for a \\
      binary outcome, or a gaussian identity link (or an {.fun lm}) for a \\
      continuous outcome."
    ),
    error_class = "propensity_ipw_family_error",
    call = call
  )
}

# Validate the propensity score and outcome links on the continuous exposure
# path, which the binary and categorical paths do not restrict this way. Both
# checks fire at spec-construction entry, before any M-estimator solve. The
# continuous propensity score design is reconstructed as the linear predictor
# X alpha, which equals the fitted mean only under an identity link; a
# non-identity gaussian propensity link would otherwise surface later as a
# misleading weights mismatch. The reported effect is labeled from the outcome
# link (identity for a slope, logit for a log odds ratio, log for a log risk
# ratio), so an outcome link outside that set has no label.
check_ipw_continuous_links <- function(
  ps_mod,
  outcome_mod,
  call = rlang::caller_env()
) {
  if (inherits(ps_mod, "glm")) {
    ps_link <- ps_mod$family$link
    if (!identical(ps_link, "identity")) {
      abort(
        c(
          "{.fun ipw} supports only an identity-link propensity score model \\
          for a continuous exposure.",
          x = "{.arg ps_mod} is a gaussian model with a {.val {ps_link}} link.",
          i = "Refit {.arg ps_mod} as an {.fun lm} or a gaussian glm with an \\
          identity link."
        ),
        error_class = "propensity_ipw_link_error",
        call = call
      )
    }
  }

  outcome_link <- if (inherits(outcome_mod, "glm")) {
    outcome_mod$family$link
  } else {
    "identity"
  }
  supported <- c("identity", "logit", "log")
  if (!outcome_link %in% supported) {
    abort(
      c(
        "{.fun ipw} does not support a {.val {outcome_link}} link for the \\
        marginal structural model of a continuous exposure.",
        x = "{.arg outcome_mod} was fit with a {.val {outcome_link}} link.",
        i = "Refit {.arg outcome_mod} with one of {.val {supported}}."
      ),
      error_class = "propensity_ipw_link_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Restrict the linearization path to an outcome model of the exposure alone. Its
# point estimates come from g-computation on the fitted weighted outcome model,
# but the influence functions in linearize_variables_for_wts and correct_for_ps
# are those of the Hajek weighted-mean estimator. The two agree only when the
# outcome model regresses on the exposure by itself; any covariate, interaction,
# or transformed term makes the reported standard error that of a different
# estimator. The exposure-only condition is that the outcome model's term labels
# are exactly the exposure name. Detect an adjusted outcome model here and direct
# the user to the mestimation path, which stacks adjusted outcome models
# correctly. Applies only to the linearization path.
check_ipw_linearization_outcome <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  term_labels <- attr(stats::terms(outcome_mod), "term.labels")

  if (!identical(term_labels, exposure_name)) {
    abort(
      c(
        "{.fun ipw} supports {.val linearization} standard errors only for an \\
        outcome model of the exposure alone.",
        x = "{.arg outcome_mod} is adjusted for terms beyond \\
        {.val {exposure_name}}.",
        i = "Use {.code se_method = \"mestimation\"} for a covariate-adjusted \\
        outcome model."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  invisible(TRUE)
}

#' @export
ipw.default <- function(
  ps_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization")
) {
  abort(
    c(
      "{.fun ipw} does not know how to handle {.arg ps_mod} of class \\
      {.cls {class(ps_mod)}}.",
      i = "{.arg ps_mod} must be a fitted propensity score model of class \\
      {.cls glm}, such as a logistic regression."
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
print.ipw <- function(x, ...) {
  cat("Inverse Probability Weight Estimator\n")
  # todo: make this more adaptive, e.g. ATE without L2FU
  cat("Estimand:", toupper(x$estimand), "\n\n")

  cat("Propensity Score Model:\n")
  cat("  Call:", format_model_call(x$ps_mod), "\n")
  cat("\n")

  cat("Outcome Model:\n")
  cat("  Call:", format_model_call(x$outcome_mod), "\n")

  cat("\n")

  cat("Estimates:\n")
  if ("comparison" %in% names(x$estimates)) {
    # A categorical exposure repeats effect labels across comparisons, so the
    # printed rows are keyed by effect and comparison together and the character
    # comparison column is dropped from the numeric matrix printCoefmat formats.
    estimates <- x$estimates[setdiff(
      names(x$estimates),
      c("effect", "comparison")
    )]
    rownames(estimates) <- paste(x$estimates$effect, x$estimates$comparison)
  } else {
    estimates <- x$estimates[-1]
    rownames(estimates) <- x$estimates$effect
  }
  printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 2:3, tst.ind = 4)

  invisible(x)
}

# Format a model's originating call for the ipw() summary. Objects that carry an
# accessible call, such as glm and lm, print the deparsed call; a weighting
# object that does not expose one, including an S7 object that cannot be
# subset, falls back to a class label so print() works for any weighting object.
format_model_call <- function(mod) {
  call <- tryCatch(stats::getCall(mod), error = function(e) NULL)
  if (is.null(call)) {
    return(paste0("<", paste(class(mod), collapse = "/"), ">"))
  }
  paste(deparse(call), collapse = "\n")
}


#' @param x An `ipw` object.
#' @param exponentiate If `TRUE`, exponentiate the log risk ratio and log odds
#'   ratio to produce risk ratios and odds ratios on their natural scale. The
#'   confidence interval bounds are also exponentiated. Standard errors, z
#'   statistics, and p-values remain on the log scale. Default is `FALSE`.
#' @param row.names,optional Passed to [base::as.data.frame()].
#' @returns `as.data.frame()` returns the `estimates` component as a data
#'   frame. When `exponentiate = TRUE`, the `log(rr)` and `log(or)` rows are
#'   transformed: point estimates and confidence limits are exponentiated and
#'   the effect labels become `"rr"` and `"or"`.
#' @export
#' @rdname ipw-methods
as.data.frame.ipw <- function(
  x,
  row.names = NULL,
  optional = NULL,
  exponentiate = FALSE,
  ...
) {
  df <- as.data.frame(
    x$estimates,
    row.names = row.names,
    optional = optional,
    ...
  )

  if (!exponentiate) {
    # Return as-is
    return(df)
  }

  # If exponentiate=TRUE, we transform the "log_rr" and "log_or" rows to the raw scale
  # by exponentiating estimate, ci.lower, and ci.upper.

  is_log_rr <- df$effect == "log(rr)"
  is_log_or <- df$effect == "log(or)"

  rows_to_expo <- is_log_rr | is_log_or

  # Exponentiate estimate, ci.lower, ci.upper
  df$estimate[rows_to_expo] <- exp(df$estimate[rows_to_expo])
  df$ci.lower[rows_to_expo] <- exp(df$ci.lower[rows_to_expo])
  df$ci.upper[rows_to_expo] <- exp(df$ci.upper[rows_to_expo])

  # Rename effect labels for clarity
  # "log_rr" -> "rr"
  # "log_or" -> "or"
  df$effect[is_log_rr] <- "rr"
  df$effect[is_log_or] <- "or"

  # note: variance, std.err, z_stat, p_value remain on the log scale.
  # significance testing is typically done on log-scale.

  df
}

calculate_estimates <- function(
  lin_vars,
  marginal_means,
  n,
  linear_regression,
  conf_level
) {
  z_val <- qnorm(1 - ((1 - conf_level) / 2))

  ### RISK DIFFERENCE (raw scale)
  # --------------------------------
  # Influence = (l1 - l0)
  rd_var <- var(lin_vars$l1 - lin_vars$l0) / n

  rd_est <- marginal_means$mu1 - marginal_means$mu0
  rd_se <- sqrt(rd_var)

  rd_ci_lower <- rd_est - z_val * rd_se
  rd_ci_upper <- rd_est + z_val * rd_se

  rd_z <- rd_est / rd_se
  rd_p <- 2 * (1 - pnorm(abs(rd_z)))

  # for continuous outcomes, only return difference
  if (linear_regression) {
    return(
      data.frame(
        effect = "diff",
        estimate = rd_est,

        # variances are on the same scale as 'estimate':
        # variance = c(rd_var, log_rr_var, log_or_var),
        std.err = rd_se,
        z = rd_z,
        ci.lower = rd_ci_lower,
        ci.upper = rd_ci_upper,
        conf.level = conf_level,
        p.value = rd_p
      )
    )
  }

  ### RISK RATIO (log scale)
  # ---------------------------
  # Risk ratio: RR = mu1 / mu0
  # We'll store 'log_rr_est = log(RR)' in the output.

  rr_raw_est <- marginal_means$mu1 / marginal_means$mu0
  log_rr_est <- log(rr_raw_est)

  # Influence function for log(RR) via the delta method:
  #   d log(mu1) - d log(mu0) = l1 / mu1 - l0 / mu0
  log_rr_inf <- lin_vars$l1 /
    marginal_means$mu1 -
    lin_vars$l0 / marginal_means$mu0

  log_rr_var <- var(log_rr_inf) / n
  log_rr_se <- sqrt(log_rr_var)

  log_rr_ci_lower <- log_rr_est - z_val * log_rr_se
  log_rr_ci_upper <- log_rr_est + z_val * log_rr_se

  rr_z <- log_rr_est / log_rr_se
  rr_p <- 2 * (1 - pnorm(abs(rr_z)))

  ### ODDS RATIO (log scale)
  # ---------------------------
  # OR = [mu1/(1 - mu1)] / [mu0/(1 - mu0)]

  or_raw_est <- (marginal_means$mu1 / (1 - marginal_means$mu1)) /
    (marginal_means$mu0 / (1 - marginal_means$mu0))
  log_or_est <- log(or_raw_est)

  # Influence function for log(OR) via the delta method:
  #   d logit(mu1) - d logit(mu0)
  #     = l1 / [mu1*(1 - mu1)] - l0 / [mu0*(1 - mu0)]
  log_or_inf <- (lin_vars$l1 /
    (marginal_means$mu1 * (1 - marginal_means$mu1))) -
    (lin_vars$l0 / (marginal_means$mu0 * (1 - marginal_means$mu0)))

  log_or_var <- var(log_or_inf) / n
  log_or_se <- sqrt(log_or_var)

  log_or_ci_lower <- log_or_est - z_val * log_or_se
  log_or_ci_upper <- log_or_est + z_val * log_or_se

  or_z <- log_or_est / log_or_se
  or_p <- 2 * (1 - pnorm(abs(or_z)))

  data.frame(
    effect = c("rd", "log(rr)", "log(or)"),
    # For RD, the estimate is raw. For RR and OR, the estimate is log-scale:
    estimate = c(rd_est, log_rr_est, log_or_est),

    # Variances are on the same scale as 'estimate'
    # but we won't return these since they can be calculated from std.err
    # variance = c(rd_var, log_rr_var, log_or_var)
    std.err = c(rd_se, log_rr_se, log_or_se),
    z = c(rd_z, rr_z, or_z),
    ci.lower = c(rd_ci_lower, log_rr_ci_lower, log_or_ci_lower),
    ci.upper = c(rd_ci_upper, log_rr_ci_upper, log_or_ci_upper),
    conf.level = conf_level,
    p.value = c(rd_p, rr_p, or_p)
  )
}

# accounts for dependence introduced by weights
# treats them as fixed, not estimated
# equivalent to usual robust sandwich estimator
linearize_variables_for_wts <- function(Z, Y, wts, marginal_means) {
  n <- length(Z)
  wts <- as.double(wts)
  l1 <- n / marginal_means$n1 * (wts * Z * (Y - marginal_means$mu1))
  l0 <- n / marginal_means$n0 * (wts * (1 - Z) * (Y - marginal_means$mu0))

  list(l1 = l1, l0 = l0)
}


# additionally adds variability for the estimation of the PS
linearize_variables_for_ps <- function(
  exposure,
  outcome,
  wts,
  ps,
  estimand,
  weight_matrix,
  marginal_means,
  uncorrected_lin_vars,
  ps_link
) {
  n <- length(exposure)
  score_factor_f <- derive_score_factor(ps_link)
  weight_derivatives <- derive_weights(
    exposure = exposure,
    ps = ps,
    weight_matrix = weight_matrix,
    ps_link = ps_link,
    estimand = estimand
  )

  correction_mat <- compute_ps_correction_matrix_inv(
    weight_matrix = weight_matrix,
    ps = ps,
    ps_link = ps_link,
    n = n
  )

  l1 <- uncorrected_lin_vars$l1 +
    correct_for_ps(
      exposure = exposure,
      outcome = outcome,
      ps = ps,
      mu = marginal_means$mu1,
      n_group = marginal_means$n1,
      weight_matrix = weight_matrix,
      weight_derivatives = weight_derivatives,
      correction_mat = correction_mat,
      score_factor_f = score_factor_f,
      n = n
    )

  l0 <- uncorrected_lin_vars$l0 +
    correct_for_ps(
      exposure = exposure,
      exposure_actual = 1 - exposure,
      outcome = outcome,
      ps = ps,
      mu = marginal_means$mu0,
      n_group = marginal_means$n0,
      weight_matrix = weight_matrix,
      weight_derivatives = weight_derivatives,
      correction_mat = correction_mat,
      score_factor_f = score_factor_f,
      n = n
    )

  list(l1 = l1, l0 = l0)
}

compute_ps_correction_matrix_inv <- function(weight_matrix, ps, ps_link, n) {
  deriv_link_f <- derive_link(ps_link)

  # row-by-row x_i x_i^T
  crossprods <- t_tcrossprod_over_rows(weight_matrix)

  # build the correction matrix
  correction_mat <- matrix(
    colSums(crossprods / (ps * (1 - ps) * (deriv_link_f(ps)^2))) / n,
    ncol(weight_matrix),
    ncol(weight_matrix)
  )

  solve(correction_mat)
}

# faster than `t(apply(mat, 1, \(x) tcrossprod(x)))`
t_tcrossprod_over_rows <- function(mat) {
  n <- nrow(mat)
  p <- ncol(mat)
  out <- matrix(0, nrow = n, ncol = p * p)
  for (i in seq_len(n)) {
    out[i, ] <- tcrossprod(mat[i, ])
  }

  out
}

correct_for_ps <- function(
  exposure,
  exposure_actual = exposure,
  outcome,
  ps,
  mu,
  n_group,
  weight_matrix,
  weight_derivatives,
  correction_mat,
  score_factor_f,
  n
) {
  # first, compute partial-derivative sums over subjects (averaged by n)
  partial_derivative_sums <- colSums(
    weight_derivatives * exposure_actual * (outcome - mu)
  ) /
    n

  # then build the transformation matrix for correction. The score of a binomial
  # GLM with link g is x_i (Z_i - p_i) / (p_i (1 - p_i) g'(p_i)); the score factor
  # 1 / (p (1 - p) g'(p)) generalizes the correction to any link. It is exactly 1
  # for the canonical logit, leaving the logit influence values unchanged.
  transformation_mat <- correction_mat %*%
    t((exposure - ps) * score_factor_f(ps) * weight_matrix)

  # and then apply the partial-derivative sums to that transformation
  correction_contrib <- rbind(partial_derivative_sums) %*% transformation_mat

  # rescale by (n / n_group)
  scaling_factor <- n / n_group
  correction_contrib <- correction_contrib * scaling_factor

  # and reduce to vector and unname
  correction_contrib |>
    drop() |>
    unname()
}

estimate_marginal_means <- function(
  outcome_mod,
  wts,
  exposure,
  exposure_name,
  .data = NULL,
  call = rlang::caller_env()
) {
  # todo: this could be generalized with split() and lapply()
  if (is.null(.data)) {
    .data <- model.frame(outcome_mod)
    check_exposure(.data, exposure_name, call = call)
  }
  # todo: make this more flexible for different values and model specs
  # maybe can optionally accept a function for g-comp
  exposure_values <- sort(unique(exposure))

  if (!isTRUE(length(exposure_values) == 2)) {
    abort(
      c(
        "{.code ipw()} currently only supports binary exposures.",
        x = "There are {length(exposure_values)} unique value{?s} of the \\
        exposure."
      ),
      call = call
    )
  }

  .data_1 <- .data
  .data_0 <- .data
  .data_1[[exposure_name]] <- exposure_values[[2]]
  .data_0[[exposure_name]] <- exposure_values[[1]]

  n1 <- sum(wts[exposure == exposure_values[[2]]])
  mu1 <- mean(predict(outcome_mod, newdata = .data_1, type = "response"))

  n0 <- sum(wts[exposure == exposure_values[[1]]])
  mu0 <- mean(predict(outcome_mod, newdata = .data_0, type = "response"))

  list(
    # exposure = 1
    n1 = n1,
    mu1 = mu1,
    # exposure = 0
    n0 = n0,
    mu0 = mu0
  )
}

derive_weights <- function(
  exposure,
  ps,
  weight_matrix,
  ps_link = c("logit", "probit", "cloglog"),
  estimand = c("ate", "att", "ato", "atm")
) {
  estimand <- rlang::arg_match(estimand)
  deriv_link_f <- derive_link(ps_link)

  deriv_vec <- switch(
    estimand,
    "ate" = ate_derivative(exposure, ps, deriv_link_f),
    "att" = att_derivative(exposure, ps, deriv_link_f),
    "ato" = ato_derivative(exposure, ps, deriv_link_f),
    "atm" = atm_derivative(exposure, ps, deriv_link_f)
  )

  weight_matrix * deriv_vec
}

ate_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    -1 / (ps^2 * deriv_link_f(ps)),
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
}

att_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    0,
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
}

ato_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    -1 / deriv_link_f(ps),
    1 / deriv_link_f(ps)
  )
}

atm_derivative <- function(exposure, ps, deriv_link_f) {
  wv <- ifelse(
    exposure == 1,
    -1 / (ps^2 * deriv_link_f(ps)),
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
  # Then set parts to zero
  wv[exposure == 1 & ps < 0.5] <- 0
  wv[exposure == 0 & ps > 0.5] <- 0
  wv
}


derive_link <- function(ps_link = c("logit", "probit", "cloglog")) {
  ps_link <- rlang::arg_match(ps_link)
  switch(
    ps_link,
    logit = function(x) 1 / (x * (1 - x)),
    probit = function(x) 1 / (dnorm(qnorm(x))),
    cloglog = function(x) -1 / ((1 - x) * log(1 - x))
  )
}

# Score factor 1 / (p (1 - p) g'(p)) that scales the residual (Z - p) into the
# binomial GLM score on the linearization ps correction. Returned in closed form
# per link so the logit factor is exactly 1 (not the round-trip
# p (1 - p) / (p (1 - p)) of derive_link), keeping the logit influence values
# bit-for-bit unchanged. Probit and cloglog match 1 / (p (1 - p) g'(p)) with
# g'(p) from derive_link().
derive_score_factor <- function(ps_link = c("logit", "probit", "cloglog")) {
  ps_link <- rlang::arg_match(ps_link)
  switch(
    ps_link,
    logit = function(p) rep(1, length(p)),
    probit = function(p) dnorm(qnorm(p)) / (p * (1 - p)),
    cloglog = function(p) -log(1 - p) / p
  )
}

extract_weights <- function(.mod) {
  .mod |>
    model.frame() |>
    model.weights()
}

# Convert an outcome response to the 0/1 numeric coding a binomial glm uses. A
# factor maps every non-first level to success (glm treats the first level as
# failure), matching outcome_mod$y; a logical or numeric response doubles
# directly. Applied to the raw outcome vector at each extraction site (a model
# frame response or a .data column) so that a factor or logical response stacks
# and linearizes as its converted value rather than as level codes.
ipw_outcome_numeric <- function(outcome) {
  if (is.factor(outcome)) {
    as.double(outcome != levels(outcome)[[1L]])
  } else {
    as.double(outcome)
  }
}

# Recode a binary exposure to 0/1 for the influence-function computations, taking
# the second sorted value as the exposed group (coded 1): a factor's second
# level, TRUE for a logical, or the larger value for a numeric. Shared by the
# M-estimation spec and the linearization path so the two recode identically. A
# numeric 0/1 exposure is returned unchanged (sort(unique(c(0, 1)))[[2]] is 1).
# An exposure with more than two unique values is returned untouched so the
# downstream two-level check still fires on the original values.
ipw_recode_binary_exposure <- function(exposure) {
  exposure_values <- sort(unique(exposure))
  if (length(exposure_values) == 2) {
    as.double(exposure == exposure_values[[2]])
  } else {
    exposure
  }
}

check_estimand <- function(wts, estimand, call = rlang::caller_env()) {
  if (is_causal_wt(wts)) {
    estimand_from_weights <- estimand(wts)
  } else {
    estimand_from_weights <- NULL
  }

  if (!is.null(estimand_from_weights) && !is.null(estimand)) {
    same_estimand <- identical(estimand_from_weights, estimand)
    if (!same_estimand) {
      .estimand <- estimand
      .estimand_from_weights <- estimand_from_weights
      abort(
        "Estimand in weights different from {.arg estimand}: \\
        {.val { .estimand_from_weights}} vs. {.val { .estimand}}",
        call = call
      )
    } else {
      return(estimand)
    }
  }

  if (is.null(estimand_from_weights) && is.null(estimand)) {
    abort(
      c(
        "Can't determine the estimand from weights.",
        i = "Please specify {.arg estimand}."
      ),
      call = call
    )
  }

  if (!is.null(estimand_from_weights)) {
    return(estimand_from_weights)
  } else {
    return(estimand)
  }
}

check_exposure <- function(.data, .exposure_name, call = rlang::caller_env()) {
  assert_class(.exposure_name, "character", .length = 1, call = call)
  if (!(.exposure_name %in% names(.data))) {
    abort(
      c(
        "{.val { .exposure_name}} not found in {.code model.frame(outcome_mod)}.",
        i = "The outcome model may have transformations in the formula.",
        i = "Please specify {.arg .data}"
      ),
      call = call,
      error_class = "propensity_columns_exist_error"
    )
  }
}

is_linear_regression <- function(outcome_mod) {
  # Check if 'outcome_mod' inherits from 'glm' first
  if (inherits(outcome_mod, "glm")) {
    # 'outcome_mod' is a glm, check if family is 'gaussian'
    if (family(outcome_mod)$family == "gaussian") {
      # It's a glm with Gaussian family => linear regression
      return(TRUE)
    } else {
      # It's a glm but not Gaussian
      return(FALSE)
    }
  }

  # Otherwise, if it inherits from 'lm' (and NOT from 'glm'),
  # this should be a plain linear model
  # since we already confirmed the class of the model
  TRUE
}
