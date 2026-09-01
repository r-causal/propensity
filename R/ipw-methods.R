#' @description
#' The `multinom` method estimates causal effects for a categorical exposure
#' from a fitted [nnet::multinom()] propensity score model and a weighted
#' outcome model. Standard errors are computed by M-estimation. Neither
#' `se_method = "linearization"` nor `se_method = "robust"` is available for a
#' categorical exposure, and both are refused with an error of class
#' `propensity_method_error`: the stacked system has a sandwich for every fit
#' this exposure type accepts, and the diagnostic sandwich describes two fitted
#' cells that a categorical result does not report.
#'
#' For a K-level exposure, the counterfactual mean at each level is reported
#' first, under the effect label `"mean"`, and then the effects for each
#' non-reference level against the reference (first) factor level. The
#' `contrast` column names the level each mean row belongs to and the pair of
#' levels each effect row compares, as it does for a binary exposure.
#'
#' Two fits are refused with an error of class
#' `propensity_model_family_error`. A [nnet::multinom()] fit to two levels
#' reports the single probability of a binary exposure rather than one
#' probability per level; fit such an exposure with [stats::glm()] and
#' `family = binomial()` and pass it to the `glm` method instead. A fit to a
#' matrix response records no levels at all, since [nnet::multinom()] reads a
#' matrix as counts; refit it with the exposure factor on the left-hand side.
#'
#' @param .focal_level For the categorical (`multinom`) method with the `att`
#'   or `atu` estimand, the focal exposure level. If `NULL`, it is taken from
#'   the `focal_category` attribute of the weights in `outcome_mod`. An explicit
#'   value overrides the attribute.
#'
#' @name ipw-methods
#' @exportS3Method causalgenerics::ipw multinom
ipw.multinom <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  .focal_level = NULL,
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)
  effects <- rlang::arg_match(effects)
  assert_class(outcome_mod, c("glm", "lm"))

  # ps_link overrides the link of a binomial glm propensity model on the binary
  # path; a multinomial propensity model has no such link, so reject a non-NULL
  # argument rather than silently ignoring it.
  check_ipw_ps_link_absent(ps_link, "multinomial")

  # The same two refusals the weight path makes of a `multinom`, in the terms
  # this method's arguments are named in. A matrix response is read as counts,
  # so the fit records no levels and the exposure would be read as the deparsed
  # response instead; a two-level fit is a binary propensity score model.
  check_multinom_response(wt_mod, arg = "wt_mod")
  check_ipw_multinom_levels(wt_mod)

  if (!identical(se_method, "mestimation")) {
    abort(
      c(
        "{.fun ipw} does not support {.val {se_method}} standard errors for \\
        categorical exposures.",
        i = "Use {.code se_method = \"mestimation\"} for a categorical \\
        exposure."
      ),
      error_class = "propensity_method_error"
    )
  }

  # Guards first, on the weights that fit the outcome model, mirroring ipw.glm:
  # the psw attributes carried by the weights detect a modified propensity score
  # before any estimand parsing that would otherwise fail obliquely.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts)

  spec <- ipw_spec_categorical(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    .focal_level = .focal_level,
    .by = .by
  )
  fit <- ipw_mestimation(spec, conf_level = conf_level)
  components <- ipw_component_models(wt_mod, outcome_mod, fit)

  new_ipw(
    estimand = spec$estimand,
    wt_mod = components$wt_mod,
    outcome_mod = components$outcome_mod,
    estimates = fit$estimates,
    se_method = "mestimation",
    fit = fit$fit,
    effects = effects
  )
}

#' @description
#' The `lm` method estimates the causal dose-response effect for a continuous
#' exposure from a fitted [stats::lm()] (or gaussian-family [stats::glm()])
#' propensity score model of the exposure and a weighted marginal structural
#' outcome model. A [MASS::rlm()] fit with one of the psi functions MASS
#' supplies and an [mgcv::gam()] of a gaussian family reach this method by
#' inheritance and are stacked at their own scores, the additive fit with its
#' smoothing parameters held at the values it reports. A gaussian model is read
#' through its link, which may be identity or log; the other gaussian links, and
#' every class the stacked system has no score for, are refused, as are weights
#' built with a `"kernel"` density, which have no closed-form standard error.
#' See **Continuous
#' propensity score models** in [ipw()] for what each refusal says, what to do
#' about it, and what the additive route conditions on. The only supported
#' estimand is `"ate"`.
#'
#' Standard errors are computed by M-estimation; neither the linearization
#' method nor the robust diagnostic is available for continuous exposures. The
#' package offers no resampling method,
#' so a fit the stacked system cannot write has no standard error here, and the
#' refusal points at a bootstrap of the whole pipeline written by hand. See
#' **Continuous propensity score models** in [ipw()].
#'
#' Every term of the marginal structural model that reads the exposure must read
#' the exposure and nothing else, and the reported effects are the coefficients
#' those terms contribute, labeled `"log(or)"` for a logit-link outcome and
#' `"log(rr)"` for a log-link one. A model with one
#' exposure coefficient keeps the eight-column contract with no contrast column
#' and labels its row `"slope"` at an identity link; a dose-response curve such
#' as `y ~ A + I(A^2)` reports one row per coefficient under `"coef"`, gaining a
#' contrast column that names them.
#'
#' A curve is also the one shape of fit that has a single reading. An exposure
#' entering the outcome model through more than one design column, whether
#' written as `y ~ A + I(A^2)` or as a basis such as `poly(A, 2)`,
#' `splines::ns(A, 3)`, or `splines::bs(A, 3)`, has no coefficient that is the
#' effect of the exposure, because the dose response has a different slope at
#' every dose. `ipw()` records the conditional reading for such a fit and
#' reports it once as a message. Asking `ipw()` itself for
#' `effects = "marginal"` is refused with an error of class
#' `propensity_ipw_effects_error`; the result declares the conditional reading
#' as the only one it supports, so [causalgenerics::as_marginal()],
#' [stats::coef()], [stats::vcov()], [stats::confint()], [generics::tidy()],
#' and `as.data.frame()` refuse it from the result class instead, with an error
#' of class `causalgenerics_unsupported_reading_marginal`. Marginalizing the
#' curve over the observed doses is a separate estimand that this package does
#' not compute; use the \pkg{marginaleffects} package on the conditional
#' result: `avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and
#' `avg_predictions()` for causal dose-response functions. A basis
#' term reading a covariate rather than the exposure is a covariate term however
#' many columns it expands to, so it contributes no row and leaves the reading
#' marginal.
#'
#' A model frame records the term rather than the variables inside it, so the
#' frame of a fit whose exposure enters through `poly()` or a spline holds no
#' exposure column and `.data` must supply one. The propensity score design is
#' then rebuilt through the model's terms object, whose `predvars` attribute
#' records the basis the model was fit with, so the rebuilt columns are the
#' fitted ones.
#'
#' @name ipw-methods
#' @exportS3Method causalgenerics::ipw lm
ipw.lm <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)

  # Read before `arg_match()` resolves the default, for the reason `ipw.glm()`
  # reads it, and reading a forwarded set of readings as naming nothing for the
  # reason it does: the reading a dose-basis fit records is the same either way,
  # and what the caller named settles whether it is announced.
  effects_named <- !missing(effects) && length(effects) == 1L
  effects <- rlang::arg_match(effects)
  assert_class(outcome_mod, c("glm", "lm"))

  check_ipw_by_method(.by, se_method)
  check_ipw_by_exposure(.by)

  ipw_continuous_estimate(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    ps_link = ps_link,
    conf_level = conf_level,
    se_method = se_method,
    effects = if (effects_named) effects
  )
}
