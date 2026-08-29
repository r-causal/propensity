#' @description
#' The `multinom` method estimates causal effects for a categorical exposure
#' from a fitted [nnet::multinom()] propensity score model and a weighted
#' outcome model. Standard errors are computed by M-estimation. Neither
#' `se_method = "linearization"` nor `se_method = "bootstrap"` is available for
#' a categorical exposure, and each is refused with an error of class
#' `propensity_method_error`: the stacked system has a sandwich for every fit
#' this exposure type accepts, and the resampling method is written for the
#' continuous fits it has no equation for. `boot_reps` and `boot_seed` describe
#' that resampling, so they are refused here as well, with class
#' `propensity_unsupported_arg_error` under the methods that resample nothing
#' and with the method refusal above when named alongside
#' `se_method = "bootstrap"`.
#'
#' For a K-level exposure, the counterfactual mean at each level is reported
#' first, under the effect label `"mean"`, and then the effects for each
#' non-reference level against the reference (first) factor level. The
#' `contrast` column names the level each mean row belongs to and the pair of
#' levels each effect row compares, as it does for a binary exposure.
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
  se_method = c("mestimation", "linearization", "bootstrap"),
  boot_reps = 500L,
  boot_seed = NULL,
  .focal_level = NULL,
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)
  check_ipw_boot_args(se_method, !missing(boot_reps), !missing(boot_seed))
  effects <- rlang::arg_match(effects)
  assert_class(outcome_mod, c("glm", "lm"))

  # ps_link overrides the link of a binomial glm propensity model on the binary
  # path; a multinomial propensity model has no such link, so reject a non-NULL
  # argument rather than silently ignoring it.
  check_ipw_ps_link_absent(ps_link, "multinomial")

  if (identical(se_method, "linearization")) {
    abort(
      c(
        "{.fun ipw} does not support {.val linearization} standard errors for \\
        categorical exposures.",
        i = "Use {.code se_method = \"mestimation\"} for a categorical \\
        exposure."
      ),
      error_class = "propensity_method_error"
    )
  }

  # The resampling method is written for the continuous route alone. A
  # categorical exposure has a sandwich for every fit it accepts, so nothing
  # here needs a method that repeats the fit, and routing one into the loop that
  # repeats it would rebuild a density ratio these weights never were.
  if (identical(se_method, "bootstrap")) {
    abort_ipw_boot_exposure("categorical")
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
#' outcome model. A [MASS::rlm()] fit with the Huber psi reaches this method by
#' inheritance and is stacked at its own score. A gaussian model is read through
#' its link, which may be identity or log; the other gaussian links, and every
#' class the stacked system has no score for, are refused, as are the two fits
#' that have no closed-form standard error, an [mgcv::gam()] propensity score
#' model and weights built with a `"kernel"` density. See **Continuous
#' propensity score models** in [ipw()] for what each refusal says and what to
#' do about it. The only supported estimand is `"ate"`.
#'
#' Standard errors are computed by M-estimation by default; the linearization
#' method is not available for continuous exposures. `se_method = "bootstrap"`
#' reports the spread of the whole fit repeated on `boot_reps` resamples of the
#' rows of `.data`, which it therefore requires. It is the method for the two
#' fits above, and it accepts every fit M-estimation accepts as well. The
#' propensity score model and an [mgcv::gam()] alike are refit on each
#' resample, the weights are rebuilt from the record the supplied weights carry,
#' and the marginal structural model is refit with them; the point estimate is
#' the one the supplied models already produced. See **Resampled standard
#' errors** in [ipw()].
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
  se_method = c("mestimation", "linearization", "bootstrap"),
  boot_reps = 500L,
  boot_seed = NULL,
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)
  check_ipw_boot_args(se_method, !missing(boot_reps), !missing(boot_seed))

  # Read before `arg_match()` resolves the default, for the reason `ipw.glm()`
  # reads it, and reading a forwarded set of readings as naming nothing for the
  # reason it does: the reading a dose-basis fit records is the same either way,
  # and what the caller named settles whether it is announced.
  effects_named <- !missing(effects) && length(effects) == 1L
  effects <- rlang::arg_match(effects)
  assert_class(outcome_mod, c("glm", "lm"))

  check_ipw_by_method(.by, se_method)
  check_ipw_by_exposure(.by)

  if (identical(se_method, "bootstrap")) {
    return(ipw_continuous_bootstrap(
      wt_mod,
      outcome_mod,
      .data = .data,
      estimand = estimand,
      ps_link = ps_link,
      conf_level = conf_level,
      boot_reps = boot_reps,
      boot_seed = boot_seed,
      effects = if (effects_named) effects
    ))
  }

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
