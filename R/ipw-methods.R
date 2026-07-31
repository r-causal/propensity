#' @description
#' The `multinom` method estimates causal effects for a categorical exposure
#' from a fitted [nnet::multinom()] propensity score model and a weighted
#' outcome model. Standard errors are computed by M-estimation; the
#' linearization method is not available for categorical exposures. For a
#' K-level exposure, effects are reported for each non-reference level against
#' the reference (first) factor level, and the estimates table gains a
#' `comparison` column identifying the contrast.
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
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization"),
  .focal_level = NULL
) {
  se_method <- rlang::arg_match(se_method)
  assert_class(outcome_mod, c("glm", "lm"))

  # ps_link overrides the link of a binomial glm propensity model on the binary
  # path; a multinomial propensity model has no such link, so reject a non-NULL
  # argument rather than silently ignoring it.
  if (!is.null(ps_link)) {
    abort(
      c(
        "{.fun ipw} does not accept {.arg ps_link} for a multinomial propensity \\
        score model.",
        x = "A multinomial propensity score model has no link for \\
        {.arg ps_link} to override.",
        i = "Omit {.arg ps_link}; it applies only to a binomial glm propensity \\
        score model."
      ),
      error_class = "propensity_ipw_link_error"
    )
  }

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
    .focal_level = .focal_level
  )
  fit <- ipw_mestimation(spec, conf_level = conf_level)

  new_ipw(
    estimand = spec$estimand,
    wt_mod = wt_mod,
    outcome_mod = outcome_mod,
    estimates = fit$estimates,
    se_method = "mestimation",
    fit = fit$fit
  )
}

#' @description
#' The `lm` method estimates the causal dose-response effect for a continuous
#' exposure from a fitted [stats::lm()] (or gaussian-family [stats::glm()])
#' propensity score model of the exposure and a weighted marginal structural
#' outcome model. The only supported estimand is `"ate"`. Standard errors are
#' computed by M-estimation; the linearization method is not available for
#' continuous exposures. The marginal structural model must contain exactly one
#' exposure term, and the reported effect is that single coefficient: `"slope"`
#' for an identity-link outcome, `"log(or)"` for a logit-link outcome, and
#' `"log(rr)"` for a log-link outcome. The estimates table keeps the eight-column
#' contract with no comparison column.
#'
#' @name ipw-methods
#' @exportS3Method causalgenerics::ipw lm
ipw.lm <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization")
) {
  se_method <- rlang::arg_match(se_method)
  assert_class(outcome_mod, c("glm", "lm"))

  ipw_continuous_estimate(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    ps_link = ps_link,
    conf_level = conf_level,
    se_method = se_method
  )
}
