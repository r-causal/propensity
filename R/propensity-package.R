#' @title propensity: A Toolkit for Calculating and Working with Propensity Scores
#'
#' @description
#' propensity provides tools for propensity score analysis in causal inference.
#' Calculate propensity score weights for a variety of causal estimands, handle
#' extreme propensity scores through trimming, truncation, and calibration, and
#' estimate causal effects with inverse probability weighting. The package
#' supports binary, categorical, and continuous exposures.
#'
#' @section Weight functions:
#' Calculate propensity score weights for different causal estimands:
#' * [wt_ate()]: Average treatment effect (ATE) weights
#' * [wt_att()]: Average treatment effect on the treated (ATT) weights
#' * [wt_atu()]: Average treatment effect on the untreated (ATU) weights
#'   (`wt_atc()` is an alias)
#' * [wt_atm()]: Average treatment effect for the evenly matchable (ATM) weights
#' * [wt_ato()]: Average treatment effect for the overlap population (ATO) weights
#' * [wt_entropy()]: Entropy balancing weights
#' * [wt_cens()]: Censoring weights
#'
#' @section Propensity score modifications:
#' Handle extreme propensity scores before calculating weights:
#' * [ps_trim()]: Trim observations with extreme propensity scores
#' * [ps_trunc()]: Truncate (winsorize) extreme propensity scores
#' * [ps_calibrate()]: Calibrate propensity scores to improve balance
#' * [ps_refit()]: Re-estimate the propensity score model after trimming
#'
#' @section Estimation:
#' * [ipw()]: Inverse probability weighted estimator with variance estimation
#'   that accounts for propensity score estimation uncertainty
#'
#' @section PSW class:
#' The [psw()] class represents propensity score weights with metadata about
#' the estimand and modifications applied:
#' * [psw()], [as_psw()], [is_psw()]: Create and test propensity score weights
#' * [estimand()]: Query the causal estimand
#' * [is_stabilized()]: Check if weights are stabilized
#'
#' @seealso
#' * `vignette("propensity")` for a getting started guide
#' * The [package website](https://r-causal.github.io/propensity/) for
#'   full documentation
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats median
#' @importFrom stats plogis
#' @importFrom stats qlogis
#' @importFrom stats quantile
#' @importFrom stats uniroot
## usethis namespace: end
NULL
