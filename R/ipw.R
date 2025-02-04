#' Inverse Probability Weights for Causal Inference
#'
#' `ipw()` is a bring-your-own-model (BYOM) inverse probability weighted
#' estimator. `ipw()` accepts a propensity score model and a weighted outcome
#' model that you have already fit. The purpose of `ipw()` is to capture the
#' uncertainty inherent to this two-step process and calculate the correct
#' standard errors for each estimate. Currently, `ipw()` supports binary
#' exposures and either binary or continuous outcomes. For binary outcomes,
#' `ipw()` calculates the marginal risk difference, log risk ratio, and log odds
#' ratio. For continuous outcomes, `ipw()` only calculates the marginal
#' difference in means.
#'
#' @param ps_mod A fitted propensity score model of class [stats::glm()],
#'   typically a logistic regression with the exposure as the dependent
#'   variable.
#' @param outcome_mod A fitted, weighted outcome model of class [stats::glm()]
#'   or [stats::lm()], with the outcome as the dependent variable.
#' @param .df A data frame containing the exposure, outcome, and covariates. If
#'   `NULL`, `ipw()` will try to extract the data from `ps_mod` and
#'   `outcome_mod`.
#' @param estimand A character string specifying the causal estimand: `ate`,
#'   `att`, `ato`, or `atm`. If `NULL`, the function attempts to infer it from
#'   existing weights in `outcome_mod`, assuming they were calculated with
#'   [wt_ate()], [wt_att()], [wt_atm()], or [wt_ato()].
#' @param ps_link A character string specifying the link function for the
#'   propensity score model: `logit`, `probit`, or `cloglog`. Defaults to
#'   whatever link was used by `ps_mod`.
#' @param conf_level Numeric. Confidence level for intervals (default `0.95`).
#'
#' @details The function constructs inverse probability weights based on the
#'   chosen `estimand`, then uses these weights (or extracts them from
#'   `outcome_mod`) to compute effect measures:
#' - `rd`: Risk difference
#' - `log(rr)`: Log risk ratio
#' - `log(or)`: Log odds ratio
#'
#'   For a linear outcome model (using [stats::lm()] or [stats::glm()] with
#'   `family = gaussian()`), only the difference in means (`diff`) is returned.
#'
#' **Variance Estimation**
#'
#'   The variance is estimated via linearization, which provide variance
#'   estimates for IPW that correctly account for the uncertainty in estimation
#'   of the propensity scores. For more details on various types of propensity
#'   score weights and their corresponding variance estimators, see:
#'
#' - *Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the inverse
#'   probability-of-treatment weighting estimator: A tutorial for different
#'   types of propensity score weights. Statistics in Medicine. 2024; 43(13):
#'   2672-2694. doi: [10.1002/sim.10078](https://doi.org/10.1002/sim.10078)*
#'
#' @return An S3 object of class `ipw` containing:
#' - `estimand`: One of `"ate"`, `"att"`, `"ato"`, or `"atm"`.
#' - `ps_mod`: The fitted propensity score model.
#' - `outcome_mod`: The fitted outcome model.
#' - `estimates`: A data frame of point estimates, standard errors, z-statistics,
#'   confidence intervals, and p-values.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' # confounder
#' x1 <- rnorm(n)
#' # exposure
#' z  <- rbinom(n, 1, plogis(0.2 * x1))
#' # binary outcome
#' y  <- rbinom(n, 1, plogis(1 + 2*z + 0.5*x1))
#'
#' dat <- data.frame(x1, z, y)
#'
#' # fit a propensity score model (exposure ~ x1)
#' ps_mod <- glm(z ~ x1, data = dat, family = binomial())
#'
#' # calculate weights for ATE
#' ps <- predict(ps_mod, type = "response")
#' wts <- wt_ate(ps, z)
#'
#' # fit an outcome model (binary y ~ z) using IPW
#' outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)
#'
#' # run IPW
#' ipw_res <- ipw(ps_mod, outcome_mod)
#'
#' ipw_res
#'
#' # Convert to a data frame with exponentiated RR/OR
#' ipw_res_df <- as.data.frame(ipw_res, exponentiate = TRUE)
#' ipw_res_df
#'
#' @seealso [stats::glm()], [stats::lm()]
#'
#' @export
#' @importFrom stats dnorm family formula model.frame model.matrix model.weights
#' @importFrom stats pnorm predict printCoefmat qnorm var
ipw <- function(ps_mod, outcome_mod, .df = NULL, estimand = NULL, ps_link = NULL, conf_level = 0.95) {
  stopifnot(inherits(ps_mod, "glm"))
  stopifnot(inherits(outcome_mod, "glm") || inherits(outcome_mod, "lm"))

  weight_matrix <- model.matrix(ps_mod)
  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  if (is.null(.df)) {
    exposure <- fmla_extract_left_vctr(ps_mod)
    outcome <- fmla_extract_left_vctr(outcome_mod)
  } else {
    exposure <- .df[[exposure_name]]
    outcome <- .df[[outcome_name]]
  }

  ps <- predict(ps_mod, type = "response", newdata = .df)

  if (is.null(ps_link)) {
    ps_link <- ps_mod$family$link
  }

  stopifnot(identical(length(exposure), length(outcome)))

  # todo: allow user to specify existing weights
  # or automatically extract weights from outcome model if they exist
  wts <- extract_weights(outcome_mod)
  estimand <- check_estimand(wts, estimand)

  marginal_means <- estimate_marginal_means(
    outcome_mod = outcome_mod,
    wts = wts,
    exposure = exposure,
    exposure_name = exposure_name,
    .df = .df
  )

  uncorrected_lin_vars <- linearize_variables_for_wts(
    exposure,
    outcome,
    wts,
    marginal_means
  )

  lin_vars <- linearize_variables_for_ps(
    exposure = exposure,
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

  structure(
    list(
      estimand = estimand,
      ps_mod = ps_mod,
      outcome_mod = outcome_mod,
      estimates = estimates
    ),
    class = "ipw"
  )
}

#' @export
print.ipw <- function(x, ...) {
  cat("Inverse Probability Weight Estimator\n")
  # todo: make this more adaptive, e.g. ATE without L2FU
  cat("Estimand:", toupper(x$estimand), "\n\n")

  cat("Propensity Score Model:\n")
  cat("  Call:", paste(deparse(x$ps_mod$call), collapse = "\n"), "\n")
  cat("\n")

  cat("Outcome Model:\n")
  cat("  Call:", paste(deparse(x$outcome_mod$call), collapse = "\n"), "\n")

  cat("\n")

  cat("Estimates:\n")
  estimates <- x$estimates[-1]
  rownames(estimates) <- x$estimates$effect
  printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 2:3, tst.ind = 4)

  invisible(x)
}


#' @param x an `ipw` object
#' @param exponentiate Logical. Should the log-RR and log-OR be exponentiated?
#' @param row.names,optional,... Passed to `as.data.frame()`.
#' @export
#' @rdname ipw
as.data.frame.ipw <- function(x, row.names = NULL, optional = NULL, exponentiate = FALSE, ...) {
  df <- as.data.frame(x$estimates, row.names = row.names, optional = optional, ...)

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

calculate_estimates <- function(lin_vars, marginal_means, n, linear_regression, conf_level) {
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

  # Influence for RR on the raw scale: (l1 / mu1 - l0 / mu0)
  rr_inf_raw <- lin_vars$l1 / marginal_means$mu1 - lin_vars$l0 / marginal_means$mu0
  # Then for log(RR), the influence is  (1 / RR) * rr_inf_raw
  log_rr_inf <- (1 / rr_raw_est) * rr_inf_raw

  log_rr_var <- var(log_rr_inf) / n
  log_rr_se <- sqrt(log_rr_var)

  log_rr_ci_lower <- log_rr_est - z_val * log_rr_se
  log_rr_ci_upper <- log_rr_est + z_val * log_rr_se

  rr_z <- log_rr_est / log_rr_se
  rr_p <- 2 * (1 - pnorm(abs(rr_z)))

  ### ODDS RATIO (log scale)
  # ---------------------------
  # OR = [mu1/(1 - mu1)] / [mu0/(1 - mu0)]

  or_raw_est <- (marginal_means$mu1 / (1 - marginal_means$mu1)) / (marginal_means$mu0 / (1 - marginal_means$mu0))
  log_or_est <- log(or_raw_est)

  # Influence for OR on the raw scale:
  #   l_or_raw = l1 / [mu1*(1 - mu1)] - l0 / [mu0*(1 - mu0)]
  or_inf_raw <- (lin_vars$l1 / (marginal_means$mu1 * (1 - marginal_means$mu1))) - (lin_vars$l0 / (marginal_means$mu0 * (1 - marginal_means$mu0)))
  # Then log(OR) influence =  (1 / OR) * or_inf_raw
  log_or_inf <- (1 / or_raw_est) * or_inf_raw

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
linearize_variables_for_ps <- function(exposure, outcome, wts, ps, estimand, weight_matrix, marginal_means, uncorrected_lin_vars, ps_link) {
  n <- length(exposure)
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

  l1 <- uncorrected_lin_vars$l1 + correct_for_ps(
    exposure = exposure,
    outcome = outcome,
    ps = ps,
    mu = marginal_means$mu1,
    n_group = marginal_means$n1,
    weight_matrix = weight_matrix,
    weight_derivatives = weight_derivatives,
    correction_mat = correction_mat,
    n = n
  )

  l0 <- uncorrected_lin_vars$l0 + correct_for_ps(
    exposure = exposure,
    exposure_actual = 1 - exposure,
    outcome = outcome,
    ps = ps,
    mu = marginal_means$mu0,
    n_group = marginal_means$n0,
    weight_matrix = weight_matrix,
    weight_derivatives = weight_derivatives,
    correction_mat = correction_mat,
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

correct_for_ps <- function(exposure, exposure_actual = exposure, outcome, ps, mu, n_group, weight_matrix, weight_derivatives, correction_mat, n) {
  drop(
    n / n_group *
      rbind(colSums(weight_derivatives * exposure_actual * (outcome - mu)) / n) %*%
      (correction_mat %*% t((exposure - ps) * weight_matrix))
  ) |> unname()
}

estimate_marginal_means <- function(outcome_mod, wts, exposure, exposure_name, .df = NULL) {
  # todo: this could be generalized with split() and lapply()
  if (is.null(.df)) {
    .df <- model.frame(outcome_mod)
    check_exposure(.df, exposure_name)
  }
  # todo: make this more flexible for different values and model specs
  # maybe can optionally accept a function for g-comp
  exposure_values <- sort(unique(exposure))
  stopifnot(isTRUE(length(exposure_values) == 2))
  .df_1 <- .df
  .df_0 <- .df
  .df_1[[exposure_name]] <- exposure_values[[2]]
  .df_0[[exposure_name]] <- exposure_values[[1]]

  n1 <- sum(wts[exposure == exposure_values[[2]]])
  mu1 <- mean(predict(outcome_mod, newdata = .df_1, type = "response"))

  n0 <- sum(wts[exposure == exposure_values[[1]]])
  mu0 <- mean(predict(outcome_mod, newdata = .df_0, type = "response"))

  list(
    # exposure = 1
    n1 = n1, mu1 = mu1,
    # exposure = 0
    n0 = n0, mu0 = mu0
  )
}

derive_weights <- function(exposure, ps, weight_matrix, ps_link = c("logit", "probit", "cloglog"), estimand = c("ate", "att", "ato", "atm")) {
  estimand <- match.arg(estimand)
  deriv_link_f <- derive_link(ps_link)

  if (estimand == "ate") {
    weights_vector <- ifelse(
      exposure == 1,
      -1 / (ps^2 * deriv_link_f(ps)),
      1 / ((1 - ps)^2 * deriv_link_f(ps))
    )
  } else if (estimand == "att") {
    weights_vector <- ifelse(
      exposure == 1,
      0,
      1 / ((1 - ps)^2 * deriv_link_f(ps))
    )
  } else if (estimand == "ato") {
    weights_vector <- ifelse(
      exposure == 1,
      -1 / (deriv_link_f(ps)),
      1 / (deriv_link_f(ps))
    )
  } else if (estimand == "atm") {
    weights_vector <- ifelse(
      exposure == 1,
      -1 / (ps^2 * deriv_link_f(ps)),
      1 / ((1 - ps)^2 * deriv_link_f(ps))
    )

    weights_vector[exposure == 1 & ps < 0.5] <- 0
    weights_vector[exposure == 0 & ps > 0.5] <- 0
  }

  weight_matrix * weights_vector
}

derive_link <- function(ps_link = c("logit", "probit", "cloglog")) {
  ps_link <- match.arg(ps_link)
  switch(ps_link,
         logit = function(x) 1 / (x * (1 - x)),
         probit = function(x) 1 / (dnorm(qnorm(x))),
         cloglog = function(x) -1 / ((1 - x) * log(1 - x))
  )
}

fmla_extract_left_vctr <- function(mod) {
  .data <- mod |>
    model.frame()

  .data[[1]]
}

fmla_extract_left_chr <- function(mod) {
  as.character(formula(mod)[[2]])
}

extract_weights <- function(.mod) {
  .mod |>
    model.frame() |>
    model.weights()
}

check_estimand <- function(wts, estimand) {
  if (is_causal_wt(wts)) {
    estimand_from_weights <- estimand(wts)
  } else {
    estimand_from_weights <- NULL
  }

  if (!is.null(estimand_from_weights) && !is.null(estimand)) {
    same_estimand <- identical(estimand_from_weights, estimand)
    if (!same_estimand) {
      stop(
        "Estimand in weights different from `estimand`: ",
        estimand_from_weights,
        " vs. ",
        estimand
      )
    } else {
      return(estimand)
    }
  }

  if (is.null(estimand_from_weights) && is.null(estimand)) {
    stop("Can't determine estimand from weights. Please specify `estimand`.")
  }

  if (!is.null(estimand_from_weights)) {
    return(estimand_from_weights)
  } else {
    return(estimand)
  }
}

check_exposure <- function(.df, .exposure_name) {
  if (!(.exposure_name %in% names(.df))) {
    stop(
      .exposure_name,
      " not found in `model.frame(outcome_mod)`. ",
      "The outcome model may have transformations in the formula. ",
      "Please specify `.df`")
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
