# M-estimation engine for the binary ipw() path. Three targets, all unexported:
# ipw_spec_binary() builds the design-contract spec from a fitted propensity
# score model and a fitted weighted outcome model, reusing the same extraction
# that ipw() uses; ipw_check_weight_consistency() is a preflight that recomputes
# the weights at the seeded init and compares them to the weights that fit the
# outcome model; ipw_mestimation() drives the deli fit from the root-valued init
# and returns the estimates table alongside the fit. The psi builder and theta
# layout live in R/ipw-psi.R.

# ---- spec constructor -------------------------------------------------------

ipw_spec_binary <- function(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL
) {
  assert_class(ps_mod, "glm")
  assert_class(outcome_mod, c("glm", "lm"))

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  if (is.null(.data)) {
    exposure <- fmla_extract_left_vctr(ps_mod)
    outcome <- fmla_extract_left_vctr(outcome_mod)
    mm_data <- model.frame(outcome_mod)
  } else {
    assert_columns_exist(.data, c(exposure_name, outcome_name))
    exposure <- .data[[exposure_name]]
    outcome <- .data[[outcome_name]]
    mm_data <- .data
  }

  if (is.null(ps_link)) {
    ps_link <- ps_mod$family$link
  }

  wts <- extract_weights(outcome_mod)
  estimand <- check_estimand(wts, estimand)

  exposure_values <- sort(unique(exposure))
  z <- as.double(exposure == exposure_values[[2]])

  out_terms <- stats::delete.response(stats::terms(outcome_mod))
  counterfactual_mm <- function(value) {
    d <- mm_data
    d[[exposure_name]] <- value
    model.matrix(out_terms, data = d)
  }
  x1 <- counterfactual_mm(exposure_values[[2]])
  x0 <- counterfactual_mm(exposure_values[[1]])

  if (is_linear_regression(outcome_mod)) {
    family <- "gaussian"
    out_link <- "identity"
    contrasts <- "diff"
  } else {
    family <- "binomial"
    out_link <- outcome_mod$family$link
    contrasts <- c("rd", "log(rr)", "log(or)")
  }

  list(
    exposure_type = "binary",
    estimand = estimand,
    n = length(z),
    exposure = z,
    ps = list(
      X = model.matrix(ps_mod),
      link = ps_link,
      coefs = stats::coef(ps_mod),
      k = 2L
    ),
    stab = list(
      stabilized = is_stabilized(wts),
      score = stabilization_score(wts)
    ),
    outcome = list(
      X = model.matrix(outcome_mod),
      y = as.double(outcome),
      family = family,
      link = out_link,
      coefs = stats::coef(outcome_mod),
      X_counterfactual = list(X1 = x1, X0 = x0),
      weights = as.double(wts)
    ),
    contrasts = contrasts,
    focal_level = NULL,
    reference_level = NULL
  )
}

# ---- weight-consistency preflight -------------------------------------------

# Recompute the weights implied by the spec at its seeded init, mirroring the
# weight computation each psi builder performs at theta = init.
ipw_weights_at_init <- function(spec, layout) {
  weight_fn <- ipw_weight_fn(spec$exposure_type, spec$estimand)
  init <- layout$init
  idx <- layout$idx
  th_ps <- init[idx$ps]
  th_stab <- init[idx$stab]
  score <- spec$stab$score

  switch(
    spec$exposure_type,
    binary = {
      inv_ps <- ipw_inv_link(spec$ps$link)
      e <- inv_ps(as.vector(spec$ps$X %*% th_ps))
      stab_prob <- if (length(th_stab)) th_stab[[1]] else NULL
      weight_fn(e, spec$exposure, list(stab_prob = stab_prob, score = score))
    },
    categorical = {
      k <- spec$ps$k
      ps_mat <- ipw_categorical_ps(spec$ps$X, th_ps, k)
      focal_idx <- if (!is.null(spec$focal_level)) {
        match(spec$focal_level, names(spec$outcome$X_counterfactual))
      } else {
        NULL
      }
      stab_probs <- if (length(th_stab)) c(1 - sum(th_stab), th_stab) else NULL
      weight_fn(
        ps_mat,
        spec$exposure,
        list(focal_idx = focal_idx, stab_probs = stab_probs, score = score)
      )
    },
    continuous = {
      n_alpha <- ncol(spec$ps$X)
      alpha <- th_ps[seq_len(n_alpha)]
      sigma2_d <- th_ps[[n_alpha + 1]]
      fitted_ps <- as.vector(spec$ps$X %*% alpha)
      mu_a <- if (length(th_stab)) th_stab[[1]] else NULL
      sigma2_a <- if (length(th_stab)) th_stab[[2]] else NULL
      weight_fn(
        fitted_ps,
        spec$exposure,
        list(
          sigma2_d = sigma2_d,
          mu_a = mu_a,
          sigma2_a = sigma2_a,
          score = score,
          stabilized = spec$stab$stabilized
        )
      )
    }
  )
}

ipw_check_weight_consistency <- function(
  spec,
  observed_wts,
  call = rlang::caller_env()
) {
  layout <- ipw_theta_layout(spec)
  recomputed <- as.double(ipw_weights_at_init(spec, layout))
  observed <- as.double(observed_wts)

  consistent <- length(recomputed) == length(observed) &&
    isTRUE(all.equal(recomputed, observed, tolerance = 1e-6))

  if (!consistent) {
    abort(
      c(
        "The weights used to fit {.arg outcome_mod} are not consistent with \\
        the propensity score model and estimand.",
        i = "The {.val {spec$estimand}} weights recomputed from {.arg ps_mod} \\
        differ from the weights supplied to {.arg outcome_mod} (compared at \\
        relative tolerance 1e-6).",
        i = "Refit {.arg outcome_mod} with weights from this propensity score \\
        model and estimand."
      ),
      error_class = "propensity_ipw_weights_mismatch_error",
      call = call
    )
  }

  invisible(TRUE)
}

# ---- engine -----------------------------------------------------------------

ipw_mestimation <- function(
  spec,
  conf_level = 0.95,
  call = rlang::caller_env()
) {
  layout <- ipw_theta_layout(spec)
  psi <- build_ipw_psi(spec, layout)

  if (!is.null(spec$outcome$weights)) {
    ipw_check_weight_consistency(spec, spec$outcome$weights, call = call)
  }

  m <- deli::MEstimator(stacked_equations = psi, init = layout$init)
  m <- deli::estimate(m)

  estimates <- ipw_mestimation_estimates(spec, m, conf_level)

  list(estimates = estimates, fit = m)
}

# Build the estimates table from the solved contrast rows of theta, mirroring
# the column layout that calculate_estimates() produces on the linearization
# path so the two SE methods return the same shape.
ipw_mestimation_estimates <- function(spec, fit, conf_level) {
  co <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  names(se) <- names(co)

  effect <- spec$contrasts
  estimate <- unname(co[effect])
  std.err <- unname(se[effect])
  z <- estimate / std.err

  z_val <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci.lower <- estimate - z_val * std.err
  ci.upper <- estimate + z_val * std.err
  p.value <- 2 * (1 - stats::pnorm(abs(z)))

  data.frame(
    effect = effect,
    estimate = estimate,
    std.err = std.err,
    z = z,
    ci.lower = ci.lower,
    ci.upper = ci.upper,
    conf.level = conf_level,
    p.value = p.value
  )
}
