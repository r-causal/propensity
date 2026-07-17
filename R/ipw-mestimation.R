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
  ps_link = NULL,
  call = rlang::caller_env()
) {
  assert_class(ps_mod, "glm")
  assert_class(outcome_mod, c("glm", "lm"))

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  if (is.null(.data)) {
    exposure <- fmla_extract_left_vctr(ps_mod)
    outcome <- fmla_extract_left_vctr(outcome_mod)
    mm_data <- model.frame(outcome_mod)
    check_exposure(mm_data, exposure_name, call = call)
  } else {
    assert_columns_exist(.data, c(exposure_name, outcome_name), call = call)
    exposure <- .data[[exposure_name]]
    outcome <- .data[[outcome_name]]
    mm_data <- .data
  }

  if (!identical(length(exposure), length(outcome))) {
    abort(
      c(
        "{.arg exposure} and {.arg outcome} must be the same length.",
        x = "{.arg exposure} is length {length(exposure)}",
        x = "{.arg outcome} is length {length(outcome)}"
      ),
      call = call
    )
  }

  if (is.null(ps_link)) {
    ps_link <- ps_mod$family$link
  }

  wts <- extract_weights(outcome_mod)
  estimand <- check_estimand(wts, estimand, call = call)

  exposure_values <- sort(unique(exposure))

  if (!isTRUE(length(exposure_values) == 2)) {
    abort(
      c(
        "{.fun ipw} currently only supports binary exposures.",
        x = "There are {length(exposure_values)} unique value{?s} of the \\
        exposure."
      ),
      call = call
    )
  }

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

ipw_spec_categorical <- function(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  .focal_level = NULL,
  call = rlang::caller_env()
) {
  assert_class(ps_mod, "multinom")
  assert_class(outcome_mod, c("glm", "lm"))

  # A multinom fit with case weights would need a weighted score in the stacked
  # system; the ee_mlogit block is unweighted, so the fitted coefficients would
  # not sit at the score root and the estimates would drift. multinom always
  # carries a length-n weights vector, unit for an unweighted fit.
  if (!is.null(ps_mod$weights) && !all(ps_mod$weights == 1)) {
    abort(
      c(
        "{.fun ipw} does not support a propensity score model fit with case \\
        weights.",
        x = "{.arg ps_mod} was fit with non-unit {.arg weights}.",
        i = "Refit {.arg ps_mod} without {.arg weights}."
      ),
      error_class = "propensity_ipw_ps_weights_error",
      call = call
    )
  }

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  # The propensity score design comes from ps_mod. nnet::multinom stores no model
  # frame (model = FALSE), so model.matrix(ps_mod) re-evaluates the fitting call
  # in the formula environment; if the data behind that call is no longer
  # reachable the extraction fails. When .data is supplied, rebuild the design
  # from the model terms against .data so no re-evaluation is needed; otherwise
  # wrap the extraction and direct the user to supply .data on failure.
  if (is.null(.data)) {
    ps_extract <- tryCatch(
      list(
        exposure = fmla_extract_left_vctr(ps_mod),
        ps_X = model.matrix(ps_mod)
      ),
      error = function(e) e
    )
    if (inherits(ps_extract, "error")) {
      cause <- conditionMessage(ps_extract)
      abort(
        c(
          "Can't reconstruct the data behind {.arg ps_mod}.",
          x = "{cause}",
          i = "Supply {.arg .data} with the exposure, outcome, and covariates."
        ),
        error_class = "propensity_ipw_data_error",
        call = call
      )
    }
    exposure <- ps_extract$exposure
    ps_X <- ps_extract$ps_X
    outcome <- fmla_extract_left_vctr(outcome_mod)
    mm_data <- model.frame(outcome_mod)
    check_exposure(mm_data, exposure_name, call = call)
  } else {
    assert_columns_exist(.data, c(exposure_name, outcome_name), call = call)
    exposure <- .data[[exposure_name]]
    outcome <- .data[[outcome_name]]
    mm_data <- .data
    ps_X <- model.matrix(
      stats::delete.response(stats::terms(ps_mod)),
      data = .data,
      xlev = ps_mod$xlevels
    )
  }

  if (!identical(length(exposure), length(outcome))) {
    abort(
      c(
        "{.arg exposure} and {.arg outcome} must be the same length.",
        x = "{.arg exposure} is length {length(exposure)}",
        x = "{.arg outcome} is length {length(outcome)}"
      ),
      call = call
    )
  }

  exposure <- as.factor(exposure)
  levs <- levels(exposure)
  k <- length(levs)

  # Reference-first indicator matrix in factor-level order, matching the column
  # order deli::ee_mlogit expects (reference level in the first column).
  z_ind <- vapply(
    levs,
    function(l) as.integer(exposure == l),
    integer(length(exposure))
  )

  ps_coefs <- as.vector(t(stats::coef(ps_mod)))

  wts <- extract_weights(outcome_mod)

  # Estimand resolution mirrors the binary path: reconcile the psw estimand
  # attribute with an explicit `estimand` argument via check_estimand(), which
  # asks for `estimand` when the weights carry no attribute and none is supplied.
  estimand <- check_estimand(wts, estimand, call = call)

  # Focal level resolution: an explicit argument overrides the focal_category
  # attribute the weights carry. The att and atu estimands require a focal
  # level; without one there is nothing to contrast the reference against.
  focal_level <- if (is.null(.focal_level)) {
    attr(wts, "focal_category")
  } else {
    .focal_level
  }
  if (!is.null(focal_level) && !focal_level %in% levs) {
    abort(
      c(
        "{.arg .focal_level} must be one of the levels of \\
        {.arg {exposure_name}}.",
        x = "{.val {focal_level}} is not an exposure level.",
        i = "Available levels: {.val {levs}}."
      ),
      error_class = "propensity_focal_level_error",
      call = call
    )
  }
  if (estimand %in% c("att", "atu") && is.null(focal_level)) {
    abort(
      c(
        "A focal level is required for the {.val {estimand}} estimand with a \\
        categorical exposure.",
        i = "Supply {.arg .focal_level}, or fit {.arg outcome_mod} with \\
        weights created with a focal level so it is recorded on the weights."
      ),
      error_class = "propensity_focal_required_error",
      call = call
    )
  }

  # One counterfactual design per level: set the exposure factor to that level,
  # preserving all levels, and rebuild the outcome model matrix.
  out_terms <- stats::delete.response(stats::terms(outcome_mod))
  x_cf <- lapply(levs, function(l) {
    d <- mm_data
    d[[exposure_name]] <- factor(l, levels = levs)
    model.matrix(out_terms, data = d)
  })
  names(x_cf) <- levs

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
    exposure_type = "categorical",
    estimand = estimand,
    n = length(exposure),
    exposure = z_ind,
    ps = list(
      X = ps_X,
      link = NULL,
      coefs = ps_coefs,
      k = k
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
      X_counterfactual = x_cf,
      weights = as.double(wts)
    ),
    contrasts = contrasts,
    focal_level = focal_level,
    reference_level = levs[[1]]
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

  estimates <- ipw_mestimation_estimates(spec, m, layout, conf_level)

  list(estimates = estimates, fit = m)
}

# Build the estimates table from the solved contrast rows of theta, mirroring
# the column layout that calculate_estimates() produces on the linearization
# path so the two SE methods return the same shape. The contrast rows are
# addressed by their theta positions, not by name: theta names are not unique
# across blocks (a ps or outcome covariate can share a contrast label such as
# "rd" or "diff"), so name subsetting could silently return the wrong row.
ipw_mestimation_estimates <- function(spec, fit, layout, conf_level) {
  co <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))

  idx <- layout$idx$contrast
  estimate <- unname(co[idx])
  std.err <- unname(se[idx])
  z <- estimate / std.err

  z_val <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci.lower <- estimate - z_val * std.err
  ci.upper <- estimate + z_val * std.err
  p.value <- 2 * (1 - stats::pnorm(abs(z)))

  out <- data.frame(
    effect = rep(spec$contrasts, times = ipw_n_comparisons(spec)),
    estimate = estimate,
    std.err = std.err,
    z = z,
    ci.lower = ci.lower,
    ci.upper = ci.upper,
    conf.level = conf_level,
    p.value = p.value
  )

  # A categorical exposure contributes one contrast per non-reference level, so
  # the table gains a comparison column, placed immediately after effect, that
  # names each contrast. Binary and continuous exposures keep the eight-column
  # contract with no comparison column.
  if (identical(spec$exposure_type, "categorical")) {
    nonref <- names(spec$outcome$X_counterfactual)[-1]
    comparison <- rep(
      paste(nonref, "vs", spec$reference_level),
      each = length(spec$contrasts)
    )
    out <- cbind(
      out["effect"],
      comparison = comparison,
      out[setdiff(names(out), "effect")]
    )
  }

  out
}

# Number of contrast blocks the estimates table holds: one per non-reference
# level for a categorical exposure, a single block otherwise.
ipw_n_comparisons <- function(spec) {
  if (identical(spec$exposure_type, "categorical")) {
    spec$ps$k - 1L
  } else {
    1L
  }
}
