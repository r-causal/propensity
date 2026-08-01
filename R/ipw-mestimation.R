# M-estimation engine for the binary ipw() path. Three targets, all unexported:
# ipw_spec_binary() builds the design-contract spec from a fitted propensity
# score model and a fitted weighted outcome model, reusing the same extraction
# that ipw() uses; ipw_check_weight_consistency() is a preflight that recomputes
# the weights at the seeded init and compares them to the weights that fit the
# outcome model; ipw_mestimation() drives the deli fit from the root-valued init
# and returns the estimates table alongside the fit. The psi builder and theta
# layout live in R/ipw-psi.R.

# ---- spec constructor -------------------------------------------------------

# Reject counterfactual designs the outcome model cannot use. The designs are
# built by setting the exposure to each level in turn, so the coding of the
# exposure decides whether a per-level marginal mean is estimable at all. Two
# ways it is not, checked in this order:
#
# A design that is identically zero. A numeric no-intercept coding such as
# `y ~ z - 1` leaves the design at the zero-coded level with no nonzero entry at
# all. The marginal mean there is then inv_link(0) for every unit, a constant
# fixed by the outcome link rather than a quantity estimated from the data (0.5
# for a logit or probit outcome link, 0 for a linear model), and every contrast
# formed against it is wrong with nothing signaled.
#
# A pair of designs that are elementwise identical. An indicator for a single
# level of a multi-level exposure, such as `y ~ as.numeric(z == "c") + x1`,
# leaves the designs at every other level agreeing in every column: the
# indicator is zero at all of them and the covariates are untouched. The model
# predicts one and the same outcome at both levels, so their contrast is zero to
# rounding error, its standard error degenerates to NaN, and the number reported
# says nothing about the exposure. The zero case is checked first because two
# all-zero designs are also identical to each other and the zero message is the
# more specific diagnosis.
#
# Neither condition is about the presence of an intercept. A saturated factor
# coding such as `y ~ 0 + zf` carries no intercept yet its dummy columns sum to
# one at every level and separate every pair, and `y ~ z + x1 - 1` still
# estimates the marginal mean from the covariate. Both are honest g-computation
# on the model as specified and are left alone.
check_ipw_counterfactual_designs <- function(
  designs,
  exposure_levels,
  exposure_name,
  call = rlang::caller_env()
) {
  degenerate <- vapply(designs, function(x) isTRUE(all(x == 0)), logical(1))

  if (any(degenerate)) {
    bad <- exposure_levels[degenerate]
    abort(
      c(
        "{.arg outcome_mod} must be able to represent the outcome at every \\
        exposure level.",
        x = "Setting {.arg {exposure_name}} to {.val {bad}} leaves the \\
        counterfactual design{?s} identically zero, which pins the marginal \\
        mean{?s} there to the outcome link's zero point instead of estimating \\
        {?it/them}.",
        i = "Include an intercept in {.arg outcome_mod}, or code the exposure \\
        as a factor, whose no-intercept coding is saturated and represents \\
        every level."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  # The designs are deterministic constructions from one terms object, so equal
  # entries are bit-for-bit equal and no tolerance is involved. Only the first
  # offending pair is reported: any fix is a recoding of the exposure, after
  # which the guard runs again over the new designs.
  bad <- first_identical_design_pair(designs)

  if (!is.null(bad)) {
    bad <- exposure_levels[bad]
    abort(
      c(
        "{.arg outcome_mod} must be able to distinguish every pair of \\
        exposure levels.",
        x = "Setting {.arg {exposure_name}} to {.val {bad}} produces identical \\
        counterfactual designs, so the model predicts the same outcome at both \\
        levels and the contrast between them is degenerate.",
        i = "Code the exposure so {.arg outcome_mod} separates every level, \\
        for example as a factor rather than an indicator for a single level."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Index pair of the first two designs that agree in every entry, or NULL if all
# of them differ. Pairs are visited in level order, so the pair reported is
# stable across calls.
first_identical_design_pair <- function(designs) {
  for (j in seq_along(designs)[-1]) {
    for (i in seq_len(j - 1)) {
      same_dim <- identical(dim(designs[[i]]), dim(designs[[j]]))
      if (same_dim && isTRUE(all(designs[[i]] == designs[[j]]))) {
        return(c(i, j))
      }
    }
  }

  NULL
}

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
  check_ipw_ps_response(ps_mod, call = call)
  check_ipw_offset(outcome_mod, call = call)
  check_ipw_outcome_family(outcome_mod, call = call)

  # A binary propensity model fit with case weights would need a weighted score
  # in the stacked system; the ee_glm ps block is unweighted, so the fitted
  # coefficients would not sit at the score root and the estimates would drift.
  # A glm records prior weights in prior.weights (all one when unweighted).
  if (!is.null(ps_mod$prior.weights) && !all(ps_mod$prior.weights == 1)) {
    abort(
      c(
        "{.fun ipw} does not support a propensity score model fit with case \\
        weights.",
        x = "The propensity score model was fit with non-unit {.arg weights}.",
        i = "Refit the propensity score model without {.arg weights}."
      ),
      error_class = "propensity_ipw_ps_weights_error",
      call = call
    )
  }

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  check_ipw_outcome_exposure(outcome_mod, exposure_name, call = call)

  extracted <- ipw_extract_ps_design(
    ps_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    xlev = ps_mod$xlevels,
    call = call
  )
  exposure <- extracted$exposure
  outcome <- extracted$outcome
  ps_X <- extracted$ps_X
  mm_data <- extracted$mm_data

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

  # Restrict the resolved link, whether it came from `ps_link` or the model's
  # family, to the three links the binary propensity score is documented and
  # derived for. Without this the link reaches ipw_inv_link(), which serves
  # outcome-link reconstruction as well and so also accepts identity and log,
  # leaving an unsupported propensity model to be estimated without complaint.
  # Membership is all that is checked: a supported link whose weights disagree
  # with the fitted model is the weight-consistency preflight's business.
  ps_link_supported <- c("logit", "probit", "cloglog")
  if (!isTRUE(ps_link %in% ps_link_supported)) {
    abort(
      c(
        "{.fun ipw} does not support the {.val {ps_link}} link for a binary \\
        propensity score model.",
        i = "Supported links: {.val {ps_link_supported}}.",
        i = "Refit {.arg wt_mod} with a supported link, or set {.arg ps_link} \\
        to one of them."
      ),
      error_class = "propensity_ipw_link_error",
      call = call
    )
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

  z <- ipw_recode_binary_exposure(exposure)

  # The rebuilt design multiplies the outcome model's coefficients positionally,
  # so it has to reproduce the coding those coefficients were fit under.
  # `contrasts.arg` carries it: setting the exposure column drops any contrasts
  # attribute it had, and a factor-free fit records NULL, which is the default.
  out_terms <- stats::delete.response(stats::terms(outcome_mod))
  counterfactual_mm <- function(value) {
    d <- mm_data
    d[[exposure_name]] <- value
    model.matrix(out_terms, data = d, contrasts.arg = outcome_mod$contrasts)
  }
  x1 <- counterfactual_mm(exposure_values[[2]])
  x0 <- counterfactual_mm(exposure_values[[1]])

  check_ipw_counterfactual_designs(
    list(x0, x1),
    as.character(exposure_values),
    exposure_name,
    call = call
  )

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
      X = ps_X,
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
      y = ipw_outcome_numeric(outcome),
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

# Extract the propensity score design and the exposure and outcome vectors,
# rebuilding from .data when it is supplied. A propensity model can lose the data
# behind its fitting call (nnet::multinom stores no model frame; an lm or glm can
# be fit in an environment that is later gone), so model.matrix(ps_mod) is wrapped
# and the user is directed to supply .data on failure. When .data is supplied the
# design is rebuilt from the model terms so no re-evaluation is needed.
ipw_extract_ps_design <- function(
  ps_mod,
  outcome_mod,
  .data,
  exposure_name,
  outcome_name,
  xlev = NULL,
  call = rlang::caller_env()
) {
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
          "Can't reconstruct the data behind {.arg wt_mod}.",
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

    # Everything downstream is sized to `.data` while the weights come from the
    # outcome model fit, so a row count that disagrees leaves the two to be
    # recycled against each other. Left alone it surfaces much later, as weights
    # that fail their consistency check, which reports a problem with how the
    # weights were built when the mistake was the data frame that was passed.
    #
    # Compare against the fitted model frame rather than `stats::nobs()`: the
    # model frame's row count is the length of the weights the outcome model was
    # fit with, which is what this reconciliation is really about, while `nobs()`
    # subtracts zero-weight observations and would reject a correct `.data` for a
    # fit that has any. The frame is already built before this point on every
    # path, by the weight extraction each `ipw()` method runs first, so asking
    # for it here cannot fail on its own.
    n_fitted <- nrow(stats::model.frame(outcome_mod))

    if (!identical(nrow(.data), n_fitted)) {
      abort(
        c(
          "{.arg .data} must have one row per observation the models were fit \\
          to.",
          x = "{.arg .data} has {nrow(.data)} rows.",
          x = "{.arg outcome_mod} was fit to {n_fitted} observations.",
          i = "Supply the data the models were fit to, or omit {.arg .data}."
        ),
        error_class = "propensity_ipw_data_error",
        call = call
      )
    }

    exposure <- .data[[exposure_name]]
    outcome <- .data[[outcome_name]]
    mm_data <- .data
    # Rebuilding a factor against `xlev` drops its contrasts attribute, so the
    # fit's own coding has to be supplied here or the design no longer matches
    # the coefficients it is multiplied by. NULL, for a factor-free fit, is the
    # default. Carrying the coding in `contrasts.arg` also makes the attribute
    # on the column redundant, so drop it first: `model.frame()` would otherwise
    # warn that re-leveling lost it, which says nothing about this rebuild. The
    # design is the same either way.
    ps_X <- model.matrix(
      stats::delete.response(stats::terms(ps_mod)),
      data = drop_contrasts_attrs(.data, names(ps_mod$contrasts)),
      xlev = xlev,
      contrasts.arg = ps_mod$contrasts
    )
  }

  list(exposure = exposure, outcome = outcome, ps_X = ps_X, mm_data = mm_data)
}

# Clear the contrasts attribute from the named columns of a data frame. Only for
# rebuilds that pass the coding separately through `contrasts.arg`, where the
# attribute duplicates it. Columns not present are skipped, so a fit whose
# covariates are absent from the data still reaches the error model.matrix
# raises rather than one from here.
drop_contrasts_attrs <- function(data, cols) {
  for (nm in intersect(cols, names(data))) {
    attr(data[[nm]], "contrasts") <- NULL
  }

  data
}

# Resolve the exposure to a factor whose levels are ordered the way the fitted
# multinomial model orders them. nnet::multinom records its training levels in
# `lev` and lays its coefficient rows out in that order, so that order governs
# the indicator matrix, the coefficient naming, the counterfactual designs, and
# the reference level. The exposure column reaching this point carries no such
# guarantee: a character column orders alphabetically and a factor can be
# releveled, either of which would silently disagree with the coefficients. A
# model that carries no `lev` leaves nothing to resolve against, so the column's
# own ordering is used, as it was before.
ipw_categorical_exposure_factor <- function(
  exposure,
  lev,
  exposure_name,
  call = rlang::caller_env()
) {
  if (is.null(lev)) {
    return(as.factor(exposure))
  }

  values <- as.character(exposure)
  unknown <- unique(values[!is.na(values) & !values %in% lev])
  if (length(unknown) > 0) {
    abort(
      c(
        "{.arg {exposure_name}} holds values the propensity score model was \\
        not fit on.",
        x = "{.val {unknown}} {?is/are} not among the model's levels.",
        i = "The model was fit on {.val {lev}}.",
        i = "Supply data whose {.arg {exposure_name}} holds only those levels."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  factor(values, levels = lev)
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
  check_ipw_offset(outcome_mod, call = call)
  check_ipw_outcome_family(outcome_mod, call = call)

  # A focal level names a single exposure level. Reject a longer argument here,
  # before the membership check reaches an `&&` on a non-scalar and raises a raw
  # coercion error. A length-1 character or factor level, or NULL (ate has no
  # focal level), is left to the resolution and membership logic below.
  if (!is.null(.focal_level) && length(.focal_level) != 1) {
    abort(
      c(
        "{.arg .focal_level} must be a single exposure level.",
        x = "{.arg .focal_level} has length {length(.focal_level)}."
      ),
      error_class = "propensity_focal_level_error",
      call = call
    )
  }

  # A multinom fit with case weights would need a weighted score in the stacked
  # system; the ee_mlogit block is unweighted, so the fitted coefficients would
  # not sit at the score root and the estimates would drift. multinom always
  # carries a length-n weights vector, unit for an unweighted fit.
  if (!is.null(ps_mod$weights) && !all(ps_mod$weights == 1)) {
    abort(
      c(
        "{.fun ipw} does not support a propensity score model fit with case \\
        weights.",
        x = "The propensity score model was fit with non-unit {.arg weights}.",
        i = "Refit the propensity score model without {.arg weights}."
      ),
      error_class = "propensity_ipw_ps_weights_error",
      call = call
    )
  }

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  check_ipw_outcome_exposure(outcome_mod, exposure_name, call = call)
  check_ipw_outcome_exposure_class(outcome_mod, exposure_name, call = call)
  check_ipw_outcome_exposure_levels(
    outcome_mod,
    exposure_name,
    ps_mod$lev,
    call = call
  )

  extracted <- ipw_extract_ps_design(
    ps_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    xlev = ps_mod$xlevels,
    call = call
  )
  exposure <- extracted$exposure
  outcome <- extracted$outcome
  ps_X <- extracted$ps_X
  mm_data <- extracted$mm_data

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

  exposure <- ipw_categorical_exposure_factor(
    exposure,
    ps_mod$lev,
    exposure_name,
    call = call
  )
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
  # preserving all levels, and rebuild the outcome model matrix. The replacement
  # is a bare factor carrying no contrasts attribute, and the design multiplies
  # the outcome model's coefficients positionally, so the fit's own coding has
  # to be supplied through `contrasts.arg`. It keys by column name and overrides
  # the column, which is what lets an ordered fit's polynomial coding be
  # reproduced here from an unordered replacement.
  out_terms <- stats::delete.response(stats::terms(outcome_mod))
  x_cf <- lapply(levs, function(l) {
    d <- mm_data
    d[[exposure_name]] <- factor(l, levels = levs)
    model.matrix(out_terms, data = d, contrasts.arg = outcome_mod$contrasts)
  })
  names(x_cf) <- levs

  check_ipw_counterfactual_designs(x_cf, levs, exposure_name, call = call)

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
      y = ipw_outcome_numeric(outcome),
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

ipw_spec_continuous <- function(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  call = rlang::caller_env()
) {
  assert_class(ps_mod, c("glm", "lm"))
  assert_class(outcome_mod, c("glm", "lm"))
  check_ipw_offset(outcome_mod, call = call)
  check_ipw_outcome_family(outcome_mod, call = call)
  check_ipw_continuous_links(ps_mod, outcome_mod, call = call)

  # A propensity model fit with prior case weights would need a weighted score in
  # the stacked system; the ee_regression ps block is unweighted, so the fitted
  # coefficients would not sit at the score root and the estimates would drift.
  # A glm records prior weights in prior.weights (all one when unweighted); an lm
  # records them in weights (NULL when unweighted).
  ps_weights <- if (inherits(ps_mod, "glm")) {
    ps_mod$prior.weights
  } else {
    ps_mod$weights
  }
  if (!is.null(ps_weights) && !all(ps_weights == 1)) {
    abort(
      c(
        "{.fun ipw} does not support a propensity score model fit with case \\
        weights.",
        x = "The propensity score model was fit with non-unit {.arg weights}.",
        i = "Refit the propensity score model without {.arg weights}."
      ),
      error_class = "propensity_ipw_ps_weights_error",
      call = call
    )
  }

  exposure_name <- fmla_extract_left_chr(ps_mod)
  outcome_name <- fmla_extract_left_chr(outcome_mod)

  extracted <- ipw_extract_ps_design(
    ps_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    xlev = ps_mod$xlevels,
    call = call
  )
  exposure <- as.double(extracted$exposure)
  outcome <- extracted$outcome
  ps_X <- extracted$ps_X

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

  wts <- extract_weights(outcome_mod)

  # A continuous exposure supports only the ate estimand. This guard fires before
  # check_estimand() so an unsupported request errors with its own class rather
  # than as an estimand-versus-weights mismatch.
  if (!is.null(estimand) && !identical(estimand, "ate")) {
    abort(
      c(
        "{.fun ipw} only supports the {.val ate} estimand for a continuous \\
        exposure.",
        x = "The {.val {estimand}} estimand is not available for a continuous \\
        exposure.",
        i = "Use {.code estimand = \"ate\"} or omit {.arg estimand}."
      ),
      error_class = "propensity_ipw_estimand_error",
      call = call
    )
  }
  estimand <- check_estimand(wts, estimand, call = call)

  # Marginal structural model term detection: exactly one term of the outcome
  # model references the exposure, contributing exactly one coefficient. A model
  # with several exposure terms has no single reported effect.
  out_X <- model.matrix(outcome_mod)
  term_labels <- attr(stats::terms(outcome_mod), "term.labels")
  is_exposure_term <- vapply(
    term_labels,
    function(l) exposure_name %in% all.vars(str2lang(l)),
    logical(1)
  )
  exposure_cols <- which(attr(out_X, "assign") %in% which(is_exposure_term))
  if (length(exposure_cols) != 1) {
    abort(
      c(
        "{.fun ipw} requires a marginal structural model with a single \\
        exposure term for a continuous exposure.",
        x = "{.arg outcome_mod} contributes {length(exposure_cols)} \\
        coefficient{?s} for {.val {exposure_name}}.",
        i = "Read the full coefficient vector from the returned {.field fit} \\
        object for a model with more than one exposure term."
      ),
      error_class = "propensity_ipw_msm_error",
      call = call
    )
  }

  if (is_linear_regression(outcome_mod)) {
    family <- "gaussian"
    out_link <- "identity"
  } else {
    family <- "binomial"
    out_link <- outcome_mod$family$link
  }

  list(
    exposure_type = "continuous",
    estimand = estimand,
    n = length(exposure),
    exposure = exposure,
    ps = list(
      X = ps_X,
      link = NULL,
      coefs = stats::coef(ps_mod),
      k = NULL
    ),
    stab = list(
      stabilized = is_stabilized(wts),
      score = stabilization_score(wts)
    ),
    outcome = list(
      X = out_X,
      y = ipw_outcome_numeric(outcome),
      family = family,
      link = out_link,
      coefs = stats::coef(outcome_mod),
      X_counterfactual = NULL,
      weights = as.double(wts),
      exposure_col = exposure_cols
    ),
    contrasts = NULL,
    focal_level = NULL,
    reference_level = NULL
  )
}

# ---- weight-consistency preflight -------------------------------------------

# Reject a propensity score model that separates. A model whose covariates
# predict the exposure without error has no finite maximum likelihood estimate,
# and `glm` returns whatever its convergence criterion stopped at, which can be
# a linear predictor in the tens of thousands. Rebuilding the propensity scores
# from those coefficients gives probabilities of exactly zero or one, and a unit
# whose own probability is zero has no finite weight.
#
# Two things this guard has to be checked against, neither visible from the code.
#
# It is not redundant with the weight comparison that follows. `glm` never hands
# back an exact zero or one, because `binomial()$linkinv` clamps its argument, so
# the weights a user builds from a separated fit are finite and the comparison
# sees a mismatch rather than a singularity. Under `ate` that produced an error
# about the weights and the focal level, neither of which was the cause; under
# the other estimands the tilt cancels the singularity, the comparison passed,
# and the solver failed later with an error from deli. This runs first so the
# diagnosis names the model.
#
# The rebuild deliberately keeps using `plogis` rather than the fitted family's
# clamped inverse link. The clamp would flatten the psi derivatives past its
# threshold and quietly corrupt the sandwich variance, so the rebuild stays
# mathematically honest and a fit it cannot represent is rejected here instead.
#
# A continuous exposure has no saturating inverse link and needs no such check.
check_ipw_ps_separation <- function(n_saturated, call = rlang::caller_env()) {
  if (n_saturated > 0) {
    abort(
      c(
        "{.arg wt_mod} must not separate the exposure.",
        x = "Rebuilding the propensity scores gives a probability of exactly \\
        0 or 1 for {n_saturated} observation{?s}, whose weight{?s} {?is/are} \\
        then undefined.",
        i = "This is separation: some covariate pattern predicts the exposure \\
        without error, so the fit has no finite maximum likelihood estimate.",
        i = "Check overlap in {.arg wt_mod} rather than the weights. Dropping \\
        or combining the covariate that separates, or penalizing the fit, \\
        gives a model with finite coefficients."
      ),
      error_class = "propensity_ipw_separation_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Recompute the weights implied by the spec at its seeded init, mirroring the
# weight computation each psi builder performs at theta = init.
ipw_weights_at_init <- function(spec, layout, call = rlang::caller_env()) {
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
      check_ipw_ps_separation(sum(e == 0 | e == 1), call = call)
      stab_prob <- if (length(th_stab)) th_stab[[1]] else NULL
      weight_fn(e, spec$exposure, list(stab_prob = stab_prob, score = score))
    },
    categorical = {
      k <- spec$ps$k
      ps_mat <- ipw_categorical_ps(spec$ps$X, th_ps, k)
      # Only the score at the level each unit was actually assigned divides the
      # weight. A softmax row can underflow to an exact zero in the columns for
      # the levels a unit was not assigned, and under ordinary separation that is
      # where the zeros land, leaving the denominator positive and the analysis
      # sound. Firing on those would reject working fits.
      check_ipw_ps_separation(
        sum(rowSums(spec$exposure * ps_mat) == 0),
        call = call
      )
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

  ipw_compare_weights(
    ipw_weights_at_init(spec, layout, call = call),
    observed_wts,
    spec$exposure_type,
    spec$estimand,
    call = call
  )
}

# Compare weights recomputed from the propensity score model against the ones
# the outcome model was actually fit with, and reject a disagreement. Shared by
# the M-estimation preflight above and the linearization path, which recomputes
# its weights from predicted propensity scores rather than from a spec but has
# to reach the same verdict and say the same thing when it fails.
ipw_compare_weights <- function(
  recomputed,
  observed_wts,
  exposure_type,
  estimand,
  call = rlang::caller_env()
) {
  recomputed <- as.double(recomputed)
  observed <- as.double(observed_wts)

  consistent <- length(recomputed) == length(observed) &&
    isTRUE(all.equal(recomputed, observed, tolerance = 1e-6))

  if (!consistent) {
    msg <- c(
      "The weights used to fit {.arg outcome_mod} are not consistent with \\
      the propensity score model and estimand.",
      i = "The {.val {estimand}} weights recomputed from {.arg wt_mod} \\
      differ from the weights supplied to {.arg outcome_mod} (compared at \\
      relative tolerance 1e-6).",
      i = "Refit {.arg outcome_mod} with weights from this propensity score \\
      model and estimand."
    )
    # The focal level is a common cause of a mismatch, but for different reasons
    # per exposure type, so the hint is only offered where it applies: a binary
    # exposure has a fixed focal convention the weights can disagree with, a
    # categorical exposure takes its focal level from the weights or the
    # argument, and a continuous exposure has no focal level at all.
    focal_hint <- switch(
      exposure_type,
      binary = "A non-default {.arg .focal_level} or {.arg .reference_level} \\
      in the weights is one cause: {.fun ipw} treats the second sorted level \\
      of a binary exposure as focal.",
      categorical = "Weights built with a different {.arg .focal_level} than \\
      the one {.fun ipw} resolved are one cause.",
      NULL
    )
    if (!is.null(focal_hint)) {
      msg <- c(msg, i = focal_hint)
    }
    abort(
      msg,
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
# Effect label for a continuous exposure, taken from the outcome-model link:
# the reported effect is the single marginal structural model coefficient.
ipw_continuous_effect_label <- function(link) {
  switch(
    link,
    identity = "slope",
    logit = "log(or)",
    log = "log(rr)",
    abort(
      "Unsupported outcome link {.val {link}} for a continuous exposure.",
      error_class = "propensity_ipw_link_error"
    )
  )
}

ipw_mestimation_estimates <- function(spec, fit, layout, conf_level) {
  co <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))

  # A continuous exposure reports the marginal structural model exposure
  # coefficient, addressed by its outcome-block theta position; binary and
  # categorical exposures report the contrast rows. Positions are used, not
  # names, since theta names are not unique across blocks.
  if (identical(spec$exposure_type, "continuous")) {
    idx <- layout$idx$out[spec$outcome$exposure_col]
    effect <- ipw_continuous_effect_label(spec$outcome$link)
  } else {
    idx <- layout$idx$contrast
    effect <- rep(spec$contrasts, times = ipw_n_comparisons(spec))
  }

  estimate <- unname(co[idx])
  std.err <- unname(se[idx])
  z <- estimate / std.err

  z_val <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci.lower <- estimate - z_val * std.err
  ci.upper <- estimate + z_val * std.err
  p.value <- 2 * (1 - stats::pnorm(abs(z)))

  out <- data.frame(
    effect = effect,
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
