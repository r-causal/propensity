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

# Require the counterfactual designs to be as wide as the design the outcome
# model was fit to. The marginal mean at each level is `X %*% beta`, a positional
# multiply that consults no names, so a design of the wrong width either fails to
# multiply at all or pairs each column with a coefficient that belongs to
# another. The width only ever moves when the designs are rebuilt from `.data`,
# because the model frame reproduces the fit exactly, which is the mirror of the
# check the propensity design gets in check_ipw_ps_design_width().
#
# A column type that disagrees with the fit is rejected before this, so what
# reaches here is a categorical column whose observed levels differ from the ones
# the fit recorded: a level that `.data` collapses away, or one it declares that
# the fit never saw. Only the count is checked, for the same reason it is there:
# a covariate recoded to the same width carries the same numbers.
check_ipw_outcome_design_width <- function(
  designs,
  outcome_mod,
  call = rlang::caller_env()
) {
  n_fitted <- length(stats::coef(outcome_mod))
  widths <- vapply(designs, ncol, integer(1))
  bad <- which(widths != n_fitted)

  if (length(bad) == 0) {
    return(invisible(TRUE))
  }

  n_rebuilt <- widths[[bad[[1]]]]

  abort(
    c(
      "{.arg .data} must rebuild the design {.arg outcome_mod} was fit to.",
      x = "The counterfactual design rebuilt from {.arg .data} has \\
      {n_rebuilt} column{?s}.",
      x = "{.arg outcome_mod} was fit with {n_fitted} coefficient{?s}.",
      i = "A column whose levels differ from the fitting data is the usual \\
      cause: a factor contributes one design column per non-reference level."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# Require the values the binary counterfactual rebuild writes to carry the coding
# the outcome model holds the exposure under. The rebuild assigns one value to
# the whole exposure column, so a factor or logical column hands the rebuild a
# value that still knows every level while a character column hands it one
# string, which `model.frame()` levels on its own and reduces to a single level.
# That died inside `model.matrix()` as "contrasts can be applied only to factors
# with 2 or more levels".
#
# The check is on the value written rather than on the column supplied, which is
# why the type check cannot stand in for it: a character column and a factor
# column are the same design coding and differ only in what survives the
# assignment. It cannot arise without `.data`, where the exposure comes out of the
# propensity model's own frame.
#
# The categorical path has no counterpart because it resolves the exposure to a
# factor of the fitted levels before it writes anything.
check_ipw_binary_exposure_coding <- function(
  outcome_mod,
  exposure_name,
  exposure_values,
  call = rlang::caller_env()
) {
  data_classes <- attr(stats::terms(outcome_mod), "dataClasses")
  fit_class <- data_classes[exposure_name]
  categorical <- !is.na(fit_class) &&
    identical(ipw_design_coding(fit_class[[1]]), "categorical")

  if (
    !categorical || is.factor(exposure_values) || is.logical(exposure_values)
  ) {
    return(invisible(TRUE))
  }

  fit_levels <- outcome_mod$xlevels[[exposure_name]]
  supplied_article <- ipw_class_articles[[stats::.MFclass(exposure_values)]]

  abort(
    c(
      "{.arg .data} must supply {.val {exposure_name}} as a column that \\
      carries its own levels.",
      x = "{.arg .data} has {.val {exposure_name}} as {supplied_article}.",
      x = "{.fun ipw} rebuilds the outcome design with {.val {exposure_name}} \\
      set to one value at a time, and that leaves {supplied_article} holding \\
      the single level it was set to rather than the levels \\
      {.val {fit_levels}} {.arg outcome_mod} was fit on.",
      i = "Supply {.val {exposure_name}} as a factor with those levels."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
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
  check_ipw_outcome_response(outcome_mod, call = call)
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

  check_ipw_outcome_exposure(outcome_mod, exposure_name, call = call)

  extracted <- ipw_extract_ps_design(
    ps_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    check_exposure_levels = TRUE,
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
      error_class = "propensity_length_error",
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
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  z <- ipw_recode_binary_exposure(exposure)

  check_ipw_binary_exposure_coding(
    outcome_mod,
    exposure_name,
    exposure_values,
    call = call
  )

  # The rebuilt design multiplies the outcome model's coefficients positionally,
  # so it has to reproduce the coding those coefficients were fit under.
  # `contrasts.arg` carries it: setting the exposure column drops any contrasts
  # attribute it had, and a factor-free fit records NULL, which is the default.
  out_terms <- stats::delete.response(stats::terms(outcome_mod))
  counterfactual_mm <- function(value) {
    d <- mm_data
    d[[exposure_name]] <- value
    ipw_counterfactual_design(
      outcome_mod,
      out_terms,
      d,
      rebuilt = !is.null(.data),
      call = call
    )
  }
  x1 <- counterfactual_mm(exposure_values[[2]])
  x0 <- counterfactual_mm(exposure_values[[1]])

  check_ipw_outcome_design_width(list(x0, x1), outcome_mod, call = call)
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
#
# `counterfactual` says whether the caller goes on to rebuild the outcome design
# at each exposure value. The continuous path does not: it reports a coefficient
# from the fitted marginal structural model, so it writes no counterfactual
# column and the rebuild check would probe every distinct dose for nothing.
#
# `check_exposure_levels` says whether the exposure column supplied in `.data`
# has to carry the level order the propensity model was fit under. The binary
# path recodes on the order the column declares and so requires it; the
# categorical path resolves the column against `ps_mod$lev` instead and passes
# `FALSE`, as does the continuous path, whose exposure has no levels.
ipw_extract_ps_design <- function(
  ps_mod,
  outcome_mod,
  .data,
  exposure_name,
  counterfactual = TRUE,
  check_exposure_levels = FALSE,
  call = rlang::caller_env()
) {
  # First, and independent of `.data`: the propensity model has to have a
  # coefficient for every column of the design that follows. Both routes below
  # produce the full design, and everything downstream multiplies the two
  # positionally.
  check_ipw_model_rank(stats::coef(ps_mod), "wt_mod", call = call)

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

    # The outcome model's frame gets the same treatment. Every public `ipw()`
    # method reads the weights out of that frame before any design work, so a
    # frame-gone outcome model is rejected there first and never arrives here;
    # this guard is for callers that reach the spec constructors directly.
    outcome_extract <- tryCatch(
      list(
        outcome = fmla_extract_left_vctr(outcome_mod),
        mm_data = model.frame(outcome_mod)
      ),
      error = function(e) e
    )
    if (inherits(outcome_extract, "error")) {
      abort_outcome_frame_gone(conditionMessage(outcome_extract), call = call)
    }
    outcome <- outcome_extract$outcome
    mm_data <- outcome_extract$mm_data
    if (counterfactual) {
      check_ipw_exposure_rebuild(
        outcome_mod,
        exposure_name,
        exposure,
        mm_data,
        rebuilt = FALSE,
        call = call
      )
      # Runs second so a term the fit recorded levels for keeps the more
      # specific report above, which names the levels it was fit with.
      check_ipw_exposure_call_terms(outcome_mod, exposure_name, call = call)
    }
    check_exposure(mm_data, exposure_name, call = call)
  } else {
    # Every column the rebuilds read has to be there before any of them runs: a
    # covariate missing from `.data` otherwise surfaces as a raw
    # object-not-found from `model.matrix` rather than as the missing column it
    # is.
    required <- ipw_required_columns(ps_mod, outcome_mod, exposure_name)
    assert_columns_exist(.data, required, call = call)

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
    outcome_frame <- stats::model.frame(outcome_mod)
    n_fitted <- nrow(outcome_frame)

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

    check_ipw_data_types(
      .data,
      ps_mod,
      outcome_mod,
      exposure_name,
      call = call
    )

    exposure <- .data[[exposure_name]]

    if (check_exposure_levels) {
      check_ipw_exposure_levels(exposure, ps_mod, exposure_name, call = call)
    }

    if (counterfactual) {
      check_ipw_exposure_rebuild(
        outcome_mod,
        exposure_name,
        exposure,
        .data,
        rebuilt = TRUE,
        call = call
      )
    }

    # Evaluated rather than looked up, for the same reason the assert above reads
    # the response's variables: the response may be a transformation of a column
    # rather than a column itself, and it has to be computed the way the fit
    # computed it.
    outcome <- fmla_eval_left(outcome_mod, .data)
    check_ipw_response_levels(
      outcome,
      levels(stats::model.response(outcome_frame)),
      fmla_extract_left_chr(outcome_mod),
      call = call
    )
    mm_data <- .data

    # Last of the `.data` guards, and immediately before the first rebuild:
    # every design below is built with `model.frame()`, which drops a row
    # holding a missing value in any column it reads, while the weights, the
    # exposure, and the outcome values keep every row of `.data`. Each guard
    # above diagnoses one column on its own terms and is unaffected by missing
    # values in it, so they keep their place and this one covers what is left.
    check_ipw_data_complete(.data, required, call = call)

    ps_X <- ipw_rebuild_design(
      ps_mod,
      stats::delete.response(stats::terms(ps_mod)),
      .data,
      call = call
    )
    check_ipw_ps_design_width(ps_X, ps_mod, call = call)
  }

  list(exposure = exposure, outcome = outcome, ps_X = ps_X, mm_data = mm_data)
}

# Terms per estimating equation in a fitted propensity model. A glm or lm has one
# equation and a coefficient vector as long as its design. nnet::multinom has one
# per non-reference level and returns a level-by-term matrix, except at two
# levels, where it returns that single row as a bare vector.
ipw_ps_terms_per_equation <- function(ps_mod) {
  coefs <- stats::coef(ps_mod)

  if (is.matrix(coefs)) ncol(coefs) else length(coefs)
}

# Variable names a model's right-hand side reads, so a rebuild from `.data` can
# be told which columns it needs before it goes looking for them.
ipw_model_covariates <- function(mod) {
  all.vars(stats::delete.response(stats::terms(mod)))
}

# The columns of `.data` the rebuilds read. Both models' covariates are needed,
# not just the exposure and the outcome: the designs are rebuilt from `.data`,
# and a covariate missing from it otherwise surfaces as a raw object-not-found
# from `model.matrix` rather than as the missing column it is.
#
# The outcome model's response contributes the variables it reads rather than
# the deparsed expression: a transformed response such as `log(y)` is computed
# from a column named `y`, and asking for one named `log` would report a column
# that could not exist under any correct `.data`.
ipw_required_columns <- function(ps_mod, outcome_mod, exposure_name) {
  unique(c(
    exposure_name,
    fmla_extract_left_vars(outcome_mod),
    ipw_model_covariates(ps_mod),
    ipw_model_covariates(outcome_mod)
  ))
}

# Reject a `.data` that is missing values in a column the rebuilds read. Every
# design is rebuilt with `model.frame()`, which drops a row holding a missing
# value in any column it reads, while the weights come from the outcome model
# fit and the exposure and outcome values are read from `.data` itself, each
# keeping every row. The design and the vectors it multiplies then differ in
# length and are recycled against each other: base R warned twice that one
# object was not a multiple of the other, and the mismatch surfaced as weights
# that failed their consistency check, which reports how the weights were built
# when what was wrong was the frame that was passed.
#
# Only the columns the rebuilds read are swept. A column `.data` happens to
# carry and neither model consults is never in a model frame and cannot drop a
# row.
#
# A model fit on data with missing values dropped those rows itself, so its
# frame is shorter than the frame it was fit from and the row-count check above
# reports that first. This guard is for the `.data` that has the right number of
# rows and holes in it.
check_ipw_data_complete <- function(
  .data,
  columns,
  call = rlang::caller_env()
) {
  present <- intersect(columns, names(.data))
  missing_by_column <- vapply(
    present,
    function(nm) sum(is.na(.data[[nm]])),
    integer(1)
  )
  incomplete <- missing_by_column[missing_by_column > 0]

  if (length(incomplete) == 0) {
    return(invisible(TRUE))
  }

  cols <- names(incomplete)
  n_rows <- sum(Reduce(`|`, lapply(.data[cols], is.na)))

  abort(
    c(
      "{.arg .data} must have no missing values in the columns the models \\
      read.",
      x = "{.val {cols}} {?has/have} missing values in {.arg .data}.",
      x = "{n_rows} row{?s} of {.arg .data} {?is/are} incomplete, and the \\
      rebuilds drop {?it/them}.",
      i = "Every design {.fun ipw} rebuilds from {.arg .data} drops a row that \\
      is missing a value, while the weights, the exposure, and the outcome \\
      values keep every row, so the two are then recycled against each other.",
      i = "Supply the data the models were fit to, or drop the incomplete rows \\
      and refit both models on what is left."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# Reject a model whose design has columns it could not separate. `lm` and `glm`
# pivot such a column out of the fit and record `NA` for its coefficient, while
# `model.matrix()` keeps returning the column, so the coefficient vector and the
# design disagree in a way only the `NA` shows.
#
# The `NA` is what reaches `ipw()`, and it is the whole signal: nothing here
# measures the design's rank or its conditioning, because a fit that kept a
# coefficient for every column carries no `NA` however badly conditioned it is,
# and one that dropped a column carries one however the dependence arose. An
# exact duplicate and a copy perturbed below the pivot tolerance are both
# dropped, so both are caught.
#
# Untreated, the `NA` propagated differently per model and per path, and none of
# the reports named a model or a column. From the propensity model it reached
# the rebuilt propensity scores, and from there the count of saturated scores
# and an `if`, which stopped with "missing value where TRUE/FALSE needed"; on
# the linearization path it reached a solve, which reported an exactly singular
# system from LAPACK; for a continuous exposure it reached the recomputed
# weights, which disagreed with the supplied ones and were reported as a
# disagreement about the estimand. From the outcome model it reached the seeded
# starting values, and the M-estimator reported that `stacked_equations`
# returned non-finite values at `init`.
#
# `nnet::multinom` optimizes rather than pivots, so it returns finite
# coefficients for a dependent column and is passed here. Its design and its
# coefficient matrix still agree, which is what this guard is about.
#
# The coefficients are taken rather than the model, because the outcome model's
# are already read into the spec by the time this runs on that side.
check_ipw_model_rank <- function(coefs, arg, call = rlang::caller_env()) {
  if (is.matrix(coefs)) {
    dropped <- colnames(coefs)[apply(!is.finite(coefs), 2, any)]
  } else {
    dropped <- names(coefs)[!is.finite(coefs)]
  }

  if (length(dropped) == 0) {
    return(invisible(TRUE))
  }

  consequence <- if (identical(arg, "wt_mod")) {
    "{.fun ipw} rebuilds the propensity scores by multiplying the fitted \\
    coefficients against that design, so a column with no coefficient leaves \\
    every score undefined."
  } else {
    "{.fun ipw} estimates the marginal means by multiplying the fitted \\
    coefficients against that design, so a column with no coefficient leaves \\
    every predicted outcome undefined."
  }

  abort(
    c(
      "{.arg {arg}} must have a coefficient for every column of its design.",
      x = "{.arg {arg}} has no fitted coefficient for {.val {dropped}}.",
      i = "A model reports that for a column its design cannot separate from \\
      the others: the column is a linear combination of them, exactly or to \\
      within the tolerance the fit pivots at, so the fit has no unique \\
      solution for it and drops it.",
      i = consequence,
      i = "Refit {.arg {arg}} without the redundant column{?s}, or combine \\
      {?it/them} with the column{?s} {?it duplicates/they duplicate}."
    ),
    error_class = "propensity_ipw_rank_error",
    call = call
  )
}

# Require the design rebuilt from `.data` to be as wide as the one the propensity
# model was fit to. Everything downstream multiplies the two together
# positionally: the propensity scores, the stacked score block, and, for a
# categorical exposure, the coefficient names built by pairing each level with
# each design column. A width that disagrees produced a raw error about a names
# attribute or a non-conformable multiply, neither of which mentions `.data`.
#
# A column whose type differs between the fitting data and `.data` is rejected
# before this by check_ipw_data_types(), which names the column and both types
# and is the more specific diagnosis. What is left for a count to catch is a term
# recorded under a call rather than a variable, which the type check has no
# column to compare and which can still rebuild to a different width.
check_ipw_ps_design_width <- function(
  ps_X,
  ps_mod,
  call = rlang::caller_env()
) {
  n_fitted <- ipw_ps_terms_per_equation(ps_mod)

  if (!identical(ncol(ps_X), n_fitted)) {
    abort(
      c(
        "{.arg .data} must rebuild the design {.arg wt_mod} was fit to.",
        x = "The design rebuilt from {.arg .data} has {ncol(ps_X)} column{?s}.",
        x = "{.arg wt_mod} was fit with {n_fitted} coefficient{?s} per \\
        equation.",
        i = "A column whose type differs from the fitting data is the usual \\
        cause: a factor expands to one column per non-reference level where a \\
        numeric takes one."
      ),
      error_class = "propensity_ipw_data_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject a `.data` column whose type contradicts the type either model was fit
# on. A fit records each variable's class, and every design rebuilt from `.data`
# is rebuilt under the coding that class implies: a numeric column takes one
# design column, while a factor, ordered factor, character, or logical column
# expands into one column per non-reference level and carries the fit's contrast
# coding. Supplying the other kind therefore either fails to multiply the
# coefficients at all or, when the two happen to be the same width, multiplies
# them by different numbers with nothing signaled. Both directions died raw: a
# factor-fit column supplied as numeric as "contrasts apply only to factors", a
# numeric-fit column supplied as a factor as a non-conformable multiply, and a
# numeric-fit column supplied as a two-level factor or as a logical not at all.
#
# Both models are swept, and the outcome model's response with them. This runs
# before the propensity design is rebuilt and before the outcome model's
# counterfactual designs are built on every path, so an outcome-only covariate is
# caught on its way to the same failure.
#
# A name that is not a column of `.data`, `factor(x)` say, is skipped: it is a
# term recorded under a call rather than a variable this can compare, and the
# rebuild reaches it with the column it reads untouched.
#
# The exposure is swept too. Its fitted type comes from the propensity model,
# whose response it is, and the counterfactual rebuild overwrites the supplied
# column, so a type that disagrees is a mistake there for the same reason it is
# anywhere else. What the written value has to carry beyond its type is the
# binary path's own check.
#
# Only the first offending column is reported: the fix is a recoding of that
# column, after which the guard runs again over the rest.
check_ipw_data_types <- function(
  .data,
  ps_mod,
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  fitted <- ipw_fitted_classes(ps_mod, outcome_mod)
  fitted <- fitted[names(fitted) %in% names(.data)]
  response <- fmla_extract_left_chr(outcome_mod)

  for (column in names(fitted)) {
    fit_class <- fitted[[column]]
    supplied <- stats::.MFclass(.data[[column]])

    mismatch <- if (identical(column, response)) {
      # The response is not a design column. It is converted by
      # ipw_outcome_numeric(), which indicator-codes a factor against its first
      # level and uses the values themselves otherwise, so what has to agree is
      # which of those two conversions applies.
      !identical(
        fit_class %in% c("factor", "ordered"),
        supplied %in% c("factor", "ordered")
      ) ||
        !supplied %in% c("factor", "ordered", "numeric", "logical")
    } else {
      !ipw_types_interchangeable(fit_class, supplied)
    }

    if (mismatch) {
      abort_ipw_type_mismatch(
        column = column,
        fit_class = fit_class,
        supplied = .data[[column]],
        supplied_class = supplied,
        ps_mod = ps_mod,
        outcome_mod = outcome_mod,
        response = identical(column, response),
        call = call
      )
    }
  }

  invisible(TRUE)
}

# The class each model recorded for every variable it reads, keyed by variable
# name. The propensity model is taken first so the exposure, which is its
# response, is described by the fit that owns it. Classes a design has no coding
# for are dropped rather than compared.
ipw_fitted_classes <- function(ps_mod, outcome_mod) {
  classes <- c(
    attr(stats::terms(ps_mod), "dataClasses"),
    attr(stats::terms(outcome_mod), "dataClasses")
  )
  classes <- classes[!duplicated(names(classes))]

  classes[classes %in% names(ipw_class_nouns)]
}

# Whether a `.data` column of class `supplied` rebuilds the design the fit
# recorded under `fit_class`.
#
# A fit that recorded levels for a column rebuilds it by re-leveling the values
# it is handed against the levels it recorded, so what a supplied column has to
# carry is values that name those levels. A factor, an ordered factor, and a
# character vector all do, and are interchangeable.
#
# A logical column does not, and it is rejected however its values line up. It
# has no levels to re-level, so `model.matrix()` codes it on its own terms,
# `FALSE` first, and the design lands on whichever mapping that gives. Against a
# two-level fitted factor the widths agree and nothing further notices: an
# aligned logical reproduces the fitted design by coincidence, while one whose
# `TRUE` marks the fit's first level silently swaps the two levels'
# coefficients and moves the estimates. The rejection is keyed on the type
# rather than on the values, because a logical column carries nothing that says
# which level its `TRUE` was meant to name.
#
# The same rejection runs the other way. A fit that recorded a logical column
# recorded no levels either, so it took `FALSE` as its reference and the design
# rebuilt from `.data` has to arrive on that same coding. A factor, an ordered
# factor, and a character vector all bring their own reference: the widths
# agree, and the design lands on whichever level sorts or declares first, which
# is the fitted coding only by coincidence. A factor declaring `TRUE` before
# `FALSE` reproduces the fitted design reversed and moves the estimates with
# nothing signaled, and a character column matches only when its values happen
# to sort the way `FALSE` and `TRUE` do.
#
# Everything else is compared by the coding its class implies, which is what
# decides how many design columns it takes.
ipw_types_interchangeable <- function(fit_class, supplied) {
  leveled <- c("factor", "ordered", "character")

  if (fit_class %in% leveled) {
    return(supplied %in% leveled)
  }

  if (identical(fit_class, "logical")) {
    return(identical(supplied, "logical"))
  }

  identical(ipw_design_coding(fit_class), ipw_design_coding(supplied))
}

# The design coding a class implies: one column, or one column per non-reference
# level.
ipw_design_coding <- function(class) {
  if (class %in% c("factor", "ordered", "character", "logical")) {
    "categorical"
  } else if (identical(class, "numeric")) {
    "numeric"
  } else {
    class
  }
}

ipw_class_nouns <- c(
  numeric = "numeric column",
  factor = "factor",
  ordered = "ordered factor",
  character = "character vector",
  logical = "logical vector"
)

ipw_class_articles <- c(
  numeric = "a numeric vector",
  factor = "a factor",
  ordered = "an ordered factor",
  character = "a character vector",
  logical = "a logical vector"
)

abort_ipw_type_mismatch <- function(
  column,
  fit_class,
  supplied,
  supplied_class,
  ps_mod,
  outcome_mod,
  response,
  call = rlang::caller_env()
) {
  fit_noun <- ipw_class_nouns[[fit_class]]
  fit_article <- ipw_class_articles[[fit_class]]
  supplied_article <- if (supplied_class %in% names(ipw_class_articles)) {
    ipw_class_articles[[supplied_class]]
  } else {
    paste("a", class(supplied)[[1]], "vector")
  }
  supplied_noun <- if (supplied_class %in% names(ipw_class_nouns)) {
    ipw_class_nouns[[supplied_class]]
  } else {
    paste(class(supplied)[[1]], "vector")
  }

  if (response) {
    abort(
      c(
        "{.arg .data} must supply the outcome {.val {column}} as the \\
        {fit_noun} {.arg outcome_mod} was fit on.",
        x = "{.arg .data} has {.val {column}} as {supplied_article}.",
        x = "{.fun ipw} reads the outcome values from {.arg .data}: a factor \\
        response becomes an indicator for its non-first levels and any other \\
        response is used as its own values, so the two are not \\
        interchangeable.",
        i = "Supply {.val {column}} as that {fit_noun}."
      ),
      error_class = "propensity_ipw_data_error",
      call = call
    )
  }

  fit_args <- c(
    if (column %in% names(attr(stats::terms(ps_mod), "dataClasses"))) "wt_mod",
    if (column %in% names(attr(stats::terms(outcome_mod), "dataClasses"))) {
      "outcome_mod"
    }
  )
  fit_levels <- c(ps_mod$xlevels, outcome_mod$xlevels)[[column]]

  recorded <- if (is.null(fit_levels)) {
    "{.arg {fit_args}} recorded {.val {column}} as {fit_article}, and the \\
    designs rebuilt from {.arg .data} use that coding."
  } else {
    "{.arg {fit_args}} recorded {.val {column}} as {fit_article} with the \\
    levels {.val {fit_levels}}, and the designs rebuilt from {.arg .data} use \\
    that coding."
  }

  abort(
    c(
      "{.arg .data} must supply {.val {column}} as the {fit_noun} the models \\
      were fit with.",
      x = "{.arg .data} has {.val {column}} as {supplied_article}.",
      x = recorded,
      i = "Supply {.val {column}} as that {fit_noun}, or refit the models on \\
      the {supplied_noun}."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
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

# Rebuild a fitted model's design from a supplied data frame, under the levels
# the fit recorded. Left to itself, `model.matrix()` takes each categorical
# column's levels from the values in front of it: a character column is ordered
# alphabetically and a factor is taken at whatever order it declares. Either one
# pairs each level with another level's coefficient, which the multiply cannot
# see, so the fit's own `xlevels` decide instead. Its contrast coding is passed
# separately, because re-leveling a column drops the attribute that carries it,
# and the attribute on the supplied column is dropped first so the two cannot
# disagree. NULL, for a factor-free fit, is `contrasts.arg`'s default.
#
# The frame is built first and re-leveled afterwards rather than handing `xlev`
# to `model.frame()`, which does the same work but reports two situations in its
# own vocabulary: a value the fit never saw stops it with a message naming
# neither `.data` nor `ipw()`, and a column that is neither a factor nor a
# character vector makes it warn that a variable is not a factor, when what is
# wrong there is the column's type. A column the fit recorded levels for reaches
# this as a factor, an ordered factor, or a character vector, because
# check_ipw_data_types() rejects every other type for such a column and names
# both types when it does.
ipw_rebuild_design <- function(
  mod,
  mod_terms,
  data,
  call = rlang::caller_env()
) {
  frame <- stats::model.frame(
    mod_terms,
    data = drop_contrasts_attrs(data, names(mod$contrasts))
  )
  frame <- ipw_relevel_frame(frame, mod$xlevels, call = call)

  model.matrix(mod_terms, data = frame, contrasts.arg = mod$contrasts)
}

# Re-level the categorical columns of a rebuilt model frame to the levels the
# fit recorded, rejecting a value the fit never saw. A column the fit recorded
# levels for arrives as a factor, an ordered factor, or a character vector;
# anything else is left alone, because its type is the disagreement rather than
# its levels.
#
# A level the supplied column declares but no observation carries is dropped
# here, which is what the fits themselves do, so a column that declares one
# rebuilds the design the model was fit to. A level an observation does carry
# has no coefficient to multiply and is rejected.
#
# The frame is keyed by the fit's own variable names, which are the deparsed
# model variables, so a term recorded under a call rather than a column, such as
# `factor(x)`, is re-leveled here too: the call has been evaluated by the time
# the frame exists.
ipw_relevel_frame <- function(frame, xlev, call = rlang::caller_env()) {
  for (nm in intersect(names(xlev), names(frame))) {
    column <- frame[[nm]]

    if (!is.factor(column) && !is.character(column)) {
      next
    }

    levs <- xlev[[nm]]
    values <- as.character(column)
    unseen <- unique(values[!is.na(values) & !values %in% levs])

    if (length(unseen) > 0) {
      abort_ipw_new_levels(nm, unseen, levs, call = call)
    }

    frame[[nm]] <- factor(values, levels = levs, ordered = is.ordered(column))
  }

  frame
}

abort_ipw_new_levels <- function(
  term,
  unseen,
  levs,
  call = rlang::caller_env()
) {
  abort(
    c(
      "{.arg .data} must hold only the levels the models were fit on.",
      x = "{.code {term}} takes the value{?s} {.val {unseen}} in {.arg .data}.",
      x = "The models were fit with {.code {term}} at {.val {levs}}, and the \\
      designs rebuilt from {.arg .data} carry one column per non-reference \\
      level of that set, so a value outside it has no coefficient to multiply.",
      i = "Supply the data the models were fit to, or refit them on data that \\
      holds every level you want represented."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# Build one counterfactual outcome design. On the `.data` route the design is
# rebuilt from the supplied frame, so every term is re-evaluated at the exposure
# value written into it and re-leveled against the fit. Without `.data` the
# design comes from the outcome model's own model frame, which carries its terms
# attribute: `model.matrix()` then reads the columns already in it and evaluates
# nothing, so a term recorded under a call keeps its fitted values whatever is
# written beside them. check_ipw_exposure_rebuild() rejects such a term on that
# route for exactly that reason.
ipw_counterfactual_design <- function(
  outcome_mod,
  mod_terms,
  data,
  rebuilt,
  call = rlang::caller_env()
) {
  if (rebuilt) {
    return(ipw_rebuild_design(outcome_mod, mod_terms, data, call = call))
  }

  model.matrix(mod_terms, data = data, contrasts.arg = outcome_mod$contrasts)
}

# Require the response column supplied in `.data` to induce the coding the
# outcome model was fit under. ipw_outcome_numeric() indicator-codes a factor
# response against its first level, following `glm`, so the level the fit treats
# as the failure has to come first here as well. A response re-leveled the other
# way is the same values under the opposite coding: on the M-estimation path
# every estimate comes back with its sign reversed, and on the linearization
# path the point estimates stay, because they are the outcome model's own
# predictions, while the standard errors move with the response values that feed
# the influence functions.
#
# Only the first level decides the conversion, so only the first level is
# compared. A column that declares further levels the fit never saw, or that
# orders the rest differently, codes every observation the way the fit did.
#
# Silent for a fit whose response carries no levels, which covers a numeric or
# logical response and a response written as a transformation, and silent for a
# supplied response that is not a factor, whose type contradicts the fit and is
# rejected before this by check_ipw_data_types().
check_ipw_response_levels <- function(
  outcome,
  fit_levels,
  response_name,
  call = rlang::caller_env()
) {
  if (is.null(fit_levels) || !is.factor(outcome)) {
    return(invisible(TRUE))
  }

  supplied <- levels(outcome)

  # A factor that declares no levels holds nothing but missing values, so it has
  # no first level to compare and codes no observation at all. Left to the
  # comparison below it stopped as a subscript out of bounds.
  if (length(supplied) == 0) {
    abort(
      c(
        "{.arg .data} must supply the outcome {.val {response_name}} on the \\
        levels {.arg outcome_mod} was fit with.",
        x = "{.arg .data} has {.val {response_name}} as a factor that declares \\
        no levels, so every observation's outcome is missing.",
        i = "{.arg outcome_mod} was fit with {.val {fit_levels}}.",
        i = "Supply {.val {response_name}} as a factor on those levels, or \\
        supply the data the models were fit to."
      ),
      error_class = "propensity_ipw_data_error",
      call = call
    )
  }

  if (identical(supplied[[1L]], fit_levels[[1L]])) {
    return(invisible(TRUE))
  }

  fit_first <- fit_levels[[1L]]
  supplied_first <- supplied[[1L]]

  abort(
    c(
      "{.arg .data} must supply the outcome {.val {response_name}} with the \\
      level {.arg outcome_mod} was fit to treat as the failure first.",
      x = "{.arg .data} declares {.val {supplied_first}} first; \\
      {.arg outcome_mod} was fit with {.val {fit_first}} first.",
      x = "{.fun ipw} codes a factor outcome as an indicator for its non-first \\
      levels, following {.fun glm}, so the two orders describe opposite \\
      outcomes and every estimate reverses.",
      i = "Re-level {.val {response_name}} so {.val {fit_first}} comes first, \\
      or supply the data the models were fit to."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# Require the exposure column supplied in `.data` to carry the level order the
# propensity score model was fit under. The binary path takes the second sorted
# value of the exposure as the exposed group, which for a factor is the second
# level it declares, and the counterfactual designs, the recoded exposure, and
# the influence functions all read that coding. A column re-leveled the other
# way therefore contrasts the levels in the opposite direction. It did error,
# at the weight-consistency preflight, because weights recomputed under the
# reversed coding no longer match the ones that fit the outcome model, but that
# message is about how the weights were built rather than about `.data`.
#
# Levels the fit never saw are tolerated wherever they sit, since no observation
# carries them and `sort(unique())` cannot reach them. What has to hold is that
# the fitted levels keep their fitted order.
#
# The exposure is the propensity model's response, so `xlevels` never records
# it and the fitted frame is the only place its levels are kept; a model that
# cannot rebuild that frame leaves nothing to compare and is passed. A supplied
# column that is not a factor carries no order of its own: a numeric one is
# rejected by check_ipw_data_types() and a character one by
# check_ipw_binary_exposure_coding().
#
# The categorical path has no counterpart. It resolves the supplied column
# against `ps_mod$lev` before anything reads it, so the order it declares says
# nothing there.
check_ipw_exposure_levels <- function(
  exposure,
  ps_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  fit_levels <- ipw_fitted_response_levels(ps_mod)

  if (is.null(fit_levels) || !is.factor(exposure)) {
    return(invisible(TRUE))
  }

  supplied <- levels(exposure)

  if (identical(intersect(supplied, fit_levels), fit_levels)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg .data} must supply {.val {exposure_name}} on the levels \\
      {.arg wt_mod} was fit with, in that order.",
      x = "{.arg .data} declares {.val {supplied}}; {.arg wt_mod} was fit on \\
      {.val {fit_levels}}.",
      x = "{.fun ipw} treats the second level of a binary exposure as the \\
      exposed group, so a different order contrasts the levels the other way \\
      round.",
      i = "Re-level {.val {exposure_name}} to the order {.arg wt_mod} was fit \\
      with, or supply the data the models were fit to."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# The levels a fitted model recorded for its response, or NULL if it has none or
# if the frame that holds them cannot be rebuilt. A response is not a design
# column, so `xlevels` never carries it.
ipw_fitted_response_levels <- function(mod) {
  frame <- tryCatch(stats::model.frame(mod), error = function(e) NULL)

  if (is.null(frame)) {
    return(NULL)
  }

  levels(stats::model.response(frame))
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
  check_ipw_outcome_response(outcome_mod, call = call)
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
      error_class = "propensity_length_error",
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
    ipw_counterfactual_design(
      outcome_mod,
      out_terms,
      d,
      rebuilt = !is.null(.data),
      call = call
    )
  })
  names(x_cf) <- levs

  check_ipw_outcome_design_width(x_cf, outcome_mod, call = call)
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
  check_ipw_ps_response(ps_mod, call = call)
  check_ipw_offset(outcome_mod, call = call)
  check_ipw_outcome_response(outcome_mod, call = call)
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

  extracted <- ipw_extract_ps_design(
    ps_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    counterfactual = FALSE,
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
      error_class = "propensity_length_error",
      call = call
    )
  }

  wts <- extract_weights(outcome_mod)

  # Membership first, so that a value naming no estimand at all is reported as
  # that on this path too. The restriction below says the estimand asked for is
  # not available here, which of a value that is not an estimand anywhere read
  # as though it were one and merely unsupported for this exposure type.
  check_estimand_known(estimand, from_weights = FALSE, call = call)

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
    abort_ipw_ps_separation(
      n_saturated,
      lead = "Rebuilding the propensity scores",
      call = call
    )
  }

  invisible(TRUE)
}

# Both standard error methods refuse a saturated propensity fit at the same
# threshold, so both raise this error and only `lead` differs: the M-estimation
# path rebuilds the scores from a parameter vector, while the linearization path
# puts the fitted linear predictors back through the link's inverse. The
# diagnosis and the remedy are the same either way, and keeping the two messages
# one function apart is what keeps them the same shape. `lead` opens the sentence
# the `x` bullet completes, and the pluralization keys off `n_saturated` because
# it is the last quantity substituted before it.
abort_ipw_ps_separation <- function(
  n_saturated,
  lead,
  call = rlang::caller_env()
) {
  abort(
    c(
      "{.arg wt_mod} must not separate the exposure.",
      x = "{lead} gives a probability of exactly 0 or 1 for {n_saturated} \\
      observation{?s}, whose weight{?s} {?is/are} then undefined.",
      i = "This is usually separation: some covariate pattern predicts the \\
      exposure without error, so the fit has no finite maximum likelihood \\
      estimate. An extreme covariate pattern can saturate the scores even \\
      where the estimate is finite.",
      i = "Check overlap in {.arg wt_mod} rather than the weights. Dropping \\
      or combining the covariate that separates, or penalizing the fit, \\
      gives a model with finite coefficients."
    ),
    error_class = "propensity_ipw_separation_error",
    call = call
  )
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
      # weight. Separation drives a softmax row toward the corner at that level:
      # the assigned level's probability rounds to exactly 1 well before the
      # columns for the other levels underflow to an exact 0, which they do only
      # under stronger separation still. A denominator of 1 is harmless, and a
      # zero in a column for a level the unit was not assigned never reaches the
      # denominator at all, so firing on either would reject working fits.
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

# `layout` defaults to the layout the spec implies, so a caller holding one
# already passes it rather than paying for a second identical partition of the
# same spec.
ipw_check_weight_consistency <- function(
  spec,
  observed_wts,
  layout = ipw_theta_layout(spec, call = call),
  call = rlang::caller_env()
) {
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
    # All this comparison establishes is that two vectors disagree. Weights that
    # are exactly right reach it, so the message reports the comparison and
    # offers the causes as possibilities rather than naming the weights as the
    # fault.
    #
    # A continuous exposure has no focal level, so its causes bullet names only
    # the estimand.
    resolved_cause <- if (exposure_type == "continuous") {
      "The estimand the weights were built for may differ from the one \\
      {.fun ipw} resolved."
    } else {
      "The estimand or the focal level the weights were built for may differ \\
      from the ones {.fun ipw} resolved."
    }
    # The focal level is a common cause of a mismatch, but for different reasons
    # per exposure type, so the specific hint is only offered where it applies: a
    # binary exposure has a fixed focal convention the weights can disagree with,
    # a categorical exposure takes its focal level from the weights or the
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
    # A continuous exposure is the only one whose weights carry a spread. The
    # stacked system estimates one pooled residual variance, so weights built on
    # observation-level standard deviations are a different function of the data
    # and cannot be reproduced here at any parameter value.
    sigma_hint <- if (exposure_type == "continuous") {
      "Weights built with an observation-level {.arg .sigma}, such as \\
      {.code influence(model)$sigma}, are one cause: {.fun ipw} models the \\
      conditional density with a single pooled residual standard deviation, \\
      which is what {.fun wt_ate} uses when no {.arg .sigma} is given."
    } else {
      NULL
    }
    msg <- c(
      "The {.val {estimand}} weights recomputed from {.arg wt_mod} differ \\
      from the weights supplied to {.arg outcome_mod} (compared at relative \\
      tolerance 1e-6).",
      i = resolved_cause
    )
    if (!is.null(focal_hint)) {
      msg <- c(msg, i = focal_hint)
    }
    if (!is.null(sigma_hint)) {
      msg <- c(msg, i = sigma_hint)
    }
    msg <- c(
      msg,
      i = "Weights trimmed, truncated, or normalized after {.arg wt_mod} was \\
      fit differ from the ones rebuilt here, which come from that model \\
      alone.",
      i = "{.arg .data} values that differ from the data the models were fit \\
      to move the recomputed weights on their own and leave the supplied \\
      weights exactly right.",
      i = "Refit {.arg outcome_mod} with weights from this propensity score \\
      model and estimand if the weights are the cause."
    )
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
  # The outcome model's turn at the rank guard the propensity model gets during
  # design extraction. It runs here rather than in the spec constructors so it
  # follows the checks those constructors make on the counterfactual designs
  # themselves: a term that is constant across the exposure levels leaves the
  # designs identical *and* leaves the fit without a coefficient for it, and
  # naming the identical designs is the more specific diagnosis of the two.
  check_ipw_model_rank(spec$outcome$coefs, "outcome_mod", call = call)

  layout <- ipw_theta_layout(spec, call = call)
  psi <- build_ipw_psi(spec, layout, call = call)

  if (!is.null(spec$outcome$weights)) {
    ipw_check_weight_consistency(
      spec,
      spec$outcome$weights,
      layout = layout,
      call = call
    )
  }

  m <- deli::MEstimator(stacked_equations = psi, init = layout$init)
  solved <- ipw_solve(m)

  # Before the estimates are read, because reading them is where a fit with no
  # variance dies: the standard errors come from the variance the solve did not
  # build.
  if (solved$no_variance) {
    abort_ipw_no_variance(spec, call = call)
  }

  estimates <- ipw_mestimation_estimates(
    spec,
    solved$fit,
    layout,
    conf_level,
    call = call
  )

  if (solved$unsolved) {
    warn_ipw_unsolved(estimates, call = call)
  } else {
    warn_ipw_degenerate_se(estimates, call = call)
  }

  # The blocks the component models are reported under. The stacked system
  # estimates the propensity score model, the outcome model, and the effects
  # together, so each model's covariance here accounts for the others having been
  # estimated from the same data and is not the one its own fitting routine
  # reports.
  vc <- stats::vcov(solved$fit)

  list(
    estimates = estimates,
    fit = solved$fit,
    ps_vcov = ipw_component_vcov(vc, layout, "ps"),
    outcome_vcov = ipw_component_vcov(vc, layout, "out")
  )
}

# The component models an M-estimation result holds: the models the user
# supplied, each carrying its block of the joint sandwich, so that `vcov()` on
# either reports the covariance the stacked system implies rather than the one
# the model's own fitting routine does. Wrapping adds a class and an attribute
# and changes nothing else, so predictions, coefficients, and model frames are
# read from the wrapped model exactly as from the model itself. The models the
# user holds are untouched.
ipw_component_models <- function(wt_mod, outcome_mod, fit) {
  list(
    wt_mod = new_ipw_model(wt_mod, fit$ps_vcov),
    outcome_mod = new_ipw_model(outcome_mod, fit$outcome_vcov)
  )
}

# Run the solve, holding back the two reports deli makes about it in its own
# vocabulary. The first is that it could not pin the solution down, signaled as
# a `deli_solver_not_converged` warning whose text talks about
# `stacked_equations`, redundant parameters, and rank deficiency in the
# estimating equations, none of which the user wrote. The second is that the
# matrix of derivatives it differentiates those equations into came back holding
# values that are not finite, signaled as a `deli_bread_na` warning about a
# bread matrix that cannot be inverted. Both are muffled here and re-raised in
# `ipw()`'s terms: the first once the estimates exist, which is what it takes to
# say which of them the report applies to, and the second immediately, since a
# fit with no variance has no inference to report and the estimates are never
# built. This follows the treatment `ipw_contrast_reporter()` already gives base
# R's unclassed "NaNs produced" from the contrast transforms: the foreign
# condition is replaced rather than chained, so the user reads one message in
# one vocabulary.
#
# Each flag records only that the condition was raised. Neither carries data
# beyond its message, so nothing more can be read off them.
#
# A `deli_bread_na` and a variance that was not built are the same event: deli
# raises the warning exactly where it abandons the sandwich, so the flag stands
# in for the missing variance and nothing here reaches into the fit to ask.
ipw_solve <- function(m) {
  unsolved <- FALSE
  no_variance <- FALSE

  fit <- withCallingHandlers(
    deli::estimate(m),
    deli_solver_not_converged = function(cnd) {
      unsolved <<- TRUE
      rlang::cnd_muffle(cnd)
    },
    deli_bread_na = function(cnd) {
      no_variance <<- TRUE
      rlang::cnd_muffle(cnd)
    }
  )

  list(fit = fit, unsolved = unsolved, no_variance = no_variance)
}

# A standard error the solve did not really produce. When the estimating
# equations have no unique root, the sandwich variance along the unidentified
# direction collapses instead of growing: the zero-events fixture reports a risk
# difference of 0.66 with a standard error of 1.5e-44, a p-value of zero, and an
# interval of no width.
#
# The signature is the size of the test statistic the two make together, which
# is free of the units the outcome is measured in. A floor on the standard error
# itself is not: an outcome measured in nanomolar units gives a healthy fit an
# estimate of 6e-10 and a standard error of 1e-10, and any absolute floor near
# machine precision reports that fit, whose z is 6, as carrying no information.
#
# Measured across the exposure types and both standard error paths, healthy fits
# reach |z| = 11, including outcomes rescaled by 1e-9; the degenerate fixtures
# sit between 4.5e16 and 1.8e44. The threshold is the reciprocal of the square
# root of machine epsilon, about 6.7e7, which the healthy fits clear by nearly
# seven orders of magnitude and the degenerate ones by nearly nine.
#
# A standard error of exactly zero is degenerate whatever the estimate is, and
# is taken on its own: it is what an outcome that is constant at zero in both
# arms reports, where the ratio is 0/0 and says nothing. An estimate of exactly
# zero beside a standard error that is not gives |z| = 0 and is an honest null.
ipw_degenerate_se <- function(estimate, std.err) {
  reported <- !is.na(estimate) & !is.na(std.err)

  z <- rep(0, length(std.err))
  scaled <- reported & std.err > 0
  z[scaled] <- abs(estimate[scaled]) / std.err[scaled]

  (reported & std.err == 0) | z > 1 / sqrt(.Machine$double.eps)
}

# The rows of an estimates table whose standard errors carry that signature,
# labeled as the table labels them: a categorical exposure names the effect and
# the comparison together, and the other exposure types name the effect alone.
ipw_degenerate_se_rows <- function(estimates) {
  degenerate <- ipw_degenerate_se(estimates$estimate, estimates$std.err)

  labels <- if (is.null(estimates$comparison)) {
    estimates$effect
  } else {
    paste(estimates$effect, "for", estimates$comparison)
  }

  labels[degenerate]
}

# Report standard errors that collapsed on a fit nothing else has anything to
# say about. The signature was consulted only where the solver had already
# reported a solve it could not pin down, which left the same collapse silent
# wherever the numbers came back cleanly: an outcome that does not vary within
# an exposure arm returns its contrast with an interval of no width and a
# p-value of zero, and said so on both standard error paths without a word.
#
# Raised from the seam each path finishes its estimates at, so it fires once per
# fit. The M-estimation path raises it only where `warn_ipw_unsolved()` does not,
# which is the more specific report of the two and already names these rows.
#
# The wording is neutral about how the numbers were arrived at, because the
# linearization path has no solver and no estimating equations to speak of.
warn_ipw_degenerate_se <- function(estimates, call = rlang::caller_env()) {
  affected <- ipw_degenerate_se_rows(estimates)

  if (length(affected) == 0) {
    return(invisible(FALSE))
  }

  warn(
    c(
      "The standard error{?s} reported for {.val {affected}} {?is/are} not \\
      meaningful.",
      x = "{cli::qty(affected)}{?It is/They are} zero, or so small beside the \\
      estimate{?s} {?it accompanies/they accompany} that the test statistic{?s} \\
      and the interval{?s} built from {?it/them} carry no information.",
      i = "An exposure group the outcome does not vary within is one cause: \\
      the contrast is then a fixed value rather than a quantity with any \\
      spread. Check the outcome within each level of the exposure.",
      i = "The estimates are reported as they were computed."
    ),
    warning_class = "propensity_ipw_degenerate_se_warning",
    call = call
  )

  invisible(TRUE)
}

# Report a solve the equations do not pin down, naming the rows whose standard
# errors came back degenerate. The estimates stay in the output: they are the
# values the solver returned and remain the g-computation point estimates, and
# it is the inference around them that is empty.
#
# The row labels come from the estimates table rather than from the spec, so a
# categorical exposure names the effect and the comparison together and the
# other exposure types name the effect alone.
warn_ipw_unsolved <- function(estimates, call = rlang::caller_env()) {
  affected <- ipw_degenerate_se_rows(estimates)

  degenerate_bullet <- if (length(affected) > 0) {
    c(
      x = "The standard error{?s} reported for {.val {affected}} {?is/are} not \\
      meaningful: {?it collapsed/they collapsed} to essentially zero, which \\
      makes the test statistic{?s} and the interval{?s} built from {?it/them} \\
      meaningless too."
    )
  } else {
    NULL
  }

  warn(
    c(
      "The estimating equations behind these estimates have no unique root at \\
      the values the solver returned.",
      degenerate_bullet,
      i = "At least one direction in the parameter space leaves the equations \\
      unchanged, so the sandwich variance along it is not identified.",
      i = "An exposure level in which every outcome is an event, or none is, \\
      is one cause: the outcome model has no finite fit within it. Check the \\
      outcome within each level of the exposure, and both designs for columns \\
      that duplicate one another.",
      i = "The estimates are reported as the solver returned them."
    ),
    warning_class = "propensity_ipw_solver_warning",
    call = call
  )
}

# The exposure arms the outcome can be looked at within, as a named list of the
# outcome values in each. A binary exposure has the two the contrasts compare, a
# categorical exposure has one per level, and a continuous exposure has none:
# its dose is not a grouping, so there is nothing to look within and this
# derives nothing.
ipw_exposure_arms <- function(spec) {
  y <- spec$outcome$y

  switch(
    spec$exposure_type,
    binary = list(
      "the unexposed group" = y[spec$exposure == 0],
      "the exposed group" = y[spec$exposure == 1]
    ),
    categorical = stats::setNames(
      lapply(
        seq_len(ncol(spec$exposure)),
        function(i) y[spec$exposure[, i] == 1]
      ),
      paste("the", colnames(spec$exposure), "level")
    ),
    NULL
  )
}

# How the outcome is degenerate within an arm, phrased for the coding it is
# read under. A binomial outcome is an indicator, so a constant one is every
# observation an event or none, and any other outcome is described by the single
# value it takes. An arm with no observations in it is degenerate for the same
# reason and is described as the empty arm it is.
ipw_degenerate_arm_label <- function(values, family) {
  if (length(values) == 0) {
    return("holds no observations")
  }

  value <- values[[1]]

  if (identical(family, "binomial") && value %in% c(0, 1)) {
    if (value == 0) {
      return("holds no events")
    }

    return("holds nothing but events")
  }

  paste0("holds the single outcome value ", format(value))
}

# Report a fit whose variance could not be built, in the terms of the fit the
# user handed over rather than of the matrix deli could not invert. Two things
# arrive here. One is an outcome that never varies within an exposure arm, which
# leaves the outcome model with no finite fit inside that arm and the equations
# without a finite derivative at the solution; the fit carries that structure
# and it is read off and named. The other is a design the fits kept every
# coefficient for but that all but repeats itself, whose derivatives are not
# finite at the solution without any arm being degenerate; nothing in the fit
# names that, so the report says what is true of both and where to look.
#
# Only the first degenerate arm is named: the fix is to drop it or to remodel
# the outcome, after which the fit runs again and finds the next.
#
# This is an error rather than a warning because there is nothing to return: the
# standard errors, the intervals, and the p-values all come from the variance,
# and `ipw()` reports estimates with inference or not at all.
abort_ipw_no_variance <- function(spec, call = rlang::caller_env()) {
  arms <- ipw_exposure_arms(spec)
  degenerate <- vapply(
    arms,
    function(values) length(unique(values)) <= 1,
    logical(1)
  )

  if (any(degenerate)) {
    arm <- names(arms)[degenerate][[1]]
    described <- ipw_degenerate_arm_label(
      arms[degenerate][[1]],
      spec$outcome$family
    )
    diagnosis <- c(
      x = "The outcome does not vary within {arm}, which {described}.",
      i = "The outcome model has no finite fit within an arm whose outcome \\
      never varies, so the stacked estimating equations have no finite \\
      derivative at the values the solver returned and there is no variance to \\
      build from them.",
      i = "Drop that arm and estimate the effect among the ones that are left, \\
      or model an outcome that varies within every arm."
    )
  } else {
    diagnosis <- c(
      i = "The stacked estimating equations have no finite derivative at the \\
      values the solver returned, so there is no variance to build from them.",
      i = "Check the outcome within each level of the exposure, and both \\
      designs for columns that all but duplicate one another: a column the fit \\
      kept a coefficient for can still leave the equations flat enough here \\
      that their derivatives are not finite."
    )
  }

  abort(
    c(
      "Can't compute a variance for this fit.",
      diagnosis,
      i = "{.fun ipw} reports estimates with inference or not at all, so no \\
      estimates are returned."
    ),
    error_class = "propensity_ipw_variance_error",
    call = call
  )
}

# Build the estimates table from the solved contrast rows of theta, mirroring
# the column layout that calculate_estimates() produces on the linearization
# path so the two SE methods return the same shape. The contrast rows are
# addressed by their theta positions, not by name: theta names are not unique
# across blocks (a ps or outcome covariate can share a contrast label such as
# "rd" or "diff"), so name subsetting could silently return the wrong row.
# Effect label for a continuous exposure, taken from the outcome-model link:
# the reported effect is the single marginal structural model coefficient.
ipw_continuous_effect_label <- function(link, call = rlang::caller_env()) {
  switch(
    link,
    identity = "slope",
    logit = "log(or)",
    log = "log(rr)",
    abort(
      "Unsupported outcome link {.val {link}} for a continuous exposure.",
      error_class = "propensity_ipw_link_error",
      call = call
    )
  )
}

ipw_mestimation_estimates <- function(
  spec,
  fit,
  layout,
  conf_level,
  call = rlang::caller_env()
) {
  co <- stats::coef(fit)
  vc <- stats::vcov(fit)
  se <- sqrt(diag(vc))

  # A continuous exposure reports the marginal structural model exposure
  # coefficient, addressed by its outcome-block theta position; binary and
  # categorical exposures report the contrast rows. Positions are used, not
  # names, since theta names are not unique across blocks.
  if (identical(spec$exposure_type, "continuous")) {
    idx <- layout$idx$out[spec$outcome$exposure_col]
    effect <- ipw_continuous_effect_label(spec$outcome$link, call = call)
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
    comparison <- rep(
      ipw_comparison_labels(spec),
      each = length(spec$contrasts)
    )
    out <- cbind(
      out["effect"],
      comparison = comparison,
      out[setdiff(names(out), "effect")]
    )
  }

  # The covariance of the reported effects, taken from the same rows and columns
  # of the sandwich the standard errors came from. The effect measures are
  # transformations of one another's marginal means, so they covary, and that
  # part of the fit is recoverable from nowhere else in the result. The block
  # arrives under the theta names, which are not unique across blocks and read as
  # seeds rather than as effects, so it is relabeled the way every surface of the
  # result labels its rows.
  attr(out, "ipw_vcov") <- ipw_vcov_block(vc, idx, ipw_effect_labels(out))

  out
}

# The labels every surface of a result keys its effects by: the effect measure
# on its own, or the effect measure and the comparison together where a
# categorical exposure reports one row per comparison.
ipw_effect_labels <- function(estimates) {
  if (is.null(estimates[["comparison"]])) {
    estimates$effect
  } else {
    paste(estimates$effect, estimates$comparison)
  }
}

# A square block of the joint sandwich, relabeled on both margins. Blocks are
# addressed by position because theta names are not unique across blocks, and
# the labels a block is read under belong to whatever reads it rather than to
# the seed the solver was started from.
ipw_vcov_block <- function(vc, idx, labels) {
  block <- vc[idx, idx, drop = FALSE]
  dimnames(block) <- list(labels, labels)
  block
}

# The block of the joint sandwich belonging to one component model, under the
# names that model calls its own parameters. Those names come from the layout,
# which seeds each block from the fitted model it stands for, rather than from
# the solved fit.
ipw_component_vcov <- function(vc, layout, block) {
  idx <- layout$idx[[block]]
  ipw_vcov_block(vc, idx, names(layout$init)[idx])
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
