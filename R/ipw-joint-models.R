# The two-model route to a joint intervention.
#
# `joint_exposure()` crosses two treatments into one categorical exposure and
# weights it with one multinomial propensity score model. This route weights the
# same intervention through the sequential factorization it really has,
# f(A | L) f(E | A, L): two treatment models, one per treatment, and the product
# of their weights. `joint_wt_models()` records the pair, `wt_joint()` builds the
# product, and everything here estimates over that container.
#
# For a pair of discrete treatments, what is reported is the surface the
# declared route reports, row for row and label for label, because it is built
# by the same machinery: the plan comes from `ipw_joint_plan()` and the rows
# from `ipw_joint_estimate_rows()`, both in R/ipw-joint.R. The crossing the plan
# is built from is constructed here from the two treatment columns rather than
# read off a declaration, which is what makes the two routes agree about which
# cells exist and what they are called without either of them owning the answer.
#
# The second treatment may instead be a dose. A dose has no cells, so there is
# no crossing to construct and nothing to set every unit to: the surface is
# coefficient-shaped rather than cell-shaped, reporting the marginal structural
# model's own causal coefficients under labels written in the same vocabulary.
# Its treatment block is a least-squares score carrying the conditional variance
# of the density its weight divides by, and the stabilizing numerator that
# density ratio needs is estimated alongside it.
#
# The stacked system carries both treatment models' score blocks, so the weights
# entering the outcome score are recomputed from both blocks of theta on every
# evaluation and the sandwich accounts for having estimated both.

# ---- the method -------------------------------------------------------------

#' @description
#' The `joint_wt_models` method estimates the effects of a joint intervention on
#' two treatments from the pair of fitted treatment models [joint_wt_models()]
#' records and a weighted outcome model that reads both treatments. Standard
#' errors are computed by M-estimation; the linearization method is not
#' available. The only supported estimand is `"ate"`, which is what the product
#' weights [wt_joint()] builds target.
#'
#' The second treatment may be a dose, in which case the surface is the marginal
#' structural model's own coefficients rather than the cells of a crossing; see
#' **Joint exposures**.
#'
#' This is the second of the two routes to a joint intervention. The other
#' declares the crossing with [causalgenerics::joint_exposure()] and weights it
#' with one multinomial model over the cells; see **Joint exposures**. Both
#' report the same surface, so the choice between them is a modeling choice
#' rather than a reporting one. Prefer this route when the two treatments call
#' for different adjustment sets, or when the dependence of the second treatment
#' on the first is what you want to model directly: it weights through the
#' sequential factorization `f(A | L) f(E | A, L)`, so each treatment gets its
#' own model and its own covariates. Prefer the declared crossing when one model
#' over the cells is the natural specification.
#'
#' The factorization is validated when the container is built rather than here:
#' [joint_wt_models()] refuses a second model that does not condition on the
#' first treatment, and a pair that condition on each other. By the time a
#' container reaches `ipw()` those questions are settled.
#'
#' @name ipw-methods
#' @exportS3Method causalgenerics::ipw joint_wt_models
ipw.joint_wt_models <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization"),
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)
  effects <- rlang::arg_match(effects)
  assert_class(outcome_mod, c("glm", "lm"))

  # The propensity score models are the container's, and each carries the link
  # it was fit with, so there is nothing here for `ps_link` to override.
  check_ipw_ps_link_absent(ps_link, "joint treatment")

  # Both refusals come before anything is read off the models, so a request this
  # route will not answer is refused on its own terms. `.by` is refused before
  # any modifier could be resolved, and so before the outcome model could be
  # diagnosed for a term a refused request would never use.
  check_ipw_joint_models_method(se_method)
  check_ipw_joint_by(
    TRUE,
    .by,
    remedy = "Drop {.arg .by} to report the joint surface, or weight the \\
    crossing of the two treatments as one plain categorical exposure to report \\
    each cell against the reference cell within each subgroup."
  )

  # Guards on the weights that fit the outcome model, mirroring the other
  # methods: the psw attributes say whether these are the product weights the
  # container names, which nothing downstream would otherwise notice.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts)
  check_ipw_joint_models_weights(wts)

  spec <- ipw_spec_joint_models(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand
  )
  fit <- ipw_mestimation(spec, conf_level = conf_level)

  new_ipw(
    estimand = spec$estimand,
    wt_mod = wt_mod,
    # The container is not a model, so only the outcome model is wrapped with
    # its block of the joint sandwich. That block is what the conditional
    # reading reports, and it accounts for both treatment models having been
    # estimated from the same data.
    outcome_mod = new_ipw_model(outcome_mod, fit$outcome_vcov),
    estimates = fit$estimates,
    se_method = "mestimation",
    fit = fit$fit,
    effects = effects
  )
}

check_ipw_joint_models_method <- function(
  se_method,
  call = rlang::caller_env()
) {
  if (!identical(se_method, "linearization")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.val linearization} standard errors for a \\
      joint treatment model.",
      x = "The cell means and every contrast built from them are parameters \\
      of the stacked estimating equations, and the linearization path solves \\
      no such system.",
      i = "Use {.code se_method = \"mestimation\"} for a joint treatment \\
      model."
    ),
    error_class = "propensity_ipw_joint_models_method_error",
    call = call
  )
}

check_ipw_joint_models_weights <- function(wts, call = rlang::caller_env()) {
  if (is_joint_wt(wts)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg outcome_mod} must be fit with the product weights the two \\
      treatment models imply.",
      x = "Its weights carry no record of being a product, so they weight one \\
      treatment rather than the crossing of the two.",
      i = "A single treatment's weights are an ordinary {.cls psw} of the \\
      right length and the right estimand, so nothing else would notice.",
      i = "Build the weights with {.fun wt_joint} from the two components, and \\
      refit {.arg outcome_mod} with them."
    ),
    error_class = "propensity_ipw_joint_models_weights_error",
    call = call
  )
}

# ---- the spec ---------------------------------------------------------------

ipw_spec_joint_models <- function(
  models,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  call = rlang::caller_env()
) {
  assert_class(outcome_mod, c("glm", "lm"))
  check_ipw_offset(outcome_mod, call = call)
  check_ipw_outcome_response(outcome_mod, call = call)
  check_ipw_outcome_family(outcome_mod, call = call)

  names <- models$names
  fits <- models$models
  types <- unname(models$exposure_type)

  # The container accepts a categorical component too, whose score block and
  # weight are neither the binomial nor the density-ratio ones this stack
  # carries.
  check_ipw_joint_models_types(models$exposure_type, call = call)
  check_ipw_joint_models_dose_link(fits, types, names, call = call)

  # Every cell of the crossing is set at once, so an outcome model reading one
  # of the two treatments has no counterfactual design for three of the four
  # cells. Both omissions are the same fault.
  for (name in names) {
    check_ipw_outcome_exposure(outcome_mod, name, call = call)
  }

  ps_X <- lapply(fits, ipw_joint_models_design, call = call)
  treatments <- lapply(fits, ipw_joint_models_treatment, call = call)
  coefs <- lapply(fits, stats::coef)

  for (i in seq_along(fits)) {
    check_ipw_model_rank(coefs[[i]], names[[i]], call = call)
  }

  mm_data <- ipw_joint_models_frame(outcome_mod, .data, call = call)
  n <- nrow(mm_data)
  ipw_joint_models_check_lengths(treatments, ps_X, names, n, call = call)

  wts <- extract_weights(outcome_mod)
  estimand <- check_estimand(wts, estimand, call = call)

  if (is_linear_regression(outcome_mod)) {
    family <- "gaussian"
    out_link <- "identity"
    contrasts <- "diff"
  } else {
    family <- "binomial"
    out_link <- outcome_mod$family$link
    contrasts <- c("rd", "log(rr)")
  }

  dose <- identical(types[[2]], "continuous")

  if (dose) {
    # A dose has no cells, so there is no crossing to build, nothing to set
    # every unit to, and no counterfactual risks to report. The surface is the
    # marginal structural model's own coefficients and this says which of them
    # carry causal content.
    cells <- NULL
    joint <- NULL
    x_cf <- NULL
    coefficients <- ipw_joint_dose_rows(
      outcome_mod,
      names,
      treatments,
      out_link,
      call = call
    )
  } else {
    # The crossing is built from the two treatment columns rather than declared,
    # and it is built with the constructor the declared route reads, so the
    # cells this route reports are the cells that route reports: same labels,
    # same order, same reference. Construction also refuses an unpopulated cell,
    # which is a positivity violation the joint effect is not identified under.
    declared <- do.call(
      causalgenerics::joint_exposure,
      stats::setNames(treatments, names)
    )
    cells <- levels(declared)
    joint <- ipw_joint_plan(declared, contrasts)
    coefficients <- NULL

    out_terms <- stats::delete.response(stats::terms(outcome_mod))
    x_cf <- ipw_joint_models_designs(
      outcome_mod,
      out_terms,
      mm_data,
      names,
      treatments,
      cells,
      rebuilt = !is.null(.data),
      call = call
    )

    check_ipw_outcome_design_width(x_cf, outcome_mod, call = call)
    check_ipw_counterfactual_designs(
      x_cf,
      cells,
      paste(names, collapse = " and "),
      call = call
    )
  }

  list(
    exposure_type = "joint_models",
    estimand = estimand,
    n = n,
    # A discrete treatment enters as the 0/1 indicator of its non-reference
    # level, which is the coding its own binomial score is written against, and
    # a dose enters as itself.
    exposure = lapply(seq_along(treatments), function(i) {
      if (identical(types[[i]], "continuous")) {
        as.double(treatments[[i]])
      } else {
        ipw_recode_binary_exposure(treatments[[i]])
      }
    }),
    ps = list(
      X = ps_X,
      link = vapply(
        seq_along(fits),
        function(i) ipw_joint_models_link(fits[[i]], types[[i]]),
        character(1)
      ),
      # Per component rather than concatenated, since the seed interleaves each
      # dose model's conditional variance after its own coefficients.
      coefs = coefs,
      # A dose model carries that variance in its block of theta, so its block
      # is one wider than the model has coefficients.
      widths = lengths(coefs) + (types == "continuous"),
      types = types,
      k = 2L
    ),
    # `wt_joint()` requires a continuous component to be stabilized, so a dose
    # brings a stabilizing numerator with it and the system estimates the two
    # marginal moments that numerator is built from. A component built with a
    # fixed stabilization score is a different function of the data, and the
    # weight-consistency preflight is what notices.
    stab = list(stabilized = dose, score = NULL),
    outcome = list(
      X = model.matrix(outcome_mod),
      y = ipw_outcome_numeric(fmla_extract_left_vctr(outcome_mod)),
      family = family,
      link = out_link,
      coefs = stats::coef(outcome_mod),
      X_counterfactual = x_cf,
      weights = as.double(wts)
    ),
    contrasts = if (!dose) contrasts,
    coefficients = coefficients,
    by = NULL,
    joint = joint,
    names = names,
    focal_level = NULL,
    reference_level = if (!dose) cells[[1]]
  )
}

# The link a treatment model's score block is written against. A dose model is
# fit by least squares, whose score is the identity-link one whether the model
# arrived as an `lm` or as a gaussian `glm`; a discrete treatment model carries
# the link its own family was fit with.
ipw_joint_models_link <- function(fit, type) {
  if (identical(type, "continuous")) {
    return("identity")
  }

  fit$family$link
}

# The pairs of treatment types this route estimates: two binary treatments, or a
# binary treatment and a dose. The position matters for the dose, which is
# supported as the second treatment, the one whose model conditions on the
# first: a dose is what the second factor of the factorization may be, and
# nothing here carries the density of a first factor that is one.
ipw_joint_models_supported <- list(
  first = "binary",
  second = c("binary", "continuous")
)

check_ipw_joint_models_types <- function(
  exposure_type,
  call = rlang::caller_env()
) {
  unsupported <- c(
    !exposure_type[[1]] %in% ipw_joint_models_supported$first,
    !exposure_type[[2]] %in% ipw_joint_models_supported$second
  )

  if (!any(unsupported)) {
    return(invisible(TRUE))
  }

  position <- c("first", "second")[unsupported][[1]]
  bad <- names(exposure_type)[unsupported][[1]]
  bad_type <- unname(exposure_type[unsupported])[[1]]

  abort(
    c(
      "{.fun ipw} currently supports a joint intervention on two binary \\
      treatments, or on a binary treatment and a dose.",
      x = "The model named {.arg {bad}} fits a {.val {bad_type}} treatment as \\
      the {position} of the two.",
      i = "The stacked system carries a binomial score for a binary treatment \\
      and a linear score with a conditional variance for a dose. A categorical \\
      treatment model sits at neither, and a dose is carried as the second \\
      treatment alone.",
      i = "Cross two binary treatments, weight a dose as the second treatment, \\
      or report the two treatments separately."
    ),
    error_class = "propensity_ipw_exposure_error",
    call = call
  )
}

# A dose model's block is the least-squares score of the exposure on its
# covariates, which is the score an identity-link fit sits at. A gaussian model
# fit through another link has the same coefficients under a different
# parameterization, and the block would put the seed off the root.
check_ipw_joint_models_dose_link <- function(
  fits,
  types,
  names,
  call = rlang::caller_env()
) {
  bad <- which(vapply(
    seq_along(fits),
    function(i) {
      identical(types[[i]], "continuous") &&
        inherits(fits[[i]], "glm") &&
        !identical(fits[[i]]$family$link, "identity")
    },
    logical(1)
  ))

  if (length(bad) == 0) {
    return(invisible(TRUE))
  }

  name <- names[[bad[[1]]]]
  link <- fits[[bad[[1]]]]$family$link

  abort(
    c(
      "{.fun ipw} supports only an identity-link model of a dose.",
      x = "The model named {.arg {name}} is a gaussian model with a \\
      {.val {link}} link.",
      i = "Refit it as an {.fun lm} or a gaussian glm with an identity link."
    ),
    error_class = "propensity_ipw_link_error",
    call = call
  )
}

# ---- the surfaces a dose reports --------------------------------------------

# Which coefficients of a joint marginal structural model carry causal content,
# and how each of them is named. A dose has no cells, so nothing is standardized
# and what is reported either way are the weighted fit's own coefficients. Which
# rows those are is one question, answered here for both surfaces, and how they
# are named is the other, answered by whichever surface the model reaches.
#
# A row is reported when the term its column belongs to reads a treatment. A
# covariate term reads neither and contributes no row however many columns it
# expands to, which is the whole of the rule and is applied once, here, so that
# neither surface can get it differently. It is applied to columns rather than
# to terms because a term is not what a row reports.
ipw_joint_dose_rows <- function(
  outcome_mod,
  names,
  treatments,
  link,
  call = rlang::caller_env()
) {
  out_X <- model.matrix(outcome_mod)
  term_of_column <- attr(out_X, "assign")
  term_labels <- attr(stats::terms(outcome_mod), "term.labels")
  term_vars <- lapply(term_labels, function(label) all.vars(str2lang(label)))

  reads <- lapply(names, function(name) {
    vapply(term_vars, function(vars) name %in% vars, logical(1))
  })
  reads_treatment <- reads[[1]] | reads[[2]]
  reads_covariate <- vapply(
    term_vars,
    function(vars) length(setdiff(vars, names)) > 0,
    logical(1)
  )

  check_ipw_joint_dose_terms(
    term_labels[reads_treatment & reads_covariate],
    names,
    call = call
  )

  # The one filter. Columns come out in the design's order, so the rows read as
  # the coefficient vector reads, and each carries the term it belongs to, which
  # is what the labels below are decided from.
  columns <- which(term_of_column %in% which(reads_treatment))
  column_terms <- term_of_column[columns]

  # A model in bare treatment terms is one the vocabulary describes, and it is
  # reported in that vocabulary. Every other model reaches the coefficient
  # surface. Nothing else decides this: a fit reports the surface its own
  # marginal structural model has a reading for.
  bare <- ipw_joint_dose_bare_terms(names)

  if (all(term_labels[reads_treatment] %in% bare$admitted)) {
    return(ipw_joint_dose_vocabulary_rows(
      out_X,
      columns,
      column_terms,
      reads,
      term_labels,
      names,
      treatments,
      link,
      call = call
    ))
  }

  ipw_joint_dose_coefficient_rows(out_X, columns, link, call = call)
}

# The vocabulary surface, which a model linear in each treatment reports. Each
# row is named the way the declared route names its rows: `contrast` names the
# treatment being varied and how, and `group` says where in the other
# treatment's range the row is evaluated. The binary treatment's coefficient is
# its effect at a dose of zero, the dose's is its slope at the binary
# treatment's reference level, and the interaction row keeps the discrete
# route's idiom of a group comparing the other treatment's levels, which for a
# dose is a one-unit step. An additive model evaluates neither row anywhere in
# particular, so both are reported as overall.
#
# A bare term of a binary treatment or of a dose contributes exactly one column,
# so the columns and the terms they belong to are in step here and each row is
# one coefficient.
ipw_joint_dose_vocabulary_rows <- function(
  out_X,
  columns,
  column_terms,
  reads,
  term_labels,
  names,
  treatments,
  link,
  call = rlang::caller_env()
) {
  level_labels <- as.character(ipw_joint_models_level_values(treatments[[1]]))

  # A bare term says which variables a column is built from and nothing about
  # how it is coded, so the columns themselves are read here. Each reported row
  # is a claim about a coefficient, and a claim about a coefficient is a claim
  # about the column it multiplies.
  check_ipw_joint_dose_coding(
    out_X,
    columns,
    ipw_joint_dose_columns(reads, column_terms, treatments),
    term_labels[column_terms],
    names,
    level_labels,
    call = call
  )

  interacted <- any(reads[[1]] & reads[[2]])
  level_contrast <- ipw_joint_contrast_label(names[[1]], level_labels, 2L)
  level_effect <- ipw_effect_label(link, "diff", call = call)
  dose_contrast <- paste0(names[[2]], ": per unit")
  dose_effect <- ipw_effect_label(link, "slope", call = call)

  named <- lapply(column_terms, function(term) {
    if (reads[[1]][[term]] && reads[[2]][[term]]) {
      return(c(
        level_effect,
        level_contrast,
        paste0(names[[2]], " + 1 vs ", names[[2]])
      ))
    }

    if (reads[[1]][[term]]) {
      return(c(
        level_effect,
        level_contrast,
        if (interacted) paste0(names[[2]], " = 0") else ipw_overall_group
      ))
    }

    c(
      dose_effect,
      dose_contrast,
      if (interacted) {
        paste0(names[[1]], " = ", level_labels[[1]])
      } else {
        ipw_overall_group
      }
    )
  })

  ipw_coefficient_rows(
    col = columns,
    effect = vapply(named, function(row) row[[1]], character(1)),
    contrast = vapply(named, function(row) row[[2]], character(1)),
    group = vapply(named, function(row) row[[3]], character(1))
  )
}

# The coefficient surface, which every other treatment-reading model reports.
# One row per treatment-reading coefficient, named by the coefficient as the fit
# names it, and no group column at all: where the vocabulary surface says which
# treatment a row varies and where it is evaluated, this one says only which
# coefficient it is, because for a curve there is no one place a row holds at.
# An interpretable dose response is built downstream from `coef()` and `vcov()`,
# which is what the covariance between these rows is for.
#
# Naming rows by the coefficient is also what lets the surface hold a model the
# vocabulary has nowhere to put, such as one carrying two distinct interactions
# between the same pair of treatments: the columns have different names, so the
# rows do.
ipw_joint_dose_coefficient_rows <- function(
  out_X,
  columns,
  link,
  call = rlang::caller_env()
) {
  ipw_coefficient_rows(
    col = columns,
    effect = rep(ipw_effect_label(link, "coef", call = call), length(columns)),
    contrast = colnames(out_X)[columns]
  )
}

# The treatment terms the vocabulary surface describes: each treatment on its
# own and their two-way interaction, written bare. The interaction is carried
# under both orders, since `terms()` labels it in the order the formula
# introduces the treatments and either order is the same term. A model whose
# treatment terms all come from this set has a reading in the vocabulary; any
# other treatment-reading model is reported by the coefficient surface, which
# is a choice of surface rather than a refusal.
ipw_joint_dose_bare_terms <- function(names) {
  list(
    admitted = c(
      names,
      paste(names, collapse = ":"),
      paste(rev(names), collapse = ":")
    )
  )
}

# Refuse a term reading a treatment and a covariate together. This is the one
# shape neither surface can report. A coefficient of such a term is the change
# in an effect per unit of the covariate, so the effect it modifies is defined
# only at a value of that covariate no row names, and naming the row after the
# column would report a quantity the table pins to nothing.
#
# What a covariate may do is enter on its own. Such a term reads no treatment,
# contributes no row, and adjusts the marginal structural model without
# appearing on either surface. The declared route carries treatment-by-covariate
# terms because it standardizes over the covariates rather than reporting
# coefficients, which is the difference between the two.
check_ipw_joint_dose_terms <- function(
  mixed_terms,
  names,
  call = rlang::caller_env()
) {
  if (length(mixed_terms) == 0) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} requires a joint marginal structural model whose treatment \\
      terms read the treatments alone.",
      x = "{.code {mixed_terms}} in {.arg outcome_mod} {?is/are} not \\
      {?a term/terms} of {.val {names}} alone.",
      i = "A coefficient of a term reading a covariate is a change in an effect \\
      per unit of that covariate, so no row could name the effect it stands for.",
      i = "A covariate entering on its own is admitted and contributes no row, \\
      so adjust {.arg outcome_mod} for it that way, or cross two discrete \\
      treatments to report a surface that standardizes over the covariates."
    ),
    error_class = "propensity_ipw_msm_error",
    call = call
  )
}

# The column each reported term has to contribute for the row it lands in to be
# true: the indicator of the first treatment's non-reference level, the dose
# itself, or their product. The indicator is built with the recode the weight
# machinery uses, so what is compared is the coding the estimator assumes
# against the coding the outcome model was fit under.
ipw_joint_dose_columns <- function(reads, column_terms, treatments) {
  indicator <- ipw_recode_binary_exposure(treatments[[1]])
  dose <- as.double(treatments[[2]])

  lapply(column_terms, function(term) {
    if (reads[[1]][[term]] && reads[[2]][[term]]) {
      return(indicator * dose)
    }

    if (reads[[1]][[term]]) {
      return(as.double(indicator))
    }

    dose
  })
}

# Refuse an outcome model whose treatment columns are coded some other way.
#
# The bare-term boundary reads term labels, which say which variables a column
# is built from and nothing about how the column is built. A factor treatment
# carries its contrasts on the vector or takes them from the session, and a
# coding other than treatment contrasts leaves the term label untouched while
# rescaling or recentering the column: an ordered two-level factor gives the
# contrast divided by the square root of two, and a sum coding gives half of it
# with the sign reversed, each with an interaction column to match. The fit runs
# to completion either way and reports numbers that are not the effects the rows
# name.
#
# Reading the columns rather than the contrast attributes is what makes this
# hold for codings not considered here, including one a user writes by hand and
# attaches to the factor. The comparison is exact: an indicator is an exact 0/1
# double and an interaction column is one multiplication of the same two
# vectors, so a column that is the one the row describes agrees to the last bit
# and anything else is a different quantity rather than a rounding of the same
# one.
check_ipw_joint_dose_coding <- function(
  out_X,
  columns,
  expected,
  term_labels,
  names,
  level_labels,
  call = rlang::caller_env()
) {
  mismatched <- !vapply(
    seq_along(columns),
    function(i) {
      identical(unname(as.double(out_X[, columns[[i]]])), expected[[i]])
    },
    logical(1)
  )

  if (!any(mismatched)) {
    return(invisible(TRUE))
  }

  bad_terms <- term_labels[mismatched]
  first <- names[[1]]
  second <- names[[2]]
  reference <- level_labels[[1]]
  focal <- level_labels[[2]]

  abort(
    c(
      "{.fun ipw} reports a joint intervention with a dose from a marginal \\
      structural model whose treatment columns are the treatments themselves.",
      x = "{.code {bad_terms}} in {.arg outcome_mod} {?contributes/contribute} \\
      a column coded some other way.",
      i = "The reported rows name the coefficients of a model in which \\
      {.val {first}} enters as 0 for {.val {reference}} and 1 for \\
      {.val {focal}}, {.val {second}} enters as itself, and their interaction \\
      is the product of the two.",
      i = "A contrast coding other than treatment contrasts rescales or \\
      recenters those columns without changing what the formula says. An \\
      ordered factor carries polynomial contrasts, and \\
      {.code options(contrasts = )} sets a coding for every factor in the \\
      session.",
      i = "Refit {.arg outcome_mod} with {.val {first}} as a 0/1 numeric, or as \\
      an unordered factor under treatment contrasts."
    ),
    error_class = "propensity_ipw_msm_error",
    call = call
  )
}

# A treatment model's design, with the guided error the other spec constructors
# give when a model has lost the data behind its fitting call.
ipw_joint_models_design <- function(fit, call = rlang::caller_env()) {
  design <- tryCatch(model.matrix(fit), error = function(e) e)

  if (inherits(design, "error")) {
    abort(
      c(
        "Can't reconstruct the data behind a treatment model.",
        x = "{conditionMessage(design)}",
        i = "Refit the treatment models where the data they were fit to is \\
        still available."
      ),
      error_class = "propensity_ipw_data_error",
      call = call
    )
  }

  design
}

ipw_joint_models_treatment <- function(fit, call = rlang::caller_env()) {
  values <- tryCatch(fmla_extract_left_vctr(fit), error = function(e) e)

  if (inherits(values, "error")) {
    abort(
      c(
        "Can't read the treatment behind a treatment model.",
        x = "{conditionMessage(values)}",
        i = "Refit the treatment models where the data they were fit to is \\
        still available."
      ),
      error_class = "propensity_ipw_data_error",
      call = call
    )
  }

  values
}

# The frame the counterfactual designs are rebuilt from: `.data` when the caller
# supplied one and the outcome model's own frame otherwise, which is where both
# treatment columns already sit.
ipw_joint_models_frame <- function(
  outcome_mod,
  .data,
  call = rlang::caller_env()
) {
  if (is.null(.data)) {
    frame <- tryCatch(model.frame(outcome_mod), error = function(e) e)
    if (inherits(frame, "error")) {
      abort_outcome_frame_gone(conditionMessage(frame), call = call)
    }
    return(frame)
  }

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

  .data
}

# Everything the stack multiplies is sized to the outcome model's observations,
# and a treatment model fit to a different subset would be recycled against them
# with nothing signaled.
ipw_joint_models_check_lengths <- function(
  treatments,
  ps_X,
  names,
  n,
  call = rlang::caller_env()
) {
  sizes <- lengths(treatments)
  rows <- vapply(ps_X, nrow, integer(1))
  bad <- which(sizes != n | rows != n)

  if (length(bad) == 0) {
    return(invisible(TRUE))
  }

  first <- bad[[1]]
  bad_name <- names[[first]]
  bad_n <- sizes[[first]]

  abort(
    c(
      "Every model must be fit to the same observations.",
      x = "The model named {.arg {bad_name}} was fit to {bad_n} \\
      observation{?s}.",
      x = "{.arg outcome_mod} was fit to {n} observations.",
      i = "Fit both treatment models and {.arg outcome_mod} to one data frame \\
      with no missing values in the columns they read."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# One counterfactual design per cell, built by setting both treatment columns at
# once. `cells` arrives in the crossing's own order, first treatment varying
# fastest, and `expand.grid()` varies its first argument fastest too, so the
# designs come back in that order.
#
# The values written are the treatments' own, in level order, rather than the
# cell label's character halves: a numeric treatment set to `"0"` would be
# levelled by `model.frame()` into a factor of one level, and a factor treatment
# has to keep every level it was fit with.
ipw_joint_models_designs <- function(
  outcome_mod,
  out_terms,
  mm_data,
  names,
  treatments,
  cells,
  rebuilt,
  call = rlang::caller_env()
) {
  values <- lapply(treatments, ipw_joint_models_level_values)
  grid <- expand.grid(
    seq_along(values[[1]]),
    seq_along(values[[2]]),
    KEEP.OUT.ATTRS = FALSE
  )

  designs <- lapply(seq_len(nrow(grid)), function(cell) {
    d <- mm_data
    d[[names[[1]]]] <- values[[1]][[grid[[1]][[cell]]]]
    d[[names[[2]]]] <- values[[2]][[grid[[2]][[cell]]]]
    ipw_counterfactual_design(
      outcome_mod,
      out_terms,
      d,
      rebuilt = rebuilt,
      call = call
    )
  })

  stats::setNames(designs, cells)
}

# A treatment's distinct values in the order the crossing levels them, carrying
# the type the outcome model was fit on. A factor keeps every level so its
# design coding is the one the coefficients were fit under.
ipw_joint_models_level_values <- function(values) {
  if (is.factor(values)) {
    return(factor(levels(values), levels = levels(values)))
  }

  sort(unique(values))
}

# ---- the stacked system -----------------------------------------------------

# The product of the two treatments' contributions, which is what `wt_joint()`
# computed at the observed fit and what this recomputes as a function of theta.
# Each factor is the weight the registry already defines for the kind of
# treatment it weights, so the product the solver sees at the seed is the
# product the caller built, to the bit.
ipw_joint_models_weight_fn <- function(
  estimand,
  components = c("binary", "binary")
) {
  factors <- lapply(components, function(type) {
    switch(
      type,
      binary = ipw_binary_weight_fn(estimand),
      continuous = ipw_continuous_weight_fn(estimand)
    )
  })

  function(ps, exposure, extras) {
    factors[[1]](ps[[1]], exposure[[1]], extras[[1]]) *
      factors[[2]](ps[[2]], exposure[[2]], extras[[2]])
  }
}

# Each treatment model's block of theta, resolved into the pieces the rest of
# the system reads: the parameters its own score sits at, the fitted quantity
# its weight factor divides by, and whatever else that factor needs.
#
# A dose model's block carries the conditional variance of its density after its
# coefficients, and its weight factor reads the stabilizing numerator's moments
# out of the stab block, neither of which a binary component has. Both routes
# read the pair from here rather than slicing theta themselves, so the preflight
# that rebuilds the weights once and the psi that rebuilds them on every
# evaluation cannot disagree about where a parameter sits.
ipw_joint_models_blocks <- function(spec, th_ps, th_stab) {
  ends <- cumsum(spec$ps$widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  lapply(seq_along(spec$ps$X), function(i) {
    th <- th_ps[starts[[i]]:ends[[i]]]
    x <- spec$ps$X[[i]]

    if (!identical(spec$ps$types[[i]], "continuous")) {
      return(list(
        type = spec$ps$types[[i]],
        coefs = th,
        ps = ipw_inv_link(spec$ps$link[[i]])(as.vector(x %*% th)),
        extras = list()
      ))
    }

    alpha <- th[seq_len(ncol(x))]
    sigma2_d <- th[[ncol(x) + 1L]]

    list(
      type = "continuous",
      coefs = alpha,
      sigma2_d = sigma2_d,
      ps = as.vector(x %*% alpha),
      extras = list(
        sigma2_d = sigma2_d,
        mu_a = if (length(th_stab)) th_stab[[1]],
        sigma2_a = if (length(th_stab)) th_stab[[2]],
        score = spec$stab$score,
        stabilized = spec$stab$stabilized
      )
    )
  })
}

# One treatment model's score rows: the unweighted score its own fit sits at. A
# dose model adds the row that estimates the conditional variance its density
# ratio divides by, which is the mean squared residual of that same fit.
ipw_joint_models_score_rows <- function(block, x, y, link) {
  if (!identical(block$type, "continuous")) {
    return(deli::ee_glm(
      block$coefs,
      X = x,
      y = y,
      distribution = "binomial",
      link = link
    ))
  }

  rbind(
    deli::ee_regression(block$coefs, X = x, y = y, model = "linear"),
    matrix((y - block$ps)^2 - block$sigma2_d, nrow = 1)
  )
}

# The stabilizing numerator's rows: the two moments of the dose's own marginal
# normal density, each the exact root of the row that estimates it.
ipw_joint_models_stab_rows <- function(th_stab, dose) {
  if (!length(th_stab)) {
    return(NULL)
  }

  rbind(
    matrix(dose - th_stab[[1]], nrow = 1),
    matrix((dose - th_stab[[1]])^2 - th_stab[[2]], nrow = 1)
  )
}

# The index of the dose among the components, or NULL where both are discrete.
# Only one component may be a dose; `check_ipw_joint_models_types()` settles
# that before anything here reads it.
ipw_joint_models_dose <- function(spec) {
  dose <- which(spec$ps$types == "continuous")

  if (length(dose)) dose[[1]] else NULL
}

# The treatment blocks' seed: each model's coefficients, and for a dose the
# conditional variance of its density, which is the mean squared residual of its
# own fit and the exact root of the row that estimates it.
ipw_init_joint_models_ps <- function(spec) {
  blocks <- lapply(seq_along(spec$ps$coefs), function(i) {
    alpha <- spec$ps$coefs[[i]]

    if (!identical(spec$ps$types[[i]], "continuous")) {
      return(alpha)
    }

    resid <- spec$exposure[[i]] - as.vector(spec$ps$X[[i]] %*% alpha)
    c(
      alpha,
      stats::setNames(mean(resid^2), paste0("sigma2_", spec$names[[i]]))
    )
  })

  unlist(unname(blocks))
}

# The stabilizing numerator's seed: the dose's own marginal moments. Blocks are
# named for the treatment they belong to, since a joint system carries two
# treatment models and a name saying only which role a parameter plays would
# not say whose.
ipw_init_joint_models_stab <- function(spec) {
  dose <- ipw_joint_models_dose(spec)

  if (is.null(dose) || !spec$stab$stabilized || !is.null(spec$stab$score)) {
    return(numeric(0))
  }

  a <- spec$exposure[[dose]]
  mu_a <- mean(a)

  stats::setNames(
    c(mu_a, mean((a - mu_a)^2)),
    paste0(c("stab_mu_", "stab_sigma2_"), spec$names[[dose]])
  )
}

ipw_init_joint_models <- function(spec, call = rlang::caller_env()) {
  beta <- spec$outcome$coefs
  ps_block <- ipw_init_joint_models_ps(spec)
  stab_block <- ipw_init_joint_models_stab(spec)

  # A dose has no cells, so the surface reports the outcome block directly and
  # the system carries no marginal means and no contrasts over them.
  if (is.null(spec$joint)) {
    return(list(
      ps = ps_block,
      stab = stab_block,
      out = beta,
      mu = numeric(0),
      contrast = numeric(0),
      by_mu = numeric(0),
      by_contrast = numeric(0)
    ))
  }

  inv_out <- ipw_inv_link(spec$outcome$link, call = call)

  # The ate tilt is the constant one, so each cell mean is seeded at the
  # ordinary mean of its counterfactual predictions, which is the exact root of
  # its psi row.
  preds <- lapply(
    spec$outcome$X_counterfactual,
    function(x) inv_out(as.vector(x %*% beta))
  )
  mu_vals <- vapply(preds, mean, numeric(1))
  mu_block <- stats::setNames(mu_vals, paste0("mu_", spec$joint$cells))

  list(
    ps = ps_block,
    stab = stab_block,
    out = beta,
    mu = mu_block,
    contrast = ipw_init_joint(spec$joint, mu_vals),
    by_mu = numeric(0),
    by_contrast = numeric(0)
  )
}

ipw_psi_joint_models <- function(
  spec,
  layout,
  weight_fn,
  call = rlang::caller_env()
) {
  idx <- layout$idx
  ps_X <- spec$ps$X
  ps_link <- spec$ps$link
  exposure <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  inv_out <- ipw_inv_link(out_link, call = call)
  x_cf <- spec$outcome$X_counterfactual
  k <- length(x_cf)
  n <- spec$n
  joint <- spec$joint
  reporters <- ipw_joint_reporters(joint)
  dose <- ipw_joint_models_dose(spec)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    # One score block per treatment model, each the unweighted score its own fit
    # sits at. Stacked in the container's order, so the second model's block
    # carries the coefficients on the first treatment and the sandwich accounts
    # for both fits.
    blocks <- ipw_joint_models_blocks(spec, th_ps, th_stab)
    ps_rows <- do.call(
      rbind,
      lapply(seq_along(blocks), function(i) {
        ipw_joint_models_score_rows(
          blocks[[i]],
          ps_X[[i]],
          exposure[[i]],
          ps_link[[i]]
        )
      })
    )

    stab_rows <- if (is.null(dose)) {
      NULL
    } else {
      ipw_joint_models_stab_rows(th_stab, exposure[[dose]])
    }

    # The weights are rebuilt from both propensity score blocks on every
    # evaluation, which is what propagates the uncertainty of having estimated
    # them into the sandwich.
    w <- weight_fn(
      lapply(blocks, function(block) block$ps),
      exposure,
      lapply(blocks, function(block) block$extras)
    )

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    # A dose has no cells, so the outcome block is the whole of what the surface
    # reports and there is nothing standardized to stack behind it.
    if (is.null(joint)) {
      return(ipw_stack(list(ps_rows, stab_rows, out_rows)))
    }

    # The ate tilt is the constant one, so each cell mean is the ordinary mean
    # of its counterfactual predictions and no tilt multiplies these rows.
    preds <- lapply(x_cf, function(x) inv_out(as.vector(x %*% th_out)))
    mu_rows <- t(vapply(
      seq_len(k),
      function(j) preds[[j]] - th_mu[[j]],
      numeric(n)
    ))

    con_rows <- ipw_joint_psi_rows(joint, th_mu, th_con, n, reporters)

    ipw_stack(list(ps_rows, stab_rows, out_rows, mu_rows, con_rows))
  }
}
