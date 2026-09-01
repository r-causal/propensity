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
# Its treatment block is the score its own model solves, read through the same
# registry the single-dose route reads, carrying the conditional variance of the
# density its weight divides by.
#
# Either component may carry a stabilizing numerator, and each one is estimated
# in a block of its own: the marginal proportion or the two marginal moments the
# default stabilizer is, or the score of a numerator model the caller fit. The
# stabilization block is therefore per component and is sliced by the widths the
# spec records, the way the treatment blocks are.
#
# The stacked system carries both treatment models' score blocks, so the weights
# entering the outcome score are recomputed from both blocks of theta on every
# evaluation and the sandwich accounts for having estimated both.

# ---- the method -------------------------------------------------------------

#' @description
#' The `joint_wt_models` method estimates the effects of a joint intervention on
#' two treatments from the pair of fitted treatment models [joint_wt_models()]
#' records and a weighted outcome model that reads both treatments. Standard
#' errors are computed by M-estimation. Neither `se_method = "linearization"`
#' nor `se_method = "robust"` is available here: neither path solves a stacked
#' system, and the cell means and every contrast built from them are parameters
#' of one.
#'
#' The only supported estimand is `"ate"`, which is what the product weights
#' [wt_joint()] builds target.
#'
#' The second treatment may be a dose, in which case the surface is the marginal
#' structural model's own coefficients rather than the cells of a crossing. The
#' dose model is read through the same registry the single-treatment route
#' reads, so an [stats::lm()], a gaussian [stats::glm()] at an identity or a log
#' link, and a [MASS::rlm()] fit with one of the psi functions MASS supplies are
#' stacked, and the dose component's weights may be any density ratio [wt_ate()]
#' builds. A dose model or a density the stacked system cannot differentiate is
#' refused rather than resampled, since this route has no resampling method; see
#' **Joint exposures**.
#'
#' Either component may be stabilized, including on a numerator model the caller
#' fit and passed to [wt_ate()]'s `stabilize`. Each component's numerator is
#' estimated in a block of its own, so the standard errors account for it having
#' been fitted; see **Joint exposures** for what each block holds and for what a
#' numerator may condition on. A numerator the caller computed and passed as
#' [wt_ate()]'s `stabilization_score` is read back off the product and held
#' fixed, contributing no block, as it is for a single treatment.
#'
#' Every score this route stacks is unweighted, so a treatment model fit with
#' prior case weights is refused by the name it was recorded under in
#' [joint_wt_models()].
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
  se_method = c("mestimation", "linearization", "robust"),
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

  # The refusals come before anything is read off the models, so a request this
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
  if (identical(se_method, "mestimation")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.val {se_method}} standard errors for a \\
      joint treatment model.",
      x = "The cell means and every contrast built from them are parameters \\
      of the stacked estimating equations, and the {se_method} path solves \\
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

  # A dose model is read through the registry the single-dose route reads, so
  # which classes, families, and links can be stacked is one answer rather than
  # two, and a fit that cannot be stacked is refused here with the registry's
  # own reason. Only the second component may be a dose, which the check above
  # has already settled.
  #
  # The refusal names the component rather than `wt_mod`, which here is the
  # container the two treatment models arrived in: a reader told to refit
  # `wt_mod` would be told to refit the wrong thing.
  dose_idx <- if (identical(types[[2]], "continuous")) 2L else NULL
  dose_model <- NULL
  if (!is.null(dose_idx)) {
    dose_model <- ipw_continuous_model(
      fits[[dose_idx]],
      hint = ipw_joint_models_dose_hint(),
      label = names[[dose_idx]],
      call = call
    )
    check_ipw_continuous_model(dose_model, call = call)
  }

  # Every score this route stacks is the unweighted one, so a component fit
  # under prior case weights is not at the root of the equation written for it
  # and the solve, seeded at its coefficients, would report a treatment model
  # nobody fit. The single-treatment routes refuse such a fit by name, and this
  # one refuses each component in turn, naming the component rather than the
  # container the models arrived in.
  for (i in seq_along(fits)) {
    check_ipw_joint_model_weights(fits[[i]], names[[i]], call = call)
  }

  # Every cell of the crossing is set at once, so an outcome model reading one
  # of the two treatments has no counterfactual design for three of the four
  # cells. Both omissions are the same fault.
  for (name in names) {
    check_ipw_outcome_exposure(outcome_mod, name, call = call)
  }

  coefs <- lapply(fits, stats::coef)

  for (i in seq_along(fits)) {
    check_ipw_model_rank(coefs[[i]], names[[i]], call = call)
  }

  # Read before the designs are built rather than with the rest of what the
  # weights record: each component's numerator model has its design rebuilt from
  # `.data` alongside them, so the columns those models read have to be among
  # the ones the rebuild asks for before it goes looking for them.
  wts <- extract_weights(outcome_mod)
  numerator_mods <- ipw_joint_models_numerator_models(wts, dose_idx)

  mm_data <- ipw_joint_models_frame(
    outcome_mod,
    .data,
    fits = fits,
    treatment_names = names,
    numerator_mods = numerator_mods,
    call = call
  )
  n <- nrow(mm_data)

  # With a `.data` the caller supplied, every design and every treatment column
  # is built from that frame, over the rows every model read. Without one they
  # are read off the fits, where an additive dose fit's entry has already
  # evaluated the smooth basis its score was checked at; a fit of another kind
  # carries no design on its entry, and one whose design could not be read at
  # all carries none either, so both are built and refused there.
  if (is.null(.data)) {
    ps_X <- lapply(seq_along(fits), function(i) {
      if (identical(i, dose_idx) && !is.null(dose_model$design)) {
        dose_model$design
      } else {
        ipw_joint_models_design(fits[[i]], call = call)
      }
    })
    treatments <- lapply(fits, ipw_joint_models_treatment, call = call)
  } else {
    ps_X <- lapply(seq_along(fits), function(i) {
      design <- ipw_rebuild_design(
        fits[[i]],
        stats::delete.response(stats::terms(fits[[i]])),
        mm_data,
        call = call
      )
      check_ipw_design_width(design, fits[[i]], names[[i]], call = call)

      design
    })
    # Each treatment is its own model's response and is named for it, which is
    # what `joint_wt_models()` requires of the pair, so the column of `.data` to
    # read is the one the component is named after.
    treatments <- lapply(names, function(name) mm_data[[name]])

    for (i in seq_along(fits)) {
      check_ipw_exposure_levels(
        treatments[[i]],
        fits[[i]],
        names[[i]],
        arg = names[[i]],
        call = call
      )
    }
  }

  ipw_joint_models_check_lengths(treatments, ps_X, names, n, call = call)

  # The response each treatment model was fit against and the residuals it left,
  # both over that model's own rows. The moments the weights were built at were
  # read from exactly these, and `.data` can leave fewer rows than a fit was
  # made over, so they are carried on the spec rather than recomputed where the
  # seeds are written.
  fit_moments <- lapply(fits, ipw_joint_models_fit_moments)

  estimand <- check_estimand(wts, estimand, call = call)

  # What ratio of densities the dose's weights are, which is what the stacked
  # system has to rebuild. A product records one density per component, so the
  # dose's record is read out of the joint record positionally; a product built
  # before the record existed has none, and its dose is the ratio the package
  # has always built.
  ratio <- if (!is.null(dose_idx)) {
    ipw_continuous_ratio_meta(
      joint_wt_dose_density(wts, dose_idx),
      stabilized = TRUE,
      hint = ipw_joint_models_dose_hint(),
      call = call
    )
  }

  # Each component's stabilization score, which is the one numerator that exists
  # nowhere but the record: a model can be estimated again and the marginal
  # numerator can be estimated from the exposure, and a score is a vector the
  # caller computed. The system holds it fixed, exactly as the single-treatment
  # routes hold theirs.
  scores <- joint_wt_stabilization_scores(joint_wt_meta(wts))

  # A score the record still holds has to describe the rows about to be
  # weighted. A model frame that drops incomplete rows re-attaches the record
  # whole, leaving a score at the length the product was built at, which is
  # multiplied against a density read over fewer rows unless it is caught first.
  #
  # An empty slot is not this: it is the drop the record marks, and what stands
  # in for it and what reports the difference are settled below.
  for (i in seq_along(scores)) {
    if (is.null(scores[[i]])) {
      next
    }

    check_ipw_stabilization_score(
      "score",
      scores[[i]],
      n,
      component = names[[i]],
      call = call
    )
  }

  # A dose whose weights record a score the record does not keep is the one
  # numerator this route cannot rebuild. A product written before the record
  # kept scores says the numerator was a score and holds no vector, so the
  # system estimates the exposure's own marginal moments instead and the
  # weight-consistency preflight is what reports the difference, naming the
  # score as the cause.
  numerator <- ratio$numerator
  if (identical(numerator, "score") && is.null(scores[[dose_idx]])) {
    numerator <- "marginal"
  }

  # Each component's stabilizing numerator is estimated in a block of its own,
  # the way the single-treatment routes estimate the one numerator they carry,
  # so a component stabilized on a fitted model is stacked rather than refused.
  # The blocks are built here, where the numerator models each go through the
  # guards their own route puts them through, and their designs are rebuilt from
  # `.data` where the caller supplied one, over the rows every other design here
  # is built over.
  stab_components <- ipw_joint_models_stab_components(
    wts,
    types = types,
    names = names,
    n = n,
    numerator = numerator,
    recorded = ratio$numerator,
    numerator_model = ratio$numerator_model,
    scores = scores,
    dose_idx = dose_idx,
    .data = if (!is.null(.data)) mm_data,
    call = call
  )

  if (is_linear_regression(outcome_mod)) {
    family <- "gaussian"
    out_link <- "identity"
    contrasts <- "diff"
  } else {
    family <- "binomial"
    out_link <- outcome_mod$family$link
    contrasts <- c("rd", "log(rr)")
  }

  dose <- !is.null(dose_idx)

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

  # A discrete treatment enters as the 0/1 indicator of its non-reference level,
  # which is the coding its own binomial score is written against, and a dose
  # enters as itself.
  exposures <- lapply(seq_along(treatments), function(i) {
    if (identical(types[[i]], "continuous")) {
      as.double(treatments[[i]])
    } else {
      ipw_recode_binary_exposure(treatments[[i]])
    }
  })

  # A dose model whose spread the system estimates carries the conditional
  # variance of its density in its block of theta, so its block is one wider
  # than the model has coefficients. A spread the caller fixed is a known
  # constant instead, which sits in no block at all.
  widths <- lengths(coefs)
  if (dose && !identical(ratio$sigma$kind, "fixed")) {
    widths[[dose_idx]] <- widths[[dose_idx]] + 1L
  }

  list(
    exposure_type = "joint_models",
    estimand = estimand,
    n = n,
    exposure = exposures,
    ps = list(
      X = ps_X,
      link = vapply(
        seq_along(fits),
        function(i) {
          if (identical(types[[i]], "continuous")) {
            dose_model$link
          } else {
            fits[[i]]$family$link
          }
        },
        character(1)
      ),
      # Per component rather than concatenated, since the seed interleaves each
      # dose model's conditional variance after its own coefficients.
      coefs = coefs,
      widths = widths,
      types = types,
      # The registry entry the dose's block is written from, which the psi and
      # the preflight both rebuild that block out of rather than reaching back
      # for the model itself.
      kind = dose_model$kind,
      psi_loss = dose_model$psi_loss,
      psi_k = dose_model$psi_k,
      penalty = dose_model$penalty,
      # Per component, and over each fit's own rows rather than over the rows
      # analyzed here, which is what the stabilization seeds are read from.
      fit_exposure = lapply(fit_moments, `[[`, "exposure"),
      fit_residuals = lapply(fit_moments, `[[`, "residuals"),
      k = 2L
    ),
    # `wt_joint()` requires a continuous component to be stabilized, so a dose
    # brings a stabilizing numerator with it, and a binary component may bring
    # one of its own. Which numerator a component carries decides what the
    # system estimates for it: a marginal dose numerator is two moments of the
    # exposure, a marginal binary one is a single proportion, a fitted model is
    # its own coefficients, and an integrated one is built from the dose block
    # and the data alone. `widths` says how wide each component's slice of the
    # stabilization block is, so the block is sliced the way `ps$widths` slices
    # the treatment blocks rather than by position.
    stab = list(
      components = stab_components,
      widths = vapply(stab_components, function(x) x$width, integer(1))
    ),
    density = if (dose) ratio$density,
    numerator = if (dose) numerator,
    sigma = if (dose) ratio$sigma,
    # The points an integrated numerator averages the conditional density over,
    # which are a function of the exposure alone and so are fixed across the
    # solve, as they were when `wt_joint()`'s dose component was built.
    grid = if (dose) ipw_numerator_grid(exposures[[dose_idx]], numerator),
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

# One component's prior case weights, refused for the same reason the
# single-treatment routes refuse them: the block written for the component is
# its unweighted score, so a fit made under weights is not at that block's root.
# The refusal names the component, since `wt_mod` here is the container the two
# models arrived in and a reader told to refit it would be told to refit the
# wrong thing.
check_ipw_joint_model_weights <- function(
  fit,
  label,
  call = rlang::caller_env()
) {
  weights <- ipw_model_prior_weights(fit)

  if (is.null(weights) || all(weights == 1)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support a treatment model fit with case weights.",
      x = "{.arg {label}} was fit with non-unit {.arg weights}, so its \\
      coefficients are not the root of the unweighted score stacked for it.",
      i = "Refit {.arg {label}} without {.arg weights}."
    ),
    error_class = "propensity_ipw_ps_weights_error",
    call = call
  )
}

# The remedy the dose refusals point to on this route. The package computes no
# standard error for a fit the stacked system cannot write, here or anywhere
# else, so the refusal says where the standard errors would have to come from
# and leaves the resampling to the user.
ipw_joint_models_dose_hint <- function() {
  "This route builds standard errors from the stacked system alone. Build the
   dose weights from a model and a density that system can differentiate, or
   bootstrap the whole joint fit yourself: resample the rows, refit both
   treatment models, rebuild the weights with {.fn wt_joint}, and refit the
   outcome model on each resample."
}

# The dose component's density record, read positionally off the product. A
# product built before the record existed carries none, which reads as the ratio
# every earlier version of the package built.
joint_wt_dose_density <- function(wts, dose) {
  density <- joint_wt_meta(wts)$density

  if (is.null(density)) {
    return(NULL)
  }

  density[[dose]]
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
      i = "A model with no intercept, written {.code - 1} or {.code + 0}, \\
      expands a factor treatment to an indicator for every level, so its first \\
      column is the reference-level indicator rather than the 0/1 indicator \\
      the rows describe. Keep the intercept, or code {.val {first}} as a 0/1 \\
      numeric.",
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

# The frame every design this route builds is rebuilt from: `.data` when the
# caller supplied one and the outcome model's own frame otherwise, which is
# where both treatment columns already sit.
#
# Every guard the single-treatment routes run over a supplied `.data` runs here,
# read for the models this route holds: the columns the rebuilds need have to be
# there before any of them goes looking, the frame has to describe the rows the
# fits analyzed, the columns have to arrive under the coding the fits recorded,
# and none of them may be missing a value, which `model.frame()` would drop a
# row for while the weights kept it.
ipw_joint_models_frame <- function(
  outcome_mod,
  .data,
  fits,
  treatment_names,
  numerator_mods = NULL,
  call = rlang::caller_env()
) {
  if (is.null(.data)) {
    frame <- tryCatch(model.frame(outcome_mod), error = function(e) e)
    if (inherits(frame, "error")) {
      abort_outcome_frame_gone(conditionMessage(frame), call = call)
    }
    return(frame)
  }

  models <- ipw_joint_models_fit_args(
    fits,
    treatment_names,
    outcome_mod,
    numerator_mods
  )
  required <- ipw_joint_models_required_columns(
    outcome_mod,
    treatment_names,
    models
  )
  assert_columns_exist(.data, required, call = call)

  n_fitted <- nrow(stats::model.frame(outcome_mod))

  # A set of fits made under `na.exclude` analyzed the rows complete over the
  # columns they read, and the frame the caller has is the frame those fits
  # were given. Restricting `.data` to those rows reconciles the two, which is
  # what the single-model routes do with the same helper; a frame that is
  # longer for any other reason has no such restriction to make and keeps the
  # report below.
  .data <- restrict_ipw_data(.data, required, n_fitted)

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

  check_ipw_data_types(.data, models, outcome_mod, call = call)
  check_ipw_data_complete(.data, required, call = call)

  .data
}

# The fits a `.data` guard reads on this route, keyed by the name the caller
# would have to look at. A treatment model arrived under the name of the
# treatment it fits, so that is what names it; a numerator model arrived under
# the `stabilize` argument of the component's own weights, and the components
# name it together, since a column two numerators read is one recoding either
# way.
ipw_joint_models_fit_args <- function(
  fits,
  treatment_names,
  outcome_mod,
  numerator_mods = NULL
) {
  numerators <- numerator_mods[!vapply(numerator_mods, is.null, logical(1))]

  c(
    stats::setNames(fits, treatment_names),
    list(outcome_mod = outcome_mod),
    stats::setNames(numerators, rep("stabilize", length(numerators)))
  )
}

# The columns of `.data` the rebuilds read: both treatments, whatever the
# outcome model's response is computed from, and every covariate any of the fits
# reads, a numerator model's among them.
ipw_joint_models_required_columns <- function(
  outcome_mod,
  treatment_names,
  models
) {
  unique(c(
    treatment_names,
    fmla_extract_left_vars(outcome_mod),
    unlist(lapply(models, ipw_model_covariates))
  ))
}

# The numerator models the components' stabilization blocks are estimated from,
# one entry per component and NULL where a component's numerator is not a fitted
# model. A discrete component's model is recorded on the product itself and a
# dose's on its density record, and the accessor spans both.
#
# What decides whether a component's model is read at all is what the record
# says stabilized that component: a score is the numerator itself, an
# unstabilized component has none, and a dose is read for the model only where
# its density record says a model is what its numerator was.
ipw_joint_models_numerator_models <- function(wts, dose_idx) {
  meta <- joint_wt_meta(wts)
  models <- joint_wt_numerator_models(meta)
  stabilized <- meta$stabilized
  scores <- joint_wt_stabilization_scores(meta)

  lapply(seq_along(models), function(i) {
    if (identical(i, dose_idx)) {
      density <- joint_wt_dose_density(wts, dose_idx)

      return(if (identical(density$numerator, "model")) models[[i]])
    }

    if (isTRUE(stabilized[[i]]) && is.null(scores[[i]])) {
      models[[i]]
    }
  })
}

# The response a treatment model was fit against and the residuals it left, both
# over the model's own rows. The moments the component's stabilization was built
# at were computed from exactly these, and the rows `.data` leaves need not be
# the rows the fit was made over, so the seeds read them from here rather than
# recomputing them over the rows the system goes on to solve over.
#
# Both are read off the fit rather than out of its model frame, which a fit made
# with `model = FALSE` in a function whose frame is gone has none of: the
# response is the conditional mean plus the residual around it, and both are
# kept whatever the frame is. Reading them on the response scale is what puts
# their sum on the treatment's own scale under a link, and it is the 0/1 coding
# for a discrete treatment, whose binomial fit takes its first level as zero.
# A fit made under `na.exclude` pads both out to the length of the frame it was
# given rather than to the rows it analyzed, so the padding is dropped.
ipw_joint_models_fit_moments <- function(fit) {
  residuals <- as.double(stats::residuals(fit, type = "response"))
  response <- as.double(stats::fitted(fit)) + residuals
  kept <- !is.na(residuals)

  list(exposure = response[kept], residuals = residuals[kept])
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
  components = c("binary", "binary"),
  call = rlang::caller_env()
) {
  factors <- lapply(components, function(type) {
    switch(
      type,
      binary = ipw_binary_weight_fn(estimand),
      continuous = ipw_continuous_weight_fn(estimand, call = call)
    )
  })

  function(ps, exposure, extras) {
    factors[[1]](ps[[1]], exposure[[1]], extras[[1]]) *
      factors[[2]](ps[[2]], exposure[[2]], extras[[2]])
  }
}

# What the system estimates for each component's stabilizing numerator, one
# entry per component and in the container's order. Each entry names the kind of
# numerator the component's weights record, the block of a fitted one, and how
# wide the component's slice of the stabilization block is; a component whose
# numerator the system estimates nothing for takes a slice of width zero.
#
# The parameters those slices hold are named `stab_<component>_<parameter>`,
# since a joint system carries two components and a name saying only which role
# a parameter plays would not say whose. A marginal binary numerator is one
# proportion, `pi`; a marginal dose numerator is the exposure's two moments,
# `mu` and `sigma2`; and a fitted numerator is one parameter per column of its
# own design, followed for a dose by `sigma2`, the spread its density is read
# at. `ipw_init_joint_models_stab()` writes those names and is the one place
# they are written.
#
# A binary component's numerator model is recorded on the product itself and a
# dose's on its density record, so the two are read from the places their own
# routes read them from and put through the guards their own routes apply. A
# score is read from the product either way, since a vector the caller computed
# belongs to no model and to no density.
ipw_joint_models_stab_components <- function(
  wts,
  types,
  names,
  n,
  numerator,
  recorded,
  numerator_model,
  scores,
  dose_idx,
  .data = NULL,
  call = rlang::caller_env()
) {
  meta <- joint_wt_meta(wts)
  stabilized <- meta$stabilized
  binary_models <- joint_wt_numerator_models(meta)
  # Whether the record could have kept a score at all. A product written before
  # the slot existed says a discrete component was stabilized without saying
  # what by, so the marginal proportion stands in for a score of the caller's
  # there and nowhere else.
  keeps_scores <- joint_wt_records_scores(meta)
  # Whether a component's score was dropped by an operation that changed the
  # length of the weights. The slot such a drop leaves is the slot a component
  # stabilized on the marginal proportion has, so without the mark the marginal
  # proportion stands in for a score of the caller's silently, and a mismatch
  # names no component.
  dropped_scores <- joint_wt_dropped_scores(meta)

  lapply(seq_along(types), function(i) {
    if (identical(i, dose_idx)) {
      return(ipw_joint_models_dose_stab(
        numerator,
        recorded,
        numerator_model,
        scores[[i]],
        names[[i]],
        n,
        .data = .data,
        call = call
      ))
    }

    # A discrete component the product records as unstabilized carries no
    # numerator, and one it records as stabilized carries a score of the
    # caller's, a model of the caller's, or the marginal proportion the default
    # stabilizer estimates.
    if (!isTRUE(stabilized[[i]])) {
      return(list(
        type = types[[i]],
        numerator = "none",
        model = NULL,
        score = NULL,
        width = 0L,
        stand_in = FALSE
      ))
    }

    # A score is the numerator itself rather than a description of one, so the
    # system multiplies by it and estimates nothing: no block, no parameters,
    # exactly as the single-treatment routes treat a score.
    if (!is.null(scores[[i]])) {
      return(list(
        type = types[[i]],
        numerator = "score",
        model = NULL,
        score = scores[[i]],
        width = 0L,
        stand_in = FALSE
      ))
    }

    model <- ipw_binary_numerator_block(
      binary_models[[i]],
      .data = .data,
      call = call
    )

    if (is.null(model)) {
      return(list(
        type = types[[i]],
        numerator = "marginal",
        model = NULL,
        score = NULL,
        width = 1L,
        # A record that keeps scores says the marginal proportion is what
        # stabilized this component, unless it also records that the component's
        # score was dropped. One written before it kept them says only that
        # something did.
        stand_in = !keeps_scores || isTRUE(dropped_scores[[i]])
      ))
    }

    check_ipw_numerator_model(
      binary_models[[i]],
      model,
      names[[i]],
      n,
      call = call
    )

    list(
      type = types[[i]],
      numerator = "model",
      model = model,
      score = NULL,
      width = ncol(model$X),
      stand_in = FALSE
    )
  })
}

# The dose component's entry. A marginal numerator is the exposure's own two
# moments, a fitted one is its coefficients and the spread its density is read
# at, and an integrated one is built from the dose block and the data alone, so
# it estimates nothing and takes no slice.
ipw_joint_models_dose_stab <- function(
  numerator,
  recorded,
  numerator_model,
  score,
  name,
  n,
  .data = NULL,
  call = rlang::caller_env()
) {
  if (identical(numerator, "marginal")) {
    return(list(
      type = "continuous",
      numerator = "marginal",
      model = NULL,
      score = NULL,
      width = 2L,
      # A dose whose weights record a score the product does not keep is
      # rebuilt from the marginal moments instead, and that substitution is a
      # cause the weights preflight reports by name.
      stand_in = identical(recorded, "score")
    ))
  }

  # A score is a known multiplier, so the ratio is rebuilt with the vector the
  # product recorded and the system estimates nothing for it.
  if (identical(numerator, "score")) {
    return(list(
      type = "continuous",
      numerator = "score",
      model = NULL,
      score = score,
      width = 0L,
      stand_in = FALSE
    ))
  }

  if (!identical(numerator, "model")) {
    return(list(
      type = "continuous",
      numerator = "none",
      model = NULL,
      score = NULL,
      width = 0L,
      stand_in = FALSE
    ))
  }

  model <- ipw_numerator_model_block(
    numerator_model,
    .data = .data,
    call = call
  )
  check_ipw_numerator_model(numerator_model, model, name, n, call = call)

  list(
    type = "continuous",
    numerator = "model",
    model = model,
    score = NULL,
    width = ncol(model$X) + 1L,
    stand_in = FALSE
  )
}

# One component's slice of the stabilization block, cut by the widths the spec
# records rather than by position, so a component that estimates nothing takes
# nothing and the component after it still reads its own parameters.
ipw_joint_models_stab_slices <- function(spec, th_stab) {
  widths <- spec$stab$widths
  ends <- cumsum(widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  lapply(seq_along(widths), function(i) {
    if (identical(widths[[i]], 0L)) {
      numeric(0)
    } else {
      th_stab[starts[[i]]:ends[[i]]]
    }
  })
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
  dose <- ipw_joint_models_dose(spec)
  ps_fns <- if (!is.null(dose)) ipw_joint_models_dose_fns(spec)
  stab_th <- ipw_joint_models_stab_slices(spec, th_stab)

  lapply(seq_along(spec$ps$X), function(i) {
    th <- th_ps[starts[[i]]:ends[[i]]]
    x <- spec$ps$X[[i]]
    stab <- spec$stab$components[[i]]

    if (!identical(spec$ps$types[[i]], "continuous")) {
      return(list(
        type = spec$ps$types[[i]],
        coefs = th,
        ps = ipw_inv_link(spec$ps$link[[i]])(as.vector(x %*% th)),
        # A discrete component's numerator is the probability of the level
        # each unit took, which the registry multiplies its unstabilized weight
        # by, or a score the caller computed, which is that numerator itself
        # and takes the probability's place. A component with no numerator
        # takes neither.
        extras = list(
          stab_prob = ipw_binary_stab_prob(stab$model, stab_th[[i]]),
          score = stab$score
        )
      ))
    }

    alpha <- th[seq_len(ncol(x))]

    # A spread the caller fixed is a known constant, and a constant sits in no
    # block of theta: the dose block is its coefficients alone and the density
    # is read at the number the weights record.
    sigma2_d <- if (identical(spec$sigma$kind, "fixed")) {
      spec$sigma$value^2
    } else {
      th[[ncol(x) + 1L]]
    }

    # A dose stabilized on a fitted model reads its own conditional mean and its
    # own spread out of its slice, where a marginal numerator reads the
    # exposure's two moments.
    th_n <- stab_th[[i]]
    marginal <- identical(stab$numerator, "marginal")
    model <- if (identical(stab$numerator, "model")) stab$model

    list(
      type = "continuous",
      coefs = alpha,
      sigma2_d = sigma2_d,
      ps = ps_fns$mean(x, alpha),
      extras = list(
        sigma2_d = sigma2_d,
        mu_a = if (marginal) th_n[[1]],
        sigma2_a = if (marginal) th_n[[2]],
        mu_n = if (!is.null(model)) {
          ipw_numerator_model_fns(model)$mean(
            model$X,
            th_n[seq_len(ncol(model$X))]
          )
        },
        sigma_n = if (!is.null(model)) sqrt(th_n[[ncol(model$X) + 1L]]),
        score = stab$score,
        stabilized = !identical(stab$numerator, "none"),
        # Which ratio the dose's weights are, so the factor rebuilt here is the
        # factor `wt_joint()` multiplied rather than the one this route used to
        # assume every dose was.
        density = spec$density,
        numerator = spec$numerator,
        grid = spec$grid
      )
    )
  })
}

# The mean and the score the dose model contributes, read off the spec's record
# of its registry entry rather than off the model itself, so the psi that
# rebuilds the block at every evaluation and the preflight that rebuilds it once
# write the same equation from the same five values.
ipw_joint_models_dose_fns <- function(spec) {
  dose <- ipw_joint_models_dose(spec)

  ipw_continuous_score_fns(
    spec$ps$kind,
    spec$ps$link[[dose]],
    spec$ps$psi_loss,
    spec$ps$psi_k,
    spec$ps$penalty
  )
}

# One treatment model's score rows: the unweighted score its own fit sits at,
# which for a dose is the equation its registry entry names rather than least
# squares whatever the fit was. A dose model adds the row that estimates the
# conditional variance its density ratio divides by, which is the mean squared
# residual read against that entry's conditional mean; a dose whose spread the
# caller fixed carries no such row, because a constant has no equation.
ipw_joint_models_score_rows <- function(
  block,
  x,
  y,
  link,
  ps_fns = NULL,
  sigma_row = TRUE,
  sigma = NULL,
  density = NULL
) {
  if (!identical(block$type, "continuous")) {
    return(deli::ee_glm(
      block$coefs,
      X = x,
      y = y,
      distribution = "binomial",
      link = link
    ))
  }

  score <- ps_fns$score(block$coefs, x, y)

  if (!sigma_row) {
    return(score)
  }

  rbind(
    score,
    ipw_continuous_sigma_row(sigma, y - block$ps, block$sigma2_d, density)
  )
}

# The stabilizing numerators' rows, one block per component and in the
# components' own order, each seeded at the exact root of the equation written
# for it. A marginal dose numerator is the two moments of the exposure's own
# density; a marginal binary one is the proportion that took the non-reference
# level; a fitted one is the score its own class solves, followed for a dose by
# the moment its spread is the root of. A component the system estimates no
# numerator for contributes no rows at all.
ipw_joint_models_stab_rows <- function(spec, th_stab, exposures) {
  slices <- ipw_joint_models_stab_slices(spec, th_stab)

  rows <- lapply(seq_along(slices), function(i) {
    th <- slices[[i]]

    if (!length(th)) {
      return(NULL)
    }

    stab <- spec$stab$components[[i]]
    a <- exposures[[i]]

    if (identical(stab$numerator, "marginal")) {
      if (identical(stab$type, "continuous")) {
        return(deli::ee_mean_variance(th, y = a))
      }

      return(matrix(a - th[[1]], nrow = 1))
    }

    model <- stab$model
    p_n <- ncol(model$X)

    if (!identical(stab$type, "continuous")) {
      return(deli::ee_glm(
        th,
        X = model$X,
        y = a,
        distribution = "binomial",
        link = model$link
      ))
    }

    fns <- ipw_numerator_model_fns(model)
    mu_n <- fns$mean(model$X, th[seq_len(p_n)])

    rbind(
      fns$score(th[seq_len(p_n)], model$X, a),
      matrix((a - mu_n)^2 - th[[p_n + 1L]], nrow = 1)
    )
  })

  do.call(rbind, rows)
}

# The index of the dose among the components, or NULL where both are discrete.
# Only one component may be a dose; `check_ipw_joint_models_types()` settles
# that before anything here reads it.
ipw_joint_models_dose <- function(spec) {
  dose <- which(spec$ps$types == "continuous")

  if (length(dose)) dose[[1]] else NULL
}

# The response a component's treatment model was fit against and the residuals
# it left, over that model's own rows, which is where the moments the weights
# were built at were read. `ipw_spec_joint_models()` carries both, because the
# rows `.data` leaves need not be the rows a fit was made over. A spec assembled
# without that record is one whose rows are the fits', so the rows it holds
# answer the same question.
ipw_joint_models_fit_exposure <- function(spec, i) {
  recorded <- spec$ps$fit_exposure

  if (!is.null(recorded)) {
    return(recorded[[i]])
  }

  spec$exposure[[i]]
}

ipw_joint_models_fit_residuals <- function(spec, i, ps_fns) {
  recorded <- spec$ps$fit_residuals

  if (!is.null(recorded)) {
    return(recorded[[i]])
  }

  spec$exposure[[i]] - ps_fns$mean(spec$ps$X[[i]], spec$ps$coefs[[i]])
}

# The treatment blocks' seed: each model's coefficients, and for a dose whose
# spread the system estimates, the conditional variance of its density, read
# from that model's own residuals and the exact root of the row that estimates
# it. A dose whose spread the caller fixed seeds its coefficients alone, since a
# constant is in no block.
#
# Those residuals are the fit's own rather than the ones the analyzed rows
# leave. The weights were built at the fit's moment, and a seed at any other
# moment rebuilds weights the caller was never given.
ipw_init_joint_models_ps <- function(spec, call = rlang::caller_env()) {
  dose <- ipw_joint_models_dose(spec)
  ps_fns <- if (!is.null(dose)) ipw_joint_models_dose_fns(spec)

  blocks <- lapply(seq_along(spec$ps$coefs), function(i) {
    alpha <- spec$ps$coefs[[i]]

    if (!identical(spec$ps$types[[i]], "continuous")) {
      return(alpha)
    }

    if (identical(spec$sigma$kind, "fixed")) {
      return(alpha)
    }

    resid <- ipw_joint_models_fit_residuals(spec, i, ps_fns)
    c(
      alpha,
      stats::setNames(
        ipw_continuous_sigma2_seed(
          spec$sigma,
          resid,
          spec$density,
          call = call
        ),
        paste0("sigma2_", spec$names[[i]])
      )
    )
  })

  unlist(unname(blocks))
}

# The stabilizing numerators' seed, one block per component and in the
# components' own order, each at the exact root of the row that estimates it.
# Parameters are named `stab_<component>_<parameter>`, since a joint system
# carries two components and a name saying only which role a parameter plays
# would not say whose. This is where that convention is written; every other
# reader of the stabilization block slices it by the widths the spec records.
#
# Every block is seeded at moments the fit it belongs to came to, taken over
# that fit's rows, since those are the moments the weights were built at. The
# rows the spec analyzes are what the equations solve the moments over, which is
# the answer rather than the starting value.
ipw_init_joint_models_stab <- function(spec) {
  blocks <- lapply(seq_along(spec$stab$components), function(i) {
    stab <- spec$stab$components[[i]]

    if (identical(stab$width, 0L)) {
      return(numeric(0))
    }

    prefix <- paste0("stab_", spec$names[[i]], "_")

    if (identical(stab$numerator, "marginal")) {
      a_fit <- ipw_joint_models_fit_exposure(spec, i)

      if (identical(stab$type, "continuous")) {
        mu_a <- mean(a_fit)
        return(stats::setNames(
          c(mu_a, mean((a_fit - mu_a)^2)),
          paste0(prefix, c("mu", "sigma2"))
        ))
      }

      return(stats::setNames(
        ipw_default_stab_seed(a_fit),
        paste0(prefix, "pi")
      ))
    }

    model <- stab$model
    coefs <- stats::setNames(model$coefs, paste0(prefix, colnames(model$X)))

    if (!identical(stab$type, "continuous")) {
      return(coefs)
    }

    c(
      coefs,
      stats::setNames(
        ipw_numerator_fit_sigma2(model, spec$exposure[[i]]),
        paste0(prefix, "sigma2")
      )
    )
  })

  out <- unlist(unname(blocks))

  if (is.null(out)) numeric(0) else out
}

ipw_init_joint_models <- function(spec, call = rlang::caller_env()) {
  beta <- spec$outcome$coefs
  ps_block <- ipw_init_joint_models_ps(spec, call = call)
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
  ps_fns <- if (!is.null(dose)) ipw_joint_models_dose_fns(spec)
  sigma_row <- !identical(spec$sigma$kind, "fixed")

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
          ps_link[[i]],
          ps_fns = ps_fns,
          sigma_row = sigma_row,
          sigma = spec$sigma,
          density = spec$density
        )
      })
    )

    stab_rows <- ipw_joint_models_stab_rows(spec, th_stab, exposure)

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
