# The two-model route to a joint intervention.
#
# `joint_exposure()` crosses two treatments into one categorical exposure and
# weights it with one multinomial propensity score model. This route weights the
# same intervention through the sequential factorization it really has,
# f(A | L) f(E | A, L): two treatment models, one per treatment, and the product
# of their weights. `joint_wt_models()` records the pair, `wt_joint()` builds the
# product, and everything here estimates over that container.
#
# What is reported is the surface the declared route reports, row for row and
# label for label, because it is built by the same machinery: the plan comes
# from `ipw_joint_plan()` and the rows from `ipw_joint_estimate_rows()`, both in
# R/ipw-joint.R. The crossing the plan is built from is constructed here from
# the two treatment columns rather than read off a declaration, which is what
# makes the two routes agree about which cells exist and what they are called
# without either of them owning the answer.
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

  # Only a pair of binary treatments is estimated here. The container accepts a
  # categorical or a continuous component, whose score block and weight are not
  # the binomial ones this stack carries.
  check_ipw_joint_models_types(models$exposure_type, call = call)

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

  # The crossing is built from the two treatment columns rather than declared,
  # and it is built with the constructor the declared route reads, so the cells
  # this route reports are the cells that route reports: same labels, same
  # order, same reference. Construction also refuses an unpopulated cell, which
  # is a positivity violation the joint effect is not identified under.
  declared <- do.call(
    causalgenerics::joint_exposure,
    stats::setNames(treatments, names)
  )
  cells <- levels(declared)

  if (is_linear_regression(outcome_mod)) {
    family <- "gaussian"
    out_link <- "identity"
    contrasts <- "diff"
  } else {
    family <- "binomial"
    out_link <- outcome_mod$family$link
    contrasts <- c("rd", "log(rr)")
  }

  joint <- ipw_joint_plan(declared, contrasts)

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

  list(
    exposure_type = "joint_models",
    estimand = estimand,
    n = n,
    # The two treatments as 0/1 indicators of their non-reference level, which
    # is the coding each model's own binomial score is written against.
    exposure = lapply(treatments, ipw_recode_binary_exposure),
    ps = list(
      X = ps_X,
      link = vapply(fits, function(fit) fit$family$link, character(1)),
      coefs = unlist(unname(coefs)),
      widths = lengths(coefs),
      k = 2L
    ),
    stab = list(stabilized = FALSE, score = NULL),
    outcome = list(
      X = model.matrix(outcome_mod),
      y = ipw_outcome_numeric(fmla_extract_left_vctr(outcome_mod)),
      family = family,
      link = out_link,
      coefs = stats::coef(outcome_mod),
      X_counterfactual = x_cf,
      weights = as.double(wts)
    ),
    contrasts = contrasts,
    by = NULL,
    joint = joint,
    names = names,
    focal_level = NULL,
    reference_level = cells[[1]]
  )
}

check_ipw_joint_models_types <- function(
  exposure_type,
  call = rlang::caller_env()
) {
  unsupported <- exposure_type != "binary"

  if (!any(unsupported)) {
    return(invisible(TRUE))
  }

  bad <- names(exposure_type)[unsupported]
  bad_type <- unname(exposure_type[unsupported])

  abort(
    c(
      "{.fun ipw} currently supports a joint intervention on two binary \\
      treatments.",
      x = "The model named {.arg {bad[[1]]}} fits a {.val {bad_type[[1]]}} \\
      treatment.",
      i = "The stacked system carries a binomial score for each treatment, \\
      which is not the score a categorical or a continuous treatment model \\
      sits at.",
      i = "Cross two binary treatments, or report the two treatments \\
      separately."
    ),
    error_class = "propensity_ipw_exposure_error",
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

# The product of the two treatments' inverse-probability contributions, which is
# what `wt_joint()` computed at the observed fit and what this recomputes as a
# function of theta. Each factor is the binary ate weight the registry already
# defines, so the product the solver sees at the seed is the product the caller
# built, to the bit.
ipw_joint_models_weight_fn <- function(estimand) {
  binary <- ipw_binary_weight_fn(estimand)

  function(ps, exposure, extras) {
    binary(ps[[1]], exposure[[1]], list()) *
      binary(ps[[2]], exposure[[2]], list())
  }
}

# The fitted probabilities each treatment model implies at its block of theta.
ipw_joint_models_ps <- function(spec, th_ps) {
  ends <- cumsum(spec$ps$widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  lapply(seq_along(spec$ps$X), function(i) {
    block <- th_ps[starts[[i]]:ends[[i]]]
    ipw_inv_link(spec$ps$link[[i]])(as.vector(spec$ps$X[[i]] %*% block))
  })
}

ipw_init_joint_models <- function(spec, call = rlang::caller_env()) {
  beta <- spec$outcome$coefs
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
    ps = spec$ps$coefs,
    stab = numeric(0),
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
  widths <- spec$ps$widths
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
  ends <- cumsum(widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    # One score block per treatment model, each the unweighted binomial score
    # its own fit sits at. Stacked in the container's order, so the second
    # model's block carries the coefficients on the first treatment and the
    # sandwich accounts for both fits.
    blocks <- lapply(seq_along(ps_X), function(i) {
      th_i <- th_ps[starts[[i]]:ends[[i]]]
      list(
        rows = deli::ee_glm(
          th_i,
          X = ps_X[[i]],
          y = exposure[[i]],
          distribution = "binomial",
          link = ps_link[[i]]
        ),
        ps = ipw_inv_link(ps_link[[i]])(as.vector(ps_X[[i]] %*% th_i))
      )
    })

    ps_rows <- do.call(rbind, lapply(blocks, function(b) b$rows))

    # The weights are rebuilt from both propensity score blocks on every
    # evaluation, which is what propagates the uncertainty of having estimated
    # them into the sandwich.
    w <- weight_fn(lapply(blocks, function(b) b$ps), exposure, list())

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    # The ate tilt is the constant one, so each cell mean is the ordinary mean
    # of its counterfactual predictions and no tilt multiplies these rows.
    preds <- lapply(x_cf, function(x) inv_out(as.vector(x %*% th_out)))
    mu_rows <- t(vapply(
      seq_len(k),
      function(j) preds[[j]] - th_mu[[j]],
      numeric(n)
    ))

    con_rows <- ipw_joint_psi_rows(joint, th_mu, th_con, n, reporters)

    ipw_stack(list(ps_rows, out_rows, mu_rows, con_rows))
  }
}
