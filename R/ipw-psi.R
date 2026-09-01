# Internal building blocks for the M-estimation path of ipw(). Three targets:
# ipw_weight_fn() (the weight registry), ipw_theta_layout() (the theta
# partition and its root-valued seed), and build_ipw_psi() (the stacked
# estimating-function builder). Everything here is unexported and follows the
# M-estimation design contract: the psi matrix stacks blocks in the fixed order
# [ps score | stabilization | outcome score | marginal means | contrasts |
# stratum means | stratum contrasts], and the weights that enter the
# outcome-score block are recomputed from the ps-parameter block of theta on
# every evaluation so the sandwich variance accounts for propensity score
# estimation.
#
# The last two blocks are the ones a `.by` request adds, and they come last for
# a reason the result depends on: no earlier equation reads a parameter of
# theirs, so the bread is block lower triangular across that partition and the
# leading block of the sandwich is the one a fit without `.by` solves. That is
# what lets a grouped fit report the same whole-sample rows, and the same
# conditional reading, as the fit it grew out of.

# Inverse link functions. Covers the propensity score links (logit, probit,
# cloglog) and the outcome links (logit, identity, log) used when reconstructing
# fitted values from a coefficient block. These are the mathematical inverses,
# not the family's clamped C routines, so a linear predictor extreme enough to
# saturate one of them returns an exact 0 or 1 here where `linkinv` would not.
#
# The registry is a list rather than a `switch()` so that a caller can ask
# whether a link is covered without triggering the error for one that is not.
# The saturation guard on the linearization path needs exactly that.
ipw_inv_links <- list(
  logit = stats::plogis,
  probit = stats::pnorm,
  cloglog = function(x) 1 - exp(-exp(x)),
  identity = function(x) x,
  log = exp
)

ipw_inv_link <- function(link, call = rlang::caller_env()) {
  inv_link <- ipw_inv_links[[link]]

  if (is.null(inv_link)) {
    abort(
      "Unsupported link {.val {link}}.",
      error_class = "propensity_ipw_link_error",
      call = call
    )
  }

  inv_link
}

# ---- weight registry --------------------------------------------------------

# Every estimand `ipw()` knows a weight for. Which of them a particular fit
# supports is narrower and depends on the exposure type, which `ipw_weight_fn()`
# resolves below; this is the set an estimand has to belong to before that
# question is worth asking, and `check_estimand()` reads it for that.
ipw_estimands <- c("ate", "att", "atu", "atm", "ato", "entropy")

ipw_weight_fn <- function(
  exposure_type,
  estimand,
  components = NULL,
  call = rlang::caller_env()
) {
  supported <- switch(
    exposure_type,
    binary = ipw_estimands,
    categorical = ipw_estimands,
    continuous = "ate",
    joint_models = "ate",
    abort(
      "Unsupported exposure type {.val {exposure_type}}.",
      error_class = "propensity_ipw_estimand_error",
      call = call
    )
  )

  if (!estimand %in% supported) {
    abort(
      c(
        "{.fun ipw} does not support the {.val {estimand}} estimand for \\
        {.val {exposure_type}} exposures.",
        i = "Supported estimands for {.val {exposure_type}} exposures: \\
        {.val {supported}}."
      ),
      error_class = "propensity_ipw_estimand_error",
      call = call
    )
  }

  switch(
    exposure_type,
    binary = ipw_binary_weight_fn(estimand),
    categorical = ipw_categorical_weight_fn(estimand),
    continuous = ipw_continuous_weight_fn(estimand, call = call),
    joint_models = ipw_joint_models_weight_fn(estimand, components, call = call)
  )
}

# Binary weight registry. `ps` is the vector of propensity scores e_i,
# `exposure` is the 0/1 exposure z, and `extras` carries the stabilization
# probability, the fixed stabilization score, and optionally a tilt already
# evaluated at `ps`. Formulas mirror the binary helpers in R/weights.R exactly.
# Every estimand is the tilt h(e) from R/ps_tilt.R over the exposure denominator
# z e + (1 - z)(1 - e); ate additionally carries the stabilization branches.
#
# A psi evaluation needs the tilt for the marginal-mean rows as well, so it
# passes the one it already holds through `extras$tilt`. Callers with no tilt in
# hand leave it out and the registry evaluates its own.
ipw_binary_weight_fn <- function(estimand) {
  if (identical(estimand, "ate")) {
    return(function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      if (!is.null(extras$score)) {
        wt <- (z / e) + ((1 - z) / (1 - e))
        wt * extras$score
      } else if (!is.null(extras$stab_prob)) {
        p1 <- extras$stab_prob
        p0 <- 1 - p1
        (z * p1 / e) + ((1 - z) * p0 / (1 - e))
      } else {
        (z / e) + ((1 - z) / (1 - e))
      }
    })
  }

  function(ps, exposure, extras) {
    z <- exposure
    e <- ps
    h_e <- if (is.null(extras$tilt)) {
      ps_tilt_binary(e, estimand)
    } else {
      extras$tilt
    }
    h_e / (z * e + (1 - z) * (1 - e))
  }
}

# Categorical weight registry. `ps` is the n-by-K propensity score matrix,
# `exposure` the n-by-K reference-first indicator, and `extras` carries the
# focal column index, the stabilization probabilities in column order, either
# length K for a numerator every unit shares or an n-by-K matrix for one a
# model gives each unit, the fixed stabilization score, and optionally a tilt
# already evaluated at `ps`, on the same terms as the binary registry. Formulas
# mirror calculate_categorical_weights() in R/weights.R.
ipw_categorical_weight_fn <- function(estimand) {
  is_ate <- identical(estimand, "ate")

  function(ps, exposure, extras) {
    e_actual <- rowSums(exposure * ps)

    # The ate tilt is the constant 1, which divides out of the weight.
    if (!is_ate) {
      h_e <- if (is.null(extras$tilt)) {
        ps_tilt_categorical(ps, estimand, extras$focal_idx)
      } else {
        extras$tilt
      }
      return(h_e / e_actual)
    }

    weights <- 1 / e_actual

    if (!is.null(extras$score)) {
      weights <- weights * extras$score
    } else if (!is.null(extras$stab_probs)) {
      # A numerator that conditions on nothing is one probability per level, the
      # same row for every unit, while a numerator model gives each unit its
      # own. The matrix is read as it stands: broadcasting it would fill each
      # row with the whole matrix read column-wise and gather every unit at
      # another unit's numerator.
      stab_row <- if (is.matrix(extras$stab_probs)) {
        extras$stab_probs
      } else {
        matrix(
          extras$stab_probs,
          nrow = nrow(ps),
          ncol = ncol(ps),
          byrow = TRUE
        )
      }
      weights <- weights * rowSums(exposure * stab_row)
    }

    weights
  }
}

# Continuous weight registry (ate only). `ps` is the fitted conditional mean,
# `exposure` the continuous A, and `extras` carries the conditional variance
# sigma2_d, the density the ratio is taken in, the numerator that stabilized it,
# the marginal moments mu_a and sigma2_a, the numerator model's own mean mu_n
# and spread sigma_n, the evaluation grid an integrated numerator marginalizes
# over, the fixed stabilization score, and the stabilized flag.
#
# The ratio itself is `continuous_density_ratio()`, the same function
# `ate_continuous()` builds the weights with, so the weights rebuilt at a value
# of theta are the weights the user was given whenever the parameters agree,
# rather than a second implementation that has to be kept in step with the
# first. The registry never emits the alert the unstabilized branch of
# `wt_ate()` does.
#
# Weights that record no density took the normal family, and their numerator is
# read off the stabilization the way `ate_continuous()` reads it: a score the
# caller fixed, the marginal density, or nothing at all.
ipw_continuous_weight_fn <- function(estimand, call = rlang::caller_env()) {
  function(ps, exposure, extras) {
    density <- extras$density
    if (is.null(density)) {
      density <- dens_normal()
    }

    numerator <- extras$numerator
    if (is.null(numerator)) {
      numerator <- ipw_continuous_numerator(extras$stabilized, extras$score)
    }

    continuous_density_ratio(
      exposure = exposure,
      mu = ps,
      sigma = sqrt(extras$sigma2_d),
      density = density,
      numerator = numerator,
      mu_a = extras$mu_a,
      sigma_a = if (!is.null(extras$sigma2_a)) sqrt(extras$sigma2_a),
      mu_n = extras$mu_n,
      sigma_n = extras$sigma_n,
      score = extras$score,
      grid = extras$grid,
      # The frame the user entered `ipw()` on. A conditional density that comes
      # out at zero part way through the solve is refused by the ratio, and
      # without this the refusal would name the frame the registry rebuilds the
      # weights in, which the caller never wrote.
      call = call
    )
  }
}

# What stabilized a set of weights that left no record of it: a score the caller
# supplied, the marginal density of the exposure, or nothing.
ipw_continuous_numerator <- function(stabilized, score) {
  if (!isTRUE(stabilized)) {
    return("none")
  }

  if (!is.null(score)) {
    return("score")
  }

  "marginal"
}

# ---- theta layout -----------------------------------------------------------

ipw_theta_layout <- function(spec, call = rlang::caller_env()) {
  blocks <- switch(
    spec$exposure_type,
    binary = ipw_init_binary(spec, call = call),
    categorical = ipw_init_categorical(spec, call = call),
    continuous = ipw_init_continuous(spec, call = call),
    joint_models = ipw_init_joint_models(spec, call = call)
  )

  block_order <- c(
    "ps",
    "stab",
    "out",
    "mu",
    "contrast",
    "by_mu",
    "by_contrast"
  )
  idx <- vector("list", length(block_order))
  names(idx) <- block_order
  pos <- 0L
  for (b in block_order) {
    n_b <- length(blocks[[b]])
    idx[[b]] <- if (n_b > 0L) pos + seq_len(n_b) else integer(0)
    pos <- pos + n_b
  }

  init <- do.call(c, unname(blocks[block_order]))

  list(idx = idx, init = init)
}

# Labels for the contrast blocks a categorical exposure carries, one per
# non-reference level in counterfactual order. The estimates table reports these
# in its contrast column and the contrast diagnostics name the contrast they
# concern, so both read them from here rather than rebuilding the pairs.
ipw_contrast_labels <- function(spec) {
  nonref <- names(spec$outcome$X_counterfactual)[-1]
  paste(nonref, "vs", spec$reference_level)
}

# The positions of the mu block in the order the estimates table reports the
# counterfactual means: level order, with the reference level first. A
# categorical spec seeds its means in that order already. A binary spec seeds
# the pair exposed-first, because the contrast rows built from it read the
# exposed mean as the first of the pair, so its two positions are swapped here
# rather than in the seed, which would move every contrast that reads them.
ipw_mu_order <- function(spec) {
  if (identical(spec$exposure_type, "binary")) {
    return(c(2L, 1L))
  }

  seq_along(spec$exposure_levels)
}

# The same order over the stratum mean block, which stores one mean per exposure
# level per stratum, stratum-major. Empty where `.by` names nothing, which is
# the block the layout leaves empty too.
ipw_by_mu_order <- function(spec) {
  if (is.null(spec$by)) {
    return(integer(0))
  }

  order <- ipw_mu_order(spec)
  offsets <- (seq_along(spec$by$labels) - 1L) * length(order)
  as.integer(unlist(lapply(offsets, function(offset) offset + order)))
}

# The transform a contrast form applies to a pair of marginal means.
ipw_contrast_transform <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo),
    "log(or)" = stats::qlogis(mu_hi) - stats::qlogis(mu_lo)
  )
}

# Whether a marginal mean lies where the form's transform is defined. A
# difference is defined everywhere, a log risk ratio needs a positive mean, and
# a log odds ratio needs a mean strictly inside (0, 1). A mean that is already
# NA carries no information about the domain and is left to the transform.
ipw_contrast_defined <- function(form, mu) {
  switch(
    form,
    rd = TRUE,
    diff = TRUE,
    "log(rr)" = is.na(mu) || mu > 0,
    "log(or)" = is.na(mu) || (mu > 0 && mu < 1)
  )
}

# The marginal means are free parameters of theta, so the solver moves them
# through values where the transform has nothing to return: a risk pushed past 1
# has no logit and a risk pushed below 0 has no logarithm. Base signals an
# unclassed "NaNs produced" there, naming neither the effect nor the contrast it
# belongs to. `reporter` replaces that with a classed warning that names both.
# The transform still runs on the same values and returns the same result, so
# the numbers the solver sees are unchanged.
ipw_contrast_value <- function(form, mu_hi, mu_lo, reporter = NULL) {
  if (ipw_contrast_defined(form, mu_hi) && ipw_contrast_defined(form, mu_lo)) {
    return(ipw_contrast_transform(form, mu_hi, mu_lo))
  }

  if (!is.null(reporter)) {
    reporter(form)
  }

  suppressWarnings(ipw_contrast_transform(form, mu_hi, mu_lo))
}

# A reporter for one block, labelled by the contrast, the stratum, or both, or
# for the single binary contrast when `label` is NULL. The solver revisits the
# same out-of-domain marginal means on every step and again on every column of
# the bread, so each effect reports once for the fit the reporter was built for.
#
# The init path builds no reporter. It seeds theta from the fitted models and
# the solver evaluates psi at that seed, so anything undefined there is reported
# from the psi block instead of twice.
ipw_contrast_reporter <- function(label = NULL) {
  reported <- new.env(parent = emptyenv())

  function(form) {
    if (isTRUE(reported[[form]])) {
      return(invisible(NULL))
    }
    reported[[form]] <- TRUE

    headline <- if (is.null(label)) {
      "The {.val {form}} effect is undefined at the marginal means the solver \\
      reached."
    } else {
      "The {.val {form}} effect for {.val {label}} is undefined at the \\
      marginal means the solver reached."
    }
    domain <- switch(
      form,
      "log(rr)" = "a positive marginal mean",
      "log(or)" = "a marginal mean strictly between 0 and 1"
    )

    warn(
      c(
        headline,
        i = "{.val {form}} needs {domain} on each side of the comparison, and \\
        at least one side is outside that range.",
        i = "An exposure level whose fitted outcomes are all events, or all \\
        non-events, drives its marginal mean to the boundary. Check the \\
        outcome within each level of the exposure.",
        i = "Estimates and standard errors from this fit are not reliable."
      ),
      warning_class = "propensity_ipw_contrast_warning",
      # The reporter runs inside the estimating function, which the solver
      # calls; the call reported here would name nothing the caller wrote.
      call = NULL
    )
  }
}

# Plug-in contrast init values, named by form (with an optional level suffix for
# categorical exposures).
ipw_init_contrasts <- function(contrasts, mu_hi, mu_lo, suffix = NULL) {
  if (is.null(contrasts) || !length(contrasts)) {
    return(numeric(0))
  }
  vals <- vapply(
    contrasts,
    function(f) ipw_contrast_value(f, mu_hi, mu_lo),
    numeric(1)
  )
  nm <- if (is.null(suffix)) contrasts else paste0(contrasts, "_", suffix)
  stats::setNames(vals, nm)
}

# The stabilizer to use for a binary exposure whose weights carry no
# per-observation stabilization score: the marginal exposure probability.
#
# Four consumers must agree on this value, and they read it from here rather
# than rebuilding it so that changing the convention changes all four at once.
# `ate_binary()` is the origin: it scales the default stabilized ATE weights by
# this probability and its complement, so the value is already baked into the
# weights a user hands to `ipw()`. The M-estimation init seeds the `stab_pi`
# parameter with it. The linearization preflight rebuilds a stabilized weight
# from it to compare against the weights that fit the outcome model, and
# `effective_stabilizer()` recovers the same factor to scale the linearization
# weight derivatives. A disagreement among them is silent in different ways: the
# preflight would reject weights the solver would have accepted, or accept
# weights it would not, while `effective_stabilizer()` would move only the
# standard errors. `exposure` is the 0/1 recode everywhere but `ate_binary()`,
# which passes the same recode under its own name.
ipw_default_stab_seed <- function(exposure) {
  mean(exposure, na.rm = TRUE)
}

# The probability the stabilization block reports at one value of theta: the
# single marginal proportion the default stabilizer estimates, or, where the
# caller stabilized on a fitted model, the probability that model gives each
# unit, rebuilt from the numerator design and the block's own parameters. The
# preflight that rebuilds the weights once and the psi that rebuilds them at
# every evaluation both read it here, so the two cannot compute a different
# numerator from the same parameters. Weights with no numerator to estimate
# leave the block empty and take no probability at all.
#
# `model` is the numerator block itself and `th_stab` the slice of theta that
# belongs to it, rather than the whole spec and the whole stabilization block:
# the joint route carries one such pair per component and would have no spec to
# hand over that named a single one of them.
ipw_binary_stab_prob <- function(model, th_stab, call = rlang::caller_env()) {
  if (!length(th_stab)) {
    return(NULL)
  }

  if (is.null(model)) {
    return(th_stab[[1]])
  }

  ipw_inv_link(model$link, call = call)(as.vector(model$X %*% th_stab))
}

ipw_init_binary <- function(spec, call = rlang::caller_env()) {
  ps_block <- spec$ps$coefs

  if (spec$stab$stabilized && is.null(spec$stab$score)) {
    # A numerator model's block is the shape the propensity score block is: one
    # parameter per coefficient, seeded at the coefficients the model was fit
    # at, which is the exact root of the score written for it. The default
    # stabilizer estimates the one marginal proportion instead.
    model <- spec$stab$model
    stab_block <- if (is.null(model)) {
      c(stab_pi = ipw_default_stab_seed(spec$exposure))
    } else {
      stats::setNames(model$coefs, paste0("stab_", colnames(model$X)))
    }
  } else {
    stab_block <- numeric(0)
  }

  beta <- spec$outcome$coefs
  inv_out <- ipw_inv_link(spec$outcome$link, call = call)
  pred1 <- inv_out(as.vector(spec$outcome$X_counterfactual$X1 %*% beta))
  pred0 <- inv_out(as.vector(spec$outcome$X_counterfactual$X0 %*% beta))
  # Seed the marginal means at the tilt-standardized weighted mean so the init
  # is the exact root of the standardized mu rows. ate keeps the ordinary mean
  # so its seed, and therefore its solve, is unchanged.
  h <- ipw_binary_seed_tilt(spec, call = call)
  if (is.null(h)) {
    mu1 <- mean(pred1)
    mu0 <- mean(pred0)
  } else {
    mu1 <- stats::weighted.mean(pred1, h)
    mu0 <- stats::weighted.mean(pred0, h)
  }
  mu_block <- c(mu1 = mu1, mu0 = mu0)

  con_block <- ipw_init_contrasts(spec$contrasts, mu1, mu0)

  by_block <- ipw_init_by(
    spec,
    preds = list(pred1, pred0),
    reference = 2L,
    comparisons = 1L,
    mu_names = names(mu_block),
    contrast_labels = NULL,
    tilt = h
  )

  list(
    ps = ps_block,
    stab = stab_block,
    out = beta,
    mu = mu_block,
    contrast = con_block,
    by_mu = by_block$mu,
    by_contrast = by_block$contrast
  )
}

# The tilt a binary estimand standardizes to at the seeded propensity score
# coefficients, or NULL for ate, whose tilt is the constant one. The
# whole-sample means and the stratum means seeded from them both standardize
# over it, so it is evaluated once and handed to both rather than rebuilt at
# each: a stratum seeded against a different tilt would not sit at the root of
# the row it seeds.
ipw_binary_seed_tilt <- function(spec, call = rlang::caller_env()) {
  if (identical(spec$estimand, "ate")) {
    return(NULL)
  }

  e_fit <- ipw_inv_link(spec$ps$link, call = call)(as.vector(
    spec$ps$X %*% spec$ps$coefs
  ))

  ps_tilt_binary(e_fit, spec$estimand)
}

# Names for a multinomial coefficient block, which is
# as.vector(t(coef(multinom))): level-major, term-minor, over the
# non-reference levels. `prefix` says which side of the weight the block is,
# the denominator writing its coefficients under their bare names and a
# numerator model's under the "stab_" every stabilization block carries. Both
# sides read the convention here rather than each writing it, because a name
# the two spelled differently would report one block under the other's terms.
ipw_name_categorical_coefs <- function(levels, terms, prefix = "") {
  as.vector(vapply(
    levels,
    function(l) paste0(prefix, l, ":", terms),
    character(length(terms))
  ))
}

# Generate readable names for the multinomial ps coefficient block.
ipw_name_categorical_ps <- function(spec) {
  levs <- names(spec$outcome$X_counterfactual)
  nm <- ipw_name_categorical_coefs(levs[-1], colnames(spec$ps$X))

  stats::setNames(spec$ps$coefs, nm)
}

ipw_init_categorical <- function(spec, call = rlang::caller_env()) {
  ps_block <- ipw_name_categorical_ps(spec)

  levs <- names(spec$outcome$X_counterfactual)

  if (spec$stab$stabilized && is.null(spec$stab$score)) {
    # A numerator model's block is the shape the propensity score block is: one
    # parameter per coefficient of the multinomial fit, seeded at the
    # coefficients it was fit at, which is the exact root of the score written
    # for it. The default stabilizer estimates the k - 1 free marginal
    # proportions instead.
    model <- spec$stab$model
    stab_block <- if (is.null(model)) {
      props <- colMeans(spec$exposure)
      stats::setNames(props[-1], paste0("stab_", levs[-1]))
    } else {
      stats::setNames(
        model$coefs,
        ipw_name_categorical_coefs(levs[-1], colnames(model$X), "stab_")
      )
    }
  } else {
    stab_block <- numeric(0)
  }

  beta <- spec$outcome$coefs
  inv_out <- ipw_inv_link(spec$outcome$link, call = call)
  k <- spec$ps$k
  # Seed each level's marginal mean at the tilt-standardized weighted mean, the
  # exact root of the standardized mu rows. ate keeps the ordinary mean so its
  # seed, and therefore its solve, is unchanged.
  h <- ipw_categorical_seed_tilt(spec, levs)
  weight_mu <- if (is.null(h)) {
    function(pred) mean(pred)
  } else {
    function(pred) stats::weighted.mean(pred, h)
  }
  preds <- lapply(
    seq_len(k),
    function(j) inv_out(as.vector(spec$outcome$X_counterfactual[[j]] %*% beta))
  )
  mu_vals <- vapply(preds, weight_mu, numeric(1))
  mu_block <- stats::setNames(mu_vals, paste0("mu_", levs))

  # A declared crossing replaces the vs-reference contrasts with the joint
  # surface's own, built over the same cell means.
  con_block <- if (!is.null(spec$joint)) {
    ipw_init_joint(spec$joint, mu_vals)
  } else {
    con_list <- lapply(seq_len(k)[-1], function(j) {
      ipw_init_contrasts(
        spec$contrasts,
        mu_vals[[j]],
        mu_vals[[1]],
        suffix = levs[[j]]
      )
    })
    if (length(con_list)) do.call(c, con_list) else numeric(0)
  }

  by_block <- ipw_init_by(
    spec,
    preds = preds,
    reference = 1L,
    comparisons = seq_len(k)[-1],
    mu_names = names(mu_block),
    contrast_labels = ipw_contrast_labels(spec),
    tilt = h
  )

  list(
    ps = ps_block,
    stab = stab_block,
    out = beta,
    mu = mu_block,
    contrast = con_block,
    by_mu = by_block$mu,
    by_contrast = by_block$contrast
  )
}

# The tilt a categorical estimand standardizes to at the seeded propensity score
# coefficients, or NULL for ate, whose tilt is the constant one. The
# whole-sample means and the stratum means seeded from them both standardize
# over it, the way `ipw_binary_seed_tilt()` serves the binary blocks.
ipw_categorical_seed_tilt <- function(spec, levs) {
  if (identical(spec$estimand, "ate")) {
    return(NULL)
  }

  ps_fit <- ipw_categorical_ps(spec$ps$X, spec$ps$coefs, spec$ps$k)
  focal_idx <- if (!is.null(spec$focal_level)) {
    match(spec$focal_level, levs)
  } else {
    NULL
  }

  ps_tilt_categorical(ps_fit, spec$estimand, focal_idx)
}

# The residuals the conditional density's spread was estimated from and the
# exposure a marginal numerator is the normal density of, both as the propensity
# score fit came to them over its own rows. `ipw_spec_continuous()` reads them
# off the model and carries them, because the rows `.data` leaves need not be
# the rows the fit was made over and the weights were built at the fit's
# reading. A spec assembled without that record is one whose rows are the fit's,
# so the rows it holds answer the same question.
ipw_continuous_fit_residuals <- function(spec) {
  if (!is.null(spec$ps$fit_residuals)) {
    return(spec$ps$fit_residuals)
  }

  spec$exposure - ipw_continuous_spec_fns(spec)$mean(spec$ps$X, spec$ps$coefs)
}

ipw_continuous_fit_exposure <- function(spec) {
  if (!is.null(spec$ps$fit_exposure)) {
    return(spec$ps$fit_exposure)
  }

  spec$exposure
}

# The same reading for a numerator model, whose block carries the spread its own
# fit's residuals came to.
ipw_numerator_fit_sigma2 <- function(spec, model) {
  if (!is.null(model$sigma2_fit)) {
    return(model$sigma2_fit)
  }

  fitted_n <- ipw_numerator_model_fns(model)$mean(model$X, model$coefs)

  mean((spec$exposure - fitted_n)^2)
}

ipw_init_continuous <- function(spec, call = rlang::caller_env()) {
  alpha <- spec$ps$coefs

  # A spread the caller fixed is a constant the weights were built with rather
  # than a quantity the data estimate, so the block is the coefficients alone
  # and nothing in the stacked system carries its uncertainty.
  #
  # The residuals the spread is read from are the propensity score fit's own,
  # carried on the spec, rather than the ones the analyzed rows leave. The
  # weights were built at the fit's moment, and a seed at any other moment
  # rebuilds weights the caller was never given.
  ps_block <- if (identical(spec$sigma$kind, "fixed")) {
    alpha
  } else {
    c(
      alpha,
      sigma2_d = ipw_continuous_sigma2_seed(
        spec$sigma,
        ipw_continuous_fit_residuals(spec),
        spec$density,
        call = call
      )
    )
  }

  # Only the marginal density of the exposure and a numerator model of the
  # caller's are read at parameters of their own. A marginalized conditional
  # density is built from the propensity score block and the data, and a fixed
  # score and unstabilized weights carry no numerator to estimate, so none of
  # the three leaves a stabilization block behind.
  #
  # A numerator model's block is the shape the propensity score block is: one
  # parameter per coefficient, seeded at the coefficients the model was fit at,
  # and one for the spread its density is read at, seeded at the second moment
  # of its residuals. Both seeds are the exact root of the row that estimates
  # them, which is what makes the stacked system carry the uncertainty of having
  # fit the model rather than move away from it.
  #
  # Both blocks are seeded at moments their own fit came to, taken over that
  # fit's rows, since those are the moments the weights were built at. The rows
  # the spec analyzes are what the equations below solve the moments over, which
  # is the answer rather than the starting value.
  if (identical(spec$numerator, "model")) {
    model <- spec$stab$model
    coefs <- model$coefs
    stab_block <- c(
      stats::setNames(coefs, paste0("stab_", colnames(model$X))),
      sigma2_n = ipw_numerator_fit_sigma2(spec, model)
    )
  } else if (identical(spec$numerator, "marginal")) {
    a_fit <- ipw_continuous_fit_exposure(spec)
    mu_a <- mean(a_fit)
    sigma2_a <- mean((a_fit - mu_a)^2)
    stab_block <- c(mu_a = mu_a, sigma2_a = sigma2_a)
  } else {
    stab_block <- numeric(0)
  }

  list(
    ps = ps_block,
    stab = stab_block,
    out = spec$outcome$coefs,
    mu = numeric(0),
    contrast = numeric(0),
    by_mu = numeric(0),
    by_contrast = numeric(0)
  )
}

# ---- psi builder ------------------------------------------------------------

build_ipw_psi <- function(spec, layout, call = rlang::caller_env()) {
  weight_fn <- ipw_weight_fn(
    spec$exposure_type,
    spec$estimand,
    spec$ps$types,
    call = call
  )
  switch(
    spec$exposure_type,
    binary = ipw_psi_binary(spec, layout, weight_fn, call = call),
    categorical = ipw_psi_categorical(spec, layout, weight_fn, call = call),
    continuous = ipw_psi_continuous(spec, layout, weight_fn, call = call),
    joint_models = ipw_psi_joint_models(spec, layout, weight_fn, call = call)
  )
}

# Deterministic contrast rows: each form's residual repeated across the n
# observation columns so the stacked matrix stays rectangular.
ipw_contrast_rows <- function(
  contrasts,
  mu_hi,
  mu_lo,
  th_con,
  n,
  reporter = NULL
) {
  if (is.null(contrasts) || !length(contrasts)) {
    return(NULL)
  }
  rows <- lapply(seq_along(contrasts), function(i) {
    val <- ipw_contrast_value(contrasts[[i]], mu_hi, mu_lo, reporter = reporter)
    matrix(val - th_con[[i]], nrow = 1, ncol = n)
  })
  do.call(rbind, rows)
}

ipw_stack <- function(blocks) {
  do.call(rbind, blocks[!vapply(blocks, is.null, logical(1))])
}

# Weighted outcome-score block. Branch on the outcome family so a gaussian
# outcome uses the weighted least-squares score and a binomial outcome uses the
# weighted GLM score with the outcome link.
ipw_outcome_rows <- function(theta, x, y, family, link, weights) {
  if (family == "binomial") {
    deli::ee_glm(
      theta,
      X = x,
      y = y,
      distribution = "binomial",
      link = link,
      weights = weights
    )
  } else {
    deli::ee_regression(
      theta,
      X = x,
      y = y,
      model = "linear",
      weights = weights
    )
  }
}

ipw_psi_binary <- function(
  spec,
  layout,
  weight_fn,
  call = rlang::caller_env()
) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  ps_link <- spec$ps$link
  inv_ps <- ipw_inv_link(ps_link, call = call)
  z <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  inv_out <- ipw_inv_link(out_link, call = call)
  x1 <- spec$outcome$X_counterfactual$X1
  x0 <- spec$outcome$X_counterfactual$X0
  score <- spec$stab$score
  contrasts <- spec$contrasts
  estimand <- spec$estimand
  # The ate tilt is the constant 1: no tilt is evaluated and none is multiplied
  # into the marginal-mean rows.
  tilted <- !identical(estimand, "ate")
  n <- spec$n
  # Built once for the fit, not once per evaluation, so the report survives the
  # solver's repeated visits to the same undefined means.
  reporter <- ipw_contrast_reporter()
  by <- spec$by
  # One reporter per stratum, on the same terms: a stratum's marginal means can
  # leave the domain of a transform while the whole sample's do not. The
  # effect-modification rows are differences of the stratum contrast parameters
  # and have no domain of their own, so they need none.
  by_reporters <- ipw_by_reporters(by, NULL)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    ps_rows <- deli::ee_glm(
      th_ps,
      X = x_ps,
      y = z,
      distribution = "binomial",
      link = ps_link
    )
    e <- inv_ps(as.vector(x_ps %*% th_ps))

    stab_prob <- ipw_binary_stab_prob(spec$stab$model, th_stab, call = call)
    # The tilt enters an evaluation twice, as the numerator of the weights and
    # as the standardization of the marginal-mean rows below. It is evaluated
    # here and handed to both.
    h <- if (tilted) ps_tilt_binary(e, estimand) else NULL
    w <- weight_fn(e, z, list(stab_prob = stab_prob, score = score, tilt = h))

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    # The numerator's own equations: the binomial score its coefficients solve,
    # written over the exposure the weights divide. The default stabilizer
    # estimates the marginal proportion instead, and a fixed score and
    # unstabilized weights estimate nothing.
    stab_rows <- if (!length(th_stab)) {
      NULL
    } else if (!is.null(spec$stab$model)) {
      deli::ee_glm(
        th_stab,
        X = spec$stab$model$X,
        y = z,
        distribution = "binomial",
        link = spec$stab$model$link
      )
    } else {
      matrix(z - th_stab[[1]], nrow = 1)
    }

    mu1 <- th_mu[[1]]
    mu0 <- th_mu[[2]]
    # Standardize the marginal means to the estimand's tilted population. The
    # tilt is recomputed from the ps block of theta on every evaluation, so the
    # sandwich variance accounts for propensity score estimation. The root of
    # sum_i h(e_i)(pred_a(x_i) - mu_a) = 0 is mu_a = weighted.mean(pred_a, h),
    # and h = 1 for ate leaves the unweighted marginal means.
    pred1 <- inv_out(as.vector(x1 %*% th_out))
    pred0 <- inv_out(as.vector(x0 %*% th_out))
    resid1 <- pred1 - mu1
    resid0 <- pred0 - mu0
    mu_rows <- if (tilted) {
      rbind(h * resid1, h * resid0)
    } else {
      rbind(resid1, resid0)
    }

    con_rows <- ipw_contrast_rows(
      contrasts,
      mu1,
      mu0,
      th_con,
      n,
      reporter = reporter
    )

    if (is.null(by)) {
      by_mu_rows <- NULL
      by_con_rows <- NULL
    } else {
      by_rows <- ipw_by_psi_rows(
        preds = list(pred1, pred0),
        reference = 2L,
        comparisons = 1L,
        tilt = h,
        indicators = by$indicators,
        forms = by$contrasts,
        th_mu = theta[idx$by_mu],
        th_con = theta[idx$by_contrast],
        reporters = by_reporters
      )
      by_mu_rows <- by_rows$mu
      by_con_rows <- by_rows$contrast
    }

    ipw_stack(list(
      ps_rows,
      stab_rows,
      out_rows,
      mu_rows,
      con_rows,
      by_mu_rows,
      by_con_rows
    ))
  }
}

# Reconstruct the n-by-K propensity score matrix from the multinomial ps block,
# matching deli::ee_mlogit's internal softmax. Column order is reference-first.
# A softmax is invariant to a common shift of its linear predictors, so each row
# is shifted by its largest predictor, the reference category's implicit 0
# included, before exponentiating. Without the shift a linear predictor above
# about 709 overflows exp() to Inf and the row normalizes to NaN, which the
# solver can reach on its own while iterating toward the root.
ipw_categorical_ps <- function(x, theta, k) {
  n <- nrow(x)
  b <- ncol(x)
  eta <- matrix(0, nrow = n, ncol = k)
  shift <- rep(0, n)
  for (j in seq_len(k - 1)) {
    idx <- ((j - 1) * b + 1):(j * b)
    eta[, j + 1] <- as.vector(x %*% theta[idx])
    shift <- pmax(shift, eta[, j + 1])
  }
  ps <- exp(eta - shift)
  ps / rowSums(ps)
}

# The probabilities the categorical stabilization block reports at one value of
# theta, in the column order the exposure indicator carries: the k - 1 free
# marginal proportions the default stabilizer estimates, completed by the
# reference level's, or, where the caller stabilized on a fitted model, the
# n-by-K matrix that model gives, rebuilt from the numerator design and the
# block's own parameters. The preflight that rebuilds the weights once and the
# psi that rebuilds them at every evaluation both read it here, so the two
# cannot compute a different numerator from the same parameters. Weights with no
# numerator to estimate leave the block empty and take no probability at all.
#
# `stab_block` is the stabilization block itself and `th_stab` the slice of
# theta that belongs to it, rather than the whole spec, on the terms
# `ipw_binary_stab_prob()` takes its own pair.
ipw_categorical_stab_probs <- function(stab_block, th_stab, k) {
  if (!length(th_stab)) {
    return(NULL)
  }

  model <- stab_block$model
  if (is.null(model)) {
    return(c(1 - sum(th_stab), th_stab))
  }

  # The same softmax the denominator is rebuilt through, so the numerator a
  # model reports is read back the way the model itself computes it, reference
  # level first.
  ipw_categorical_ps(model$X, th_stab, k)
}

ipw_psi_categorical <- function(
  spec,
  layout,
  weight_fn,
  call = rlang::caller_env()
) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  k <- spec$ps$k
  z_ind <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  inv_out <- ipw_inv_link(out_link, call = call)
  x_cf <- spec$outcome$X_counterfactual
  levs <- names(x_cf)
  stab <- spec$stab
  score <- stab$score
  focal_idx <- if (!is.null(spec$focal_level)) {
    match(spec$focal_level, levs)
  } else {
    NULL
  }
  contrasts <- spec$contrasts
  estimand <- spec$estimand
  tilted <- !identical(estimand, "ate")
  n <- spec$n
  # One reporter per contrast, built once for the fit rather than once per
  # evaluation, so each undefined effect reports once however often the solver
  # revisits it. The labels are the ones the estimates table reports.
  contrast_labels <- ipw_contrast_labels(spec)
  reporters <- lapply(contrast_labels, ipw_contrast_reporter)
  by <- spec$by
  by_reporters <- ipw_by_reporters(by, contrast_labels)
  joint <- spec$joint
  joint_reporters <- ipw_joint_reporters(joint)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    ps_rows <- deli::ee_mlogit(th_ps, X = x_ps, y = z_ind)
    ps_mat <- ipw_categorical_ps(x_ps, th_ps, k)

    stab_probs <- ipw_categorical_stab_probs(stab, th_stab, k)
    # Evaluated once and handed to both the weights and the marginal-mean rows.
    h <- if (tilted) {
      ps_tilt_categorical(ps_mat, estimand, focal_idx)
    } else {
      NULL
    }
    w <- weight_fn(
      ps_mat,
      z_ind,
      list(
        focal_idx = focal_idx,
        stab_probs = stab_probs,
        score = score,
        tilt = h
      )
    )

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    stab_rows <- if (!length(th_stab)) {
      NULL
    } else if (!is.null(stab$model)) {
      # The numerator model's own score, which is the multinomial score the
      # denominator's block solves, written over the numerator's design and the
      # same reference-first indicator matrix.
      deli::ee_mlogit(th_stab, X = stab$model$X, y = z_ind)
    } else {
      do.call(
        rbind,
        lapply(seq_len(k - 1), function(j) {
          matrix(z_ind[, j + 1] - th_stab[[j]], nrow = 1)
        })
      )
    }

    # Standardize each level's marginal mean to the estimand's tilted
    # population, recomputing the tilt from the ps block on every evaluation.
    # Row j of the block is h_i (pred_j(x_i) - mu_j) across the n columns, which
    # makes the root mu_j = weighted.mean(pred_j, h). Each row is built at its
    # own length n, so the tilt is never expanded to the full k-by-n matrix.
    preds <- lapply(
      seq_len(k),
      function(j) inv_out(as.vector(x_cf[[j]] %*% th_out))
    )
    mu_rows <- t(vapply(
      seq_len(k),
      function(j) {
        resid <- preds[[j]] - th_mu[[j]]
        if (tilted) h * resid else resid
      },
      numeric(n)
    ))

    con_rows <- if (!is.null(joint)) {
      ipw_joint_psi_rows(
        joint,
        th_mu,
        theta[idx$contrast],
        n,
        joint_reporters
      )
    } else if (!is.null(contrasts) && length(contrasts)) {
      mu_ref <- th_mu[[1]]
      con_index <- 0L
      row_list <- list()
      for (j in seq_len(k)[-1]) {
        mu_j <- th_mu[[j]]
        reporter <- reporters[[j - 1L]]
        for (form in contrasts) {
          con_index <- con_index + 1L
          val <- ipw_contrast_value(form, mu_j, mu_ref, reporter = reporter)
          row_list[[con_index]] <- matrix(
            val - th_con[[con_index]],
            nrow = 1,
            ncol = n
          )
        }
      }
      do.call(rbind, row_list)
    } else {
      NULL
    }

    if (is.null(by)) {
      by_mu_rows <- NULL
      by_con_rows <- NULL
    } else {
      by_rows <- ipw_by_psi_rows(
        preds = preds,
        reference = 1L,
        comparisons = seq_len(k)[-1],
        tilt = h,
        indicators = by$indicators,
        forms = by$contrasts,
        th_mu = theta[idx$by_mu],
        th_con = theta[idx$by_contrast],
        reporters = by_reporters
      )
      by_mu_rows <- by_rows$mu
      by_con_rows <- by_rows$contrast
    }

    ipw_stack(list(
      ps_rows,
      stab_rows,
      out_rows,
      mu_rows,
      con_rows,
      by_mu_rows,
      by_con_rows
    ))
  }
}

ipw_psi_continuous <- function(
  spec,
  layout,
  weight_fn,
  call = rlang::caller_env()
) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  a <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  ps_fns <- ipw_continuous_spec_fns(spec)

  # The grid an integrated numerator marginalizes over is fixed by the exposure
  # rather than by theta, so it is built once here and read at every evaluation.
  # Rebuilding it inside the closure would leave it unchanged and cost the
  # sandwich a pass over the data for each of its finite differences.
  grid <- ipw_continuous_grid(spec)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]

    inputs <- ipw_continuous_inputs(spec, th_ps, th_stab, grid = grid)

    ps_score <- ps_fns$score(inputs$alpha, x_ps, a)

    # A spread the caller fixed is a constant, and a constant has no equation.
    # The conditional variance is estimated only where the weights estimated it,
    # by the equation the spread they were built with is the root of.
    ps_rows <- if (identical(spec$sigma$kind, "fixed")) {
      ps_score
    } else {
      rbind(
        ps_score,
        ipw_continuous_sigma_row(
          spec$sigma,
          a - inputs$mu,
          inputs$extras$sigma2_d,
          spec$density
        )
      )
    }

    # The numerator's own equations: the score its coefficients solve, and the
    # second moment of its residuals, which is the spread its density is read
    # at. The marginal numerator estimates the exposure's own two moments
    # instead, and every other numerator estimates nothing.
    stab_rows <- if (!length(th_stab)) {
      NULL
    } else if (identical(spec$numerator, "model")) {
      model <- spec$stab$model
      p_n <- ncol(model$X)

      rbind(
        ipw_numerator_model_fns(model)$score(
          th_stab[seq_len(p_n)],
          model$X,
          a
        ),
        matrix((a - inputs$extras$mu_n)^2 - th_stab[[p_n + 1L]], nrow = 1)
      )
    } else {
      deli::ee_mean_variance(th_stab, y = a)
    }

    w <- weight_fn(inputs$mu, a, inputs$extras)

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    ipw_stack(list(ps_rows, stab_rows, out_rows))
  }
}

# What the continuous weight function reads at one value of theta: the
# conditional mean the propensity score block implies under its link, and the
# rest of the density ratio's inputs. The preflight that rebuilds the weights
# once and the psi that rebuilds them at every evaluation both come through
# here, so the two cannot compute a different weight from the same parameters.
ipw_continuous_inputs <- function(
  spec,
  th_ps,
  th_stab,
  grid = ipw_continuous_grid(spec)
) {
  n_alpha <- ncol(spec$ps$X)
  alpha <- th_ps[seq_len(n_alpha)]
  sigma2_d <- if (identical(spec$sigma$kind, "fixed")) {
    spec$sigma$value^2
  } else {
    th_ps[[n_alpha + 1L]]
  }

  # A numerator model reads its own conditional mean and its own spread out of
  # the stabilization block, where the marginal numerator reads the exposure's
  # two moments.
  model <- spec$stab$model
  numerator_model <- identical(spec$numerator, "model") && !is.null(model)
  mu_n <- NULL
  sigma_n <- NULL
  if (numerator_model) {
    p_n <- ncol(model$X)
    mu_n <- ipw_numerator_model_fns(model)$mean(model$X, th_stab[seq_len(p_n)])
    sigma_n <- sqrt(th_stab[[p_n + 1L]])
  }

  list(
    alpha = alpha,
    mu = ipw_continuous_spec_fns(spec)$mean(spec$ps$X, alpha),
    extras = list(
      sigma2_d = sigma2_d,
      mu_a = if (length(th_stab) && !numerator_model) th_stab[[1]],
      sigma2_a = if (length(th_stab) && !numerator_model) th_stab[[2]],
      mu_n = mu_n,
      sigma_n = sigma_n,
      score = spec$stab$score,
      stabilized = spec$stab$stabilized,
      density = spec$density,
      numerator = spec$numerator,
      grid = grid
    )
  )
}

# The points an integrated numerator averages the conditional density over,
# which `wt_ate()` builds from the exposure of the units it weights. The grid is
# a function of the data alone, so holding it fixed across the solve is what
# makes the weights the sandwich differentiates the weights the user was given.
ipw_continuous_grid <- function(spec) {
  ipw_numerator_grid(spec$exposure, spec$numerator)
}

# The same grid built from an exposure vector, which the joint route holds one
# of per component rather than one of outright.
ipw_numerator_grid <- function(exposure, numerator) {
  if (!identical(numerator, "integrated")) {
    return(NULL)
  }

  present <- exposure[!is.na(exposure)]

  seq(min(present), max(present), length.out = continuous_grid_n)
}

# The equation the stacked system estimates the conditional variance of the
# density by, which is the equation the spread the weights were built with is
# the root of. The pooled spread is the root mean square of the residuals, so
# its row is the uncentered second moment; a scale fit by maximum likelihood is
# the root of the score of the t itself for the scale, so its row is that score,
# multiplied through by the scale so that each residual enters through a bounded
# term.
ipw_continuous_sigma_row <- function(sigma, residuals, sigma2_d, density) {
  if (identical(sigma$kind, "mle")) {
    df <- density$params$df

    return(matrix(
      (df + 1) * residuals^2 / (df * sigma2_d + residuals^2) - 1,
      nrow = 1
    ))
  }

  matrix(residuals^2 - sigma2_d, nrow = 1)
}

# The seed for that parameter: the exact root of whichever row estimates it, so
# that the weights the system rebuilds at its starting value are the weights the
# user was given. That holds only where `residuals` are the ones the spread was
# estimated from, which are the fitted model's own over the rows it was fit
# over; the rows the system goes on to solve over are a separate question the
# row itself answers.
ipw_continuous_sigma2_seed <- function(
  sigma,
  residuals,
  density,
  call = rlang::caller_env()
) {
  if (identical(sigma$kind, "mle")) {
    return(t_sigma_mle(residuals, density$params$df, call = call)^2)
  }

  mean(residuals^2)
}
