# Resampled standard errors for a continuous exposure.
#
# The stacked routes describe a fit as a system of estimating equations and read
# the sandwich variance off it. Two continuous fits have no such description: an
# additive propensity score model chooses how much to smooth from the data, and
# a kernel density chooses its bandwidth from the residuals, so neither is a
# smooth function of a parameter vector. Resampling asks for neither. It repeats
# the whole fit on a resample of the rows and reads the spread of what comes
# back, so a step it cannot differentiate is simply a step it runs again.
#
# What is resampled is the pipeline the user wrote: the rows of `.data`, the
# propensity score model refit on them, the weights rebuilt from the record the
# supplied weights carry, and the marginal structural model refit with those
# weights. What is not resampled is the point estimate. That is the coefficient
# the supplied outcome model already reports, exactly as on the stacked route,
# so the two routes differ in the interval they put around a number and not in
# the number.

# The builder, which assembles a result from replicates rather than from a
# solved system. It runs the guards that belong to the estimator through
# `ipw_spec_continuous()` and lifts the two that belong to the sandwich.
ipw_continuous_bootstrap <- function(
  wt_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  boot_reps = 500L,
  boot_seed = NULL,
  effects = NULL,
  call = rlang::caller_env()
) {
  check_ipw_boot_data(.data, call = call)
  boot_reps <- check_boot_reps(boot_reps, call = call)
  check_boot_seed(boot_seed, call = call)
  check_ipw_ps_link_absent(ps_link, "continuous", call = call)

  # Guard the weights that fit the outcome model before building the spec, so a
  # modified propensity score is detected before any estimand parsing.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts, call = call)

  spec <- ipw_spec_continuous(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    stacked = FALSE,
    call = call
  )

  # The readings the result declares and the one it records, settled before any
  # replicate runs so that a request for a reading the fit does not have is
  # refused rather than answered after five hundred refits.
  readings <- c("marginal", "conditional")
  if (ipw_is_dose_basis(spec)) {
    check_ipw_dose_basis_effects(effects, call = call)
    if (is.null(effects)) {
      ipw_dose_basis_announce()
    }
    effects <- "conditional"
    readings <- "conditional"
  } else if (is.null(effects)) {
    effects <- "marginal"
  }

  reported <- colnames(spec$outcome$X)[spec$coefficients$col]
  context <- ipw_boot_context(wt_mod, outcome_mod, .data, spec)

  replicates <- ipw_boot_replicates(
    context,
    n_reps = boot_reps,
    boot_seed = boot_seed,
    reported = reported
  )

  ipw_boot_report(replicates, boot_reps, call = call)

  # The reported rows of the replicates, and the table their spread describes.
  # The table is built first because the labels its rows carry are the labels
  # the replicates are keyed by, and those labels are read off the table.
  reps <- replicates$coefs[, reported, drop = FALSE]
  estimates <- ipw_boot_estimates(spec, reps, conf_level = conf_level)
  colnames(reps) <- ipw_effect_labels(estimates)

  fit <- ipw_boot_fit(
    reps,
    n_reps = boot_reps,
    n_failed = replicates$n_failed,
    boot_seed = boot_seed,
    conf_level = conf_level
  )
  attr(estimates, "ipw_vcov") <- fit$vcov

  # The propensity score model carries no block of anything, because resampling
  # solves no system the two models share, so it is stored as it was supplied.
  # The outcome model carries the covariance of its own coefficients across the
  # replicates, which is what the conditional reading reports.
  new_ipw(
    estimand = spec$estimand,
    wt_mod = wt_mod,
    outcome_mod = new_ipw_model(outcome_mod, ipw_boot_vcov(replicates$coefs)),
    estimates = estimates,
    se_method = "bootstrap",
    fit = fit,
    effects = effects,
    readings = readings
  )
}

# Everything one replicate needs, read off the supplied fit once rather than
# once per replicate. The weights are rebuilt by calling `wt_ate()` the way the
# user called it, so what travels here is the record they left on the weights,
# translated back into the arguments that produced them: `"none"` is an
# unstabilized ratio, `"score"` is one the caller's own score stabilized, and
# the other two name the numerator `wt_ate()` builds.
ipw_boot_context <- function(wt_mod, outcome_mod, .data, spec) {
  numerator <- spec$numerator

  list(
    data = .data,
    ps_mod = wt_mod,
    outcome_mod = outcome_mod,
    exposure_name = fmla_extract_left_chr(wt_mod),
    density = spec$density,
    numerator = if (identical(numerator, "integrated")) {
      "integrated"
    } else {
      "marginal"
    },
    stabilize = !identical(numerator, "none"),
    score = if (identical(numerator, "score")) spec$stab$score,
    sigma = spec$sigma$value
  )
}

# One replicate: refit both models on the resampled rows and report the
# coefficients the outcome model comes back with. The caller draws the row
# indices, so this is a pure function of them and can be stood in for.
#
# The propensity score model is refit through its own call, evaluated in the
# environment its formula carries. That is where an unqualified `rlm()` or
# `gam()` was written, and where a formula built in a local variable still
# lives. The weights are then rebuilt from the record rather than reused: a
# kernel bandwidth is re-estimated on the replicate's residuals and an additive
# model re-selects its smoothing parameters, which is the point of resampling
# rather than stacking.
ipw_boot_once <- function(idx, context) {
  resampled <- context$data[idx, , drop = FALSE]

  ps_fit <- ipw_boot_refit(context$ps_mod, resampled)
  mu <- as.double(extract_continuous_ps(ps_fit)$mu)

  wts <- wt_ate(
    mu,
    as.double(resampled[[context$exposure_name]]),
    .sigma = context$sigma,
    exposure_type = "continuous",
    stabilize = context$stabilize,
    stabilization_score = ipw_boot_score(context$score, idx),
    .density = context$density,
    numerator = context$numerator
  )

  stats::coef(ipw_boot_refit(context$outcome_mod, resampled, wts))
}

# A stabilizing score the caller supplied travels with the rows it describes, so
# one value per observation is resampled alongside them. A single value
# describes every observation and is replayed as it stands.
ipw_boot_score <- function(score, idx) {
  if (is.null(score) || length(score) == 1L) {
    return(score)
  }

  score[idx]
}

# Refit a model on the resampled rows, and on the resampled weights where there
# are any. The call is the one the model was fit by, with its data and its
# weights redirected at bindings this function makes, so that everything else
# the user wrote, a family, a control list, a tuning constant, is carried over
# untouched. The bindings live in a child of the formula's environment, which
# leaves every name the call reads resolving where it resolved when the model
# was fit: an unqualified `rlm()` or `gam()`, a formula held in a local
# variable, a covariate the formula transforms.
#
# The call is evaluated there and the formula is put back into it carrying that
# environment as well. Both are needed. A fitting function evaluates its own
# call in the caller's frame, which is what finds the data, but builds the model
# frame in the formula's environment, which is what finds the weights, and the
# two are the same frame only if the formula is told so.
ipw_boot_refit <- function(model, data, weights = NULL) {
  env <- new.env(parent = environment(stats::formula(model)))
  env$boot_data <- data

  fml <- stats::formula(model)
  environment(fml) <- env

  model_call <- stats::update(model, evaluate = FALSE)
  model_call$formula <- fml
  model_call$data <- quote(boot_data)

  if (!is.null(weights)) {
    model_call$weights <- quote(boot_weights)
    env$boot_weights <- weights
  }

  eval(model_call, env)
}

# Run the replicates and report what survived. The rows every replicate
# resamples are drawn up front, so the draw is a function of the seed alone and
# a refit that consumes random numbers of its own cannot move it. The refits
# announce nothing: a rebuild of the weights is the same rebuild five hundred
# times over, and its alerts describe the fit the user already has.
ipw_boot_replicates <- function(context, n_reps, boot_seed, reported) {
  draws <- ipw_boot_draw(nrow(context$data), n_reps, boot_seed)

  quiet <- options(propensity.quiet = TRUE)
  on.exit(options(quiet), add = TRUE)

  results <- lapply(draws, function(idx) {
    tryCatch(ipw_boot_once(idx, context), error = function(e) NULL)
  })

  coefs <- ipw_boot_bind(results, reported)

  list(coefs = coefs, n_failed = n_reps - ipw_boot_kept(coefs))
}

# The row indices each replicate resamples. A seed of its own is set for the
# draw and the session's random state is put back afterwards, so a bootstrap
# leaves the stream where it found it and a caller who set no seed of their own
# gets a result their own seed reproduces.
ipw_boot_draw <- function(n, n_reps, boot_seed) {
  if (!is.null(boot_seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", state, envir = globalenv()), add = TRUE)
    } else {
      on.exit(
        suppressWarnings(rm(".Random.seed", envir = globalenv())),
        add = TRUE
      )
    }

    set.seed(boot_seed)
  }

  lapply(seq_len(n_reps), function(i) sample.int(n, n, replace = TRUE))
}

# The replicates that count, as a matrix of one row each. A replicate counts
# when it returned a coefficient vector that is finite throughout, holds every
# coefficient the table reports, and is the same vector the replicates before it
# returned. The last of those is what a rank-deficient refit fails: it drops a
# column the others carry, and a column the replicates disagree about has no
# covariance to report.
ipw_boot_bind <- function(results, reported) {
  kept <- list()
  reference <- NULL

  for (value in results) {
    if (!is.numeric(value) || !all(is.finite(value))) {
      next
    }

    if (!all(reported %in% names(value))) {
      next
    }

    if (is.null(reference)) {
      reference <- names(value)
    } else if (!identical(names(value), reference)) {
      next
    }

    kept[[length(kept) + 1L]] <- unname(value)
  }

  if (length(kept) == 0L) {
    return(NULL)
  }

  matrix(
    unlist(kept),
    nrow = length(kept),
    byrow = TRUE,
    dimnames = list(NULL, reference)
  )
}

ipw_boot_kept <- function(coefs) {
  if (is.null(coefs)) {
    return(0L)
  }

  nrow(coefs)
}

# What the count of dropped replicates decides. Too few that survived and there
# is no standard error to report, however well the rest went; a few and the
# result says how much of the request it was built from. The floor is checked
# first, because a request that cannot clear it is refused rather than warned
# about and then refused.
ipw_boot_floor <- 50L

ipw_boot_report <- function(replicates, n_reps, call = rlang::caller_env()) {
  kept <- ipw_boot_kept(replicates$coefs)

  if (kept < ipw_boot_floor) {
    abort(
      c(
        "{.fun ipw} needs at least {ipw_boot_floor} successful bootstrap \\
        replicates.",
        x = "{kept} of {n_reps} replicate{?s} succeeded, and \\
        {replicates$n_failed} failed.",
        i = "Raise {.arg boot_reps}, or refit the models on data every \\
        resample of them can be fit to."
      ),
      error_class = "propensity_ipw_bootstrap_error",
      call = call
    )
  }

  if (replicates$n_failed / n_reps > 0.05) {
    warn(
      c(
        "{replicates$n_failed} of {n_reps} bootstrap replicate{?s} failed and \\
        {?was/were} dropped.",
        x = "A replicate fails when a refit errors or leaves a coefficient \\
        that is not finite, so the spread is the spread of the resamples the \\
        models could be fit to.",
        i = "Check whether either model is fit to something a resample of the \\
        rows can break, such as a factor covariate with a level so rare that a \\
        resample can leave it out."
      ),
      warning_class = "propensity_ipw_bootstrap_warning",
      call = call
    )
  }

  invisible(TRUE)
}

# The fit a resampled result holds: the replicates of the reported coefficients,
# the seed they were drawn under, how much of the request they are, their
# covariance, and the percentile bounds of each of them. The percentile bounds
# are reported alongside the normal-approximation interval rather than in place
# of it, so that a reader who wants the interval the replicates draw themselves
# has it without the accessors reporting two different intervals.
ipw_boot_fit <- function(reps, n_reps, n_failed, boot_seed, conf_level) {
  alpha <- (1 - conf_level) / 2
  percentile <- t(apply(
    reps,
    2,
    stats::quantile,
    probs = c(alpha, 1 - alpha)
  ))

  structure(
    list(
      reps = reps,
      seed = boot_seed,
      n_reps = n_reps,
      n_failed = n_failed,
      vcov = ipw_boot_vcov(reps),
      percentile = percentile
    ),
    class = "ipw_bootstrap"
  )
}

# The covariance of a set of replicates, which is the whole of what resampling
# estimates. `stats::cov()` of a one-column matrix is the one-by-one matrix a
# single reported effect needs, so no case is special.
ipw_boot_vcov <- function(reps) {
  stats::cov(reps)
}

# The estimates table, built from the coefficients the supplied outcome model
# reports and the spread of the replicates around them. The point estimate is
# not resampled, so what the replicates contribute is the standard error alone.
ipw_boot_estimates <- function(spec, reps, conf_level) {
  ipw_estimates_frame(
    effect = spec$coefficients$effect,
    contrast = spec$coefficients$contrast,
    group = spec$coefficients$group,
    estimate = unname(spec$outcome$coefs[spec$coefficients$col]),
    std.err = unname(apply(reps, 2, stats::sd)),
    conf_level = conf_level
  )
}

# ---- what the bootstrap needs -----------------------------------------------

# Resampling refits both models on rows it selects, and the rows have to come
# from somewhere. A fitted model frame is not enough: the frame records the
# terms rather than the columns behind them, and the propensity score model is
# refit from its own call, which reads a data frame.
check_ipw_boot_data <- function(.data, call = rlang::caller_env()) {
  if (!is.null(.data)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} requires {.arg .data} for {.val bootstrap} standard errors.",
      x = "Each replicate refits both models on a resample of the rows, and \\
      there are no rows to resample without {.arg .data}.",
      i = "Supply the data the models were fit to in {.arg .data}."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# The two arguments that describe the resampling, refused under a standard error
# method that resamples nothing. They are not ignored there: a caller who set a
# replicate count under the sandwich asked for something, and the answer is that
# this method has no replicates rather than five hundred silent ones.
check_ipw_boot_args <- function(
  se_method,
  boot_reps_named,
  boot_seed_named,
  call = rlang::caller_env()
) {
  if (identical(se_method, "bootstrap")) {
    return(invisible(TRUE))
  }

  supplied <- c("boot_reps" = boot_reps_named, "boot_seed" = boot_seed_named)

  if (!any(supplied)) {
    return(invisible(TRUE))
  }

  named <- names(supplied)[supplied]

  abort(
    c(
      "{.arg {named}} {?is/are} not supported with {.val {se_method}} \\
      standard errors.",
      x = "{.arg {named}} describe{?s/} the resampling, and \\
      {.val {se_method}} resamples nothing.",
      i = "Use {.code se_method = \"bootstrap\"}, or drop {.arg {named}}."
    ),
    error_class = "propensity_unsupported_arg_error",
    call = call
  )
}

# How many replicates to run, which is a count and so is read as one. The floor
# on successful replicates is checked after they run rather than here: a request
# below it is refused for what it produced, since a request above it can fail
# the same way.
check_boot_reps <- function(boot_reps, call = rlang::caller_env()) {
  whole <- rlang::is_scalar_integerish(boot_reps) && !is.na(boot_reps)

  if (!whole || boot_reps < 1) {
    abort(
      c(
        "{.arg boot_reps} must be a single positive whole number.",
        x = "{.arg boot_reps} is {.val {boot_reps}}.",
        i = "Use {.code boot_reps = 500} to run five hundred replicates."
      ),
      error_class = "propensity_ipw_bootstrap_error",
      call = call
    )
  }

  as.integer(boot_reps)
}

# The seed the replicates are drawn under, or `NULL` for the session's own
# random state. A vector of numbers is not a seed here: `set.seed()` takes one.
check_boot_seed <- function(boot_seed, call = rlang::caller_env()) {
  if (is.null(boot_seed)) {
    return(invisible(boot_seed))
  }

  single <- (rlang::is_scalar_double(boot_seed) ||
    rlang::is_scalar_integer(boot_seed)) &&
    !is.na(boot_seed)

  if (!single) {
    abort(
      c(
        "{.arg boot_seed} must be a single number or {.code NULL}.",
        x = "{.arg boot_seed} is {.obj_type_friendly {boot_seed}}.",
        i = "Omit {.arg boot_seed} to draw from the session's random state."
      ),
      error_class = "propensity_ipw_bootstrap_error",
      call = call
    )
  }

  invisible(boot_seed)
}

# Resampling is written for the continuous route alone. A binary or categorical
# exposure has a sandwich for every fit it accepts, so nothing there needs a
# method that repeats the fit, and routing one into a loop written for a density
# ratio would rebuild weights it never had.
abort_ipw_boot_exposure <- function(
  exposure_type,
  call = rlang::caller_env()
) {
  abort(
    c(
      "{.fun ipw} supports {.val bootstrap} standard errors only for a \\
      continuous exposure.",
      x = "{.arg wt_mod} is a propensity score model of a {exposure_type} \\
      exposure.",
      i = "Use {.code se_method = \"mestimation\"}, which builds a sandwich \\
      variance for every fit this exposure type accepts."
    ),
    error_class = "propensity_method_error",
    call = call
  )
}

# The two-model joint route is refused for the same reason, said in its own
# terms. Its weights are the product of two treatment models' weights, and the
# loop above rebuilds one density ratio from the record one set of weights
# carries, so there is no product for it to repeat.
abort_ipw_boot_joint <- function(call = rlang::caller_env()) {
  abort(
    c(
      "{.fun ipw} supports {.val bootstrap} standard errors only for a \\
      continuous exposure weighted by one propensity score model.",
      x = "{.arg wt_mod} is a pair of treatment models, and each replicate \\
      rebuilds the weights of one model rather than the product of two.",
      i = "Use {.code se_method = \"mestimation\"}, which builds a sandwich \\
      variance for every fit a joint treatment model accepts."
    ),
    error_class = "propensity_method_error",
    call = call
  )
}
