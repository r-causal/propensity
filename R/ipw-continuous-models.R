# The propensity score models a continuous exposure's sandwich can be built
# from, and the refusals for the ones it cannot.
#
# Weights for a continuous exposure are a ratio of densities centered on the
# conditional mean the propensity score model fits, so the stacked system has to
# carry that model's own score: the equation its coefficients solve. A class
# whose coefficients solve a different equation would drift to the root of the
# one stacked here, and the estimates would describe a model the user never fit.
# The registry is therefore a list of the classes whose score can be written
# down, and every class is named rather than reached by inheritance.
#
# What the registry returns is data rather than a decision. `stackable` says
# whether the score exists here, `reason` says why not when it does not, and the
# caller that needs a stacked system asks for the refusal by calling
# `check_ipw_continuous_model()`. A class the registry has never heard of, and a
# link no propensity score model here is written for, are refused on the spot
# instead: neither describes a model this path could stack under any standard
# error method.
#
# Lookup is in class order rather than by inheritance, so a `gam` is read as a
# `gam` rather than as the `glm` it inherits from, and an `rlm` as an `rlm`
# rather than as the `lm` it inherits from.
ipw_continuous_model <- function(ps_mod, call = rlang::caller_env()) {
  ps_class <- class(ps_mod)

  if (inherits(ps_mod, "gam")) {
    return(ipw_continuous_entry(
      kind = "gam",
      classes = ps_class,
      link = ps_mod$family$link,
      stackable = FALSE,
      error_class = c(
        "propensity_ipw_se_method_unavailable_error",
        "propensity_method_error"
      ),
      reason = c(
        "{.fun ipw} cannot build a sandwich variance for a \\
        {.cls {entry$classes}} propensity score model of a continuous \\
        exposure.",
        x = "An additive model chooses how much to smooth by REML, and no \\
        estimating equation stacked here reproduces that choice, so the \\
        stacked system would describe a different fit.",
        i = ipw_continuous_bootstrap_hint()
      )
    ))
  }

  if (inherits(ps_mod, "rlm")) {
    return(ipw_continuous_rlm_entry(ps_mod, ps_class))
  }

  if (inherits(ps_mod, "glm")) {
    # The stacked score is the one deli writes for a normal model, so a family
    # whose spread changes with its fitted values is refused before its link is
    # read: there is no single conditional density for the weights to divide by
    # and no equation here that its coefficients solve.
    check_ipw_continuous_ps_family(ps_mod, call = call)

    link <- ps_mod$family$link
    check_ipw_continuous_ps_link(link, call = call)

    return(ipw_continuous_entry(kind = "glm", link = link, stackable = TRUE))
  }

  if (identical(ps_class, "lm")) {
    return(ipw_continuous_entry(
      kind = "lm",
      link = "identity",
      stackable = TRUE
    ))
  }

  abort(
    c(
      "{.fun ipw} supports only {.fun stats::lm}, gaussian \\
      {.fun stats::glm}, or {.fun MASS::rlm} propensity score models for a \\
      continuous exposure.",
      x = "{.arg wt_mod} has class {.cls {ps_class}}.",
      i = "A {.cls gam} is recognized and refused on its own terms; every \\
      other class reaches this refusal.",
      i = "Refit {.arg wt_mod} with {.fun stats::lm} or \\
      {.code stats::glm(family = gaussian())}."
    ),
    error_class = "propensity_class_error",
    call = call
  )
}

# The entry a robust fit gets, which is the one place a class carries a constant
# of its own into the stacked system. `MASS::rlm()` minimizes a loss rather than
# the sum of squares, so its coefficients are the root of the psi score read at
# the scale the fit settled on, and that score is stacked here as deli's Huber
# robust regression.
#
# The two clip the same residual on different scales: `rlm` clips the residual
# divided by its scale estimate at the psi's own constant, and deli clips the
# raw residual at `k`, so the equations agree when `k` is the product of the
# two. The constant is read out of the psi's formals rather than out of `k2`,
# because that is where `rlm` writes a caller's `k`; `k2` tunes the scale
# estimator instead and is unchanged when the psi's threshold moves.
#
# The scale itself enters as a known constant, and the system carries none of
# its uncertainty. That is a choice rather than an oversight: the MAD scale
# solves an equation of its own that this stack does not write, and a fit whose
# scale is uncertain enough for that to matter is better served by resampling.
ipw_continuous_rlm_entry <- function(ps_mod, classes = class(ps_mod)) {
  psi_name <- ipw_rlm_psi_name(ps_mod)
  huber <- identical(psi_name, "MASS::psi.huber")
  mm <- identical(ipw_rlm_method(ps_mod), "MM")

  if (mm || !huber) {
    found <- if (mm && !is.null(psi_name)) {
      "{.arg wt_mod} was fit with {.code method = \"MM\"}, which starts from a
       high-breakdown fit and finishes on {.fun {entry$psi}}."
    } else if (mm) {
      "{.arg wt_mod} was fit with {.code method = \"MM\"}, which starts from a
       high-breakdown fit and finishes on a redescending psi."
    } else if (is.null(psi_name)) {
      "{.arg wt_mod} was fit with a psi function this path cannot recognize."
    } else {
      "{.arg wt_mod} was fit with {.fun {entry$psi}}."
    }

    return(ipw_continuous_entry(
      kind = "rlm",
      classes = classes,
      link = "identity",
      stackable = FALSE,
      psi = psi_name,
      error_class = "propensity_ipw_robust_psi_error",
      reason = c(
        "{.fun ipw} stacks only the Huber score of a {.cls {entry$classes}} \\
        propensity score model of a continuous exposure.",
        x = found,
        i = "Refit {.arg wt_mod} with {.fun MASS::psi.huber}, the default, \\
        whose threshold {.fun ipw} reads off the fit.",
        i = ipw_continuous_bootstrap_hint()
      )
    ))
  }

  # A fit that stopped short of its own tolerance is not at the root of the
  # score stacked here, so a system seeded from its coefficients would move away
  # from them and report a propensity score model the user never saw.
  if (!isTRUE(ps_mod$converged)) {
    return(ipw_continuous_entry(
      kind = "rlm",
      classes = classes,
      link = "identity",
      stackable = FALSE,
      psi = psi_name,
      error_class = "propensity_ipw_convergence_error",
      reason = c(
        "{.fun ipw} cannot stack a {.cls {entry$classes}} propensity score \\
        model that did not converge.",
        x = "{.arg wt_mod} reports {.code converged = FALSE}, so its \\
        coefficients are not the root of the score stacked here.",
        i = "Refit {.arg wt_mod} with a larger {.arg maxit}, or a looser \\
        {.arg acc}, until it converges."
      )
    ))
  }

  ipw_continuous_entry(
    kind = "rlm",
    classes = classes,
    link = "identity",
    stackable = TRUE,
    psi = psi_name,
    huber_k = as.numeric(formals(ps_mod$psi)$k) * ps_mod$s
  )
}

# The method a robust fit was made by, which `rlm` records only in its call. A
# caller who named it through a variable left an unevaluated symbol there, so
# the call is evaluated in the environment the model's formula carries, which is
# where that variable was written. A value that cannot be evaluated, or is not a
# string, names no method here and reads as `NULL` rather than stopping the
# lookup.
ipw_rlm_method <- function(ps_mod) {
  method <- ps_mod$call$method
  if (is.null(method)) {
    return(NULL)
  }

  value <- tryCatch(
    eval(method, environment(stats::formula(ps_mod))),
    error = function(e) NULL
  )

  if (!is.character(value) || length(value) != 1L) {
    return(NULL)
  }

  value
}

# Which of MASS's psi functions a robust fit was given, or `NULL` for one this
# path has no name for. `rlm` tunes a psi by rewriting its formals, so a fit
# whose caller passed `k` no longer carries the function it was given and
# `identical()` on the function itself reports the wrong answer. The body is
# what the rewriting leaves alone, so that is what is compared.
ipw_rlm_psi_name <- function(ps_mod) {
  psi <- ps_mod$psi
  if (!is.function(psi)) {
    return(NULL)
  }

  known <- c(
    "MASS::psi.huber" = MASS::psi.huber,
    "MASS::psi.bisquare" = MASS::psi.bisquare,
    "MASS::psi.hampel" = MASS::psi.hampel
  )

  match <- vapply(
    known,
    function(fn) identical(body(psi), body(fn)),
    logical(1)
  )

  if (!any(match)) {
    return(NULL)
  }

  names(known)[[which(match)[[1]]]]
}

# One entry of the registry. `mean` reconstructs the conditional mean from the
# propensity score block of theta, and `score` writes the equation the model's
# coefficients solve; both are the identity's least squares unless the entry
# carries another link or another kind. A model that cannot be stacked carries
# neither, so a caller reaching for one has already had to pass the refusal.
#
# `psi` and `huber_k` are the two things a robust fit carries that no other
# entry has: the name of the psi it was fit with, which its refusal reads, and
# the threshold its score clips the raw residual at, which the spec carries so
# that the psi builder writes the same equation the registry did.
ipw_continuous_entry <- function(
  kind,
  link,
  stackable,
  classes = kind,
  error_class = NULL,
  reason = NULL,
  psi = NULL,
  huber_k = NULL
) {
  entry <- list(
    kind = kind,
    classes = classes,
    link = link,
    stackable = stackable,
    error_class = error_class,
    reason = reason,
    psi = psi,
    huber_k = huber_k,
    mean = NULL,
    score = NULL
  )

  if (!stackable) {
    return(entry)
  }

  utils::modifyList(entry, ipw_continuous_score_fns(kind, link, huber_k))
}

# The mean and the score a stacked propensity score model contributes, keyed on
# the registry entry's kind as well as its link. Every least squares fit is read
# through its link alone, and a robust fit is the one kind whose score is a
# different equation at the same identity link, so the kind is what separates
# them.
ipw_continuous_score_fns <- function(kind, link, huber_k = NULL) {
  if (!identical(kind, "rlm")) {
    return(ipw_continuous_link_fns(link))
  }

  list(
    mean = function(X, alpha) as.vector(X %*% alpha),
    score = function(alpha, X, y) {
      deli::ee_robust_regression(
        alpha,
        X = X,
        y = y,
        model = "linear",
        k = huber_k,
        loss = "huber"
      )
    }
  )
}

# The same pair read off a spec, which is what the psi builder and the preflight
# that has to agree with it both hold. The spec records the kind, the link, and
# the robust threshold rather than the model they came from, so the two rebuild
# the propensity score block from the same three values.
ipw_continuous_spec_fns <- function(spec) {
  ipw_continuous_score_fns(spec$ps$kind, spec$ps$link, spec$ps$huber_k)
}

# The two things a link decides: how the conditional mean is read back from the
# propensity score coefficients, and which equation those coefficients solve. An
# identity link is ordinary least squares, and every other link is the score
# deli writes for a normal model fit through it.
#
# A spec records the link and the kind rather than the model they came from, so
# the psi builder and the preflight that has to agree with it read these from
# here as well, which is what makes the two the same arithmetic.
ipw_continuous_link_fns <- function(link) {
  inv_link <- ipw_inv_link(link)

  score <- if (identical(link, "identity")) {
    function(alpha, X, y) {
      deli::ee_regression(alpha, X = X, y = y, model = "linear")
    }
  } else {
    function(alpha, X, y) {
      deli::ee_glm(alpha, X = X, y = y, distribution = "normal", link = link)
    }
  }

  list(
    mean = function(X, alpha) inv_link(as.vector(X %*% alpha)),
    score = score
  )
}

# The families whose score the stacked system writes. A gaussian model is the
# one it writes, and the weights are a ratio of densities with a single spread,
# which is the same model read from the other side.
check_ipw_continuous_ps_family <- function(
  ps_mod,
  call = rlang::caller_env()
) {
  family <- ps_mod$family
  if (identical(family$family, "gaussian")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} supports only a gaussian propensity score model for a \\
      continuous exposure.",
      x = "{.arg wt_mod} was fit with \\
      {.code {model_family_label(family)}}.",
      i = "Refit {.arg wt_mod} with {.fun stats::lm} or \\
      {.code stats::glm(family = gaussian())}."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# The links a gaussian propensity score model can be stacked on. An identity
# link is least squares and a log link is the score deli writes for a normal
# model of `exp(X alpha)`. The remaining gaussian links are refused by name
# rather than stacked: deli evaluates the score at the coefficients the IRLS
# iteration stopped at, and for those links that point is not close enough to
# the root to seed the solve from.
ipw_continuous_ps_links <- c("identity", "log")

check_ipw_continuous_ps_link <- function(link, call = rlang::caller_env()) {
  if (link %in% ipw_continuous_ps_links) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support the {.val {link}} link for the propensity \\
      score model of a continuous exposure.",
      x = "{.arg wt_mod} is a gaussian model fit with the {.val {link}} \\
      link.",
      i = "Refit {.arg wt_mod} with one of {.val {ipw_continuous_ps_links}}, \\
      or as an {.fun lm}."
    ),
    error_class = "propensity_ipw_link_error",
    call = call
  )
}

# Refuse a registry entry whose score the stacked system cannot write. The
# refusal carries the entry's own class and reason, so a model that is wrong for
# this path and a model that is right but has no sandwich are distinguishable by
# the condition they raise.
#
# A reason is glued where it is raised rather than where it is written, so one
# that names the model's classes reads them out of `entry` here.
check_ipw_continuous_model <- function(entry, call = rlang::caller_env()) {
  if (entry$stackable) {
    return(invisible(TRUE))
  }

  abort(entry$reason, error_class = entry$error_class, call = call)
}

# What ratio of densities a set of continuous weights is: the family, the
# numerator that stabilized it, and the spread the conditional density was read
# at. The weights record all three, and the record is what the stacked system
# rebuilds them from.
#
# Weights that carry no record were either written by hand or built before the
# record existed. Both are the ratio the package has always built: a normal
# density spread by the pooled residual standard deviation, over whatever their
# stabilization says stabilized them.
ipw_continuous_ratio <- function(wts, call = rlang::caller_env()) {
  meta <- density_meta(wts)

  if (is.null(meta)) {
    return(list(
      density = dens_normal(),
      numerator = ipw_continuous_numerator(
        is_stabilized(wts),
        stabilization_score(wts)
      ),
      sigma = list(kind = "pooled", value = NULL)
    ))
  }

  check_ipw_continuous_density(meta$density, call = call)

  list(
    density = meta$density,
    numerator = meta$numerator,
    sigma = ipw_continuous_spread(meta, call = call)
  )
}

# The spread the stacked system reads the conditional density at. A pooled
# spread is the residual moment the system estimates alongside the coefficients;
# a single spread the caller supplied is a constant it holds fixed, carrying
# none of its uncertainty, which is what fixing it says. A spread supplied for
# each observation is neither: it is a function of the data that no parameter
# value here reproduces, so it is refused before anything is solved rather than
# reported afterwards as two vectors that disagree.
ipw_continuous_spread <- function(meta, call = rlang::caller_env()) {
  if (identical(meta$sigma, "pooled")) {
    return(list(kind = "pooled", value = NULL))
  }

  if (!is.null(meta$sigma_value)) {
    return(list(kind = "fixed", value = meta$sigma_value))
  }

  abort(
    c(
      "{.fun ipw} does not support weights built with an observation-level \\
      {.arg .sigma}.",
      x = "The weights supplied to {.arg outcome_mod} record a spread for \\
      each observation, which has no counterpart in the stacked system: it \\
      estimates one conditional spread, or holds one fixed.",
      i = "Rebuild the weights with the pooled default, by leaving \\
      {.arg .sigma} unset, or with a single {.arg .sigma}, which \\
      {.fun ipw} takes as a known constant."
    ),
    error_class = "propensity_ipw_sigma_error",
    call = call
  )
}

# A kernel estimate is fit to the residuals of the propensity score model, and
# its bandwidth is chosen from them, so the density is not a function of the
# parameters that the sandwich could differentiate. The refusal is the one a
# `gam` gets: what is unavailable is the standard error method rather than the
# weights, which are exactly what they claim to be.
check_ipw_continuous_density <- function(density, call = rlang::caller_env()) {
  if (!identical(density$family, "kernel")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} cannot build a sandwich variance for weights built with a \\
      {.val kernel} density.",
      x = "The bandwidth of a kernel estimate is chosen from the residuals of \\
      the propensity score model, so the weights are not a differentiable \\
      function of that model's parameters.",
      i = ipw_continuous_bootstrap_hint()
    ),
    error_class = c(
      "propensity_ipw_se_method_unavailable_error",
      "propensity_method_error"
    ),
    call = call
  )
}

# The remedy both unavailable-method refusals point to. Resampling asks the
# weights for nothing but their value, so it reaches the fits a stacked score
# cannot describe.
ipw_continuous_bootstrap_hint <- function() {
  "Use {.code se_method = \"bootstrap\"}, passing the data the models were fit
   to in {.arg .data}, which resamples the whole fit instead of stacking it."
}
