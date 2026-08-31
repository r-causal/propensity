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
#
# `hint` is the remedy every refusal here ends on, which the route asking for
# the entry supplies. The package computes no standard error for these fits, so
# the remedy is a bootstrap of the whole pipeline written by hand; a route where
# even that does not apply supplies a hint of its own.
#
# `label` is what a refusal calls the model it refused. The single-dose route
# reads it out of `wt_mod`, which is the argument the model arrived in; the
# joint route reads a dose out of a container of two models, where `wt_mod`
# names the container and the component has a name of its own.
ipw_continuous_model <- function(
  ps_mod,
  hint = ipw_continuous_resample_hint(),
  label = "wt_mod",
  role = "propensity score model",
  call = rlang::caller_env()
) {
  ps_class <- class(ps_mod)

  if (inherits(ps_mod, "gam")) {
    return(ipw_continuous_gam_entry(
      ps_mod,
      ps_class,
      hint = hint,
      label = label,
      role = role,
      call = call
    ))
  }

  if (inherits(ps_mod, "rlm")) {
    return(ipw_continuous_rlm_entry(
      ps_mod,
      ps_class,
      hint = hint,
      label = label,
      role = role
    ))
  }

  if (inherits(ps_mod, "glm")) {
    # The stacked score is the one deli writes for a normal model, so a family
    # whose spread changes with its fitted values is refused before its link is
    # read: there is no single conditional density for the weights to divide by
    # and no equation here that its coefficients solve.
    check_ipw_continuous_ps_family(
      ps_mod,
      label = label,
      role = role,
      call = call
    )

    link <- ps_mod$family$link
    check_ipw_continuous_ps_link(
      link,
      label = label,
      role = role,
      call = call
    )

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
      {.fun stats::glm}, {.fun MASS::rlm}, or {.fun mgcv::gam} as the {role} \\
      of a continuous exposure.",
      x = "{.arg {label}} has class {.cls {ps_class}}.",
      i = "Each of those is read as the class it is rather than by what it \\
      inherits from, so a subclass of one of them reaches this refusal too.",
      i = "Refit {.arg {label}} with {.fun stats::lm} or \\
      {.code stats::glm(family = gaussian())}."
    ),
    error_class = "propensity_class_error",
    call = call
  )
}

# The entry an additive fit gets. `mgcv::gam()` is read as the penalized least
# squares it is, with its smoothing parameters held at the values the fit
# reports: its design is the smooth basis it was built on rather than the
# columns its formula names, and the equation its coefficients solve is the
# least-squares score of that basis less the penalty those smoothing parameters
# define.
#
# Holding the smoothing parameters fixed is what makes the score writable at
# all, and it is also the approximation. The system conditions on the smoothing
# choice rather than estimating it, so the variance it reports carries none of
# that choice's uncertainty, and what it reports for the coefficients is the
# frequentist covariance of a penalized fit rather than the Bayesian one mgcv
# reports by default. The family and link envelope is the one every other
# gaussian fit here has, and it is checked in the same place.
ipw_continuous_gam_entry <- function(
  ps_mod,
  classes = class(ps_mod),
  hint = ipw_continuous_resample_hint(),
  label = "wt_mod",
  role = "propensity score model",
  call = rlang::caller_env()
) {
  check_ipw_continuous_ps_family(
    ps_mod,
    label = label,
    role = role,
    call = call
  )

  link <- ps_mod$family$link
  check_ipw_continuous_ps_link(
    link,
    label = label,
    role = role,
    call = call
  )

  smoothing <- ipw_gam_smoothing_parameters(ps_mod)

  # A fit holding a penalty this path cannot place is refused rather than
  # stacked without it: a penalty left out of the score moves its root, and the
  # solve would walk the coefficients off the fit the user has.
  if (is.null(smoothing)) {
    return(ipw_continuous_entry(
      kind = "gam",
      classes = classes,
      link = link,
      stackable = FALSE,
      label = label,
      role = role,
      error_class = c(
        "propensity_ipw_se_method_unavailable_error",
        "propensity_method_error"
      ),
      reason = c(
        "{.fun ipw} reads the penalty of a {.cls {entry$classes}} \\
        {entry$role} off the smooth terms it was fit with.",
        x = "{.arg {entry$label}} records more smoothing parameters than its \\
        smooth terms account for, which is what a penalty on a parametric \\
        term, such as one from {.arg paraPen}, adds.",
        i = "Refit {.arg {entry$label}} so that every penalty in it belongs \\
        to a smooth term.",
        i = hint
      )
    ))
  }

  penalty <- ipw_gam_penalty(ps_mod, smoothing)

  # A fit whose own coefficients are not at the root of the score built from its
  # record is refused for the same reason, one step later and on the arithmetic
  # rather than on a count. `min.sp` is the shape that reaches here: it puts a
  # floor under the smoothing parameters, adds that floor to the penalty, and
  # records the smoothing parameters as excluding it, so the penalty rebuilt
  # from the record is not the one the fit was made under and nothing in the
  # fitted object says so. The score is the check that no field is.
  #
  # The tolerance is on the score divided by the number of observations, so it
  # reads the same at any sample size. The fits this route stacks sit five or
  # more orders of magnitude below it, and a floor small enough to move the
  # smoothing parameters in their fourth significant digit already sits above
  # it at the sample size the fixtures use. That last margin is the one thing
  # here that is not a property of the check: what passes is a floor whose
  # effect on the stacked fit is proportionally below the tolerance, and at a
  # larger sample size a floor of the same relative size can fall under it.
  gap <- ipw_gam_score_gap(ps_mod, link, penalty)

  if (!isTRUE(gap < ipw_gam_score_tolerance)) {
    return(ipw_continuous_entry(
      kind = "gam",
      classes = classes,
      link = link,
      stackable = FALSE,
      label = label,
      role = role,
      error_class = c(
        "propensity_ipw_se_method_unavailable_error",
        "propensity_method_error"
      ),
      reason = c(
        "{.fun ipw} stacks a {.cls {entry$classes}} {entry$role} at the \\
        penalized score its coefficients solve, which it rebuilds from the \\
        penalty the fit records.",
        x = "The penalty {.arg {entry$label}} records does not reproduce the \\
        score {.arg {entry$label}} is at, so a system seeded at its \\
        coefficients would settle away from them.",
        i = "A smoothing floor from {.arg min.sp} is one cause: it is added to \\
        the penalty and left out of the smoothing parameters the fit reports. \\
        Refit {.arg {entry$label}} without {.arg min.sp}.",
        i = hint
      )
    ))
  }

  ipw_continuous_entry(
    kind = "gam",
    link = link,
    stackable = TRUE,
    penalty = penalty
  )
}

# How far a fitted additive model's coefficients sit from the root of the
# penalized score rebuilt from its record, per observation:
#
#   max |X'W(a - mu(X alpha)) - S_lambda alpha| / n.
#
# The data half is the score this route stacks, read at the fit's own design and
# response, and the penalty half is the one built from the fit's smoothing
# parameters. A fit whose record describes the penalty it was made under leaves
# the two at the same place, and the difference is the arithmetic of the solve
# rather than a quantity.
#
# The prior weights are the fit's own, so what is asked here is whether the
# record reproduces the fit's score rather than whether the fit is one this
# route can stack in every other respect: a fit made under case weights is
# refused by name further along, where the remedy is its own.
#
# A fit that cannot be read at all, which is one keeping neither its design nor
# its response, returns a gap of `NA` and is refused with the rest: a score that
# cannot be checked is not one to stack on.
ipw_gam_score_gap <- function(fit, link, penalty) {
  alpha <- stats::coef(fit)

  gap <- tryCatch(
    {
      X <- stats::model.matrix(fit)
      y <- as.numeric(fit$y)
      weights <- fit$prior.weights
      if (is.null(weights)) {
        weights <- rep(1, length(y))
      }

      score <- ipw_continuous_link_fns(link)$score
      total <- as.vector(score(alpha, X = X, y = y) %*% weights)
      max(abs(total - as.vector(penalty %*% alpha))) / length(y)
    },
    error = function(e) NA_real_
  )

  if (!is.numeric(gap) || !is.finite(gap)) {
    return(NA_real_)
  }

  gap
}

# How close to the root of its own penalized score a fitted additive model's
# coefficients have to sit, per observation, to be stacked. The fits this route
# reads solve their score to rounding and the shapes it refuses miss it by
# orders of magnitude, so the tolerance is placed between the two with room on
# both sides rather than at either edge.
#
# Room on both sides is not the same as a bright line: what passes is a floor
# whose effect on the stacked fit is proportionally below the tolerance, which
# depends on the sample size the gap is read over as well as on the floor.
ipw_gam_score_tolerance <- 1e-6

# The smoothing parameters of a fitted additive model, one for each penalty
# matrix its smooth terms carry, or `NULL` for a fit whose record of them does
# not line up with those terms.
#
# A fit records every penalty it holds in `full.sp` and only the ones it chose
# in `sp`, so the two agree until a smoothing parameter is handed to the fit and
# `full.sp` is the field that describes the whole penalty. It is therefore the
# one read wherever the fit keeps it, and a count of it that the smooth terms do
# not account for is what a penalty on a parametric term, such as one from
# `paraPen`, leaves behind: `sp` can still be as long as those terms are, and
# reading it there would place the wrong parameters on the smooths and drop the
# parametric penalty from the score. `sp` is read only for a fit with no full
# record at all, which is one carrying no penalty.
#
# A smooth fit with `fx = TRUE` carries no penalty matrix and no smoothing
# parameter, and one selected out carries several of each, so the count comes
# from the terms rather than from the number of them.
ipw_gam_smoothing_parameters <- function(fit) {
  n_penalties <- sum(vapply(
    fit$smooth,
    function(smooth) length(smooth$S),
    integer(1)
  ))

  smoothing <- if (length(fit$full.sp) > 0) fit$full.sp else fit$sp

  if (identical(length(smoothing), n_penalties)) {
    return(as.numeric(smoothing))
  }

  NULL
}

# The total penalty of a fitted additive model, in the parameterization its
# coefficients are reported in: each smooth term's penalty matrices scaled by
# their own smoothing parameters and placed at the coefficients that term owns.
# The smoothing parameters are passed in rather than read here, because a fit
# whose record of them does not line up with its smooth terms has no penalty
# this can place and has already been refused by the time this is reached.
ipw_gam_penalty <- function(fit, smoothing) {
  p <- length(stats::coef(fit))
  penalty <- matrix(0, nrow = p, ncol = p)
  k <- 0L

  for (smooth in fit$smooth) {
    idx <- smooth$first.para:smooth$last.para
    for (S in smooth$S) {
      k <- k + 1L
      penalty[idx, idx] <- penalty[idx, idx] + smoothing[[k]] * S
    }
  }

  penalty
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
ipw_continuous_rlm_entry <- function(
  ps_mod,
  classes = class(ps_mod),
  hint = ipw_continuous_resample_hint(),
  label = "wt_mod",
  role = "propensity score model"
) {
  psi_name <- ipw_rlm_psi_name(ps_mod)
  huber <- identical(psi_name, "MASS::psi.huber")
  mm <- identical(ipw_rlm_method(ps_mod), "MM")

  if (mm || !huber) {
    found <- if (mm && !is.null(psi_name)) {
      "{.arg {entry$label}} was fit with {.code method = \"MM\"}, which starts from a
       high-breakdown fit and finishes on {.fun {entry$psi}}."
    } else if (mm) {
      "{.arg {entry$label}} was fit with {.code method = \"MM\"}, which starts from a
       high-breakdown fit and finishes on a redescending psi."
    } else if (is.null(psi_name)) {
      "{.arg {entry$label}} was fit with a psi function this path cannot recognize."
    } else {
      "{.arg {entry$label}} was fit with {.fun {entry$psi}}."
    }

    return(ipw_continuous_entry(
      kind = "rlm",
      classes = classes,
      link = "identity",
      stackable = FALSE,
      psi = psi_name,
      label = label,
      role = role,
      error_class = "propensity_ipw_robust_psi_error",
      reason = c(
        "{.fun ipw} stacks only the Huber score of a {.cls {entry$classes}} \\
        {entry$role} of a continuous exposure.",
        x = found,
        i = "Refit {.arg {entry$label}} with {.fun MASS::psi.huber}, the default, \\
        whose threshold {.fun ipw} reads off the fit.",
        i = hint
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
      label = label,
      role = role,
      error_class = "propensity_ipw_convergence_error",
      reason = c(
        "{.fun ipw} cannot stack a {.cls {entry$classes}} {entry$role} that \\
        did not converge.",
        x = "{.arg {entry$label}} reports {.code converged = FALSE}, so its \\
        coefficients are not the root of the score stacked here.",
        i = "Refit {.arg {entry$label}} with a larger {.arg maxit}, or a looser \\
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
# that the psi builder writes the same equation the registry did. `penalty` is
# the one an additive fit carries, and is read the same way: the fixed matrix
# its score subtracts, held in the spec rather than recomputed from the model.
ipw_continuous_entry <- function(
  kind,
  link,
  stackable,
  classes = kind,
  error_class = NULL,
  reason = NULL,
  psi = NULL,
  huber_k = NULL,
  penalty = NULL,
  label = "wt_mod",
  role = "propensity score model"
) {
  entry <- list(
    kind = kind,
    classes = classes,
    label = label,
    role = role,
    link = link,
    stackable = stackable,
    error_class = error_class,
    reason = reason,
    psi = psi,
    huber_k = huber_k,
    penalty = penalty,
    mean = NULL,
    score = NULL
  )

  if (!stackable) {
    return(entry)
  }

  utils::modifyList(
    entry,
    ipw_continuous_score_fns(kind, link, huber_k, penalty)
  )
}

# The mean and the score a stacked propensity score model contributes, keyed on
# the registry entry's kind as well as its link. Every least squares fit is read
# through its link alone, and a robust fit is the one kind whose score is a
# different equation at the same identity link, so the kind is what separates
# them. A penalty is the one thing that changes a least squares score without
# changing its link, so it is carried beside them rather than keyed on.
ipw_continuous_score_fns <- function(
  kind,
  link,
  huber_k = NULL,
  penalty = NULL
) {
  if (identical(kind, "rlm")) {
    return(list(
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
    ))
  }

  fns <- ipw_continuous_link_fns(link)

  if (is.null(penalty)) {
    return(fns)
  }

  ipw_penalized_score_fns(fns, penalty)
}

# The same pair with a fixed penalty subtracted from the score, which is the
# equation an additive fit's coefficients are the root of. The penalty is one
# term of the whole equation rather than a sum over observations, so it is
# spread evenly across the columns of the psi matrix: the column sums are then
# the equation the fit solves, and the bread picks up the penalized
# cross-product divided by the sample size, as every other block is.
ipw_penalized_score_fns <- function(fns, penalty) {
  unpenalized <- fns$score

  fns$score <- function(alpha, X, y) {
    unpenalized(alpha, X, y) - as.vector(penalty %*% alpha) / length(y)
  }

  fns
}

# The same pair read off a spec, which is what the psi builder and the preflight
# that has to agree with it both hold. The spec records the kind, the link, the
# robust threshold, and the penalty rather than the model they came from, so the
# two rebuild the propensity score block from the same four values.
ipw_continuous_spec_fns <- function(spec) {
  ipw_continuous_score_fns(
    spec$ps$kind,
    spec$ps$link,
    spec$ps$huber_k,
    spec$ps$penalty
  )
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
  label = "wt_mod",
  role = "propensity score model",
  call = rlang::caller_env()
) {
  family <- ps_mod$family
  if (identical(family$family, "gaussian")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} supports only a gaussian {role} for a continuous exposure.",
      x = "{.arg {label}} was fit with \\
      {.code {model_family_label(family)}}.",
      i = "Refit {.arg {label}} with {.fun stats::lm} or \\
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

check_ipw_continuous_ps_link <- function(
  link,
  label = "wt_mod",
  role = "propensity score model",
  call = rlang::caller_env()
) {
  if (link %in% ipw_continuous_ps_links) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support the {.val {link}} link for the {role} of a \\
      continuous exposure.",
      x = "{.arg {label}} is a gaussian model fit with the {.val {link}} \\
      link.",
      i = "Refit {.arg {label}} with one of {.val {ipw_continuous_ps_links}}, \\
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
# density spread by the pooled residual root mean square, over whatever their
# stabilization says stabilized them.
#
# `stacked` says whether the caller is going to differentiate the ratio. A
# kernel density is refused only for the caller that is, since a route that
# rebuilds the weights by calling `wt_ate()` again asks the density for nothing
# a kernel cannot give.
ipw_continuous_ratio <- function(
  wts,
  stacked = TRUE,
  call = rlang::caller_env()
) {
  ipw_continuous_ratio_meta(
    density_meta(wts),
    stabilized = is_stabilized(wts),
    score = stabilization_score(wts),
    stacked = stacked,
    call = call
  )
}

# The same reading from the record itself, for the route whose weights carry
# theirs somewhere other than on the vector. A product weight records one
# density per component, so the dose's record is read out of `joint_wt_meta()`
# and the stabilization that stands in for a missing one is read from there too.
ipw_continuous_ratio_meta <- function(
  meta,
  stabilized,
  score = NULL,
  stacked = TRUE,
  hint = ipw_continuous_resample_hint(),
  call = rlang::caller_env()
) {
  if (is.null(meta)) {
    return(list(
      density = dens_normal(),
      numerator = ipw_continuous_numerator(stabilized, score),
      sigma = list(kind = "pooled", value = NULL)
    ))
  }

  if (stacked) {
    check_ipw_continuous_density(meta$density, hint = hint, call = call)
  }

  list(
    density = meta$density,
    numerator = meta$numerator,
    # The model a numerator was estimated with, which the stacked system
    # estimates alongside everything else rather than reading the numerator it
    # evaluated to.
    numerator_model = meta$numerator_model,
    sigma = ipw_continuous_spread(meta, call = call)
  )
}

# The stabilization block a numerator model contributes: the design its
# coefficients multiply, the score they solve, and the link its mean is read
# back through. The model goes through the registry the propensity score model
# goes through, so a class whose score this stack cannot write is refused with
# that registry's own reason, naming the argument the model arrived in.
ipw_numerator_model_block <- function(
  numerator_model,
  call = rlang::caller_env()
) {
  if (is.null(numerator_model)) {
    return(NULL)
  }

  entry <- ipw_continuous_model(
    numerator_model,
    label = "stabilize",
    role = "numerator model",
    call = call
  )
  check_ipw_continuous_model(entry, call = call)

  # The block below multiplies the fitted coefficients against the design
  # positionally, as the propensity score block does, so the model needs a
  # coefficient for every column of its design. Without this the numerator comes
  # back missing at every value of theta, and the first report of it is the
  # weights the system rebuilds failing to match the ones the caller built.
  check_ipw_model_rank(stats::coef(numerator_model), "stabilize", call = call)

  list(
    # An additive fit's design is the smooth basis it reports rather than the
    # columns its formula names, which `model.matrix()` returns for it as well.
    X = stats::model.matrix(numerator_model),
    kind = entry$kind,
    link = entry$link,
    huber_k = entry$huber_k,
    penalty = entry$penalty,
    coefs = stats::coef(numerator_model)
  )
}

# The mean and the score the numerator model's block contributes, read off the
# block the way the propensity score block's are read off the spec.
ipw_numerator_model_fns <- function(model) {
  ipw_continuous_score_fns(
    model$kind,
    model$link,
    model$huber_k,
    model$penalty
  )
}

# What the stacked system needs of a numerator model beyond a score it can
# write: that the score is the score of the exposure being weighted, and that it
# is written over the same observations as everything else the system stacks.
# The block's own equations read the exposure the propensity score model reads,
# so a model of another response would sit away from the root of the row seeded
# for it and the solve would move it, reporting a numerator the user never fit.
#
# What is checked is the name of the response and how many observations it was
# read over, not that those observations are the same rows in the same order.
# Row identity is what the preflight at initialization covers: weights rebuilt
# over misaligned rows do not match the ones the caller supplied, and that
# comparison reports the misalignment.
check_ipw_numerator_model <- function(
  numerator_model,
  block,
  exposure_name,
  n,
  call = rlang::caller_env()
) {
  response <- tryCatch(
    fmla_extract_left_chr(numerator_model),
    error = function(e) NA_character_
  )

  if (!identical(response, exposure_name)) {
    abort(
      c(
        "The model supplied to {.arg stabilize} must model the exposure.",
        x = "It models {.val {response}} and {.arg wt_mod} models
             {.val {exposure_name}}.",
        i = "The numerator of the weights is the density of the exposure given
             what the numerator model reads, so both models describe the same
             response."
      ),
      error_class = "propensity_ipw_numerator_error",
      call = call
    )
  }

  if (!identical(nrow(block$X), as.integer(n))) {
    abort(
      c(
        "The model supplied to {.arg stabilize} must be fit to the observations
         the other models were fit to.",
        x = "It was fit to {nrow(block$X)} observation{?s} and
             {.arg outcome_mod} to {n}.",
        i = "Refit the numerator model on the data the other models were fit
             to, and rebuild the weights from it."
      ),
      error_class = "propensity_ipw_numerator_error",
      call = call
    )
  }

  invisible(TRUE)
}

# The spread the stacked system reads the conditional density at. A pooled
# spread is the residual moment the system estimates alongside the coefficients,
# and a scale fit by maximum likelihood is the score of the t estimated in the
# same place; a single spread the caller supplied is a constant it holds fixed,
# carrying none of its uncertainty, which is what fixing it says. A spread
# supplied for each observation is neither: it is a function of the data that no
# parameter value here reproduces, so it is refused before anything is solved
# rather than reported afterwards as two vectors that disagree.
ipw_continuous_spread <- function(meta, call = rlang::caller_env()) {
  if (identical(meta$sigma, "pooled")) {
    return(list(kind = "pooled", value = NULL))
  }

  if (identical(meta$sigma, "mle")) {
    return(list(kind = "mle", value = NULL))
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
# parameters that the sandwich could differentiate. What is unavailable is the
# standard error method rather than the weights, which are exactly what they
# claim to be, so the refusal says so and points at a bootstrap.
check_ipw_continuous_density <- function(
  density,
  hint = ipw_continuous_resample_hint(),
  call = rlang::caller_env()
) {
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
      i = hint
    ),
    error_class = c(
      "propensity_ipw_se_method_unavailable_error",
      "propensity_method_error"
    ),
    call = call
  )
}

# The remedy both unavailable-method refusals point to. The package computes no
# standard error for these fits, but resampling asks the weights for nothing but
# their value, so a bootstrap written by hand reaches what a stacked score
# cannot describe.
ipw_continuous_resample_hint <- function() {
  "propensity has no resampling method; bootstrap the whole fit yourself:
   resample the rows, refit the propensity score model, rebuild the weights with
   {.fn wt_ate}, and refit the outcome model on each resample."
}
