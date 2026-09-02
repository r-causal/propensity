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

  # The design is the smooth basis the fit has to evaluate, which is the larger
  # part of what reading an additive fit costs, and it is the same matrix
  # wherever it is read: the score checked below and the psi function's own
  # multiply are both written against it. It is therefore built once here and
  # travels on the entry to every reader.
  design <- ipw_gam_design(ps_mod)

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
  gap <- ipw_gam_score_gap(ps_mod, link, penalty, X = design)

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
    penalty = penalty,
    design = design
  )
}

# The design a fitted additive model is read at: the smooth basis it was built
# on, evaluated over the rows it was fit to. `mgcv` returns it for
# `model.matrix()` as for `predict(type = "lpmatrix")`, since they are the same
# call.
#
# A fit that cannot be read at all, which is one keeping neither its design nor
# the frame to rebuild it from, has no design here and is refused further along
# by the score it leaves unreadable.
ipw_gam_design <- function(fit) {
  tryCatch(stats::model.matrix(fit), error = function(e) NULL)
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
# The design is the one the entry built, since the score is read at the same
# basis the psi function multiplies later.
#
# A fit that cannot be read at all, which is one keeping neither its design nor
# its response, returns a gap of `NA` and is refused with the rest: a score that
# cannot be checked is not one to stack on. That branch is defensive rather than
# reachable from any fitting call known here: mgcv keeps the model frame even
# under `model = FALSE`, so `model.matrix()` and the response are there for
# every fit this route sees.
ipw_gam_score_gap <- function(fit, link, penalty, X = ipw_gam_design(fit)) {
  if (is.null(X)) {
    return(NA_real_)
  }

  alpha <- stats::coef(fit)

  gap <- tryCatch(
    {
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

# The entry a robust fit gets, which is the one place a class carries constants
# of its own into the stacked system. `MASS::rlm()` minimizes a loss rather than
# the sum of squares, so its coefficients are the root of the psi score read at
# the scale the fit settled on, and that score is stacked here as the deli
# robust regression loss the psi corresponds to.
#
# Each of the three psi functions MASS supplies has such a loss: `psi.huber` is
# deli's `"huber"`, `psi.bisquare` is its `"tukey"`, and `psi.hampel` is its
# `"hampel"`. The two clip the same residual on different scales: `rlm` clips
# the residual divided by its scale estimate at the psi's own constants, and
# deli clips the raw residual at `k`, so the equations agree when each constant
# is multiplied by the scale the fit settled on. How many constants there are
# differs by psi, one for Huber and for the bisquare and three for Hampel, so
# what the entry carries is a vector rather than a number.
#
# The constants are read out of the psi's formals rather than out of `k2`,
# because that is where `rlm` writes a caller's choice; `k2` tunes the scale
# estimator instead and is unchanged when the psi's thresholds move. Which
# formals hold them differs by psi as well: `k` for Huber, `c` for the bisquare,
# and `a`, `b`, and `c` for Hampel.
#
# The bisquare and Hampel psi functions redescend, so the equation stacked for
# them has more than one root and the solve reports whichever one it is seeded
# at. The seed is the fit's own coefficients, so what the sandwich describes is
# the covariance of that root, read locally, and it carries no claim about the
# others.
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
  mm <- identical(ipw_rlm_method(ps_mod), "MM")

  if (mm || is.null(psi_name)) {
    found <- if (mm && !is.null(psi_name)) {
      "{.arg {entry$label}} was fit with {.code method = \"MM\"}, whose
       high-breakdown start decides which root of {.fun {entry$psi}} the fit
       finishes at and supplies the scale it clips at, and neither of those is
       an equation this system writes, so a sandwich read at those coefficients
       would not describe how the fit behaves across samples."
    } else if (mm) {
      "{.arg {entry$label}} was fit with {.code method = \"MM\"}, whose
       high-breakdown start decides which root the fit finishes at and supplies
       the scale it clips at, and neither of those is an equation this system
       writes."
    } else {
      "{.arg {entry$label}} was fit with a psi function this path cannot recognize."
    }

    remedy <- if (mm) {
      "Refit {.arg {entry$label}} with the default {.code method = \"M\"}, whose
       psi score {.fun ipw} writes."
    } else {
      "Refit {.arg {entry$label}} with {.fun MASS::psi.huber},
       {.fun MASS::psi.bisquare}, or {.fun MASS::psi.hampel}, whose constants
       {.fun ipw} reads off the fit."
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
        "{.fun ipw} cannot write the equation this {.cls {entry$classes}} \\
        {entry$role} of a continuous exposure is the root of.",
        x = found,
        i = remedy,
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
    psi_loss = ipw_rlm_psi_losses[[psi_name]],
    psi_k = ipw_rlm_psi_constants(ps_mod, psi_name)
  )
}

# Which deli robust regression loss each of MASS's psi functions is, keyed on
# the name `ipw_rlm_psi_name()` reports. A psi with no entry here is one whose
# score this path does not write.
ipw_rlm_psi_losses <- c(
  "MASS::psi.huber" = "huber",
  "MASS::psi.bisquare" = "tukey",
  "MASS::psi.hampel" = "hampel"
)

# The constants a robust fit clips the raw residual at, which is what deli's
# loss takes: the constants the psi clips the standardized residual at, each
# multiplied by the scale the fit settled on. Which formals hold them is the
# psi's own naming, and `rlm` records a caller's choice by rewriting exactly
# those formals, so the fit's psi is where they are read from.
ipw_rlm_psi_constants <- function(ps_mod, psi_name) {
  psi_formals <- formals(ps_mod$psi)

  constants <- switch(
    psi_name,
    "MASS::psi.huber" = psi_formals$k,
    "MASS::psi.bisquare" = psi_formals$c,
    "MASS::psi.hampel" = c(psi_formals$a, psi_formals$b, psi_formals$c)
  )

  as.numeric(constants) * ps_mod$s
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
# `psi`, `psi_loss`, and `psi_k` are the three things a robust fit carries that
# no other entry has: the name of the psi it was fit with, which its refusal
# reads, the deli loss that psi's score is written as, and the constants that
# loss clips the raw residual at. The spec carries the last two so that the psi
# builder writes the same equation the registry did. `penalty` is the one an
# additive fit carries, and is read the same way: the fixed matrix its score
# subtracts, held in the spec rather than recomputed from the model.
#
# `design` is the second thing an additive fit carries, for the reason it
# carries the penalty and one more: evaluating a smooth basis is expensive
# enough that every reader building its own is most of what the route costs.
# The entry is what those readers have in common, so the design it built to
# check the fit's score travels on it, and a reader holding an entry holds the
# design as well. A fit of any other kind has a design that is a lookup, so
# only the additive one carries one here.
ipw_continuous_entry <- function(
  kind,
  link,
  stackable,
  classes = kind,
  error_class = NULL,
  reason = NULL,
  psi = NULL,
  psi_loss = NULL,
  psi_k = NULL,
  penalty = NULL,
  design = NULL,
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
    psi_loss = psi_loss,
    psi_k = psi_k,
    penalty = penalty,
    design = design,
    mean = NULL,
    score = NULL
  )

  if (!stackable) {
    return(entry)
  }

  utils::modifyList(
    entry,
    ipw_continuous_score_fns(kind, link, psi_loss, psi_k, penalty)
  )
}

# The mean and the score a stacked propensity score model contributes, keyed on
# the registry entry's kind as well as its link. Every least squares fit is read
# through its link alone, and a robust fit is the one kind whose score is a
# different equation at the same identity link, so the kind is what separates
# them. Which equation that is depends on the psi the fit descended, so the
# loss and its constants travel with the kind. A penalty is the one thing that
# changes a least squares score without changing its link, so it is carried
# beside them rather than keyed on.
ipw_continuous_score_fns <- function(
  kind,
  link,
  psi_loss = NULL,
  psi_k = NULL,
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
          k = psi_k,
          loss = psi_loss
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
# robust loss and its constants, and the penalty rather than the model they came
# from, so the two rebuild the propensity score block from the same five values.
ipw_continuous_spec_fns <- function(spec) {
  ipw_continuous_score_fns(
    spec$ps$kind,
    spec$ps$link,
    spec$ps$psi_loss,
    spec$ps$psi_k,
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
#
# `n` is the number of observations the ratio is about to be read over, which is
# what the score the record names is checked against.
ipw_continuous_ratio <- function(
  wts,
  n,
  stacked = TRUE,
  call = rlang::caller_env()
) {
  score <- stabilization_score(wts)

  ratio <- ipw_continuous_ratio_meta(
    density_meta(wts),
    stabilized = is_stabilized(wts),
    score = score,
    stacked = stacked,
    call = call
  )

  check_ipw_stabilization_score(ratio$numerator, score, n, call = call)

  ratio
}

# Refuse weights whose record names a score as their numerator when the score in
# hand is not one value for each observation being weighted. The record
# describes how the weights were built rather than what survived reaching here,
# and two things break the agreement between the two.
#
# Subsetting the weights drops a per-observation score outright, since nothing
# rebuilding a psw is given the indices behind a length change; the record still
# says a score stabilized them and there is no vector left to read.
# `stats::model.frame()` instead re-attaches the original attributes to the
# shortened vector, so the score arrives at the length the weights were recorded
# at rather than the length they now are.
#
# Either way the ratio would divide a density read over the analyzed rows by a
# numerator describing other rows: zero-length in the first case, recycling in
# the second. Neither is arithmetic the caller wrote, and the report that
# followed named the estimand, the spread, and `.data` rather than the score, so
# the score is checked here, before any ratio is built.
#
# A scalar score is exempt for the reason it is the remedy: one number scales
# every weight whatever the rows are, so nothing about it is indexed and nothing
# drops it.
#
# `component` names the component of a product weight the score belongs to, and
# is `NULL` for the single-treatment route, whose weights carry the one score
# they have on the vector itself.
check_ipw_stabilization_score <- function(
  numerator,
  score,
  n,
  component = NULL,
  call = rlang::caller_env()
) {
  if (!identical(numerator, "score")) {
    return(invisible(TRUE))
  }

  # A score with no values in it is refused alongside a missing one.
  # `stabilization_score_aligns()` passes a zero-length score, because it is
  # written for the prototypes and empty subsets a `psw` restore passes through,
  # where a score lines up with observations that have not arrived yet. Here the
  # observations have arrived and the record names a score as their numerator,
  # so a score holding nothing is a numerator with nothing to read.
  if (length(score) > 0 && stabilization_score_aligns(score, n)) {
    return(invisible(TRUE))
  }

  # Both halves are cli templates rather than values spliced into one, because
  # cli reads the message once: a brace inside a value is text by the time it
  # arrives.
  headline <- if (is.null(component)) {
    "{.fun ipw} can't read the {.arg stabilization_score} the weights \\
    supplied to {.arg outcome_mod} record as their numerator."
  } else {
    "{.fun ipw} can't read the {.arg stabilization_score} the \\
    {.arg {component}} component of the weights supplied to \\
    {.arg outcome_mod} records as its numerator."
  }

  held <- if (is.null(score)) {
    "hold no score at all."
  } else {
    "hold {length(score)} of them."
  }

  abort(
    c(
      headline,
      x = paste("They weight {n} observation{?s} and", held),
      i = "A per-observation score is one value for each unit, so it does not \\
      survive the rows being restricted: subsetting the weights drops it, and \\
      a model frame that drops incomplete rows leaves it at the length the \\
      weights were built at.",
      i = "Rebuild the weights on the rows being analyzed, or stabilize on a \\
      single {.arg stabilization_score}, which scales every weight and \\
      survives any restriction."
    ),
    error_class = "propensity_ipw_stabilization_score_error",
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
#
# `.data` is the frame every other design in the stack is rebuilt from, over the
# rows every model read, and this design is rebuilt from it alongside them.
# Without one the design is read off the fit itself.
#
# `component` names the component whose weights the numerator was built for on
# the joint route, where naming the argument alone would leave a caller with two
# stabilized components unable to tell which fit is refused.
ipw_numerator_model_block <- function(
  numerator_model,
  .data = NULL,
  component = NULL,
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
  check_ipw_numerator_model_weights(
    numerator_model,
    component = component,
    call = call
  )

  # The block below multiplies the fitted coefficients against the design
  # positionally, as the propensity score block does, so the model needs a
  # coefficient for every column of its design. Without this the numerator comes
  # back missing at every value of theta, and the first report of it is the
  # weights the system rebuilds failing to match the ones the caller built.
  check_ipw_model_rank(
    stats::coef(numerator_model),
    "stabilize",
    component = component,
    call = call
  )

  numerator_design <- ipw_numerator_model_design(
    numerator_model,
    entry,
    .data = .data,
    component = component,
    call = call
  )

  list(
    # The numerator model's own design, whichever way it was arrived at: a block
    # multiplied by another model's design would rebuild a numerator nobody fit.
    X = numerator_design,
    kind = entry$kind,
    link = entry$link,
    psi_loss = entry$psi_loss,
    psi_k = entry$psi_k,
    penalty = entry$penalty,
    coefs = stats::coef(numerator_model),
    # The spread the numerator's density was read at, taken the way
    # `numerator_model_moments()` takes it: the mean square of the fit's own
    # response-scale residuals, over the fit's own rows. The design above is
    # rebuilt over the rows `.data` leaves, and the same moment read over those
    # rows is a different spread from the one the weights carry.
    sigma2_fit = mean(
      as.numeric(stats::residuals(numerator_model, type = "response"))^2,
      na.rm = TRUE
    )
  )
}

# The design a numerator model's block multiplies its coefficients against.
#
# With a `.data` the caller supplied it is rebuilt from that frame under the
# fit's own terms, contrasts, and levels, which is what puts it over the rows
# every other design in the stack is built over. Without one it is read off the
# fit: an additive fit's entry has already evaluated its smooth basis, and every
# other fit is asked for its own model matrix.
#
# An `lm` or `glm` usually keeps its model frame, but one fit with
# `model = FALSE` keeps none and rebuilds it by re-evaluating the fitting call,
# which a fit made inside a function whose frame is gone cannot do. That is the
# denominator's recovery and the failure it can meet, so it is reported the way
# the denominator's is, named for the argument this model arrived in and, on the
# joint route, for the component it was built for.
ipw_numerator_model_design <- function(
  numerator_model,
  entry,
  .data = NULL,
  component = NULL,
  call = rlang::caller_env()
) {
  if (!is.null(.data)) {
    rebuilt <- ipw_rebuild_design(
      numerator_model,
      stats::delete.response(stats::terms(numerator_model)),
      .data,
      call = call
    )
    check_ipw_design_width(
      rebuilt,
      numerator_model,
      "stabilize",
      component = component,
      call = call
    )

    return(rebuilt)
  }

  if (!is.null(entry$design)) {
    return(entry$design)
  }

  recovered <- tryCatch(
    stats::model.matrix(numerator_model),
    error = function(e) e
  )

  if (inherits(recovered, "error")) {
    abort_ipw_numerator_frame_gone(
      conditionMessage(recovered),
      evaluates_to = "numerator",
      component = component,
      call = call
    )
  }

  recovered
}

# The design a model's own rows enter through, which is what an integrated
# numerator's average has to be read over: that numerator is the conditional
# density averaged over the units the model was fit to, and a `.data` the caller
# supplied can leave the rest of the stack standing on fewer of them. An
# additive fit's entry has already evaluated its smooth basis over those rows,
# and every other fit is asked for its own model matrix.
#
# An `lm` or a `glm` fit with `model = FALSE` keeps no model frame and rebuilds
# one by re-evaluating its fitting call, which a fit made inside a function whose
# frame is gone cannot do. Unlike every other design here, `.data` cannot stand
# in for it, because `.data` describes the rows about to be weighted rather than
# the rows the numerator was read over.
ipw_integrated_fit_design <- function(
  mod,
  entry,
  label,
  call = rlang::caller_env()
) {
  if (!is.null(entry$design)) {
    return(entry$design)
  }

  recovered <- tryCatch(stats::model.matrix(mod), error = function(e) e)

  if (inherits(recovered, "error")) {
    abort_ipw_integrated_frame_gone(
      conditionMessage(recovered),
      label,
      call = call
    )
  }

  recovered
}

# The numerator model's prior case weights, refused for the reason the
# propensity score routes refuse their own: the block written here is the
# model's unweighted score, so a fit made under weights is not at that block's
# root and the solve would move it, reporting a numerator nobody fit. The
# refusal names `stabilize`, since that is the argument the model arrived in
# and a reader told to refit the propensity score model would be told to refit
# the wrong thing.
#
# `component` names the component the numerator was built for where there is
# one to name, which is what the joint route has and the single-treatment routes
# do not.
check_ipw_numerator_model_weights <- function(
  numerator_model,
  component = NULL,
  call = rlang::caller_env()
) {
  weights <- ipw_model_prior_weights(numerator_model)

  if (is.null(weights) || all(weights == 1)) {
    return(invisible(TRUE))
  }

  labels <- ipw_numerator_labels(component)

  abort(
    c(
      "{.fun ipw} does not support a numerator model fit with case weights.",
      x = paste0(
        labels$numerator,
        " was fit with non-unit {.arg weights}, so its coefficients are not \\
        the root of the unweighted score stacked for it."
      ),
      i = paste0("Refit ", labels$numerator, " without {.arg weights}.")
    ),
    error_class = "propensity_ipw_ps_weights_error",
    call = call
  )
}

# The mean and the score the numerator model's block contributes, read off the
# block the way the propensity score block's are read off the spec.
ipw_numerator_model_fns <- function(model) {
  ipw_continuous_score_fns(
    model$kind,
    model$link,
    model$psi_loss,
    model$psi_k,
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
#
# `component` names the component the numerator was built for where there is
# one to name, which is what the joint route has and the single-treatment routes
# do not.
check_ipw_numerator_model <- function(
  numerator_model,
  block,
  exposure_name,
  n,
  component = NULL,
  call = rlang::caller_env()
) {
  labels <- ipw_numerator_labels(component)

  response <- tryCatch(
    fmla_extract_left_chr(numerator_model),
    error = function(e) NA_character_
  )

  if (!identical(response, exposure_name)) {
    abort(
      c(
        paste0(
          "The model supplied to ",
          labels$numerator,
          " must model the exposure."
        ),
        x = paste0(
          "It models {.val {response}} and ",
          labels$model,
          " models {.val {exposure_name}}."
        ),
        i = "The numerator of the weights is what the numerator model reports
             about the exposure given what it reads, so both models describe
             the same response."
      ),
      error_class = "propensity_ipw_numerator_error",
      call = call
    )
  }

  if (!identical(nrow(block$X), as.integer(n))) {
    abort(
      c(
        paste0(
          "The model supplied to ",
          labels$numerator,
          " must be fit to the observations the other models were fit to."
        ),
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
