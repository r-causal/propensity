abort_no_method <- function(.propensity, call = rlang::caller_env()) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    call = call,
    error_class = "propensity_method_error"
  )
}

# The propensity score modifiers hold their vector route in a default method, so
# an object with no method of its own arrives there. The shapes that route reads
# are a vector of scores and, where the modifier accepts one, a data frame of
# them; anything else, a fitted model of a class the package has no method for
# among it, carries no scores to read and is reported as such rather than as a
# value out of range or a type the modifier cannot use.
check_ps_method <- function(.propensity, call = rlang::caller_env()) {
  if (!is.atomic(.propensity) && !is.data.frame(.propensity)) {
    abort_no_method(.propensity, call = call)
  }

  invisible(NULL)
}

# lifecycle decides whether a deprecation belongs to the caller or to the
# package that raised it from `user_env`, and says so: one it judges indirect
# carries a bullet naming this package and asking the reader to report an issue
# against it. A helper that signals a deprecation on another function's behalf
# therefore has to name the frame that function was called from, which is one
# further out than its own caller. That is the frame the default names, which is
# also the default lifecycle's own signalling functions carry. A route that
# reaches this helper from deeper than the function whose argument was
# deprecated names the frame it was entered on, the way it already names `call`.
handle_focal_deprecation <- function(
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  fn_name,
  user_env = rlang::caller_env(2)
) {
  # Handle deprecation warnings and parameter mapping
  if (!is.null(.treated)) {
    lifecycle::deprecate_warn(
      "0.1.0",
      paste0(fn_name, "(.treated)"),
      paste0(fn_name, "(.focal_level)"),
      user_env = user_env
    )
    if (is.null(.focal_level)) {
      .focal_level <- .treated
    }
  }

  if (!is.null(.untreated)) {
    lifecycle::deprecate_warn(
      "0.1.0",
      paste0(fn_name, "(.untreated)"),
      paste0(fn_name, "(.reference_level)"),
      user_env = user_env
    )
    if (is.null(.reference_level)) {
      .reference_level <- .untreated
    }
  }

  list(.focal_level = .focal_level, .reference_level = .reference_level)
}

# The categorical path reads one column of propensity scores per exposure level
# and never resolves a focal level, so an argument naming one describes a coding
# it does not use. Dropping such an argument silently leaves the caller believing
# the modification honored it.
check_no_focal_levels <- function(
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  call = rlang::caller_env()
) {
  supplied <- c(
    ".focal_level" = !is.null(.focal_level),
    ".reference_level" = !is.null(.reference_level),
    ".treated" = !is.null(.treated),
    ".untreated" = !is.null(.untreated)
  )

  if (!any(supplied)) {
    return(invisible(NULL))
  }

  named <- names(supplied)[supplied]
  abort(
    c(
      "{.arg {named}} {?is/are} not supported for categorical exposures.",
      x = "A categorical exposure is described by one column of scores per
           level, so no level is treated as focal and none as reference.",
      i = "Drop {.arg {named}}, or supply a binary exposure and a single column
           of scores."
    ),
    error_class = "propensity_unsupported_arg_error",
    call = call
  )
}

# The propensity score argument was named `ps` in ps_trim(), ps_trunc(),
# ps_calibrate(), and ps_tilt() before it was standardized as `.propensity`, the
# name the weight functions already read it under. The two spellings are
# resolved here: a call that supplies both is refused, the old one signals its
# deprecation, and the scores come back either way.
handle_propensity_deprecation <- function(
  .propensity,
  ps,
  fn_name,
  call = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  if (!lifecycle::is_present(ps)) {
    # Neither name was supplied, so there are no scores to modify. Dispatch
    # would report the formal it could not evaluate, which names an argument the
    # caller may have written under the other spelling.
    if (rlang::is_missing(.propensity)) {
      abort(
        c(
          "{.arg .propensity} must be supplied.",
          i = "It holds the propensity scores {.fun {fn_name}} modifies."
        ),
        error_class = "propensity_missing_arg_error",
        call = call
      )
    }

    return(.propensity)
  }

  if (!rlang::is_missing(.propensity)) {
    abort(
      c(
        "The propensity scores must be supplied once.",
        x = "They arrived under both {.arg .propensity} and {.arg ps}.",
        i = "{.arg ps} is deprecated; supply the scores as {.arg .propensity}.",
        i = "An argument given by position binds to {.arg .propensity}, so a \\
        call that names {.arg ps} must name the arguments after it as well."
      ),
      error_class = "propensity_duplicate_arg_error",
      call = call
    )
  }

  lifecycle::deprecate_warn(
    "0.2.0",
    paste0(fn_name, "(ps)"),
    paste0(fn_name, "(.propensity)"),
    user_env = user_env
  )

  ps
}

# `UseMethod()` matches the arguments of the original call against the formals
# of the method, so scores the generic read out of the deprecated `ps` do not
# reach the method under `.propensity`. Each method reads the old name for
# itself, and silently: the generic has already signaled the deprecation for
# this call and refused one that names the scores twice.
read_method_propensity <- function(.propensity, ps) {
  if (lifecycle::is_present(ps)) {
    return(ps)
  }

  .propensity
}

# A named level is one level. Both the binary coding and the categorical column
# lookup compare the exposure against the level elementwise, so a level holding
# more than one value recycles across the observations and sorts them into
# groups nobody named: on the binary route two levels alternate down the
# exposure, and on the categorical route the comparison is answered for the
# first value alone. Checked before either route reads the levels, so both give
# the same refusal.
check_focal_levels <- function(
  .focal_level,
  .reference_level,
  call = rlang::caller_env()
) {
  check_focal_level_scalar(.focal_level, ".focal_level", call = call)
  check_focal_level_scalar(.reference_level, ".reference_level", call = call)

  invisible(TRUE)
}

check_focal_level_scalar <- function(level, arg, call = rlang::caller_env()) {
  if (is.null(level) || length(level) == 1) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg {arg}} must name a single level of {.arg .exposure}.",
      x = "It holds {length(level)} value{?s}.",
      i = "Supply one level; the units on the other side of the split are
           whatever {.arg .exposure} takes besides it."
    ),
    call = call,
    error_class = "propensity_focal_level_error"
  )
}

transform_exposure_binary <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  # A matrix or data frame survives every branch below unchanged, because
  # comparison and coercion are elementwise, and only fails much later where the
  # weights are given their class. Refuse it here, where the shape is still
  # attached to the argument that carries it.
  if (!is.null(dim(.exposure))) {
    abort(
      c(
        "{.arg .exposure} must be a vector.",
        x = "It is {.cls {class(.exposure)[[1]]}} with \\
        {length(dim(.exposure))} dimension{?s}.",
        i = "Supply one value per observation."
      ),
      call = call,
      error_class = "propensity_binary_transform_error"
    )
  }

  if (!is.null(.focal_level)) {
    check_named_binary_level(
      .exposure,
      .focal_level,
      arg = ".focal_level",
      call = call
    )
    return(ifelse(.exposure == .focal_level, 1, 0))
  }

  if (!is.null(.reference_level)) {
    check_named_binary_level(
      .exposure,
      .reference_level,
      arg = ".reference_level",
      call = call
    )
    return(ifelse(.exposure != .reference_level, 1, 0))
  }

  # With no level named, an exposure that already carries the 0/1 coding, or a
  # logical one, is taken at face value and needs no announcement.
  if (is_binary(.exposure)) {
    return(.exposure)
  }

  if (is.logical(.exposure)) {
    return(as.numeric(.exposure))
  }

  if (causalgenerics::has_two_levels(.exposure)) {
    levels <- binary_exposure_levels(.exposure)
    alert_info("Setting focal level to {.val {levels[[2]]}}")
    return(ifelse(.exposure == levels[[2]], 1, 0))
  } else {
    abort(
      c(
        "Don't know how to transform `.exposure` to 0/1 binary variable.",
        i = "Specify `.focal_level` and `.reference_level`."
      ),
      call = call,
      error_class = "propensity_binary_transform_error"
    )
  }
}

# A level named as focal or as reference has to be one the exposure takes.
# Naming one it never takes sorts every unit into a single group without a word:
# a focal level nobody holds leaves every unit reference, and a reference level
# nobody holds leaves every unit focal. The categorical path already refuses
# this, and the binary path owes the same refusal whatever the exposure's
# storage. Membership is tested with the comparison the coding itself makes, so
# a level is accepted exactly when it would select the units the caller means.
check_named_binary_level <- function(
  .exposure,
  level,
  arg,
  call = rlang::caller_env()
) {
  if (any(.exposure == level, na.rm = TRUE)) {
    return(invisible(TRUE))
  }

  # An exposure that is missing everywhere takes no levels, so the hint that
  # lists the levels present has nothing to list. Report the emptiness itself,
  # which is the thing the caller has to correct.
  levels_present <- binary_exposure_levels(.exposure)
  if (length(levels_present) == 0) {
    abort(
      c(
        "{.arg {arg}} must be a level that {.arg .exposure} takes.",
        x = "{.arg .exposure} has no observed values, so it takes no levels.",
        i = "Every observation is missing. Supply an exposure with observed
             values."
      ),
      call = call,
      error_class = "propensity_focal_level_error"
    )
  }

  abort(
    c(
      "{.arg {arg}} must be a level that {.arg .exposure} takes.",
      x = "No observation takes the value {.val {level}}.",
      i = "Levels present: {.val {levels_present}}."
    ),
    call = call,
    error_class = "propensity_focal_level_error"
  )
}

# The exposure level `transform_exposure_binary()` actually codes as 1. A named
# level always wins; only when neither is named does the coding fall back to the
# exposure's own ordering, which is 1 for a 0/1 exposure and `TRUE` for a
# logical one. Mirrors that function's branching in the same order.
effective_binary_focal_level <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  if (!is.null(.focal_level)) {
    return(.focal_level)
  }

  if (!is.null(.reference_level)) {
    # Every level other than the reference is coded as focal, so a single focal
    # level is named only when the exposure takes exactly one other value.
    # Missing values are dropped first: `!=` returns NA rather than FALSE for
    # them, so an NA would survive the subscript and count as a second remaining
    # value, leaving a two-level exposure with no focal level recorded.
    observed <- .exposure[!is.na(.exposure)]
    remaining <- unique(observed[observed != .reference_level])
    if (length(remaining) != 1) {
      return(NULL)
    }
    return(as_focal_label(remaining))
  }

  if (is_binary(.exposure)) {
    return(1)
  }

  if (is.logical(.exposure)) {
    return(TRUE)
  }

  if (!causalgenerics::has_two_levels(.exposure)) {
    return(NULL)
  }

  binary_exposure_levels(.exposure)[[2]]
}

# Subsetting a factor yields a factor; record its label, which is what the
# categorical path stores and what comparisons against exposure values expect.
as_focal_label <- function(x) {
  if (is.factor(x)) as.character(x) else x
}

# Whether the named levels leave the focal level where the default coding puts
# it. Both sides come from `effective_binary_focal_level()`, so this agrees with
# `transform_exposure_binary()` by construction. A level that cannot be resolved
# counts as the default, leaving the exposure to be handled downstream as it is
# today. That covers a reference level with more than one level left beside it,
# and an exposure that does not take exactly two values. Missing values are not
# values, so an exposure that carries them is answered on the levels it takes.
resolved_focal_is_default <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  resolved <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )
  default <- effective_binary_focal_level(.exposure)

  if (is.null(resolved) || is.null(default)) {
    return(TRUE)
  }

  # `==` rather than `identical()`, matching how `transform_exposure_binary()`
  # compares the exposure against a named level: 1 and 1L are the same level.
  isTRUE(resolved == default)
}

# Binary att and atu weights are mirror images of one another: att weights built
# on one exposure level are numerically identical to atu weights built on the
# other, so nothing in the values records which level they target. Store the
# level under the attribute name the categorical path uses, but only when the
# caller named one; with no level named there is no intent to check against.
record_binary_focal_level <- function(
  psw_obj,
  .exposure,
  exposure_type,
  .focal_level = NULL,
  .reference_level = NULL
) {
  level_supplied <- !is.null(.focal_level) || !is.null(.reference_level)
  if (!identical(exposure_type, "binary") || !level_supplied) {
    return(psw_obj)
  }

  attr(psw_obj, "focal_category") <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  psw_obj
}

is_binary <- function(.exposure) {
  is.numeric(.exposure) &&
    identical(sort(unique(as.double(.exposure))), c(0, 1))
}

check_refit <- function(.propensity, call = rlang::caller_env()) {
  if (!is_refit(.propensity)) {
    warn(
      c(
        "It appears you trimmed your propensity score but did not refit the model.",
        i = "Use {.code ps_refit()} for more accurate re-estimation."
      ),
      warning_class = "propensity_no_refit_warning",
      call = call
    )
  }
}

# `call` is plumbing. Each route into a weight function, and into the
# categorical trimming and truncation methods, binds the frame it was dispatched
# into and hands it down, so that a condition raised several frames below the
# surface still names the function the user called. The generics pass their dots
# to their methods, so the argument is reachable from user code, and a value the
# condition system cannot read as a call would turn the next guard that fires
# into a report of rlang's internals. Checked where the value arrives, before
# anything is done with it.
check_call_arg <- function(call, error_call = rlang::caller_env()) {
  if (is.null(call) || rlang::is_environment(call) || rlang::is_call(call)) {
    return(invisible(call))
  }

  abort(
    c(
      "{.arg call} must be a call or an environment.",
      x = "It has class {.cls {class(call)}}.",
      i = "{.arg call} names the frame a condition is attributed to. Leave it
           unset unless you are wrapping a function from this package."
    ),
    error_class = "propensity_call_arg_error",
    call = error_call
  )
}

# The refusal names `.propensity`, which is what every function that reads
# propensity scores calls the argument they arrive in.
check_ps_range <- function(ps, call = rlang::caller_env()) {
  if (is.matrix(ps) || is.data.frame(ps)) {
    # For matrices/data frames, check all values
    ps_vals <- as.numeric(as.matrix(ps))
    # Check only non-NA values
    non_na_vals <- ps_vals[!is.na(ps_vals)]
    if (
      length(non_na_vals) > 0 &&
        any(non_na_vals <= 0 | non_na_vals >= 1 | !is.finite(non_na_vals))
    ) {
      abort(
        c(
          "All propensity scores must be between 0 and 1.",
          i = "The range of values in {.arg .propensity} is \\
        {format(range(ps_vals, na.rm = TRUE), nsmall = 1, digits = 1)}"
        ),
        call = call,
        error_class = "propensity_range_error"
      )
    }
  } else {
    ps <- as.numeric(ps)
    # Check only non-NA values
    non_na_vals <- ps[!is.na(ps)]
    if (
      length(non_na_vals) > 0 &&
        any(non_na_vals <= 0 | non_na_vals >= 1 | !is.finite(non_na_vals))
    ) {
      abort(
        c(
          "The propensity score must be between 0 and 1.",
          i = "The range of {.arg .propensity} is \\
        {format(range(ps, na.rm = TRUE), nsmall = 1, digits = 1)}"
        ),
        call = call,
        error_class = "propensity_range_error"
      )
    }
  }

  invisible(TRUE)
}

# The default for `stabilize` is read from the exposure type rather than fixed
# once. A continuous exposure is weighted by a density ratio whose denominator
# carries the exposure's own units, so only the stabilized form of those weights
# is recommended, and it is what an unnamed `stabilize` asks for. A binary or
# categorical exposure is weighted by a probability ratio that needs no such
# numerator, so it stays unstabilized, as it always has been. An explicit `TRUE`
# or `FALSE` is passed through untouched.
resolve_stabilize <- function(stabilize, exposure_type) {
  if (!is.null(stabilize)) {
    return(stabilize)
  }

  identical(exposure_type, "continuous")
}

# The conditional spread of a continuous exposure, checked where the exposure
# type it belongs to is known. `.sigma` is read straight into a normal density,
# so a value that is not a number reaches the density as though it were one, and
# a length that neither matches nor divides the exposure recycles into weights
# nobody asked for. `.sigma` also sits in the third position, which is where an
# exposure type supplied without a name arrives; a binary or categorical
# exposure has no conditional spread to set, so a score that reaches one is
# reported rather than absorbed.
check_sigma <- function(
  .sigma,
  exposure_type,
  n,
  call = rlang::caller_env()
) {
  if (is.null(.sigma)) {
    return(invisible(NULL))
  }

  if (!is.numeric(.sigma)) {
    abort(
      c(
        "{.arg .sigma} must be numeric or {.code NULL}.",
        x = "It has class {.cls {class(.sigma)}}.",
        i = "{.arg .sigma} is the third argument, so a value meant for
             {.arg exposure_type} arrives here when it is supplied by position.
             Name the argument you mean."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  # A standard deviation scales the density one observation at a time, so a
  # shape the exposure does not have is refused rather than flattened.
  if (!is.null(dim(.sigma))) {
    shape <- paste(dim(.sigma), collapse = " x ")
    abort(
      c(
        "{.arg .sigma} must not have dimensions.",
        x = "It is {shape}.",
        i = "Supply a single standard deviation, or one for each observation."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  if (length(.sigma) != 1 && length(.sigma) != n) {
    abort(
      c(
        "{.arg .sigma} must hold one value or one value per observation.",
        x = "It holds {length(.sigma)} value{?s}.",
        x = "{.arg .exposure} has {n} observation{?s}.",
        i = "Supply a single standard deviation to spread every conditional
             density alike, or one for each observation."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  if (!identical(exposure_type, "continuous")) {
    abort(
      c(
        "{.arg .sigma} applies only to continuous exposures.",
        x = "{.arg .exposure} is being treated as {exposure_type}.",
        i = "{.arg .sigma} spreads the conditional density of a continuous
             exposure. Leave it unset for a {exposure_type} exposure."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  # The density reads the spread as it arrives, so a value that cannot scale one
  # is reported rather than evaluated: a negative spread comes back as `NaN`
  # under a base R warning about them, a zero comes back infinite, and an
  # unbounded one comes back as no density at all. A spread that is missing for
  # one observation leaves that observation's weight missing, the way a missing
  # score does; one missing for every observation describes no density anywhere.
  absent <- is.na(.sigma) & !is.nan(.sigma)
  if (all(absent)) {
    abort(
      c(
        "{.arg .sigma} must hold at least one standard deviation.",
        x = "Every value is missing.",
        i = "{.arg .sigma} spreads the conditional density of the exposure,
             which needs a width to be evaluated at."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  supplied <- .sigma[!absent]
  unusable <- supplied[!is.finite(supplied) | supplied <= 0]
  if (length(unusable) > 0) {
    abort(
      c(
        "{.arg .sigma} must hold positive, finite standard deviations.",
        x = "It holds {.val {unique(unusable)}}.",
        i = "{.arg .sigma} spreads the conditional density of the exposure, and
             a width at or below zero, or without a bound, spreads nothing."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  invisible(.sigma)
}

# A matrix of one column per exposure level holds what the equivalent data frame
# holds, so a binary exposure reads it the same way: the column carrying the
# probability of the resolved focal level. A matrix left as it is survives every
# check elementwise, because comparison and coercion are elementwise, and only
# fails where the weights are given their class, several frames below the
# argument that carries the shape.
resolve_binary_propensity <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  if (!is.matrix(.propensity)) {
    return(.propensity)
  }

  extract_propensity_from_df(
    as.data.frame(.propensity),
    .propensity_col_quo = rlang::quo(NULL),
    .exposure = .exposure,
    exposure_type = "binary",
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )
}

# `lower` and `upper` are read into comparisons that decide which units to keep
# or pin, and a missing bound answers neither `TRUE` nor `FALSE`, so the first
# comparison it reaches raises a bare error about a missing value and names
# nothing the caller wrote.
check_bounds_not_missing <- function(
  lower,
  upper,
  call = rlang::caller_env()
) {
  missing_bounds <- c("lower", "upper")[c(anyNA(lower), anyNA(upper))]

  if (length(missing_bounds) == 0) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg {missing_bounds}} must not be missing.",
      i = "Each bound is read into the comparison that decides what happens to a
           propensity score, and a missing bound decides nothing.",
      i = "Supply a value, or leave the argument unset to take the default for
           this method."
    ),
    call = call,
    error_class = "propensity_missing_value_error"
  )
}

check_lower_upper <- function(lower, upper, call = rlang::caller_env()) {
  check_bounds_not_missing(lower, upper, call = call)

  if (lower >= upper) {
    abort(
      c(
        "{.arg lower} must be smaller than {.arg upper}",
        x = "{.arg lower} is {lower} and {.arg upper} is {upper}"
      ),
      call = call,
      error_class = "propensity_range_error"
    )
  }

  invisible(TRUE)
}

# The percentile methods read their cutoffs off `quantile()`, whose `probs`
# argument is not one `ps_trim()` or `ps_trunc()` takes. A bound outside the unit
# interval is refused here so the caller is told about the argument they wrote
# rather than about one they cannot reach.
check_quantile_probs <- function(lower, upper, call = rlang::caller_env()) {
  check_bounds_not_missing(lower, upper, call = call)

  in_unit_interval <- function(x) x >= 0 && x <= 1

  if (in_unit_interval(lower) && in_unit_interval(upper)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "For {.code method = 'pctl'}, {.arg lower} and {.arg upper} must be
       between 0 and 1.",
      x = "{.arg lower} is {lower} and {.arg upper} is {upper}.",
      i = "The bounds are quantile probabilities rather than propensity scores,
           so each must lie in [0, 1]."
    ),
    call = call,
    error_class = "propensity_range_error"
  )
}

# A cutoff read off the exposure groups is undefined for a unit that belongs to
# neither, and the rule then comes out missing rather than absent: every score
# compares as outside a missing bound, or as inside one, so the whole sample is
# trimmed or the recorded bounds are never applied. The missing exposure is
# reported instead.
check_exposure_complete <- function(
  .exposure,
  method,
  call = rlang::caller_env()
) {
  n_missing <- sum(is.na(.exposure))

  if (n_missing == 0) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg .exposure} must not have missing values for method {.val {method}}.",
      x = "{n_missing} exposure value{?s} {?is/are} missing.",
      i = "The cutoffs are read off the exposure groups, and a unit that belongs
           to neither leaves them undefined.",
      i = "Remove or impute the missing exposure values."
    ),
    call = call,
    error_class = "propensity_missing_value_error"
  )
}

# The common range runs from the lowest score among the focal units to the
# highest among the reference units, so each bound needs at least one observed
# score in its own group. Over none of them `min()` returns `Inf` and `max()`
# returns `-Inf`, both under a base R warning, and the bounds then cross: the
# group with no scores to read a bound off is reported as distributions that do
# not overlap, or leaves a range no unit falls inside.
check_cr_groups_observed <- function(
  focal_scores,
  reference_scores,
  call = rlang::caller_env()
) {
  empty <- c(
    "focal" = length(focal_scores) == 0,
    "reference" = length(reference_scores) == 0
  )

  if (!any(empty)) {
    return(invisible(TRUE))
  }

  groups <- names(empty)[empty]
  abort(
    c(
      "For {.code method = 'cr'}, every exposure group must hold at least one
       observed propensity score.",
      x = "No score is observed among the {groups} units.",
      i = "The bounds are the lowest score among the focal units and the highest
           among the reference units, and a group with no observed score has
           neither to give.",
      i = "Remove or impute the missing propensity scores, or bound the scores
           with {.code method = 'ps'} or {.code method = 'pctl'}."
    ),
    call = call,
    error_class = "propensity_no_data_error"
  )
}

# The shape a continuous exposure's `.propensity` has to have. The weights are
# a ratio of densities evaluated at one conditional mean for each unit, and the
# arithmetic that evaluates them vectorizes over whatever it is handed, so a
# matrix of means yields a weight for every cell rather than a weight for every
# unit. Anything carrying a dimension is refused here, rather than left to
# produce a result of the wrong length.
check_continuous_ps_shape <- function(
  .propensity,
  call = rlang::caller_env()
) {
  dims <- dim(.propensity)

  if (is.null(dims)) {
    return(invisible(TRUE))
  }

  shape <- paste(dims, collapse = " by ")

  abort(
    c(
      "Weights for a continuous exposure need one conditional mean for each
       unit.",
      x = "{.arg .propensity} is {shape} and of class
           {.cls {class(.propensity)[[1]]}}.",
      i = "Pass a numeric vector of conditional means, such as the single
           column of {.arg .propensity} that holds the mean of this exposure."
    ),
    call = call,
    error_class = "propensity_ps_shape_error"
  )
}

check_lengths_match <- function(
  .propensity,
  .exposure,
  call = rlang::caller_env()
) {
  # Handle matrix/data.frame inputs
  if (is.matrix(.propensity) || is.data.frame(.propensity)) {
    len_prop <- nrow(.propensity)
  } else {
    len_prop <- length(.propensity)
  }

  len_exp <- length(.exposure)

  if (len_prop != len_exp) {
    abort(
      c(
        "{.arg .propensity} and {.arg .exposure} must have the same length.",
        i = "{.arg .propensity} has {if (is.matrix(.propensity) || is.data.frame(.propensity)) 'rows' else 'length'} {len_prop}",
        i = "{.arg .exposure} has length {len_exp}"
      ),
      call = call,
      error_class = "propensity_length_error"
    )
  }

  invisible(TRUE)
}

transform_exposure_categorical <- function(
  .exposure,
  .focal_level = NULL,
  call = rlang::caller_env()
) {
  # Convert to factor if not already
  if (!is.factor(.exposure)) {
    .exposure <- as.factor(.exposure)
  }

  # Check if we have more than 2 levels
  n_levels <- nlevels(.exposure)
  if (n_levels <= 2) {
    abort(
      c(
        "Categorical exposure must have more than 2 levels.",
        i = "Found {n_levels} levels.",
        i = "Use binary exposure methods for 2-level exposures."
      ),
      call = call,
      error_class = "propensity_categorical_levels_error"
    )
  }

  # Validate focal category if provided
  if (!is.null(.focal_level)) {
    if (!.focal_level %in% levels(.exposure)) {
      abort(
        c(
          "Focal category must be one of the exposure levels.",
          i = "Focal category: {.val {(.focal_level)}}",
          i = "Available levels: {.val {levels(.exposure)}}"
        ),
        call = call,
        error_class = "propensity_focal_category_error"
      )
    }
  }

  .exposure
}

check_ps_matrix <- function(
  ps_matrix,
  .exposure,
  call = rlang::caller_env()
) {
  # Convert to matrix if data frame first
  if (is.data.frame(ps_matrix)) {
    ps_matrix <- as.matrix(ps_matrix)
  }

  # Check if it's a matrix
  if (!is.matrix(ps_matrix)) {
    abort(
      "For categorical exposures, `.propensity` must be a matrix or data frame.",
      call = call,
      error_class = "propensity_matrix_type_error"
    )
  }

  # Check dimensions
  n_obs <- length(.exposure)
  n_cats <- nlevels(.exposure)

  if (nrow(ps_matrix) != n_obs) {
    abort(
      c(
        "Number of rows in propensity score matrix must match number of observations.",
        i = "Matrix rows: {nrow(ps_matrix)}",
        i = "Observations: {n_obs}"
      ),
      call = call,
      error_class = "propensity_matrix_dims_error"
    )
  }

  if (ncol(ps_matrix) != n_cats) {
    abort(
      c(
        "Number of columns in propensity score matrix must match number of exposure categories.",
        i = "Matrix columns: {ncol(ps_matrix)}",
        i = "Categories: {n_cats}"
      ),
      call = call,
      error_class = "propensity_matrix_dims_error"
    )
  }

  check_ps_matrix_rowsums(ps_matrix, call = call)
  check_ps_matrix_range(ps_matrix, call = call)

  # Ensure columns are in the same order as factor levels
  # This is critical for correct weight calculation
  exp_levels <- levels(.exposure)

  # Check if columns have names
  if (!is.null(colnames(ps_matrix))) {
    # Try to match column names to factor levels
    # Handle both plain names (A, B, C) and parsnip-style names (.pred_A, .pred_B, .pred_C)
    col_names <- colnames(ps_matrix)

    # Remove common prefixes like ".pred_"
    clean_names <- gsub("^\\.pred_", "", col_names)

    # Check if clean names match factor levels
    if (setequal(clean_names, exp_levels)) {
      # Reorder columns to match factor levels
      if (!identical(clean_names, exp_levels)) {
        col_order <- match(exp_levels, clean_names)
        ps_matrix <- ps_matrix[, col_order, drop = FALSE]
        # Update column names to match
        colnames(ps_matrix) <- col_names[col_order]
      }
    } else {
      # Column names don't match factor levels
      abort(
        c(
          "Column names of propensity score matrix must match exposure levels.",
          i = "Expected levels: {.val {exp_levels}}",
          i = "Found columns: {.val {clean_names}}"
        ),
        call = call,
        error_class = "propensity_matrix_names_error"
      )
    }
  } else {
    # No column names - assume they're in factor level order
    # Issue a warning as this is risky
    warn(
      c(
        "Propensity score matrix has no column names.",
        i = "Assuming columns are in factor level order: {.val {exp_levels}}",
        i = "This may lead to incorrect results if columns are misaligned."
      ),
      warning_class = "propensity_matrix_no_names_warning",
      call = call
    )
  }

  ps_matrix
}

# Each row of a categorical propensity score matrix is a probability vector.
# Shared with ps_tilt(), which reads a matrix with no exposure beside it and so
# cannot run the rest of check_ps_matrix().
check_ps_matrix_rowsums <- function(ps_matrix, call = rlang::caller_env()) {
  # Only check non-NA rows
  row_sums <- rowSums(ps_matrix, na.rm = FALSE)
  ROW_SUM_TOLERANCE <- 1e-6 # Tolerance for floating point comparison
  non_na_rows <- !is.na(row_sums)

  if (any(non_na_rows)) {
    # Check only the rows that don't have NA values
    if (any(abs(row_sums[non_na_rows] - 1) > ROW_SUM_TOLERANCE)) {
      bad_rows <- which(abs(row_sums - 1) > ROW_SUM_TOLERANCE & non_na_rows)
      abort(
        c(
          "Propensity score matrix rows must sum to 1.",
          i = "Problem rows: {bad_rows[1:min(5, length(bad_rows))]}{if (length(bad_rows) > 5) ' ...' else ''}"
        ),
        call = call,
        error_class = "propensity_matrix_sum_error"
      )
    }
  }

  invisible(TRUE)
}

# The bounds are the same open interval that check_ps_range() enforces on the
# binary path. A score of exactly 0 for a unit's observed level divides by zero
# in the weight calculation, so the endpoints are rejected rather than carried
# through as infinite weights.
#
# The refusal names `.propensity`, the argument every function that reads
# propensity scores takes them in; reporting a bare range says nothing about
# what the caller has to correct.
check_ps_matrix_range <- function(ps_matrix, call = rlang::caller_env()) {
  # `min()` and `max()` compare strings rather than numbers, so a matrix that
  # is not made of numbers has to be refused before the bounds are read.
  if (!is.numeric(ps_matrix) && !is.logical(ps_matrix)) {
    abort(
      "{.arg .propensity} must be numeric.",
      call = call,
      error_class = "propensity_matrix_type_error"
    )
  }

  # Missing values are not scores, so they neither decide the bounds nor are
  # refused. `min()` and `max()` skip them without copying the matrix, but they
  # have no bounds to report when there is nothing left to read, so that case is
  # answered first.
  if (length(ps_matrix) == 0 || (anyNA(ps_matrix) && all(is.na(ps_matrix)))) {
    return(invisible(TRUE))
  }

  lower <- min(ps_matrix, na.rm = TRUE)
  upper <- max(ps_matrix, na.rm = TRUE)

  # An infinite score is out of bounds on the side it is infinite on, so the
  # bounds alone catch it.
  if (lower <= 0 || upper >= 1) {
    abort(
      c(
        "All propensity scores must be between 0 and 1.",
        i = "The bounds are exclusive: a score of exactly 0 or 1 leaves the \\
        weight undefined.",
        i = "The range of values in {.arg .propensity} is \\
        {format(as.numeric(c(lower, upper)), nsmall = 1, digits = 1)}"
      ),
      call = call,
      error_class = "propensity_range_error"
    )
  }

  invisible(TRUE)
}

# Helper for ps_trim and ps_trunc methods
calculate_weight_from_modified_ps <- function(
  .propensity,
  .exposure,
  weight_fn,
  modification_type = c("trim", "trunc", "calib"),
  ...,
  call = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  # `weight_fn` is a numeric method invoked through a local symbol, so its own
  # frame deparses as `weight_fn()`. The default is the dispatching method's
  # frame, which is handed down so that conditions raised inside the weight
  # machinery name the weight function the user called. `user_env` is handed
  # down for the same reason: the numeric method is reached from here rather
  # than from the caller who wrote the deprecated argument, so its own caller
  # would attribute the deprecation to this package.
  check_call_arg(call, error_call = rlang::caller_env())

  modification_type <- rlang::arg_match(modification_type)

  # Only check refit for trim
  if (modification_type == "trim") {
    check_refit(.propensity, call = call)
  }

  # Handle matrix or vector propensity scores
  if (inherits(.propensity, c("ps_trim_matrix", "ps_trunc_matrix"))) {
    # For matrix propensity scores, pass them directly
    # The weight function should handle the matrix appropriately
    base_wt <- weight_fn(
      .propensity,
      .exposure = .exposure,
      call = call,
      user_env = user_env,
      ...
    )
  } else {
    # Convert to numeric for vector propensity scores
    numeric_ps <- as.numeric(.propensity)

    # Call the weight function with the numeric propensity scores
    base_wt <- weight_fn(
      numeric_ps,
      .exposure = .exposure,
      call = call,
      user_env = user_env,
      ...
    )
  }

  # Update estimand
  if (modification_type == "trim") {
    old_est <- estimand(base_wt)
    estimand(base_wt) <- paste0(old_est, "; trimmed")
    attr(base_wt, "trimmed") <- TRUE
    attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")
  } else if (modification_type == "trunc") {
    estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")
    attr(base_wt, "truncated") <- TRUE
    attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)
  } else {
    # calib
    estimand(base_wt) <- paste0(estimand(base_wt), "; calibrated")
    attr(base_wt, "calibrated") <- TRUE
    attr(base_wt, "ps_calib_meta") <- ps_calib_meta(.propensity)
  }

  base_wt
}

# What the weights record about the exposure they were built for. The type is
# the one the weight function resolved, and the density record is read off the
# weights the formula returned: only the continuous ATE formula writes one, and
# only there is a weight a ratio of densities.
record_exposure_attrs <- function(psw_obj, wts, exposure_type) {
  attr(psw_obj, "exposure_type") <- exposure_type
  attr(psw_obj, "density_meta") <- attr(wts, "density_meta")

  psw_obj
}

# Helper to preserve categorical attributes on psw objects
preserve_categorical_attrs <- function(psw_obj, wts, exposure_type) {
  if (exposure_type == "categorical") {
    attr(psw_obj, "n_categories") <- attr(wts, "n_categories")
    attr(psw_obj, "category_names") <- attr(wts, "category_names")
    # focal_category might not always exist
    if (!is.null(attr(wts, "focal_category"))) {
      attr(psw_obj, "focal_category") <- attr(wts, "focal_category")
    }
  }
  psw_obj
}

# The exposure levels the binary coding resolves against, in the order it
# resolves them: the levels a factor takes, or the sorted values anything else
# takes. A factor answers for every level it declares, whether or not any
# observation holds it, so a declared level nobody holds would otherwise sit
# between the two real ones and be coded as focal, returning weights for an
# exposure nobody has. `sort()` and `droplevels()` both leave missing values
# out, which is what keeps an NA from counting as a level of its own.
binary_exposure_levels <- function(.exposure) {
  if (is.factor(.exposure)) {
    levels(droplevels(.exposure))
  } else {
    sort(unique(.exposure))
  }
}

# Position of the column holding the probability of the resolved focal level,
# or NULL when the column names cannot answer the question. The names must
# cover every exposure level before any of them is trusted, so that a lone
# coincidental match cannot redirect the selection. A reference level that
# leaves more than one candidate resolves to no focal level at all, which is
# likewise unanswerable.
match_focal_level_column <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  exposure_levels <- as.character(binary_exposure_levels(.exposure))
  if (!all(exposure_levels %in% names(.propensity))) {
    return(NULL)
  }

  focal_level <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )
  if (is.null(focal_level)) {
    return(NULL)
  }

  col_pos <- match(as.character(focal_level), names(.propensity))
  if (is.na(col_pos)) {
    return(NULL)
  }

  warn_ambiguous_focal_column(
    .propensity,
    focal_level = focal_level,
    col_pos = col_pos,
    call = call
  )

  col_pos
}

# `match()` takes the first column of a given name, and a data frame is allowed
# more than one column of the same name. The selection is then made on nothing
# the caller expressed, between columns that may hold different numbers, and the
# one it lands on is as likely to be the wrong one. That the columns might also
# be copies of each other is not worth checking for: the ambiguity is in the
# names, and a caller who meant one of two identically named columns is owed the
# same report as one whose columns differ. Only the selected name is at issue: a
# name repeated elsewhere in the frame leaves this column unambiguous.
warn_ambiguous_focal_column <- function(
  .propensity,
  focal_level,
  col_pos,
  call = rlang::caller_env()
) {
  focal_name <- as.character(focal_level)
  # A data frame is allowed a missing column name, which compares to the focal
  # name as `NA` rather than as no match. Left in, that `NA` reaches the
  # comparison below and stops the call on a frame this function has nothing to
  # say about.
  n_matching <- sum(names(.propensity) == focal_name, na.rm = TRUE)
  if (n_matching < 2) {
    return(invisible(FALSE))
  }

  warn(
    c(
      "{.arg .propensity} has {n_matching} columns named {.val {focal_name}}, \\
      the level resolved as focal.",
      i = "Read column {col_pos}, the first of them.",
      i = "Give the columns distinct names, or set {.arg .propensity_col}, to \\
      select the column you mean."
    ),
    warning_class = "propensity_df_duplicate_column_warning",
    call = call
  )

  invisible(TRUE)
}

# Helper function to extract propensity scores from data frames
# This consolidates the logic used across multiple weight functions
extract_propensity_from_df <- function(
  .propensity,
  .propensity_col_quo = NULL,
  .exposure = NULL,
  exposure_type = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  if (!rlang::quo_is_null(.propensity_col_quo)) {
    col_pos <- tryCatch(
      tidyselect::eval_select(
        .propensity_col_quo,
        data = .propensity
      ),
      error = function(e) {
        abort(
          paste0("Column selection failed: ", e$message),
          call = call,
          error_class = "propensity_df_column_error"
        )
      }
    )

    if (length(col_pos) != 1) {
      abort(
        "`.propensity_col` must select exactly one column.",
        call = call,
        error_class = "propensity_df_column_error"
      )
    }

    return(.propensity[[col_pos]])
  }

  if (ncol(.propensity) == 0) {
    abort(
      "`.propensity` data frame must have at least one column.",
      call = call,
      error_class = "propensity_df_ncol_error"
    )
  }

  is_binary_exposure <- identical(exposure_type, "binary")

  if (is_binary_exposure) {
    col_pos <- match_focal_level_column(
      .propensity,
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    if (!is.null(col_pos)) {
      return(.propensity[[col_pos]])
    }
  }

  # Default behavior: use second column if available, otherwise first
  col_pos <- if (ncol(.propensity) >= 2) 2L else 1L

  # A single column is the caller supplying the probability of the focal level
  # directly, so there is nothing to choose between and nothing to report.
  level_named <- !is.null(.focal_level) || !is.null(.reference_level)
  if (is_binary_exposure && level_named && ncol(.propensity) > 1) {
    warn(
      c(
        "Can't tell which column of {.arg .propensity} holds the probability of the focal level.",
        i = "Selected {.val {names(.propensity)[[col_pos]]}} by position.",
        i = "Name the columns after the levels of {.arg .exposure}, or set {.arg .propensity_col}, to select the column by name."
      ),
      warning_class = "propensity_df_column_warning",
      call = call
    )
  }

  .propensity[[col_pos]]
}

# Where a fitted model keeps the probability of a binary exposure. Most classes
# report it through `predict(type = "response")`; a class whose `predict()`
# method offers no such type, such as `nnet::multinom()`, has a method of its
# own. The family check runs first, so a method here is only ever handed a
# model that fits probabilities.
extract_binary_ps <- function(model, call = rlang::caller_env()) {
  UseMethod("extract_binary_ps")
}

#' @export
extract_binary_ps.default <- function(model, call = rlang::caller_env()) {
  stats::predict(model, type = "response")
}

# The propensity score a fitted model reports, read according to the exposure it
# is a model of. A binary exposure needs the probability of the exposure, which
# only the binomial families fit; a continuous exposure needs a conditional mean
# with a single spread, which `extract_continuous_ps()` reads from the classes
# that fit one; a categorical exposure needs a probability for every level,
# which `extract_categorical_ps()` reads from the classes that fit one per
# level.
#
# `exposure_levels` carries the levels of the exposure being weighted, so that a
# categorical model fit to some other set of them is reported against the model
# rather than against the shape of what it returns.
extract_model_propensity <- function(
  model,
  exposure_type,
  exposure_levels = NULL,
  call = rlang::caller_env()
) {
  if (identical(exposure_type, "binary")) {
    check_binary_model_family(model, call = call)

    return(extract_binary_ps(model, call = call))
  }

  if (identical(exposure_type, "continuous")) {
    return(extract_continuous_ps(model, call = call)$mu)
  }

  extract_categorical_ps(model, exposure_levels, call = call)
}

# Helper function to handle common data frame method pattern
# This encapsulates the logic used across all weight function data.frame methods
#
# The deprecated arguments are resolved here rather than downstream, for the
# same reason as in `prepare_model_weight_args()`: the resolved focal level
# decides which column carries `.propensity`, and mapping it twice would emit
# the deprecation warning twice, so the numeric method receives the mapped
# levels and no deprecated arguments.
handle_data_frame_weight_calculation <- function(
  weight_fn_numeric,
  .propensity,
  .exposure,
  exposure_type,
  valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
  .propensity_col_quo,
  .focal_level = NULL,
  .reference_level = NULL,
  .treated = NULL,
  .untreated = NULL,
  fn_name,
  ...,
  call = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  # The default is the frame of the data frame method this was called from,
  # which is the frame the generic dispatched into and so the one that names the
  # weight function the user called. A caller may name another frame, and the
  # generics' dots deliver it here rather than to the numeric method, so that
  # the whole route is attributed to the same place. A value that cannot be
  # read as a call is reported against the method rather than against itself.
  check_call_arg(call, error_call = rlang::caller_env())

  # Validate inputs
  if (!is.data.frame(.propensity)) {
    abort(
      "`.propensity` must be a data frame.",
      call = call,
      error_class = "propensity_matrix_type_error"
    )
  }

  # Check exposure type
  exposure_type_check <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    valid_exposure_types,
    announce = !be_quiet(),
    call = call
  )

  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    fn_name,
    user_env = user_env
  )
  check_focal_levels(
    focal_params$.focal_level,
    focal_params$.reference_level,
    call = call
  )

  # The resolved type is what the numeric method is handed, not the argument it
  # was resolved from. Resolving it twice makes the same decision twice and
  # announces it twice, once for a call that made one decision.
  if (exposure_type_check == "categorical") {
    # For categorical exposures, pass the whole data frame
    return(weight_fn_numeric(
      .propensity = .propensity,
      .exposure = .exposure,
      exposure_type = exposure_type_check,
      .focal_level = focal_params$.focal_level,
      .reference_level = focal_params$.reference_level,
      call = call,
      ...
    ))
  }

  # For non-categorical exposures, extract single column
  ps_vec <- extract_propensity_from_df(
    .propensity,
    .propensity_col_quo,
    .exposure = .exposure,
    exposure_type = exposure_type_check,
    .focal_level = focal_params$.focal_level,
    .reference_level = focal_params$.reference_level,
    call = call
  )

  # The numeric method is reached by value here, so its own frame deparses as
  # `weight_fn_numeric()`. Hand it the call the route was entered on.
  weight_fn_numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type_check,
    .focal_level = focal_params$.focal_level,
    .reference_level = focal_params$.reference_level,
    call = call,
    ...
  )
}

# Helper functions for extracting information from GLM formulas
# (moved from ipw.R to be shared across the package)
fmla_extract_left_vctr <- function(mod) {
  .data <- mod |>
    model.frame()

  .data[[1]]
}

fmla_extract_left_chr <- function(mod) {
  as.character(formula(mod)[[2]])
}

# Variables a model's left-hand side reads. `fmla_extract_left_chr()` deparses
# the whole expression, so a transformed response such as `log(y)` gives back the
# function name alongside its argument; a rebuild from a supplied data frame needs
# the columns the response is computed from, which is what `all.vars()` gives and
# what the user can actually supply.
fmla_extract_left_vars <- function(mod) {
  all.vars(formula(mod)[[2]])
}

# Evaluate a model's left-hand side against a data frame. A response written as a
# transformation has to be computed rather than looked up by name, and the
# functions and constants it reads resolve from the formula's environment, the
# same way they did when the model was fit.
#
# The formula's environment is the enclosure, so a column absent from `data` is
# taken from there instead of reported. Every caller asserts the columns the
# response reads are present before reaching this, which is what keeps a
# `.data` missing the response from being answered out of the enclosure. That
# ordering is the whole protection and is pinned as such.
#
# `scale()` and friends return a one-column matrix; drop that dimension so the
# result is an ordinary vector. The model-frame route reads its response through
# `fmla_extract_left_vctr()` and keeps the matrix, so the two shapes are not the
# same here; they converge at `ipw_outcome_numeric()`, whose `as.double()`
# flattens either.
fmla_eval_left <- function(mod, data) {
  fmla <- formula(mod)
  left <- eval(fmla[[2]], data, environment(fmla))

  if (is.matrix(left) && ncol(left) == 1L) {
    left <- drop(left)
  }

  left
}

# The exposure a fitted model was fit to. Most classes keep it in their model
# frame, which is what the default method reads; a class that has to rebuild its
# frame from the fitting call, and so cannot when that call is no longer
# reachable, has a method of its own.
extract_model_exposure <- function(model) {
  UseMethod("extract_model_exposure")
}

#' @export
extract_model_exposure.default <- function(model) {
  fmla_extract_left_vctr(model)
}

# Helper function to handle optional exposure in the model methods
extract_exposure_from_model <- function(
  model,
  .exposure = NULL
) {
  if (is.null(.exposure)) {
    # Extract the exposure from the model's response
    .exposure <- extract_model_exposure(model)
    exposure_name <- fmla_extract_left_chr(model)
    alert_info(
      "Using exposure variable {.val {exposure_name}} from the propensity score model"
    )
  }
  .exposure
}

# Shared preparation for the model methods of the weight functions.
#
# Fitted values report the probability of the level the response's default
# coding treats as focal: the second factor level, or the second sorted value.
# The weight formulas read `.propensity` as the probability of the level
# actually resolved as focal, so naming the other level means the fitted values
# must be inverted before the numeric method sees them.
#
# The deprecated arguments are resolved here rather than downstream. The
# resolved focal level is needed to decide whether to invert, and mapping it
# twice would emit the deprecation warning twice, so the numeric method
# receives the mapped levels and no deprecated arguments.
prepare_model_weight_args <- function(
  model,
  .exposure,
  exposure_type,
  valid_exposure_types,
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  fn_name,
  call = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  .exposure <- extract_exposure_from_model(model, .exposure)
  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    valid_exposure_types,
    announce = !be_quiet(),
    call = call
  )

  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    fn_name,
    user_env = user_env
  )
  check_focal_levels(
    focal_params$.focal_level,
    focal_params$.reference_level,
    call = call
  )

  # The levels a categorical model has to report a probability for are the
  # levels of the exposure it is weighting, read the way the categorical path
  # reads them. The other exposure types have no levels to compare.
  exposure_levels <- if (identical(exposure_type, "categorical")) {
    levels(as.factor(.exposure))
  }

  ps_vec <- extract_model_propensity(
    model,
    exposure_type,
    exposure_levels = exposure_levels,
    call = call
  )

  invert <- identical(exposure_type, "binary") &&
    !resolved_focal_is_default(
      .exposure,
      .focal_level = focal_params$.focal_level,
      .reference_level = focal_params$.reference_level
    )

  if (invert) {
    ps_vec <- 1 - ps_vec
  }

  list(
    propensity = ps_vec,
    exposure = .exposure,
    exposure_type = exposure_type,
    focal_level = focal_params$.focal_level,
    reference_level = focal_params$.reference_level
  )
}

# Shared preparation for the model methods of the propensity score modifiers.
#
# A modifier reads a fitted model the way the weight functions read one, and
# the extraction is the same: the scores come off the fit, and so does the
# exposure when the modification needs one. What differs is that a modifier has
# no exposure type to resolve. Which route a fit takes is a property of the fit
# itself, read by `model_fits_levels()`, so a call that needs no exposure never
# reads one off the model and is not told about an exposure it makes no use of.
#
# `needs_exposure` is what the modification itself requires: trimming to the
# preference scale or to a common range reads the exposure, trimming a matrix of
# scores reads it to hold the columns against the levels, and calibration always
# does, while a tilt and a trimming to a fixed threshold read none. A named focal
# level requires it too, because the fitted values report the probability of the
# level the response's default coding treats as focal, so naming the other level
# means inverting them, exactly as `prepare_model_weight_args()` does.
prepare_model_ps <- function(
  model,
  .exposure = NULL,
  needs_exposure = FALSE,
  .focal_level = NULL,
  .reference_level = NULL,
  .treated = NULL,
  .untreated = NULL,
  fn_name,
  call = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  # The deprecated arguments are resolved here rather than downstream, for the
  # same reason as in `prepare_model_weight_args()`: the resolved focal level
  # decides whether the fitted values are inverted, and mapping it twice would
  # emit the deprecation warning twice, so the route below receives the mapped
  # levels and no deprecated arguments.
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    fn_name,
    user_env = user_env
  )
  check_focal_levels(
    focal_params$.focal_level,
    focal_params$.reference_level,
    call = call
  )

  if (model_fits_levels(model)) {
    # The levels a categorical model has to report a probability for are the
    # levels of the exposure being modified, read the way the categorical path
    # reads them. A modification that reads no exposure has nothing to hold the
    # columns against, and recovering one anyway would announce a reading that
    # never happened, so the levels are taken off the fit, which records them.
    if (needs_exposure || !is.null(.exposure)) {
      .exposure <- extract_exposure_from_model(model, .exposure)
      exposure_levels <- levels(as.factor(.exposure))
    } else {
      exposure_levels <- model_levels(model)
    }

    return(list(
      propensity = extract_model_propensity(
        model,
        "categorical",
        exposure_levels = exposure_levels,
        call = call
      ),
      exposure = .exposure,
      exposure_type = "categorical",
      focal_level = focal_params$.focal_level,
      reference_level = focal_params$.reference_level
    ))
  }

  # Read before the exposure, so that a fit whose scores cannot be read at all
  # is reported as such rather than after an announcement about an exposure
  # nothing goes on to use.
  ps_vec <- extract_model_propensity(model, "binary", call = call)

  focal_named <- !is.null(focal_params$.focal_level) ||
    !is.null(focal_params$.reference_level)

  if (needs_exposure || focal_named) {
    .exposure <- extract_exposure_from_model(model, .exposure)
  }

  invert <- !is.null(.exposure) &&
    !resolved_focal_is_default(
      .exposure,
      .focal_level = focal_params$.focal_level,
      .reference_level = focal_params$.reference_level
    )

  if (invert) {
    ps_vec <- 1 - ps_vec
  }

  list(
    propensity = ps_vec,
    exposure = .exposure,
    exposure_type = "binary",
    focal_level = focal_params$.focal_level,
    reference_level = focal_params$.reference_level
  )
}
