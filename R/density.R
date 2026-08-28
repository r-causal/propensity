# The conditional density family for continuous-exposure weights. A
# specification records the family, the parameters that identify it, and, for
# every family but the kernel, the function that evaluates it. The kernel has no
# function until the residuals it is fit on are known, so it carries its
# `stats::density()` controls and is fit by `density_eval()`.

# The bandwidth selectors `stats::density()` accepts by name.
density_bw_rules <- c("nrd0", "nrd", "ucv", "bcv", "SJ", "SJ-ste", "SJ-dpi")

# The kernels `stats::density()` accepts, read from the function itself so that
# the two stay in step.
density_kernels <- function() {
  eval(formals(stats::density.default)$kernel)
}

new_density_spec <- function(family, params = list(), fn = NULL) {
  structure(
    list(family = family, params = params, fn = fn),
    class = "propensity_density"
  )
}

#' Density specifications for continuous exposures
#'
#' @description
#' Weights for a continuous exposure are a ratio of densities. Both densities
#' are evaluated on the standardized residual
#' \eqn{z_i = (A_i - \mu_i) / \sigma}, where \eqn{A_i} is the exposure,
#' \eqn{\mu_i} the fitted conditional mean, and \eqn{\sigma} the residual
#' spread, so a specification describes a density on a standardized scale rather
#' than on the scale of the exposure itself.
#'
#' These constructors name the family and record the parameters that identify
#' it. Pass one to the `.density` argument of [wt_ate()] or [wt_cens()], which
#' also accept the strings `"normal"`, `"laplace"`, and `"kernel"` for the
#' families that need no parameters, and a bare function of one argument, which
#' is wrapped with `dens_fn()`.
#'
#' * `dens_normal()` is the standard normal density, the default.
#' * `dens_laplace()` is the standard Laplace density,
#'   \eqn{\exp(-|z|) / 2}, which puts more mass in the tails than the normal.
#' * `dens_t()` is Student's t density with `df` degrees of freedom, heavier
#'   tailed still, and heavier the smaller `df` is.
#' * `dens_kernel()` is a kernel density estimate of the standardized
#'   residuals, fit with [stats::density()] and interpolated to each
#'   observation. It assumes no family at all, at the cost of a density that
#'   is not a smooth function of the model's parameters: weights built from it
#'   have no closed-form standard error.
#' * `dens_fn()` is a density you write yourself.
#'
#' @param df Degrees of freedom for Student's t, a single positive, finite
#'   number.
#' @param bw The bandwidth passed to [stats::density()]: a single positive
#'   number, or the name of one of its selection rules
#'   (`"nrd0"`, `"nrd"`, `"ucv"`, `"bcv"`, `"SJ"`, `"SJ-ste"`, or `"SJ-dpi"`).
#' @param adjust A single positive number the bandwidth is multiplied by, as in
#'   [stats::density()]. Values above 1 smooth the estimate further.
#' @param kernel The smoothing kernel, one of the kernels [stats::density()]
#'   accepts.
#' @param n The number of grid points [stats::density()] evaluates, at least 2.
#'   The estimate is interpolated between them, so a larger `n` follows the
#'   shape of the residuals more closely.
#' @param f A function of one argument, the standardized residual, returning one
#'   non-negative, finite density value for each element it is given.
#'
#' @return An object of class `propensity_density`: a list with the elements
#'   `family`, `params`, and `fn`. `fn` is the function that evaluates the
#'   density, and is `NULL` for `dens_kernel()`.
#'
#' @examples
#' dens_normal()
#'
#' dens_t(df = 4)
#'
#' dens_kernel(adjust = 1.5)
#'
#' dens_fn(function(z) stats::dt(z, df = 4))
#'
#' @name dens_normal
#' @export
dens_normal <- function() {
  new_density_spec("normal", fn = function(z) stats::dnorm(z))
}

#' @rdname dens_normal
#' @export
dens_laplace <- function() {
  new_density_spec("laplace", fn = function(z) exp(-abs(z)) / 2)
}

#' @rdname dens_normal
#' @export
dens_t <- function(df) {
  check_density_number(df, "df")
  force(df)

  new_density_spec(
    "t",
    params = list(df = df),
    fn = function(z) stats::dt(z, df = df)
  )
}

#' @rdname dens_normal
#' @export
dens_kernel <- function(
  bw = "nrd0",
  adjust = 1,
  kernel = "gaussian",
  n = 512
) {
  bw <- check_density_bw(bw)
  check_density_number(adjust, "adjust")
  kernel <- rlang::arg_match0(kernel, density_kernels(), arg_nm = "kernel")
  check_density_n(n)

  new_density_spec(
    "kernel",
    params = list(bw = bw, adjust = adjust, kernel = kernel, n = n)
  )
}

#' @rdname dens_normal
#' @export
dens_fn <- function(f) {
  check_density_fn(f, arg = "f")

  new_density_spec("function", fn = f)
}

#' @export
format.propensity_density <- function(x, ...) {
  if (length(x$params) == 0) {
    return(x$family)
  }

  arguments <- paste0(
    names(x$params),
    " = ",
    vapply(x$params, format_density_param, character(1))
  )

  paste0(x$family, "(", paste(arguments, collapse = ", "), ")")
}

#' @export
print.propensity_density <- function(x, ...) {
  cat("<density: ", format(x), ">\n", sep = "")

  invisible(x)
}

# A parameter reads back the way it was written: a rule or kernel name in
# quotes, a number bare.
format_density_param <- function(value) {
  if (is.character(value)) {
    encodeString(value, quote = "\"")
  } else {
    format(value)
  }
}

# Coerce whatever the user supplied to `.density` into a specification. A
# specification is handed back as it came, so the object the user built is the
# object that is recorded.
as_density_spec <- function(
  x,
  arg = ".density",
  call = rlang::caller_env()
) {
  if (inherits(x, "propensity_density")) {
    return(x)
  }

  if (is.character(x)) {
    # Checked ahead of the match, which describes a longer vector in terms of
    # its own argument rather than the one the user wrote.
    if (length(x) != 1) {
      abort(
        c(
          "{.arg {arg}} must be a single string.",
          x = "It has length {length(x)}.",
          i = "Supply {.val normal}, {.val laplace}, or {.val kernel}."
        ),
        error_class = "propensity_density_error",
        call = call
      )
    }

    family <- rlang::arg_match0(
      x,
      c("normal", "laplace", "kernel"),
      arg_nm = arg,
      error_call = call
    )

    return(switch(
      family,
      normal = dens_normal(),
      laplace = dens_laplace(),
      kernel = dens_kernel()
    ))
  }

  if (is.function(x)) {
    # Checked here as well as in `dens_fn()` so that the refusal names the
    # argument the user wrote rather than the constructor it routes through.
    check_density_fn(x, arg = arg, call = call)

    return(dens_fn(x))
  }

  abort(
    c(
      "{.arg {arg}} must be a density specification, a function, or a string.",
      x = "It has class {.cls {class(x)}}.",
      i = "Supply {.val normal}, {.val laplace}, or {.val kernel}, one of
           {.fun dens_normal}, {.fun dens_laplace}, {.fun dens_t},
           {.fun dens_kernel}, or {.fun dens_fn}, or a function of the
           standardized residual."
    ),
    error_class = "propensity_density_error",
    call = call
  )
}

# Evaluate a specification at the standardized residuals `z`. A kernel is fit
# once on `fit_on`, bounded by `range`, and interpolated to `z`, so `range` must
# cover every point the density is asked for: a point outside it interpolates to
# a missing value, which the output check refuses.
density_eval <- function(
  spec,
  z,
  fit_on = z,
  # Qualified because the argument's own name would otherwise shadow the
  # function the default calls.
  range = base::range(z),
  call = rlang::caller_env()
) {
  values <- if (identical(spec$family, "kernel")) {
    check_density_kernel_fit(fit_on, range, call = call)

    fit <- stats::density(
      fit_on,
      bw = spec$params$bw,
      adjust = spec$params$adjust,
      kernel = spec$params$kernel,
      n = spec$params$n,
      from = range[1],
      to = range[2]
    )

    stats::approxfun(fit$x, fit$y)(z)
  } else {
    # A density written with matrix arithmetic returns its values with a
    # dimension attached. The values are what is wanted, one per observation, so
    # the shape is dropped rather than held against the function.
    as.vector(spec$fn(z))
  }

  check_density_values(values, length(z), z = z, call = call)

  values
}

# What a kernel needs before it can be fit at all. `stats::density()` and
# `stats::approxfun()` each refuse some of these on their own terms, in messages
# written about their own arguments rather than about the exposure the residuals
# came from, and a range of zero width reaches `approxfun()` as a warning before
# it reaches it as an error. They are refused here instead, once, in the terms
# the caller supplied them in.
check_density_kernel_fit <- function(
  fit_on,
  range,
  call = rlang::caller_env()
) {
  if (anyNA(fit_on)) {
    n_missing <- sum(is.na(fit_on))
    abort(
      c(
        "{.arg .density} cannot fit a kernel on missing standardized
         residuals.",
        x = "{n_missing} of them {?is/are} missing.",
        i = "A missing exposure or fitted value leaves a missing residual.
             Drop those observations before weighting."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (anyNA(range)) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel over a missing range.",
        x = "The range a kernel is fit over is read from the standardized
             residuals, and at least one of them is missing.",
        i = "A missing exposure or fitted value leaves a missing residual.
             Drop those observations before weighting."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (sum(is.finite(fit_on)) < 2) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel on fewer than two standardized
         residuals.",
        x = "It was given {sum(is.finite(fit_on))} to fit on.",
        i = "The bandwidth is chosen from the spread of the residuals, which
             takes at least two of them. Use a parametric family, such as
             {.fun dens_normal}, on a sample this small."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (range[[2]] <= range[[1]]) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel on standardized residuals that do
         not vary.",
        x = "They are every one {.val {range[[1]]}}, so the estimate has no
             range to be fit over.",
        i = "A model that fits the exposure exactly, or an exposure that is
             constant, leaves residuals a kernel cannot smooth. Use a
             parametric family, such as {.fun dens_normal}."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  invisible(TRUE)
}

# The output check every density passes, whoever wrote it. `z` is the
# standardized residuals the density was evaluated at, and is used only to say
# how extreme the exposure was where the density failed.
check_density_values <- function(
  p,
  n,
  what = ".density",
  z = NULL,
  call = rlang::caller_env()
) {
  if (!is.numeric(p)) {
    abort(
      c(
        "{.arg {what}} must return a numeric vector.",
        x = "It returned {.cls {class(p)}}.",
        i = "A density returns one number for each standardized residual it is
             given."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (length(p) != n) {
    abort(
      c(
        "{.arg {what}} must return one value for each observation.",
        x = "It returned {length(p)} value{?s} for {n} observation{?s}.",
        i = "The density is evaluated on the whole vector of standardized
             residuals at once, so it must be vectorized over its argument."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (anyNA(p)) {
    n_missing <- sum(is.na(p))
    abort(
      c(
        "{.arg {what}} must not return missing values.",
        x = "It returned {n_missing} missing value{?s}.",
        i = density_range_hint(z),
        i = "The density must be defined at every standardized residual, and a
             kernel density is defined only over the range it was fit on."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  unusable <- !is.finite(p) | p < 0
  if (any(unusable)) {
    n_unusable <- sum(unusable)
    abort(
      c(
        "{.arg {what}} must return finite, non-negative values.",
        x = "It returned {n_unusable} value{?s} that {?is/are} negative or not
             finite.",
        i = density_range_hint(z, unusable),
        i = "A density is non-negative everywhere. Check that {.arg {what}}
             returns a density rather than a log density or a distribution
             function."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (length(p) > 0 && all(p == 0)) {
    abort(
      c(
        "{.arg {what}} must not be zero everywhere.",
        x = "Every one of the {length(p)} value{?s} it returned is zero.",
        i = density_range_hint(z),
        i = "Every weight built from a density that is zero at every
             observation is infinite. Check that the family matches the scale
             of the standardized residuals."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  invisible(TRUE)
}

# How extreme the standardized residuals were, over all of them or over the
# subset the density failed at. Written as plain text so that it can be passed
# to `abort()` as a bullet without being interpolated a second time.
density_range_hint <- function(z, at = NULL) {
  if (!is.numeric(z) || length(z) == 0) {
    return(NULL)
  }

  if (!is.null(at)) {
    if (length(at) != length(z)) {
      return(NULL)
    }

    z <- z[at]
  }

  z <- z[is.finite(z)]
  if (length(z) == 0) {
    return(NULL)
  }

  lowest <- signif(min(z), 3)
  highest <- signif(max(z), 3)

  if (lowest == highest) {
    return(paste0(
      if (is.null(at)) {
        "The standardized residual is "
      } else {
        "It failed at the standardized residual "
      },
      lowest,
      "."
    ))
  }

  paste0(
    if (is.null(at)) {
      "The standardized residuals run from "
    } else {
      "It failed at standardized residuals from "
    },
    lowest,
    " to ",
    highest,
    "."
  )
}

# A single positive, finite number: the shape both `df` and `adjust` take.
check_density_number <- function(x, arg, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    abort(
      c(
        "{.arg {arg}} must be a single positive, finite number.",
        x = density_value_bullet(x),
        i = "Supply one number, greater than zero."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  invisible(TRUE)
}

# The grid size, hardwired to `dens_kernel()`'s `n` because that is the only
# argument that takes this shape. The value keeps the name the bullet helper
# writes its templates against.
check_density_n <- function(x, call = rlang::caller_env()) {
  usable <- is.numeric(x) &&
    length(x) == 1 &&
    is.finite(x) &&
    x == trunc(x) &&
    x >= 2

  if (!usable) {
    abort(
      c(
        "{.arg n} must be a single whole number of at least 2.",
        x = density_value_bullet(x),
        i = "{.arg n} is the number of points {.fun stats::density} evaluates
             the estimate at before it is interpolated."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  invisible(TRUE)
}

# A bandwidth is either a number or the name of a selection rule. The names are
# matched the way an argument is, so a near miss is corrected; a number is held
# to the same shape as the other numeric controls.
check_density_bw <- function(bw, call = rlang::caller_env()) {
  if (is.character(bw)) {
    return(rlang::arg_match0(
      bw,
      density_bw_rules,
      arg_nm = "bw",
      error_call = call
    ))
  }

  check_density_number(bw, "bw", call = call)

  bw
}

check_density_fn <- function(x, arg, call = rlang::caller_env()) {
  if (!is.function(x)) {
    abort(
      c(
        "{.arg {arg}} must be a function.",
        x = density_value_bullet(x),
        i = "Supply a function of the standardized residual, such as
             {.code function(z) dnorm(z)}."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  # `args()` gives a primitive its formals, which it does not otherwise have.
  if (length(formals(args(x))) == 0) {
    abort(
      c(
        "{.arg {arg}} must be a function of at least one argument.",
        x = "It takes no arguments.",
        i = "The density is called on the standardized residuals, so its first
             argument is the vector it is evaluated at."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  invisible(TRUE)
}

# What was supplied, described the way it went wrong. The bullet is a template
# interpolated in the frame that raises the error, so every validator holds the
# value it is checking in a local named `x`.
density_value_bullet <- function(x) {
  if (!is.numeric(x)) {
    "It has class {.cls {class(x)}}."
  } else if (length(x) != 1) {
    "It has length {length(x)}."
  } else {
    "It is {.val {x}}."
  }
}
