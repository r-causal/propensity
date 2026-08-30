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
#'   have no closed-form standard error, and [ipw()] reports none for them at
#'   all. Bootstrap the whole fit by hand to put an interval around such an
#'   estimate.
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

# Two specifications describe the same density when they name the same family
# with the same parameters. The function is left out of the comparison for the
# families that build one: each constructor writes a fresh closure, so two
# specifications built the same way in two calls hold functions that evaluate
# alike and are not `identical()`. A density the user wrote has no parameters to
# be told apart by, and the function on it is the object the user supplied
# rather than one a constructor built, so there it is the comparison.
density_specs_agree <- function(x, y) {
  if (!identical(x$family, y$family) || !identical(x$params, y$params)) {
    return(FALSE)
  }

  !identical(x$family, "function") || identical(x$fn, y$fn)
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

# A density other than the default describes a continuous exposure and has
# nothing to say about any other kind, so it is refused there the way
# `check_sigma()` refuses a spread. The default is accepted for every type and
# ignored outside the continuous route, so that a caller who writes out the
# family they were getting anyway is not refused for saying so.
check_density_arg <- function(
  spec,
  exposure_type,
  call = rlang::caller_env()
) {
  if (
    identical(exposure_type, "continuous") || identical(spec$family, "normal")
  ) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.arg .density} applies only to continuous exposures.",
      x = "{.arg .exposure} is being treated as {exposure_type}.",
      i = "{.arg .density} chooses the family of the conditional density a
           continuous exposure's weights are a ratio of. A {exposure_type}
           exposure has a probability rather than a density, so leave
           {.arg .density} unset for one."
    ),
    error_class = "propensity_density_error",
    call = call
  )
}

# What a numerator has to be able to be. The default is accepted for every
# exposure type and ignored outside the continuous route, so a caller who writes
# out the numerator the weights were getting anyway is not refused for saying
# so. An integrated numerator marginalizes the conditional density a continuous
# exposure's weights divide by, and each of the settings refused below leaves it
# with no such density to marginalize.
check_numerator <- function(
  numerator,
  exposure_type,
  stabilize,
  stabilization_score,
  .sigma,
  call = rlang::caller_env()
) {
  if (!identical(numerator, "integrated")) {
    return(invisible(NULL))
  }

  if (!identical(exposure_type, "continuous")) {
    abort(
      c(
        "{.arg numerator} = {.val integrated} applies only to continuous
         exposures.",
        x = "{.arg .exposure} is being treated as {exposure_type}.",
        i = "An integrated numerator averages the conditional density of a
             continuous exposure over the units. A {exposure_type} exposure has
             a probability rather than a density, so leave {.arg numerator}
             unset for one."
      ),
      error_class = "propensity_numerator_error",
      call = call
    )
  }

  if (!isTRUE(stabilize)) {
    abort(
      c(
        "{.arg numerator} = {.val integrated} needs stabilized weights.",
        x = "{.arg stabilize} is {.code FALSE}, so the weights carry no
             numerator at all.",
        i = "Set {.code stabilize = TRUE} to stabilize the weights on the
             marginalized conditional density, or leave {.arg numerator} unset."
      ),
      error_class = "propensity_numerator_error",
      call = call
    )
  }

  if (!is.null(stabilization_score)) {
    abort(
      c(
        "{.arg numerator} = {.val integrated} cannot be used with
         {.arg stabilization_score}.",
        x = "A score you supply is itself the numerator of the weights.",
        i = "Drop {.arg stabilization_score} to stabilize on the marginalized
             conditional density, or leave {.arg numerator} unset to keep the
             numerator you wrote."
      ),
      error_class = "propensity_numerator_error",
      call = call
    )
  }

  if (!is.null(.sigma)) {
    abort(
      c(
        "{.arg numerator} = {.val integrated} cannot be used with
         {.arg .sigma}.",
        x = "The marginalization reads the conditional density the propensity
             model estimated at every unit's fitted mean, and {.arg .sigma}
             replaces the spread of that density with one of your own.",
        i = "Leave {.arg .sigma} unset to spread the conditional density by the
             pooled residual root mean square, or use {.arg numerator} =
             {.val marginal}, which takes a spread you supply."
      ),
      error_class = "propensity_numerator_error",
      call = call
    )
  }

  invisible(NULL)
}

# The spread of the conditional density: the one the caller supplied, or the
# pooled uncentered root mean square of the residuals. It is uncentered because
# the residuals of a fitted model already average to zero, and because the
# estimating equation `ipw()` solves for the same quantity is the uncentered
# moment.
continuous_sigma <- function(exposure, mu, .sigma = NULL) {
  if (!is.null(.sigma)) {
    return(.sigma)
  }

  sqrt(mean((exposure - mu)^2, na.rm = TRUE))
}

# The density ratio a continuous exposure's weights are. Both densities are
# evaluated on a standardized residual and divided by the spread that
# standardized it, which is the Jacobian that puts each back on the exposure's
# own scale, so the normal family returns exactly what a normal density in the
# exposure's units returns.
#
# The numerator is the marginal density of the exposure, the conditional density
# marginalized over the units, a stabilization score the caller supplied, or
# nothing at all. `grid` is the evaluation grid an integrated numerator averages
# the conditional density over, and is built from the exposure when it is not
# supplied.
continuous_density_ratio <- function(
  exposure,
  mu,
  sigma,
  density,
  numerator = c("marginal", "integrated", "none", "score"),
  mu_a = NULL,
  sigma_a = NULL,
  score = NULL,
  grid = NULL,
  call = rlang::caller_env()
) {
  numerator <- rlang::arg_match(numerator, error_call = call)

  if (identical(numerator, "integrated")) {
    return(continuous_integrated_ratio(
      exposure = exposure,
      mu = mu,
      sigma = sigma,
      density = density,
      grid = grid,
      call = call
    ))
  }

  z <- (exposure - mu) / sigma
  f_den <- density_eval_present(density, z, call = call) / sigma

  if (identical(numerator, "none")) {
    return(1 / f_den)
  }

  if (identical(numerator, "score")) {
    return(score / f_den)
  }

  # The marginal density is the same family read at the exposure's own center
  # and spread, so it is standardized by those rather than by the conditional
  # spread. A kernel numerator is therefore fit on its own values.
  z_a <- (exposure - mu_a) / sigma_a
  f_num <- density_eval_present(density, z_a, call = call) / sigma_a

  f_num / f_den
}

# The number of points an integrated numerator marginalizes the conditional
# density over. It is the grid WeightIt uses, which is what the agreement with
# it is written against. A grid four times as long removed the negative
# interpolated densities the simulation study behind the default saw in its
# heavy-tailed cell, and changed nothing else it measured.
continuous_grid_n <- 50L

# Whether a vector holds one number, to within the arithmetic that produced it.
# Fitted values that ought to be identical seldom are: an intercept-only model
# reaches each of them through a decomposition, and the last few bits differ, so
# counting distinct doubles would find several where the model means one.
#
# `scale` is the size the spread of `x` is read against, and is the spread of
# the residuals wherever this is called from. A density ratio is the same
# number whatever unit the exposure is measured in, so the test that decides
# which ratio to return has to be the same number too: a fixed floor of one
# would read an exposure in nanograms as a single value and return weights of
# one for a model that conditions on covariates. `abs(mean(x))` is kept in the
# comparison so that values which are large and constant are still read as
# constant, as `all.equal()` reads them.
is_constant <- function(x, scale) {
  diff(range(x)) <= sqrt(.Machine$double.eps) * max(scale, abs(mean(x)))
}

# The density ratio under an integrated numerator, which is the conditional
# density in both places: the denominator reads it at each unit's own fitted
# mean, and the numerator reads the average of it over the units.
#
# A unit whose exposure or fitted mean is missing has no standardized residual,
# so it is not read by the average, its exposure does not reach the ends of the
# grid, and its weight is missing, exactly as it is under any other numerator.
continuous_integrated_ratio <- function(
  exposure,
  mu,
  sigma,
  density,
  grid = NULL,
  call = rlang::caller_env()
) {
  z <- (exposure - mu) / sigma
  present <- !is.na(z)

  exposure <- exposure[present]
  mu <- mu[present]
  z <- z[present]

  wt <- rep(NA_real_, length(present))

  # With nothing to condition on, every unit's conditional density is the same
  # density, the average of it over the units is that density again, and the
  # ratio is one. Sending it through the grid and the interpolation instead
  # returns values near but not equal to one, so the case is answered directly
  # rather than to within the grid's error.
  #
  # The fitted values of an intercept-only model are constant to within the
  # arithmetic that produced them rather than exactly, so the case is read from
  # the spread of the fitted means rather than from how many distinct doubles
  # they are. Fitted means that vary by that little leave weights of one either
  # way; this returns the ones the model means. How little that is depends on
  # the units the exposure is measured in, so the spread of the fitted means is
  # read against the spread of the residuals rather than against one.
  if (length(mu) == 0 || is_constant(mu, scale = sigma)) {
    wt[present] <- 1
    return(wt)
  }

  if (is.null(grid)) {
    grid <- seq(min(exposure), max(exposure), length.out = continuous_grid_n)
  }

  # An exposure that takes one value leaves the grid no width, and the
  # interpolation nowhere to run. `stats::spline()` reaches that as a base
  # warning about collapsing tied points and then returns the density at the one
  # point for every unit; it is refused here instead, in terms of the exposure
  # the grid was built from.
  if (grid[[length(grid)]] <= grid[[1]]) {
    abort(
      c(
        "{.arg numerator} = {.val integrated} cannot marginalize over an
         exposure that does not vary.",
        x = "Every exposure with a fitted conditional mean is
             {.val {grid[[1]]}}, so the grid the conditional density is
             averaged over has no width.",
        i = "The integrated numerator reads the conditional density at points
             spanning {.arg .exposure}. Use {.arg numerator} = {.val marginal}
             for an exposure that takes one value."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  # The conditional density of each unit read at each grid point, on the
  # standardized scale the family is written on. The spread is a single number
  # here, which is why an integrated numerator refuses a `.sigma` of one per
  # observation.
  standardized <- outer(grid, mu, "-") / sigma

  # A kernel is one estimate answering for both densities, so it is fit on the
  # standardized residuals over a range that covers them and the whole
  # standardized grid, which reaches further than they do. The parametric
  # families have no fit and ignore both, so the range is read only for a
  # kernel: the grid holds one value per unit per grid point, and concatenating
  # it with the residuals to take their range would copy all of it.
  span <- if (identical(density$family, "kernel")) {
    range(range(z), range(standardized))
  } else {
    NULL
  }

  f_den <- density_eval(
    density,
    z,
    fit_on = z,
    range = span,
    call = call
  ) /
    sigma

  f_num <- continuous_numerator_integrated(
    exposure = exposure,
    z = z,
    grid = grid,
    standardized = standardized,
    sigma = sigma,
    density = density,
    span = span,
    call = call
  )

  wt[present] <- f_num / f_den
  wt
}

# The integrated numerator: the conditional density averaged over the units at
# each grid point, \eqn{f_A(t_j) = n^{-1} \sum_i g((t_j - \mu_i) / \sigma) /
# \sigma}, interpolated back to each observed exposure with a cubic spline.
#
# The interpolation is the reason the result is checked rather than trusted. A
# spline through points that come close to zero can undershoot below it, which
# is a property of the interpolation and not of the family, so the values are
# held to the same rules any density is held to before they become weights.
continuous_numerator_integrated <- function(
  exposure,
  z,
  grid,
  standardized,
  sigma,
  density,
  span,
  call = rlang::caller_env()
) {
  on_grid <- rowMeans(matrix(
    density_eval(
      density,
      as.vector(standardized),
      fit_on = z,
      range = span,
      call = call
    ),
    nrow = length(grid)
  )) /
    sigma

  values <- stats::spline(grid, on_grid, xout = exposure, method = "fmm")$y

  check_density_values(
    values,
    length(exposure),
    what = "The integrated numerator",
    what_is_arg = FALSE,
    z = z,
    remedy = "The marginalized density is interpolated back to the exposure
              with a cubic spline, which can dip below zero where the density
              on the grid comes close to it. Use {.arg numerator} =
              {.val marginal}, or a density with a heavier tail, to stabilize
              on a density that is positive everywhere.",
    call = call
  )

  values
}

# A density is asked only about the observations it can answer for. A missing
# exposure, fitted value, or spread leaves a missing standardized residual, and
# the weight there is missing for that reason rather than for anything the
# density did: a propensity score trimmed away is the ordinary case. Asking only
# about the residuals that are present also lets a kernel be fit on the
# observations that have one.
density_eval_present <- function(spec, z, call = rlang::caller_env()) {
  present <- !is.na(z)

  if (all(present)) {
    return(density_eval(spec, z, call = call))
  }

  values <- rep(NA_real_, length(z))

  # Nothing is present to ask about, so nothing is asked. Every family returns
  # a missing value at a missing residual, and a kernel asked to fit on no
  # residuals at all would instead be refused for the sample size, describing a
  # density that was never the reason the weights are missing.
  if (!any(present)) {
    return(values)
  }

  values[present] <- density_eval(spec, z[present], call = call)

  values
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
  # The sample size is read first because the range a kernel is fit over
  # defaults to the range of the residuals, and `base::range()` of nothing at
  # all warns about its missing arguments twice before it returns a pair of
  # infinities. Nothing is learned from that pair, and the count below says
  # everything there is to say about it.
  if (length(fit_on) < 2) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel on fewer than two standardized
         residuals.",
        x = "It was given {length(fit_on)} to fit on.",
        i = "The bandwidth is chosen from the spread of the residuals, which
             takes at least two of them. Use a parametric family, such as
             {.fun dens_normal}, on a sample this small."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  # A residual that is infinite is refused alongside one that is missing. An
  # infinite residual among the values the ends of the estimate are read from
  # reaches `stats::density()` as an error about its own `from` and `to`; one
  # that arrives only in `fit_on` reaches nothing at all, and leaves a kernel
  # fit over a range its residuals overflow, which integrates to less than one
  # with nothing said.
  if (!all(is.finite(fit_on))) {
    n_unusable <- sum(!is.finite(fit_on))
    abort(
      c(
        "{.arg .density} cannot fit a kernel on missing or infinite
         standardized residuals.",
        x = "{n_unusable} of them {?is/are} missing or infinite.",
        i = "A missing exposure or fitted value leaves a missing residual, and
             an infinite exposure or a spread of zero leaves an infinite one.
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

  # Ends given the wrong way round are a range the caller wrote backwards, not
  # residuals that do not vary, and are said to be that.
  if (range[[2]] < range[[1]]) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel over a reversed range.",
        x = "The range it was given runs from {.val {range[[1]]}} down to
             {.val {range[[2]]}}.",
        i = "A kernel is fit from the lower end of the range to the upper one.
             Give the ends in that order."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (range[[2]] == range[[1]]) {
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

# The output check every density passes, whoever wrote it. `what` names the
# thing the values came from: an argument the caller wrote, styled as one, or,
# under `what_is_arg = FALSE`, a noun phrase for a density the package computed
# itself. `z` is the standardized residuals the density was evaluated at, and is
# used only to say how extreme the exposure was where the density failed.
# `remedy` is how to fix a density that came back negative, which a family the
# caller chose and a density read off an interpolation are fixed in different
# ways.
check_density_values <- function(
  p,
  n,
  what = ".density",
  what_is_arg = TRUE,
  z = NULL,
  remedy = "A density is non-negative everywhere. Check that {.arg {what}}
            returns a density rather than a log density or a distribution
            function.",
  call = rlang::caller_env()
) {
  # Styled once, here, rather than by each message below: a template cannot
  # interpolate markup that arrives in a value, so the subject reaches the
  # messages already written the way it reads. The styling is deferred until a
  # message is built, because the check almost always succeeds and formatting
  # the subject costs more than every test below it put together.
  delayedAssign(
    "subject",
    if (what_is_arg) cli::format_inline("{.arg {what}}") else what
  )

  if (!is.numeric(p)) {
    abort(
      c(
        "{subject} must return a numeric vector.",
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
        "{subject} must return one value for each observation.",
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
        "{subject} must not return missing values.",
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
        "{subject} must return finite, non-negative values.",
        x = "It returned {n_unusable} value{?s} that {?is/are} negative or not
             finite.",
        i = density_range_hint(z, unusable),
        i = remedy
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  if (length(p) > 0 && all(p == 0)) {
    abort(
      c(
        "{subject} must not be zero everywhere.",
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
