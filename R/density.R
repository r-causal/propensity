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

new_density_spec <- function(
  family,
  params = list(),
  fn = NULL,
  sigma_method = NULL
) {
  spec <- list(family = family, params = params, fn = fn)

  # How the spread of the standardized residuals is estimated is not a parameter
  # of the density: `t(df = 4)` is the same density however its scale was
  # arrived at. It sits beside the parameters, and only for a family that offers
  # a choice, so a family with one estimator carries no field saying which.
  if (!is.null(sigma_method)) {
    spec$sigma_method <- sigma_method
  }

  structure(spec, class = "propensity_density")
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
#'   tailed still, and heavier the smaller `df` is. Its scale is estimated by
#'   the root mean square of the residuals, as every other family's spread is,
#'   or by maximum likelihood under the t itself.
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
#' @param sigma_method How the spread of the conditional density is estimated
#'   from the residuals of the propensity score model, either `"rms"` or
#'   `"mle"`. See the section below.
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
#'   `family`, the name of the family; `params`, the parameters that identify
#'   it, which is an empty list for a family that takes none; and `fn`, the
#'   function that evaluates the density, which is `NULL` for `dens_kernel()`.
#'   `dens_t()` carries a fourth element, `sigma_method`, the name of the
#'   estimator its scale is read with.
#'
#' @section The spread of a t density:
#'
#' Both densities are evaluated on a residual standardized by a spread, and for
#' the conditional density that spread is estimated from the residuals of the
#' propensity score model. Every family takes the root mean square of those
#' residuals, which is `sigma_method = "rms"`, the default, and is the maximum
#' likelihood estimator of the spread of a normal density.
#'
#' It is not the maximum likelihood estimator of the scale of a t. The scale of
#' a t is smaller than its standard deviation by a factor that depends on `df`,
#' and the root mean square is pulled outward by the large residuals a heavy
#' tail produces, which is what the family was chosen to accommodate.
#' `sigma_method = "mle"` estimates the scale under the t itself, as the root of
#' \eqn{\sum_i \left[(\nu + 1) r_i^2 / (\nu \sigma^2 + r_i^2)\right] = n},
#' where \eqn{r_i} is the residual and \eqn{\nu} is `df`. Each residual enters
#' that sum through a bounded term, so a residual far out in the tail moves the
#' estimate by less than it moves the root mean square. Prefer it when the
#' residuals are heavy tailed, which is the case the t family is for; the two
#' estimators answer the same question, and answer it alike, as `df` grows and
#' the t approaches the normal.
#'
#' The choice describes the conditional density alone. The marginal density
#' that stabilizes the weights is the exposure's own, read at the exposure's
#' mean and root mean square, whatever the conditional spread was estimated by.
#' A scale estimated by maximum likelihood is recorded by [density_meta()] as
#' `sigma = "mle"`, and [ipw()] estimates it alongside the propensity score
#' model's coefficients, solving the equation above as part of its stacked
#' system so that the standard errors account for it. Supplying a `.sigma` says
#' the spread is a number of your own rather than one estimated from the
#' residuals, so the two cannot be given together.
#'
#' @examples
#' dens_normal()
#'
#' dens_t(df = 4)
#'
#' dens_t(df = 4, sigma_method = "mle")
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
dens_t <- function(df, sigma_method = c("rms", "mle")) {
  check_density_number(df, "df")
  sigma_method <- rlang::arg_match(sigma_method)
  force(df)

  new_density_spec(
    "t",
    params = list(df = df),
    fn = function(z) stats::dt(z, df = df),
    sigma_method = sigma_method
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
  if (
    !identical(x$family, y$family) ||
      !identical(x$params, y$params) ||
      !identical(x$sigma_method, y$sigma_method)
  ) {
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

# An infinite exposure or fitted value, refused where it arrives. A missing
# value is a unit with nothing to weight and leaves that unit's weight missing,
# which is a local answer to a local gap. An infinite one is not local: the
# pooled spread computed from it is infinite, so every weight comes back
# missing however few units carry the infinity, and under a kernel the
# residuals it leaves are reported as residuals that do not vary, which they
# do. Neither report names the value that caused it, so the value is named
# here, before anything is computed from it.
check_continuous_finite <- function(x, arg, call = rlang::caller_env()) {
  unusable <- !is.na(x) & !is.finite(x)

  if (!any(unusable)) {
    return(invisible(NULL))
  }

  abort(
    c(
      "Weights for a continuous exposure cannot be computed from an infinite
       value.",
      x = "{sum(unusable)} value{?s} of {.arg {arg}} {?is/are} infinite.",
      i = "An infinite value leaves the spread of the conditional density
           infinite, so every weight is missing rather than that unit's alone.
           Drop those observations before weighting."
    ),
    error_class = "propensity_density_error",
    call = call
  )
}

# The spread of the conditional density: the one the caller supplied, the scale
# of a t density estimated by maximum likelihood, or the pooled uncentered root
# mean square of the residuals. The root mean square is uncentered because the
# residuals of a fitted model already average to zero, and because the
# estimating equation `ipw()` solves for the same quantity is the uncentered
# moment.
continuous_sigma <- function(
  exposure,
  mu,
  .sigma = NULL,
  density = NULL,
  call = rlang::caller_env()
) {
  if (!is.null(.sigma)) {
    return(.sigma)
  }

  residuals <- exposure - mu

  if (density_sigma_is_mle(density)) {
    return(t_sigma_mle(residuals, density$params$df, call = call))
  }

  sqrt(mean(residuals^2, na.rm = TRUE))
}

# Whether a density asks for its scale to be estimated under itself rather than
# by the root mean square every family otherwise takes. Only `dens_t()` records
# the field, so every other specification answers no.
density_sigma_is_mle <- function(density) {
  identical(density$sigma_method, "mle")
}

# Where the spread of the conditional density came from, as the record keeps it.
# A spread the caller supplied is the caller's whatever the family would have
# estimated, so it is read first, and it is the pairing `check_sigma_method()`
# refuses.
density_sigma_source <- function(.sigma, density) {
  if (!is.null(.sigma)) {
    return("supplied")
  }

  if (density_sigma_is_mle(density)) {
    return("mle")
  }

  "pooled"
}

# The maximum likelihood scale of a Student's t with `df` degrees of freedom fit
# to the residuals that are there. A missing exposure or fitted value leaves a
# unit with no residual for the likelihood to read, exactly as it leaves that
# unit out of the root mean square.
#
# The score for the scale, multiplied through by the scale so that it is a mean
# of bounded terms, is the function below. It falls from `df` at a scale of
# nothing to -1 as the scale grows without bound, so the root is bracketed
# either side of the root mean square and found rather than searched for. The
# upper end is the root mean square scaled by enough to put the score below zero
# for any degrees of freedom: the mean is at most `(df + 1) / df` times the
# squared root mean square over the squared scale.
t_sigma_mle <- function(residuals, df, call = rlang::caller_env()) {
  residuals <- residuals[!is.na(residuals)]
  rms <- sqrt(mean(residuals^2))

  # Nothing to fit on, or residuals a model reproduced exactly. Neither leaves a
  # likelihood with a maximum to find, and both are the number the root mean
  # square returns for them, so the two estimators agree on the degenerate cases
  # and whatever the ratio does with such a spread it does either way.
  if (!is.finite(rms) || rms == 0) {
    return(rms)
  }

  score <- function(sigma) {
    mean((df + 1) * residuals^2 / (df * sigma^2 + residuals^2)) - 1
  }

  lower <- rms / 1e4
  upper <- rms * max(50, 2 * sqrt((df + 1) / df))

  # The likelihood has no maximum at a positive scale when so many residuals are
  # exactly zero that the score is negative however small the scale is: the fit
  # is degenerate rather than merely tight, and the root finder would report it
  # as a bracket that does not straddle a root.
  if (score(lower) <= 0) {
    n_zero <- sum(residuals == 0)
    abort(
      c(
        "The scale of a {.val t} density cannot be estimated by maximum
         likelihood from these residuals.",
        x = "{n_zero} of the {length(residuals)} residuals of the propensity
             score model {?is/are} exactly zero, so the likelihood grows
             without bound as the scale falls to zero.",
        i = "Use {.code sigma_method = \"rms\"}, or supply a spread with
             {.arg .sigma}."
      ),
      error_class = "propensity_density_error",
      call = call
    )
  }

  stats::uniroot(
    score,
    interval = c(lower, upper),
    tol = .Machine$double.eps^0.75
  )$root
}

# A spread the caller supplied and a spread estimated under the density are two
# instructions about the same quantity, and the second is not a family the
# weights take but an estimator that is then not run. They are refused together
# the way an integrated numerator refuses a `.sigma`, and for the same reason.
check_sigma_method <- function(.sigma, density, call = rlang::caller_env()) {
  if (is.null(.sigma) || !density_sigma_is_mle(density)) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.code sigma_method = \"mle\"} cannot be used with {.arg .sigma}.",
      x = "{.code sigma_method = \"mle\"} estimates the scale of the
           conditional density from the residuals of the propensity score
           model, and {.arg .sigma} is a spread of your own that replaces it.",
      i = "Drop {.arg .sigma} to estimate the scale under the t density, or
           build the density with {.code sigma_method = \"rms\"} to spread the
           one you supplied."
    ),
    error_class = "propensity_density_error",
    call = call
  )
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
  check_density_denominator(f_den, z, call = call)

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
# one for a model that conditions on covariates.
#
# `abs(mean(x))` is a second floor, for values that are large and constant: an
# offset of a billion leaves the last bits of the arithmetic differing by about
# 1e-6, which is a spread of nothing at that offset but far more than the
# residuals' own floor allows. It carries a much smaller coefficient than the
# scale term, because it is measuring the arithmetic rather than the model: at
# the coefficient the scale term uses, fitted means that differ by a unit or
# two at that offset read as one number, and a model that conditions on a
# covariate silently returns a weight of one for every unit.
is_constant <- function(x, scale) {
  diff(range(x)) <= max(sqrt(.Machine$double.eps) * scale, 1e-12 * abs(mean(x)))
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

  check_density_denominator(f_den, z, index = which(present), call = call)

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

  # The two ends are read one at a time below, so there have to be two of them.
  # A range of some other length was indexed for an end it does not have.
  if (length(range) != 2L) {
    abort(
      c(
        "{.arg .density} fits a kernel between the two ends of a range.",
        x = "It was given {length(range)} end{?s} rather than two.",
        i = "A kernel is fit from the lower end of the range to the upper one.
             Give both ends, in that order."
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

  # An infinite end passed the missingness test and reached `stats::density()`
  # as an error about its own `from` and `to`. A kernel is fit between two
  # numbers, so it is refused here in terms of the range it was given.
  if (!all(is.finite(range))) {
    abort(
      c(
        "{.arg .density} cannot fit a kernel over an infinite range.",
        x = "The range it was given runs from {.val {range[[1]]}} to
             {.val {range[[2]]}}.",
        i = "A kernel is fit between two finite ends. An infinite end comes
             from an infinite exposure or fitted value; drop those
             observations before weighting."
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

# Stabilized normal density-ratio weights have a finite second moment only while
# the marginal variance of the exposure stays below twice the variance of the
# conditional density. Past that boundary the weights have infinite variance,
# and estimates built from them are erratic however large the sample is, so the
# boundary is reported rather than left for the reader to work out.
#
# It is a property of the normal family read against the exposure's own marginal
# density, so it is reported for that configuration alone: the normal family, a
# marginal numerator, and one conditional spread describing every unit. A
# heavier-tailed family sits at a different boundary, an integrated numerator
# and a supplied stabilization score are not the marginal density, and a spread
# for each observation is not a single conditional variance the marginal
# variance can be read against.
#
# The report is a warning rather than a refusal because the weights are still
# the weights the estimand asks for; what fails is the precision of anything
# built from them.
check_density_variance <- function(
  density,
  sigma,
  sigma_a,
  numerator,
  call = rlang::caller_env()
) {
  reportable <- identical(numerator, "marginal") &&
    identical(density$family, "normal") &&
    length(sigma) == 1L &&
    length(sigma_a) == 1L &&
    !is.na(sigma) &&
    !is.na(sigma_a)

  if (!reportable || sigma_a^2 < 2 * sigma^2) {
    return(invisible(TRUE))
  }

  var_a <- signif(sigma_a^2, 3)
  var_d <- signif(sigma^2, 3)

  warn(
    c(
      "Stabilized normal weights have no finite variance for this exposure.",
      x = "The marginal variance of the exposure is {var_a}, which is at least
           twice the conditional variance {var_d}.",
      i = "The second moment of a stabilized normal density ratio exists only
           while the marginal variance stays below twice the conditional one.
           The weights are returned, but estimates built from them are erratic
           however large the sample is.",
      i = "A model that explains more of the exposure lowers the conditional
           variance, and a family with heavier tails, such as {.fun dens_t},
           has a different boundary."
    ),
    warning_class = "propensity_density_variance_warning",
    call = call
  )

  invisible(TRUE)
}

# The conditional density is what a continuous exposure's weights divide by, so
# a unit whose own exposure falls where that density is exactly zero has a
# weight of infinity. It is the continuous form of a propensity score of exactly
# 0 or 1, and is refused for the same reason: an infinite weight leaves every
# estimate built from it undefined.
#
# Only the denominator is held to this. A zero in the numerator is a marginal
# density that gives the unit no weight, which is a legitimate weight of zero.
#
# `index` maps the positions of `f_den` back to the observations they came from,
# for a caller that dropped the missing units before evaluating the density.
check_density_denominator <- function(
  f_den,
  z,
  index = NULL,
  call = rlang::caller_env()
) {
  zero <- !is.na(f_den) & f_den == 0
  if (!any(zero)) {
    return(invisible(TRUE))
  }

  units <- if (is.null(index)) which(zero) else index[zero]

  abort(
    c(
      "The conditional density is zero where the weight would divide by it.",
      x = "It is zero at {length(units)} observation{?s}: {units}.",
      i = density_range_hint(z, zero),
      i = "A weight built from a conditional density of zero is infinite. Use a
           family with heavier tails, such as {.fun dens_t} or
           {.fun dens_kernel}."
    ),
    error_class = "propensity_density_error",
    call = call
  )
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
