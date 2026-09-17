#' Trim Propensity Scores
#'
#' @description
#' Trim observations with extreme propensity scores by replacing them with `NA`,
#' effectively removing those units from downstream analyses. The returned object
#' has the same length (or dimensions) as the input, with trimmed entries set to
#' `NA`. After trimming, refit the propensity score model on the retained
#' observations with [ps_refit()].
#'
#' @param .propensity A numeric vector of propensity scores in (0, 1) for binary
#'   exposures, a matrix / data frame where each column gives the propensity
#'   score for one level of a categorical exposure, or a fitted model (see
#'   **Fitted models** and **Dose models** in Details). A data frame trimmed
#'   for a binary exposure is reduced to a single column: the second column of
#'   a two column data frame, which is the probability of the second level in
#'   the layout model predictions come in, and the first column otherwise. The
#'   column taken is announced; `options(propensity.quiet = TRUE)` silences the
#'   announcement. A matrix is held to the same open interval as a vector, so a
#'   score of exactly 0 or 1 in any cell is refused and a separated multinomial
#'   fit cannot be repaired by trimming it: setting an extreme score to missing
#'   gains nothing from a score already at an endpoint. [ps_trunc()] reads the
#'   closed interval for a categorical matrix and is the repair; see
#'   **Propensity scores at 0 and 1** in [wt_ate()]. A data frame holding a
#'   `.pred_class` column, which a fitted tidymodels classification model
#'   returns when no prediction type is named, carries predicted levels rather
#'   than probabilities and is refused with an error of class
#'   `propensity_df_class_column_error`.
#' @param method Trimming method. One of:
#'
#'   * **`"ps"`** (default): Fixed threshold. Observations with propensity
#'     scores outside `[lower, upper]` are trimmed. For categorical exposures,
#'     observations where *any* column falls below `lower` (the symmetric
#'     threshold delta) are trimmed.
#'   * **`"adaptive"`**: Data-driven threshold that minimizes the asymptotic
#'     variance of the IPW estimator (Crump et al., 2009). The `lower` and
#'     `upper` arguments are ignored.
#'   * **`"pctl"`**: Quantile-based. Observations outside the `[lower, upper]`
#'     quantiles of the propensity score distribution are trimmed. Defaults:
#'     `lower = 0.05`, `upper = 0.95`.
#'   * **`"pref"`**: Preference score trimming. Transforms propensity scores
#'     to the preference scale (Walker et al., 2013) and trims outside
#'     `[lower, upper]`. Requires `.exposure`. Binary exposures only. Defaults:
#'     `lower = 0.3`, `upper = 0.7`.
#'   * **`"cr"`**: Common range (clinical equipoise). Trims to the overlap
#'     region of the propensity score distributions across exposure groups.
#'     Requires `.exposure`. Binary exposures only. The `lower` and `upper`
#'     arguments are ignored. When the two distributions do not overlap at all,
#'     so that the lowest score among the focal units sits above the highest
#'     among the reference units, the overlap region is empty and every observed
#'     unit is trimmed. That is a truthful record of an empty region rather than
#'     an error, and it differs deliberately from [ps_trunc()], which refuses
#'     the same data with an error of class `propensity_no_overlap_error`
#'     because there is no range left to bound the scores to.
#'   * **`"optimal"`**: Multi-category optimal trimming (Yang et al., 2016).
#'     Categorical exposures only. Requires `.exposure`.
#'   * **`"density"`**: Quantile floor on the conditional density of a
#'     continuous exposure. Units whose conditional density at their observed
#'     dose falls below its `lower` quantile are trimmed. Default:
#'     `lower = 0.01`.
#'   * **`"resid"`**: Bound on the absolute standardized residual of a
#'     continuous exposure. Units more than `upper` spreads from their
#'     predicted dose are trimmed. `upper` has no default.
#'
#'   `"density"` and `"resid"` need the residuals and the family of the model
#'   that fit the exposure's conditional mean, so they accept only a dose model
#'   (see **Dose models** in Details). They refuse a vector or matrix of values,
#'   a binomial or quasibinomial `glm`, and a `multinom` with an error of class
#'   `propensity_method_error`. A dose model accepts only these two methods, and
#'   refuses every other one, including the default `"ps"` when `method` is not
#'   supplied, with the same class.
#'
#'   For categorical exposures, only `"ps"` and `"optimal"` are supported.
#' @param lower,upper Numeric thresholds whose interpretation depends on
#'   `method`:
#'
#'   * `"ps"`: absolute propensity score bounds (defaults: 0.1, 0.9). For
#'     categorical exposures, only `lower` is used, as the symmetric threshold
#'     delta, and it defaults to 0.1, the default [ps_trunc()] uses as well.
#'     With `k` exposure levels, a threshold of `1/k` or larger cannot be met by
#'     every column of a row that sums to one, and is an error.
#'   * `"pctl"`: quantile probabilities (defaults: 0.05, 0.95).
#'   * `"pref"`: preference score bounds (defaults: 0.3, 0.7).
#'   * `"adaptive"`, `"cr"`, `"optimal"`: ignored (thresholds are data-driven).
#'   * `"density"`: `lower` is the quantile probability of the floor, in
#'     (0, 0.5) (default 0.01). `upper` is refused with an error of class
#'     `propensity_unsupported_arg_error`.
#'   * `"resid"`: `upper` is the bound on the absolute standardized residual, a
#'     single positive finite number, and is required (an error of class
#'     `propensity_missing_arg_error` without it). `lower` is refused with an
#'     error of class `propensity_unsupported_arg_error`.
#'
#'   A `"density"` `lower` or `"resid"` `upper` out of range is refused with an
#'   error of class `propensity_range_error`.
#' @param .exposure An exposure variable. Required for `"pref"`, `"cr"` (binary
#'   vector), and `"optimal"` (factor or character). For `"density"` and
#'   `"resid"`, the dose is read from the model's response unless supplied here.
#'   Not required for other methods.
#' @inheritParams wt_ate
#' @param .focal_level The value of `.exposure` representing the focal
#'   (treated) group, used by `"pref"` and `"cr"`. Every binary coding honors
#'   it: 0/1 numeric, logical, two-level factor, and two-level character
#'   exposures are all coded with the named level as focal, and a level the
#'   exposure never takes is an error. With no level named, a binary exposure
#'   defaults to its higher level, which is `1` for a 0/1 exposure, `TRUE` for
#'   a logical one, and the second of the two levels a factor or character
#'   exposure takes. Levels a factor declares but never takes are not
#'   candidates. Naming any other level reverses the coding, so `.propensity`
#'   must then hold the probability of the named level.
#' @param .reference_level The value of `.exposure` representing the reference
#'   (control) group. Naming it makes the exposure's other level focal, with
#'   the same consequence for `.propensity`, and a level the exposure never
#'   takes is an error. Automatically detected if not supplied.
#' @param ... Additional arguments passed to methods.
#' @param .sigma For `"density"` and `"resid"`, a single residual spread to
#'   read the conditional density at, recorded as `"supplied"`. With none
#'   supplied, the spread is the one the density family estimates from the
#'   residuals, as in [wt_ate()]: the root mean square, unless the family
#'   estimates a scale of its own. A spread for each unit is refused with an
#'   error of class `propensity_sigma_error`, and a spread supplied with a
#'   family that estimates its own scale is refused with an error of class
#'   `propensity_density_error`. The other methods refuse any `.sigma` with an
#'   error of class `propensity_sigma_error`.
#' @param .density For `"density"` and `"resid"`, the family of the conditional
#'   density, in any form [wt_ate()] accepts: `"normal"` (the default),
#'   `"laplace"`, `"kernel"`, a specification such as [dens_t()], or a function
#'   of the standardized residual. `"resid"` accepts only [dens_normal()],
#'   [dens_t()], and [dens_laplace()], and refuses a kernel or user-written
#'   density with an error of class `propensity_density_error`. The other
#'   methods refuse any family but the normal with the same class.
#' @param ps `r lifecycle::badge("deprecated")` Use `.propensity` instead. A
#'   call that names `ps` must name the arguments after it as well, since a
#'   positional argument binds to `.propensity`.
#'
#' @details
#' ## How trimming works
#'
#' Trimming identifies observations with extreme (near 0 or 1) propensity
#' scores and sets them to `NA`. These observations are excluded from
#' subsequent weight calculations and effect estimation. The goal is to
#' remove units that lack sufficient overlap between exposure groups, which
#' would otherwise receive extreme weights and destabilize estimates.
#'
#' ## Choosing a method
#'
#' * Use `"ps"` when you have a specific threshold in mind or want a simple
#'   default.
#' * Use `"adaptive"` for a principled, data-driven cutoff that targets
#'   variance reduction.
#' * Use `"pctl"` to trim a fixed percentage of extreme values from each tail.
#' * Use `"pref"` when you want to restrict to the region of clinical
#'   equipoise based on the preference score.
#' * Use `"cr"` to restrict to the common support region where both exposure
#'   groups have observed propensity scores.
#' * Use `"optimal"` for multi-category (3+) exposures; this is the only
#'   data-driven method available for categorical treatments.
#' * Use `"density"` or `"resid"` for a continuous exposure, on its dose model.
#'
#' ## Typical workflow
#'
#' 1. Fit a propensity score model
#' 2. Apply `ps_trim()` to flag extreme values
#' 3. Call [ps_refit()] to re-estimate propensity scores on the retained sample
#' 4. Compute weights with [wt_ate()] or another weight function
#'
#' ## Fitted models
#'
#' `.propensity` can be the fitted propensity score model instead of the scores
#' it reports. A binomial [stats::glm()] and a two-level `nnet::multinom()` are
#' read as one score per unit; a `nnet::multinom()` of three or more levels is
#' read as one column per level. Those are the shapes `predict(fit, type =
#' "response")` and `fitted(fit)` give, and trimming a fit trims exactly what
#' trimming those values would.
#'
#' The methods that read an exposure (`"pref"`, `"cr"`, and every method on the
#' categorical route) take it from the model when `.exposure` is not supplied,
#' announcing the variable they read; `options(propensity.quiet = TRUE)`
#' silences the announcement. An `.exposure` you supply is used instead, and a
#' categorical model's columns are matched to its levels by name, so an exposure
#' whose levels are ordered differently is still trimmed against the right
#' column.
#'
#' ## Dose models
#'
#' For a continuous exposure, `.propensity` can be the model of the exposure's
#' conditional mean: a [stats::lm()], a `glm` of the `gaussian()` family (or
#' another family whose variance is constant), a `MASS::rlm()`, or an
#' `mgcv::gam()`, the models [wt_ate()] reads a dose from. A `glm` whose spread
#' changes with its mean, such as `poisson()`, is refused with an error of
#' class `propensity_model_family_error`. Only `"density"` and `"resid"` apply.
#'
#' Both read the fitted conditional mean `mu`, the spread `sigma` of the
#' residuals under the family in `.density`, the standardized residual
#' `z = (a - mu) / sigma`, and the conditional density `f = g(z) / sigma`,
#' where `g` is the family's standardized density. `"density"` keeps the units
#' whose `f` is at least its `lower` quantile. `"resid"` keeps those whose
#' `|z|` is at most `upper`; for the symmetric unimodal families it accepts,
#' that is the floor `g(upper) / sigma` on `f`, which is the threshold it
#' records.
#'
#' The retained values are the conditional means, and the trimmed ones are
#' `NA`. A unit with a missing mean or dose takes no part in the trim. The
#' population a trimmed analysis describes is the units whose observed dose
#' was plausible under the model, not the full sample. `.focal_level`,
#' `.reference_level`, and their deprecated forms are refused with an error of
#' class `propensity_focal_level_error`, since a dose has no levels.
#'
#' ## Object behavior
#'
#' Arithmetic operations on `ps_trim` objects return plain numeric vectors,
#' since transformed propensity scores (e.g., `1/ps`) are no longer propensity
#' scores. Trimmed values propagate as `NA` in calculations; use `na.rm = TRUE`
#' where appropriate.
#'
#' When combining `ps_trim` objects with [c()], trimming parameters must match.
#' Mismatched parameters trigger a warning and return a numeric vector.
#'
#' Use [ps_trim_meta()] to inspect the trimming metadata, including the method,
#' cutoffs, and which observations were retained or trimmed.
#'
#' ## Missing values
#'
#' A propensity score that arrives missing is not one this function removed, so
#' it takes no part in the trimming record: it joins neither the retained nor
#' the trimmed positions, and [is_unit_trimmed()] reports `FALSE` for it. The
#' value propagates as `NA` into the result. For a matrix of categorical
#' propensity scores, a row with a missing cell comes back exactly as it
#' arrived, since there is no complete probability vector to place against a
#' threshold.
#'
#' The methods that read a cutoff off the scores read it off the scores they
#' have. `"adaptive"`, `"pctl"`, `"cr"`, and `"optimal"` all work their cutoffs
#' out from the complete scores or rows, so the cutoff is the one the same call
#' would produce with the missing observations dropped.
#'
#' `"pref"` centers its preference scores on the proportion exposed across the
#' whole sample, which is a fact about the exposure rather than about the
#' scores, so a unit whose propensity score is missing still counts toward it
#' and the cutoff is not the one the shorter call would produce. That unit's own
#' preference score is missing, which leaves it outside both the retained and
#' the trimmed positions like any other missing score.
#'
#' A missing exposure is a different matter, and `"pref"` and `"cr"` refuse one
#' with an error of class `propensity_missing_value_error`. Their cutoffs come
#' from the exposure groups, and a unit that belongs to neither leaves them
#' undefined. Remove or impute the missing exposure values first.
#'
#' ## The trimming record
#'
#' A `ps_trim` records which units were trimmed as positions among the
#' observations it was written for, along with how many observations that was.
#' Operations that hand this package the subscript re-index those positions onto
#' the result: subsetting with `[`, [sort()], [rep()], and [na.omit()] all
#' return a record written for what they return, and a subscript naming a
#' position more than once reports that unit at every place it now holds.
#' [unique()] does the same when it can, as described below.
#'
#' Operations that change how many observations there are without supplying a
#' subscript cannot re-index the record, and it is dropped rather than worked
#' out from the values, since reading membership back from the `NA` pattern
#' would report a propensity score that arrived missing as one this package
#' removed. [vctrs::vec_slice()], which is how filtering, joining, and grouped
#' verbs in dplyr reach a column, is the usual route, and combining two or more
#' vectors with [c()] is another, because concatenation appends one set of
#' observations to another. The record is dropped without comment on every
#' route: most of these length changes build vectors the caller never holds,
#' such as the pieces a grouped verb slices a column into or the rows a tibble
#' slices off to print, so a warning would mostly describe something that is not
#' the result. The values, the class, and the method and its cutoffs are
#' untouched.
#'
#' [unique()] keeps one element for each distinct value, and that element
#' stands for every unit holding the value. The record is re-indexed onto the
#' result when all of those units share one status, and dropped otherwise: a
#' trimmed score and one that arrived missing are both `NA`, so a vector holding
#' both returns a single `NA` that neither status describes.
#'
#' A record can also outlive the observations it describes, because it travels
#' by routes vctrs does not see: growing a `ps_trim` by subassignment carries it
#' across the length change. [is_unit_trimmed()] and [ps_refit()] therefore
#' check that the record covers the object they are given and raise an error of
#' class `propensity_missing_meta_error` when it does not, rather than name
#' trimmed units at stale positions.
#'
#' That check compares how many observations the record was written for against
#' how many the object holds, which a reordering does not change. An operation
#' that reorders the observations through vctrs, rather than through `[`,
#' therefore keeps a record written for the order they used to be in:
#' `vctrs::vec_slice(x, 5:1)` and `dplyr::arrange()` both return the values in
#' a new order under positions still naming the old one. [is_unit_trimmed()]
#' answers from those positions and names the wrong units, and [ps_refit()]
#' refits on the wrong rows. Subsetting with `[` is handed the subscript and
#' re-indexes, so reorder with `[`, or put the propensity scores in the order
#' you want before trimming them.
#'
#' Casting a numeric vector into a `ps_trim` with [vctrs::vec_cast()] is a type
#' operation and not a trimming. The result is described by the method and
#' cutoffs of the target and records that none of the arriving values was
#' trimmed, so it can hold scores outside those cutoffs, including 0 and 1.
#' Call `ps_trim()` on the scores to trim them.
#'
#' @return A **`ps_trim`** object (a numeric vector with class `"ps_trim"`, or a
#'   matrix with class `"ps_trim_matrix"`). Trimmed observations are `NA`.
#'   Metadata is stored in the `"ps_trim_meta"` attribute and can be accessed
#'   with [ps_trim_meta()]. Key fields include:
#'
#'   * `method`: the trimming method used
#'   * `keep_idx`: integer indices of retained observations
#'   * `trimmed_idx`: integer indices of trimmed (NA) observations
#'   * `n_obs`: the number of observations those indices describe
#'   * Method-specific fields such as `cutoff` (adaptive), `q_lower`/`q_upper`
#'     (pctl), `cr_lower`/`cr_upper` (cr), `delta` (categorical ps),
#'     or `lambda` (optimal)
#'   * `focal_inverted` (vector scores only): `TRUE` when the scores are one
#'     minus the probability a fitted model reports, because the model was
#'     trimmed with its first level named as focal, and `FALSE` otherwise,
#'     including for every vector of scores supplied directly
#'
#'   A trim of a dose model records `method`, `lower` (`"density"`, else
#'   `NULL`), `upper` (`"resid"`, else `NULL`), `threshold` (the realized floor
#'   on the conditional density), `sigma` (the spread it was read at),
#'   `sigma_kind` (`"pooled"`, `"mle"`, or `"supplied"`), `density` (the
#'   density specification), `keep_idx`, `trimmed_idx`, and `n_obs`.
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#' with limited overlap in estimation of average treatment effects.
#' *Biometrika*, 96(1), 187--199.
#'
#' Walker, A. M., Patrick, A. R., Lauer, M. S., et al. (2013). A tool for
#' assessing the feasibility of comparative effectiveness research.
#' *Comparative Effectiveness Research*, 3, 11--20.
#'
#' Yang, S., Imbens, G. W., Cui, Z., Faries, D. E., & Kadziola, Z. (2016).
#' Propensity score matching and subclassification in observational studies
#' with multi-level treatments. *Biometrics*, 72(4), 1055--1065.
#'
#' @seealso [ps_trunc()] for bounding (winsorizing) instead of discarding,
#'   [ps_refit()] to re-estimate propensity scores after trimming,
#'   [ps_calibrate()] for calibration-based adjustment,
#'   [ps_trim_meta()] to inspect trimming metadata,
#'   [is_ps_trimmed()] and [is_unit_trimmed()] for logical queries.
#'
#' @examples
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.3 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # Fixed threshold trimming (default)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#' trimmed
#'
#' # How many observations were trimmed?
#' sum(is_unit_trimmed(trimmed))
#'
#' # Data-driven adaptive trimming
#' ps_trim(ps, method = "adaptive")
#'
#' # Quantile-based trimming at 5th and 95th percentiles
#' ps_trim(ps, method = "pctl")
#'
#' # Refit after trimming, then compute weights
#' trimmed <- ps_trim(ps, method = "adaptive")
#' refitted <- ps_refit(trimmed, fit)
#' wt_ate(refitted, .exposure = z)
#'
#' # Trim the scores a fitted model reports, reading the exposure off the model
#' ps_trim(fit, method = "cr")
#'
#' # Trim a dose model on the scale of its conditional density
#' dose <- 1 + 0.5 * x + rt(n, df = 3)
#' dose_fit <- lm(dose ~ x)
#' ps_trim(dose_fit, method = "density", lower = 0.05)
#' ps_trim(dose_fit, method = "resid", upper = 3, .density = dens_t(4))
#'
#' if (rlang::is_installed("nnet")) {
#'   trt <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#'   multinomial_fit <- nnet::multinom(trt ~ x, trace = FALSE)
#'   ps_trim(multinomial_fit, method = "optimal")
#' }
#'
#' @export
ps_trim <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated()
) {
  .propensity <- handle_propensity_deprecation(
    rlang::maybe_missing(.propensity),
    ps,
    "ps_trim"
  )

  UseMethod("ps_trim", .propensity)
}

#' @export
ps_trim.default <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  # Whether `.propensity` is one minus the probability a fitted model reports,
  # which the model route decides when the caller names the model's first
  # level as focal. Supplied scores are taken as given.
  focal_inverted = FALSE,
  # Two frames arrive here because two condition systems read them.
  # `user_env` is the frame lifecycle reports a deprecation from, which decides
  # whether the reader is told to change their own call or to report an issue.
  # `call` is the frame every other condition is attributed to, which rlang
  # reads to name the function in the report. A route that reaches this method
  # by a call rather than by dispatch has to supply both.
  user_env = rlang::caller_env(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  # The dots reach every route but are read by none of them, so a misspelled
  # bound such as `lowr` would trim at the default the caller believed they had
  # replaced. The model methods forward their dots here, so the refusal covers
  # a fit as well as a vector.
  rlang::check_dots_empty(call = call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_ps_method(.propensity, call = call)

  # A calibrated score is a propensity score, so this reads the scores it holds
  # the way the weight functions do rather than refusing a vector of numbers
  # for the class attached to it. What comes back is a trimmed score rather
  # than a calibrated one, so the class is dropped and what it recorded is kept
  # in the trim metadata, where `is_ps_calibrated()` reads it.
  calibrated <- inherits(.propensity, "ps_calib")
  if (calibrated) {
    .propensity <- as.numeric(.propensity)
  }

  method <- rlang::arg_match(method, error_call = call)

  # Optimal trimming is defined over the rows of a propensity score matrix, so
  # there is nothing for it to do with a vector of scores.
  if (method == "optimal") {
    abort(
      c(
        "Method {.val optimal} is only supported for categorical exposures.",
        i = "Supply the propensity scores as a matrix or data frame with one column per exposure level."
      ),
      error_class = "propensity_wt_not_supported_error",
      call = call
    )
  }

  # The density methods read the conditional density of a continuous exposure,
  # which takes the residuals and the family of the model that fit its mean. A
  # vector of fitted means carries neither. The refusal comes before the range
  # check, so that means outside the unit interval are refused for the method
  # asked of them rather than for their values.
  if (method %in% c("density", "resid")) {
    abort(
      c(
        "Method {.val {method}} cannot trim a vector of values.",
        x = "It reads the conditional density of a continuous exposure, which
             needs the residuals and the family of the model that fit its
             conditional mean, and {.arg .propensity} carries neither.",
        i = "Supply the fitted model of the exposure itself as
             {.arg .propensity}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  check_score_trim_density(.sigma, .density, call = call)
  check_ps_range(.propensity, call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trim",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  if (method == "ps") {
    if (is.null(lower)) {
      lower <- 0.1
    }
    if (is.null(upper)) {
      upper <- 0.9
    }
    check_lower_upper(lower, upper, call = call)
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'adaptive'}, {.code lower} and {.code upper} are ignored.",
        call = call
      )
    }
  } else if (method == "pctl") {
    if (is.null(lower)) {
      lower <- 0.05
    }
    if (is.null(upper)) {
      upper <- 0.95
    }
    check_quantile_probs(lower, upper, call = call)
    check_lower_upper(lower, upper, call = call)
  } else if (method == "pref") {
    if (is.null(lower)) {
      lower <- 0.3
    }
    if (is.null(upper)) {
      upper <- 0.7
    }
    check_lower_upper(lower, upper, call = call)
  } else if (method == "cr") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'cr'}, {.code lower} and {.code upper} are ignored.",
        call = call
      )
    }
  }

  n <- length(.propensity)

  # A bound the method never reads describes nothing about the trimming it did,
  # and recording it would leave two trimmings that ran the same rule to the same
  # cutoff describing themselves differently, which is enough for a combine to
  # report different parameters and fall back to numeric.
  meta_list <- if (method %in% c("adaptive", "cr")) {
    list(method = method)
  } else {
    list(method = method, lower = lower, upper = upper)
  }

  if (calibrated) {
    meta_list$calibrated <- TRUE
  }

  # A refit predicts the probability the fitted model reports, and this says
  # whether the retained scores are that probability or its complement.
  meta_list$focal_inverted <- focal_inverted

  # A score that arrived missing is not one this function can place against a
  # cutoff, so it takes no part in working the cutoff out and no part in the
  # record. Every rule below compares scores with `which()`, which leaves a
  # missing comparison out of the retained positions on its own; the trimmed
  # positions are then everything else that was observed.
  observed <- !is.na(.propensity)

  # Decide which indices are kept
  if (method == "ps") {
    keep_idx <- which(.propensity >= lower & .propensity <= upper)
  } else if (method == "adaptive") {
    sum_wt <- 1 / (.propensity[observed] * (1 - .propensity[observed]))
    k <- 2 * mean(sum_wt) - max(sum_wt)

    if (k >= 0) {
      cutoff <- 0
    } else {
      trim_fun <- function(x) {
        sum_wt_trim <- sum_wt[sum_wt <= x]
        x - 2 * mean(sum_wt_trim)
      }
      rng <- range(sum_wt)
      lambda <- uniroot(trim_fun, lower = rng[1], upper = rng[2])$root
      cutoff <- 0.5 - sqrt(0.25 - 1 / lambda)
    }
    meta_list$cutoff <- cutoff
    keep_idx <- which(pmin(.propensity, 1 - .propensity) > cutoff)
  } else if (method == "pctl") {
    # `quantile()` names its result for the probability it was asked for, which
    # says nothing about the cutoff and reappears wherever the cutoff is printed
    # or compared.
    q_lower <- unname(quantile(.propensity, probs = lower, na.rm = TRUE))
    q_upper <- unname(quantile(.propensity, probs = upper, na.rm = TRUE))
    meta_list$q_lower <- q_lower
    meta_list$q_upper <- q_upper
    keep_idx <- which(.propensity >= q_lower & .propensity <= q_upper)
  } else if (method == "pref") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'pref'}, must supply {.arg .exposure}.",
        error_class = "propensity_missing_arg_error",
        call = call
      )
    }
    check_exposure_complete(.exposure, method, observed = observed, call = call)
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    prop_exposure <- mean(.exposure, na.rm = TRUE)
    pref_score <- plogis(qlogis(.propensity) - qlogis(prop_exposure))
    meta_list$P <- prop_exposure
    keep_idx <- which(pref_score >= lower & pref_score <= upper)
  } else if (method == "cr") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'cr'}, must supply {.arg .exposure}.",
        error_class = "propensity_missing_arg_error",
        call = call
      )
    }
    check_exposure_complete(.exposure, method, observed = observed, call = call)
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    ps_treat <- .propensity[observed & .exposure == 1]
    ps_untrt <- .propensity[observed & .exposure == 0]
    check_cr_groups_observed(ps_treat, ps_untrt, call = call)
    cr_lower <- min(ps_treat)
    cr_upper <- max(ps_untrt)
    meta_list$cr_lower <- cr_lower
    meta_list$cr_upper <- cr_upper

    keep_idx <- which(.propensity >= cr_lower & .propensity <= cr_upper)
  }

  trimmed_idx <- setdiff(seq_len(n), c(keep_idx, which(!observed)))

  # Replace trimmed entries with NA
  ps_na <- .propensity
  ps_na[trimmed_idx] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx,
        n_obs = n
      )
    )
  )
}

#' @export
ps_trim.matrix <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  # The generic offers every method, so match against that full set and then
  # reject the ones the categorical path does not define.
  method <- rlang::arg_match(
    method,
    values = c(
      "ps",
      "adaptive",
      "pctl",
      "pref",
      "cr",
      "optimal",
      "density",
      "resid"
    ),
    error_call = call
  )
  if (!method %in% c("ps", "optimal")) {
    abort(
      c(
        "Method {.val {method}} is not supported for categorical exposures.",
        i = "Use {.val ps} or {.val optimal}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  check_score_trim_density(.sigma, .density, call = call)

  check_no_focal_levels(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    call = call
  )

  # Validate exposure for categorical
  if (is.null(.exposure)) {
    abort(
      "`.exposure` must be provided for categorical propensity score trimming.",
      error_class = "propensity_missing_arg_error",
      call = call
    )
  }

  # Transform to factor and validate
  .exposure <- transform_exposure_categorical(.exposure, call = call)

  # Validate matrix
  .propensity <- check_ps_matrix(
    .propensity,
    .exposure,
    call = call
  )

  n <- nrow(.propensity)
  k <- ncol(.propensity)

  # A row with a missing score has no complete probability vector to place
  # against a threshold, so it takes no part in working the threshold out and no
  # part in the record. Both rules below compare rows with `which()`, which
  # leaves a missing comparison out of the retained positions on its own, and the
  # group-preservation reset falls back to the complete rows for the same reason.
  incomplete_rows <- which(apply(.propensity, 1, anyNA))
  complete_rows <- setdiff(seq_len(n), incomplete_rows)

  # Initialize metadata
  meta_list <- list(method = method, is_matrix = TRUE)

  if (method == "ps") {
    # Symmetric trimming
    if (is.null(lower)) {
      lower <- 0.1
    }
    delta <- lower # Use lower as delta for consistency

    # A threshold at or above 1/k cannot be met by every column of a row that
    # sums to one, so there is no trimming rule left to apply. Both numbers are
    # facts about this call: the threshold supplied and the limit the width of
    # the matrix imposes on it.
    if (delta >= 1 / k) {
      limit <- format(1 / k, digits = 7)
      abort(
        c(
          "The trimming threshold must fall below 1/k, for k columns of
           propensity scores.",
          x = "{.arg lower} is {.val {delta}}, and 1/k is {limit} for the {k}
               column{?s} the scores hold.",
          i = "No row summing to one can hold every score above 1/k, so a
               threshold there leaves no rule to apply."
        ),
        error_class = "propensity_range_error",
        call = call
      )
    }

    # Apply symmetric trimming rule: keep if min(propensity scores) > delta
    keep_idx <- which(apply(.propensity, 1, function(x) min(x) > delta))

    # Check if all treatment groups are preserved
    if (length(unique(.exposure[keep_idx])) < k) {
      warn(
        "One or more groups removed after trimming; returning original data",
        warning_class = "propensity_no_data_warning",
        call = call
      )
      keep_idx <- complete_rows
    }

    meta_list$delta <- delta
  } else {
    # optimal
    # Multi-category optimal trimming (Yang et al., 2016)
    # Calculate sum of inverse propensity scores
    sum_inv_ps <- rowSums(1 / .propensity)
    sum_inv_complete <- sum_inv_ps[complete_rows]

    # Define trimming function
    trim_fun <- function(x) {
      sum_trim <- sum_inv_complete[sum_inv_complete <= x]
      if (length(sum_trim) == 0) {
        return(x)
      }
      x - 2 * mean(sum_trim) / mean(sum_inv_complete <= x)
    }

    # Check if trimming is needed
    if (trim_fun(max(sum_inv_complete)) < 0) {
      # No valid solution, use maximum + 1
      lambda <- max(sum_inv_complete) + 1
      keep_idx <- complete_rows # Keep all
    } else {
      # Find optimal lambda
      result <- tryCatch(
        {
          uniroot(
            trim_fun,
            lower = min(sum_inv_complete),
            upper = max(sum_inv_complete)
          )$root
        },
        error = function(e) {
          warn(
            "Could not find optimal trimming threshold; using no trimming",
            warning_class = "propensity_convergence_warning",
            call = call
          )
          NULL
        }
      )

      if (!is.null(result)) {
        lambda <- result
        keep_idx <- which(sum_inv_ps <= lambda)

        # Check if all treatment groups are preserved
        if (length(unique(.exposure[keep_idx])) < k) {
          warn(
            "One or more groups removed after trimming; returning original data",
            warning_class = "propensity_no_data_warning",
            call = call
          )
          keep_idx <- complete_rows
          lambda <- NULL
        }
      } else {
        keep_idx <- complete_rows
        lambda <- NULL
      }
    }

    meta_list$lambda <- lambda
  }

  trimmed_idx <- setdiff(seq_len(n), c(keep_idx, incomplete_rows))

  # Replace trimmed entries with NA
  ps_na <- .propensity
  ps_na[trimmed_idx, ] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx,
        n_obs = n
      )
    )
  )
}

#' @export
ps_trim.data.frame <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_predicted_class_column(.propensity, call = call)

  # For categorical exposures, convert to matrix and call matrix method
  if (!is.null(.exposure)) {
    exposure_type <- causalgenerics::detect_exposure_type(
      .exposure,
      announce = !be_quiet(),
      call = call
    )
    if (exposure_type == "categorical") {
      ps_matrix <- as.matrix(.propensity)
      # The focal arguments travel on so that the matrix method refuses them.
      # Dropping them here would leave the caller believing the trimming honored
      # a coding the categorical path never reads.
      #
      # The matrix method is reached here by a call rather than by dispatch, so
      # it is handed a frame to report against. Left to its own, a refusal on
      # this route would name the method the caller never wrote.
      return(ps_trim.matrix(
        .propensity = ps_matrix,
        method = method,
        lower = lower,
        upper = upper,
        .exposure = .exposure,
        .focal_level = .focal_level,
        .reference_level = .reference_level,
        ...,
        .sigma = .sigma,
        .density = .density,
        .treated = .treated,
        .untreated = .untreated,
        call = call
      ))
    }
  }

  ps_vec <- binary_ps_column(.propensity, "ps_trim")

  # The default method is reached here by a call rather than by dispatch, so it
  # is handed a frame to report against. Left to its own, a refusal on this
  # route would name the method the caller never wrote.
  ps_trim.default(
    .propensity = ps_vec,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .sigma = .sigma,
    .density = .density,
    .treated = .treated,
    .untreated = .untreated,
    user_env = rlang::caller_env(),
    call = call
  )
}

# The fitted models trimming reads. A `glm` of a binomial family fits the
# probability of a binary exposure and a `multinom` fits a probability for every
# level, so both are trimmed on the scale of that probability. A `lm`, and a
# `glm` of any other family, fits the conditional mean of a dose, which is
# trimmed on the scale of its conditional density.
#' @export
ps_trim.glm <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  trim <- if (is_binomial_family(.propensity[["family"]])) {
    ps_trim_from_model
  } else {
    ps_trim_dose
  }

  trim(
    .propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .sigma = .sigma,
    .density = .density,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# A `lm` has no family, and its fitted values are the conditional mean of a
# dose. `MASS::rlm()` fits and other subclasses arrive here by inheritance.
#' @export
ps_trim.lm <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  ps_trim_dose(
    .propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .sigma = .sigma,
    .density = .density,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# A `multinom` is `c("multinom", "nnet")` and inherits from neither `glm` nor
# `lm`, so it reaches a method of its own. What it was fit to is checked before
# anything reads it, so that a fit the route cannot read is reported as such
# rather than as whatever its fitted values are made of.
#' @export
ps_trim.multinom <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_multinom_response(.propensity, call = call)

  ps_trim_from_model(
    .propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .sigma = .sigma,
    .density = .density,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# The model methods of the score route differ only in the class they are
# registered for and in what that class needs checked before it is read, so
# both route through here. The method is resolved first because it decides
# whether the trimming reads an exposure at all, and a score method the vector
# route does not define is left to the route that does not define it, so that
# the refusal is the same one a vector of scores gets.
#
# The vector and matrix methods are reached by a call rather than by dispatch,
# so they are handed the frame the route was entered on. Left to their own, a
# refusal here would name a method the caller never wrote.
ps_trim_from_model <- function(
  model,
  method,
  lower,
  upper,
  .exposure,
  .focal_level,
  .reference_level,
  ...,
  .sigma,
  .density,
  .treated,
  .untreated,
  call,
  user_env
) {
  method <- rlang::arg_match(
    method,
    values = ps_trim_methods,
    error_call = call
  )

  # Checked before anything is read off the fit, so that a refusal is never
  # preceded by an announcement about an exposure nothing goes on to use.
  if (method %in% ps_trim_dose_methods) {
    abort(
      c(
        "Method {.val {method}} cannot trim a model of a probability.",
        x = "{.arg .propensity} is {.cls {class(model)[[1]]}}, whose fitted
             values are probabilities of the exposure, and a model of a
             probability has no conditional density to trim.",
        i = "Trim its propensity scores with a score method such as
             {.val ps}, or supply the model of a continuous exposure's
             conditional mean to trim with {.val {method}}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }
  check_score_trim_density(.sigma, .density, call = call)

  args <- prepare_model_ps(
    model,
    .exposure,
    # Trimming a matrix of scores reads the exposure whatever the method, both
    # to hold the columns against its levels and to check that no group is
    # trimmed away, so a fit that reports a column for every level is read for
    # its exposure too.
    needs_exposure = method %in% c("pref", "cr") || model_fits_levels(model),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "ps_trim",
    call = call,
    user_env = user_env
  )

  if (identical(args$exposure_type, "categorical")) {
    return(ps_trim.matrix(
      .propensity = args$propensity,
      method = method,
      lower = lower,
      upper = upper,
      .exposure = args$exposure,
      .focal_level = args$focal_level,
      .reference_level = args$reference_level,
      ...,
      call = call
    ))
  }

  ps_trim.default(
    .propensity = args$propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = args$exposure,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    ...,
    focal_inverted = args$focal_inverted,
    user_env = user_env,
    call = call
  )
}

ps_trim_score_methods <- c("ps", "adaptive", "pctl", "pref", "cr", "optimal")
ps_trim_dose_methods <- c("density", "resid")
ps_trim_methods <- c(ps_trim_score_methods, ps_trim_dose_methods)

# The score routes read no density, so a family or a spread handed to them
# would be an instruction silently ignored. The default family is accepted, so
# that a caller who writes out the family they were getting anyway is not
# refused for saying so, which is how the weight functions treat a binary
# exposure.
check_score_trim_density <- function(
  .sigma,
  .density,
  call = rlang::caller_env()
) {
  if (!is.null(.sigma)) {
    abort(
      c(
        "{.arg .sigma} applies only to trimming a dose model.",
        x = "The trimming method bounds propensity scores, which have no
             spread to read.",
        i = "Leave {.arg .sigma} unset, or supply the model of a continuous
             exposure with {.code method = \"density\"} or
             {.code method = \"resid\"}."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  spec <- as_density_spec(.density, arg = ".density", call = call)
  if (identical(spec$family, "normal")) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.arg .density} applies only to trimming a dose model.",
      x = "The trimming method bounds propensity scores, which have no
           conditional density to read.",
      i = "Leave {.arg .density} unset, or supply the model of a continuous
           exposure with {.code method = \"density\"} or
           {.code method = \"resid\"}."
    ),
    error_class = "propensity_density_error",
    call = call
  )
}

# Trimming a dose model on the scale of its conditional density. The fitted
# values are the conditional mean of the dose, the residuals are spread by the
# estimator the density family asks for, as the weights spread them, and each
# unit's conditional density at its own dose decides whether it is kept.
#
# Every argument is checked before the exposure is read, so a refusal is never
# preceded by an announcement about an exposure nothing goes on to use.
ps_trim_dose <- function(
  model,
  method,
  lower,
  upper,
  .exposure,
  .focal_level,
  .reference_level,
  ...,
  .sigma,
  .density,
  .treated,
  .untreated,
  call,
  user_env
) {
  rlang::check_dots_empty(call = call)
  method <- rlang::arg_match(
    method,
    values = ps_trim_methods,
    error_call = call
  )

  check_continuous_model_family(
    model,
    problem = "Trimming a dose model needs a model of its conditional mean
               with a single spread.",
    remedy = "Fit the dose model with {.fun gaussian}, {.fun lm},
              {.fun mgcv::gam}, or {.fun MASS::rlm}.",
    call = call
  )
  check_continuous_model_response(
    model,
    problem = "Trimming a dose model needs a model of one conditional mean for
               each unit.",
    remedy = "Fit the dose on its own and trim that model.",
    call = call
  )

  if (!method %in% ps_trim_dose_methods) {
    abort(
      c(
        "Method {.val {method}} cannot trim a dose model.",
        x = "{.arg .propensity} is {.cls {class(model)[[1]]}}, a model of a
             continuous exposure's conditional mean, which has no propensity
             score in (0, 1) for {.val {method}} to bound.",
        i = "Trim it on the scale of its conditional density with
             {.code method = \"density\"} or {.code method = \"resid\"}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trim",
    user_env = user_env
  )
  check_dose_no_levels(
    focal_params$.focal_level,
    focal_params$.reference_level,
    call = call
  )

  density <- as_density_spec(.density, arg = ".density", call = call)
  check_dose_trim_sigma(.sigma, call = call)
  check_sigma_method(.sigma, density, call = call)

  if (method == "density") {
    if (!is.null(upper)) {
      abort(
        c(
          "{.arg upper} is not used by {.code method = \"density\"}.",
          x = "The method sets aside the units whose conditional density falls
               below a quantile floor, and a large conditional density is a
               plausible dose, so there is no upper tail to trim.",
          i = "Drop {.arg upper}, and set the floor with {.arg lower}."
        ),
        error_class = "propensity_unsupported_arg_error",
        call = call
      )
    }
    if (is.null(lower)) {
      lower <- 0.01
    }
    check_density_trim_lower(lower, call = call)
  } else {
    if (!is.null(lower)) {
      abort(
        c(
          "{.arg lower} is not used by {.code method = \"resid\"}.",
          x = "The method bounds the absolute standardized residual, which has
               no lower end to trim.",
          i = "Drop {.arg lower}, and set the bound with {.arg upper}."
        ),
        error_class = "propensity_unsupported_arg_error",
        call = call
      )
    }
    if (is.null(upper)) {
      abort(
        c(
          "{.code method = \"resid\"} needs {.arg upper}.",
          x = "The bound on the absolute standardized residual has no
               default.",
          i = "Supply the number of spreads from the predicted dose beyond
               which a unit is trimmed, such as {.code upper = 3}."
        ),
        error_class = "propensity_missing_arg_error",
        call = call
      )
    }
    check_resid_trim_upper(upper, call = call)
    check_resid_trim_density(density, call = call)
  }

  mu <- as.numeric(extract_continuous_ps(model, call = call)$mu)
  exposure <- extract_exposure_from_model(model, .exposure)
  check_dose_trim_exposure(exposure, length(mu), call = call)

  n <- length(mu)
  sigma <- continuous_sigma(exposure, mu, .sigma, density, call = call)
  z <- (exposure - mu) / sigma
  f <- density_eval_present(density, z, call = call) / sigma

  # A unit with a missing mean or dose has no density to place against the
  # floor, so it takes no part in working the floor out and no part in the
  # record, as a missing score takes none in a trim of scores.
  observed <- !is.na(f)

  if (method == "density") {
    # `quantile()` names its result for the probability it was asked for, which
    # says nothing about the floor.
    threshold <- unname(stats::quantile(f, probs = lower, na.rm = TRUE))
    keep_idx <- which(f >= threshold)
  } else {
    # For a symmetric unimodal density the bound on `|z|` is a floor on the
    # density written in residual units, and the floor is what is recorded, so
    # that both methods describe the same kind of cut.
    threshold <- density_eval(density, upper, call = call) / sigma
    keep_idx <- which(abs(z) <= upper)
  }

  trimmed_idx <- setdiff(which(observed), keep_idx)

  values <- mu
  values[trimmed_idx] <- NA_real_

  # Every field is named whichever method made the record, and the bound a
  # method does not read is `NULL`, so that two records compare field by field.
  new_trimmed_ps(
    x = values,
    ps_trim_meta = list(
      method = method,
      lower = if (method == "density") lower,
      upper = if (method == "resid") upper,
      threshold = threshold,
      sigma = sigma,
      sigma_kind = density_sigma_source(.sigma, density),
      density = density,
      keep_idx = keep_idx,
      trimmed_idx = trimmed_idx,
      n_obs = n
    )
  )
}

# A dose has no levels, so an argument naming one describes a coding the dose
# route never reads.
check_dose_no_levels <- function(
  .focal_level,
  .reference_level,
  call = rlang::caller_env()
) {
  supplied <- c(
    ".focal_level" = !is.null(.focal_level),
    ".reference_level" = !is.null(.reference_level)
  )

  if (!any(supplied)) {
    return(invisible(NULL))
  }

  named <- names(supplied)[supplied]
  abort(
    c(
      "{.arg {named}} cannot be used when trimming a dose model.",
      x = "A continuous exposure has no levels, so none is focal and none is
           reference.",
      i = "Drop {.arg {named}}."
    ),
    error_class = "propensity_focal_level_error",
    call = call
  )
}

# A dose trim is recorded with the single spread it was read at, which a refit
# can keep. A spread for each unit describes a density no refit could
# reproduce, so it is refused here rather than accepted as the weights accept
# it.
check_dose_trim_sigma <- function(.sigma, call = rlang::caller_env()) {
  if (is.null(.sigma)) {
    return(invisible(NULL))
  }

  if (is.numeric(.sigma) && is.null(dim(.sigma)) && length(.sigma) != 1) {
    abort(
      c(
        "{.arg .sigma} must be a single spread when trimming a dose model.",
        x = "It holds {length(.sigma)} value{?s}.",
        i = "A trim is recorded at one spread, which a refit of the dose model
             can keep; a spread for each unit cannot be refit. Supply one
             value, or leave {.arg .sigma} unset to estimate it."
      ),
      error_class = "propensity_sigma_error",
      call = call
    )
  }

  check_sigma(.sigma, "continuous", n = 1, call = call)
}

check_density_trim_lower <- function(lower, call = rlang::caller_env()) {
  if (
    is.numeric(lower) &&
      length(lower) == 1 &&
      !is.na(lower) &&
      lower > 0 &&
      lower < 0.5
  ) {
    return(invisible(NULL))
  }

  abort(
    c(
      "For {.code method = \"density\"}, {.arg lower} must be a single
       probability between 0 and 0.5.",
      x = "{.arg lower} is {.val {lower}}.",
      i = "{.arg lower} is the quantile of the conditional density below which
           units are trimmed, such as the default {.code lower = 0.01}."
    ),
    error_class = "propensity_range_error",
    call = call
  )
}

check_resid_trim_upper <- function(upper, call = rlang::caller_env()) {
  if (
    is.numeric(upper) &&
      length(upper) == 1 &&
      is.finite(upper) &&
      upper > 0
  ) {
    return(invisible(NULL))
  }

  abort(
    c(
      "For {.code method = \"resid\"}, {.arg upper} must be a single positive,
       finite number.",
      x = "{.arg upper} is {.val {upper}}.",
      i = "{.arg upper} is the number of spreads from the predicted dose
           beyond which a unit is trimmed, such as {.code upper = 3}."
    ),
    error_class = "propensity_range_error",
    call = call
  )
}

# A bound on the absolute residual is a floor on the density only for a
# symmetric unimodal family. A kernel estimate or a function of the caller's
# need be neither, and the bound and the floor would then keep different units.
check_resid_trim_density <- function(density, call = rlang::caller_env()) {
  if (density$family %in% c("normal", "t", "laplace")) {
    return(invisible(NULL))
  }

  family <- density$family
  abort(
    c(
      "{.code method = \"resid\"} needs a symmetric unimodal density.",
      x = "{.arg .density} is a {.val {family}} density, which need not be
           symmetric or unimodal, so a bound on the absolute residual is not a
           floor on it.",
      i = "Use {.code method = \"density\"}, which trims on the density itself
           under any family, or a {.fun dens_normal}, {.fun dens_t}, or
           {.fun dens_laplace} density."
    ),
    error_class = "propensity_density_error",
    call = call
  )
}

# The dose is read against the model's fitted means one unit at a time, so it
# has to be a number for each of them.
check_dose_trim_exposure <- function(exposure, n, call = rlang::caller_env()) {
  if (!is.numeric(exposure) || !is.null(dim(exposure))) {
    abort(
      c(
        "{.arg .exposure} must be a numeric vector when trimming a dose model.",
        x = "It is {.obj_type_friendly {exposure}}.",
        i = "Supply the dose the model was fit to, one value per unit."
      ),
      error_class = "propensity_type_error",
      call = call
    )
  }

  if (length(exposure) != n) {
    abort(
      c(
        "{.arg .exposure} must have one value for each fitted mean.",
        x = "It has {length(exposure)} value{?s}, and the model reports
             {n} fitted mean{?s}.",
        i = "Supply the dose for the rows the model was fit to, or refit it
             with {.code na.action = na.exclude}, which reports a mean for
             every row of its data."
      ),
      error_class = "propensity_length_error",
      call = call
    )
  }

  invisible(NULL)
}

#' @export
ps_trim.ps_trim <- function(
  .propensity,
  method = c(
    "ps",
    "adaptive",
    "pctl",
    "pref",
    "cr",
    "optimal",
    "density",
    "resid"
  ),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .sigma = NULL,
  .density = "normal",
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated()
) {
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  warn(
    "Propensity scores have already been trimmed. Returning original object.",
    warning_class = "propensity_already_modified_warning"
  )
  .propensity
}

new_trimmed_ps <- function(x, ps_trim_meta = list()) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_trim_meta = ps_trim_meta,
      class = c("ps_trim_matrix", "ps_trim", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_trim_meta = ps_trim_meta,
      class = "ps_trim",
      inherit_base_type = TRUE
    )
  }
}

# The positional half of a trimming record. The rest of it, the method and its
# cutoffs, describes the trimming rather than the units and means the same thing
# at any length.
reindex_trim_record <- function(meta, i) {
  meta$keep_idx <- reindex_positions(meta$keep_idx, i)
  meta$trimmed_idx <- reindex_positions(meta$trimmed_idx, i)
  meta$n_obs <- length(i)

  meta
}

# A record that no longer describes the observations in front of it is dropped
# rather than guessed at. Nothing in the values says which units a shorter or a
# longer vector once trimmed, and reading membership back from the `NA` pattern
# would report a propensity score that arrived missing as one this package
# removed.
#
# The drop is silent. Most of the length changes that reach it are vectors
# built along the way to a result, such as the pieces a grouped dplyr verb
# slices a column into or the rows a tibble slices off to print, so a warning
# would mostly describe something the caller never holds. A positional query on
# the result refuses to answer from a record that is gone, which is where the
# loss matters.
drop_trim_record <- function(meta) {
  meta[["keep_idx"]] <- NULL
  meta[["trimmed_idx"]] <- NULL
  meta[["n_obs"]] <- NULL

  meta
}

# Each unit's standing in a record that covers it: retained, trimmed, or
# neither, which is a score that arrived missing.
trim_unit_status <- function(meta) {
  status <- rep("missing", meta$n_obs)
  status[meta$keep_idx] <- "kept"
  status[meta$trimmed_idx] <- "trimmed"
  status
}

# A positional query reads its answer out of the record, so a record that does
# not cover the object in front of it has no answer to give: reporting every
# unit as retained would be a wrong answer rather than a missing one.
check_trim_record <- function(meta, n, fn, call = rlang::caller_env()) {
  if (record_covers(meta, n)) {
    return(invisible(meta))
  }

  recorded <- meta$n_obs
  problem <- if (is.null(recorded)) {
    "These propensity scores carry no record of which units were trimmed."
  } else {
    "The record describes {recorded} observation{?s} and these propensity
     scores have {n}, so its positions do not describe them."
  }

  abort(
    c(
      "{.code {fn}} has no usable trimming record for these propensity scores.",
      x = problem,
      i = "An operation that changed the number of observations dropped it.
           Trim the propensity scores you want to work with, or query the
           object the record was written for."
    ),
    error_class = "propensity_missing_meta_error",
    call = call
  )
}

#' Test whether propensity scores have been trimmed
#'
#' @description `is_ps_trimmed()` returns `TRUE` if `x` is a `ps_trim` object
#'   or a `psw` object created from trimmed propensity scores, and `FALSE`
#'   otherwise. This tests whether the *object* carries trimming information, not
#'   which individual units were trimmed; see [is_unit_trimmed()] for that.
#'
#' @param x An object to test.
#' @return A logical scalar (`TRUE` or `FALSE`).
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_unit_trimmed()] to
#'   identify which units were trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_ps_trimmed(trimmed)
#' is_ps_trimmed(ps)
#'
#' @export
is_ps_trimmed <- function(x) {
  UseMethod("is_ps_trimmed")
}

#' @export
is_ps_trimmed.default <- function(x) {
  FALSE
}

#' @export
is_ps_trimmed.ps_trim <- function(x) {
  TRUE
}

#' @export
is_ps_trimmed.ps_trim_matrix <- function(x) {
  TRUE
}

#' Identify which units were trimmed
#'
#' @description `is_unit_trimmed()` returns a logical vector indicating which
#'   observations were removed by trimming. This is a per-unit query, as opposed
#'   to [is_ps_trimmed()], which tests whether the object has been trimmed at
#'   all.
#'
#'   The answer comes from the trimming record, which is written for a fixed set
#'   of observations and can both be lost and outlive them. On a `ps_trim`,
#'   [vctrs::vec_slice()] and [c()] drop it, and subassignment that grows the
#'   vector carries it across a length change; see [ps_trim()] for the whole
#'   contract. On a [psw] vector built from trimmed propensity scores, a subset
#'   drops it, while `model.frame()` re-attaches it to the shortened weights
#'   column of an outcome model fit on these weights.
#'
#'   `is_unit_trimmed()` therefore checks that the record covers the object it
#'   is given, and raises an error of class `propensity_missing_meta_error` when
#'   it does not, or when an object marked as trimmed carries no record at all,
#'   rather than name trimmed units at stale positions. Query the `ps_trim`
#'   object the record was written for instead.
#'
#'   That check counts observations, which a reordering does not change, so it
#'   does not catch one. A `ps_trim` reordered through vctrs rather than through
#'   `[`, by `vctrs::vec_slice(x, 5:1)` or `dplyr::arrange()`, keeps a record
#'   written for the old order, and a `psw` keeps one through any same-length
#'   operation, a reordering included. `is_unit_trimmed()` answers from those
#'   positions and names the wrong units. See [ps_trim()] and [psw] for the
#'   whole contract.
#'
#' @param x A `ps_trim` object created by [ps_trim()], or a [psw] vector built
#'   from one.
#' @return A logical vector the same length as `x`, where `TRUE` marks a
#'   trimmed unit.
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] to
#'   test whether an object has been trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_unit_trimmed(trimmed)
#'
#' # Use to subset data to retained observations
#' kept <- !is_unit_trimmed(trimmed)
#' ps[kept]
#'
#' @export
is_unit_trimmed <- function(x) {
  UseMethod("is_unit_trimmed")
}

#' @export
is_unit_trimmed.default <- function(x) {
  abort(
    "{.code is_unit_trimmed()} not supported for class {.val {class(x)}}",
    error_class = "propensity_method_error"
  )
}

#' @export
is_unit_trimmed.ps_trim <- function(x) {
  # No observations, no answers. A record kept on an empty vector describes
  # observations it does not have, and indexing an empty logical by the
  # positions it names would grow one padded with `NA`.
  if (length(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trim_meta(x)
  check_trim_record(meta, length(x), "is_unit_trimmed()")

  out <- vector("logical", length = length(x))
  out[meta$trimmed_idx] <- TRUE

  out
}

#' @export
is_unit_trimmed.ps_trim_matrix <- function(x) {
  if (nrow(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trim_meta(x)
  check_trim_record(meta, nrow(x), "is_unit_trimmed()")

  out <- vector("logical", length = nrow(x))
  out[meta$trimmed_idx] <- TRUE

  out
}


#' @export
`[.ps_trim_matrix` <- function(x, i, j, ..., drop = TRUE) {
  # Get metadata
  meta <- ps_trim_meta(x)

  # Handle single index (matrix as vector) - bypass ps_trim method
  if (nargs() == 2) {
    return(unclass(x)[i])
  }

  # Handle missing i (all rows) - bypass ps_trim method
  if (missing(i)) {
    return(unclass(x)[, j, ..., drop = drop])
  }

  # For single element extraction, call base method directly to avoid ps_trim method
  # Check if this will result in a single element
  if (!missing(j) && length(i) == 1 && length(j) == 1) {
    return(as.numeric(unclass(x)[i, j, drop = TRUE]))
  }

  # Perform subsetting with base matrix method to avoid calling [.ps_trim
  result <- unclass(x)[i, j, ..., drop = drop]

  # If not a matrix anymore or dropped dimensions, return as-is (no metadata)
  if (!is.matrix(result)) {
    return(result)
  }

  # If result is a single element, return as numeric (no metadata)
  if (length(result) == 1) {
    return(as.numeric(result))
  }

  # The rows the result was built from, as positions in `x`, which is what
  # re-indexing the record takes. These are worked out from the subscript rather
  # than read off the result, so a form they do not reproduce, seen as a count
  # that does not match the rows base returned, degrades the record rather than
  # certifying it against rows it does not describe.
  rows <- subscript_row_positions(i, x)
  new_meta <- if (length(rows) == nrow(result)) {
    carry_trim_record(meta, nrow(x), rows)
  } else {
    drop_trim_record(meta)
  }

  attr(result, "ps_trim_meta") <- new_meta
  class(result) <- c("ps_trim_matrix", "ps_trim", "matrix")
  result
}


# Print methods for ps_trim_matrix

#' @export
print.ps_trim_matrix <- function(x, ..., n = NULL) {
  meta <- ps_trim_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)
  n_trimmed <- length(meta$trimmed_idx)

  # Create header
  cat(sprintf(
    "<ps_trim_matrix[%d x %d]; trimmed %d of %d; method=%s>\n",
    n_rows,
    k,
    n_trimmed,
    n_rows,
    meta$method
  ))

  # Determine how many rows to show
  if (is.null(n)) {
    n_print <- getOption("propensity.print_max", default = 10)
  } else {
    n_print <- n
  }

  # The header summarizes the record, and the scores are what is left to show.
  # The record itself is a set of index vectors as long as the data, so printed
  # after the matrix it is a second and longer object no one asked to see.
  x_print <- unclass(x)
  attr(x_print, "ps_trim_meta") <- NULL

  # Show all rows if n is Inf or very large
  if (is.infinite(n_print) || n_print >= n_rows) {
    print(x_print)
  } else {
    # Show first n_print rows
    n_show <- min(n_print, n_rows)
    print(x_print[seq_len(n_show), , drop = FALSE])

    if (n_rows > n_show) {
      cat("# ... with", n_rows - n_show, "more rows\n")
    }
  }

  invisible(x)
}

# vctrs machinery for ps_trim

#' @export
vec_ptype_abbr.ps_trim <- function(x, ...) "ps_trim"

#' @export
vec_ptype_full.ps_trim <- function(x, ...) {
  meta <- ps_trim_meta(x)

  # Without a record there is no count to report, and reporting none trimmed
  # would say something about the units the record no longer speaks for.
  if (is.null(meta$n_obs)) {
    return("ps_trim; record dropped")
  }

  # A count of trimmed units means something against the number of units there
  # were, and both numbers are read off the record so that they describe the
  # same set of observations.
  paste0("ps_trim; trimmed ", length(meta$trimmed_idx), " of ", meta$n_obs)
}

#' @export
#' @method vec_arith ps_trim
vec_arith.ps_trim <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trim", y)
}

#' @export
#' @method vec_arith.ps_trim default
vec_arith.ps_trim.default <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim ps_trim
vec_arith.ps_trim.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim MISSING
vec_arith.ps_trim.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_trim numeric
vec_arith.ps_trim.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trim
vec_arith.numeric.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim integer
vec_arith.ps_trim.integer <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim list
vec_arith.ps_trim.list <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

# How a trimming is described, as opposed to which units it touched: the
# method, the bounds it was given, and the cutoffs it settled on. A cutoff is
# read off the scores, so two objects that agree on the method and the bounds
# can still have been trimmed at different places. A trim of a dose model is
# also described by the spread and the density family it was read at, and a
# vector of scores by whether they are the complement of what the model reports,
# since scores for different levels are not scores of one trimming.
trim_parameters <- function(meta) {
  fields <- c(
    "method",
    "lower",
    "upper",
    "cutoff",
    "q_lower",
    "q_upper",
    "P",
    "cr_lower",
    "cr_upper",
    "threshold",
    "sigma",
    "sigma_kind",
    "density",
    "focal_inverted"
  )

  params <- rlang::set_names(
    lapply(fields, function(field) meta[[field]]),
    fields
  )
  params["density"] <- list(trim_density_parameters(meta$density))
  # A record written without the field describes scores as the model reports
  # them, which is what `FALSE` says.
  params$focal_inverted <- isTRUE(meta$focal_inverted)

  params
}

# What identifies a density specification, as `density_specs_agree()` reads
# it. Each constructor writes a fresh function, so two specifications built the
# same way hold functions that are not `identical()`; only a function the caller
# wrote is the specification's identity.
trim_density_parameters <- function(density) {
  if (is.null(density)) {
    return(NULL)
  }

  params <- list(
    family = density$family,
    params = density$params,
    sigma_method = density$sigma_method
  )
  if (identical(density$family, "function")) {
    params$fn <- density$fn
  }

  params
}

#' @export
vec_ptype2.ps_trim.ps_trim <- function(x, y, ...) {
  x_meta <- ps_trim_meta(x)
  y_meta <- ps_trim_meta(y)

  # Check if trim parameters match
  if (!identical(trim_parameters(x_meta), trim_parameters(y_meta))) {
    warn_incompatible_metadata(
      x,
      y,
      "different trimming parameters"
    )
    return(double())
  }

  # Check if refit status matches
  if (!identical(x_meta$refit, y_meta$refit)) {
    warn_incompatible_metadata(
      x,
      y,
      "different refit status"
    )
    return(double())
  }

  # The prototype carries the description of the trimming across and holds no
  # observations, so it names no positions: it is shared by inputs whose
  # observations are appended one after another, and the positions either record
  # names would describe units from the other input. Trimming an empty vector
  # again would work the cutoffs out from scores that are not there, and the
  # rules defined against the exposure have none to be handed.
  new_trimmed_ps(double(), ps_trim_meta = drop_trim_record(x_meta))
}
#' @export
vec_ptype2.ps_trim.double <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}
#' @export
vec_ptype2.double.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}

# A cast returns the values it was handed in the type it was handed, and a
# `ps_trim`'s type is the whole description of the trimming. That is the
# comparison `vec_ptype2()` already makes when it refuses to find a common type,
# so the cast makes it too: a cast comparing less than the combine does hands
# `x` back describing itself under the target's name. The positional half of the
# record describes the values arriving rather than the type they arrive in, so
# it is left out of the comparison, which is also what lets a prototype built by
# `drop_trim_record()` be cast to.
#
# `vec_ptype_full()` names none of what is compared, so the two types render
# identically and the refusal would read as a type that cannot be converted to
# itself. What disagrees is named alongside them, the way the combine names it.
#' @export
vec_cast.ps_trim.ps_trim <- function(x, to, ...) {
  x_meta <- ps_trim_meta(x)
  to_meta <- ps_trim_meta(to)

  problem <- if (
    !identical(trim_parameters(x_meta), trim_parameters(to_meta))
  ) {
    "different trimming parameters"
  } else if (!identical(x_meta$refit, to_meta$refit)) {
    "different refit status"
  }

  if (!is.null(problem)) {
    vctrs::stop_incompatible_cast(
      x,
      to,
      x_arg = "",
      to_arg = "",
      details = problem
    )
  }

  x
}

#' @export
vec_cast.double.ps_trim <- function(x, to, ...) {
  # degrade to numeric with NAs
  vec_data(x)
}

# A cast returns the values it was handed in the type it was handed, and the
# trimming is part of that type, so the description comes from the target. The
# values are the ones arriving, so the positions are written for them: none of
# them was removed on the way here.
trim_meta_for_cast <- function(to, x) {
  c(
    drop_trim_record(ps_trim_meta(to)),
    list(
      keep_idx = seq_along(x),
      trimmed_idx = integer(0),
      n_obs = length(x)
    )
  )
}

#' @export
vec_cast.ps_trim.double <- function(x, to, ...) {
  new_trimmed_ps(x, ps_trim_meta = trim_meta_for_cast(to, x))
}

#' @export
vec_ptype2.psw.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_trim"))
  double()
}

#' @export
vec_ptype2.ps_trim.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "psw"))
  double()
}

#' @export
vec_cast.character.ps_trim <- function(x, to, ...) as.character(vec_data(x))

# Propensity scores lie strictly between 0 and 1, so meeting an integer in the
# integers would round every one of them away. The common type is the one that
# holds both sets of values.
#' @export
vec_ptype2.ps_trim.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}
#' @export
vec_ptype2.integer.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}

#' @export
vec_ptype2.ps_trim.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "ps_trunc"))
  double()
}
#' @export
vec_ptype2.ps_trunc.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "ps_trim"))
  double()
}

#' @export
vec_cast.integer.ps_trim <- function(x, to, ...) {
  # A propensity score has no integer to be, so vctrs' own check reports the
  # loss rather than silently rounding it away.
  vec_cast(vec_data(x), integer(), x_arg = "ps_trim")
}

#' @export
vec_cast.ps_trim.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_trimmed_ps(xx, ps_trim_meta = trim_meta_for_cast(to, xx))
}

#' @export
vec_cast.ps_trim.ps_trunc <- function(x, to, ...) {
  # Convert ps_trunc to ps_trim (no trimming, just convert)
  new_trimmed_ps(
    vec_data(x),
    ps_trim_meta = list(
      method = "from_trunc",
      keep_idx = seq_along(x),
      trimmed_idx = integer(0),
      n_obs = length(x)
    )
  )
}

#' @export
vec_cast.ps_trunc.ps_trim <- function(x, to, ...) {
  # Convert ps_trim to ps_trunc (ignore NAs)
  cast_to_unbounded_ps_trunc(vec_data(x))
}

#' @export
vec_math.ps_trim <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_trim <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_trim <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}


#' @export
`[.ps_trim` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get original metadata
  meta <- ps_trim_meta(x)

  # Convert i to positive integer indices if needed
  # For logical indices, only allow length 1 or length n
  if (is.logical(i) && length(i) != 1 && length(i) != length(x)) {
    abort(
      "Logical subscript `i` must be size 1 or {length(x)}, not {length(i)}.",
      error_class = "propensity_length_error"
    )
  }
  i <- vec_as_location(i, n = length(x))

  # The subscript is in hand here, which is what re-indexing a record takes and
  # what the restore behind `NextMethod()` lacks: that restore would see a
  # record written for a different number of observations and drop it. Building
  # the result here keeps ordinary subsetting describing what it returns.
  new_trimmed_ps(
    vec_data(x)[i],
    ps_trim_meta = carry_trim_record(meta, length(x), i)
  )
}

# `i` holds the old positions the result is built from, in the order it holds
# them, among the `n_obs` observations it was taken from. A record that covers
# those observations re-indexes onto `i`; one that does not cannot be placed
# against them at all.
carry_trim_record <- function(meta, n_obs, i) {
  if (record_covers(meta, n_obs)) {
    reindex_trim_record(meta, i)
  } else {
    drop_trim_record(meta)
  }
}

#' @export
sort.ps_trim <- function(x, decreasing = FALSE, na.last = NA, ...) {
  meta <- ps_trim_meta(x)
  x_data <- vec_data(x)

  # `order()` returns the old positions in their new order, and drops the `NA`
  # positions entirely when `na.last = NA`, so it is the subscript the result is
  # built from and the record follows it.
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  new_trimmed_ps(
    x_data[ord],
    ps_trim_meta = carry_trim_record(meta, length(x), ord)
  )
}

#' @export
unique.ps_trim <- function(x, incomparables = FALSE, ...) {
  check_incomparables(incomparables, "ps_trim")

  # `vec_unique_loc()` names the position each retained value came from, which
  # is the subscript re-indexing the record takes. Without this the restore
  # behind vctrs' own method sees only a shorter vector and drops the record.
  loc <- vec_unique_loc(x)
  meta <- ps_trim_meta(x)

  # A retained value stands for every unit holding it, and the record can speak
  # for it only when all of those units share one standing. A trimmed score and
  # one that arrived missing are both `NA`, so they merge into one element that
  # neither status describes, and the record is dropped rather than left to
  # report that element as one of them.
  if (
    !is.matrix(x) &&
      record_covers(meta, length(x)) &&
      !merged_units_agree(vec_data(x), trim_unit_status(meta))
  ) {
    return(new_trimmed_ps(vec_data(x)[loc], drop_trim_record(meta)))
  }

  x[loc]
}

#' @export
rep.ps_trim <- function(x, ...) {
  # The positions `rep()` produces are the subscript the result is built from,
  # so the record follows them and a repeated unit is reported at each place it
  # now holds.
  x[rep(seq_along(x), ...)]
}

#' @export
summary.ps_trim <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @importFrom stats na.omit
#' @export
na.omit.ps_trim <- function(object, ...) {
  # Get metadata
  meta <- ps_trim_meta(object)

  # Get non-NA values and their positions
  not_na <- !is.na(object)
  kept_positions <- which(not_na)

  # Get the clean data
  clean_data <- vec_data(object)[not_na]

  # Create the result with na.action attribute
  result <- new_trimmed_ps(
    clean_data,
    ps_trim_meta = carry_trim_record(meta, length(object), kept_positions)
  )

  # Add na.action attribute as base na.omit does
  attr(result, "na.action") <- which(!not_na)
  class(attr(result, "na.action")) <- "omit"

  result
}

# Combining a single `ps_trim` hands back that vector, so its record still
# describes it. Every combine through vctrs restores against a zero-length
# prototype that cannot be told apart from one a caller supplied, and drops the
# positions. Anything other than one unnamed input with the default
# `recursive` and `use.names` goes to vctrs unchanged, and a matrix of scores,
# which is not a vctrs vector, keeps base `c()`'s flattening.
#' @export
c.ps_trim <- function(..., recursive = FALSE, use.names = TRUE) {
  if (
    ...length() == 1L &&
      !is.matrix(..1) &&
      is.null(...names()) &&
      isFALSE(recursive) &&
      isTRUE(use.names)
  ) {
    return(..1)
  }

  NextMethod()
}

#' @export
vec_restore.ps_trim <- function(x, to, ...) {
  # vec_data in case x is already a vctr
  data <- vec_data(x)
  meta <- ps_trim_meta(to)

  # Nothing rebuilding a `ps_trim` is handed the subscript behind a length
  # change, so a record written for a different number of observations cannot be
  # re-indexed onto the data arriving here. Zero-length data is exempt: a
  # prototype or an empty slice holds no observations, so no position in the
  # record contradicts it.
  #
  # Data restored against a zero-length `to` did not come from `to`: a
  # prototype, whether sliced from the one input of a combine or supplied by the
  # caller, is restored onto observations it was never written for, in an order
  # nothing here is told. Its positions are dropped without comment, as a
  # combine of several inputs drops them, and base `c()` of a single vector keeps
  # them by returning that vector.
  if (length(data) > 0 && length(to) == 0) {
    meta <- drop_trim_record(meta)
  } else if (length(data) > 0 && !record_covers(meta, length(data))) {
    meta <- drop_trim_record(meta)
  }

  new_trimmed_ps(data, ps_trim_meta = meta)
}

#' @export
quantile.ps_trim <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
anyDuplicated.ps_trim <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trim <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}

#' Refit a Propensity Score Model on Retained Observations
#'
#' @description
#' Re-estimates a propensity score model using only the observations retained
#' after trimming. This is the recommended intermediate step between
#' [ps_trim()] and weight calculation (e.g. [wt_ate()]):
#'
#' **`ps_trim()` -> `ps_refit()` -> `wt_*()`**
#'
#' Trimming changes the target population by removing observations with extreme
#' propensity scores. Refitting the model on the retained subset produces
#' propensity scores that better reflect this population, improving both model
#' fit and downstream weight estimation. Weight functions warn if a trimmed
#' propensity score has not been refit.
#'
#' @param trimmed_ps A `ps_trim` object returned by [ps_trim()]. Refitting reads
#'   the retained positions out of the trimming record, so an object whose
#'   record was dropped or no longer covers it raises an error of class
#'   `propensity_missing_meta_error`; see [ps_trim()].
#' @param model The original fitted model used to estimate the propensity
#'   scores (e.g. a [glm][stats::glm] or [multinom][nnet::multinom] object).
#'   The model is refit via [update()][stats::update] on the retained subset.
#'   Trimmed propensity scores are refit only by a model of the probability of
#'   the exposure; a model of a conditional mean, such as a [stats::lm()] fit
#'   or a gaussian [stats::glm()], never produced them and raises an error of
#'   class `propensity_model_family_error`. A dose model trimmed with
#'   `method = "density"` or `method = "resid"` is refit only by a model of
#'   the dose's conditional mean, such as a gaussian [stats::glm()],
#'   [stats::lm()], [mgcv::gam()], or [MASS::rlm()] fit, and any other model,
#'   such as one of a probability, raises the same error.
#' @param .data A data frame with one row per observation in `trimmed_ps`, in
#'   the same order. If `NULL` (the default), the data are recovered from
#'   `model`: its [model.frame()][stats::model.frame] when that already holds
#'   every variable the refit reads, and otherwise the data the model names,
#'   restricted by row name to the rows the model analyzed. A model fit without
#'   a data argument names none, and its variables are read out of the formula's
#'   environment instead. A formula that transforms a term, such as `z ~ log(x)`
#'   or a spline basis, stores that term already computed, so only the
#'   underlying variables let the transformation be recomputed from the retained
#'   rows. Pass `.data` when the data the model was fit on can no longer be
#'   reached.
#' @param ... Additional arguments passed to [update()][stats::update], such
#'   as `subset` or `weights`. Each is evaluated against the retained rows, the
#'   data the refit is handed: a name is read from their columns first and
#'   otherwise from the environment `ps_refit()` was called from, so
#'   `subset = x1 > 0` refits on the retained rows where `x1` is positive. A
#'   logical or index vector supplied instead indexes the retained rows, not
#'   the full data. The `.data` and `.env` pronouns are available to tell a
#'   column from a variable of the same name. Without `.data`, the retained
#'   rows hold only the variables `model` reads, so an argument naming any
#'   other column needs the data frame passed to `.data`.
#'
#' @details
#' ## Composing with a `subset`
#'
#' A `subset` in the original call has already chosen the sample the propensity
#' scores are about, and the trimming record indexes that sample rather than
#' every row the data carry. Refitting narrows that sample further, to the
#' retained rows, so the original `subset` is dropped from the call rather than
#' put to work a second time on rows it was never about. A `subset` passed
#' through `...` is an instruction of its own and is honored.
#'
#' ## Which level a refit predicts
#'
#' For a vector of scores, the refit predicts the probability of the level the
#' scores describe: when [ps_trim()] was given a model with its first level
#' named as focal, the refit reports one minus the model's prediction, and for
#' scores supplied directly it reports the level the model predicts by default.
#'
#' ## Arguments read from outside the formula
#'
#' `weights`, `offset`, and `na.action` in the original call are re-evaluated
#' against the retained rows. A `weights` or `offset` naming a column of the
#' data the model was fit on is read from that column and follows the retained
#' rows, whether the data are recovered from `model` or passed to `.data`. A
#' vector held outside the data cannot follow them: it keeps the length it had
#' and raises an error about differing variable lengths.
#'
#' Scores predicted from a fit with `na.action = na.exclude` are padded back to
#' the full length of the data, so they describe more observations than the fit
#' read and the trimming record indexes a sample the model never analyzed.
#' `ps_refit()` refuses such scores. Trim scores from a fit whose `na.action`
#' drops those rows instead.
#'
#' @return A `ps_trim` object with re-estimated propensity scores for retained
#'   observations and `NA` for trimmed observations. Use [is_refit()] to
#'   confirm refitting was applied. Refitting replaces the values a calibrated
#'   score held with predictions from the refit model, so a score trimmed after
#'   [ps_calibrate()] is no longer calibrated and [is_ps_calibrated()] answers
#'   `FALSE` for the result.
#'
#'   For a trimmed dose model, the retained values are the refit model's
#'   conditional means, and the spread in the record is re-estimated from the
#'   retained residuals under the recorded density family, unless the trim was
#'   made at a `.sigma` the caller supplied, which is kept. The recorded
#'   threshold describes the cut that was made and is left as it was, so the
#'   refit model can place a retained unit's conditional density below it;
#'   nothing is trimmed again.
#'
#' @seealso [ps_trim()] for the trimming step, [is_refit()] to check refit
#'   status, [wt_ate()] and other weight functions for the next step in the
#'   pipeline.
#'
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#'
#' # fit a propensity score model
#' ps_model <- glm(z ~ x, family = binomial)
#' ps <- predict(ps_model, type = "response")
#'
#' # trim -> refit -> weight pipeline
#' trimmed <- ps_trim(ps, lower = 0.1, upper = 0.9)
#' refit <- ps_refit(trimmed, ps_model)
#' wts <- wt_ate(refit, .exposure = z)
#'
#' is_refit(refit)
#'
#' @export
ps_refit <- function(trimmed_ps, model, .data = NULL, ...) {
  assert_class(trimmed_ps, "ps_trim")
  meta <- ps_trim_meta(trimmed_ps)

  # Get the number of observations
  n_obs <- if (is.matrix(trimmed_ps)) nrow(trimmed_ps) else length(trimmed_ps)

  # Refitting subsets the data by the retained positions, so a record that does
  # not cover these scores would name rows of a sample it was never about.
  check_trim_record(meta, n_obs, "ps_refit()")

  if (length(meta$keep_idx) == 0) {
    abort(
      "No retained rows to refit on (all were trimmed).",
      error_class = "propensity_no_data_error"
    )
  }

  # Checked before any data are recovered or any model is refit, so that a
  # model that could never have produced the trimmed values is refused for that
  # rather than for whatever refitting it would run into.
  check_refit_model(meta, model, is.matrix(trimmed_ps))
  density_record <- is_density_trim_record(meta)

  from_model <- is.null(.data)
  if (from_model) {
    .data <- refit_data(model)
  }

  if (nrow(.data) != n_obs) {
    # Data recovered from the model are the rows the model analyzed, so scores
    # that outnumber those rows were read over a longer sample than the fit did,
    # and an `na.action` that pads rather than drops is what does that. Scores
    # that fall short of them describe some narrower sample, which no padding
    # could account for.
    padding_hint <- if (from_model && n_obs > nrow(.data)) {
      c(
        i = "Scores predicted from a fit with {.code na.action = na.exclude} are padded back to the full length of the data, so they outnumber the rows the model analyzed.",
        i = "Trim scores from a fit whose {.arg na.action} drops those rows, such as {.fun stats::na.omit}."
      )
    }

    abort(
      c(
        "{.arg .data} must have the same number of rows as observations in {.arg trimmed_ps}.",
        x = "{.arg .data} has {nrow(.data)} row{?s}.",
        x = "{.arg trimmed_ps} has {n_obs} observation{?s}.",
        padding_hint
      ),
      error_class = "propensity_length_error"
    )
  }

  # refit on untrimmed rows. The retained rows are the sample to fit on, so a
  # `subset` the original call carried has already chosen its rows: putting it
  # to work again would choose among rows it was never about. It is dropped
  # unless the caller names one, which is an instruction of its own.
  data_sub <- .data[meta$keep_idx, , drop = FALSE]
  refit_args <- refit_extra_args(rlang::enquos(...), data_sub, from_model)
  if (!"subset" %in% names(refit_args)) {
    refit_args <- c(refit_args, list(subset = NULL))
  }
  refit_call <- rlang::inject(
    stats::update(model, data = data_sub, !!!refit_args, evaluate = FALSE)
  )
  refit_model <- eval(refit_call_function(refit_call, model))

  if (density_record) {
    return(refit_dose_trim(meta, refit_model, data_sub, n_obs))
  }

  # predict new PS for all rows
  if (is.matrix(trimmed_ps)) {
    # For matrix propensity scores (categorical exposures)
    new_ps <- matrix(NA_real_, nrow = n_obs, ncol = ncol(trimmed_ps))
    colnames(new_ps) <- colnames(trimmed_ps)

    # Predict probabilities for retained observations
    if (inherits(refit_model, "multinom")) {
      # For multinomial models from nnet
      pred_probs <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "probs"
      )
      # Ensure it's a matrix
      if (!is.matrix(pred_probs)) {
        pred_probs <- matrix(pred_probs, nrow = 1)
      }
      new_ps[meta$keep_idx, ] <- pred_probs
    } else {
      # Generic prediction
      new_ps[meta$keep_idx, ] <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "response"
      )
    }
  } else {
    # For vector propensity scores (binary exposures)
    new_ps <- rep(NA_real_, n_obs)
    refit_ps <- predict_binary_ps(refit_model, data_sub)
    if (isTRUE(meta$focal_inverted)) {
      refit_ps <- 1 - refit_ps
    }
    new_ps[meta$keep_idx] <- refit_ps
  }

  meta$refit <- TRUE

  # Every retained score is now a prediction from the model fit on the retained
  # rows, and nothing calibrates those predictions. A calibration the trimming
  # recorded describes scores that are no longer here, so it is dropped rather
  # than carried onto scores it was never about.
  meta$calibrated <- NULL

  new_trimmed_ps(
    x = new_ps,
    ps_trim_meta = meta
  )
}

# The arguments a caller passes to `ps_refit()` for `update()`, evaluated
# before the refit call is built. The call is evaluated in the frame of
# `ps_refit()`, where an expression such as `subset = x1 > 0` names nothing, so
# each argument is read here the way the fitting function would read it on the
# retained rows: a column of `data` first, and otherwise a name in the frame the
# argument was written in. `formula.` is a formula to update the model's own
# with rather than something to read from the rows, so it is read without them.
#
# Data recovered from the model hold only the variables the model itself reads,
# so an argument naming any other column fails there, and the error says how to
# make that column available. Only a failure to find a name is relabeled: any
# other error the argument raises is the caller's own and passes through.
refit_extra_args <- function(
  dots,
  data,
  from_model,
  call = rlang::caller_env()
) {
  arg_names <- names(dots)
  values <- lapply(seq_along(dots), function(i) {
    if (identical(arg_names[[i]], "formula.")) {
      return(rlang::eval_tidy(dots[[i]]))
    }
    if (!from_model) {
      return(rlang::eval_tidy(dots[[i]], data = data))
    }
    tryCatch(
      rlang::eval_tidy(dots[[i]], data = data),
      error = function(cnd) {
        if (!is_unfound_name_error(cnd, dots[[i]], data)) {
          rlang::cnd_signal(cnd)
        }
        abort(
          c(
            "Can't evaluate {.arg {arg_names[[i]]}} against the retained rows.",
            i = "Without {.arg .data}, only the variables {.arg model} reads
                 are available as columns. Pass the data frame to
                 {.arg .data} to use any other column."
          ),
          error_class = "propensity_no_data_error",
          call = call,
          parent = cnd
        )
      }
    )
  })
  names(values) <- arg_names
  values
}

# Whether evaluating `quo` against `data` failed because a name it reads is
# found neither among the columns nor in the environment it was written in,
# which is the failure passing `.data` can repair. Base R gives that failure no
# class of its own, so it is recognized by its message, compared in the
# session's language against each such name. A `.data$` lookup of a missing
# column raises rlang's own class. A `.env$` lookup is never a column, so the
# names it reads are left out and its failure passes through.
is_unfound_name_error <- function(cnd, quo, data) {
  if (inherits(cnd, "rlang_error_data_pronoun_not_found")) {
    return(TRUE)
  }

  env <- rlang::quo_get_env(quo)
  read <- unique(column_candidate_names(rlang::quo_get_expr(quo)))
  unfound <- read[
    !read %in% names(data) &
      !vapply(read, exists, logical(1), envir = env)
  ]
  if (length(unfound) == 0) {
    return(FALSE)
  }

  messages <- gettextf("object '%s' not found", unfound, domain = "R")
  conditionMessage(cnd) %in% messages
}

# The bare names an expression reads as variables, which are the names that
# could be columns. A function's name, the field after `$` or `@`, and a
# lookup through either pronoun are not.
column_candidate_names <- function(expr) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %in% c("", ".data", ".env")) {
      return(character())
    }
    return(name)
  }
  if (!is.call(expr)) {
    return(character())
  }

  head <- expr[[1]]
  args <- as.list(expr)[-1]
  if (identical(head, quote(`$`)) || identical(head, quote(`@`))) {
    args <- args[1]
  }
  if (
    (identical(head, quote(`$`)) || identical(head, quote(`[[`))) &&
      identical(args[[1]], quote(.env))
  ) {
    args <- args[-1]
    if (identical(head, quote(`$`))) {
      args <- list()
    }
  }
  if (!is.symbol(head)) {
    args <- c(list(head), args)
  }

  unlist(lapply(args, column_candidate_names), use.names = FALSE)
}

# Whether `model` can refit what the trimming record was made from. A record of
# trimmed propensity scores needs a model that fits a probability: one that fits
# a probability for every level, or one whose family fits the probability of a
# binary exposure, which a two-level `multinom` also answers. A model of a
# conditional mean never produced those scores, however close to the unit
# interval its fitted values fall.
#
# A record of a trimmed dose model needs the other kind. A `multinom` is refused
# by its class before the family is read, because the continuous family check
# reads a model that carries no family as a least squares fit, and a two-level
# `multinom` carries none and does not answer to `model_fits_levels()` either.
#
# The shape of the trimmed scores decides which kind of probability model is
# needed: a matrix has a column for every level and a vector holds one
# probability, so `matrix_record` is read from `trimmed_ps` itself.
check_refit_model <- function(
  meta,
  model,
  matrix_record,
  call = rlang::caller_env()
) {
  density_record <- is_density_trim_record(meta)

  if (!density_record) {
    if (matrix_record) {
      if (model_fits_levels(model)) {
        return(invisible(NULL))
      }

      abort(
        c(
          "A matrix of trimmed propensity scores can only be refit with a
           model of the probability of every exposure level.",
          x = "{.arg trimmed_ps} holds one column per level, and
               {.arg model} is {.cls {class(model)[[1]]}}, which fits a single
               probability.",
          i = "Refit with the model the scores were read from, such as a
               {.fun nnet::multinom} fit to all of the exposure's levels."
        ),
        error_class = "propensity_model_family_error",
        call = call
      )
    }

    # A vector of scores is the probability of one level, so a model that
    # reports a probability for every level of three or more has nothing to
    # refit it with, and the binary family check refuses such a model.
    check_binary_model_family(
      model,
      arg = "model",
      problem = "Trimmed propensity scores can only be refit with a model of
                 the probability of the exposure.",
      remedy = "{.arg trimmed_ps} holds propensity scores this model never
                produced. To set aside the units whose dose is implausible
                under a model of a continuous exposure, trim that model itself
                with {.code ps_trim(method = \"density\")}.",
      call = call
    )
    return(invisible(NULL))
  }

  dose_problem <- "A trimmed dose model can only be refit with a model of the
                   exposure's conditional mean."

  if (inherits(model, "multinom") || model_fits_levels(model)) {
    abort(
      c(
        dose_problem,
        x = "The trimming record holds a dose model's conditional means, and
             {.arg model} is {.cls {class(model)[[1]]}}, a model of the
             probabilities of a discrete exposure's levels.",
        i = "Refit with the model of the continuous exposure the trimming was
             made from, such as one fit with {.fun gaussian}, {.fun lm},
             {.fun mgcv::gam}, or {.fun MASS::rlm}."
      ),
      error_class = "propensity_model_family_error",
      call = call
    )
  }

  check_continuous_model_family(
    model,
    arg = "model",
    problem = dose_problem,
    remedy = "The trimming record holds a dose model's conditional means, so
              refit with the model of the continuous exposure the trimming was
              made from, such as one fit with {.fun gaussian}, {.fun lm},
              {.fun mgcv::gam}, or {.fun MASS::rlm}.",
    call = call
  )
}

# Whether a trimming record was made by trimming a dose model, whose values are
# conditional means read under a density family rather than propensity scores.
is_density_trim_record <- function(meta) {
  isTRUE(meta$method %in% ps_trim_dose_methods)
}

# A refit call names the fitting function the way the original fit recorded it.
# `MASS::rlm()` records itself as a bare `rlm`, which cannot be found unless its
# package is attached, so a bare name that does not resolve is qualified with
# the namespace that defines the methods for the model's class.
refit_call_function <- function(refit_call, model) {
  # `predict()` is the generic looked up, because every model class a refit
  # accepts has to answer it.
  fn <- refit_call[[1]]
  if (!is.symbol(fn) || exists(as.character(fn), mode = "function")) {
    return(refit_call)
  }

  fn_name <- as.character(fn)
  for (cls in class(model)) {
    method <- utils::getS3method("predict", cls, optional = TRUE)
    home <- if (is.function(method)) environment(method)
    if (
      isNamespace(home) &&
        exists(fn_name, envir = home, mode = "function", inherits = FALSE)
    ) {
      refit_call[[1]] <- call("::", as.symbol(getNamespaceName(home)), fn)
      return(refit_call)
    }
  }

  refit_call
}

# The refit of a trimmed dose model. The retained values are the refit model's
# conditional means, and the spread the density is read at is re-estimated from
# the retained residuals under the recorded family, because the density ratio
# takes the spread as an argument of its own. A spread the caller supplied is
# theirs whatever the residuals would estimate, and is kept. The threshold
# describes the cut that was made, so it is left as recorded even where the
# refit model puts a retained unit's density below it.
refit_dose_trim <- function(meta, refit_model, data_sub, n_obs) {
  keep_idx <- meta$keep_idx
  new_mu <- rep(NA_real_, n_obs)
  new_mu[keep_idx] <- stats::predict(
    refit_model,
    newdata = data_sub,
    type = "response"
  )

  # The spread is read from the refit model's own residuals, which cover
  # exactly the rows it analyzed. A `subset` or a covariate missing on a kept
  # row leaves fewer of those than there are kept rows.
  if (!identical(meta$sigma_kind, "supplied")) {
    residuals <- as.numeric(stats::residuals(refit_model, type = "response"))
    meta$sigma <- density_scale_estimate(
      residuals[!is.na(residuals)],
      meta$density
    )
  }

  meta$refit <- TRUE

  new_trimmed_ps(x = new_mu, ps_trim_meta = meta)
}

# The variables the refit call reads: the ones the formula names, plus any named
# by `weights` or `offset`, which are evaluated against whatever data the refit
# is handed. The original `subset` is dropped from that call, so the variables it
# names are not among them.
refit_call_vars <- function(model, model_formula) {
  model_call <- stats::getCall(model)
  extras <- lapply(
    list(model_call$weights, model_call$offset),
    function(arg) if (is.null(arg)) character() else all.vars(arg)
  )

  unique(c(all.vars(model_formula), unlist(extras, use.names = FALSE)))
}

# The data `ps_refit()` refits on when the caller names none. The model frame
# holds each term as the formula computed it, so `z ~ log(x)` stores a column of
# logged values and nothing to recompute them from; refitting on the retained
# rows has to recompute the transformation from those rows alone, which only the
# raw variables allow. It also stores `weights` and `offset` under fixed names
# rather than the ones the call reads. The model frame is enough whenever it
# already carries every variable the refit reads, and it costs nothing to read.
#
# Otherwise the data the model names are read back and cut down to the rows the
# model analyzed. The whole frame is used, not just the variables the formula
# names: it is the frame the original call read, so `weights`, `offset`, and a
# `.` in the formula mean against it what they meant when the model was fit. A
# fit that named no data has no frame to read back, and its variables come from
# the formula's environment instead.
refit_data <- function(model, call = rlang::caller_env()) {
  model_formula <- stats::formula(model)
  model_frame <- stats::model.frame(model)

  if (all(refit_call_vars(model, model_formula) %in% names(model_frame))) {
    return(model_frame)
  }

  raw <- tryCatch(
    {
      model_data <- eval(stats::getCall(model)$data, environment(model_formula))
      if (is.data.frame(model_data)) {
        model_data
      } else {
        stats::get_all_vars(model_formula, data = model_data)
      }
    },
    error = function(cnd) cnd
  )

  if (inherits(raw, "condition")) {
    abort(
      c(
        "Can't recover the data {.arg model} was fit on.",
        x = conditionMessage(raw),
        i = "Pass the data frame to {.arg .data}."
      ),
      error_class = "propensity_no_data_error",
      call = call
    )
  }

  # The recovered data carry every row they were read from, while the propensity
  # scores, and with them the trimming record, are about the rows the model
  # analyzed. Row names line the two up across the rows a fit dropped as missing
  # or a `subset` excluded.
  rows <- match(rownames(model_frame), rownames(raw))

  if (anyNA(rows)) {
    abort(
      c(
        "Can't match the rows {.arg model} was fit on to the data it read.",
        i = "Pass the data frame to {.arg .data}."
      ),
      error_class = "propensity_no_data_error",
      call = call
    )
  }

  raw[rows, , drop = FALSE]
}

#' Check if propensity scores have been refit
#'
#' `is_refit()` tests whether `x` is a `ps_trim` object whose propensity
#' model has been refit on the retained (non-trimmed) observations via
#' [ps_refit()].
#'
#' On a [psw] vector built from trimmed propensity scores, the answer comes from
#' the trimming record carried on the weights. `is_refit()` reads a single flag
#' out of that record rather than a position, so unlike [is_unit_trimmed()] it
#' answers from any record present, whatever length the weights have since taken
#' on. It raises an error of class `propensity_missing_meta_error` only when
#' weights marked as trimmed carry no record at all, which is what a subset of
#' such weights leaves behind. Weights that were never trimmed have nothing
#' missing and return `FALSE`.
#'
#' @param x An object to test (typically a [ps_trim] vector).
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_refit()] to refit a propensity model after trimming,
#'   [ps_trim()] to trim propensity scores.
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' trimmed <- ps_trim(ps, lower = 0.2, upper = 0.8)
#' is_refit(trimmed)
#'
#' refit <- ps_refit(trimmed, fit)
#' is_refit(refit)
#'
#' @export
is_refit <- function(x) {
  UseMethod("is_refit")
}

#' @export
is_refit.default <- function(x) {
  FALSE
}

#' @export
is_refit.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  isTRUE(meta$refit)
}

#' Extract trimming metadata from a `ps_trim` object
#'
#' @description `ps_trim_meta()` returns the metadata list attached to a
#'   `ps_trim` object by [ps_trim()].
#'
#' @param x A `ps_trim` object.
#' @return A named list with elements:
#' \describe{
#'   \item{`method`}{Character string indicating the trimming method used.}
#'   \item{`keep_idx`}{Integer vector of retained observation indices.}
#'   \item{`trimmed_idx`}{Integer vector of trimmed observation indices.}
#'   \item{`n_obs`}{The number of observations those indices describe.}
#'   \item{`lower`, `upper`}{Numeric cutoffs, when applicable.}
#'   \item{`refit`}{Logical, `TRUE` if the model was refit via [ps_refit()].}
#'   \item{`focal_inverted`}{For a vector of scores, logical, `TRUE` if the
#'     scores are one minus the probability the fitted model reports, which
#'     [ps_refit()] reads to predict the same level.}
#' }
#' Additional method-specific elements (e.g. `cutoff`, `delta`, `lambda`) may
#' also be present.
#'
#' `keep_idx`, `trimmed_idx`, and `n_obs` are absent when an operation that
#' changed the number of observations dropped the record; see [ps_trim()].
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] and
#'   [is_unit_trimmed()] for predicate queries.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' ps_trim_meta(trimmed)
#'
#' @export
ps_trim_meta <- function(x) {
  attr(x, "ps_trim_meta")
}
