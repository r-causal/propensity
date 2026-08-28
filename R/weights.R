#' Calculate propensity score weights
#'
#' @description
#' Compute inverse probability weights for causal inference under different
#' estimands. Each function targets a different population:
#'
#' - `wt_ate()`: **Average Treatment Effect** -- the full population.
#' - `wt_att()`: **Average Treatment Effect on the Treated** -- the treated
#'   (focal) group.
#' - `wt_atu()`: **Average Treatment Effect on the Untreated** -- the
#'   untreated (reference) group. `wt_atc()` is an alias.
#' - `wt_atm()`: **Average Treatment Effect for the Evenly Matchable** --
#'   units with the most overlap.
#' - `wt_ato()`: **Average Treatment Effect for the Overlap Population** --
#'   weights proportional to overlap.
#' - `wt_entropy()`: **Entropy-weighted Average Treatment Effect** --
#'   the entropy-tilted population.
#' - `wt_cens()`: **Inverse probability of censoring weights** -- uses the
#'   same formula as `wt_ate()` but labels the estimand `"uncensored"`. Use
#'   these to adjust for censoring in survival analysis, not for treatment
#'   weighting. Remaining uncensored is a two-level event, so `wt_cens()`
#'   supports binary and continuous exposures only.
#'
#' `.propensity` accepts a numeric vector of predicted probabilities, a
#' `data.frame` or matrix of per-level probabilities, a fitted `glm` object, or
#' a modified propensity score created by [ps_trim()], [ps_trunc()],
#' [ps_refit()], or [ps_calibrate()].
#'
#' All functions return a [`psw`] object -- a numeric vector that tracks the
#' estimand, stabilization status, and any trimming or truncation applied.
#'
#' @details
#' ## Exposure types
#'
#' All weight functions support binary exposures. `wt_ate()` and `wt_cens()`
#' also support continuous exposures. All except `wt_cens()` support
#' categorical exposures; naming a categorical exposure to `wt_cens()`, or
#' handing it one detection reads as categorical, is an error of class
#' `causalgenerics_unsupported_exposure_type`.
#'
#' - **Binary**: `.exposure` is a two-level vector (e.g., 0/1, logical, or a
#'   two-level factor). `.propensity` is a numeric vector of P(treatment | X),
#'   or a matrix or data frame holding it in one of its columns.
#' - **Categorical**: `.exposure` is a factor or character vector with 3+
#'   levels. `.propensity` must be a matrix or data frame with one column per
#'   level, where rows sum to 1.
#' - **Continuous**: `.exposure` is a numeric vector. `.propensity` is a
#'   vector of conditional means (fitted values). Weights are a ratio of
#'   densities, whose family is chosen by `.density`; stabilization is strongly
#'   recommended.
#' - **Auto** (default): Detects the exposure type from `.exposure`.
#'
#' ## Stabilization
#'
#' Setting `stabilize = TRUE` multiplies the base weight by an estimate of
#' P(A) (binary) or f_A(A) (continuous), reducing variance. When no
#' `stabilization_score` is supplied, that estimate is the marginal mean of
#' `.exposure` for a binary or categorical exposure, and for a continuous
#' exposure the marginal density of the family `.density` names, evaluated at
#' the population mean and standard deviation of `.exposure`. Stabilization is
#' supported for ATE and censoring weights (`wt_ate()` and `wt_cens()`) and is
#' strongly recommended for continuous exposures.
#'
#' For a continuous exposure, `numerator` chooses how that marginal density is
#' arrived at. `"marginal"`, the default, reads the family `.density` names at
#' the population mean and standard deviation of `.exposure`. Those two moments
#' are parameters of the weights, and [ipw()] estimates them alongside the rest
#' of its parameter vector, so the standard errors account for the numerator
#' having been estimated. `"integrated"` marginalizes the conditional density
#' numerically instead: it averages \eqn{f_{A|X}(t \mid X_i)} over the units at
#' each of 50 points spanning `.exposure`, then interpolates that average back
#' to each observed exposure with a cubic spline. It estimates no parameters of
#' its own, and it is what WeightIt has done since version 2.0.0.
#'
#' The two agree whenever the family is normal and the fitted conditional means
#' are themselves normal, since an average of normal densities centered on
#' normally distributed means is again a normal density. They separate as
#' those means grow skewed or bimodal, or as the family's tails grow heavier
#' than the normal's.
#'
#' The marginal numerator is the default because it is the published convention
#' (Robins et al., 2000; Naimi et al., 2014) and because it held up: in a
#' simulation study run for this package, over linear, skewed, heteroskedastic,
#' heavy-tailed, and bimodal designs, marginalizing never improved bias, root
#' mean squared error, or weighted covariate balance by more than Monte Carlo
#' noise, including in the scenarios drawn to favor it, and in the heavy-tailed
#' cells it was worse: it inflated the root mean squared error, and in some
#' replicates the interpolated density came back negative. The integrated
#' numerator is here for
#' agreement with WeightIt, and for a conditional density whose marginal has no
#' closed form in the same family. Being read off an interpolation rather than
#' a formula, it can dip below zero where the density on the grid comes close
#' to it, which is refused with an error of class `propensity_density_error`
#' rather than turned into a negative weight.
#'
#' The default numerator conditions on nothing, and `stabilization_score`
#' replaces it with one of your own, which leaves `numerator` nothing to
#' choose: the two cannot be supplied together. The case that pays is a
#' numerator conditioning on a variable the model the estimates are read from
#' also reads, such as an effect modifier. Fit that numerator and evaluate it at
#' the exposure each unit actually took, which is the shape the default already
#' has:
#'
#' ```
#' num <- glm(A ~ V, data = dat, family = binomial())
#' p <- fitted(num)
#' wt_ate(ps_mod, stabilize = TRUE,
#'        stabilization_score = ifelse(dat$A == 1, p, 1 - p))
#' ```
#'
#' Passing `fitted(num)` on its own is P(A = 1 | V) for every unit, which is not
#' the probability of the exposure an untreated unit took, so the untreated
#' weights would carry the wrong numerator. Whether a numerator may condition on
#' `V` at all is a question about the reported model rather than about the
#' weights; see **Effect modification** in [ipw()].
#'
#' ## Handling extreme weights
#'
#' Extreme weights signal positivity violations, poor model fit, or limited
#' overlap. You can address them by:
#'
#' - Choosing an overlap-focused estimand (`wt_ato()`, `wt_atm()`,
#'   `wt_entropy()`), which down-weight units in regions of poor overlap.
#' - Trimming ([ps_trim()]) or truncating ([ps_trunc()]) propensity scores
#'   before computing weights.
#' - Calibrating weights with [ps_calibrate()].
#' - Stabilizing ATE weights (`stabilize = TRUE`).
#'
#' See the [halfmoon](https://CRAN.R-project.org/package=halfmoon) package
#' for weight diagnostics and visualization.
#'
#' ## Propensity scores at 0 and 1
#'
#' Propensity scores must lie strictly inside (0, 1), since a score at either
#' endpoint leaves the weight undefined. A score of exactly 0 or 1 is refused
#' with an error of class `propensity_range_error`. For a categorical exposure
#' that rule reads every cell of the matrix rather than only the column holding
#' each unit's observed level, which is the same open interval the binary path
#' enforces on the single score and its complement.
#'
#' [ipw()] applies a narrower rule to the propensity scores it rebuilds from its
#' own propensity score model. For a categorical exposure it refuses only a
#' probability of exactly zero at the level a unit was actually assigned, since
#' that alone leaves the unit's own weight undefined, and a column for a level
#' the unit was not assigned may reach zero without the fit being refused. The
#' difference is deliberate: the weight functions can refuse a matrix that holds
#' nothing [ipw()] would have objected to. The weights still have to be built
#' before [ipw()] can be reached, so the refusal here is the one to resolve.
#'
#' [nnet::multinom()] reaches the endpoints readily. Under separation the
#' softmax puts the probability at a unit's assigned level at exactly 1 in
#' double precision, and the columns for the other levels can underflow to
#' exactly 0. Trimming is no way around the refusal: [ps_trim()] and
#' [ps_trunc()] validate a categorical matrix under the same open interval and
#' refuse it before either reaches a threshold. Bound the fitted probabilities
#' away from 0 and 1 yourself, renormalizing each row to sum to 1, or refit the
#' propensity score model so that it does not separate.
#'
#' ## Weight formulas
#'
#' ### Binary exposures
#'
#' For binary treatments (\eqn{A \in \{0, 1\}}), with propensity score
#' \eqn{e(X) = P(A=1 \mid X)}:
#'
#' - **ATE**: \eqn{w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}}
#' - **ATT**: \eqn{w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}}
#' - **ATU**: \eqn{w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)}
#' - **ATM**: \eqn{w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}
#' - **ATO**: \eqn{w = A \cdot (1-e(X)) + (1-A) \cdot e(X)}
#' - **Entropy**: \eqn{w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot (1-e(X))}}, where \eqn{h(e) = -[e \log(e) + (1-e) \log(1-e)]}
#'
#' The entropy weight tilts the propensity score by \eqn{h(e)}, the entropy of
#' the score itself, in the sense of Zhou, Matsouaka, and Thomas (2020). It is
#' not entropy balancing, which solves for weights that satisfy exact covariate
#' moment constraints rather than tilting a fitted propensity score.
#'
#' ### Continuous exposures
#'
#' Weights use the density ratio
#' \eqn{w = f_A(A) / f_{A|X}(A \mid X)}, where \eqn{f_A} is the marginal
#' density of the exposure and \eqn{f_{A|X}} is the conditional density the
#' propensity model estimates. Only `wt_ate()` and `wt_cens()` support
#' continuous exposures.
#'
#' Both densities are read on a standardized residual. The conditional density
#' is evaluated at \eqn{z_i = (A_i - \hat{A}_i) / \sigma}, where
#' \eqn{\hat{A}_i} is the fitted conditional mean in `.propensity` and
#' \eqn{\sigma} is the residual spread; the marginal density is evaluated at
#' \eqn{z^A_i = (A_i - \bar{A}) / s_A}, where \eqn{\bar{A}} and \eqn{s_A} are
#' the population mean and standard deviation of `.exposure`. Each density is
#' then divided by the spread that standardized it, the Jacobian of that change
#' of variable, which returns both to the exposure's own units so that each
#' integrates to one:
#'
#' \deqn{w_i = \frac{g(z^A_i) / s_A}{g(z_i) / \sigma}}
#'
#' `.density` chooses the family \eqn{g}:
#'
#' | `.density` | \eqn{g(z)} |
#' | --- | --- |
#' | `"normal"`, [dens_normal()] | the standard normal density, the default |
#' | `"laplace"`, [dens_laplace()] | \eqn{\exp(-\lvert z \rvert) / 2} |
#' | [dens_t()] | Student's t with `df` degrees of freedom |
#' | `"kernel"`, [dens_kernel()] | a kernel density estimate of the standardized residuals |
#' | [dens_fn()], or a bare function | a density you write yourself |
#'
#' The numerator and the denominator take the same family, so the choice
#' describes the shape of the residuals rather than the estimand. A heavier tail
#' than the normal's, such as [dens_t()] or `"laplace"`, holds down the weight
#' of a unit whose exposure the model fits poorly, and `"kernel"` assumes no
#' shape at all, at the cost of a density that is not a smooth function of the
#' model's parameters. A `stabilization_score` replaces the numerator outright
#' and `stabilize = FALSE` leaves it out, so under either the family describes
#' the denominator alone. Under `numerator = "integrated"` it describes both
#' again, the numerator there being the same conditional density averaged over
#' the units rather than the same family read at the exposure's own moments;
#' see **Stabilization**.
#'
#' A unit whose exposure or fitted conditional mean is missing has no
#' standardized residual, and so no weight: the density is asked only about the
#' units that have one, and a kernel is fit on those. A propensity score
#' [ps_trim()] set aside is the ordinary way this arises.
#'
#' \eqn{\sigma} is the pooled residual standard deviation
#' \eqn{\sqrt{\mathrm{mean}((A - \hat{A})^2)}} unless `.sigma` supplies a single
#' standard deviation, or one for each unit, which the `glm` methods do not do
#' on their own. The marginal density is spread by \eqn{s_A} whatever `.sigma`
#' says, the spread of the exposure itself being no business of the conditional
#' model's.
#'
#' [ipw()] models the conditional density with a single pooled residual
#' variance, estimated jointly with the rest of the parameter vector. Weights
#' built with an observation-level `.sigma` therefore have no counterpart in
#' that estimating equation, and `ipw()` rejects them in its weight consistency
#' check. Build weights with the pooled default when the outcome model is
#' headed for `ipw()`.
#'
#' ### Categorical exposures
#'
#' For \eqn{K}-level treatments, weights take the tilting-function form
#' \eqn{w_i = h(\mathbf{e}_i) / e_{i,Z_i}}, where \eqn{e_{i,Z_i}} is the
#' propensity for unit \eqn{i}'s observed level and \eqn{h(\cdot)} depends
#' on the estimand:
#'
#' - **ATE**: \eqn{h(\mathbf{e}) = 1}
#' - **ATT**: \eqn{h(\mathbf{e}) = e_{\text{focal}}}
#' - **ATU**: \eqn{h(\mathbf{e}) = 1 - e_{\text{focal}}}
#' - **ATM**: \eqn{h(\mathbf{e}) = \min(e_1, \ldots, e_K)}
#' - **ATO**: \eqn{h(\mathbf{e}) = \bigl(\sum_k 1/e_k\bigr)^{-1}}
#' - **Entropy**: \eqn{h(\mathbf{e}) = -\sum_k e_k \log(e_k)}
#'
#' @param .propensity Propensity scores in one of several forms:
#'   * A **numeric vector** of predicted probabilities (binary/continuous).
#'   * A **data frame** or **matrix** with one column per exposure level. Both
#'     shapes are read the same way, for categorical exposures and for binary
#'     ones alike, including the default choice of column described under
#'     `.propensity_col`. That argument is itself a formal of the data frame
#'     methods only, so selecting a column by name means passing a data frame.
#'   * A fitted **`glm`** object -- fitted values are extracted automatically.
#'   * A modified propensity score created by [ps_trim()], [ps_trunc()],
#'     [ps_refit()], or [ps_calibrate()].
#'
#'   For a binary exposure, `.propensity` is the probability of the focal
#'   level, so a numeric vector must be supplied on the scale of whichever
#'   level `.focal_level` or `.reference_level` resolves to. The `glm` methods
#'   derive it from the fitted values, which give the probability of the
#'   response's second level, and subtract them from one when the resolved
#'   focal level is the first level instead. A data frame or matrix reduces to
#'   a single column, read on the same scale; see `.propensity_col` for how
#'   that column is chosen. A matrix reduces by those same defaults, but
#'   `.propensity_col` belongs to the data frame methods alone, so a matrix
#'   whose column you want to name has to be converted with `as.data.frame()`
#'   first.
#' @param .exposure The exposure (treatment) variable. For binary exposures, a
#'   numeric 0/1 vector, logical, or two-level factor. For categorical
#'   exposures, a factor or character vector. For continuous exposures, a
#'   numeric vector. Optional when `.propensity` is a `glm` object (extracted
#'   from the model). Missing values are not counted as a level of their own,
#'   and are carried through to the weights as missing.
#' @param exposure_type Type of exposure: `"auto"` (default), `"binary"`,
#'   `"categorical"`, or `"continuous"`. `"auto"` detects the type from
#'   `.exposure`. Not every weight function answers every type; see **Exposure
#'   types** in Details for which does which. Naming a type a function does not
#'   support, or supplying an exposure detection reads as one, is an error of
#'   class `causalgenerics_unsupported_exposure_type`.
#' @param .sigma Observation-level standard deviations of the conditional
#'   density for continuous exposures, one per observation (e.g.,
#'   `influence(model)$sigma`). Optional: with none supplied, including when
#'   `.propensity` is a `glm` object, the conditional density uses the pooled
#'   residual standard deviation of `.exposure` around `.propensity`. Weights
#'   built with an observation-level `.sigma` cannot be used with [ipw()]; see
#'   **Continuous exposures** in Details.
#'
#'   Must be numeric, holding either a single standard deviation or one per
#'   observation, and applies only to continuous exposures. `.sigma` sits in
#'   the third position, which is where a value meant for `exposure_type`
#'   arrives when it is supplied without a name, so anything else is refused
#'   with an error of class `propensity_sigma_error`.
#' @param .density The family of the conditional density for continuous
#'   exposures, described under **Continuous exposures** in Details. One of the
#'   strings `"normal"` (the default), `"laplace"`, or `"kernel"`; a
#'   specification built by [dens_normal()], [dens_laplace()], [dens_t()],
#'   [dens_kernel()], or [dens_fn()]; or a bare function of the standardized
#'   residual, which is wrapped with [dens_fn()].
#'
#'   A function is called on the whole vector of standardized residuals at once
#'   rather than on one residual at a time, and is called a second time on the
#'   standardized exposure when the marginal numerator stabilizes the weights.
#'   Each call must return one finite, non-negative value for each value it was
#'   given, not every one of them zero. Anything else is refused with an error of
#'   class `propensity_density_error`, which is also what a family other than the
#'   default gets for a binary or categorical exposure, whose weights are not a
#'   ratio of densities. The default is accepted for every exposure type and
#'   ignored outside a continuous one.
#'
#'   `.density` sits after `...` and so can only be supplied by name. The family
#'   the weights were built from is recorded on the result and read back with
#'   [density_meta()].
#' @param numerator How the marginal density that stabilizes a continuous
#'   exposure's weights is obtained, described under **Stabilization** in
#'   Details. Either `"marginal"`, the default, which reads the family
#'   `.density` names at the population mean and standard deviation of
#'   `.exposure`, or `"integrated"`, which averages the conditional density
#'   over the units on a grid spanning `.exposure` and interpolates the result
#'   back to each observed exposure.
#'
#'   `"integrated"` describes a continuous exposure alone, and needs a
#'   conditional density to marginalize: it is refused with an error of class
#'   `propensity_numerator_error` for a binary or categorical exposure, with
#'   `stabilize = FALSE`, with a `stabilization_score`, and with any `.sigma`.
#'   The default is accepted for every exposure type and ignored outside a
#'   continuous one.
#'
#'   `numerator` sits after `...` and so can only be supplied by name. The
#'   numerator the weights were built from is recorded on the result and read
#'   back with [density_meta()].
#' @param .treated `r lifecycle::badge("deprecated")` Use `.focal_level` instead.
#' @param .untreated `r lifecycle::badge("deprecated")` Use `.reference_level` instead.
#' @param .focal_level The value of `.exposure` representing the focal
#'   (treated) group. Every binary coding honors it: 0/1 numeric, logical,
#'   two-level factor, and two-level character exposures are all coded with the
#'   named level as focal, and a level the exposure never takes is an error.
#'   With no level named, a binary exposure defaults to its higher level, which
#'   is `1` for a 0/1 exposure, `TRUE` for a logical one, and the second of the
#'   two levels a factor or character exposure takes. Levels a factor declares
#'   but never takes are not candidates. Naming any other level reverses the
#'   coding, so `.propensity` must then hold the probability of the named
#'   level. Required for `wt_att()` and `wt_atu()` with categorical exposures.
#' @param .reference_level The value of `.exposure` representing the reference
#'   (control) group. For a binary exposure, naming it makes the exposure's
#'   other level focal, with the same consequence for `.propensity`, and a
#'   level the exposure never takes is an error. Automatically detected if not
#'   supplied.
#' @param ... These dots are for future extensions and must be empty.
#' @param stabilize If `TRUE`, multiply weights by an estimate of the marginal
#'   treatment probability (binary) or density (continuous). Only supported by
#'   `wt_ate()` and `wt_cens()`. See **Stabilization** in Details.
#' @param stabilization_score Optional stabilization multiplier to use instead
#'   of the default described under **Stabilization**: the marginal mean of
#'   `.exposure`, or its marginal normal density for a continuous exposure.
#'   Either a single value
#'   applied to every weight or a numeric vector holding one value per
#'   observation, which is multiplied into the weights observation by
#'   observation. Every value must be positive and finite, and any other length
#'   or value is refused with an error of class
#'   `propensity_stabilization_score_error`. A per-observation score is recorded
#'   on the result and is dropped, with a warning, by any operation that changes
#'   the length of the weight vector, since it cannot be re-indexed for the new
#'   length. Ignored when `stabilize = FALSE`.
#' @param .propensity_col Column to use when `.propensity` is a data frame
#'   with a binary exposure. Accepts a column name (quoted or unquoted) or
#'   numeric index. Whichever column is selected is read as the probability of
#'   the resolved focal level.
#'
#'   With no column named, the exposure's levels drive the choice. When the
#'   data frame has a column named for every level of `.exposure`, the column
#'   named for the resolved focal level is used, wherever it sits in the frame.
#'   Otherwise the second column is used, or the only column when the frame has
#'   just one. Falling back to position after `.focal_level` or
#'   `.reference_level` was supplied warns and reports the column used, since
#'   the named level could not be matched to a column; the fallback is silent
#'   when no level was named and when the frame has a single column.
#'
#'   Ignored for categorical exposures, where all columns are used.
#'
#' @return A [`psw`] vector (a double vector with class `psw`) carrying
#'   these attributes:
#'   - `estimand`: character, e.g. `"ate"`, `"att"`, `"uncensored"`.
#'   - `stabilized`: logical, whether stabilization was applied.
#'   - `stabilization_score`: double, the user-supplied stabilization score
#'     (one value, or one per observation), or `NULL` if none was supplied.
#'   - `trimmed`: logical, whether the propensity scores were trimmed.
#'   - `truncated`: logical, whether the propensity scores were truncated.
#'   - `calibrated`: logical, whether the propensity scores were calibrated.
#'   - `exposure_type`: character, the type of exposure the weights were built
#'     for, read back with [exposure_type()].
#'   - `density_meta`: for a continuous exposure, whose weights are a ratio of
#'     densities, the record of that ratio: the density family, the numerator
#'     that stabilized the weights, and where the residual spread came from.
#'     Read back with [density_meta()], and `NULL` for a binary or categorical
#'     exposure.
#'
#' @examples
#' # -- Binary exposure, numeric propensity scores ----------------------
#' set.seed(123)
#' ps <- runif(100, 0.1, 0.9)
#' trt <- rbinom(100, 1, ps)
#'
#' wt_ate(ps, trt)
#' wt_att(ps, trt)
#' wt_atu(ps, trt)
#' wt_atm(ps, trt)
#' wt_ato(ps, trt)
#' wt_entropy(ps, trt)
#'
#' # Stabilized ATE weights (reduces variance)
#' wt_ate(ps, trt, stabilize = TRUE)
#'
#' # Inspect the result
#' w <- wt_ate(ps, trt)
#' estimand(w)
#' summary(w)
#'
#' # -- Overlap-focused estimands handle extreme PS better --------------
#' ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
#' trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)
#'
#' max(wt_ate(ps_extreme, trt_extreme))
#' max(wt_ato(ps_extreme, trt_extreme))
#'
#' # -- From a fitted GLM -----------------------------------------------
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' trt2 <- rbinom(100, 1, plogis(0.5 * x1 + 0.3 * x2))
#' ps_model <- glm(trt2 ~ x1 + x2, family = binomial)
#'
#' # Exposure is extracted from the model automatically
#' wt_ate(ps_model)
#'
#' # -- Data frame input ------------------------------------------------
#' ps_df <- data.frame(
#'   control = c(0.9, 0.7, 0.3, 0.1),
#'   treated = c(0.1, 0.3, 0.7, 0.9)
#' )
#' exposure <- c(0, 0, 1, 1)
#' wt_ate(ps_df, exposure)
#' wt_ate(ps_df, exposure, .propensity_col = "treated")
#'
#' # -- Censoring weights -----------------------------------------------
#' cens_ps <- runif(50, 0.6, 0.95)
#' cens_ind <- rbinom(50, 1, cens_ps)
#' wt_cens(cens_ps, cens_ind)
#' estimand(wt_cens(cens_ps, cens_ind))  # "uncensored"
#'
#' @references
#' Barrett, M., D'Agostino McGowan, L., & Gerke, T. *Causal Inference in R*.
#' \url{https://www.r-causal.org/}
#'
#' Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity
#' score in observational studies for causal effects. *Biometrika*, 70(1),
#' 41--55.
#'
#' Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
#' propensity score analysis. *The International Journal of Biostatistics*,
#' 9(2), 215--234. (ATM weights)
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. *Journal of the American Statistical
#' Association*, 113(521), 390--400. (ATO weights)
#'
#' Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score weighting
#' under limited overlap and model misspecification. *Statistical Methods in
#' Medical Research*, 29(12), 3721--3756. (Entropy weights)
#'
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. In *Applied Bayesian Modeling and Causal Inference from
#' Incomplete-Data Perspectives* (pp. 73--84).
#'
#' Robins, J. M., Hernán, M. A., & Brumback, B. (2000). Marginal structural
#' models and causal inference in epidemiology. *Epidemiology*, 11(5),
#' 550--560.
#'
#' Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
#' Constructing inverse probability weights for continuous exposures: a
#' comparison of methods. *Epidemiology*, 25(2), 292--299.
#'
#' Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
#' using inverse probability of treatment weighting (IPTW). *Statistics in
#' Medicine*, 34(28), 3661--3679.
#'
#' @seealso
#' * [psw()] for the returned weight vector class.
#' * [ps_trim()], [ps_trunc()], [ps_refit()], and [ps_calibrate()] for
#'   modifying propensity scores before weighting.
#' * [ipw()] for inverse-probability-weighted estimation of causal effects.
#' * [ps_tilt()] for the tilting function each weight divides by the propensity
#'   score of the received exposure level, which also standardizes a
#'   g-computation estimate to the same target population.
#'
#' @export
wt_ate <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_ate")
}

#' @export
wt_ate.default <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_ate.numeric <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_ate",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    announce = !be_quiet(),
    call = call
  )

  check_sigma(.sigma, exposure_type, length(.exposure), call = call)

  # The density is resolved once, here, for every route: the model, data frame,
  # and modified-score methods all funnel through this one, and `wt_cens()`
  # reaches the continuous formula through it as well. The numerator is resolved
  # in the same place and for the same reason.
  .density <- as_density_spec(.density, arg = ".density", call = call)
  check_density_arg(.density, exposure_type, call = call)

  numerator <- rlang::arg_match(numerator, error_call = call)
  check_numerator(
    numerator,
    exposure_type,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .sigma = .sigma,
    call = call
  )

  # The exposure supplies the number of observations the score is checked
  # against, so a `.propensity` of a different length has to be caught first.
  # Otherwise a score that matches `.propensity` is reported against the
  # exposure's count, which describes the score rather than the mismatch that
  # made it wrong. The categorical route checks the two lengths against each
  # other in `calculate_categorical_weights()`.
  if (exposure_type %in% c("binary", "continuous")) {
    check_lengths_match(.propensity, .exposure, call = call)
  }

  # The score is checked again by `psw()`, which is the choke point for every
  # route that records one. It is checked here as well because the formulas
  # below multiply it into the weights first, and a length that neither matches
  # nor divides them recycles under a base R warning before the constructor is
  # reached. A score is only checked when it is used; `stabilize = FALSE`
  # ignores it, as documented.
  if (isTRUE(stabilize)) {
    stabilization_score <- check_stabilization_score(
      stabilization_score,
      length(.exposure),
      call = call
    )
  }

  if (exposure_type == "binary") {
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- ate_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      stabilize = stabilize,
      stabilization_score = stabilization_score,
      call = call
    )
  } else if (exposure_type == "continuous") {
    wts <- ate_continuous(
      .propensity = .propensity,
      .exposure = .exposure,
      .sigma = .sigma,
      .density = .density,
      numerator = numerator,
      stabilize = stabilize,
      stabilization_score = stabilization_score,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    # including the more specific matrix checks
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "ate",
      .focal_level = NULL,
      stabilize = stabilize,
      stabilization_score = stabilization_score,
      call = call
    )
  }

  # Create psw object with appropriate attributes. The stabilization score is
  # only recorded when the user supplied one and stabilization is on; the
  # default marginal stabilizer stores nothing so downstream code can tell the
  # cases apart.
  psw_obj <- psw(
    wts,
    "ate",
    stabilized = isTRUE(stabilize),
    stabilization_score = if (isTRUE(stabilize)) stabilization_score else NULL
  )

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_ate.data.frame <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_ate.numeric,
    fn_name = "wt_ate",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
    .propensity_col_quo = col_quo,
    .sigma = .sigma,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ate.glm <- function(
  .propensity,
  .exposure = NULL,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_ate",
    call = call
  )

  # Call the numeric method
  wt_ate.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    .sigma = .sigma,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    call = call,
    ...
  )
}


ate_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  if (isTRUE(stabilize) && is.null(stabilization_score)) {
    p1 <- ipw_default_stab_seed(.exposure)
    p0 <- 1 - p1

    wt <- (.exposure * p1 / .propensity) +
      ((1 - .exposure) * p0 / (1 - .propensity))
  } else if (isTRUE(stabilize) && !is.null(stabilization_score)) {
    wt <- (.exposure / .propensity) + ((1 - .exposure) / (1 - .propensity))
    wt <- wt * stabilization_score
  } else {
    wt <- (.exposure / .propensity) + ((1 - .exposure) / (1 - .propensity))
  }

  wt
}

ate_continuous <- function(
  .propensity,
  .exposure,
  .sigma,
  .density = dens_normal(),
  numerator = "marginal",
  stabilize = FALSE,
  stabilization_score = NULL,
  call = rlang::caller_env()
) {
  # The conditional density f_{A|X}(A_i | X_i) is the one the propensity model
  # estimates, so its spread is that model's residual standard deviation: the
  # observation-level `.sigma` when the caller supplies one, the pooled residual
  # standard deviation otherwise.
  sigma_i <- continuous_sigma(.exposure, .propensity, .sigma)

  # What divides the conditional density: nothing at all when the weights are
  # not stabilized, a score the caller supplied when there is one, and otherwise
  # the numerator asked for. The combinations `check_numerator()` refuses cannot
  # reach here, so the two that answer for themselves are read first.
  numerator <- if (!isTRUE(stabilize)) {
    "none"
  } else if (!is.null(stabilization_score)) {
    "score"
  } else {
    numerator
  }

  # The marginal density f_A(A_i) is read at the exposure's population moments,
  # which no other numerator needs and a constant exposure has no spread for.
  mu_a <- NULL
  sigma_a <- NULL
  if (identical(numerator, "marginal")) {
    mu_a <- mean(.exposure, na.rm = TRUE)
    sigma_a <- sqrt(mean((.exposure - mu_a)^2, na.rm = TRUE))
  }

  wt <- continuous_density_ratio(
    exposure = .exposure,
    mu = .propensity,
    sigma = sigma_i,
    density = .density,
    numerator = numerator,
    mu_a = mu_a,
    sigma_a = sigma_a,
    score = stabilization_score,
    call = call
  )

  if (identical(numerator, "none")) {
    alert_info(
      "Using unstabilized weights for continuous exposures is not recommended."
    )
  }

  # What the weights are a ratio of, recorded on them for a reader to inspect
  # and for `ipw()` to rebuild the same weights from. The record travels back
  # with the weights rather than beside them so that it reaches the `psw` the
  # same way the categorical attributes do.
  attr(wt, "density_meta") <- new_density_meta(
    density = .density,
    numerator = numerator,
    sigma = if (is.null(.sigma)) "pooled" else "supplied"
  )

  wt
}


#' @export
#' @rdname wt_ate
wt_att <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_att")
}

#' @export
wt_att.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_att.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_att",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical"),
    announce = !be_quiet(),
    call = call
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure, call = call)
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- att_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "att",
      .focal_level = .focal_level,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = call
    )
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "att")

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  psw_obj <- record_binary_focal_level(
    psw_obj,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_att.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_att.numeric,
    fn_name = "wt_att",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_att.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_att",
    call = call
  )

  # Call the numeric method
  wt_att.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    call = call,
    ...
  )
}

att_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  tilted_binary_weights(.propensity, .exposure, "att")
}

#' @export
#' @rdname wt_ate
wt_atu <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_atu")
}

#' @export
wt_atu.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_atu.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_atu",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical"),
    announce = !be_quiet(),
    call = call
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure, call = call)
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- atu_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atu",
      .focal_level = .focal_level,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = call
    )
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "atu")

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  psw_obj <- record_binary_focal_level(
    psw_obj,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_atu.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atu.numeric,
    fn_name = "wt_atu",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_atu",
    call = call
  )

  # Call the numeric method
  wt_atu.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    call = call,
    ...
  )
}

atu_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  tilted_binary_weights(.propensity, .exposure, "atu")
}

#' @export
#' @rdname wt_ate
wt_atm <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_atm")
}

#' @export
wt_atm.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_atm.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_atm",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical"),
    announce = !be_quiet(),
    call = call
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure, call = call)
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- atm_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "atm",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = call
    )
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "atm")

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_atm.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_atm.numeric,
    fn_name = "wt_atm",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_atm",
    call = call
  )

  # Call the numeric method
  wt_atm.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    call = call,
    ...
  )
}

atm_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  tilted_binary_weights(.propensity, .exposure, "atm")
}


#' @export
#' @rdname wt_ate
wt_ato <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_ato")
}

#' @export
wt_ato.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_ato.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_ato",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical"),
    announce = !be_quiet(),
    call = call
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure, call = call)
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- ato_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "ato",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = call
    )
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "ato")

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_ato.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_ato.numeric,
    fn_name = "wt_ato",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_ato",
    call = call
  )

  # Call the numeric method
  wt_ato.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    call = call,
    ...
  )
}

ato_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  tilted_binary_weights(.propensity, .exposure, "ato")
}

#' @export
#' @rdname wt_ate
wt_entropy <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_entropy")
}

#' @export
wt_entropy.default <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_entropy.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_entropy",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "categorical"),
    announce = !be_quiet(),
    call = call
  )

  if (exposure_type == "binary") {
    check_lengths_match(.propensity, .exposure, call = call)
    .propensity <- resolve_binary_propensity(
      .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    check_ps_range(.propensity, call = call)
    wts <- entropy_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
  } else {
    # For categorical, let calculate_categorical_weights handle all validation
    wts <- calculate_categorical_weights(
      ps_matrix = .propensity,
      .exposure = .exposure,
      estimand = "entropy",
      .focal_level = NULL,
      stabilize = FALSE,
      stabilization_score = NULL,
      call = call
    )
  }

  # Create psw object with appropriate attributes
  psw_obj <- psw(wts, "entropy")

  # Preserve categorical attributes if they exist
  psw_obj <- preserve_categorical_attrs(psw_obj, wts, exposure_type)

  record_exposure_attrs(psw_obj, wts, exposure_type)
}

#' @export
#' @rdname wt_ate
wt_entropy.data.frame <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_entropy.numeric,
    fn_name = "wt_entropy",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .propensity_col_quo = col_quo,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.glm <- function(
  .propensity,
  .exposure = NULL,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "categorical"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_entropy",
    call = call
  )

  # Call the numeric method
  wt_entropy.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    call = call,
    ...
  )
}

entropy_binary <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  tilted_binary_weights(.propensity, .exposure, "entropy")
}

# Every binary weight but the ate's is the estimand's tilt h(e) divided by the
# propensity score of the exposure level the unit received: h(e)/e for the focal
# level, h(e)/(1 - e) for the reference level. The ate's tilt is the constant 1,
# and its formula carries the stabilization arms instead.
tilted_binary_weights <- function(.propensity, .exposure, estimand) {
  ps_tilt_binary(.propensity, estimand) /
    (.exposure * .propensity + (1 - .exposure) * (1 - .propensity))
}

# --------------------------------------------------------------------
#  wt_atc() - alias for wt_atu()
# --------------------------------------------------------------------

#' @export
#' @rdname wt_ate
wt_atc <- wt_atu

# --------------------------------------------------------------------
#  wt_cens() - censoring weights
# --------------------------------------------------------------------

#' @export
#' @rdname wt_ate
wt_cens <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("wt_cens")
}

#' @export
wt_cens.default <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  abort_no_method(.propensity)
}

#' @export
wt_cens.numeric <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env(),
  user_env = rlang::caller_env()
) {
  check_call_arg(call)

  # Resolved here rather than left to the ATE machinery, which would name
  # `wt_ate()` in the deprecation warning for an argument the user passed to
  # `wt_cens()`.
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "wt_cens",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  # Censoring is a two-level event, so the exposure type is resolved against the
  # types censoring weights answer rather than against the wider set the ATE
  # machinery answers.
  exposure_type <- causalgenerics::match_exposure_type(
    exposure_type,
    .exposure,
    c("auto", "binary", "continuous"),
    announce = !be_quiet(),
    call = call
  )

  # Get weights using ATE formula. `call` carries the frame this method was
  # dispatched into, so a rejection raised inside the ATE machinery still names
  # `wt_cens()`.
  wts <- wt_ate.numeric(
    .propensity = .propensity,
    .exposure = .exposure,
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = NULL,
    .untreated = NULL,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    call = call,
    ...
  )

  # Update estimand to "uncensored"
  estimand(wts) <- "uncensored"
  wts
}

#' @export
#' @rdname wt_ate
wt_cens.data.frame <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  col_quo <- rlang::enquo(.propensity_col)
  handle_data_frame_weight_calculation(
    weight_fn_numeric = wt_cens.numeric,
    fn_name = "wt_cens",
    .propensity = .propensity,
    .exposure = .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "continuous"),
    .propensity_col_quo = col_quo,
    .sigma = .sigma,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

#' @export
wt_cens.glm <- function(
  .propensity,
  .exposure = NULL,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL,
  call = rlang::current_env()
) {
  check_call_arg(call)

  args <- prepare_glm_weight_args(
    .propensity,
    .exposure,
    exposure_type = exposure_type,
    valid_exposure_types = c("auto", "binary", "continuous"),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "wt_cens",
    call = call
  )

  # Call the numeric method
  wt_cens.numeric(
    .propensity = args$propensity,
    .exposure = args$exposure,
    .sigma = .sigma,
    exposure_type = args$exposure_type,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    call = call,
    ...
  )
}

# --------------------------------------------------------------------
#  methods for `ps_trim()`
# --------------------------------------------------------------------

#' @export
wt_ate.ps_trim <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trim",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

#' @export
wt_att.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ate.ps_trunc <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "trunc",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

#' @export
wt_att.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trim",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "trunc",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atc.ps_trim <- wt_atu.ps_trim

#' @export
wt_atc.ps_trunc <- wt_atu.ps_trunc

#' @export
wt_cens.ps_trim <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trim",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

#' @export
wt_cens.ps_trunc <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "trunc",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

# --------------------------------------------------------------------
#  Categorical exposure weight calculations
# --------------------------------------------------------------------

calculate_categorical_weights <- function(
  ps_matrix,
  .exposure,
  estimand,
  .focal_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  call = rlang::caller_env()
) {
  # Ensure exposure is a factor
  .exposure <- transform_exposure_categorical(
    .exposure,
    .focal_level,
    call = call
  )

  # Validate propensity score matrix
  ps_matrix <- check_ps_matrix(ps_matrix, .exposure, call = call)

  # Get dimensions
  n <- length(.exposure)
  k <- nlevels(.exposure)
  levels_exp <- levels(.exposure)

  # Create indicator matrix for exposures
  # Each row i has a 1 in column j if observation i has exposure level j
  Z <- matrix(0, nrow = n, ncol = k)
  for (j in 1:k) {
    Z[.exposure == levels_exp[j], j] <- 1
  }

  # Extract the propensity score for each unit's actual exposure
  # e_{i,Z_i} = sum over j of Z_{ij} * e_{ij}
  e_actual <- rowSums(Z * ps_matrix)

  # A unit with no observed level has an all-zero indicator row, so e_{i,Z_i}
  # comes back as 0. There is no weight to give such a unit, so it is missing
  # here just as it is on the binary and continuous paths.
  missing_exposure <- is.na(.exposure)
  e_actual[missing_exposure] <- NA_real_

  if (!estimand %in% ipw_estimands) {
    abort(
      "Unknown estimand: {estimand}",
      error_class = "propensity_unknown_estimand_error",
      call = call
    )
  }

  # The tilting function h(e_i) reads the focal column for ATT and ATU and
  # treats every level alike otherwise.
  focal_idx <- NULL
  if (estimand %in% c("att", "atu")) {
    check_focal_level_required(.focal_level, estimand, call = call)
    focal_idx <- which(levels_exp == .focal_level)
  }
  h_e <- ps_tilt_categorical(ps_matrix, estimand, focal_idx)

  # Calculate weights: w_i = h(e_i) / e_{i,Z_i}
  weights <- h_e / e_actual

  # Apply stabilization only for ATE
  if (isTRUE(stabilize)) {
    if (estimand != "ate") {
      abort(
        "Stabilization is only supported for ATE with categorical exposures.",
        error_class = "propensity_stabilize_categorical_error",
        call = call
      )
    }

    if (!is.null(stabilization_score)) {
      weights <- weights * stabilization_score
    } else {
      # Every marginal is a share of the units with an observed level, so with
      # none observed each one is 0 / 0. Weights that are missing everywhere are
      # not an answer to the call that was made, so the degenerate stabilizer is
      # reported rather than returned.
      if (all(missing_exposure)) {
        abort(
          c(
            "Can't stabilize categorical weights when {.arg .exposure} has no
             observed values.",
            x = "Every observation is missing, so every marginal probability
                 is undefined.",
            i = "Supply {.arg stabilization_score}, or leave
                 {.code stabilize = FALSE}."
          ),
          error_class = "propensity_stabilize_categorical_error",
          call = call
        )
      }

      # For categorical, use marginal probabilities. `table()` counts only the
      # units with an observed level, so the denominator has to match or the
      # marginals sum to less than 1 whenever the exposure is missing.
      p_marginal <- table(.exposure) / sum(!missing_exposure)

      # Create stabilization weights based on marginal probabilities
      stab_wts <- rep(NA_real_, n)
      for (j in 1:k) {
        stab_wts[.exposure == levels_exp[j]] <- p_marginal[j]
      }

      weights <- weights * stab_wts
    }
  }

  # Add attributes for categorical weights
  attr(weights, "n_categories") <- k
  attr(weights, "category_names") <- levels_exp
  if (!is.null(.focal_level)) {
    attr(weights, "focal_category") <- .focal_level
  }

  weights
}

# ps_calib methods ----

#' @export
wt_ate.ps_calib <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ate.numeric,
    modification_type = "calib",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}

#' @export
wt_att.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_att.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atu.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atu.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atm.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_atm.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_ato.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_ato.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_entropy.ps_calib <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_entropy.numeric,
    modification_type = "calib",
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
wt_atc.ps_calib <- wt_atu.ps_calib

#' @export
wt_cens.ps_calib <- function(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
) {
  calculate_weight_from_modified_ps(
    .propensity = .propensity,
    .exposure = .exposure,
    weight_fn = wt_cens.numeric,
    modification_type = "calib",
    .sigma = .sigma,
    exposure_type = exposure_type,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density,
    numerator = numerator,
    ...
  )
}
