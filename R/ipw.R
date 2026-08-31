#' Inverse Probability Weighted Estimation
#'
#' @description
#' `ipw()` is a bring-your-own-model (BYOM) inverse probability weighted
#' estimator for causal inference. You supply a fitted propensity score model
#' and a fitted weighted outcome model; `ipw()` computes causal effect estimates
#' with standard errors that correctly account for the two-step estimation
#' process.
#'
#' `ipw()` supports binary, categorical, and continuous exposures. For a binary
#' or categorical exposure, a binary outcome returns the risk difference, log
#' risk ratio, and log odds ratio, and a continuous outcome returns the
#' difference in means. A continuous exposure is supplied through an
#' [stats::lm()] or gaussian [stats::glm()] propensity score model with a
#' weighted marginal structural outcome model whose link is identity, logit, or
#' log, and `ipw()` reports the single exposure coefficient of that model. A
#' subclass of either propensity score model, such as a robust or an additive
#' fit, errors.
#'
#' @param wt_mod The weighting object: a fitted propensity score model that
#'   produced the weights, typically a logistic regression of class
#'   [stats::glm()] with the exposure as the left-hand side of the formula.
#'   `ipw()` is an S3 generic that dispatches on this object. The left-hand side
#'   must be the exposure column itself, for every exposure type: a matrix
#'   response such as `cbind(successes, failures)` and a transformed response
#'   such as `log(a)` or `factor(z)` both error, because `ipw()` reads the
#'   exposure's name from that expression. Compute a transformed exposure into
#'   its own column and fit `wt_mod` on that column.
#' @param outcome_mod A fitted weighted outcome model of class [stats::glm()]
#'   or [stats::lm()], with the outcome as the dependent variable and
#'   propensity score weights supplied via the `weights` argument. The weights
#'   should be created with a propensity weight function such as [wt_ate()].
#'   Supported outcome models are an [stats::lm()], a gaussian
#'   [stats::glm()] with an identity link, and a binomial or quasibinomial
#'   [stats::glm()] with a logit, probit, cloglog, log, or identity link; any
#'   other family (such as poisson or Gamma), an unsupported link (such as
#'   cauchit), or a non-identity gaussian link errors. A factor or logical
#'   outcome response is
#'   converted to `0`/`1` following glm's coding (the first factor level is the
#'   failure, every other level is the success).
#' @param .data A data frame containing the exposure, outcome, and covariates.
#'   If `NULL` (the default), `ipw()` extracts data from the model objects.
#'   Supply `.data` explicitly if the outcome model formula contains
#'   transformations that prevent extraction of the exposure variable from
#'   [stats::model.frame()], or if the propensity score model cannot reconstruct
#'   its design (for example, fit with `model = FALSE` with the fitting data no
#'   longer available); `ipw()` then rebuilds the propensity design from `.data`.
#'   `.data` must be the data the models were fit to. It must have one row per
#'   observation the models were fit to, and each column must carry the type its
#'   model was fit on: a column fit as a factor may also be supplied as character
#'   or as an ordered factor, since every design rebuilt from `.data` is rebuilt
#'   under the levels the fit recorded, but a factor supplied as numeric, or a
#'   numeric supplied as a factor or a logical, errors. A row count or a column
#'   type that disagrees with the fitted models errors rather than rebuilding a
#'   design the models were not fit to.
#'
#'   The levels a categorical column declares are its fit's to decide, so a
#'   column re-leveled after fitting is rebuilt in the fitted order, and one that
#'   declares a level no observation carries rebuilds the fitted design, since
#'   the fits drop unused levels themselves. A value the fit never saw has no
#'   coefficient to multiply and errors. The exposure column and the outcome
#'   model's response are read rather than rebuilt, and both are rejected when
#'   the order they declare contradicts the fit: `ipw()` treats the second level
#'   of a binary exposure as the exposed group, and codes a factor response as an
#'   indicator for its non-first levels, so either one re-leveled the other way
#'   describes the opposite contrast.
#'
#'   Supplying `.data` also decides how a term that transforms the exposure is
#'   handled. The counterfactual designs are built by setting the exposure column
#'   to each level in turn and rebuilding the outcome design from `.data`, so a
#'   term such as `as.numeric(a == "c")`, `I(z^2)`, or `factor(z)` is
#'   re-evaluated at the level being set rather than held at its observed value,
#'   and such a model is g-computation on the model as specified. Without
#'   `.data` there is nothing to re-evaluate: the designs come from the outcome
#'   model's own model frame, which holds each such term at the values it was fit
#'   on. A term that has to derive the exposure's levels from the values it sees,
#'   such as `factor(a)`, is therefore rejected on that route, with `.data` named
#'   as the remedy.
#' @param .by A single unquoted variable naming a modifier to report effects
#'   within, or `NULL` (the default) to report the effect for the sample as a
#'   whole and nothing else. The variable is looked up in `.data` when you
#'   supply one and in the outcome model's own model frame otherwise, and it
#'   must be a factor or a character vector with no missing values. The result
#'   then reports the effects it reports without `.by`, keyed by the group
#'   `"overall"`, followed by one set within each level of the modifier and one
#'   set contrasting each level against the first, and the `estimates` frame
#'   gains a `group` column naming the subgroup each row belongs to. See
#'   **Effect modification** for what those rows estimate and what they require.
#' @param estimand A character string specifying the causal estimand: one of
#'   `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`. The available
#'   estimands depend on the exposure type: a binary or categorical exposure
#'   supports all six, while a continuous exposure supports only `"ate"`. For a
#'   categorical exposure, `"att"` and `"atu"` require a focal level (see
#'   `.focal_level`). See the **Estimands** section for the full support matrix.
#'   If `NULL`, the estimand is inferred from the weights in `outcome_mod`, which
#'   requires weights created with a propensity weight function such as
#'   [wt_ate()].
#' @param ps_link `r lifecycle::badge("deprecated")` A character string
#'   specifying the link function used in the propensity score model: `"logit"`,
#'   `"probit"`, or `"cloglog"`. The link is always read from `wt_mod`, so
#'   `ps_link` can only restate what the fitted model already supplies and can
#'   simply be dropped. Supplying it warns.
#'
#'   The argument applies only to a binomial [stats::glm()] propensity score
#'   model on the binary path. A multinomial or continuous propensity score
#'   model has no link for `ps_link` to override, so a non-`NULL` value errors
#'   on both of those paths; leave it `NULL` there. `ps_link` cannot name a link
#'   other than the one `wt_mod` was fit with: `se_method = "linearization"`
#'   and `se_method = "robust"` error, and `se_method = "mestimation"` reports
#'   the resulting weights as inconsistent with the propensity score model.
#' @param conf_level Confidence level for intervals. Default is `0.95`.
#' @param se_method Method for standard error estimation. `"mestimation"` (the
#'   default) stacks the propensity score, outcome, and estimand estimating
#'   equations and returns the empirical sandwich variance. `"linearization"`
#'   uses the influence-function linearization of Kostouraki et al. (2024).
#'   Both account for the uncertainty of estimating the propensity scores.
#'
#'   `"robust"` does not, and is a diagnostic rather than an inference. It is
#'   the linearization route with the correction for having estimated the
#'   propensity score dropped, so it reports the sandwich the weighted outcome
#'   model computes for itself: the delta-method reading of
#'   `sandwich::vcovHC(outcome_mod, type = "HC0")`, exactly. Ignoring that the
#'   weights were estimated generally understates the variance, and understates
#'   it most where the weights are least stable, so this is offered for reading
#'   beside one of the other two rather than in place of it. A result fit this
#'   way says so when printed and marks the tables it coerces to; see
#'   **Standard errors as a diagnostic** below.
#'
#'   Which exposure types support which: `"mestimation"` supports every one of
#'   them; `"linearization"` and `"robust"` support a binary exposure alone. A
#'   binary or a categorical exposure has a sandwich for every fit it accepts.
#'   A few
#'   continuous fits have no sandwich, and the package computes no standard
#'   error for them: an [mgcv::gam()] propensity score model and weights built
#'   with a `"kernel"` density are refused, and so is a [MASS::rlm()] fit the
#'   stacked system cannot write. See **Continuous propensity score models**
#'   below.
#'
#'   `"linearization"` and `"robust"` support only an outcome model of the
#'   exposure alone, fit with an intercept; a covariate-adjusted outcome model
#'   requires `"mestimation"`. For a binary or categorical exposure, every
#'   method requires an outcome model that can represent a baseline, so a
#'   numeric no-intercept coding such as `y ~ z - 1` errors on any of them. The
#'   standard errors of the conditional reading described under `effects`
#'   require `"mestimation"` as well: the linearization and robust routes store
#'   the outcome model unwrapped, so [stats::vcov()], [stats::confint()], and
#'   [`tidy()`][tidy.ipw()] refuse `effects = "conditional"` there with a
#'   classed error rather than reporting the covariance the outcome model
#'   computed for itself.
#' @param effects The presentation mode the result records, either `"marginal"`
#'   (the default) or `"conditional"`. The marginal reading reports the
#'   population-averaged causal contrasts; the conditional reading reports the
#'   outcome model's coefficient surface. Both surfaces are computed whichever
#'   mode is named, so the argument settles which one the result presents and
#'   nothing else. [causalgenerics::as_marginal()] and
#'   [causalgenerics::as_conditional()] move a result between the two readings
#'   afterwards, and the accessors take an `effects` argument of their own for a
#'   single call.
#'
#'   One kind of fit has a single reading rather than two. A continuous exposure
#'   entering `outcome_mod` through more than one design column, such as a
#'   `poly()` or spline dose response, has no coefficient that is the effect of
#'   the exposure, because a curve has a different slope at every dose. Such a
#'   result records the conditional reading and says so once. `ipw()` refuses
#'   `effects = "marginal"` itself with an error of class
#'   `propensity_ipw_effects_error`, and the result declares the conditional
#'   reading as the only one it supports, so the accessors refuse it from the
#'   result class with an error of class
#'   `causalgenerics_unsupported_reading_marginal`. See the continuous exposure
#'   section below.
#'
#'   The covariance the conditional reading reports is the outcome block of the
#'   jointly estimated sandwich, which every route that stacks estimating
#'   equations attaches to the outcome model it stores:
#'   `se_method = "mestimation"` for a binary exposure, the categorical and
#'   joint routes, which run on M-estimation alone, and the continuous route
#'   under that method. A linearization or robust fit stacks no system and
#'   records no covariance of either kind, so [stats::vcov()] and
#'   [stats::confint()] error on its conditional reading, and printing it writes
#'   the coefficients under a note saying the standard errors are not reported
#'   rather than beside the ones the outcome model computed for itself.
#' @param ... These dots are for future extensions and must be empty. They
#'   separate the two models from the remaining arguments, which must therefore
#'   be supplied by name. Anything left in them, such as a `.data` argument
#'   given by position or a misspelled argument name, errors.
#'
#' @details
#' # Workflow
#'
#' `ipw()` is designed around a three-step workflow:
#'
#' 1. Fit a propensity score model (e.g., logistic regression of exposure on
#'    confounders).
#' 2. Calculate propensity score weights for your estimand (e.g., [wt_ate()])
#'    and fit a weighted outcome model.
#' 3. Pass both models to `ipw()` to obtain causal effect estimates with
#'    correct standard errors.
#'
#' You are responsible for specifying and fitting both models. `ipw()` handles
#' the variance estimation.
#'
#' # Estimands
#'
#' The available estimands depend on the exposure type:
#'
#' | Estimand    | Binary | Categorical             | Continuous |
#' |-------------|:------:|:-----------------------:|:----------:|
#' | `"ate"`     | yes    | yes                     | yes        |
#' | `"att"`     | yes    | yes (needs focal level) | no         |
#' | `"atu"`     | yes    | yes (needs focal level) | no         |
#' | `"atm"`     | yes    | yes                     | no         |
#' | `"ato"`     | yes    | yes                     | no         |
#' | `"entropy"` | yes    | yes                     | no         |
#'
#' A continuous exposure supports only `"ate"`; any other request errors. For a
#' categorical exposure, `"att"` and `"atu"` require a focal level, supplied
#' through `.focal_level` or taken from the weights (see the `multinom` method).
#'
#' A binary exposure takes no focal level: `ipw()` always treats the second
#' sorted exposure level as the exposed group, which is the second level of a
#' factor, `TRUE` for a logical, and the larger of two numeric values. Because
#' `"att"` and `"atu"` weights are mirror images of each other, weights built
#' with [wt_att()] or [wt_atu()] on the other level carry values that are
#' correct for the reversed roles, and `ipw()` rejects them rather than
#' correcting them as though the roles matched. To target the other level,
#' relevel the exposure so that it sorts second and refit both models.
#'
#' The estimand a weighted outcome model was fit for is not necessarily the one a
#' tool that averages its per-row predictions reports, so see **Standardizing
#' model predictions** in [ps_tilt()] before taking such an average of a model fit
#' for anything other than `"ate"`.
#'
#' # Effect measures
#'
#' A binary or categorical exposure reports the counterfactual mean at each
#' exposure level before the contrasts built from those means. Those rows carry
#' the effect label `"mean"` and lead the table, one per level in level order
#' with the reference level first. Each is the outcome model's prediction with
#' every unit set to that level, averaged over the sample under the estimand's
#' tilt, so a binomial outcome reports the marginal risk under each level and a
#' continuous one reports the marginal mean. Reading a risk difference beside
#' the two risks it is a difference of is what those rows are for: a difference
#' of 0.05 is a different finding at risks of 0.02 and 0.07 than at 0.50 and
#' 0.55.
#'
#' For a binary exposure, the reported measures depend on the outcome model. A
#' binary outcome ([stats::glm()] with `family = binomial()`) returns three
#' measures:
#' - `rd`: Risk difference (marginal risk in exposed minus unexposed)
#' - `log(rr)`: Log risk ratio
#' - `log(or)`: Log odds ratio
#'
#' A continuous outcome ([stats::lm()] or [stats::glm()] with
#' `family = gaussian()`) returns only the difference in means (`diff`).
#'
#' For a categorical exposure, the same measures are reported for each
#' non-reference level against the reference (first) factor level, so a K-level
#' exposure produces one block of measures per non-reference level.
#'
#' The estimates table of either exposure carries a `contrast` column naming
#' what each row is about: the level a mean row belongs to, and the pair of
#' levels a contrast row compares, written `"1 vs 0"` or `"b vs a"` in the
#' levels the models were fit on. A binary exposure names its one pair the same
#' way a categorical one names each of its several, so the two tables are read
#' alike.
#'
#' For a continuous exposure, `ipw()` reports the exposure coefficients of the
#' weighted marginal structural outcome model. A model with one exposure
#' coefficient reports one row and no `contrast` column, and a dose-response
#' curve such as `y ~ A + I(A^2)` reports one row per coefficient, with the
#' `contrast` column naming the coefficient each row belongs to.
#'
#' The `effect` label names the scale: `log(or)` for a logit link and `log(rr)`
#' for a log link, since a coefficient there is a log odds ratio or a log risk
#' ratio per unit of whatever column it multiplies. At an identity link the word
#' depends on what the row claims. A model with one exposure coefficient reports
#' `slope`, because that coefficient is the slope of the dose response
#' everywhere. A model with several reports `coef`, because a curve has a
#' different slope at every dose and no one of its coefficients is that slope;
#' the row says which coefficient it is and leaves the dose response to be built
#' downstream from `coef()` and `vcov()`.
#'
#' Every term of the marginal structural model that reads the exposure must read
#' the exposure and nothing else. The requirement is on what a term reads rather
#' than on how it is written, so `A + I(A^2)` and `A + sin(A)` are admitted and
#' `A * x1` is not: a term reading a covariate as well contributes a coefficient
#' that depends on that covariate, and no row could name the effect it stands
#' for. Read the whole coefficient vector from the returned `fit` object for a
#' model this surface cannot report.
#'
#' A basis matrix reads the exposure alone, so `poly(A, 2)`,
#' `splines::ns(A, 3)`, and `splines::bs(A, 3)` are admitted on the same
#' footing. One such term contributes one coefficient per basis column and each
#' of them is a row, reported under `coef` at an identity link like any other
#' multi-coefficient surface. The `contrast` column names the coefficient rather
#' than the term, a distinction a basis makes and `A + I(A^2)` does not: the
#' single term `poly(A, 2)` reports the contrasts `poly(A, 2)1` and
#' `poly(A, 2)2`. A basis reading a covariate is a covariate term and
#' contributes no row, and one reading both, such as `poly(A, 2):x1`, is refused
#' with every other mixed term.
#'
#' A marginal structural model built on a basis of the exposure requires
#' `.data`. A model frame records the term rather than the variables inside
#' it, so the frame of such a fit carries the basis matrix and no column
#' holding the exposure, which is where `ipw()` reads the exposure from when
#' `.data` is absent. Omitting it errors, and the message names the exposure
#' and directs to `.data`. Supplying it does not refit the basis: the outcome
#' design is read from the fit, and a design rebuilt from `.data` is rebuilt
#' through the fitted terms object, whose `predvars` attribute records the
#' basis that term was fit with.
#'
#' Use [`as.data.frame()`][causalgenerics::new_ipw()] with
#' `exponentiate = TRUE` to obtain risk ratios and odds ratios on their natural
#' scale.
#'
#' # Effect modification
#'
#' `.by` names a modifier and asks for the effect within each of its levels
#' alongside the effect for the sample as a whole. The reported rows come in
#' four blocks, in this order:
#'
#' - the overall rows, unchanged from the fit without `.by`, the counterfactual
#'   means and then the contrasts, under the group label `"overall"`;
#' - the counterfactual mean at each exposure level within each level of the
#'   modifier, under a `"var = value"` label such as `"sex = female"`;
#' - one block of contrasts per level of the modifier, under the same labels;
#' - one block of contrasts per non-reference level against the reference level,
#'   which is the modifier's first level, under a
#'   `"var = value vs var = value"` label.
#'
#' The last block reports no means of its own. Its rows compare the effects in
#' two strata, and a difference of two effects belongs to neither stratum's
#' population, so there is no level whose mean it would be.
#'
#' A stratum's effect is g-computation restricted to that stratum: the outcome
#' model predicts each unit's outcome at each exposure level, and the
#' predictions are averaged over the units of the stratum, weighted by the
#' estimand's tilt. What a stratum row estimates is therefore the effect in that
#' stratum's own covariate distribution, and the last block is a difference of
#' two such effects. Those two populations differ whenever the covariates are
#' distributed differently across the strata, so a difference between strata is
#' a difference between two conditional effects and not, on its own, evidence
#' that the modifier changes the effect. It is the estimand's population that
#' the tilt fixes, and the modifier's stratum that the indicator fixes, and the
#' row reports the effect where those two meet: an `att` fit reports the effect
#' among the exposed units of a stratum, not among all of its units. For a
#' categorical exposure the tilt an `att` fit standardizes over is the focal
#' level's own propensity score, so a stratum row there reports the effect
#' among the units of that stratum who look like the ones assigned the focal
#' level.
#'
#' A categorical exposure already reports one contrast per non-reference level,
#' and `.by` crosses those contrasts with the strata: every block holds each
#' contrast on each reported measure, so a K-level exposure with S strata
#' reports the whole-sample block and then S + (S - 1) blocks of (K - 1)
#' contrasts, beside S blocks of K means. A row is named by three things, the
#' measure, the contrast, and the subgroup, and the `estimates` frame carries
#' `contrast` and `group` in that order after `effect`. Within a contrast block
#' the ordering is the one an ungrouped fit already uses, contrast-major and
#' measure-minor, so a contrast's measures sit together; the blocks themselves
#' run group-major.
#'
#' The stratum rows and the stratum contrasts report the risk difference and the
#' log risk ratio for a binary outcome, and the difference in means for a
#' continuous one. They report no odds ratio. An odds ratio is noncollapsible,
#' so a difference of two of them moves with the outcome distribution in each
#' stratum whether or not the effect does, and reporting it as effect
#' modification would invite exactly the reading it does not support. The
#' overall rows keep all three measures.
#'
#' Three requirements apply. The modifier must be a factor or a character vector
#' with no missing values: a missing value names no subgroup, and a continuous
#' variable names no fixed set of them, so cut it into groups first. A level no
#' unit carries is dropped rather than refused, since the fits drop it too and
#' an empty stratum has no mean to standardize. Every one of the modifier's
#' remaining levels must hold every exposure level, since a stratum in which
#' nobody took one of them identifies no contrast against it and the outcome
#' model would extrapolate one from the strata that do; the remedy is a coarser
#' modifier rather than a refit. Finally, `.by` requires
#' `se_method = "mestimation"` and a binary or categorical exposure: the stratum
#' means and their contrasts are parameters of the stacked estimating equations,
#' which is what puts them in the same sandwich as everything else, and a
#' continuous exposure reports its marginal structural model's own coefficient
#' rather than a set of standardized means to modify. Each of those errors.
#'
#' An outcome model with no term reading the exposure and the modifier together
#' warns rather than errors. The stratum effects are g-computation on the model
#' as specified, so such a model forces the same effect on every stratum up to
#' the covariate distribution each one has, and the modification it reports is a
#' property of the model rather than of the data. The fit still returns.
#'
#' Stabilizing the weights on the modifier is worth doing here. A stabilized
#' weight divides by the propensity score and multiplies by a numerator, and the
#' numerator may condition on any variable the outcome model the effects are
#' computed from already reads, without changing what is estimated. The modifier
#' is such a variable by construction, so a numerator conditioning on it leaves
#' the estimator consistent for the same stratum effects and changes only how
#' variable the weights are. That is where it pays: the part of the exposure's
#' variation the modifier explains cancels between the numerator and the
#' propensity score, and the weights tighten by more than the default numerator
#' tightens them.
#'
#' The converse is the constraint, and it is the reason this is a statement
#' about the model rather than a free choice. A numerator conditioning on a
#' variable the outcome model does not read changes the estimand rather than the
#' precision, because the pseudo-population it builds is one in which that
#' variable still predicts the exposure. Condition the numerator on the modifier
#' and on nothing the model does not hold. The same rule governs an
#' exposure-by-modifier interaction reported without `.by`, since what matters
#' is what the reported model reads and not which argument asked for the rows.
#'
#' Build the numerator with [wt_ate()]'s `stabilization_score`, evaluated at the
#' exposure each unit actually took; see **Stabilization** in [wt_ate()] for the
#' recipe. The sandwich treats a supplied score as a known constant, so the
#' standard errors do not account for the numerator model having been fitted.
#' The default stabilizer is estimated in the stacked system instead, which is
#' one parameter wider as a result. A continuous exposure has a third way: pass
#' the fitted numerator model to `stabilize` itself, and its own equations join
#' the stack, which is what a supplied score has none of.
#'
#' A stabilized fit reporting effects within the strata of a modifier, and
#' carrying the default numerator, warns with class
#' `propensity_ipw_by_stabilizer_warning` and returns: the default numerator
#' conditions on nothing, which is consistent for the same stratum effects and
#' leaves precision on the table. A numerator the caller supplied, whether as a
#' score or as a model, says nothing here, since which variables it conditions
#' on is a statement about the reported model that the weights do not carry.
#'
#' Every reported row is a parameter of one stacked system solved on one sample,
#' so the covariance couples them: [stats::vcov()] reports nonzero entries
#' between the overall rows, the stratum rows, and the stratum contrasts, and
#' the standard error of a stratum contrast accounts for the two strata having
#' been estimated together. Fitting each stratum on its own subset would report
#' the same stratum effects with no covariance between them, leaving the
#' difference between two strata untestable.
#'
#' # Joint exposures
#'
#' [causalgenerics::joint_exposure()] crosses two discrete treatments into one
#' categorical exposure and records the crossing on the vector. The result is a
#' factor, so `ipw()` fits it exactly as it fits any exposure over the same
#' cells. What the declaration changes is what gets reported: instead of each
#' cell against the reference cell, the result is written in the two treatments.
#'
#' - One row per cell, under the effect label `"mean"`, carrying the
#'   counterfactual risk (or mean) had everyone been set to that cell. The cell
#'   is named in the `contrast` column and the group is `"overall"`.
#' - The simple effects: each treatment's effect within a fixed level of the
#'   other, with the `contrast` naming the treatment and the level it is set to
#'   (`"a: 1 vs 0"`) and the `group` naming the level the other treatment is
#'   held at (`"e = 1"`). These include the comparisons that are not against the
#'   reference cell, which the vs-reference reporting cannot express at all.
#' - The interaction: the difference between two of the first treatment's simple
#'   effects, keyed by the two levels of the second being compared
#'   (`"e = 1 vs e = 0"`). On the risk difference scale this is the additive
#'   interaction, and on the log risk ratio scale the multiplicative one, which
#'   is the log of the ratio of the two risk ratios.
#'
#' A 2-by-2 crossing therefore reports fourteen rows for a binary outcome and
#' nine for a continuous one. The interaction is reported once, under the first
#' treatment's framing, rather than twice: interaction is symmetric in the two
#' treatments, so the difference between the first treatment's two simple
#' effects is the difference between the second treatment's two, and reporting
#' both would be one quantity under two labels.
#'
#' No odds ratio is reported anywhere on this surface. It is noncollapsible, so
#' neither a simple effect reported beside one nor a difference of two of them
#' says what it appears to.
#'
#' Every row is a parameter of the same stacked system, so the covariance
#' couples them: the interaction covaries with each simple effect it is built
#' from, and with each cell mean underneath those. The cell means are the
#' marginal-mean block the categorical path already estimates, and the simple
#' effects and interaction rows are contrasts over that block, so an interaction
#' equals the corresponding double difference of the cell means exactly rather
#' than to within a solver residual.
#'
#' Two restrictions apply. A declared crossing is reported for the `"ate"`
#' estimand alone: every cell mean on this surface standardizes to one
#' population, and a focal estimand standardizes each of them to a population
#' the simple effects and the interaction are not defined over, so weights built
#' for anything else error. And `.by` cannot be combined with a declared
#' crossing: effect modification of a joint intervention is a three-way
#' question, the interaction between two treatments within the levels of a
#' third variable, which this surface neither answers nor sensibly projects.
#' Either restriction is escaped the same way, by dropping the declaration with
#' `factor(x)` and reporting each cell against the reference cell.
#'
#' ## Two treatment models instead of a declared crossing
#'
#' A joint intervention can also be weighted through the sequential
#' factorization it really has, \eqn{f(A | L) f(E | A, L)}: one treatment model
#' per treatment, recorded with [joint_wt_models()], and the product of their
#' weights, built with [wt_joint()]. Pass the container as `wt_mod` and an
#' outcome model that reads both treatments as separate columns.
#'
#' The reported surface is the one above, row for row and label for label: the
#' same cell means, simple effects, and interaction, under the same
#' `effect`, `contrast`, and `group` values. The two routes estimate the same
#' estimand through different models of the same propensity score, so they agree
#' as estimators rather than to the last bit, and they coincide when the two
#' parameterizations are saturated in the same covariates.
#'
#' Prefer the two-model route when the two treatments call for different
#' adjustment sets, or when the dependence of the second treatment on the first
#' is what you want to model directly, since each treatment then gets its own
#' model and its own covariates. Prefer the declared crossing when one model over
#' the cells is the natural specification. The factorization itself is validated
#' when the container is built rather than here: [joint_wt_models()] refuses a
#' second model that does not condition on the first treatment, and a pair that
#' condition on each other.
#'
#' Both treatment models are stacked alongside the outcome model, so the weights
#' entering the outcome score are rebuilt from both propensity score blocks on
#' every evaluation and the reported standard errors account for having
#' estimated both. The same two restrictions apply: the `"ate"` estimand only,
#' and no `.by`. Standard errors come from M-estimation alone.
#'
#' ## A joint intervention with a dose
#'
#' The second treatment may be a dose rather than a second discrete treatment,
#' recorded with any of the classes the single-treatment route reads a
#' conditional density from, and weighted with a stabilized density ratio, which
#' [wt_joint()] requires of a continuous component. A dose has no cells, so
#' there is no crossing to report counterfactual risks over: the surface is the
#' marginal structural model's own coefficients, which for an identity-link
#' model are exactly the weighted fit's coefficients.
#'
#' Which of two surfaces a fit reports is decided by the marginal structural
#' model alone. A model written in bare treatment terms reports the vocabulary
#' surface below, whose rows name the treatment each one varies and where it is
#' evaluated. Every other treatment-reading model reports the coefficient
#' surface, whose rows are named after the coefficients they report.
#'
#' ### The vocabulary surface
#'
#' For `y ~ a * e` with `a` binary and `e` a dose, three coefficients carry
#' causal content and each is a row:
#'
#' - the binary treatment's effect at a dose of zero, with the `contrast`
#'   naming the treatment and its levels (`"a: 1 vs 0"`) and the `group`
#'   naming where it is evaluated (`"e = 0"`);
#' - the dose's slope at the binary treatment's reference level
#'   (`"e: per unit"`, `"a = 0"`);
#' - the interaction, the change in the binary treatment's effect per unit of
#'   the dose, keyed by the one-unit step of the other treatment
#'   (`"a: 1 vs 0"`, `"e + 1 vs e"`).
#'
#' The `effect` column names the scale, keeping the vocabulary each treatment
#' type already uses: a level contrast is a `diff` and a per-unit change in the
#' dose is a `slope` under an identity link, and both are `log(or)` under a
#' logit link and `log(rr)` under a log link. An additive model, `y ~ a + e`,
#' forces one effect for the binary treatment at every dose and one slope at
#' either of its levels, so it reports two rows under the group `"overall"` and
#' no interaction.
#'
#' Those three readings hold of a model that is linear in each treatment and of
#' no other, which is why such a model reports this surface and no other model
#' does. The columns are read as well as the formula: a factor treatment under a
#' coding other than treatment contrasts leaves the terms bare while rescaling
#' or recentering the column, so its coefficients are no longer the effects
#' these rows name, and such a fit is refused rather than reported on either
#' surface. Refit the outcome model with the treatment as a 0/1 numeric or as an
#' unordered factor under treatment contrasts.
#'
#' ### The coefficient surface
#'
#' A marginal structural model carrying a transformed or a basis treatment term,
#' such as `y ~ a * sin(e)`, `y ~ a + e + I(e^2) + a:e`, or
#' `y ~ a * splines::ns(e, 3)`, reports one row per treatment-reading
#' coefficient. The `contrast` column carries the coefficient name exactly as
#' the fit writes it, and there is no `group` column at all: the surface makes
#' no claim about where a row is evaluated, because for a curve there is no one
#' place. Build an interpretable dose response downstream from `coef()` and
#' `vcov()`, whose covariance between the rows is what that takes.
#'
#' The `effect` label follows the outcome link, `coef` at an identity link and
#' the same `log(or)` and `log(rr)` as everywhere else. `coef` is a word of its
#' own rather than `diff` or `slope`, both of which carry an evaluated-at claim
#' on the vocabulary surface.
#'
#' ### Covariates and mixed terms
#'
#' A covariate entering on its own is admitted on either surface and contributes
#' no row, however many columns it expands to: it adjusts the marginal
#' structural model and the surface has nothing to say about its coefficient. A
#' term reading a treatment and a covariate together is refused on both, since
#' its coefficient is a change in an effect per unit of that covariate and no
#' row could name the effect it stands for. Nothing here is standardized over
#' the covariates, which is what the declared crossing does instead.
#'
#' Unlike the single-treatment route, this one needs no `.data` for a basis
#' marginal structural model. Each treatment is read off its own treatment
#' model, so the outcome model frame is never asked for a column it does not
#' hold.
#'
#' The dose model contributes its coefficients and the conditional variance its
#' density ratio divides by, and a marginal stabilizing numerator contributes
#' the dose's two marginal moments, all as parameters of the same stacked
#' system.
#'
#' ### The dose model and the ratio its weights are
#'
#' The dose model goes through the registry described under **Continuous
#' propensity score models**, so the classes and links this route stacks are the
#' ones the single-treatment route stacks: an [stats::lm()], a gaussian
#' [stats::glm()] fit with an identity or a log link, and a [MASS::rlm()] fit
#' with the Huber psi. A gaussian dose model fit with another link errors with
#' class `propensity_ipw_link_error`, and a robust fit whose score the system
#' cannot write raises the same conditions it raises for a single dose.
#'
#' The dose component's weights may be any ratio [wt_ate()] builds: a
#' heavier-tailed `.density`, `numerator = "integrated"`, or a single `.sigma`
#' the caller fixed. [wt_joint()] records each component's ratio on the product,
#' and the stacked system rebuilds the dose factor from that record, so what is
#' differentiated is the weight the outcome model was fit with. An integrated
#' numerator is built from the dose block and the data alone, so it adds no
#' stabilization parameters; a fixed scalar `.sigma` is a known constant, so the
#' dose block is its coefficients alone and none of that number's uncertainty is
#' carried. A spread supplied for each observation has no counterpart in the
#' system and errors with class `propensity_ipw_sigma_error`.
#'
#' A dose component built with a `stabilization_score` is the one numerator this
#' route cannot reproduce. The product records that the numerator was a score
#' without recording the vector it was, so the system rebuilds the dose weights
#' from the exposure's own marginal moments and the weight-consistency check
#' reports the difference. Weight the dose on its own to use a numerator
#' computed by hand. A dose stabilized on a fitted model is refused outright,
#' with class `propensity_ipw_numerator_error`: this route has no stabilization
#' block to estimate such a model in, where the single-dose route has one.
#'
#' An [mgcv::gam()] dose model, and dose weights built with a `"kernel"`
#' density, are refused with class
#' `propensity_ipw_se_method_unavailable_error`. Neither is a function of the
#' parameters that the sandwich could differentiate, and the package has no
#' second method to fall back on, so the refusal points at building the dose
#' weights from a model and a density the stacked system can write, or at
#' bootstrapping the whole fit by hand.
#'
#' # Variance estimation
#'
#' By default (`se_method = "mestimation"`), standard errors are computed by
#' M-estimation. The propensity score model, the weighted outcome model, the
#' marginal means, and the effect contrasts are stacked as a single system of
#' estimating equations, and the empirical sandwich variance of that joint
#' system is used. The stacked marginal means are standardized to the estimand's
#' tilted target population, so a non-`ate` estimand with a covariate-adjusted
#' outcome model reports the contrast for that target population rather than the
#' full-sample average. Stacking the models accounts for the uncertainty
#' introduced by estimating the propensity scores, avoiding the underestimated
#' standard errors that arise from treating estimated weights as fixed. See
#' Stefanski and Boos (2002) for the M-estimation framework.
#'
#' M-estimation standard errors are available for all exposure types: binary
#' (from a [stats::glm()] propensity score model), categorical (from a
#' [nnet::multinom()] model), and continuous (from an [stats::lm()], a gaussian
#' [stats::glm()], or a [MASS::rlm()] model; see **Continuous propensity score
#' models** below). The `atm` weight
#' `pmin(e, 1 - e)` is not differentiable at a propensity score of `0.5`; deli's
#' central finite difference straddles the kink there and averages the one-sided
#' slopes. Its effect on the variance is negligible unless many observations sit
#' at exactly `0.5`.
#'
#' M-estimation standard errors for a categorical exposure allocate memory
#' roughly linearly in the number of observations, on the order of 70 to 90
#' kilobytes per observation for a single fit, so a sample of 10,000
#' observations allocates on the order of hundreds of megabytes. The cost comes
#' from the stacked estimating-equation machinery rather than any single term,
#' and for very large samples you can expect long, garbage-collection-heavy
#' fits. This is expected behavior.
#'
#' ## Continuous propensity score models
#'
#' A continuous exposure's weights are a ratio of densities centered on the
#' conditional mean the propensity score model fits, so the stacked system
#' carries that model's own score. Three classes supply one: an [stats::lm()],
#' a gaussian [stats::glm()] fit with an identity or a log link, and a
#' [MASS::rlm()] fit with the Huber psi. A gaussian model fit with any other
#' link is refused by name, because the coefficients its iteration stops at are
#' not a tight enough root to seed the stacked solve from; refit it with one of
#' the two links or as an [stats::lm()]. Any other class errors, since solving
#' it here would report a propensity score model that was never fit.
#'
#' A [MASS::rlm()] is stacked as the Huber score its coefficients are the root
#' of, clipped where the fit itself clipped: at the psi's own constant, which is
#' the constant a caller passed as `k` when they passed one, times the scale the
#' fit settled on. That scale enters the location score as a known constant. It
#' solves an equation of its own that the stacked system does not write, so none
#' of its sampling variability is propagated, and the standard errors are those
#' of a Huber score read at a scale treated as fixed. The spread of the
#' conditional density is a separate quantity, and is the pooled residual
#' spread there as it is for every other class.
#'
#' Weights built with `.sigma = fit$s`, the robust scale, are read as the fixed
#' spread they are, but that combination is usually a poor choice: the robust
#' scale resists the very residuals the density is meant to spread over, so on
#' contaminated data the conditional density is far narrower than the residuals
#' are, and the ratio can collapse onto an effective sample size of about one.
#' Leave `.sigma` unset unless you have a reason to hold the spread at a
#' particular value.
#'
#' Three robust fits are refused. One fit with a psi that is not Huber, and one
#' fit with `method = "MM"`, which starts from a high-breakdown fit and finishes
#' on a redescending psi, are the roots of equations the stacked system does not
#' write; both error with class `propensity_ipw_robust_psi_error`. One whose
#' iteration stopped short of its own tolerance is the root of nothing, and
#' errors with class `propensity_ipw_convergence_error`. The first two point to
#' a Huber refit or to a bootstrap of the whole fit written by hand; the last
#' points to a larger `maxit` or a looser `acc`, which is what a fit needs to
#' reach its own tolerance.
#'
#' Every density [wt_ate()] builds a ratio in is stacked as it was built: the
#' normal, Laplace, and Student t families, and a density you wrote yourself,
#' under either the marginal or the integrated numerator. `ipw()` reads the
#' family, the numerator, and the spread off the record the weights carry and
#' rebuilds the same ratio at each value of the parameter vector, so the weights
#' it differentiates are the weights the outcome model was fit with. A marginal
#' numerator contributes the exposure's two marginal moments to the parameter
#' vector; an integrated one is built from the propensity score block and the
#' data alone and contributes none. Weights carrying no such record, including
#' any written with [psw()] by hand, are read as the normal ratio the package
#' built before the record existed.
#'
#' A numerator estimated by a model the caller passed to [wt_ate()]'s
#' `stabilize` contributes that model: one parameter per coefficient, named with
#' a `stab_` prefix so the terms it shares with the propensity score model stay
#' apart from it, and one for the spread its density is read at, `sigma2_n`. Its
#' score and the moment its spread is the root of are stacked alongside the
#' propensity score model's, so the standard errors account for the numerator
#' having been fitted, which is what separates this from a
#' `stabilization_score`: a score the caller computed is carried as a known
#' constant and contributes no parameter at all. The numerator model is read
#' through the registry above, so a class, family, or link whose score this
#' system cannot write is refused there in the terms that registry refuses a
#' propensity score model in, naming `stabilize`. A model of a response other
#' than the exposure, or one fit to a different set of observations, errors with
#' class `propensity_ipw_numerator_error`.
#'
#' Two fits are refused for the standard error rather than for the model. An
#' [mgcv::gam()] chooses how much to smooth by REML, and a `"kernel"` density
#' chooses its bandwidth from the residuals of the propensity score model;
#' neither is a function of the parameters that the sandwich could
#' differentiate, so weights built from either have no closed-form standard
#' error. Both raise an error of class
#' `propensity_ipw_se_method_unavailable_error`. The package offers no
#' resampling method of its own, so the remedy the refusal names is a bootstrap
#' of the whole pipeline written by hand: resample the rows, refit the
#' propensity score model, rebuild the weights with [wt_ate()], and refit the
#' outcome model on each resample, then read the spread of the coefficient
#' across the resamples. Doing it that way rather than inside `ipw()` keeps what
#' each replicate re-estimates, such as a kernel bandwidth or a smoothing
#' parameter, under the writer's control.
#'
#' ## Standard errors by linearization
#'
#' Setting `se_method = "linearization"` instead uses the influence-function
#' linearization of Kostouraki et al. (2024). It is available only for a binary
#' exposure and only for the `ate`, `att`, `ato`, and `atm` estimands; a
#' categorical or continuous exposure, or the `atu` or `entropy` estimand,
#' requires `se_method = "mestimation"`. The linearization also supports only an
#' outcome model whose formula contains the exposure alone. Its point estimates
#' come from g-computation on the fitted outcome model, while its influence
#' functions are those of the Hajek weighted-mean estimator; the two agree only
#' for an exposure-only outcome model. A covariate-adjusted outcome model (any
#' term beyond the exposure, including covariates, interactions, or transformed
#' terms) errors on the linearization path and requires
#' `se_method = "mestimation"`, which stacks adjusted outcome models correctly.
#' The conditional reading is restricted the same way: linearization stacks no
#' system, so the outcome model it stores carries no block of one, and the
#' standard errors that reading reports are refused with a classed error rather
#' than replaced by the covariance the outcome model computed for itself, which
#' treats the estimated weights as fixed.
#'
#' The linearization outcome model must also carry an intercept, which is a
#' stricter requirement than the M-estimation path imposes. Under a numeric
#' coding such as `y ~ z - 1` the linear predictor under no exposure is fixed at
#' the link's zero point rather than estimated, so the g-computation marginal
#' means no longer match the Hajek means the influence functions describe. The
#' rejection is broader than that mechanism: a saturated factor coding such as
#' `y ~ 0 + zf` does estimate both cell means, and the M-estimation path accepts
#' it, but the linearization influence functions are derived for the intercept
#' parameterization, so every no-intercept outcome model errors here. See
#' **Model requirements** for the baseline contract both methods impose.
#'
#' ## Standard errors as a diagnostic
#'
#' Setting `se_method = "robust"` reports the sandwich the weighted outcome
#' model computes for itself, with the weights entering as known constants. It
#' is the linearization route with the correction for having estimated the
#' propensity score dropped, so the point estimates are the linearization ones
#' and only the standard errors move. What it reports for the two counterfactual
#' means and the contrasts built from them is the delta-method reading of
#' `sandwich::vcovHC(outcome_mod, type = "HC0")`, exactly: the meat is divided
#' by `n` rather than by `n - 1`, so what comes back is the HC0 sandwich itself
#' rather than a finite-sample rescaling of it.
#'
#' This is offered as a diagnostic and not as an inference. It ignores that the
#' weights were estimated, and the correction it drops is generally what a
#' propensity score model contributes to the variance, so the standard errors it
#' reports generally understate it, and understate it most where the weights are
#' least stable. Use `"mestimation"` or `"linearization"` to report a standard
#' error. Reach for the diagnostic to see what the propensity correction is
#' worth on a given fit, by reading it beside one of them.
#'
#' A result fit this way carries the class `"ipw_diagnostic_se"` ahead of
#' `"ipw"`, and printing it writes one line naming the method, so the numbers
#' are not read as accounting for the propensity score model. `as.data.frame()`
#' and [`tidy()`][tidy.ipw()] carry the same mark as an `ipw_se_diagnostic`
#' attribute on the table. Every requirement of the linearization route applies
#' here: a binary exposure, the `ate`, `att`, `ato`, and `atm` estimands, an
#' outcome model of the exposure alone fit with an intercept, no `.by`, and no
#' conditional reading. The categorical, continuous, and joint treatment routes
#' refuse it as they refuse linearization.
#'
#' The printed mark does not survive [pool_ipw()]: a pooled object is a summary
#' of results rather than a result, and what says its standard errors are the
#' diagnostic ones is the method it records, `se_method = "robust"`.
#'
#' # Multiple imputation
#'
#' With missing data, fit the whole analysis once per imputed dataset and pool
#' the results. Everything the analysis needs is rebuilt inside a single
#' expression, so each imputation gets its own propensity model, its own weights,
#' and its own outcome model:
#'
#' ```r
#' imp <- mice::mice(dat, m = 20, print = FALSE)
#' fits <- with(imp, {
#'   ps <- glm(z ~ x1 + x2, family = binomial())
#'   w <- wt_ate(ps)
#'   om <- glm(y ~ z, family = quasibinomial(), weights = w)
#'   ipw(ps, om)
#' })
#' pool_ipw(fits)
#' ```
#'
#' [pool_ipw()] is the recommended verb. It reads the results themselves rather
#' than a tidied table of them, so the effect labels survive, a categorical
#' result keeps the contrast each row reports, and the complete-data degrees of
#' freedom fall back to the outcome models whenever a result records none of its
#' own. That fallback is written for the condition rather than for one route: a
#' fit under `se_method = "linearization"` or `se_method = "robust"` records
#' none, and so does a result another package built on the same class from
#' estimating equations that report no residual degrees of freedom. It also
#' takes the smallest degrees of freedom across the pooled results, which is the
#' conservative choice when they differ.
#' [`tidy()`][tidy.ipw_pooled()] and [`glance()`][glance.ipw_pooled()] report
#' what it returns.
#'
#' `mice::pool()` also works, for every exposure type, and is the right choice
#' when a result has to travel through the same pipeline as other analyses. Four
#' things are worth knowing about that route. It reaches the result through
#' [`tidy()`][tidy.ipw()], which pools a categorical result grouped by `term` and
#' `contrast` together, a grouping mice has supported since 3.15.0.
#' `mice::pool()` groups correctly but its `summary()` prints only `term`, so a
#' categorical pooled table shows each effect measure repeated with no contrast
#' label beside it, and the labels have to be read off `pooled$pooled` instead.
#' A result that records no complete-data degrees of freedom of its own, a fit
#' under `se_method = "linearization"` or `se_method = "robust"` among them,
#' leaves the pooled degrees of freedom missing unless `dfcom` is passed
#' explicitly, as in
#' `mice::pool(fits, dfcom = df.residual(fits$analyses[[1]]$outcome_mod))`. Every
#' result of that kind has to be told, where [pool_ipw()] reads the outcome
#' models for all of them unasked. That remedy is what the package's requirement
#' of mice 3.18.0 or later is for: mice 3.17.0 introduced a regression in the
#' `dfcom` argument of `pool()`, and 3.18.0 is the version that repairs it. And
#' the `exponentiate` argument of
#' `summary()` on a pooled `mipo` is not the one these methods take: a `mipo`
#' records no scale for the rows it holds, so it exponentiates every one of them,
#' returning a risk difference as its exponential and still labeling the row
#' `rd`. [`tidy()`][tidy.ipw_pooled()] and the pooled result's own frame
#' exponentiate the rows reported on the log scale alone and relabel those.
#'
#' To pool the conditional reading, record it inside the same expression with
#' [as_conditional()], and the pooled result reports the outcome models'
#' coefficients rather than the causal contrasts.
#'
#' Recording the reading beforehand is one route to it rather than the only one.
#' [pool_ipw()] pools both readings of the results it is given, storing the one
#' they record and keeping the other beside it, so [as_marginal()] and
#' [as_conditional()] move a pooled result between them after the pooling and
#' report a table the pooling already built rather than pooling a second time.
#' [`tidy()`][tidy.ipw_pooled()], [stats::coef()], [stats::vcov()],
#' [stats::confint()], and the pooled result's own frame each take an `effects`
#' argument naming a reading for the one call, which leaves the pooled result as
#' it is. `mice::pool()` offers neither: it reaches a result through
#' [`tidy()`][tidy.ipw()] and returns a `mipo` holding the one table it asked
#' for, so on that route the reading still has to be recorded before the pooling.
#'
#' A reading the pooling could not compute is refused rather than reported under
#' the other one's name. The conditional reading needs the covariance the joint
#' estimation of the weights and the outcome implies, and a fit under
#' `se_method = "linearization"` or `se_method = "robust"` stacks no such
#' system and records no such block, so a set of those fits pools the marginal
#' reading alone. The pooled result
#' records why the other reading is missing, and asking for it, whether by
#' flipping the result or by naming it for one call, errors with that recorded
#' explanation under the classes
#' `causalgenerics_pool_missing_surface_conditional` and
#' `causalgenerics_pool_missing_surface`.
#'
#' The analysis belongs inside `with()` rather than outside it because the
#' propensity score model has to be estimated once per imputation. Weights built
#' from propensity scores averaged across imputations are identical in every
#' imputation, which leaves no between-imputation variance for the weights to
#' contribute and no per-imputation estimation for the corrected standard errors
#' to account for; both components of the variance are lost. Fitting within each
#' imputation keeps the uncertainty of having estimated the weights inside each
#' result, and the pooling adds the uncertainty the imputation itself
#' contributed. Leyrat et al. (2019) compare the approaches and recommend this
#' one.
#'
#' # Model requirements
#'
#' `ipw()` cannot yet account for propensity scores that were trimmed
#' ([ps_trim()]), truncated ([ps_trunc()]), or calibrated ([ps_calibrate()])
#' before weighting. The stacked estimating equations rebuild the weights from
#' the propensity score model on every evaluation, so a weight that is no longer
#' a deterministic function of that model breaks the sandwich variance. Supplying
#' weights built from a modified score errors on either standard error method;
#' refit the weights from the unmodified propensity score model. An outcome model
#' fit without weights also errors on either method. The outcome model must not
#' carry an offset term on either method, since neither the stacked outcome score
#' nor the linearization influence functions thread an offset; supplying one
#' errors.
#'
#' For a binary or categorical exposure, the outcome model formula must contain
#' the exposure. The counterfactual designs are built by setting the exposure to
#' each level in turn, so a model fit on covariates alone gives one identical
#' design per level; it errors on either standard error method.
#'
#' For those two exposure types, the outcome model must also be able to represent
#' a baseline at every level, for example through an intercept or through a
#' saturated factor coding of the exposure. A numeric no-intercept coding such as
#' `y ~ z - 1` errors on both standard error methods, for different reasons. On
#' the M-estimation path the counterfactual design at the zero-coded level is
#' identically zero, so the marginal mean there is fixed by the outcome link
#' (`0.5` under a logit or probit link, `0` for a linear model) rather than
#' estimated from the data. On the linearization and robust paths the intercept
#' is required outright, because without it the g-computation means stop
#' matching the Hajek means the influence functions describe. A saturated factor
#' coding such as `y ~ 0 + zf` is a reparameterization whose designs are the
#' level indicators, so the M-estimation path accepts it and reproduces the
#' with-intercept fit; the linearization and robust paths still require the
#' intercept. A no-intercept model that keeps a covariate, such as
#' `y ~ z + x1 - 1`, still estimates the marginal mean from that covariate and
#' runs on the M-estimation path. A continuous exposure
#' has no counterfactual designs, since `ipw()` reports the marginal structural
#' model's own exposure coefficient, and is unaffected.
#'
#' Three further requirements apply to the weights and the propensity score
#' model. The outcome model weights must match the values implied by the
#' propensity score model; a mismatch errors, on both standard error methods.
#'
#' For a continuous exposure that requirement bears on the spread of the
#' conditional density. `ipw()` stacks a single pooled residual variance
#' alongside the propensity score coefficients, so the weights it rebuilds are
#' the ones [wt_ate()] produces with its pooled default. A single `.sigma` is
#' taken instead as a known constant: the weights are rebuilt at the number that
#' was supplied, and the stacked system carries none of that number's
#' uncertainty, which is what fixing a spread says. Weights built with an
#' observation-level `.sigma`, such as `influence(model)$sigma`, are a different
#' function of the data with no counterpart in the stacked system, and are
#' refused before anything is solved, with an error of class
#' `propensity_ipw_sigma_error`. Rebuild the weights with the pooled default or
#' with one number to use `ipw()`.
#'
#' The propensity score model must not separate the exposure. A model whose
#' covariates predict the exposure without error has no finite maximum likelihood
#' estimate, and the propensity scores its coefficients imply reach exactly zero
#' or one, leaving the corresponding weights undefined. For a binary exposure,
#' whose propensity score model is a binomial `glm()`, both standard error
#' methods reject such a fit at the same threshold: the fitted linear predictors,
#' put through the link's inverse, saturate for at least one observation. Nothing
#' short of saturation is rejected. A fitted model's `predict()` cannot show this
#' on its own, because the inverse link `glm` uses is bounded away from zero and
#' one, so a separated fit otherwise yields finite weights and a standard error
#' that looks ordinary. For a categorical exposure, which runs on the
#' M-estimation path alone, the threshold is narrower: only a probability of
#' exactly zero at the level a unit was actually assigned is rejected, since that
#' alone leaves the unit's own weight undefined, and a softmax column for a level
#' the unit was not assigned may reach zero without the fit being refused. A
#' continuous exposure has no saturating inverse link and is not checked.
#'
#' The weight functions apply a stricter rule to the propensity scores they are
#' handed: a categorical matrix holding an exact 0 or 1 anywhere is refused,
#' whichever level the cell belongs to. A separated [nnet::multinom()] fit is
#' therefore stopped where its weights are built rather than here. [ps_trunc()]
#' repairs such a matrix, its categorical matrix method reading the closed
#' interval so that a cell at an endpoint is bounded rather than refused, while
#' [ps_trim()] holds the open interval as the weight functions do. Weights built
#' from a truncated score are their own refusal here, since [ipw()] cannot yet
#' account for a modified propensity score, so a separated fit is repaired for
#' the weight functions rather than for this one. See **Propensity scores at 0
#' and 1** in [wt_ate()] for the remedy.
#'
#' The propensity score model must also be fit without case weights, since the
#' stacked propensity score equations are unweighted and a weighted fit would not
#' sit at the score root; that requirement applies to `se_method = "mestimation"`
#' alone, because neither the linearization path nor the robust one restacks the
#' propensity score model. Linearization still corrects for the uncertainty of
#' estimating the propensity scores; it does so through the influence functions
#' rather than through a stacked score. The robust diagnostic drops that
#' correction, which is what makes it a diagnostic.
#'
#' @references
#' Stefanski LA, Boos DD. The calculus of M-estimation. *The American
#' Statistician*. 2002;56(1):29--38. \doi{10.1198/000313002753631330}
#'
#' Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
#' inverse probability-of-treatment weighting estimator: A tutorial for
#' different types of propensity score weights. *Statistics in Medicine*.
#' 2024;43(13):2672--2694. \doi{10.1002/sim.10078}
#'
#' Leyrat C, Seaman SR, White IR, et al. Propensity score analysis with partially
#' observed covariates: How should multiple imputation be used? *Statistical
#' Methods in Medical Research*. 2019;28(1):3--19.
#' \doi{10.1177/0962280217713032}
#'
#' @return Methods of `ipw()` return an S3 object of class `ipw`. That result
#'   class is shared across packages and its components, its `print()` method,
#'   and its `as.data.frame()` method are documented at
#'   [causalgenerics::new_ipw()]. propensity's methods fill every component;
#'   three of them take propensity-specific values:
#' \describe{
#'   \item{`estimand`}{One of `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or
#'     `"entropy"`.}
#'   \item{`se_method`}{One of `"mestimation"`, `"linearization"`, or
#'     `"robust"`.}
#'   \item{`fit`}{The fitted M-estimator object when `se_method` is
#'     `"mestimation"`. Calling [stats::vcov()] or
#'     [generics::tidy()] on this object exposes the full stacked system of
#'     estimating equations, including the propensity score, outcome, and
#'     estimand parameters. It is `NULL` under `"linearization"` and
#'     `"robust"`, neither of which solves a system.}
#' }
#'
#'   The result answers the standard model accessors, which causalgenerics
#'   registers for the shared class: [stats::coef()] and [stats::confint()] for
#'   the reported effects, [stats::vcov()] for their covariance,
#'   [stats::nobs()] and [stats::df.residual()] for the counts describing the
#'   fit, and [stats::weights()] for the [psw()] vector the outcome model was
#'   fit with. Coefficients are named for the effect measure and the contrast
#'   together, the contrast being the level a counterfactual mean belongs to or
#'   the pair of levels an effect measure compares, and for the subgroup as well
#'   where `.by` reports one row per subgroup. Which surface [stats::coef()],
#'   [stats::vcov()], and [stats::confint()] report follows the presentation
#'   mode the result records, described under `effects` above.
#'
#'   Under `se_method = "mestimation"` the stored `wt_mod` and `outcome_mod` are
#'   the models supplied, carrying their block of the joint sandwich. Calling
#'   [stats::vcov()] on either reports a covariance that accounts for the whole
#'   system having been estimated from the same data, rather than the one the
#'   model's own fitting routine reports. The models supplied are left as they
#'   were fit. Linearization solves no such system, so its component models are
#'   stored unchanged.
#'
#' @examples
#' # Simulate data with a confounder, binary exposure, and binary outcome
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.5 * x1))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
#' dat <- data.frame(x1, z, y)
#'
#' # Step 1: Fit a propensity score model
#' ps_mod <- glm(z ~ x1, data = dat, family = binomial())
#'
#' # Step 2: Calculate ATE weights and fit a weighted outcome model
#' wts <- wt_ate(ps_mod)
#' outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
#'
#' # Step 3: Estimate causal effects with correct standard errors
#' result <- ipw(ps_mod, outcome_mod)
#' result
#'
#' # Exponentiate log-RR and log-OR to get RR and OR
#' as.data.frame(result, exponentiate = TRUE)
#'
#' # Continuous outcome example
#' y_cont <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
#' dat$y_cont <- y_cont
#' outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
#' ipw(ps_mod, outcome_cont)
#'
#' # Effect modification: the effect within each level of a factor, and the
#' # difference between them. The outcome model carries the exposure-by-modifier
#' # term that lets the effect differ across the strata.
#' dat$grp <- factor(rep(c("a", "b"), length.out = n))
#' ps_grp <- glm(z ~ x1 + grp, data = dat, family = binomial())
#' wts_grp <- wt_ate(ps_grp)
#' outcome_grp <- glm(
#'   y ~ z * grp + x1,
#'   data = dat,
#'   family = quasibinomial(),
#'   weights = wts_grp
#' )
#' ipw(ps_grp, outcome_grp, .by = grp)
#'
#' # The same fit with the weights stabilized on the modifier. The numerator is
#' # the probability of the exposure each unit took given `grp`, which is the
#' # shape the default stabilizer has and conditions only on a variable the
#' # outcome model reads, so it tightens the weights without moving the estimand.
#' num_grp <- glm(z ~ grp, data = dat, family = binomial())
#' p_grp <- fitted(num_grp)
#' wts_stab <- wt_ate(
#'   ps_grp,
#'   stabilize = TRUE,
#'   stabilization_score = ifelse(dat$z == 1, p_grp, 1 - p_grp)
#' )
#' sd(wts_grp)
#' sd(wts_stab)
#' outcome_stab <- glm(
#'   y ~ z * grp + x1,
#'   data = dat,
#'   family = quasibinomial(),
#'   weights = wts_stab
#' )
#' ipw(ps_grp, outcome_stab, .by = grp)
#'
#' # Continuous exposure: an lm propensity model of the dose on covariates,
#' # weights, which are stabilized by default for a continuous exposure, and a
#' # weighted marginal structural outcome model
#' a <- 0.5 + 0.8 * x1 + rnorm(n)
#' y_dose <- 1 + 0.6 * a + 0.3 * x1 + rnorm(n)
#' dat$a <- a
#' dat$y_dose <- y_dose
#' ps_cont <- lm(a ~ x1, data = dat)
#' wts_cont <- wt_ate(ps_cont)
#' msm <- lm(y_dose ~ a, data = dat, weights = wts_cont)
#' ipw(ps_cont, msm)
#'
#' # A dose-response curve written as a basis matrix reports one row per basis
#' # coefficient. The model frame of such a fit holds the basis rather than the
#' # dose, so `.data` supplies the exposure.
#' msm_curve <- lm(y_dose ~ poly(a, 2), data = dat, weights = wts_cont)
#' ipw(ps_cont, msm_curve, .data = dat)
#'
#' @examplesIf requireNamespace("nnet", quietly = TRUE)
#' # Categorical exposure: a multinomial propensity model and per-level contrasts
#' set.seed(2)
#' n <- 300
#' x1 <- rnorm(n)
#' a <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.3 * x1))
#' dat_cat <- data.frame(x1, a, y)
#'
#' ps_cat <- nnet::multinom(a ~ x1, data = dat_cat, trace = FALSE)
#' ps_mat <- predict(ps_cat, type = "probs")
#' wts_cat <- wt_ate(ps_mat, dat_cat$a, exposure_type = "categorical")
#' outcome_cat <- glm(y ~ a, data = dat_cat, family = quasibinomial(), weights = wts_cat)
#' ipw(ps_cat, outcome_cat)
#'
#' @seealso
#' [wt_ate()], [wt_att()], [wt_atu()], [wt_atm()], [wt_ato()], [wt_entropy()]
#'   for calculating propensity score weights.
#'
#' [ps_trim()], [ps_trunc()] for handling extreme propensity scores before
#'   weighting.
#'
#' @name ipw-methods
#' @exportS3Method causalgenerics::ipw glm
#' @importFrom causalgenerics new_ipw new_ipw_model
#' @importFrom stats cov dnorm family formula model.frame model.matrix
#' @importFrom stats model.weights pnorm predict qnorm var
ipw.glm <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  .by <- rlang::enquo(.by)
  se_method <- rlang::arg_match(se_method)

  # Whether the caller named a reading, read before `arg_match()` resolves the
  # default and loses the distinction. The continuous route needs it: a fit
  # whose exposure enters the outcome model through several columns records the
  # conditional reading either way, and what changes is whether that is a
  # default being explained or a reading the caller asked for.
  #
  # A whole set of readings is one of them no more than a missing argument is.
  # It is what a wrapper carrying this argument's default forwards, and a caller
  # of that wrapper who named nothing has named nothing here either.
  effects_named <- !missing(effects) && length(effects) == 1L
  effects <- rlang::arg_match(effects)
  assert_class(wt_mod, "glm")
  assert_class(outcome_mod, c("glm", "lm"))
  check_ipw_ps_response(wt_mod)

  # The two things about a `.by` request that can be settled before the modifier
  # is read, in the order a request fails in: the standard error method, then
  # the exposure type. The exposure refusal sits inside the continuous branch
  # below, so a request naming both is refused for the method it cannot use
  # rather than for the exposure it cannot use either.
  check_ipw_by_method(.by, se_method)

  # A gaussian-family propensity model indicates a continuous exposure. Route it
  # to the shared continuous path, which an lm propensity model also uses; the two
  # share fitted values and so produce identical estimates and standard errors.
  if (identical(wt_mod$family$family, "gaussian")) {
    check_ipw_by_exposure(.by)

    return(ipw_continuous_estimate(
      wt_mod,
      outcome_mod,
      .data = .data,
      estimand = estimand,
      ps_link = ps_link,
      conf_level = conf_level,
      se_method = se_method,
      effects = if (effects_named) effects
    ))
  }

  # Both standard error methods accept only the link `wt_mod` was fit with, so
  # the one value `ps_link` can carry without changing the result is the value
  # the default already resolves to. That leaves the argument as pure
  # redundancy. The warning comes after the continuous route above, so a
  # propensity model with no link to override is rejected outright rather than
  # deprecated and then rejected.
  if (!is.null(ps_link)) {
    lifecycle::deprecate_warn(
      "0.1.0",
      "ipw(ps_link)",
      details = "The link is always read from `wt_mod`, so `ps_link` can be \\
      dropped."
    )
  }

  # Guards first, on the weights that fit the outcome model. These carry the
  # psw attributes, so a modified propensity score is detected here before any
  # length reconciliation or estimand parsing that would otherwise fail
  # obliquely. Both standard error methods share these guards.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts)
  check_ipw_binary_focal(wts, wt_mod, .data = .data, estimand = estimand)

  if (identical(se_method, "mestimation")) {
    spec <- ipw_spec_binary(
      wt_mod,
      outcome_mod,
      .data = .data,
      estimand = estimand,
      ps_link = ps_link,
      .by = .by
    )
    fit <- ipw_mestimation(spec, conf_level = conf_level)
    components <- ipw_component_models(wt_mod, outcome_mod, fit)
    return(new_ipw(
      estimand = spec$estimand,
      wt_mod = components$wt_mod,
      outcome_mod = components$outcome_mod,
      estimates = fit$estimates,
      se_method = "mestimation",
      fit = fit$fit,
      effects = effects
    ))
  }

  exposure_name <- fmla_extract_left_chr(wt_mod)

  # Guards on the model structure run before the data is reconstructed, so a
  # model this path cannot support is rejected on its own terms rather than
  # through whatever the reconstruction happens to fail on first.
  check_ipw_offset(outcome_mod)
  check_ipw_outcome_response(outcome_mod)
  check_ipw_linearization_outcome(outcome_mod, exposure_name, se_method)
  check_ipw_outcome_family(outcome_mod)

  # Shared with the mestimation specs, so a propensity model that has lost the
  # data behind its fitting call raises the same guided error here and rebuilds
  # its design from `.data` the same way. The shared helper also reconciles a
  # `.data` whose row count disagrees with the fits, which this path needs
  # because everything below is sized to `.data` while the weights come from the
  # outcome model, and an unreconciled mismatch recycles the two against each
  # other and shrinks the standard errors with nothing signaled.
  extracted <- ipw_extract_ps_design(
    wt_mod,
    outcome_mod,
    .data = .data,
    exposure_name = exposure_name,
    check_exposure_levels = TRUE
  )
  exposure <- extracted$exposure
  outcome <- extracted$outcome
  weight_matrix <- extracted$ps_X

  # Convert a factor or logical response to its 0/1 coding so the influence
  # values compute Y - mu on numeric values rather than factor level codes.
  # Covers both extraction routes above.
  outcome <- ipw_outcome_numeric(outcome)

  # `predict.lm` applies the fit's own contrasts to `newdata`, so the scores are
  # right either way, but it re-levels against `xlevels` first and warns that
  # doing so dropped a contrasts attribute the column happened to carry. That
  # attribute is redundant to the one predict is about to apply, so clear it.
  #
  # The linear predictors are taken and the scores derived from them, rather than
  # both being predicted, because the saturation guard below needs the same
  # linear predictors and predicting twice invites the two to be computed from
  # different data. This is what `predict(type = "response")` does: it puts the
  # linear predictors through the family's own inverse link, and the two forms
  # agree to the bit.
  eta <- predict(
    wt_mod,
    newdata = drop_contrasts_attrs(.data, names(wt_mod$contrasts))
  )
  ps <- wt_mod$family$linkinv(eta)

  # The score factor, the weight derivatives, and the correction matrix are all
  # derived from `ps_link` and none of them consults the fitted model, so naming
  # a link other than the one `wt_mod` was fit with scales the estimation
  # correction by the wrong derivative. Only the standard errors move: the
  # estimates come from the outcome model and its weights and are untouched, and
  # the weight-consistency preflight cannot see it either, because the weights do
  # not depend on the link. That leaves nothing to notice, so reject it here.
  if (is.null(ps_link)) {
    ps_link <- wt_mod$family$link
  } else if (!identical(ps_link, wt_mod$family$link)) {
    fitted_link <- wt_mod$family$link
    abort(
      c(
        "{.arg ps_link} must match the link {.arg wt_mod} was fit with for \\
        {.val {se_method}} standard errors.",
        x = "{.arg ps_link} is {.val {ps_link}}; {.arg wt_mod} was fit with a \\
        {.val {fitted_link}} link.",
        i = "Omit {.arg ps_link} to use the model's own link, or refit \\
        {.arg wt_mod} with the link you intend."
      ),
      error_class = "propensity_ipw_link_error"
    )
  }

  # Runs before the weight-consistency preflight and before the estimand is
  # resolved, so a fit with no overlap at all is diagnosed as what it is rather
  # than as whatever the weights or the estimand happen to trip over first. The
  # estimand restriction below is the one worth naming: an estimand this path
  # has no influence functions for is refused with a pointer to the M-estimation
  # path, which on a separated fit would only meet the same separation there.
  check_ipw_ps_saturation(eta, link = wt_mod$family$link)

  if (!identical(length(exposure), length(outcome))) {
    abort(
      c(
        "{.arg exposure} and {.arg outcome} must be the same length.",
        x = "{.arg exposure} is length {length(exposure)}",
        x = "{.arg outcome} is length {length(outcome)}"
      ),
      error_class = "propensity_length_error"
    )
  }

  estimand <- check_estimand(wts, estimand)

  # The linearization influence functions are derived only for ate, att, ato, and
  # atm; atu and entropy require the mestimation path. Reject them here with the
  # documented classed error rather than letting the request reach derive_weights,
  # whose internal arg_match raises a bare, misleading error.
  #
  # The robust diagnostic uses none of those derivatives, but it is the same
  # route with the propensity correction dropped and it reports the same means,
  # so it accepts the estimands this route accepts and no others. That leaves
  # one set of requirements to learn rather than two.
  if (!estimand %in% c("ate", "att", "ato", "atm")) {
    abort(
      c(
        "{.fun ipw} does not support {.val {se_method}} standard errors for \\
        the {.val {estimand}} estimand.",
        i = "Use {.code se_method = \"mestimation\"} for the {.val {estimand}} \\
        estimand."
      ),
      error_class = "propensity_method_error"
    )
  }

  # Recode the exposure to 0/1 for the influence functions, in parity with the
  # M-estimation recode. The marginal means below keep the original exposure so
  # they set the counterfactual levels the outcome model expects; the influence
  # values need the numeric indicator (wts * Z, exposure == 1, Z - ps).
  exposure_binary <- ipw_recode_binary_exposure(exposure)

  # The propensity scores above and the weights that fit the outcome model are
  # taken from two different places and nothing so far requires them to describe
  # the same analysis. Recompute the weights the scores imply and reject a
  # disagreement, the same check and the same message the M-estimation path
  # applies. It matters here because a mismatch moves only the standard errors:
  # the estimates come from the outcome model and its weights, so the output
  # looks untouched beside a variance that is quietly wrong.
  #
  # The stabilizer mirrors the M-estimation init seed exactly: a per-observation
  # score is used when the weights carry one, and otherwise a stabilized weight
  # is rebuilt from the group-constant default both paths share,
  # `ipw_default_stab_seed()`. That equivalence is what lets weights whose score
  # was dropped in a length-changing operation be caught here rather than
  # silently rebuilt from the wrong stabilizer.
  stab_score <- stabilization_score(wts)
  recomputed_wts <- ipw_weight_fn("binary", estimand)(
    ps,
    exposure_binary,
    list(
      stab_prob = if (is_stabilized(wts) && is.null(stab_score)) {
        ipw_default_stab_seed(exposure_binary)
      },
      score = stab_score
    )
  )
  ipw_compare_weights(recomputed_wts, wts, "binary", estimand)

  marginal_means <- estimate_marginal_means(
    outcome_mod = outcome_mod,
    wts = wts,
    exposure = exposure,
    exposure_name = exposure_name,
    .data = .data
  )

  uncorrected_lin_vars <- linearize_variables_for_wts(
    exposure_binary,
    outcome,
    wts,
    marginal_means
  )

  # The robust diagnostic stops here. Its influence functions are those of the
  # Hajek means with the weights entering as known constants, which is what
  # `linearize_variables_for_wts()` returns; the correction below is the whole
  # of what accounts for having estimated the propensity score, and dropping it
  # is what makes the reported variance the outcome model's own sandwich.
  lin_vars <- if (identical(se_method, "robust")) {
    uncorrected_lin_vars
  } else {
    linearize_variables_for_ps(
      exposure = exposure_binary,
      outcome = outcome,
      wts = wts,
      ps = ps,
      estimand = estimand,
      weight_matrix = weight_matrix,
      marginal_means = marginal_means,
      uncorrected_lin_vars = uncorrected_lin_vars,
      ps_link = ps_link
    )
  }

  estimates <- calculate_estimates(
    lin_vars = lin_vars,
    marginal_means = marginal_means,
    n = length(outcome),
    linear_regression = is_linear_regression(outcome_mod),
    conf_level = conf_level,
    exposure_levels = marginal_means$levels,
    se_method = se_method
  )

  # The same collapse the M-estimation path reports, from the seam this path
  # finishes its estimates at. A degenerate fit is a property of the data and
  # the models rather than of how the variance was computed, and this path
  # produces the same shape: an outcome constant within each exposure arm
  # returns its contrast with an interval of no width.
  warn_ipw_degenerate_se(estimates)

  result <- new_ipw(
    estimand = estimand,
    wt_mod = wt_mod,
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = se_method,
    fit = NULL,
    effects = effects
  )

  ipw_mark_diagnostic_se(result)
}

# Shared continuous-exposure route for ipw.lm and the gaussian-family branch of
# ipw.glm. An lm and a gaussian glm propensity model share fitted values, so both
# reach the same M-estimation fit. Standard errors come only from M-estimation
# for a continuous exposure; the linearization path is not available.
#
# `effects` is `NULL` when the caller named no reading, which the methods above
# read off `missing()` before `arg_match()` resolves the default. The two are
# the same request for every fit but one: an exposure entering the outcome
# model through several columns records the conditional reading whichever way
# it was asked for, and only the unnamed request is announced.
ipw_continuous_estimate <- function(
  wt_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = "mestimation",
  effects = NULL,
  call = rlang::caller_env()
) {
  # Before the class guard below, which reports a matrix-response fit as the
  # `mlm` its class vector says it is and sends the user to refit with the
  # `lm()` they already used. The response guard names the response instead.
  check_ipw_ps_response(wt_mod, call = call)

  # The stacked system carries the propensity score model's own score and
  # rebuilds its design from the model formula, so the registry decides what can
  # be here: a class whose coefficients are not the root of any equation it
  # writes, such as a robust fit descending a redescending loss, would drift
  # silently to the root of the one it does write, and a class whose design is
  # not the formula's, such as an additive fit's smooth basis, would be
  # reconstructed as something else.
  # Either way the estimates would describe a propensity score model the user
  # never fit. The refusal runs before every other guard here so that a
  # recognized model with no sandwich is reported as the method it lacks rather
  # than as the class it is.
  check_ipw_continuous_model(
    ipw_continuous_model(wt_mod, call = call),
    call = call
  )

  check_ipw_ps_link_absent(ps_link, "continuous", call = call)

  if (!identical(se_method, "mestimation")) {
    abort(
      c(
        "{.fun ipw} does not support {.val {se_method}} standard errors for \\
        continuous exposures.",
        i = "Use {.code se_method = \"mestimation\"} for a continuous exposure."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  # Guard the weights that fit the outcome model before building the spec, so a
  # modified propensity score is detected before any estimand parsing.
  wts <- extract_weights(outcome_mod)
  check_ipw_weights(wts, call = call)

  spec <- ipw_spec_continuous(
    wt_mod,
    outcome_mod,
    .data = .data,
    estimand = estimand,
    call = call
  )

  # The readings the result declares and the one it records, settled before
  # anything is solved so that a request for a reading the fit does not have is
  # refused rather than answered after the work of building it.
  dose_basis <- ipw_is_dose_basis(spec)
  readings <- c("marginal", "conditional")
  if (dose_basis) {
    check_ipw_dose_basis_effects(effects, call = call)
    if (is.null(effects)) {
      ipw_dose_basis_announce()
    }
    effects <- "conditional"
    readings <- "conditional"
  } else if (is.null(effects)) {
    effects <- "marginal"
  }

  fit <- ipw_mestimation(spec, conf_level = conf_level, call = call)
  components <- ipw_component_models(wt_mod, outcome_mod, fit)

  # A dose-basis fit declares the one reading it has, so the result class
  # refuses the other wherever it is asked for and this package adds no methods
  # of its own to carry that refusal.
  new_ipw(
    estimand = spec$estimand,
    wt_mod = components$wt_mod,
    outcome_mod = components$outcome_mod,
    estimates = fit$estimates,
    se_method = "mestimation",
    fit = fit$fit,
    effects = effects,
    readings = readings
  )
}

# Guard the weights that fit the outcome model. ipw() cannot yet account for
# propensity scores that were trimmed, truncated, or calibrated before
# weighting, so detect them here and direct the user to refit from the
# unmodified model. An outcome model fit without weights cannot yield an IPW
# estimate at all.
check_ipw_weights <- function(wts, call = rlang::caller_env()) {
  if (is.null(wts)) {
    abort(
      c(
        "{.arg outcome_mod} must be fit with propensity score weights.",
        x = "No {.arg weights} were found in {.arg outcome_mod}.",
        i = "Fit {.arg outcome_mod} with weights from a propensity weight \\
        function such as {.fun wt_ate}."
      ),
      error_class = "propensity_ipw_weights_missing_error",
      call = call
    )
  }

  if (is_ps_trimmed(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with trimmed propensity score weights.",
        x = "{.fun ipw} cannot yet account for trimmed propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_trimmed_error",
      call = call
    )
  }

  if (is_ps_truncated(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with truncated propensity score weights.",
        x = "{.fun ipw} cannot yet account for truncated propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_truncated_error",
      call = call
    )
  }

  if (is_ps_calibrated(wts)) {
    abort(
      c(
        "{.arg outcome_mod} was fit with calibrated propensity score weights.",
        x = "{.fun ipw} cannot yet account for calibrated propensity scores.",
        i = "Refit the weights from the unmodified propensity score model; \\
        support for modified weights is planned."
      ),
      error_class = "propensity_ipw_calibrated_error",
      call = call
    )
  }

  invisible(wts)
}

# Reject a propensity fit whose scores have saturated, at the same threshold the
# M-estimation path applies. Nothing on the linearization path breaks on such a
# fit, which is exactly the problem: `predict(type = "response")` goes through
# the family's inverse link, a clamped C routine that never returns an exact 0 or
# 1, so a model with no finite maximum likelihood estimate hands back interior
# scores, finite weights, and a confident standard error beside an estimate that
# no data support. Recomputing the fitted linear predictors through the link's
# true inverse is what makes it visible, and it is the same quantity the
# M-estimation guard counts, so the two paths refuse the same fits.
#
# `eta` is the fitted linear predictors the analysis itself runs on, so the guard
# covers the model-frame and the `.data` routes alike and cannot measure a
# different fit than the one being analyzed. `link` is the link `wt_mod` was fit
# with rather than the caller's `ps_link`, because saturation is a property of
# the fit and not of an argument.
#
# A fit with a link outside `ipw_inv_links`, cauchit for instance, has no
# unclamped inverse to recompute through, so the guard steps aside rather than
# aborting: it must never be the thing that stops an analysis it cannot measure.
# An `NA` score is missing rather than saturated and is likewise not this guard's
# business; leaving it in the count would put an `NA` into the comparison below
# and stop the call with base R's report about a missing value instead. No route
# through `ipw()` produces one today: a `.data` carrying a missing covariate is
# refused before the scores are predicted, and the model-frame route predicts
# over the design the fit kept, which has the incomplete rows already dropped.
# The exclusion is held against a change to either.
check_ipw_ps_saturation <- function(
  eta,
  link,
  call = rlang::caller_env()
) {
  inv_link <- ipw_inv_links[[link]]
  if (is.null(inv_link)) {
    return(invisible(TRUE))
  }

  ps <- inv_link(as.vector(eta))
  n_saturated <- sum(ps == 0 | ps == 1, na.rm = TRUE)

  if (n_saturated > 0) {
    abort_ipw_ps_separation(
      n_saturated,
      lead = "Putting the fitted linear predictors through the link's inverse",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject att or atu weights built with a focal level other than the one the
# binary path treats as focal. That path has no focal level of its own: it codes
# the second sorted exposure level as 1, and both the estimand derivatives and
# the weight recomputation assume that coding. att and atu weights are mirror
# images of each other, so weights targeting the other level carry values that
# are correct for the roles reversed and would be corrected as though they were
# not, silently reporting a different estimand than the one requested.
#
# Only weights that record a focal level can be checked, and only against an
# exposure whose own values contain it: a level this exposure never takes says
# nothing about this analysis, since equivalent codings of the same exposure
# carry different labels. Those fall to the weight-consistency preflight.
check_ipw_binary_focal <- function(
  wts,
  wt_mod,
  .data = NULL,
  estimand = NULL,
  call = rlang::caller_env()
) {
  focal <- attr(wts, "focal_category")
  if (is.null(focal)) {
    return(invisible(TRUE))
  }

  # Only att and atu weights ever record a focal level, so the weights' own
  # estimand answers this whenever it is present.
  wt_estimand <- if (is_causal_wt(wts)) estimand(wts) else estimand
  if (!isTRUE(wt_estimand %in% c("att", "atu"))) {
    return(invisible(TRUE))
  }

  exposure_name <- fmla_extract_left_chr(wt_mod)
  exposure <- if (!is.null(.data) && exposure_name %in% names(.data)) {
    .data[[exposure_name]]
  } else {
    # A propensity model that has lost the data behind its fitting call raises
    # its own guided error further down; say nothing here.
    tryCatch(fmla_extract_left_vctr(wt_mod), error = function(e) NULL)
  }
  if (is.null(exposure)) {
    return(invisible(TRUE))
  }

  # The same ordering ipw_recode_binary_exposure() uses, so the coded-1 value is
  # exactly the level this path treats as focal. Levels are compared as
  # characters: a focal level of 1 recorded against a 0/1 numeric exposure is a
  # match, not a mismatch.
  exposure_values <- as.character(sort(unique(exposure)))
  if (length(exposure_values) != 2) {
    return(invisible(TRUE))
  }
  focal <- as.character(focal)
  coded_focal <- exposure_values[[2]]

  if (
    length(focal) != 1 ||
      !focal %in% exposure_values ||
      identical(focal, coded_focal)
  ) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "The weights that fit {.arg outcome_mod} target a different exposure \\
      level than {.fun ipw} does.",
      x = "The weights record {.val {focal}} as the focal level, but \\
      {.fun ipw} treats {.val {coded_focal}} as focal: it takes the second \\
      sorted level of a binary exposure as the exposed group.",
      i = "Relevel {.arg {exposure_name}} so that {.val {focal}} sorts \\
      second, then refit both {.arg wt_mod} and {.arg outcome_mod}."
    ),
    error_class = "propensity_focal_level_error",
    call = call
  )
}

# Reject a non-NULL `ps_link` for a propensity model whose link it cannot name.
# `ps_link` restates the link of a binomial glm on the binary path, and takes
# only that model's links; a multinomial propensity model has no link at all and
# a continuous one is read through a link of its own that this argument cannot
# name, so accepting it there would silently ignore it. `model_kind` names the
# model in the message and is the only thing that differs between the two
# callers.
check_ipw_ps_link_absent <- function(
  ps_link,
  model_kind,
  call = rlang::caller_env()
) {
  if (!is.null(ps_link)) {
    abort(
      c(
        "{.fun ipw} does not accept {.arg ps_link} for a {model_kind} \\
        propensity score model.",
        x = "{.arg ps_link} names the link of a binomial glm, which a \\
        {model_kind} propensity score model is not.",
        i = "Omit {.arg ps_link}; it applies only to a binomial glm propensity \\
        score model."
      ),
      error_class = "propensity_ipw_link_error",
      call = call
    )
  }

  invisible(TRUE)
}

# A `multinom` fit to two levels is a binary propensity score model written with
# a multinomial fitter: it reports one probability, the shape a binomial glm
# reports, and the weight functions read such a fit on the binary path. The
# categorical stacked system has no multinomial coefficient block to build from
# it, so the estimator refuses the pair of models here, before deli is asked for
# a multinomial score it cannot write and reports the shortfall in terms of its
# own arguments.
check_ipw_multinom_levels <- function(wt_mod, call = rlang::caller_env()) {
  if (!identical(length(wt_mod$lev), 2L)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} needs a propensity score model of three or more exposure \\
      levels for a categorical exposure.",
      x = "{.arg wt_mod} was fit to 2 levels ({.val {wt_mod$lev}}), so it \\
      fits the single probability of a binary exposure.",
      i = "For a binary exposure, fit the propensity score model with \\
      {.fun stats::glm} and {.code family = binomial()}; for a categorical \\
      one, fit a {.fun nnet::multinom} to its three or more levels."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# Require the propensity score model's response to be the exposure column
# itself. Everything downstream names the exposure by deparsing that left-hand
# side, so a left-hand side that is not a single symbol yields a vector of names
# and is then indexed as though it were one. Two shapes do this and the guard
# rejects both: a matrix response, `cbind(successes, failures)`, which is not one
# exposure at all, and a transformed or otherwise call-form response, `log(a)` or
# `factor(z)`, which is one column but not one the data holds under that name.
#
# The matrix case is detected two ways: the model frame's response is a matrix,
# and, as a frame-independent fallback for a model fit without a stored frame,
# the left-hand side deparses to more than one element. Which of the two shapes
# a rejected model has is read from the frame where there is one and from the
# left-hand side where there is not, since a model fit with `model = FALSE`
# whose data are gone can still be seen to write `cbind()`. Any other way of
# building a matrix response, a variable that already holds one for instance,
# reads as call-form once the frame is unavailable. The call-form case is what
# remains once neither says matrix.
#
# Shared by the ipw.glm entry, ipw_spec_binary, and the continuous route, whose
# lm and gaussian glm propensity models share fitted values and therefore have to
# reject the same shapes with the same message. Left unguarded, the continuous
# route died on an internal length assertion about `.exposure_name` without
# `.data` and asked for a column named after the transforming function with it.
#
# The error class says what is wrong with the model rather than which model it
# is: `propensity_ipw_response_error` belongs to this guard and to the outcome
# model's counterpart below, and to nothing else. It is deliberately not
# `propensity_ipw_exposure_error`, which is reserved for the guards about what
# the exposure means rather than what shape a response has: an outcome model
# that omits the exposure, one that codes it numerically, one whose factor
# levels are ordered against the propensity model's, the counterfactual design
# guards, and an exposure value the propensity model never saw.
#
# `arg` names the argument the fit arrived as, the way it does in
# `check_multinom_response()`. Every route into this guard takes the propensity
# score model as `wt_mod`, which is the default.
check_ipw_ps_response <- function(
  ps_mod,
  arg = "wt_mod",
  call = rlang::caller_env()
) {
  # A fit made through a matrix interface records no formula and no terms, so
  # there is no left-hand side to read the exposure off. Reading the formula
  # defensively keeps that a refusal of this package's rather than a subscript
  # error raised on a formula that was never written.
  fmla <- tryCatch(stats::formula(ps_mod), error = function(e) NULL)
  if (!rlang::is_formula(fmla) || length(fmla) < 3L) {
    abort(
      c(
        "{.fun ipw} needs a propensity score model fit through the formula \\
        interface.",
        x = "{.arg {arg}} records no formula with a response, so the exposure \\
        cannot be read off it.",
        i = "Refit {.arg {arg}} from a formula whose left-hand side is the \\
        exposure, as in {.code exposure ~ x}, rather than from a design \\
        matrix and a response vector."
      ),
      error_class = "propensity_ipw_response_error",
      call = call
    )
  }

  lhs <- fmla[[2]]
  lhs_not_single <- length(fmla_extract_left_chr(ps_mod)) != 1L
  frame_response_is_matrix <- tryCatch(
    is.matrix(stats::model.response(stats::model.frame(ps_mod))),
    error = function(e) FALSE
  )
  response_is_matrix <- frame_response_is_matrix ||
    (is.call(lhs) && identical(lhs[[1]], quote(cbind)))

  if (!lhs_not_single && !response_is_matrix) {
    return(invisible(TRUE))
  }

  # A matrix response and a call-form response are different mistakes, and the
  # cbind wording sent the writer of `factor(z) ~ x` looking for a matrix they
  # never wrote.
  problem <- if (response_is_matrix) {
    c(
      "{.fun ipw} does not support a matrix response in the propensity score \\
      model.",
      x = "{.arg {arg}} has a matrix response, such as \\
      {.code cbind(successes, failures)}; the exposure must be a \\
      single-column response."
    )
  } else {
    lhs_text <- deparse1(lhs)
    c(
      "{.fun ipw} does not support a transformed response in the propensity \\
      score model.",
      x = "{.arg {arg}} reads the exposure through {.code {lhs_text}}, an \\
      expression rather than a single column."
    )
  }

  abort(
    c(
      problem,
      i = "Fit {.arg {arg}} with the exposure itself as the response, adding \\
      it to the data as its own column first if it has to be computed."
    ),
    error_class = "propensity_ipw_response_error",
    call = call
  )
}

# Reject an outcome model with a matrix response, the counterpart of the guard
# above. A cbind(successes, failures) response leaves the model frame twice as
# long as the exposure, and what the user was told depended on the route:
# without `.data` the length reconciliation reported an outcome twice the length
# of the exposure, and with `.data` the extraction asked for a column named
# "cbind", which cannot exist. Neither named the response, and the second sent
# the user to add a column named after a function.
#
# The prongs are not the ps guard's. There, a left-hand side that is not a bare
# symbol is rejected, which is right for an exposure. Here it would be wrong: a
# transformed single-column response such as `log(y)`, `sqrt(y)`, or `I(y / 10)`
# is an ordinary outcome model that estimates marginal means correctly, and "not
# a single symbol" is true of every one of them. So the frame-independent prong
# is narrowed to a literal `cbind` call, which is a matrix response by
# construction and cannot match a transformation, leaving the general test to
# `is.matrix` on the built frame.
#
# That narrowed prong is redundant on every public route, since each `ipw()`
# method builds the outcome frame to read its weights before any spec is
# constructed. It is kept for callers that reach the spec constructors directly,
# where this guard runs before that frame is ever demanded.
check_ipw_outcome_response <- function(
  outcome_mod,
  call = rlang::caller_env()
) {
  lhs <- fmla_extract_left_chr(outcome_mod)
  lhs_is_cbind <- length(lhs) > 1L && identical(lhs[[1]], "cbind")
  response_is_matrix <- tryCatch(
    is.matrix(stats::model.response(stats::model.frame(outcome_mod))),
    error = function(e) FALSE
  )

  if (lhs_is_cbind || response_is_matrix) {
    abort(
      c(
        "{.fun ipw} does not support a matrix response in the outcome model.",
        x = "{.arg outcome_mod} has a matrix response, such as \\
        {.code cbind(successes, failures)}; the marginal means are estimated \\
        from a single-column response.",
        i = "Refit {.arg outcome_mod} on one row per observation, weighted by \\
        the propensity score weights. Aggregated counts have to be expanded \\
        first, because the outcome model must carry exactly those weights."
      ),
      error_class = "propensity_ipw_response_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject an outcome model that carries an offset. Neither standard error path
# accounts for one: the stacked estimating equations rebuild the weighted
# outcome score without threading an offset, and the linearization influence
# functions are derived for the same unoffset weighted fit, so the reported
# results would disagree with the model that was fit. Detect an offset supplied
# through the model formula (the terms offset attribute) or through the `offset`
# argument (the model frame) and direct the user to refit without it. Called
# from every mestimation spec path and from the linearization branch.
check_ipw_offset <- function(outcome_mod, call = rlang::caller_env()) {
  has_offset <- !is.null(attr(stats::terms(outcome_mod), "offset")) ||
    !is.null(stats::model.offset(stats::model.frame(outcome_mod)))

  if (has_offset) {
    abort(
      c(
        "{.arg outcome_mod} was fit with an offset term.",
        x = "{.fun ipw} does not support offset terms in the outcome model; \\
        the stacked estimating equations do not thread an offset.",
        i = "Refit {.arg outcome_mod} without the offset."
      ),
      error_class = "propensity_ipw_offset_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject an outcome model that never references the exposure. The counterfactual
# designs are built by setting the exposure in the outcome model's data to each
# level in turn, so a model whose formula omits the exposure yields identical
# designs: every contrast collapses to zero with a degenerate standard error
# rather than failing. Read the model terms rather than the model frame so the
# guard fires on both extraction routes. Supplying `.data` bypasses the
# frame-based check_exposure() entirely, and without `.data` check_exposure()
# would direct the user to supply it, which drops them into that same silent
# degenerate path. A transformed exposure such as `I(z^2)` still names the
# exposure among its variables, so this guard passes it through: with `.data` the
# design is rebuilt from the supplied column and the fit runs, and without it the
# term is unreadable from the model frame and check_exposure() reports that.
check_ipw_outcome_exposure <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  out_terms <- stats::delete.response(stats::terms(outcome_mod))

  if (!exposure_name %in% all.vars(out_terms)) {
    abort(
      c(
        "{.arg outcome_mod} must contain the exposure {.val {exposure_name}}.",
        x = "{.val {exposure_name}} does not appear in the formula of \\
        {.arg outcome_mod}.",
        i = "Refit {.arg outcome_mod} with {.val {exposure_name}} on the \\
        right-hand side of the formula."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject a categorical outcome model that holds the exposure as a numeric term.
# The counterfactual designs are built by replacing the exposure column with a
# factor of the fitted levels, so `model.matrix` expands it into one dummy per
# non-reference level. A numeric term occupies a single coefficient slot, so the
# design gains columns the coefficient vector does not have, and the marginal
# mean is `X %*% beta` with no names consulted: the multiply is positional, and
# here it is not even conformable.
#
# The two codings mean different things in any case. A numeric term fits one
# slope across the levels and assumes they are equally spaced, which is not the
# model the counterfactual designs represent, so reconciling the two would
# silently answer a different question than the user asked.
#
# Only "numeric" is rejected. `.MFclass` folds integer columns into "numeric",
# so that single value covers both storage types. A character column is left
# alone because the model frame keeps it as character and `model.matrix`
# expands it into the same dummy columns an explicit factor would give, for the
# same estimates. A logical column
# cannot reach here at all: it carries at most two levels, and the weight layer
# rejects a two-level categorical exposure before `ipw()` is ever called.
#
# The guard is silent when the exposure name is absent from `dataClasses`, which
# happens when the formula transforms it, as in `y ~ factor(a) + x1`. The term
# label is the call rather than the name there, and the levels and the coding
# belong to the label, which the counterfactual rebuild reproduces by level
# rather than by position. check_ipw_exposure_rebuild() decides that case: it is
# rebuilt from `.data` and rejected without it.
check_ipw_outcome_exposure_class <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  data_classes <- attr(stats::terms(outcome_mod), "dataClasses")

  if (!exposure_name %in% names(data_classes)) {
    return(invisible(TRUE))
  }

  if (identical(data_classes[[exposure_name]], "numeric")) {
    abort(
      c(
        "{.arg outcome_mod} must hold a categorical exposure as a factor.",
        x = "{.val {exposure_name}} enters {.arg outcome_mod} as a numeric \\
        term, which gives it one coefficient instead of one per level.",
        i = "Refit {.arg outcome_mod} after converting {.val {exposure_name}} \\
        to a factor in the data, rather than wrapping it in the formula."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject an outcome model whose exposure term cannot be rebuilt at a
# counterfactual value. The marginal means are estimated by setting the exposure
# column to one value at a time and rebuilding the outcome design, so a term that
# reads the exposure through a call has to be re-evaluated at that constant
# column and land on the coding the fit recorded.
#
# What that takes depends on the route, so `rebuilt` says which one is in play.
#
# With `.data` the design is rebuilt from the supplied frame, which evaluates
# every term afresh and re-levels the result against the fit. A term such as
# `factor(z)` has only the one level at a constant column, but the fitted levels
# are recorded under the term's own label and put back, so the design is the
# fitted design and the marginal means are the ones the pre-converted factor
# model gives. What is left to reject there is a term that cannot be evaluated
# at the value being set, or that evaluates to something outside the levels it
# was fit with, neither of which any re-leveling can repair.
#
# Without `.data` the design comes from the outcome model's own model frame,
# which carries its terms attribute. `model.matrix()` then reads the columns
# already in that frame and evaluates nothing, so the counterfactual column
# written beside them is dropped and every level's design is the fitted one:
# each contrast is zero with nothing signaled. Any term whose levels the fit had
# to work out from the values it saw is rejected there, and the remedy is to
# supply `.data`.
#
# The check evaluates rather than reads the spelling, because the spelling does
# not decide it. `cut(z, breaks)` takes its levels from the call and carries all
# of them at a constant column; `factor(z)` does not. A numeric transformation
# such as `I(z^2)` is not considered at all: it has no levels to lose.
#
# Only the first offending term is reported: the fix is a refit or a `.data`,
# after which the guard runs again over the rest.
check_ipw_exposure_rebuild <- function(
  outcome_mod,
  exposure_name,
  exposure,
  data,
  rebuilt,
  call = rlang::caller_env()
) {
  xlevels <- outcome_mod$xlevels
  data_classes <- attr(stats::terms(outcome_mod), "dataClasses")
  candidates <- names(data_classes)[
    data_classes %in% c("factor", "ordered", "character")
  ]
  candidates <- intersect(setdiff(candidates, exposure_name), names(xlevels))
  candidates <- candidates[
    vapply(candidates, ipw_term_reads, logical(1), exposure_name)
  ]

  if (length(candidates) == 0) {
    return(invisible(TRUE))
  }

  values <- sort(unique(exposure))
  columns <- as.list(data)
  columns[[exposure_name]] <- NULL
  # Terms are evaluated where the formula was written, so a call to a function
  # the user defined resolves the same way the fit resolved it.
  enclos <- environment(stats::formula(outcome_mod))
  if (!is.environment(enclos)) {
    enclos <- baseenv()
  }

  for (term in candidates) {
    expr <- str2lang(term)
    fit_levels <- xlevels[[term]]

    for (value in values) {
      probe <- tryCatch(
        eval(
          expr,
          c(
            columns,
            stats::setNames(list(rep(value, length(exposure))), exposure_name)
          ),
          enclos
        ),
        error = function(e) e
      )
      rebuilds <- if (inherits(probe, "error")) {
        FALSE
      } else if (rebuilt) {
        # Carrying some of the fitted levels is enough here, and a constant
        # column carries exactly one: the rebuild re-levels the term's value
        # against the whole set, so the design it is read on is the fitted one.
        all(levels(as.factor(probe)) %in% fit_levels)
      } else {
        identical(levels(as.factor(probe)), fit_levels)
      }

      if (!rebuilds) {
        abort_ipw_exposure_rebuild(
          exposure_name,
          term,
          fit_levels,
          rebuilt = rebuilt,
          call = call
        )
      }
    }
  }

  invisible(TRUE)
}

# Reject a model variable recorded under a call that reads the exposure, on the
# route that cannot recompute it. The check above settles every variable that
# carries levels, by evaluating it and comparing what comes back. One with no
# levels has nothing to compare, and without `.data` it does not need comparing:
# the design comes from the outcome model's own model frame, `model.matrix()`
# reads the column already sitting there, and the counterfactual value is
# written into a column the term never consults again. `y ~ z + x1 + I(z * x1)`
# fit on a numerically coded exposure then contributes its fitted interaction to
# both counterfactual designs, and the marginal means come back as if the
# interaction were not in the model. It is a movement in the estimates rather
# than an error, and only on this route: with `.data` the term is recomputed at
# each value and gives what `y ~ z * x1` gives.
#
# What decides this is the model frame's variables rather than the formula's
# terms. A frame column is what the counterfactual write can reach, and a
# variable that is a bare name is such a column: `model.matrix()` reads it back
# at whatever value the write left there. A variable recorded under a call is
# not, whatever term it goes on to build. `z * x1` and `z:x1` contribute the
# variables `z` and `x1`, both plain columns, and their product is formed at
# design time from the column being set, so both routes agree on them and the
# ordinary way of writing an effect modification is left alone. `I(z * x1):x2`
# contributes the variable `I(z * x1)`, which is frozen in the frame at the
# values it was fit on however deep in an interaction its term sits.
#
# So every variable that parses to a call reading the exposure is rejected,
# whatever the call computes and whatever order the terms built from it carry,
# because nothing on this route recomputes it. The response is skipped: it is
# not a design column, and the outcome values are read once from the fit rather
# than rebuilt at each counterfactual value, so the write never reaches it.
check_ipw_exposure_call_terms <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  mod_terms <- stats::terms(outcome_mod)
  variables <- as.list(attr(mod_terms, "variables"))[-1]
  response <- attr(mod_terms, "response")

  if (isTRUE(response > 0)) {
    variables <- variables[-response]
  }

  for (variable in variables) {
    if (is.name(variable) || !exposure_name %in% all.vars(variable)) {
      next
    }

    abort_ipw_exposure_rebuild(
      exposure_name,
      deparse1(variable),
      fit_levels = NULL,
      rebuilt = FALSE,
      call = call
    )
  }

  invisible(TRUE)
}

# The halves of the rebuild rejection. The diagnosis is the same throughout, a
# term that reads the exposure and cannot be put back the way it was fit, and
# what differs is why the rebuild cannot do it and what fixes it. `fit_levels`
# is NULL for a term the fit recorded no levels for, which only the no-`.data`
# route rejects.
abort_ipw_exposure_rebuild <- function(
  exposure_name,
  term,
  fit_levels,
  rebuilt,
  call = rlang::caller_env()
) {
  route <- if (rebuilt) {
    c(
      x = "{.fun ipw} estimates the marginal means by setting \\
      {.val {exposure_name}} to one value at a time and rebuilding the outcome \\
      design from {.arg .data}, and {.code {term}} does not evaluate to the \\
      levels {.val {fit_levels}} it was fit with once the column is constant.",
      i = "Refit {.arg outcome_mod} on the plain {.val {exposure_name}} \\
      column, converting it to a factor in the data first if you want a factor \\
      coding."
    )
  } else if (!is.null(fit_levels)) {
    c(
      x = "Without {.arg .data} the designs come from {.arg outcome_mod}'s own \\
      model frame, which holds {.code {term}} at the values it was fit on, so \\
      the counterfactual value {.fun ipw} sets is ignored and every level is \\
      given the fitted design.",
      i = "Supply {.arg .data}, which rebuilds {.code {term}} at each value \\
      under the levels {.val {fit_levels}} it was fit with, or refit \\
      {.arg outcome_mod} on the plain {.val {exposure_name}} column."
    )
  } else {
    # The note about `*` and `:` is for a term that mixes the exposure with
    # something else, `I(z * x1)`, which is an interaction written the long way
    # and has a short way that needs no `.data`. Beside a term of the exposure
    # alone, `I(z^2)`, there is no interaction to write either way, and the note
    # left the reader looking for one they had not written.
    covariates <- setdiff(ipw_term_vars(term), exposure_name)

    interaction_hint <- if (length(covariates) > 0) {
      c(
        i = "An interaction written with {.code *} or {.code :} is formed from \\
        the frame's own columns and is rebuilt on either route, so it needs no \\
        {.arg .data}."
      )
    } else {
      NULL
    }

    c(
      x = "Without {.arg .data} the designs come from {.arg outcome_mod}'s own \\
      model frame, which holds {.code {term}} at the values it was fit on, so \\
      the counterfactual value {.fun ipw} sets never reaches it and the term \\
      contributes its fitted values to every level's design.",
      i = "Supply {.arg .data}, which recomputes {.code {term}} at each value \\
      {.fun ipw} sets, or refit {.arg outcome_mod} on the plain \\
      {.val {exposure_name}} column.",
      interaction_hint
    )
  }

  abort(
    c(
      "{.arg outcome_mod} must read the exposure from a column {.fun ipw} can \\
      set to counterfactual values.",
      x = "{.arg outcome_mod} reads {.val {exposure_name}} through the term \\
      {.code {term}}.",
      route
    ),
    error_class = "propensity_ipw_exposure_error",
    call = call
  )
}

# The variables a deparsed model term reads. A term that does not parse reads
# none this can name.
ipw_term_vars <- function(term) {
  parsed <- tryCatch(str2lang(term), error = function(e) NULL)

  if (is.null(parsed)) {
    return(character(0))
  }

  all.vars(parsed)
}

# Whether a `dataClasses` name, which is a deparsed model variable, reads the
# exposure. A name that does not parse is not a term this can judge.
ipw_term_reads <- function(term, exposure_name) {
  exposure_name %in% ipw_term_vars(term)
}

# Require both models to order the exposure's levels the same way. The
# counterfactual designs set the exposure to a factor of the propensity score
# model's levels, so their dummy columns follow that order, while the
# coefficients they multiply follow the outcome model's own coding. `X %*% beta`
# is positional and consults no names, so an outcome model fit on a releveled
# copy of the same exposure pairs each design column with the wrong coefficient
# and returns estimates that are wrong, and can be sign-flipped, with nothing
# signaled. Releveling is not the only way in: a character exposure column is
# factored alphabetically by `model.frame`, which disagrees with any fitted
# order that is not alphabetical.
#
# Read the levels from the fitted `xlevels` and never from `.data`. R drops
# levels that appear in the data but never in the fit, so an exposure declared
# with an extra unobserved level is analyzed correctly and its `xlevels` match;
# a check against the column's own `levels()` would see the extra level and
# reject that correct analysis.
#
# Silent when the outcome model records no levels under the exposure's name, as
# when the formula transforms it, and when the propensity score model carries no
# `lev` to compare against. Both are handled elsewhere.
check_ipw_outcome_exposure_levels <- function(
  outcome_mod,
  exposure_name,
  ps_lev,
  call = rlang::caller_env()
) {
  xlevels <- outcome_mod$xlevels

  if (is.null(ps_lev) || !exposure_name %in% names(xlevels)) {
    return(invisible(TRUE))
  }

  out_lev <- xlevels[[exposure_name]]

  if (!identical(out_lev, ps_lev)) {
    abort(
      c(
        "{.arg outcome_mod} and {.arg wt_mod} must code the exposure on the \\
        same levels in the same level order.",
        x = "{.arg outcome_mod} was fit on {.val {out_lev}}; {.arg wt_mod} was \\
        fit on {.val {ps_lev}}.",
        i = "Refit {.arg outcome_mod} with {.val {exposure_name}} factored in \\
        the propensity score model's order. A character column is factored \\
        alphabetically, so convert it to a factor with that order first."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Reject an outcome model whose family or link the estimating equations cannot
# stack. Both standard error paths classify the outcome model as either a
# gaussian identity linear score or a binomial score with the model's link, so
# only three shapes are supported: an lm, a gaussian glm with an identity link,
# and a binomial or quasibinomial glm with a link the outcome score reconstructs
# (the links ipw_inv_link handles). Any other family (poisson, quasipoisson,
# Gamma, inverse gaussian, ...) would otherwise be silently stacked as a binomial
# score, solving to a different root than the fitted model; a non-identity
# gaussian link would be silently treated as identity. Shared by every mestimation
# spec path and the linearization entry.
check_ipw_outcome_family <- function(outcome_mod, call = rlang::caller_env()) {
  # An lm (not a glm) is a gaussian identity linear model.
  if (!inherits(outcome_mod, "glm")) {
    return(invisible(TRUE))
  }

  family <- outcome_mod$family$family
  link <- outcome_mod$family$link
  supported_links <- c("logit", "probit", "cloglog", "log", "identity")

  if (identical(family, "gaussian")) {
    if (!identical(link, "identity")) {
      abort(
        c(
          "{.fun ipw} supports only an identity link for a gaussian outcome \\
          model.",
          x = "{.arg outcome_mod} is a gaussian model with a {.val {link}} link.",
          i = "Refit {.arg outcome_mod} with a gaussian identity link, or use a \\
          binomial or quasibinomial outcome model."
        ),
        error_class = "propensity_ipw_family_error",
        call = call
      )
    }
    return(invisible(TRUE))
  }

  if (family %in% c("binomial", "quasibinomial")) {
    if (!link %in% supported_links) {
      abort(
        c(
          "{.fun ipw} does not support a {.val {link}} link for a \\
          {.val {family}} outcome model.",
          x = "{.arg outcome_mod} was fit with a {.val {link}} link.",
          i = "Refit {.arg outcome_mod} with one of {.val {supported_links}}."
        ),
        error_class = "propensity_ipw_family_error",
        call = call
      )
    }
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.val {family}} outcome models.",
      x = "{.arg outcome_mod} was fit with the {.val {family}} family.",
      i = "Fit {.arg outcome_mod} with a binomial or quasibinomial family for a \\
      binary outcome, or a gaussian identity link (or an {.fun lm}) for a \\
      continuous outcome."
    ),
    error_class = "propensity_ipw_family_error",
    call = call
  )
}

# Validate the marginal structural model's link on the continuous exposure path,
# which the binary and categorical paths do not restrict this way. The check
# fires at spec-construction entry, before any M-estimator solve. The reported
# effect is labeled from the outcome link (identity for a slope, logit for a log
# odds ratio, log for a log risk ratio), so an outcome link outside that set has
# no label.
#
# The propensity score model's own link is the registry's business
# (R/ipw-continuous-models.R): it decides which score the stacked system carries
# rather than whether the fitted mean can be read back at all.
check_ipw_continuous_outcome_link <- function(
  outcome_mod,
  call = rlang::caller_env()
) {
  outcome_link <- if (inherits(outcome_mod, "glm")) {
    outcome_mod$family$link
  } else {
    "identity"
  }
  supported <- c("identity", "logit", "log")
  if (!outcome_link %in% supported) {
    abort(
      c(
        "{.fun ipw} does not support a {.val {outcome_link}} link for the \\
        marginal structural model of a continuous exposure.",
        x = "{.arg outcome_mod} was fit with a {.val {outcome_link}} link.",
        i = "Refit {.arg outcome_mod} with one of {.val {supported}}."
      ),
      error_class = "propensity_ipw_link_error",
      call = call
    )
  }

  invisible(TRUE)
}

# Restrict the linearization path to an outcome model of the exposure alone. Its
# point estimates come from g-computation on the fitted weighted outcome model,
# but the influence functions in linearize_variables_for_wts and correct_for_ps
# are those of the Hajek weighted-mean estimator. The two agree only when the
# outcome model regresses on the exposure by itself; any covariate, interaction,
# or transformed term makes the reported standard error that of a different
# estimator. The exposure-only condition is that the outcome model's term labels
# are exactly the exposure name. Detect an adjusted outcome model here and direct
# the user to the mestimation path, which stacks adjusted outcome models
# correctly. Applies only to the linearization path.
#
# The intercept is load-bearing for the same reason and is checked separately,
# because dropping it leaves the term labels untouched. Under a numeric coding
# such as `y ~ z - 1` the linear predictor at the unexposed level is zero for
# every unit, so the g-computation mean under no exposure is pinned at the link's
# zero point (0.5 under a logit or probit outcome link, 0 for a linear model)
# instead of being estimated. The reported estimates are then no longer the Hajek
# means the influence functions describe, and the standard errors describe an
# estimator that was never fit.
#
# The rejection is broader than that mechanism. A saturated factor coding such as
# `y ~ 0 + zf` does estimate both cell means, and the M-estimation path accepts
# it, but the linearization influence functions are derived for the intercept
# parameterization and this path requires it regardless. The error message is
# therefore worded as the fact of the missing intercept rather than as the
# numeric coding's consequence. The check is on the terms object, so it applies
# to lm and glm outcome models alike.
# `se_method` names the method the caller asked for, since the robust diagnostic
# reaches this guard through the same route and has the same requirement: it is
# the linearization route with the correction for having estimated the weights
# dropped, so it reports the same Hajek means and needs the same outcome model
# to report them from.
check_ipw_linearization_outcome <- function(
  outcome_mod,
  exposure_name,
  se_method = "linearization",
  call = rlang::caller_env()
) {
  outcome_terms <- stats::terms(outcome_mod)
  term_labels <- attr(outcome_terms, "term.labels")

  if (!identical(term_labels, exposure_name)) {
    # Two opposite models reach this guard: one that carries the exposure and
    # more, and one that carries the exposure not at all. `y ~ 1` has no terms
    # whatever, so reporting it as adjusted described a model the user did not
    # fit. A term that reads the exposure without being it, `I(z^2)` say, is
    # adjusted for the purposes of this path, since the influence functions are
    # derived for the bare exposure alone.
    reads_exposure <- vapply(
      term_labels,
      ipw_term_reads,
      logical(1),
      exposure_name = exposure_name
    )

    problem <- if (any(reads_exposure)) {
      "{.arg outcome_mod} is adjusted for terms beyond {.val {exposure_name}}."
    } else {
      "{.arg outcome_mod} does not include the exposure {.val {exposure_name}}."
    }

    abort(
      c(
        "{.fun ipw} supports {.val {se_method}} standard errors only for an \\
        outcome model of the exposure alone.",
        x = problem,
        i = "Use {.code se_method = \"mestimation\"} for a covariate-adjusted \\
        outcome model."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  if (!identical(attr(outcome_terms, "intercept"), 1L)) {
    abort(
      c(
        "{.fun ipw} supports {.val {se_method}} standard errors only for an \\
        outcome model with an intercept.",
        x = "{.arg outcome_mod} was fit without an intercept.",
        i = "This covers every no-intercept coding, including a saturated \\
        factor coding of the exposure, which the M-estimation path accepts.",
        i = "Include an intercept in {.arg outcome_mod}, or use \\
        {.code se_method = \"mestimation\"}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  invisible(TRUE)
}

#' @exportS3Method causalgenerics::ipw default
ipw.default <- function(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
) {
  rlang::check_dots_empty()
  abort(
    c(
      "{.fun ipw} does not know how to handle {.arg wt_mod} of class \\
      {.cls {class(wt_mod)}}.",
      i = "{.arg wt_mod} must be a fitted propensity score model: a \\
      {.cls glm} for a binary exposure, an {.cls lm} or gaussian \\
      {.cls glm} for a continuous exposure, or a {.cls multinom} for a \\
      categorical exposure."
    ),
    error_class = "propensity_method_error"
  )
}

# The reported rows of the linearization path, in the order the M-estimation
# path reports them: the counterfactual mean at each exposure level, reference
# level first, then the contrasts built from that pair. `exposure_levels` names
# the two levels, and every row is keyed by the level it belongs to or the pair
# it compares.
#
# The mean rows come from the same influence functions the contrast rows do:
# `l1` and `l0` are the influence functions of `mu1` and `mu0` themselves, so
# each mean's variance is the variance of its own influence function over n,
# which is the calculation the contrasts already apply to differences of the two.
#
# `se_method` decides only how an influence function is read into a variance,
# which the two helpers below settle. Everything else about the calculation is
# shared, the diagnostic being the same influence functions with the propensity
# correction dropped.
calculate_estimates <- function(
  lin_vars,
  marginal_means,
  n,
  linear_regression,
  conf_level,
  exposure_levels,
  se_method = "linearization"
) {
  z_val <- qnorm(1 - ((1 - conf_level) / 2))
  exact <- identical(se_method, "robust")

  mu_est <- c(marginal_means$mu0, marginal_means$mu1)
  mu_inf <- cbind(lin_vars$l0, lin_vars$l1)
  mu_se <- sqrt(apply(mu_inf, 2, influence_variance, n = n, exact = exact))
  mu_z <- mu_est / mu_se
  mu_contrast <- as.character(exposure_levels)
  contrast_label <- paste(mu_contrast[[2]], "vs", mu_contrast[[1]])

  ### RISK DIFFERENCE (raw scale)
  # --------------------------------
  # Influence = (l1 - l0)
  rd_inf <- lin_vars$l1 - lin_vars$l0
  rd_var <- influence_variance(rd_inf, n, exact)

  rd_est <- marginal_means$mu1 - marginal_means$mu0
  rd_se <- sqrt(rd_var)

  rd_ci_lower <- rd_est - z_val * rd_se
  rd_ci_upper <- rd_est + z_val * rd_se

  rd_z <- rd_est / rd_se

  # for continuous outcomes, only return difference
  if (linear_regression) {
    estimates <- ipw_linearization_estimates(
      effect = c("mean", "mean", "diff"),
      contrast = c(mu_contrast, contrast_label),
      estimate = c(mu_est, rd_est),
      std.err = c(mu_se, rd_se),
      z = c(mu_z, rd_z),
      ci.lower = c(mu_est - z_val * mu_se, rd_ci_lower),
      ci.upper = c(mu_est + z_val * mu_se, rd_ci_upper),
      conf_level = conf_level
    )

    return(ipw_attach_influence_vcov(
      estimates,
      cbind(mu_inf, rd_inf),
      n,
      exact
    ))
  }

  ### RISK RATIO (log scale)
  # ---------------------------
  # Risk ratio: RR = mu1 / mu0
  # We'll store 'log_rr_est = log(RR)' in the output.

  rr_raw_est <- marginal_means$mu1 / marginal_means$mu0
  log_rr_est <- log(rr_raw_est)

  # Influence function for log(RR) via the delta method:
  #   d log(mu1) - d log(mu0) = l1 / mu1 - l0 / mu0
  log_rr_inf <- lin_vars$l1 /
    marginal_means$mu1 -
    lin_vars$l0 / marginal_means$mu0

  log_rr_var <- influence_variance(log_rr_inf, n, exact)
  log_rr_se <- sqrt(log_rr_var)

  log_rr_ci_lower <- log_rr_est - z_val * log_rr_se
  log_rr_ci_upper <- log_rr_est + z_val * log_rr_se

  rr_z <- log_rr_est / log_rr_se

  ### ODDS RATIO (log scale)
  # ---------------------------
  # OR = [mu1/(1 - mu1)] / [mu0/(1 - mu0)]

  or_raw_est <- (marginal_means$mu1 / (1 - marginal_means$mu1)) /
    (marginal_means$mu0 / (1 - marginal_means$mu0))
  log_or_est <- log(or_raw_est)

  # Influence function for log(OR) via the delta method:
  #   d logit(mu1) - d logit(mu0)
  #     = l1 / [mu1*(1 - mu1)] - l0 / [mu0*(1 - mu0)]
  log_or_inf <- (lin_vars$l1 /
    (marginal_means$mu1 * (1 - marginal_means$mu1))) -
    (lin_vars$l0 / (marginal_means$mu0 * (1 - marginal_means$mu0)))

  log_or_var <- influence_variance(log_or_inf, n, exact)
  log_or_se <- sqrt(log_or_var)

  log_or_ci_lower <- log_or_est - z_val * log_or_se
  log_or_ci_upper <- log_or_est + z_val * log_or_se

  or_z <- log_or_est / log_or_se

  estimates <- ipw_linearization_estimates(
    effect = c("mean", "mean", "rd", "log(rr)", "log(or)"),
    contrast = c(mu_contrast, rep(contrast_label, 3)),
    # The mean rows and the risk difference are on the outcome's own scale; the
    # ratio measures are reported as logs of themselves.
    estimate = c(mu_est, rd_est, log_rr_est, log_or_est),
    # Each standard error is on the scale of the estimate beside it.
    std.err = c(mu_se, rd_se, log_rr_se, log_or_se),
    z = c(mu_z, rd_z, rr_z, or_z),
    ci.lower = c(
      mu_est - z_val * mu_se,
      rd_ci_lower,
      log_rr_ci_lower,
      log_or_ci_lower
    ),
    ci.upper = c(
      mu_est + z_val * mu_se,
      rd_ci_upper,
      log_rr_ci_upper,
      log_or_ci_upper
    ),
    conf_level = conf_level
  )

  ipw_attach_influence_vcov(
    estimates,
    cbind(mu_inf, rd_inf, log_rr_inf, log_or_inf),
    n,
    exact
  )
}

# The estimates frame of the linearization path, in the column order the
# M-estimation path builds: the columns that name a row lead, then the numbers.
# The p-value is derived here rather than passed, since it is a function of the
# test statistic alone and deriving it once keeps the two branches from
# disagreeing about which tail it comes from.
ipw_linearization_estimates <- function(
  effect,
  contrast,
  estimate,
  std.err,
  z,
  ci.lower,
  ci.upper,
  conf_level
) {
  data.frame(
    effect = effect,
    contrast = contrast,
    estimate = estimate,
    std.err = std.err,
    z = z,
    ci.lower = ci.lower,
    ci.upper = ci.upper,
    conf.level = conf_level,
    p.value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

# The covariance of the reported effects on the linearization path. There is no
# stacked system here, so it comes from the influence functions the standard
# errors themselves come from: each reported row's variance is the variance of
# its influence function over n, and the covariance of two of them is the
# covariance of theirs over the same n. The diagonal therefore reproduces the
# reported standard errors, and the off-diagonal records that the effect
# measures are transformations of the same pair of marginal means, which the
# table reports rows of their own.
#
# The block is labeled the way every surface of the result labels its rows,
# read off the frame it belongs to rather than off the influence matrix, whose
# columns are the influence functions and not the rows they describe.
ipw_attach_influence_vcov <- function(estimates, influence, n, exact = FALSE) {
  covariance <- influence_covariance(influence, n, exact)
  labels <- ipw_effect_labels(estimates)
  dimnames(covariance) <- list(labels, labels)
  attr(estimates, "ipw_vcov") <- covariance
  estimates
}

# How an influence function is read into the variance of the estimate it
# belongs to. The linearization route takes the sample variance of the influence
# values over n, the estimator's asymptotic variance read off the sample. The
# robust diagnostic reports the HC0 sandwich of the weighted outcome model
# exactly, whose meat is the plain sum of squares over n. The influence values
# already sum to zero over the sample, so the two readings differ by the factor
# n / (n - 1) and nothing else, and dividing by n is what makes the diagnostic
# the sandwich itself rather than a finite-sample rescaling of it.
influence_variance <- function(x, n, exact) {
  if (exact) {
    return(sum(x^2) / n^2)
  }

  var(x) / n
}

# The same two readings for a set of influence functions at once, giving the
# covariance the reported rows are read at.
influence_covariance <- function(influence, n, exact) {
  if (exact) {
    return(crossprod(influence) / n^2)
  }

  cov(influence) / n
}

# Mark a result whose standard errors are the diagnostic ones. The class leads
# the vector, so the result is an `ipw` to everything that does not care and the
# few surfaces that do can say so; a result fit by a method that accounts for
# the estimated weights is left exactly as it was built.
ipw_mark_diagnostic_se <- function(result) {
  if (!identical(result$se_method, "robust")) {
    return(result)
  }

  class(result) <- c("ipw_diagnostic_se", class(result))
  result
}

#' Print a result whose standard errors are a diagnostic
#'
#' `se_method = "robust"` reports the sandwich the weighted outcome model
#' computes for itself, which treats the estimated weights as known. Printing
#' such a result writes what [print()][base::print] writes for any [ipw()]
#' result and then one line naming the method, so that the number beside each
#' estimate is not read as one accounting for the propensity score model. See
#' **Standard errors as a diagnostic** in [ipw()] for what the diagnostic does
#' and does not cover.
#'
#' @param x An `ipw` result fit with `se_method = "robust"`.
#' @param ... Passed to the next method.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' y <- rbinom(n, 1, plogis(-0.3 + 0.7 * z + 0.5 * x))
#' dat <- data.frame(x, z, y)
#'
#' ps_mod <- glm(z ~ x, data = dat, family = binomial())
#' wts <- wt_ate(ps_mod)
#' outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
#'
#' print(ipw(ps_mod, outcome_mod, se_method = "robust"))
#'
#' @seealso [ipw()] for the estimator and the standard error methods.
#'
#' @export
print.ipw_diagnostic_se <- function(x, ...) {
  NextMethod()
  cat("\n")
  cat(paste0(
    "Standard errors: ",
    x$se_method,
    ", a diagnostic that treats the weights as known\n"
  ))

  invisible(x)
}

#' @rdname print.ipw_diagnostic_se
#'
#' @param row.names,optional Passed to the next method.
#'
#' @return `as.data.frame()` returns what it returns for any [ipw()] result,
#'   carrying an `ipw_se_diagnostic` attribute naming the method.
#'
#' @export
as.data.frame.ipw_diagnostic_se <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  # Forced here rather than passed on as an unevaluated argument: a promise
  # holding `NextMethod()` runs in whatever frame first forces it, and a
  # refusal raised inside it would then name that frame rather than the call
  # the user made.
  frame <- NextMethod()

  mark_ipw_se_diagnostic(frame, x)
}

#' @rdname print.ipw_diagnostic_se
#'
#' @return `tidy()` returns what it returns for any [ipw()] result, carrying an
#'   `ipw_se_diagnostic` attribute naming the method.
#'
#' @exportS3Method generics::tidy ipw_diagnostic_se
tidy.ipw_diagnostic_se <- function(x, ...) {
  # Forced for the reason the coercion above forces it.
  tidied <- NextMethod()

  # The tidier reaches the table through the coercion above, which marks it, but
  # a tidied table is a tibble built from that frame and an attribute is not one
  # of its columns. The mark is therefore reattached here rather than relied on
  # to survive the coercion.
  mark_ipw_se_diagnostic(tidied, x)
}

# The mark a table of a diagnostic result carries: the name of the method,
# rather than a flag, so that a caller reading it learns which diagnostic it
# was. It is an attribute rather than a column, the table being the one every
# other route reports.
mark_ipw_se_diagnostic <- function(table, result) {
  attr(table, "ipw_se_diagnostic") <- result$se_method
  table
}

# accounts for dependence introduced by weights
# treats them as fixed, not estimated
# equivalent to usual robust sandwich estimator
linearize_variables_for_wts <- function(Z, Y, wts, marginal_means) {
  n <- length(Z)
  wts <- as.double(wts)
  l1 <- n / marginal_means$n1 * (wts * Z * (Y - marginal_means$mu1))
  l0 <- n / marginal_means$n0 * (wts * (1 - Z) * (Y - marginal_means$mu0))

  list(l1 = l1, l0 = l0)
}


# additionally adds variability for the estimation of the PS
linearize_variables_for_ps <- function(
  exposure,
  outcome,
  wts,
  ps,
  estimand,
  weight_matrix,
  marginal_means,
  uncorrected_lin_vars,
  ps_link
) {
  n <- length(exposure)
  score_factor_f <- derive_score_factor(ps_link)

  # derive_weights() differentiates the unstabilized weight. A stabilizer
  # multiplies each weight by a factor that does not depend on the propensity
  # coefficients, so the derivative of the weight actually in use is the
  # unstabilized derivative times that same factor, observation by observation.
  weight_derivatives <- derive_weights(
    exposure = exposure,
    ps = ps,
    weight_matrix = weight_matrix,
    ps_link = ps_link,
    estimand = estimand
  ) *
    effective_stabilizer(wts, exposure, estimand)

  correction_mat <- compute_ps_correction_matrix_inv(
    weight_matrix = weight_matrix,
    ps = ps,
    ps_link = ps_link,
    n = n
  )

  l1 <- uncorrected_lin_vars$l1 +
    correct_for_ps(
      exposure = exposure,
      outcome = outcome,
      ps = ps,
      mu = marginal_means$mu1,
      n_group = marginal_means$n1,
      weight_matrix = weight_matrix,
      weight_derivatives = weight_derivatives,
      correction_mat = correction_mat,
      score_factor_f = score_factor_f,
      n = n
    )

  l0 <- uncorrected_lin_vars$l0 +
    correct_for_ps(
      exposure = exposure,
      exposure_actual = 1 - exposure,
      outcome = outcome,
      ps = ps,
      mu = marginal_means$mu0,
      n_group = marginal_means$n0,
      weight_matrix = weight_matrix,
      weight_derivatives = weight_derivatives,
      correction_mat = correction_mat,
      score_factor_f = score_factor_f,
      n = n
    )

  list(l1 = l1, l0 = l0)
}

compute_ps_correction_matrix_inv <- function(weight_matrix, ps, ps_link, n) {
  deriv_link_f <- derive_link(ps_link)

  # row-by-row x_i x_i^T
  crossprods <- t_tcrossprod_over_rows(weight_matrix)

  # build the correction matrix
  correction_mat <- matrix(
    colSums(crossprods / (ps * (1 - ps) * (deriv_link_f(ps)^2))) / n,
    ncol(weight_matrix),
    ncol(weight_matrix)
  )

  solve(correction_mat)
}

# faster than `t(apply(mat, 1, \(x) tcrossprod(x)))`
t_tcrossprod_over_rows <- function(mat) {
  n <- nrow(mat)
  p <- ncol(mat)
  out <- matrix(0, nrow = n, ncol = p * p)
  for (i in seq_len(n)) {
    out[i, ] <- tcrossprod(mat[i, ])
  }

  out
}

correct_for_ps <- function(
  exposure,
  exposure_actual = exposure,
  outcome,
  ps,
  mu,
  n_group,
  weight_matrix,
  weight_derivatives,
  correction_mat,
  score_factor_f,
  n
) {
  # first, compute partial-derivative sums over subjects (averaged by n)
  partial_derivative_sums <- colSums(
    weight_derivatives * exposure_actual * (outcome - mu)
  ) /
    n

  # then build the transformation matrix for correction. The score of a binomial
  # GLM with link g is x_i (Z_i - p_i) / (p_i (1 - p_i) g'(p_i)); the score factor
  # 1 / (p (1 - p) g'(p)) generalizes the correction to any link. It is exactly 1
  # for the canonical logit, leaving the logit influence values unchanged.
  transformation_mat <- correction_mat %*%
    t((exposure - ps) * score_factor_f(ps) * weight_matrix)

  # and then apply the partial-derivative sums to that transformation
  correction_contrib <- rbind(partial_derivative_sums) %*% transformation_mat

  # rescale by (n / n_group)
  scaling_factor <- n / n_group
  correction_contrib <- correction_contrib * scaling_factor

  # and reduce to vector and unname
  correction_contrib |>
    drop() |>
    unname()
}

estimate_marginal_means <- function(
  outcome_mod,
  wts,
  exposure,
  exposure_name,
  .data = NULL,
  call = rlang::caller_env()
) {
  # todo: this could be generalized with split() and lapply()
  if (is.null(.data)) {
    .data <- model.frame(outcome_mod)
    check_exposure(.data, exposure_name, call = call)
  }
  # todo: make this more flexible for different values and model specs
  # maybe can optionally accept a function for g-comp
  exposure_values <- sort(unique(exposure))

  if (!isTRUE(length(exposure_values) == 2)) {
    abort(
      c(
        "{.code ipw()} currently only supports binary exposures.",
        x = "There are {length(exposure_values)} unique value{?s} of the \\
        exposure."
      ),
      error_class = "propensity_ipw_exposure_error",
      call = call
    )
  }

  .data_1 <- .data
  .data_0 <- .data
  .data_1[[exposure_name]] <- exposure_values[[2]]
  .data_0[[exposure_name]] <- exposure_values[[1]]

  # Group totals are sums of the underlying numbers, so work from the bare
  # doubles. Masking the weights themselves would build a shorter psw that only
  # ever feeds sum(), and any metadata that cannot follow that subset would be
  # reported against an internal temporary the caller never sees.
  wt_values <- as.double(wts)

  n1 <- sum(wt_values[exposure == exposure_values[[2]]])
  mu1 <- mean(predict(outcome_mod, newdata = .data_1, type = "response"))

  n0 <- sum(wt_values[exposure == exposure_values[[1]]])
  mu0 <- mean(predict(outcome_mod, newdata = .data_0, type = "response"))

  list(
    # exposure = 1
    n1 = n1,
    mu1 = mu1,
    # exposure = 0
    n0 = n0,
    mu0 = mu0,
    # The levels the means are keyed to: `mu0` is the mean at the first and
    # `mu1` the mean at the second. Reporting them beside the means is what
    # lets the estimates be keyed by the pair the means were computed at,
    # rather than by a second sort of the exposure taken on trust to agree.
    levels = exposure_values
  )
}

derive_weights <- function(
  exposure,
  ps,
  weight_matrix,
  ps_link = c("logit", "probit", "cloglog"),
  estimand = c("ate", "att", "ato", "atm")
) {
  estimand <- rlang::arg_match(estimand)
  deriv_link_f <- derive_link(ps_link)

  deriv_vec <- switch(
    estimand,
    "ate" = ate_derivative(exposure, ps, deriv_link_f),
    "att" = att_derivative(exposure, ps, deriv_link_f),
    "ato" = ato_derivative(exposure, ps, deriv_link_f),
    "atm" = atm_derivative(exposure, ps, deriv_link_f)
  )

  weight_matrix * deriv_vec
}

ate_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    -1 / (ps^2 * deriv_link_f(ps)),
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
}

att_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    0,
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
}

ato_derivative <- function(exposure, ps, deriv_link_f) {
  ifelse(
    exposure == 1,
    -1 / deriv_link_f(ps),
    1 / deriv_link_f(ps)
  )
}

atm_derivative <- function(exposure, ps, deriv_link_f) {
  wv <- ifelse(
    exposure == 1,
    -1 / (ps^2 * deriv_link_f(ps)),
    1 / ((1 - ps)^2 * deriv_link_f(ps))
  )
  # Then set parts to zero
  wv[exposure == 1 & ps < 0.5] <- 0
  wv[exposure == 0 & ps > 0.5] <- 0
  wv
}


derive_link <- function(ps_link = c("logit", "probit", "cloglog")) {
  ps_link <- rlang::arg_match(ps_link)
  switch(
    ps_link,
    logit = function(x) 1 / (x * (1 - x)),
    probit = function(x) 1 / (dnorm(qnorm(x))),
    cloglog = function(x) -1 / ((1 - x) * log(1 - x))
  )
}

# Score factor 1 / (p (1 - p) g'(p)) that scales the residual (Z - p) into the
# binomial GLM score on the linearization ps correction. Returned in closed form
# per link so the logit factor is exactly 1 (not the round-trip
# p (1 - p) / (p (1 - p)) of derive_link), keeping the logit influence values
# bit-for-bit unchanged. Probit and cloglog match 1 / (p (1 - p) g'(p)) with
# g'(p) from derive_link().
derive_score_factor <- function(ps_link = c("logit", "probit", "cloglog")) {
  ps_link <- rlang::arg_match(ps_link)
  switch(
    ps_link,
    logit = function(p) rep(1, length(p)),
    probit = function(p) dnorm(qnorm(p)) / (p * (1 - p)),
    cloglog = function(p) -log(1 - p) / p
  )
}

# The per-observation factor a stabilizer multiplies the weights by, recovered
# from the weights themselves. Two forms reach here. An explicit
# `stabilization_score` is recorded on the weights and used as recorded; a scalar
# recycles. The default stabilizer records no score and instead scales the
# treated weights by the sample exposure mean and the untreated weights by its
# complement, so it is reconstructed from the same recoded exposure through
# `ipw_default_stab_seed()`, the seed `ate_binary()` builds those weights with.
# Treating that sample mean as fixed is exact rather than an approximation: each
# group constant scales that group's weights and that group's weight total
# identically, so it cancels from the Hajek ratio and contributes no influence
# term of its own. Weights carrying neither form take the identity factor, which
# covers every unstabilized analysis.
effective_stabilizer <- function(wts, exposure, estimand) {
  score <- stabilization_score(wts)
  if (!is.null(score)) {
    return(score)
  }

  if (is_stabilized(wts) && identical(estimand, "ate")) {
    p1 <- ipw_default_stab_seed(exposure)
    return(ifelse(exposure == 1, p1, 1 - p1))
  }

  1
}

# Report an outcome model whose fitting data can no longer be reached, the
# counterpart of the guard `ipw_extract_ps_design()` applies to the propensity
# model. The remedy is different and the difference is the point: `.data` rebuilds
# a propensity design, but it cannot stand in for the outcome model's frame,
# because the weights are read from that frame and `.data` carries covariates
# rather than weights. Offering `.data` here would send the user somewhere that
# cannot work, so the message says to refit and says why.
abort_outcome_frame_gone <- function(cause, call = rlang::caller_env()) {
  abort(
    c(
      "Can't reconstruct the data behind {.arg outcome_mod}.",
      x = "{cause}",
      i = "Refit {.arg outcome_mod} where its data is available, or fit it \\
      with {.code model = TRUE} so the model frame is kept.",
      i = "{.arg .data} cannot stand in here: the weights are read from \\
      {.arg outcome_mod}'s own model frame."
    ),
    error_class = "propensity_ipw_data_error",
    call = call
  )
}

# Every caller supplies the outcome model, so the guard above names it directly.
extract_weights <- function(.mod, call = rlang::caller_env()) {
  mf <- tryCatch(model.frame(.mod), error = function(e) e)

  if (inherits(mf, "error")) {
    abort_outcome_frame_gone(conditionMessage(mf), call = call)
  }

  model.weights(mf)
}

# Convert an outcome response to the 0/1 numeric coding a binomial glm uses. A
# factor maps every non-first level to success (glm treats the first level as
# failure), matching outcome_mod$y; a logical or numeric response doubles
# directly. Applied to the raw outcome vector at each extraction site (a model
# frame response or a .data column) so that a factor or logical response stacks
# and linearizes as its converted value rather than as level codes.
ipw_outcome_numeric <- function(outcome) {
  if (is.factor(outcome)) {
    as.double(outcome != levels(outcome)[[1L]])
  } else {
    as.double(outcome)
  }
}

# Recode a binary exposure to 0/1 for the influence-function computations, taking
# the second sorted value as the exposed group (coded 1): a factor's second
# level, TRUE for a logical, or the larger value for a numeric. Shared by the
# M-estimation spec and the linearization path so the two recode identically. A
# numeric 0/1 exposure is returned unchanged (sort(unique(c(0, 1)))[[2]] is 1).
# An exposure with more than two unique values is returned untouched so the
# downstream two-level check still fires on the original values.
ipw_recode_binary_exposure <- function(exposure) {
  exposure_values <- sort(unique(exposure))
  if (length(exposure_values) == 2) {
    as.double(exposure == exposure_values[[2]])
  } else {
    exposure
  }
}

check_estimand <- function(wts, estimand, call = rlang::caller_env()) {
  # Membership before the comparison below, and for both sources, so that the
  # value this returns is always one `ipw()` has a weight and a tilt for. A value
  # that names no estimand at all was reported by whatever noticed it first, and
  # no report said so: an argument checked against weights that record an
  # estimand came back as a disagreement between two estimands when only one of
  # them exists, and one checked against plain numeric weights was not compared
  # at all and reached the weighted means, where base R reported that `x` and `w`
  # had different lengths. An estimand read off the weights was checked against
  # nothing whatever, since `psw()` records the estimand it is given, so the same
  # base R report came back from a hand-built `psw()`.
  check_estimand_known(estimand, from_weights = FALSE, call = call)

  if (is_causal_wt(wts)) {
    estimand_from_weights <- estimand(wts)
  } else {
    estimand_from_weights <- NULL
  }

  check_estimand_known(estimand_from_weights, from_weights = TRUE, call = call)

  if (!is.null(estimand_from_weights) && !is.null(estimand)) {
    same_estimand <- identical(estimand_from_weights, estimand)
    if (!same_estimand) {
      .estimand <- estimand
      .estimand_from_weights <- estimand_from_weights
      abort(
        "Estimand in weights different from {.arg estimand}: \\
        {.val { .estimand_from_weights}} vs. {.val { .estimand}}",
        error_class = "propensity_ipw_estimand_error",
        call = call
      )
    } else {
      return(estimand)
    }
  }

  if (is.null(estimand_from_weights) && is.null(estimand)) {
    abort(
      c(
        "Can't determine the estimand from weights.",
        i = "Please specify {.arg estimand}."
      ),
      error_class = "propensity_ipw_estimand_error",
      call = call
    )
  }

  if (!is.null(estimand_from_weights)) {
    return(estimand_from_weights)
  } else {
    return(estimand)
  }
}

# Membership in the set of estimands `ipw()` knows a weight for, reported
# against the source the value came from. The two sources are different
# mistakes: an `estimand` argument is a typo the caller can correct in place,
# while an estimand recorded in the weights is a property of an object built
# earlier and is corrected where that object is built. `NULL` is not a failure
# here; it means the source is silent and the other one is asked.
check_estimand_known <- function(
  estimand,
  from_weights,
  call = rlang::caller_env()
) {
  if (is.null(estimand)) {
    return(invisible(TRUE))
  }

  # The type before the membership. `%in%` reads its left side through
  # `as.character()`, so a value that is not a string matches the name it prints
  # as and travels on in its own type: `list("ate")` reached the marginal means
  # as a list, where base R reported that `x` and `w` had different lengths, and
  # `factor("att")` reached the tilt, which is a `switch()` and reads a factor as
  # its integer level code, selecting a branch by position and returning some
  # other estimand's weights with only base R's note about the coercion to say
  # so. A value of the wrong type is reported as one, since a list holding
  # "ate" spells an estimand that exists and being told it is not one of them
  # would be false.
  spelled <- is.character(estimand) && length(estimand) == 1L

  if (spelled && estimand %in% ipw_estimands) {
    return(invisible(TRUE))
  }

  problem <- if (from_weights) {
    found <- if (spelled) {
      "Their estimand is {.val {estimand}}."
    } else {
      "Their estimand has class {.cls {class(estimand)}} and length \\
      {length(estimand)}; an estimand is a single string."
    }

    c(
      "The weights supplied to {.arg outcome_mod} record an estimand \\
      {.fun ipw} does not know.",
      x = found,
      i = "Valid estimands: {.val {ipw_estimands}}.",
      i = "Rebuild the weights with a weight function such as {.fun wt_ate}, \\
      or record one of those estimands in {.fun psw}."
    )
  } else {
    found <- if (spelled) {
      "{.val {estimand}} is not one of them."
    } else {
      "{.arg estimand} has class {.cls {class(estimand)}} and length \\
      {length(estimand)}; an estimand is a single string."
    }

    c(
      "{.arg estimand} must name an estimand {.fun ipw} knows.",
      x = found,
      i = "Valid estimands: {.val {ipw_estimands}}."
    )
  }

  abort(
    problem,
    error_class = "propensity_ipw_estimand_error",
    call = call
  )
}

check_exposure <- function(.data, .exposure_name, call = rlang::caller_env()) {
  assert_class(.exposure_name, "character", .length = 1, call = call)

  if (.exposure_name %in% names(.data)) {
    return(invisible(TRUE))
  }

  # A model frame is keyed by the terms a model was fit with rather than by the
  # variables inside them, so a transformed or a basis exposure sits under the
  # term's own name and no column holds the exposure. Naming that column is what
  # separates a model that transforms the exposure, which `.data` fixes, from
  # one that never reads it, which nothing here can.
  #
  # The name is rebound without its dot, since cli reads a substitution opening
  # with one as inline markup rather than as a value.
  exposure <- .exposure_name
  holding <- columns_reading(names(.data), exposure)

  problem <- if (length(holding)) {
    "The frame holds {.code {holding}} instead, which {?is a term/are terms} \\
    built from it."
  } else {
    "No column of the frame is built from it either."
  }

  abort(
    c(
      "{.val {exposure}} not found in {.code model.frame(outcome_mod)}.",
      x = problem,
      i = "The outcome model may have transformations in the formula, and a \\
      model frame records the term rather than the variables inside it.",
      i = "Please specify {.arg .data}"
    ),
    call = call,
    error_class = "propensity_columns_exist_error"
  )
}

# The column names of a frame that are built from a variable. Names are parsed
# rather than matched as text, so a name mentioning the variable inside a longer
# one is not mistaken for a term reading it. A name that does not parse is not a
# term and reads nothing.
columns_reading <- function(column_names, variable) {
  reads <- vapply(
    column_names,
    function(name) {
      parsed <- tryCatch(str2lang(name), error = function(e) NULL)
      !is.null(parsed) && variable %in% all.vars(parsed)
    },
    logical(1)
  )

  column_names[reads]
}

is_linear_regression <- function(outcome_mod) {
  # Check if 'outcome_mod' inherits from 'glm' first
  if (inherits(outcome_mod, "glm")) {
    # 'outcome_mod' is a glm, check if family is 'gaussian'
    if (family(outcome_mod)$family == "gaussian") {
      # It's a glm with Gaussian family => linear regression
      return(TRUE)
    } else {
      # It's a glm but not Gaussian
      return(FALSE)
    }
  }

  # Otherwise, if it inherits from 'lm' (and NOT from 'glm'),
  # this should be a plain linear model
  # since we already confirmed the class of the model
  TRUE
}
