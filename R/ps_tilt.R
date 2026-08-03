# The estimand tilting function h and its exported interface. Everywhere the
# package needs a tilt it reads it from here: the binary and categorical weight
# formulas in R/weights.R, and the weight closures and marginal-mean rows of the
# M-estimation psi builders in R/ipw-psi.R.
#
# `ps_tilt_binary()` and `ps_tilt_categorical()` are the arithmetic and validate
# nothing. The psi builders evaluate them once per solver step on propensity
# scores the solver is free to push to exactly zero or one, so the checking
# `ps_tilt()` performs for a user belongs at the exported surface rather than in
# the arithmetic.

#' Propensity score tilting functions
#'
#' @description
#' Every estimand this package supports targets a population whose covariate
#' distribution is the study population reweighted by a tilting function
#' \eqn{h} of the propensity score. `ps_tilt()` evaluates \eqn{h} at each
#' observation's propensity score.
#'
#' The tilt is the numerator of every weight: a weight is \eqn{h} divided by the
#' propensity score of the exposure level the unit actually received. It is also
#' what standardizes a g-computation estimate to a target population, which for
#' `"atm"`, `"ato"`, and `"entropy"` is the only route to the estimand, since
#' those populations are not subsets of the rows and cannot be reached by
#' filtering.
#'
#' @details
#' # Tilting functions
#'
#' For a binary exposure with propensity score \eqn{e = P(Z = \text{focal} \mid
#' X)}:
#'
#' | estimand | \eqn{h(e)} | target population |
#' | --- | --- | --- |
#' | `"ate"` | \eqn{1} | everyone |
#' | `"att"` | \eqn{e} | the focal group |
#' | `"atu"` | \eqn{1 - e} | the reference group |
#' | `"atm"` | \eqn{\min(e, 1 - e)} | the evenly matchable |
#' | `"ato"` | \eqn{e(1 - e)} | the overlap population |
#' | `"entropy"` | \eqn{-e \log(e) - (1 - e) \log(1 - e)} | the entropy-tilted population |
#'
#' For a categorical exposure with propensity score vector \eqn{(e_1, \ldots,
#' e_K)} and focal level \eqn{f}: `"ate"` is \eqn{1}, `"att"` is \eqn{e_f},
#' `"atu"` is \eqn{1 - e_f}, `"atm"` is \eqn{\min_k e_k}, `"ato"` is
#' \eqn{(\sum_k 1 / e_k)^{-1}}, and `"entropy"` is \eqn{-\sum_k e_k \log(e_k)}.
#'
#' Stabilization and censoring weights are not tilts. Stabilization multiplies a
#' weight by a marginal quantity that does not depend on the covariates, and
#' `wt_cens()` reuses the `"ate"` formula, so neither has an entry here.
#'
#' # Propensity score range
#'
#' Every propensity score in `ps` must lie strictly inside \eqn{(0, 1)}, the
#' same requirement [wt_ate()] and the rest of the weight family impose before
#' they tilt. The bound holds for a matrix or data frame `ps` entry by entry,
#' and each row must sum to one on top of it. A fitted model that separates the
#' exposure can return a probability of exactly zero or one; those scores have
#' no weight to divide and are rejected here rather than tilted.
#'
#' A missing propensity score gives a missing tilt under every estimand,
#' `"ate"` included, so an observation whose propensity score is unknown never
#' counts toward a tilted mean. A numeric `ps` propagates `NA` position by
#' position, and a matrix `ps` gives `NA` for any row holding one: a probability
#' vector with a missing entry is not one this can tilt on, whichever level the
#' tilt reads.
#'
#' # Modified propensity scores
#'
#' `ps_tilt()` takes plain propensity scores. A score modified by [ps_trim()],
#' [ps_trunc()], or [ps_calibrate()] carries a class of its own and has no
#' method here; pass the scores underneath it, with `as.numeric(x)` for a binary
#' exposure or `as.matrix(x)` for a categorical one. Units [ps_trim()] set to
#' `NA` stay `NA` through that extraction and take an `NA` tilt.
#'
#' @param ps Propensity scores. A numeric vector of \eqn{P(Z = \text{focal} \mid
#'   X)} for a binary exposure, or a matrix or data frame with one column per
#'   exposure level, named for that level, for a categorical exposure.
#' @param estimand One of `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or
#'   `"entropy"`.
#' @param ... These dots are for future extensions and must be empty.
#' @param .focal_level The exposure level the `"att"` and `"atu"` tilts target,
#'   matched against the column names of `ps`. A column named for the level
#'   exactly wins; failing that, a `.pred_` prefix is stripped, so `.pred_a`
#'   matches the level `a`. Required for those two estimands with a categorical
#'   `ps`, and accepted nowhere else: a numeric `ps` is already the probability
#'   of the focal level, and the remaining tilts treat every level alike.
#'
#' @return A plain double vector, unnamed, with one element per observation: the
#'   length of `ps` for the numeric method, and the number of rows of `ps` for
#'   the matrix and data frame methods.
#'
#' @seealso [wt_ate()] and the rest of the weight family, which divide the tilt
#'   by the propensity score of the received exposure level.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.6 * x))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + x))
#' sim <- data.frame(x = x, z = z, y = y)
#' ps <- unname(
#'   predict(glm(z ~ x, data = sim, family = binomial()), type = "response")
#' )
#'
#' # a weight is the tilt over the propensity of the received exposure level
#' received <- z * ps + (1 - z) * (1 - ps)
#' all.equal(
#'   as.numeric(wt_ato(ps, z, exposure_type = "binary")),
#'   ps_tilt(ps, "ato") / received
#' )
#'
#' # tilted g-computation: standardize counterfactual predictions to the
#' # overlap population with an h-weighted mean, no weights in the outcome model
#' fit <- glm(y ~ z * x, data = sim, family = binomial())
#' m1 <- predict(fit, transform(sim, z = 1), type = "response")
#' m0 <- predict(fit, transform(sim, z = 0), type = "response")
#'
#' h <- ps_tilt(ps, "ato")
#' weighted.mean(m1, h) - weighted.mean(m0, h)
#'
#' # the same predictions standardized to everyone give the ATE instead
#' mean(m1) - mean(m0)
#'
#' @export
ps_tilt <- function(ps, estimand, ..., .focal_level = NULL) {
  rlang::check_dots_empty()

  if (missing(estimand)) {
    abort(
      c(
        "{.arg estimand} must be supplied.",
        i = "It is one of {.val {ipw_estimands}}."
      ),
      error_class = "propensity_missing_arg_error"
    )
  }

  UseMethod("ps_tilt")
}

#' @export
ps_tilt.default <- function(ps, estimand, ..., .focal_level = NULL) {
  abort(
    c(
      "No method for objects of class {.cls {class(ps)}}.",
      i = "{.fun ps_tilt} takes plain propensity scores. A score modified by \\
      {.fun ps_trim}, {.fun ps_trunc}, or {.fun ps_calibrate} carries a class \\
      of its own; pass the scores underneath it, with {.code as.numeric(x)} \\
      for a binary exposure or {.code as.matrix(x)} for a categorical one."
    ),
    error_class = "propensity_method_error"
  )
}

#' @export
ps_tilt.numeric <- function(ps, estimand, ..., .focal_level = NULL) {
  estimand <- check_tilt_estimand(estimand)
  check_tilt_focal_unused(.focal_level, estimand, matrix_ps = FALSE)
  check_ps_range(ps)

  ps <- as.double(ps)

  mask_missing_tilt(ps_tilt_binary(ps, estimand), is.na(ps))
}

#' @export
ps_tilt.matrix <- function(ps, estimand, ..., .focal_level = NULL) {
  tilt_from_matrix(ps, estimand, .focal_level, call = rlang::current_env())
}

#' @export
ps_tilt.data.frame <- function(ps, estimand, ..., .focal_level = NULL) {
  tilt_from_matrix(
    as.matrix(ps),
    estimand,
    .focal_level,
    call = rlang::current_env()
  )
}

# The matrix and data frame methods differ only in the coercion, and both report
# their errors against the `ps_tilt()` call the user wrote.
tilt_from_matrix <- function(ps, estimand, .focal_level, call) {
  estimand <- check_tilt_estimand(estimand, call = call)
  ps <- check_tilt_ps_matrix(ps, call = call)
  focal_idx <- resolve_tilt_focal(ps, estimand, .focal_level, call = call)

  tilt <- ps_tilt_categorical(ps, estimand, focal_idx)

  mask_missing_tilt(tilt, rowSums(is.na(ps)) > 0)
}

# The seam every exported method returns through. It strips the names the
# arithmetic picks up from a row-named matrix or a named vector, so a tilt is a
# plain double whatever it was computed from, and it makes a missing propensity
# score give a missing tilt under every estimand. ate is what needs the second
# part: its tilt is the constant 1 and would otherwise report a value for an
# observation whose propensity score is unknown, which for a score trimmed to NA
# would silently count a unit the analysis dropped. The tilts that read the
# scores mostly carry the NA through on their own; masking covers att and atu on
# a matrix as well, where a missing entry outside the focal column still leaves
# a row that is not a probability vector.
mask_missing_tilt <- function(tilt, missing) {
  tilt <- unname(as.double(tilt))
  tilt[missing] <- NA_real_

  tilt
}

# ---- the tilting functions --------------------------------------------------

ps_tilt_binary <- function(e, estimand) {
  switch(
    estimand,
    ate = rep(1, length(e)),
    att = e,
    atu = 1 - e,
    atm = pmin(e, 1 - e),
    ato = e * (1 - e),
    entropy = {
      # The inverse link saturates to exactly 0 or 1 for extreme linear
      # predictors, where the tilt would evaluate 0 * log(0). Nudge the
      # saturated scores off the boundary as the categorical tilt does.
      e_safe <- e
      e_safe[e == 0] <- .Machine$double.eps
      e_safe[e == 1] <- 1 - .Machine$double.eps
      -e_safe * log(e_safe) - (1 - e_safe) * log(1 - e_safe)
    }
  )
}

# `ps` is the n-by-K propensity score matrix and `focal_idx` the focal column
# for att and atu.
ps_tilt_categorical <- function(ps, estimand, focal_idx = NULL) {
  switch(
    estimand,
    ate = rep(1, nrow(ps)),
    att = ps[, focal_idx],
    atu = 1 - ps[, focal_idx],
    ato = 1 / rowSums(1 / ps),
    atm = do.call(pmin, lapply(seq_len(ncol(ps)), function(j) ps[, j])),
    entropy = {
      ps_safe <- ps
      ps_safe[ps == 0] <- .Machine$double.eps
      -rowSums(ps_safe * log(ps_safe))
    }
  )
}

# ---- validation -------------------------------------------------------------

check_tilt_estimand <- function(estimand, call = rlang::caller_env()) {
  if (!rlang::is_string(estimand) || !estimand %in% ipw_estimands) {
    abort(
      c(
        "{.arg estimand} must be one of {.val {ipw_estimands}}.",
        i = "Stabilization and censoring weights are not tilts of the \\
        propensity score."
      ),
      error_class = "propensity_unknown_estimand_error",
      call = call
    )
  }

  estimand
}

# The exposure-free half of check_ps_matrix(): a categorical propensity score
# matrix has to be a numeric matrix of at least two columns whose rows are
# probability vectors, whatever exposure it goes on to describe.
check_tilt_ps_matrix <- function(ps, call = rlang::caller_env()) {
  if (!is.numeric(ps)) {
    abort(
      "{.arg ps} must be a numeric matrix or data frame.",
      call = call,
      error_class = "propensity_matrix_type_error"
    )
  }

  if (ncol(ps) < 2) {
    abort(
      c(
        "{.arg ps} must have one column per exposure level.",
        i = "It has {ncol(ps)} column{?s}."
      ),
      call = call,
      error_class = "propensity_matrix_dims_error"
    )
  }

  check_ps_matrix_rowsums(ps, call = call)
  check_ps_matrix_range(ps, arg = "ps", call = call)

  ps
}

# The column of `ps` the att and atu tilts read, or NULL for the tilts that
# treat every level alike.
resolve_tilt_focal <- function(
  ps,
  estimand,
  .focal_level,
  call = rlang::caller_env()
) {
  if (!estimand %in% c("att", "atu")) {
    check_tilt_focal_unused(
      .focal_level,
      estimand,
      matrix_ps = TRUE,
      call = call
    )
    return(NULL)
  }

  check_focal_level_required(.focal_level, estimand, call = call)

  columns <- colnames(ps)
  if (is.null(columns)) {
    abort(
      c(
        "{.arg ps} must have column names to resolve {.arg .focal_level}.",
        i = "Name each column for the exposure level whose probability it \\
        holds."
      ),
      call = call,
      error_class = "propensity_matrix_names_error"
    )
  }

  focal <- .focal_level
  if (length(focal) != 1) {
    abort(
      c(
        "{.arg .focal_level} must name a single column of {.arg ps}.",
        x = "It is length {length(focal)}."
      ),
      call = call,
      error_class = "propensity_tilt_focal_error"
    )
  }

  # A column named for the level exactly wins over one that only matches after
  # the prefix comes off, so a matrix carrying both `.pred_a` and `a` resolves
  # to `a` rather than to whichever sits first.
  matches <- which(columns == as.character(focal))
  if (!length(matches)) {
    matches <- which(sub("^\\.pred_", "", columns) == as.character(focal))
  }

  if (length(matches) > 1) {
    tied <- columns[matches]
    abort(
      c(
        "{.arg .focal_level} must name a single column of {.arg ps}.",
        x = "{.val {focal}} matches {.val {tied}}."
      ),
      call = call,
      error_class = "propensity_tilt_focal_error"
    )
  }

  if (!length(matches)) {
    abort(
      c(
        "{.arg .focal_level} must name a single column of {.arg ps}.",
        x = "{.val {focal}} is not one of {.val {columns}}."
      ),
      call = call,
      error_class = "propensity_tilt_focal_error"
    )
  }

  matches
}

check_focal_level_required <- function(
  .focal_level,
  estimand,
  call = rlang::caller_env()
) {
  if (is.null(.focal_level)) {
    abort(
      "Focal category must be specified for {toupper(estimand)} with \\
      categorical exposures.",
      error_class = "propensity_focal_required_error",
      call = call
    )
  }

  invisible(.focal_level)
}

check_tilt_focal_unused <- function(
  .focal_level,
  estimand,
  matrix_ps,
  call = rlang::caller_env()
) {
  if (is.null(.focal_level)) {
    return(invisible(NULL))
  }

  reason <- if (matrix_ps) {
    "The {.val {estimand}} tilt treats every exposure level alike, so it has \\
    no focal level."
  } else {
    "A numeric {.arg ps} is already the probability of the focal level, so no \\
    binary tilt takes one."
  }

  abort(
    c(
      "{.arg .focal_level} applies to the {.val att} and {.val atu} tilts of a \\
      categorical propensity score matrix.",
      x = reason
    ),
    error_class = "propensity_tilt_focal_error",
    call = call
  )
}
