# Effect modification for the M-estimation path of ipw(). `.by` names a
# modifier, and the result reports the effects it already reported for the whole
# sample, then one set of effects within each stratum of that modifier, then one
# set contrasting each non-reference stratum against the reference one.
# Everything here is unexported and splits three ways: the refusals that decide
# whether a request can be answered at all, the resolution of the modifier into
# the strata the stacked blocks in R/ipw-psi.R are built from, and the labels the
# estimates table keys its rows by.

# The group label the rows estimated over the whole sample carry. Every other
# row names a stratum, so the whole-sample rows need a name of their own for the
# column to name a subgroup in every row.
ipw_overall_group <- "overall"

# Whether a captured `.by` names nothing. That is the default, and every path
# outside this file has to behave exactly as it did before `.by` existed when it
# holds.
ipw_by_absent <- function(.by) {
  is.null(.by) || rlang::quo_is_null(.by)
}

# ---- refusals that do not need the modifier ---------------------------------

# Refuse `.by` on the linearization path. The stratum means and the contrasts
# built from them are parameters of the stacked system, which is how their
# standard errors come out of the same sandwich as everything else; the
# linearization path solves no such system and its influence functions are
# derived for the whole-sample Hajek means alone.
check_ipw_by_method <- function(.by, se_method, call = rlang::caller_env()) {
  if (ipw_by_absent(.by) || !identical(se_method, "linearization")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} with {.val linearization} \\
      standard errors.",
      x = "The stratum effects and their contrasts are parameters of the \\
      stacked estimating equations, and the linearization path solves no such \\
      system.",
      i = "Use {.code se_method = \"mestimation\"} to report effects by \\
      {.arg .by}."
    ),
    error_class = "propensity_ipw_by_method_error",
    call = call
  )
}

# Refuse `.by` for an exposure whose path reports no counterfactual contrast to
# modify. A continuous exposure reports the marginal structural model's own
# coefficient rather than a pair of standardized means, and a categorical
# exposure reports one contrast per non-reference level, neither of which the
# stratum blocks are built for.
check_ipw_by_exposure <- function(
  .by,
  exposure_type,
  call = rlang::caller_env()
) {
  if (ipw_by_absent(.by)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} for {.val {exposure_type}} \\
      exposures.",
      x = "Effect modification is reported for a binary exposure alone.",
      i = "Omit {.arg .by} to report the overall effect.",
      i = "Fitting each subgroup on its own subset reports the same stratum \\
      effects, but no covariance between them, so the difference between two \\
      subgroups cannot be tested from those fits."
    ),
    error_class = "propensity_ipw_by_exposure_error",
    call = call
  )
}

# ---- resolving the modifier -------------------------------------------------

# The strata a `.by` request names, or NULL when it names nothing. `data` is the
# frame the modifier is read from: `.data` when the caller supplied one, and the
# outcome model's own model frame otherwise. `exposure` is the 0/1 recode, which
# the estimability check reads, and `contrasts` is the set of effect measures
# the whole-sample rows report, which the stratum measures are taken from.
#
# The refusals run in the order a request fails in rather than in the order the
# pieces are computed: what the modifier is has to be settled before what its
# strata hold, and both before the outcome model is inspected for a term letting
# the effect differ across them. That last one is a warning, so a request this
# function will not answer is refused rather than diagnosed first.
ipw_resolve_by <- function(
  .by,
  data,
  exposure,
  exposure_name,
  outcome_mod,
  contrasts,
  call = rlang::caller_env()
) {
  if (ipw_by_absent(.by)) {
    return(NULL)
  }

  name <- ipw_by_column(.by, data, call = call)
  values <- data[[name]]

  check_ipw_by_missing(values, name, call = call)
  check_ipw_by_type(values, name, call = call)

  # A level no unit carries names an empty stratum, which has no mean to
  # standardize and no contrast to report. The fits drop unused levels
  # themselves, so dropping them here is the same reading of the modifier the
  # models were fit under.
  values <- droplevels(as.factor(values))
  levels <- levels(values)
  labels <- paste0(name, " = ", levels)

  indicators <- matrix(
    0,
    nrow = length(values),
    ncol = length(levels),
    dimnames = list(NULL, levels)
  )
  for (s in seq_along(levels)) {
    indicators[, s] <- as.numeric(values == levels[[s]])
  }

  check_ipw_by_levels(
    indicators,
    exposure,
    exposure_name,
    labels,
    call = call
  )
  check_ipw_by_interaction(outcome_mod, exposure_name, name, call = call)

  list(
    name = name,
    levels = levels,
    labels = labels,
    # Each non-reference stratum against the reference one, which is the first
    # level of the modifier.
    em_labels = if (length(labels) > 1L) {
      paste(labels[-1], "vs", labels[[1]])
    } else {
      character(0)
    },
    indicators = indicators,
    # An odds ratio is noncollapsible, so a difference of two of them is not the
    # difference in effect it reads as. The whole-sample rows keep it and the
    # stratum rows report the collapsible measures alone.
    contrasts = setdiff(contrasts, "log(or)")
  )
}

# The one column of `data` a `.by` request names. Selection is tidyselect's, so
# a bare name, a string, and a name held in a variable all reach the same
# column, and a column that is not there is reported by tidyselect in the terms
# the caller wrote.
ipw_by_column <- function(.by, data, call = rlang::caller_env()) {
  loc <- tidyselect::eval_select(
    .by,
    data,
    allow_rename = FALSE,
    error_call = call
  )

  if (!identical(length(loc), 1L)) {
    abort(
      c(
        "{.arg .by} must name exactly one modifier.",
        x = "It names {length(loc)} column{?s}.",
        i = "Effects are reported within the levels of a single variable. \\
        Cross two variables into one column, with {.fun interaction}, and name \\
        that column instead."
      ),
      error_class = "propensity_ipw_by_arg_error",
      call = call
    )
  }

  names(loc)
}

check_ipw_by_missing <- function(values, name, call = rlang::caller_env()) {
  n_missing <- sum(is.na(values))

  if (n_missing == 0L) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg .by} must name a modifier with no missing values.",
      x = "{.val {name}} has {n_missing} missing value{?s}.",
      i = "A missing value names no subgroup, so the units carrying one belong \\
      to none of the strata the effects would be reported within.",
      i = "Drop those rows and refit both models, or recode the missing values \\
      as a level of their own."
    ),
    error_class = "propensity_ipw_by_missing_error",
    call = call
  )
}

check_ipw_by_type <- function(values, name, call = rlang::caller_env()) {
  if (is.factor(values) || is.character(values)) {
    return(invisible(TRUE))
  }

  supplied <- ipw_class_article(values)

  abort(
    c(
      "{.arg .by} must name a factor or a character modifier.",
      x = "{.val {name}} is {supplied}.",
      i = "The effects are reported within the levels of the modifier, so it \\
      has to name a fixed set of subgroups.",
      i = "Cut a continuous variable into groups with {.fun cut}, or convert a \\
      logical or a numeric code to a factor, and refit both models on that \\
      column."
    ),
    error_class = "propensity_ipw_by_type_error",
    call = call
  )
}

# The article-bearing description of a column's type, for a message naming the
# type that was supplied. `ipw_class_articles` covers the model-frame classes;
# anything else is named by its class.
ipw_class_article <- function(x) {
  class_name <- stats::.MFclass(x)

  if (class_name %in% names(ipw_class_articles)) {
    return(ipw_class_articles[[class_name]])
  }

  paste("an object of class", class(x)[[1L]])
}

# Require every stratum to hold both exposure levels. The effect within a
# stratum is the contrast of the two counterfactual means taken over that
# stratum's units, and a stratum in which nobody took one of the levels
# identifies neither the mean there nor the contrast against it. The outcome
# model still predicts at both levels, so the fit would return a number: an
# extrapolation from the strata that do hold both, reported as though it had
# been estimated. Neither model can be refit to supply what is not there, so the
# remedy is a coarser modifier.
check_ipw_by_levels <- function(
  indicators,
  exposure,
  exposure_name,
  labels,
  call = rlang::caller_env()
) {
  incomplete <- vapply(
    seq_along(labels),
    function(s) length(unique(exposure[indicators[, s] == 1])) < 2L,
    logical(1)
  )

  if (!any(incomplete)) {
    return(invisible(TRUE))
  }

  bad <- labels[incomplete]

  abort(
    c(
      "{.arg .by} must name a modifier whose subgroups each hold both \\
      exposure levels.",
      x = "Every unit in {.val {bad}} has the same value of \\
      {.val {exposure_name}}.",
      i = "An effect within a subgroup contrasts the exposure levels inside \\
      it, so a subgroup holding one of them has no contrast to report.",
      i = "Use a coarser modifier, one whose subgroups each hold both exposure \\
      levels. Refitting either model does not help: the data hold no \\
      comparison there."
    ),
    error_class = "propensity_ipw_by_levels_error",
    call = call
  )
}

# Report an outcome model with no term reading the exposure and the modifier
# together. The stratum effects are g-computation on the model as specified, so
# such a model predicts the same contrast in every stratum up to the covariate
# distribution each one has, and the effect modification it reports is a
# property of the model rather than of the data. It is a warning rather than an
# error because the fit is still an honest reading of the model that was
# supplied.
check_ipw_by_interaction <- function(
  outcome_mod,
  exposure_name,
  name,
  call = rlang::caller_env()
) {
  term_labels <- attr(stats::terms(outcome_mod), "term.labels")
  interacts <- any(vapply(
    term_labels,
    function(term) {
      vars <- ipw_term_vars(term)
      exposure_name %in% vars && name %in% vars
    },
    logical(1)
  ))

  if (interacts) {
    return(invisible(TRUE))
  }

  term <- paste0(exposure_name, ":", name)

  warn(
    c(
      "{.arg outcome_mod} has no term reading both {.val {exposure_name}} and \\
      {.val {name}}.",
      i = "The stratum effects are g-computation on {.arg outcome_mod} as it \\
      was specified, so a model with no such term forces one and the same \\
      effect on every subgroup.",
      i = "Add {.code {term}} to {.arg outcome_mod} and refit it to let the \\
      effect differ across the levels of {.val {name}}."
    ),
    warning_class = "propensity_ipw_by_interaction_warning",
    call = call
  )

  invisible(FALSE)
}

# ---- the rows a resolved `.by` contributes ----------------------------------

# The weights one stratum's marginal means standardize over: the estimand's tilt
# restricted to that stratum. `tilt` is NULL for ate, whose tilt is the constant
# one, leaving the stratum indicator on its own. The seed and the psi rows both
# read it from here, because a seed that weighted the predictions differently
# from the row it seeds would not sit at that row's root.
ipw_by_stratum_weight <- function(indicators, s, tilt) {
  if (is.null(tilt)) indicators[, s] else tilt * indicators[, s]
}

# Plug-in init values for the stratum blocks, seeded the way the whole-sample
# blocks are: each stratum mean at the tilt-weighted mean of the counterfactual
# predictions over that stratum, which is the exact root of its psi row, and
# each contrast at the value those means imply.
ipw_init_by <- function(spec, pred1, pred0, tilt) {
  by <- spec$by

  if (is.null(by)) {
    return(list(mu = numeric(0), contrast = numeric(0)))
  }

  weights <- lapply(
    seq_along(by$labels),
    function(s) ipw_by_stratum_weight(by$indicators, s, tilt)
  )
  mu1 <- vapply(weights, function(w) stats::weighted.mean(pred1, w), numeric(1))
  mu0 <- vapply(weights, function(w) stats::weighted.mean(pred0, w), numeric(1))

  # Interleaved so that the pair belonging to one stratum sits together, which
  # is the order the psi rows and the layout index read them in.
  mu_block <- stats::setNames(
    as.vector(rbind(mu1, mu0)),
    as.vector(rbind(paste0("mu1_", by$labels), paste0("mu0_", by$labels)))
  )

  stratum <- lapply(seq_along(by$labels), function(s) {
    ipw_init_contrasts(
      by$contrasts,
      mu1[[s]],
      mu0[[s]],
      suffix = by$labels[[s]]
    )
  })
  em <- lapply(seq_along(by$em_labels), function(j) {
    stats::setNames(
      stratum[[j + 1L]] - stratum[[1L]],
      paste0(by$contrasts, "_", by$em_labels[[j]])
    )
  })

  list(
    mu = mu_block,
    contrast = c(unlist(stratum), unlist(em, use.names = TRUE))
  )
}

# The stratum blocks of one psi evaluation. The marginal-mean rows are the
# whole-sample rows restricted to a stratum: the root of
# sum_i h(e_i) I(v_i = s) (pred_a(x_i) - mu_a_s) = 0 is the tilt-weighted mean of
# the counterfactual predictions over the units of stratum s.
#
# The contrast rows come in two groups. The stratum contrasts transform that
# stratum's pair of means the way the whole-sample contrasts transform theirs.
# The effect-modification rows are read off the stratum contrast parameters
# rather than recomputed from the means, so each one is the difference of two
# parameters the system already carries and its derivative is exact.
ipw_by_psi_rows <- function(
  pred1,
  pred0,
  tilt,
  indicators,
  forms,
  th_mu,
  th_con,
  reporters
) {
  n <- length(pred1)
  n_strata <- ncol(indicators)
  n_forms <- length(forms)

  mu_rows <- matrix(0, nrow = 2L * n_strata, ncol = n)
  con_rows <- matrix(0, nrow = length(th_con), ncol = n)

  for (s in seq_len(n_strata)) {
    w <- ipw_by_stratum_weight(indicators, s, tilt)
    mu1 <- th_mu[[2L * s - 1L]]
    mu0 <- th_mu[[2L * s]]
    mu_rows[2L * s - 1L, ] <- w * (pred1 - mu1)
    mu_rows[2L * s, ] <- w * (pred0 - mu0)

    for (f in seq_len(n_forms)) {
      row <- (s - 1L) * n_forms + f
      value <- ipw_contrast_value(
        forms[[f]],
        mu1,
        mu0,
        reporter = reporters[[s]]
      )
      con_rows[row, ] <- value - th_con[[row]]
    }
  }

  for (j in seq_len(n_strata - 1L)) {
    for (f in seq_len(n_forms)) {
      row <- (n_strata + j - 1L) * n_forms + f
      value <- th_con[[j * n_forms + f]] - th_con[[f]]
      con_rows[row, ] <- value - th_con[[row]]
    }
  }

  list(mu = mu_rows, contrast = con_rows)
}

# The effect measure and the subgroup each reported row carries. Without `.by`
# the reported rows are the contrast block and no subgroup is named at all,
# which is the frame every result reported before `.by` existed. With it the
# block is followed by one set of measures within each stratum and one set for
# each non-reference stratum against the reference one, and every row names the
# subgroup it belongs to.
ipw_estimate_rows <- function(spec) {
  effect <- rep(spec$contrasts, times = ipw_n_contrasts(spec))

  if (is.null(spec$by)) {
    return(list(effect = effect, group = NULL))
  }

  groups <- c(spec$by$labels, spec$by$em_labels)

  list(
    effect = c(effect, rep(spec$by$contrasts, times = length(groups))),
    group = c(
      rep(ipw_overall_group, length(effect)),
      rep(groups, each = length(spec$by$contrasts))
    )
  )
}
