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

# Refuse `.by` for a continuous exposure, the one exposure type whose path
# reports no counterfactual contrast to modify. A binary or a categorical fit
# reports contrasts of standardized means, and a stratum's effect is those means
# taken over the stratum's units. A continuous exposure has no such means:
# `ipw()` reports the marginal structural model's own exposure coefficient,
# which is a single number for the whole sample with no per-stratum counterpart
# the stacked blocks could build.
check_ipw_by_exposure <- function(.by, call = rlang::caller_env()) {
  if (ipw_by_absent(.by)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} for continuous exposures.",
      x = "A continuous exposure reports the marginal structural model's own \\
      exposure coefficient rather than a set of standardized means, so there \\
      is no effect within a subgroup for it to report.",
      i = "Omit {.arg .by} to report the overall effect.",
      i = "Fitting each subgroup on its own subset reports a coefficient per \\
      subgroup, but no covariance between them, so the difference between two \\
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
  exposure_levels,
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
    exposure_levels,
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

# Require every stratum to hold every exposure level. The effects within a
# stratum are contrasts of the counterfactual means taken over that stratum's
# units, and a stratum in which nobody took some level identifies neither the
# mean there nor any contrast against it. The outcome model still predicts at
# that level, so the fit would return a number: an extrapolation from the strata
# that do hold it, reported as though it had been estimated.
#
# Counting the distinct values a stratum holds is not enough once the exposure
# has more than two levels, since a stratum can hold two of three and still be
# missing a contrast. What is required is membership, level by level.
#
# Only the first incomplete stratum is named, the way
# `first_identical_design_pair()` reports the first offending pair: every fix is
# a coarser modifier, after which the check runs again over the new strata.
# Neither model can be refit to supply a comparison the data do not hold.
check_ipw_by_levels <- function(
  indicators,
  exposure,
  exposure_levels,
  exposure_name,
  labels,
  call = rlang::caller_env()
) {
  exposure <- as.character(exposure)
  exposure_levels <- as.character(exposure_levels)

  absent <- lapply(
    seq_along(labels),
    function(s) setdiff(exposure_levels, exposure[indicators[, s] == 1])
  )
  incomplete <- lengths(absent) > 0L

  if (!any(incomplete)) {
    return(invisible(TRUE))
  }

  first <- which(incomplete)[[1L]]
  stratum <- labels[[first]]
  missing_levels <- absent[[first]]

  abort(
    c(
      "{.arg .by} must name a modifier whose subgroups each hold every \\
      exposure level.",
      x = "{.val {stratum}} holds no unit with {.val {exposure_name}} set to \\
      {.val {missing_levels}}.",
      i = "An effect within a subgroup contrasts the exposure levels inside \\
      it, so a subgroup missing one of them has no contrast to report there.",
      i = "Use a coarser modifier, one whose subgroups each hold every \\
      exposure level. Refitting either model does not help: the data hold no \\
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

# The stratum blocks are written against an exposure of any arity. `preds` is
# one counterfactual prediction vector per exposure level, in the order the
# stratum mu block stores them; `reference` is the position in it of the level
# every contrast is taken against, and `comparisons` the positions of the levels
# compared to it, in the order the estimates table reports them.
#
# A binary spec stores its pair exposed-first, matching the whole-sample mu
# block it grew out of, so its reference is the second element and its single
# comparison the first. A categorical spec stores its levels reference-first, so
# its reference is the first element and its comparisons the rest. Stating both
# as positions is what lets one set of blocks serve the two.

# The suffix each contrast seed is named by, and the label its out-of-domain
# reports are made under: a subgroup on its own where the exposure names no
# contrasts, and the contrast within the subgroup where it does.
ipw_by_contrast_suffixes <- function(contrast_labels, group) {
  if (is.null(contrast_labels)) {
    return(group)
  }

  paste(contrast_labels, group)
}

# Plug-in init values for the stratum blocks, seeded the way the whole-sample
# blocks are: each stratum mean at the tilt-weighted mean of the counterfactual
# predictions over that stratum, which is the exact root of its psi row, and
# each contrast at the value those means imply. `mu_names` is the whole-sample
# mu block's naming, which each stratum's block reuses under its own label.
ipw_init_by <- function(
  spec,
  preds,
  reference,
  comparisons,
  mu_names,
  contrast_labels,
  tilt
) {
  by <- spec$by

  if (is.null(by)) {
    return(list(mu = numeric(0), contrast = numeric(0)))
  }

  # One mean per exposure level per stratum, stratum-major, so the pair or the
  # tuple belonging to one stratum sits together, which is the order the psi
  # rows and the layout index read them in.
  mu <- lapply(seq_along(by$labels), function(s) {
    w <- ipw_by_stratum_weight(by$indicators, s, tilt)
    vapply(preds, function(pred) stats::weighted.mean(pred, w), numeric(1))
  })

  mu_block <- stats::setNames(
    unlist(mu, use.names = FALSE),
    unlist(lapply(by$labels, function(label) paste0(mu_names, "_", label)))
  )

  stratum <- lapply(seq_along(by$labels), function(s) {
    suffixes <- ipw_by_contrast_suffixes(contrast_labels, by$labels[[s]])
    unlist(lapply(seq_along(comparisons), function(j) {
      ipw_init_contrasts(
        by$contrasts,
        mu[[s]][[comparisons[[j]]]],
        mu[[s]][[reference]],
        suffix = suffixes[[j]]
      )
    }))
  })

  em <- lapply(seq_along(by$em_labels), function(m) {
    suffixes <- ipw_by_contrast_suffixes(contrast_labels, by$em_labels[[m]])
    stats::setNames(
      stratum[[m + 1L]] - stratum[[1L]],
      unlist(lapply(suffixes, function(s) paste0(by$contrasts, "_", s)))
    )
  })

  list(mu = mu_block, contrast = c(unlist(stratum), unlist(em)))
}

# One out-of-domain reporter per stratum and comparison, in the order the
# contrast rows visit them, so a stratum whose marginal means leave the domain
# of a transform names the subgroup, and the contrast within it, that did.
# Built once for the fit rather than once per evaluation, so each undefined
# effect reports once however often the solver revisits it.
#
# The effect-modification rows are differences of the stratum contrast
# parameters and have no domain of their own, so they need none.
ipw_by_reporters <- function(by, contrast_labels) {
  if (is.null(by)) {
    return(NULL)
  }

  labels <- unlist(lapply(by$labels, function(group) {
    ipw_by_contrast_suffixes(contrast_labels, group)
  }))

  lapply(labels, ipw_contrast_reporter)
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
  preds,
  reference,
  comparisons,
  tilt,
  indicators,
  forms,
  th_mu,
  th_con,
  reporters
) {
  n <- length(preds[[1L]])
  n_levels <- length(preds)
  n_strata <- ncol(indicators)
  n_comparisons <- length(comparisons)
  n_forms <- length(forms)
  # The contrast rows one group contributes: contrast-major and effect-minor,
  # the order an ungrouped fit already reports its whole-sample block in.
  block <- n_comparisons * n_forms

  mu_rows <- matrix(0, nrow = n_levels * n_strata, ncol = n)
  con_rows <- matrix(0, nrow = length(th_con), ncol = n)

  for (s in seq_len(n_strata)) {
    w <- ipw_by_stratum_weight(indicators, s, tilt)
    mu_at <- (s - 1L) * n_levels

    for (level in seq_len(n_levels)) {
      mu_rows[mu_at + level, ] <- w * (preds[[level]] - th_mu[[mu_at + level]])
    }

    for (j in seq_len(n_comparisons)) {
      mu_hi <- th_mu[[mu_at + comparisons[[j]]]]
      mu_lo <- th_mu[[mu_at + reference]]
      reporter <- reporters[[(s - 1L) * n_comparisons + j]]

      for (f in seq_len(n_forms)) {
        row <- (s - 1L) * block + (j - 1L) * n_forms + f
        value <- ipw_contrast_value(
          forms[[f]],
          mu_hi,
          mu_lo,
          reporter = reporter
        )
        con_rows[row, ] <- value - th_con[[row]]
      }
    }
  }

  for (m in seq_len(n_strata - 1L)) {
    for (within in seq_len(block)) {
      row <- (n_strata + m - 1L) * block + within
      value <- th_con[[m * block + within]] - th_con[[within]]
      con_rows[row, ] <- value - th_con[[row]]
    }
  }

  list(mu = mu_rows, contrast = con_rows)
}

# The identity columns of the reported rows: the effect measure, the level or
# the pair of levels the contrast column names, and the subgroup where `.by`
# names one. The counterfactual mean at each exposure level leads, one row per
# level in level order with the reference level first, and the contrasts built
# from those means follow. Without `.by` no subgroup is named at all.
#
# With `.by` the whole-sample block is followed by the means within each
# stratum, then one contrast block per stratum and one per non-reference
# stratum against the reference one, each repeating the whole-sample block's
# contrast-major, effect-minor order over the measures a subgroup reports. A
# stratum contrast compares two strata's effects and has no mean of its own, so
# the mean rows stop at the strata themselves.
ipw_estimate_rows <- function(spec) {
  if (!is.null(spec$joint)) {
    return(ipw_joint_estimate_rows(spec$joint))
  }

  levels <- spec$exposure_levels
  means <- rep("mean", length(levels))
  contrast_labels <- ipw_estimate_contrast_labels(spec)
  effect <- rep(spec$contrasts, times = ipw_n_contrasts(spec))
  contrast <- rep(contrast_labels, each = length(spec$contrasts))

  if (is.null(spec$by)) {
    return(list(
      effect = c(means, effect),
      contrast = c(levels, contrast),
      group = NULL
    ))
  }

  strata <- spec$by$labels
  groups <- c(strata, spec$by$em_labels)
  by_effect <- rep(spec$by$contrasts, times = ipw_n_contrasts(spec))
  by_contrast <- rep(contrast_labels, each = length(spec$by$contrasts))

  list(
    effect = c(
      means,
      effect,
      rep(means, times = length(strata)),
      rep(by_effect, times = length(groups))
    ),
    contrast = c(
      levels,
      contrast,
      rep(levels, times = length(strata)),
      rep(by_contrast, times = length(groups))
    ),
    group = c(
      rep(ipw_overall_group, length(means) + length(effect)),
      rep(strata, each = length(means)),
      rep(groups, each = length(by_effect))
    )
  )
}

# The pair of levels each contrast row compares. A binary exposure compares one
# pair and a categorical exposure one pair per non-reference level, and both are
# written the same way, so a reader of either table finds the comparison named
# rather than left to the reader of a binary one to infer.
ipw_estimate_contrast_labels <- function(spec) {
  if (identical(spec$exposure_type, "categorical")) {
    return(ipw_contrast_labels(spec))
  }

  levels <- spec$exposure_levels
  paste(levels[[2]], "vs", levels[[1]])
}
