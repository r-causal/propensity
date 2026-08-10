# Reporting for a joint exposure: two treatments intervened on at once.
#
# `causalgenerics::joint_exposure()` crosses two discrete treatments into one
# categorical exposure and records the crossing on the vector. The result is a
# factor, so the whole categorical path already fits it and would report what it
# reports for any exposure over those cells: each cell against the reference
# cell. Those rows are correct and they are not the question a joint
# intervention asks. What the declaration buys is a surface written in the two
# treatments rather than in the cells, which is what everything here builds:
#
#   * the counterfactual risk of each cell as a row of its own, under the effect
#     label "mean", keyed by the cell it belongs to;
#   * the simple effects, each treatment's effect within a fixed level of the
#     other, including the comparisons that are not against the reference cell
#     and which the vs-reference reporting cannot express at all;
#   * the interaction, the difference between two of the first treatment's
#     simple effects, on each collapsible scale.
#
# The cell means are already parameters of the stacked system, as the
# categorical mu block, so the surface adds no means: it reports that block and
# replaces the vs-reference contrast block with contrasts written over the same
# means. The simple effects are contrasts of two mu parameters, and the
# interaction rows are differences of two simple-effect parameters, so the
# double difference is exact by construction rather than by arithmetic that has
# to agree to the last bits.

# The plan a declared exposure is reported under, or NULL when the exposure
# declares no crossing. `forms` is the set of effect measures the contrast rows
# report, which is the whole-sample set with the odds ratio removed: an odds
# ratio is noncollapsible, so neither a simple effect reported beside one nor a
# difference of two of them says what it appears to.
#
# The plan is positions and labels and nothing numeric, so the seed, the psi
# rows, and the estimates table all read the one description of the surface
# rather than three that have to agree.
ipw_joint_plan <- function(exposure, forms) {
  if (!causalgenerics::is_joint_exposure(exposure)) {
    return(NULL)
  }

  components <- causalgenerics::joint_components(exposure)
  component_names <- names(components)
  component_levels <- unname(components)
  sizes <- lengths(component_levels)

  # The crossing varies the first component fastest, which is what puts the
  # reference cell first and is the order `joint_cell_labels()` builds upstream.
  cell_at <- function(first, second) first + (second - 1L) * sizes[[1L]]

  simple <- list()
  for (component in 1:2) {
    other <- 3L - component
    for (held in seq_len(sizes[[other]])) {
      for (level in seq_len(sizes[[component]])[-1]) {
        cells <- if (component == 1L) {
          c(cell_at(level, held), cell_at(1L, held))
        } else {
          c(cell_at(held, level), cell_at(held, 1L))
        }
        simple[[length(simple) + 1L]] <- list(
          contrast = ipw_joint_contrast_label(
            component_names[[component]],
            component_levels[[component]],
            level
          ),
          group = paste0(
            component_names[[other]],
            " = ",
            component_levels[[other]][[held]]
          ),
          hi = cells[[1L]],
          lo = cells[[2L]]
        )
      }
    }
  }

  # The interaction is reported once, under the first component's framing,
  # rather than twice. Interaction is symmetric in the two treatments, so the
  # difference between the first's simple effects is the difference between the
  # second's, and reporting both would be one quantity under two labels.
  #
  # Each row differences two entries of the first component's block, which is
  # ordered with the held level outer and the compared level inner.
  per_held <- sizes[[1L]] - 1L
  interaction <- list()
  for (held in seq_len(sizes[[2L]])[-1]) {
    for (level in seq_len(sizes[[1L]])[-1]) {
      interaction[[length(interaction) + 1L]] <- list(
        contrast = ipw_joint_contrast_label(
          component_names[[1L]],
          component_levels[[1L]],
          level
        ),
        group = paste0(
          component_names[[2L]],
          " = ",
          component_levels[[2L]][[held]],
          " vs ",
          component_names[[2L]],
          " = ",
          component_levels[[2L]][[1L]]
        ),
        hi = (held - 1L) * per_held + (level - 1L),
        lo = level - 1L
      )
    }
  }

  list(
    cells = levels(exposure),
    forms = forms,
    simple = simple,
    interaction = interaction
  )
}

# How a treatment's own contrast is written: the treatment, then the level it is
# set to against its reference level. The cells carry `"name = level"` labels
# already, so the colon is what keeps a contrast over one treatment from reading
# as a cell of the crossing.
ipw_joint_contrast_label <- function(name, levels, level) {
  paste0(name, ": ", levels[[level]], " vs ", levels[[1L]])
}

# The label each simple effect's out-of-domain reports are made under. The
# interaction rows difference two parameters the system already carries and have
# no domain of their own, so they need none.
ipw_joint_reporters <- function(joint) {
  if (is.null(joint)) {
    return(NULL)
  }

  lapply(
    vapply(
      joint$simple,
      function(entry) paste(entry$contrast, entry$group),
      character(1)
    ),
    ipw_contrast_reporter
  )
}

# Plug-in init values for the joint contrast block: each simple effect at the
# value the seeded cell means imply, and each interaction at the difference of
# the two simple effects it is built from, which is the exact root of its row.
ipw_init_joint <- function(joint, mu) {
  simple <- lapply(joint$simple, function(entry) {
    ipw_init_contrasts(
      joint$forms,
      mu[[entry$hi]],
      mu[[entry$lo]],
      suffix = paste(entry$contrast, entry$group)
    )
  })

  interaction <- lapply(joint$interaction, function(entry) {
    stats::setNames(
      simple[[entry$hi]] - simple[[entry$lo]],
      paste0(joint$forms, "_", paste(entry$contrast, entry$group))
    )
  })

  c(unlist(simple), unlist(interaction))
}

# The joint contrast rows of one psi evaluation. A simple effect transforms the
# pair of cell means it compares, the way the vs-reference contrasts transform
# theirs. An interaction row is read off the two simple-effect parameters rather
# than recomputed from the four means, so it is the difference of two parameters
# the system already carries: its derivative is exact, and the row equals the
# double difference by construction rather than to within a solver residual.
ipw_joint_psi_rows <- function(joint, th_mu, th_con, n, reporters) {
  n_forms <- length(joint$forms)
  rows <- matrix(0, nrow = length(th_con), ncol = n)

  for (s in seq_along(joint$simple)) {
    entry <- joint$simple[[s]]
    for (f in seq_len(n_forms)) {
      row <- (s - 1L) * n_forms + f
      value <- ipw_contrast_value(
        joint$forms[[f]],
        th_mu[[entry$hi]],
        th_mu[[entry$lo]],
        reporter = reporters[[s]]
      )
      rows[row, ] <- value - th_con[[row]]
    }
  }

  simple_rows <- length(joint$simple) * n_forms
  for (i in seq_along(joint$interaction)) {
    entry <- joint$interaction[[i]]
    for (f in seq_len(n_forms)) {
      row <- simple_rows + (i - 1L) * n_forms + f
      value <- th_con[[(entry$hi - 1L) * n_forms + f]] -
        th_con[[(entry$lo - 1L) * n_forms + f]]
      rows[row, ] <- value - th_con[[row]]
    }
  }

  rows
}

# The identity columns of the rows a joint fit reports: the cell means first,
# each keyed by its cell and belonging to no subgroup, then the simple effects
# keyed by the treatment they contrast and the level the other is held at, then
# the interaction rows keyed by the same treatment and the two levels of the
# other that are being compared.
ipw_joint_estimate_rows <- function(joint) {
  entries <- c(joint$simple, joint$interaction)
  n_cells <- length(joint$cells)

  list(
    effect = c(
      rep("mean", n_cells),
      rep(joint$forms, times = length(entries))
    ),
    contrast = c(
      joint$cells,
      rep(
        vapply(entries, function(e) e$contrast, character(1)),
        each = length(joint$forms)
      )
    ),
    group = c(
      rep(ipw_overall_group, n_cells),
      rep(
        vapply(entries, function(e) e$group, character(1)),
        each = length(joint$forms)
      )
    )
  )
}

# ---- refusals ---------------------------------------------------------------

# Refuse a focal estimand on a declared joint exposure. The surface standardizes
# every cell mean to one population, and a focal estimand names one cell as that
# population. The cell means would then be the risks among the units treated at
# one corner of the crossing, which is a coherent quantity and not the one the
# simple effects and the interaction rows are built to combine, so reporting it
# under these labels would say something the numbers do not.
check_ipw_joint_estimand <- function(
  exposure,
  estimand,
  call = rlang::caller_env()
) {
  declared <- causalgenerics::is_joint_exposure(exposure)

  if (!declared || identical(estimand, "ate")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} reports a joint exposure for the {.val {'ate'}} estimand \\
      only.",
      x = "The weights were built for the {.val {estimand}} estimand.",
      i = "Every cell mean on the joint surface standardizes to one \\
      population, and a tilted estimand standardizes each of them to a \\
      population the simple effects and the interaction are not defined over.",
      i = "Refit the weights for {.val {'ate'}}, or drop the declaration with \\
      {.code factor(x)} to report each cell against the reference cell under \\
      the estimand you have."
    ),
    error_class = "propensity_ipw_joint_estimand_error",
    call = call
  )
}

# Refuse `.by` on a declared joint exposure. Effect modification of a joint
# intervention is a three-way question, the interaction between two treatments
# within the levels of a third variable, and this surface reports neither that
# nor a projection of it that could stand in.
#
# `joint` says whether the fit is over a crossing at all, as a flag rather than
# as the thing the flag was read off: a declared exposure carries the crossing
# on the vector, and a `joint_wt_models()` container is the crossing without one
# to read.
#
# `remedy` is the one bullet the two routes do not share. Both offer the same
# way out, reporting the crossing as a plain categorical exposure, but they
# reach it differently: a declared exposure is a factor already and gives up its
# declaration when it is asked to be one, while a container has no declaration
# to drop and the crossing has to be built and weighted as a factor instead.
check_ipw_joint_by <- function(joint, .by, remedy, call = rlang::caller_env()) {
  if (!joint || ipw_by_absent(.by)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} for a joint exposure.",
      x = "A joint exposure already reports an interaction between two \\
      treatments, and reporting it within the levels of a modifier is a \\
      three-way question this surface does not answer.",
      i = remedy
    ),
    error_class = "propensity_ipw_joint_by_error",
    call = call
  )
}
