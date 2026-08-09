# .by refuses linearization standard errors

    Code
      expr
    Condition <propensity_ipw_by_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` with "linearization" standard errors.
      x The stratum effects and their contrasts are parameters of the stacked estimating equations, and the linearization path solves no such system.
      i Use `se_method = "mestimation"` to report effects by `.by`.

# .by refuses a continuous exposure

    Code
      expr
    Condition <propensity_ipw_by_exposure_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for continuous exposures.
      x A continuous exposure reports the marginal structural model's own exposure coefficient rather than a set of standardized means, so there is no effect within a subgroup for it to report.
      i Omit `.by` to report the overall effect.
      i Fitting each subgroup on its own subset reports a coefficient per subgroup, but no covariance between them, so the difference between two subgroups cannot be tested from those fits.

# .by refuses a selection that does not name exactly one modifier

    Code
      expr
    Condition <propensity_ipw_by_arg_error>
      Error in `ipw()`:
      ! `.by` must name exactly one modifier.
      x It names 2 columns.
      i Effects are reported within the levels of a single variable. Cross two variables into one column, with `interaction()`, and name that column instead.

# .by refuses a modifier with missing values

    Code
      expr
    Condition <propensity_ipw_by_missing_error>
      Error in `ipw()`:
      ! `.by` must name a modifier with no missing values.
      x "w" has 1 missing value.
      i A missing value names no subgroup, so the units carrying one belong to none of the strata the effects would be reported within.
      i Drop those rows and refit both models, or recode the missing values as a level of their own.

# .by refuses a numeric modifier

    Code
      expr
    Condition <propensity_ipw_by_type_error>
      Error in `ipw()`:
      ! `.by` must name a factor or a character modifier.
      x "x1" is a numeric vector.
      i The effects are reported within the levels of the modifier, so it has to name a fixed set of subgroups.
      i Cut a continuous variable into groups with `cut()`, or convert a logical or a numeric code to a factor, and refit both models on that column.

# .by refuses a stratum missing an exposure level

    Code
      expr
    Condition <propensity_ipw_by_levels_error>
      Error in `ipw()`:
      ! `.by` must name a modifier whose subgroups each hold every exposure level.
      x "g = b" holds no unit with "z" set to "0".
      i An effect within a subgroup contrasts the exposure levels inside it, so a subgroup missing one of them has no contrast to report there.
      i Use a coarser modifier, one whose subgroups each hold every exposure level. Refitting either model does not help: the data hold no comparison there.

# .by warns when the outcome model has no exposure-by-modifier term

    Code
      out <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
    Condition <propensity_ipw_by_interaction_warning>
      Warning in `ipw()`:
      `outcome_mod` has no term reading both "z" and "v".
      i The stratum effects are g-computation on `outcome_mod` as it was specified, so a model with no such term forces one and the same effect on every subgroup.
      i Add `z:v` to `outcome_mod` and refit it to let the effect differ across the levels of "v".

