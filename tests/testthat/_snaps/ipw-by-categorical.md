# .by refuses a stratum missing one of a categorical exposure's levels

    Code
      expr
    Condition <propensity_ipw_by_levels_error>
      Error in `ipw()`:
      ! `.by` must name a modifier whose subgroups each hold every exposure level.
      x "g = b" holds no unit with "a" set to "c".
      i An effect within a subgroup contrasts the exposure levels inside it, so a subgroup missing one of them has no contrast to report there.
      i Use a coarser modifier, one whose subgroups each hold every exposure level. Refitting either model does not help: the data hold no comparison there.

# .by warns when a categorical outcome model has no exposure-by-modifier term

    Code
      out <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
    Condition <propensity_ipw_by_interaction_warning>
      Warning in `ipw()`:
      `outcome_mod` has no term reading both "a" and "v".
      i The stratum effects are g-computation on `outcome_mod` as it was specified, so a model with no such term forces one and the same effect on every subgroup.
      i Add `a:v` to `outcome_mod` and refit it to let the effect differ across the levels of "v".

