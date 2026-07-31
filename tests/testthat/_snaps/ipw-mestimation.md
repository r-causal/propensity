# ipw_spec_binary error names the unsupported outcome family

    Code
      ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    Condition
      Error:
      ! `ipw()` does not support "poisson" outcome models.
      x `outcome_mod` was fit with the "poisson" family.
      i Fit `outcome_mod` with a binomial or quasibinomial family for a binary outcome, or a gaussian identity link (or an `lm()`) for a continuous outcome.

# the missing-exposure error names the exposure and directs the user to refit

    Code
      ipw(mods$ps_mod, no_exposure, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must contain the exposure "z".
      x "z" does not appear in the formula of `outcome_mod`.
      i Refit `outcome_mod` with "z" on the right-hand side of the formula.

