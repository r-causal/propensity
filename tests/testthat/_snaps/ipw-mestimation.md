# ipw_spec_binary error names the unsupported outcome family

    Code
      ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    Condition
      Error:
      ! `ipw()` does not support "poisson" outcome models.
      x `outcome_mod` was fit with the "poisson" family.
      i Fit `outcome_mod` with a binomial or quasibinomial family for a binary outcome, or a gaussian identity link (or an `lm()`) for a continuous outcome.

