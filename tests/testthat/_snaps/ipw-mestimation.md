# ipw_spec_binary error names the unsupported outcome family

    Code
      ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    Condition
      Error:
      ! `ipw()` does not support "poisson" outcome models.
      x `outcome_mod` was fit with the "poisson" family.
      i Fit `outcome_mod` with a binomial or quasibinomial family for a binary outcome, or a gaussian identity link (or an `lm()`) for a continuous outcome.

# the propensity score link error names the link and the supported set

    Code
      ipw(mods$ps_mod, mods$outcome_mod)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support the "cauchit" link for a binary propensity score model.
      i Supported links: "logit", "probit", and "cloglog".
      i Refit `wt_mod` with a supported link, or set `ps_link` to one of them.

# the missing-exposure error names the exposure and directs the user to refit

    Code
      ipw(mods$ps_mod, no_exposure, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must contain the exposure "z".
      x "z" does not appear in the formula of `outcome_mod`.
      i Refit `outcome_mod` with "z" on the right-hand side of the formula.

# the degenerate-design error names the pinned level and the remedy

    Code
      ipw(mods$ps_mod, no_intercept)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must be able to represent the outcome at every exposure level.
      x Setting `z` to "0" leaves the counterfactual design identically zero, which pins the marginal mean there to the outcome link's zero point instead of estimating it.
      i Include an intercept in `outcome_mod`, or code the exposure as a factor, whose no-intercept coding is saturated and represents every level.

