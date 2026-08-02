# the weight-consistency error message is stable on the binary mestimation route

    Code
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
    Condition <propensity_ipw_weights_mismatch_error>
      Error in `ipw()`:
      ! The "ate" weights recomputed from `wt_mod` differ from the weights supplied to `outcome_mod` (compared at relative tolerance 1e-6).
      i The estimand or the focal level the weights were built for may differ from the ones `ipw()` resolved.
      i A non-default `.focal_level` or `.reference_level` in the weights is one cause: `ipw()` treats the second sorted level of a binary exposure as focal.
      i Weights trimmed, truncated, or normalized after `wt_mod` was fit differ from the ones rebuilt here, which come from that model alone.
      i `.data` values that differ from the data the models were fit to move the recomputed weights on their own and leave the supplied weights exactly right.
      i Refit `outcome_mod` with weights from this propensity score model and estimand if the weights are the cause.

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

# the ps frame-gone error names the propensity model and the remedy

    Code
      ipw(ps_gone, out, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! Can't reconstruct the data behind `wt_mod`.
      x object 'dat' not found
      i Supply `.data` with the exposure, outcome, and covariates.

# the frame-gone outcome error explains why .data cannot help

    Code
      ipw(ps_mod, outcome_gone, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! Can't reconstruct the data behind `outcome_mod`.
      x object 'd_local' not found
      i Refit `outcome_mod` where its data is available, or fit it with `model = TRUE` so the model frame is kept.
      i `.data` cannot stand in here: the weights are read from `outcome_mod`'s own model frame.

# the fitted-factor .data error names the column and how the fit recorded it

    Code
      ipw(ps_mod, out, .data = dat_num, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! `.data` must supply "cov" as the factor the models were fit with.
      x `.data` has "cov" as a numeric vector.
      x `wt_mod` recorded "cov" as a factor with the levels "a", "b", and "c", and the designs rebuilt from `.data` use that coding.
      i Supply "cov" as that factor, or refit the models on the numeric column.

# the separation error names the count and points at the model

    Code
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! `wt_mod` must not separate the exposure.
      x Rebuilding the propensity scores gives a probability of exactly 0 or 1 for 389 observations, whose weights are then undefined.
      i This is usually separation: some covariate pattern predicts the exposure without error, so the fit has no finite maximum likelihood estimate. An extreme covariate pattern can saturate the scores even where the estimate is finite.
      i Check overlap in `wt_mod` rather than the weights. Dropping or combining the covariate that separates, or penalizing the fit, gives a model with finite coefficients.

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

