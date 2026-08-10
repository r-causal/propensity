# the refusal without .data names the exposure and the remedy

    Code
      expr
    Condition <propensity_columns_exist_error>
      Error in `ipw()`:
      ! "e" not found in `model.frame(outcome_mod)`.
      i The outcome model may have transformations in the formula.
      i Please specify `.data`

# the mixed-term refusal names the term and the way out

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a marginal structural model whose exposure terms read the exposure alone.
      x `poly(e, 2):x1` is not a term of "e" alone.
      i A term reading anything else contributes a coefficient that depends on what it reads, so there is no one effect for a row to report.
      i A term reading the exposure alone is admitted however it is written, so a curve such as `e + I(e^2)` reports one row per coefficient.
      i Read the full coefficient vector from the returned fit object for a model this surface cannot report.

