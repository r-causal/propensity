# a basis fit announces the reading it reports

    Code
      res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
    Message
      i `ipw()` reports only the conditional reading when the exposure enters `outcome_mod` through several columns.
      i The coefficient surface is the outcome model's own, and no single row of it is a causal effect.
      i Marginalizing over the dose is left to the marginaleffects package: call `avg_slopes()` or `avg_comparisons()` on the conditional result. See <https://marginaleffects.com/chapters/interactions.html>.
      i Set `effects = "conditional"` to silence this message.

# a basis fit refuses the marginal reading

    Code
      expr
    Condition <propensity_ipw_effects_error>
      Error in `ipw()`:
      ! `ipw()` reports only the conditional reading when the exposure enters `outcome_mod` through several columns.
      x `effects = "marginal"` asks for a reading this model has none of.
      i The coefficient surface is the outcome model's own, and no single row of it is a causal effect.
      i Marginalizing over the dose is left to the marginaleffects package: call `avg_slopes()` or `avg_comparisons()` on the conditional result. See <https://marginaleffects.com/chapters/interactions.html>.
      i Use `effects = "conditional"` or omit `effects`.

# the refusal without .data names the exposure and the remedy

    Code
      expr
    Condition <propensity_columns_exist_error>
      Error in `ipw()`:
      ! "e" not found in `model.frame(outcome_mod)`.
      x The frame holds `poly(e, 2)` instead, which is a term built from it.
      i The outcome model may have transformations in the formula, and a model frame records the term rather than the variables inside it.
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

