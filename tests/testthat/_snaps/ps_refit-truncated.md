# ps_refit() refuses a truncated score vector and says why

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_refit()`:
      ! `ps_refit()` cannot refit truncated propensity scores.
      x Truncation bounds the scores and keeps every unit, so a refit on the kept units is the original model, and it would return the scores before they were bounded.
      i To refit the model, refit it before truncating, then call `ps_trunc()` on the refit model's scores.

# ps_refit() refuses a truncated score matrix and says why

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_refit()`:
      ! `ps_refit()` cannot refit truncated propensity scores.
      x Truncation bounds the scores and keeps every unit, so a refit on the kept units is the original model, and it would return the scores before they were bounded.
      i To refit the model, refit it before truncating, then call `ps_trunc()` on the refit model's scores.

