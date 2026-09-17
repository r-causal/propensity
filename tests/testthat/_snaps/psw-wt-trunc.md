# is_unit_wt_truncated() refuses flagged weights with no record

    Code
      expr
    Condition <propensity_missing_meta_error>
      Error in `is_unit_wt_truncated()`:
      ! `is_unit_wt_truncated()` has no usable weight truncation record for these weights.
      x These weights are marked as truncated but carry no record of which weights were truncated.
      i Call `is_unit_wt_truncated()` on the weights the truncation returned, before subsetting or combining them.

# growing a weight-truncated psw leaves a record that no longer covers it

    Code
      expr
    Condition <propensity_missing_meta_error>
      Error in `is_unit_wt_truncated()`:
      ! `is_unit_wt_truncated()` has no usable weight truncation record for these weights.
      x The record covers 5 observations and these weights have 7, so its positions do not describe them.
      i Call `is_unit_wt_truncated()` on the weights the truncation returned, before subsetting or combining them.

# is_unit_wt_truncated() refuses anything that is not a psw

    Code
      expr
    Condition <propensity_method_error>
      Error in `is_unit_wt_truncated()`:
      ! `is_unit_wt_truncated()` not supported for class "numeric"

