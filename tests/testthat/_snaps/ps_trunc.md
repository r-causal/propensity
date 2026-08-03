# ps_trunc() errors on invalid usage or .exposure

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trunc()`:
      ! For `method = 'cr'`, must supply `.exposure`.

---

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_trunc()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.focal_level` and `.reference_level`.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! `lower` must be smaller than `upper`
      x `lower` is 0.8 and `upper` is 0.3

# ps_trunc warns when combining objects with different parameters

    Code
      out <- c(ps_trunc1, ps_trunc2)
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trunc.ps_trunc()`:
      Converting ps_trunc to numeric: different truncation parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# ps_trunc names `.exposure` when the common range method requires one

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trunc()`:
      ! For `method = 'cr'`, must supply `.exposure`.

# ps_trunc() refuses an exposure with missing values

    Code
      expr
    Condition <propensity_missing_value_error>
      Error in `ps_trunc()`:
      ! `.exposure` must not have missing values for method "cr".
      x 1 exposure value is missing.
      i The cutoffs are read off the exposure groups, and a unit that belongs to neither leaves them undefined.
      i Remove or impute the missing exposure values.

# ps_trunc() requires lower below upper for the pctl method

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! `lower` must be smaller than `upper`
      x `lower` is 0.95 and `upper` is 0.05

# ps_trunc() refuses percentile bounds outside the unit interval

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! For `method = 'pctl'`, `lower` and `upper` must be between 0 and 1.
      x `lower` is -0.1 and `upper` is 0.95.
      i The bounds are quantile probabilities rather than propensity scores, so each must lie in [0, 1].

# ps_trunc() refuses a bound that is missing

    Code
      expr
    Condition <propensity_missing_value_error>
      Error in `ps_trunc()`:
      ! `lower` must not be missing.
      i Each bound is read into the comparison that decides what happens to a propensity score, and a missing bound decides nothing.
      i Supply a value, or leave the argument unset to take the default for this method.

# ps_trunc() refuses a common range the exposure groups do not share

    Code
      expr
    Condition <propensity_no_overlap_error>
      Error in `ps_trunc()`:
      ! The exposure groups' propensity score distributions do not overlap.
      x The lowest score among the focal units is 0.7, above the highest among the reference units, 0.6.
      i `method = 'cr'` bounds the scores to the region both groups reach, and these groups reach none in common.
      i Refit the propensity score model, or truncate with `method = 'ps'` or `method = 'pctl'`.

# combining a ps_trunc with an integer keeps the propensity scores

    Code
      out <- vec_c(x, 1L)
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.ps_trunc.integer()`:
      Converting ps_trunc to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# the column ps_trunc() announces is named in a full sentence

    Code
      truncated <- ps_trunc(fixture$ps, method = "ps", lower = 0.3, upper = 0.7)
    Message
      i Using the ".pred_1" column as the propensity score for `ps_trunc()`.

