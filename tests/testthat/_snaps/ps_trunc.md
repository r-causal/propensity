# ps_trunc() errors on invalid usage or .exposure

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_trunc()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.focal_level` and `.reference_level`.

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
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trunc.ps_trunc()`:
      Converting ps_trunc to numeric: different truncation parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

