# tidyr::pivot_longer works with propensity classes

    Code
      out <- tidyr::pivot_longer(df, cols = c(ate_wts, att_wts), names_to = "weight_type",
      values_to = "weight")
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# tidyr::pivot_longer works with mixed propensity classes

    Code
      out <- tidyr::pivot_longer(df, cols = c(psw_col, trim_col, trunc_col),
      names_to = "type", values_to = "value")
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.ps_trim()`:
      Converting psw and ps_trim to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning
      Warning in `vec_ptype2.double.ps_trunc()`:
      Converting ps_trunc to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      out <- expect_propensity_warning(tidyr::pivot_longer(df, cols = c(psw_col,
        trim_col, trunc_col), names_to = "type", values_to = "value"))

# c() works as expected with warnings

    Code
      out <- c(x, y)
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

---

    Code
      out <- c(x, z)
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# rbind and data frame operations work

    Code
      out <- vctrs::vec_rbind(df1, df2)
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# tidyr works with stabilized weights

    Code
      out <- tidyr::pivot_longer(df, cols = c(wt_stab, wt_unstab), names_to = "stab_type",
      values_to = "weight")
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: different stabilization status
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

