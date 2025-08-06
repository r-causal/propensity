# c() with psw objects of different estimands warns and returns correct numeric values

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.double.psw()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr

# c() with psw and numeric values warns and returns numeric

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# c() with ps_trim objects of different parameters warns and returns numeric

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trim.ps_trim()`:
      Converting ps_trim to numeric: different trimming parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# c() with ps_trunc objects behaves correctly

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.ps_trunc.ps_trunc()`:
      Converting ps_trunc to numeric: different truncation parameters
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# c() with mixed propensity classes warns and returns numeric

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.ps_trim()`:
      Converting psw and ps_trim to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.ps_trunc()`:
      Converting psw and ps_trunc to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.ps_trim.ps_trunc()`:
      Converting ps_trim and ps_trunc to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
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
      expr

# c() with empty vectors works correctly

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# c() with single values works correctly

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'ato'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# c() with different stabilization status warns

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: different stabilization status
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# subsetting operations work correctly

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# append() works like c()

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# data.frame operations work as expected

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# vctrs vec_ptype2 returns appropriate prototypes

    Code
      expr
    Condition <propensity_coercion_warning>
      Warning in `vec_ptype2.psw.psw()`:
      Converting psw to numeric: incompatible estimands 'ate' and 'att'
      i Metadata cannot be preserved when combining incompatible objects
      i Use identical objects or explicitly cast to numeric to avoid this warning

# comparison operations warn about class downgrade

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# c() ordering matters for warnings

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

