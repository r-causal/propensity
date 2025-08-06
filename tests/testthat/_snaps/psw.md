# vec_cast works for psw and integer with precision checks

    Code
      expr
    Condition <vctrs_error_cast_lossy>
      Error in `vec_cast.integer.psw()`:
      ! Can't convert from `psw` <double> to <integer> due to loss of precision.
      * Locations: 1, 2, 3

# vec_ptype2 combines psw and other types correctly

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
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.double()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.double.psw()`:
      Converting psw to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.psw.integer()`:
      Converting psw to integer
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

---

    Code
      expr
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.integer.psw()`:
      Converting psw to integer
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# vec_arith performs arithmetic on psw objects

    Code
      expr
    Condition <vctrs_error_incompatible_op>
      Error in `vec_arith()`:
      ! <psw{estimand = ate}> * <thing> is not permitted

---

    Code
      expr
    Condition <vctrs_error_incompatible_op>
      Error in `vec_arith()`:
      ! <psw{estimand = ate}> * <list> is not permitted

