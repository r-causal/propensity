# errors for non-numeric ps

    Code
      expr
    Condition <propensity_type_error>
      Error in `ps_calibrate()`:
      ! `.propensity` must be a numeric vector.

# errors for out-of-range ps

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_calibrate()`:
      ! `.propensity` values must be between 0 and 1.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_calibrate()`:
      ! `.propensity` values must be between 0 and 1.

# errors when ps and .exposure have different lengths

    Code
      expr
    Condition <propensity_length_error>
      Error in `ps_calibrate()`:
      ! Propensity score vector `.propensity` must be the same length as `.exposure`.

# error handling for ambiguous treatment coding

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_calibrate()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.focal_level` and `.reference_level`.

# errors when trying to calibrate already calibrated ps

    Code
      expr
    Condition <propensity_already_calibrated_error>
      Error in `ps_calibrate()`:
      ! `.propensity` is already calibrated. Cannot calibrate already calibrated propensity scores.

# method parameter validation works

    Code
      expr
    Condition <rlang_error>
      Error in `ps_calibrate()`:
      ! `method` must be one of "logistic" or "isoreg", not "invalid".

# combining a ps_calib with an integer keeps the calibrated scores

    Code
      out <- vec_c(cal, 1L)
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.ps_calib.integer()`:
      Converting ps_calib to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# combining an integer with a ps_calib keeps the calibrated scores

    Code
      out <- vec_c(1L, cal)
    Condition <propensity_class_downgrade_warning>
      Warning in `vec_ptype2.integer.ps_calib()`:
      Converting ps_calib to numeric
      i Class-specific attributes and metadata have been dropped
      i Use explicit casting to numeric to avoid this warning

# the fallback message is recorded

    Code
      calibrated <- ps_calibrate(fixture$ps, fixture$exposure, smooth = TRUE)
    Message
      i 4 distinct propensity scores are too few to place the knots of a spline, so calibration falls back to logistic regression without one.

# ps_calibrate reports a matrix of propensity scores informatively

    Code
      expr
    Condition <propensity_type_error>
      Error in `ps_calibrate()`:
      ! `.propensity` must be a numeric vector.
      x It is <matrix> with 2 dimensions.
      i Calibration reads one propensity score per observation. Calibrate the columns of a matrix of scores one at a time.

