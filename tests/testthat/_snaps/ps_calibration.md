# errors for non-numeric ps

    Code
      expr
    Condition <propensity_type_error>
      Error in `ps_calibrate()`:
      ! `ps` must be a numeric vector.

# errors for out-of-range ps

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_calibrate()`:
      ! `ps` values must be between 0 and 1.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_calibrate()`:
      ! `ps` values must be between 0 and 1.

# errors when ps and .exposure have different lengths

    Code
      expr
    Message <cliMessage>
      i Setting treatment to `1`
    Condition <propensity_length_error>
      Error in `ps_calibrate()`:
      ! Propensity score vector `ps` must be the same length as `.exposure`.

# error handling for ambiguous treatment coding

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `ps_calibrate()`:
      ! Don't know how to transform `.exposure` to 0/1 binary variable.
      i Specify `.treated` and `.untreated.`

# errors when trying to calibrate already calibrated ps

    Code
      expr
    Condition <propensity_already_calibrated_error>
      Error in `ps_calibrate()`:
      ! `ps` is already calibrated. Cannot calibrate already calibrated propensity scores.

# method parameter validation works

    Code
      expr
    Condition <rlang_error>
      Error in `ps_calibrate()`:
      ! `method` must be one of "logistic" or "isoreg", not "invalid".

