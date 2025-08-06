# ps_trim subsetting with [ preserves class and updates indices

    Code
      expr
    Condition <propensity_length_error>
      Error in `x[c(TRUE, FALSE)]`:
      ! Logical subscript `i` must be size 1 or 10, not 2.

# ps_trim rejects infinite values

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trim()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is 0.1 and Inf

