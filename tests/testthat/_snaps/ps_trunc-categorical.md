# ps_trunc validates delta < 1/k

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! Invalid truncation threshold (delta >= 1/k).

# ps_trunc errors for unsupported methods with categorical

    Code
      expr
    Condition <rlang_error>
      Error in `ps_trunc()`:
      ! `method` must be one of "ps" or "pctl", not "cr".

# ps_trunc requires exposure for categorical

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trunc()`:
      ! `.exposure` must be provided for categorical propensity score truncation.

# ps_trunc warns when no column names provided

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.

# ps_trunc.ps_trunc warns about already truncated scores

    Code
      expr
    Condition <propensity_already_modified_warning>
      Warning in `ps_trunc()`:
      Propensity scores have already been truncated. Returning original object.

