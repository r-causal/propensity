# ps_trunc validates delta < 1/k

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! Invalid truncation threshold (delta >= 1/k).

# ps_trunc errors for unsupported methods with categorical

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trunc()`:
      ! Method "cr" is not supported for categorical exposures.
      i Use "ps" or "pctl".

# ps_trunc requires exposure for categorical

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trunc()`:
      ! `.exposure` must be provided for categorical propensity score truncation.

# ps_trunc warns when no column names provided

    Code
      out <- ps_trunc(ps_matrix, .exposure = exposure, method = "ps", lower = 0.05)
    Condition <propensity_matrix_no_names_warning>
      Warning in `ps_trunc()`:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.

# ps_trunc.ps_trunc warns about already truncated scores

    Code
      out <- ps_trunc(truncated_once, .exposure = exposure, method = "ps", lower = 0.1)
    Condition <propensity_already_modified_warning>
      Warning in `ps_trunc()`:
      Propensity scores have already been truncated. Returning original object.

