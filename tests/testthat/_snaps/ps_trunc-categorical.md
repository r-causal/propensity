# ps_trunc validates delta < 1/k

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! The truncation threshold must fall below 1/k, for k columns of propensity scores.
      x `lower` is 0.35, and 1/k is 0.3333333 for the 3 columns the scores hold.
      i No row summing to one can hold every score above 1/k, so a threshold there leaves no rule to apply.

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

# a truncation whose own bounds sit at an endpoint is refused

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! Truncation left 5 propensity scores at 0 or 1.
      x The bounds computed from the scores run from 0 to 0.966.
      i A bound at an endpoint pins a score to where it already sat.
      i Supply `lower` and `upper` explicitly, inside the open interval, so that every score is bounded away from 0 and 1.

# the closed interval a categorical truncation reads is described

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! All propensity scores must be between 0 and 1.
      i The bounds are inclusive here: a score of exactly 0 or 1 is one truncation bounds, but a score outside the interval is not a probability.
      i The range of values in `.propensity` is -0.1 and 1.0

