# ps_trim preserves all treatment groups

    Code
      out <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_trim validates delta < 1/k

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trim()`:
      ! The trimming threshold must fall below 1/k, for k columns of propensity scores.
      x `lower` is 0.35, and 1/k is 0.3333333 for the 3 columns the scores hold.
      i No row summing to one can hold every score above 1/k, so a threshold there leaves no rule to apply.

# ps_trim errors for unsupported methods with categorical

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "adaptive" is not supported for categorical exposures.
      i Use "ps" or "optimal".

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "pctl" is not supported for categorical exposures.
      i Use "ps" or "optimal".

# ps_trim requires exposure for categorical

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! `.exposure` must be provided for categorical propensity score trimming.

# is_unit_trimmed works for matrix objects

    Code
      out <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.25)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_trim warns when no column names provided

    Code
      out <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
    Condition <propensity_matrix_no_names_warning>
      Warning in `ps_trim()`:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.

# ps_trim.ps_trim warns about already trimmed scores

    Code
      out <- ps_trim(trimmed_once, .exposure = exposure, method = "ps", lower = 0.2)
    Condition <propensity_already_modified_warning>
      Warning in `ps_trim()`:
      Propensity scores have already been trimmed. Returning original object.

# ps_trim handles edge cases consistently with PSweight

    Code
      out <- ps_trim(.propensity = ps_matrix, .exposure = trt, method = "ps", lower = 0.06)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_refit errors when all observations are trimmed for categorical

    Code
      out <- ps_trim(.propensity = ps_matrix, .exposure = exposure, method = "ps",
        lower = 0.3)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

---

    Code
      out <- ps_trim(.propensity = ps_matrix_extreme, .exposure = exposure, method = "ps",
        lower = 0.02)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

---

    Code
      expr
    Condition <propensity_no_data_error>
      Error in `ps_refit()`:
      ! No retained rows to refit on (all were trimmed).

# ps_refit handles minimal data for categorical exposures

    Code
      out <- ps_trim(.propensity = ps_matrix, .exposure = test_data$trt, method = "ps",
      lower = 0.25)
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# wt_ate warns when using trimmed but not refitted categorical PS

    Code
      out <- wt_ate(trimmed_ps, .exposure = exposure)
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

