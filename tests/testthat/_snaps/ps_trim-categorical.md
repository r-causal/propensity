# ps_trim preserves all treatment groups

    Code
      expr
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_trim validates delta < 1/k

    Code
      expr
    Condition <propensity_range_warning>
      Warning in `ps_trim()`:
      Invalid trimming threshold (delta >= 1/k); returning original data

# ps_trim errors for unsupported methods with categorical

    Code
      expr
    Condition <rlang_error>
      Error in `ps_trim()`:
      ! `method` must be one of "ps" or "optimal", not "adaptive".

---

    Code
      expr
    Condition <rlang_error>
      Error in `ps_trim()`:
      ! `method` must be one of "ps" or "optimal", not "pctl".

# ps_trim requires exposure for categorical

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! `.exposure` must be provided for categorical propensity score trimming.

# is_unit_trimmed works for matrix objects

    Code
      expr
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_trim warns when no column names provided

    Code
      expr
    Condition <propensity_matrix_no_names_warning>
      Warning:
      Propensity score matrix has no column names.
      i Assuming columns are in factor level order: "A", "B", and "C"
      i This may lead to incorrect results if columns are misaligned.

# ps_trim.ps_trim warns about already trimmed scores

    Code
      expr
    Condition <propensity_already_modified_warning>
      Warning in `ps_trim()`:
      Propensity scores have already been trimmed. Returning original object.

# ps_trim handles edge cases consistently with PSweight

    Code
      expr
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# ps_refit errors when all observations are trimmed for categorical

    Code
      expr
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

---

    Code
      expr
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
      expr
    Condition <propensity_no_data_warning>
      Warning in `ps_trim()`:
      One or more groups removed after trimming; returning original data

# wt_ate warns when using trimmed but not refitted categorical PS

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.
    Message <cliMessage>
      i Treating `.exposure` as categorical

