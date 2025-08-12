# ATE works for binary cases

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_ate()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is -0.1 and 3.3

# ATE errors appropriately for categorical with vector propensity scores

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

# wt_ate() with ps_trim issues refit warning if not refit, no warning if refit

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# Other estimands (att, atu, etc.) with ps_trim or ps_trunc

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_att()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_atu.ps_trim triggers refit check, sets 'atu; trimmed'

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_atu()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_atm.ps_trim triggers refit check, sets 'atm; trimmed'

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_atm()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_ato.ps_trim triggers refit check, sets 'ato; trimmed'

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_ato()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_entropy works for binary cases

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_entropy()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is -0.1 and 3.3

# wt_entropy works with ps_trim objects

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_entropy()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# entropy weights error on unsupported exposure types

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

---

    Code
      expr
    Condition <propensity_wt_not_supported_error>
      Error in `wt_entropy()`:
      ! Exposure type "continuous" not currently supported for entropy

# wt_ate works with data frames

    Code
      expr
    Condition <propensity_df_ncol_error>
      Error:
      ! `.propensity` data frame must have at least one column.

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error:
      ! Column selection failed:

# GLM methods error on non-GLM objects

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_ate()`:
      ! No method for objects of class character

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_att()`:
      ! No method for objects of class list

# wt_* functions error appropriately on invalid inputs

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_ate()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is -0.1 and 1.1

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_att()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is 0.0 and 1.0

---

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate()`:
      ! `.propensity` and `.exposure` must have the same length.
      i `.propensity` has length 2
      i `.exposure` has length 3

---

    Code
      expr
    Condition <rlang_error>
      Error in `match_exposure_type()`:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "invalid".

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_ate()`:
      ! No method for objects of class character

---

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

# data frame methods error appropriately

    Code
      expr
    Condition <propensity_df_ncol_error>
      Error:
      ! `.propensity` data frame must have at least one column.

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error:
      ! Column selection failed:

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error:
      ! Column selection failed:

---

    Code
      expr
    Condition <simpleWarning>
      Warning in `check_ps_range()`:
      NAs introduced by coercion
    Condition <simpleError>
      Error in `.exposure / .propensity`:
      ! non-numeric argument to binary operator

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `weight_fn_numeric()`:
      ! The propensity score must be between 0 and 1.
      i The range of `ps` is 0.5 and 1.5

# GLM methods error appropriately

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_ate()`:
      ! No method for objects of class lm

---

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate.numeric()`:
      ! `.propensity` and `.exposure` must have the same length.
      i `.propensity` has length 2
      i `.exposure` has length 4

# default methods provide informative errors

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_ate()`:
      ! No method for objects of class my_custom_class

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_att()`:
      ! No method for objects of class my_custom_class

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_atu()`:
      ! No method for objects of class my_custom_class

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_atm()`:
      ! No method for objects of class my_custom_class

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_ato()`:
      ! No method for objects of class my_custom_class

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `wt_entropy()`:
      ! No method for objects of class my_custom_class

# GLM methods handle non-binomial families appropriately

    Code
      expr
    Condition <rlang_error>
      Error in `match_exposure_type()`:
      ! `exposure_type` must be one of "auto", "binary", or "categorical", not "continuous".

# all methods handle NAs appropriately

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate.numeric()`:
      ! `.propensity` and `.exposure` must have the same length.
      i `.propensity` has length 18
      i `.exposure` has length 20

