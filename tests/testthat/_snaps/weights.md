# ATE works for binary cases

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_ate()`:
      ! The propensity score must be between 0 and 1.
      i The range of `.propensity` is -0.1 and 3.3

# the refusal of an absent level names the levels the exposure takes

    Code
      expr
    Condition <propensity_focal_level_error>
      Error in `wt_att()`:
      ! `.focal_level` must be a level that `.exposure` takes.
      x No observation takes the value "Treated".
      i Levels present: "control" and "treated".

# ATE errors appropriately for categorical with vector propensity scores

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error in `wt_ate()`:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

# wt_ate() with ps_trim issues refit warning if not refit, no warning if refit

    Code
      out <- wt_ate(trimmed_ps, .exposure = z, exposure_type = "binary",
        .focal_level = 1)
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# Other estimands (att, atu, etc.) with ps_trim or ps_trunc

    Code
      out <- wt_att(trimmed_ps, .exposure = z, exposure_type = "binary",
        .focal_level = 1)
    Condition <propensity_no_refit_warning>
      Warning in `wt_att()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_atu.ps_trim triggers refit check, sets 'atu; trimmed'

    Code
      out <- wt_atu(trimmed_obj, .exposure = z, exposure_type = "binary",
        .focal_level = 1)
    Condition <propensity_no_refit_warning>
      Warning in `wt_atu()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_atm.ps_trim triggers refit check, sets 'atm; trimmed'

    Code
      out <- wt_atm(trimmed_obj, .exposure = z, exposure_type = "binary",
        .focal_level = 1)
    Condition <propensity_no_refit_warning>
      Warning in `wt_atm()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# wt_ato.ps_trim triggers refit check, sets 'ato; trimmed'

    Code
      out <- wt_ato(trimmed_obj, .exposure = z, exposure_type = "binary",
        .focal_level = 1)
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
      i The range of `.propensity` is -0.1 and 3.3

# wt_entropy works with ps_trim objects

    Code
      out <- wt_entropy(ps_trimmed, .exposure = c(0, 0, 1, 0), exposure_type = "binary")
    Condition <propensity_no_refit_warning>
      Warning in `wt_entropy()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# entropy weights error on unsupported exposure types

    Code
      expr
    Condition <propensity_matrix_type_error>
      Error in `wt_entropy()`:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

---

    Code
      expr
    Condition <causalgenerics_unsupported_exposure_type>
      Error in `wt_entropy()`:
      ! Exposure type "continuous" is not supported.
      i Supported exposure types: "binary" and "categorical".

# wt_ate works with data frames

    Code
      expr
    Condition <propensity_df_ncol_error>
      Error in `wt_ate()`:
      ! `.propensity` data frame must have at least one column.

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error in `wt_ate()`:
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
      i The range of `.propensity` is -0.1 and 1.1

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_att()`:
      ! The propensity score must be between 0 and 1.
      i The range of `.propensity` is 0.0 and 1.0

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
      Error in `wt_ate()`:
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
      Error in `wt_att()`:
      ! For categorical exposures, `.propensity` must be a matrix or data frame.

# an invalid exposure_type on the numeric method names wt_ate()

    Code
      expr
    Condition <rlang_error>
      Error in `wt_ate()`:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "wrong".

# an invalid exposure_type on the data frame method names wt_att()

    Code
      expr
    Condition <rlang_error>
      Error in `wt_att()`:
      ! `exposure_type` must be one of "auto", "binary", or "categorical", not "wrong".

# an invalid exposure_type on the numeric method names wt_ato()

    Code
      expr
    Condition <rlang_error>
      Error in `wt_ato()`:
      ! `exposure_type` must be one of "auto", "binary", or "categorical", not "wrong".

# an invalid exposure_type on the glm method names wt_ate()

    Code
      expr
    Condition <rlang_error>
      Error in `wt_ate()`:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "wrong".

# an invalid exposure_type on the wt_cens route names wt_cens()

    Code
      expr
    Condition <rlang_error>
      Error in `wt_cens()`:
      ! `exposure_type` must be one of "auto", "binary", or "continuous", not "wrong".

# an invalid exposure_type on a trimmed propensity score names wt_ate()

    Code
      expr
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.
    Condition <rlang_error>
      Error in `wt_ate()`:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "wrong".

# a one-dimensional exposure counts its dimension in the singular

    Code
      expr
    Condition <propensity_binary_transform_error>
      Error in `wt_ate()`:
      ! `.exposure` must be a vector.
      x It is <array> with 1 dimension.
      i Supply one value per observation.

# data frame methods error appropriately

    Code
      expr
    Condition <propensity_df_ncol_error>
      Error in `wt_ate()`:
      ! `.propensity` data frame must have at least one column.

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error in `wt_ate()`:
      ! Column selection failed:

---

    Code
      expr
    Condition <propensity_df_column_error>
      Error in `wt_ate()`:
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
      Error in `wt_ate()`:
      ! The propensity score must be between 0 and 1.
      i The range of `.propensity` is 0.5 and 1.5

# GLM methods error appropriately

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! A binary propensity score needs a model of the probability of the exposure.
      x `.propensity` is <lm>, whose fitted values are conditional means rather than probabilities.
      i Fit the propensity score model with `binomial()`, or pass fitted probabilities to `.propensity` directly.

---

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate()`:
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
    Condition <causalgenerics_unsupported_exposure_type>
      Error in `wt_att()`:
      ! Exposure type "continuous" is not supported.
      i Supported exposure types: "binary" and "categorical".

# all methods handle NAs appropriately

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate()`:
      ! `.propensity` and `.exposure` must have the same length.
      i `.propensity` has length 18
      i `.exposure` has length 20

