# the frame refusals name the column to set and the type to predict

    Code
      expr
    Condition <propensity_df_ambiguous_column_error>
      Error in `wt_ate()`:
      ! Can't tell which column of `.propensity` holds the conditional mean of the exposure.
      x `.propensity` has 2 columns, and a continuous exposure has no levels to read them against.
      i Set `.propensity_col` to the column of conditional means.

---

    Code
      expr
    Condition <propensity_df_class_column_error>
      Error in `wt_ate()`:
      ! `.propensity` holds predicted classes rather than propensity scores.
      x The ".pred_class" column holds the level predicted for each unit.
      i Predict probabilities instead, with `predict(fit, new_data, type = "prob")`, and pass the columns that returns.

---

    Code
      expr
    Condition <propensity_df_class_column_error>
      Error in `ps_trim()`:
      ! `.propensity` holds predicted classes rather than propensity scores.
      x The ".pred_class" column holds the level predicted for each unit.
      i Predict probabilities instead, with `predict(fit, new_data, type = "prob")`, and pass the columns that returns.

---

    Code
      expr
    Condition <propensity_df_class_column_error>
      Error in `ps_calibrate()`:
      ! `.propensity` holds predicted classes rather than propensity scores.
      x The ".pred_class" column holds the level predicted for each unit.
      i Predict probabilities instead, with `predict(fit, new_data, type = "prob")`, and pass the columns that returns.

