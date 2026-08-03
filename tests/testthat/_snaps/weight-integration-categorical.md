# weight functions work with trimmed categorical propensity scores

    Code
      out <- wt_ate(trimmed_ps, .exposure = exposure)
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# weight functions work with data.frame propensity scores for categorical

    Code
      out <- wt_ate(trimmed_ps, .exposure = exposure)
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# ATT weights work with categorical trimmed propensity scores

    Code
      out <- wt_att(trimmed_ps, .exposure = exposure, .focal_level = "Treat1")
    Condition <propensity_no_refit_warning>
      Warning in `wt_att()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

