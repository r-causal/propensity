# ps_refit() refuses a model of level probabilities for a density trim

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_refit()`:
      ! A trimmed dose model can only be refit with a model of the exposure's conditional mean.
      x The trimming record holds a dose model's conditional means, and `model` fits a probability for each of its levels.
      i Refit with the model of the continuous exposure the trimming was made from.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_refit()`:
      ! A trimmed dose model can only be refit with a model of the exposure's conditional mean.
      x The trimming record holds a dose model's conditional means, and `model` fits a probability for each of its levels.
      i Refit with the model of the continuous exposure the trimming was made from.

# an unrefit density trim warns that it was not refit

    Code
      out <- wt_ate(fx$trimmed, fx$dat$a, exposure_type = "continuous")
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# weights from a density trim are refused for an exposure that is not a dose

    Code
      expr
    Condition <propensity_modified_continuous_error>
      Error in `wt_ate()`:
      ! Weights for a binary exposure cannot be built from a density-trimmed dose model.
      x The trimming record holds conditional means and a density family, not propensity scores.
      i Pass the dose as `.exposure`, or trim a propensity score model for the binary exposure instead.

---

    Code
      expr
    Condition <propensity_modified_continuous_error>
      Error in `wt_ate()`:
      ! Weights for a binary exposure cannot be built from a density-trimmed dose model.
      x The trimming record holds conditional means and a density family, not propensity scores.
      i Pass the dose as `.exposure`, or trim a propensity score model for the binary exposure instead.

