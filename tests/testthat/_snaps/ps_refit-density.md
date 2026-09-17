# ps_refit() refuses a model of level probabilities for a density trim

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_refit()`:
      ! A trimmed dose model can only be refit with a model of the exposure's conditional mean.
      x The trimming record holds a dose model's conditional means, and `model` is <multinom>, a model of the probabilities of a discrete exposure's levels.
      i Refit with the model of the continuous exposure the trimming was made from, such as one fit with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_refit()`:
      ! A trimmed dose model can only be refit with a model of the exposure's conditional mean.
      x The trimming record holds a dose model's conditional means, and `model` is <multinom>, a model of the probabilities of a discrete exposure's levels.
      i Refit with the model of the continuous exposure the trimming was made from, such as one fit with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`.

# ps_refit() refuses a model of a probability for a density trim and names every dose model

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_refit()`:
      ! A trimmed dose model can only be refit with a model of the exposure's conditional mean.
      x `model` was fit with `binomial()`, whose spread changes with its fitted values.
      i The trimming record holds a dose model's conditional means, so refit with the model of the continuous exposure the trimming was made from, such as one fit with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`.

# an unrefit density trim warns that it was not refit

    Code
      out <- wt_ate(fx$trimmed, fx$dat$a, exposure_type = "continuous")
    Condition <propensity_no_refit_warning>
      Warning in `wt_ate()`:
      It appears you trimmed your propensity score but did not refit the model.
      i Use `ps_refit()` for more accurate re-estimation.

# weights from a density trim refuse a density the record does not hold

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `.density` must be the density the dose model was trimmed under.
      x The trimming record holds `dens_normal()`, and `.density` is `dens_laplace(sigma_method = "mle")`.
      i The units were kept by how plausible their dose is under the recorded density, so the weights are read under it too. Drop `.density`, or trim the dose model again under the density you want.

# weights from a density trim refuse a spread of the caller's

    Code
      expr
    Condition <propensity_sigma_error>
      Error in `wt_ate()`:
      ! `.sigma` cannot be used with a density-trimmed dose model.
      x The trimming record holds the spread the conditional density was read at, and `.sigma` would be a second spread for the same density.
      i Drop `.sigma`. To read the density at a spread of your own, pass it to `ps_trim()` when trimming, and `ps_refit()` keeps it.

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

