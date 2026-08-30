# a model of the wrong kind is refused by what it was fit with

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
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! A binary propensity score needs a model of the probability of the exposure.
      x `.propensity` was fit with `gaussian()`, whose fitted values are conditional means rather than probabilities.
      i Fit the propensity score model with `binomial()`, or pass fitted probabilities to `.propensity` directly.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure need a model of its conditional mean with a single spread.
      x `.propensity` was fit with `poisson()`, whose spread changes with its fitted values.
      i Fit the propensity score model with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`, or pass fitted conditional means to `.propensity` directly.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure need a model of its conditional mean with a single spread.
      x `.propensity` was fit with `quasi(variance = "mu")`, whose spread changes with its fitted values.
      i Fit the propensity score model with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`, or pass fitted conditional means to `.propensity` directly.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error:
      ! A binary propensity score needs a model of the probability of the exposure.
      x `.propensity` was fit with an unnamed family, whose fitted values are conditional means rather than probabilities.
      i Fit the propensity score model with `binomial()`, or pass fitted probabilities to `.propensity` directly.

