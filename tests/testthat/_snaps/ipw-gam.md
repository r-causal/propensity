# the registry refuses an additive fit whose penalty it cannot place

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` reads the penalty of a <gam/glm/lm> propensity score model off the smooth terms it was fit with.
      x `wt_mod` records more smoothing parameters than its smooth terms account for, which is what a penalty on a parametric term, such as one from `paraPen`, adds.
      i Refit `wt_mod` so that every penalty in it belongs to a smooth term.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# the registry refuses a gam holding a smoothing floor

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` stacks a <gam/glm/lm> propensity score model at the penalized score its coefficients solve, which it rebuilds from the penalty the fit records.
      x The penalty `wt_mod` records does not reproduce the score `wt_mod` is at, so a system seeded at its coefficients would settle away from them.
      i A smoothing floor from `min.sp` is one cause: it is added to the penalty and left out of the smoothing parameters the fit reports. Refit `wt_mod` without `min.sp`.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

