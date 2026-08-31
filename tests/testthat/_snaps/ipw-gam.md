# the registry refuses an additive fit whose penalty it cannot place

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` reads the penalty of a <gam/glm/lm> propensity score model off the smooth terms it was fit with.
      x `wt_mod` records more smoothing parameters than its smooth terms account for, which is what a penalty on a parametric term, such as one from `paraPen`, adds.
      i Refit `wt_mod` so that every penalty in it belongs to a smooth term.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

