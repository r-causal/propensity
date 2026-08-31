# ipw() continuous print output is stable

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: lm(formula = A ~ x1 + x2, data = dat) 
      
      Outcome Model:
        Call: lm(formula = msm_fmla, data = dat, weights = wts) 
      
      Marginal estimates:
            estimate  std.err      z ci.lower ci.upper conf.level   p.value    
      slope 0.679126 0.037942 17.899  0.60476  0.75349       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the continuous ps_link error explains why the argument does not apply

    Code
      ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit")
    Condition
      Error in `ipw()`:
      ! `ipw()` does not accept `ps_link` for a continuous propensity score model.
      x `ps_link` names the link of a binomial glm, which a continuous propensity score model is not.
      i Omit `ps_link`; it applies only to a binomial glm propensity score model.

# the observation-level spread refusal names both spreads it takes

    Code
      expr
    Condition <propensity_ipw_sigma_error>
      Error in `ipw()`:
      ! `ipw()` does not support weights built with an observation-level `.sigma`.
      x The weights supplied to `outcome_mod` record a spread for each observation, which has no counterpart in the stacked system: it estimates one conditional spread, or holds one fixed.
      i Rebuild the weights with the pooled default, by leaving `.sigma` unset, or with a single `.sigma`, which `ipw()` takes as a known constant.

# the transformed continuous propensity response error reads in the user's terms

    Code
      ipw(ps_mod, outcome_mod, .data = dat)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support a transformed response in the propensity score model.
      x `wt_mod` reads the exposure through `log(Apos)`, an expression rather than a single column.
      i Fit `wt_mod` with the exposure itself as the response, adding it to the data as its own column first if it has to be computed.

# the robust psi refusal names the psi it found and the remedy

    Code
      expr
    Condition <propensity_ipw_robust_psi_error>
      Error in `ipw()`:
      ! `ipw()` stacks only the Huber score of a <rlm/lm> propensity score model of a continuous exposure.
      x `wt_mod` was fit with `MASS::psi.bisquare()`.
      i Refit `wt_mod` with `MASS::psi.huber()`, the default, whose threshold `ipw()` reads off the fit.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# the MM refusal names the method and the psi it finishes on

    Code
      expr
    Condition <propensity_ipw_robust_psi_error>
      Error in `ipw()`:
      ! `ipw()` stacks only the Huber score of a <rlm/lm> propensity score model of a continuous exposure.
      x `wt_mod` was fit with `method = "MM"`, which starts from a high-breakdown fit and finishes on `MASS::psi.bisquare()`.
      i Refit `wt_mod` with `MASS::psi.huber()`, the default, whose threshold `ipw()` reads off the fit.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# the robust convergence refusal names the arguments that fix it

    Code
      expr
    Condition <propensity_ipw_convergence_error>
      Error in `ipw()`:
      ! `ipw()` cannot stack a <rlm/lm> propensity score model that did not converge.
      x `wt_mod` reports `converged = FALSE`, so its coefficients are not the root of the score stacked here.
      i Refit `wt_mod` with a larger `maxit`, or a looser `acc`, until it converges.

# the unknown-subclass error names the classes that are supported

    Code
      expr
    Condition <propensity_class_error>
      Error in `ipw()`:
      ! `ipw()` supports only `stats::lm()`, gaussian `stats::glm()`, or `MASS::rlm()` as the propensity score model of a continuous exposure.
      x `wt_mod` has class <mymodel/lm>.
      i A <gam> is recognized and refused on its own terms; every other class reaches this refusal.
      i Refit `wt_mod` with `stats::lm()` or `stats::glm(family = gaussian())`.

# ipw() refuses a propensity model fit without a formula

    Code
      expr
    Condition <propensity_ipw_response_error>
      Error in `ipw()`:
      ! `ipw()` needs a propensity score model fit through the formula interface.
      x `wt_mod` records no formula with a response, so the exposure cannot be read off it.
      i Refit `wt_mod` from a formula whose left-hand side is the exposure, as in `exposure ~ x`, rather than from a design matrix and a response vector.

# the continuous propensity-link error names the link and the remedy

    Code
      expr
    Condition <propensity_ipw_link_error>
      Error in `ipw()`:
      ! `ipw()` does not support the "inverse" link for the propensity score model of a continuous exposure.
      x `wt_mod` is a gaussian model fit with the "inverse" link.
      i Refit `wt_mod` with one of "identity" and "log", or as an `lm()`.

# the unavailable-method error explains what the sandwich cannot do

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` cannot build a sandwich variance for a <gam/glm/lm> propensity score model of a continuous exposure.
      x An additive model chooses how much to smooth by REML, and no estimating equation stacked here reproduces that choice, so the stacked system would describe a different fit.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# the kernel-density refusal names the bandwidth as the reason

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` cannot build a sandwich variance for weights built with a "kernel" density.
      x The bandwidth of a kernel estimate is chosen from the residuals of the propensity score model, so the weights are not a differentiable function of that model's parameters.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# ipw() refuses a numerator model whose score it cannot write

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` cannot build a sandwich variance for a <gam/glm/lm> numerator model of a continuous exposure.
      x An additive model chooses how much to smooth by REML, and no estimating equation stacked here reproduces that choice, so the stacked system would describe a different fit.
      i propensity has no resampling method; bootstrap the whole fit yourself: resample the rows, refit the propensity score model, rebuild the weights with `wt_ate()`, and refit the outcome model on each resample.

# ipw() refuses a numerator model of another response

    Code
      expr
    Condition <propensity_ipw_numerator_error>
      Error in `ipw()`:
      ! The model supplied to `stabilize` must model the exposure.
      x It models "yc" and `wt_mod` models "A".
      i The numerator of the weights is the density of the exposure given what the numerator model reads, so both models describe the same response.

# ipw() refuses a numerator model fit to other observations

    Code
      expr
    Condition <propensity_ipw_numerator_error>
      Error in `ipw()`:
      ! The model supplied to `stabilize` must be fit to the observations the other models were fit to.
      x It was fit to 200 observations and `outcome_mod` to 197.
      i Refit the numerator model on the data the other models were fit to, and rebuild the weights from it.

# ipw() refuses a numerator model with a dropped coefficient

    Code
      expr
    Condition <propensity_ipw_rank_error>
      Error in `ipw()`:
      ! `stabilize` must have a coefficient for every column of its design.
      x `stabilize` has no fitted coefficient for "x1_again".
      i A model reports that for a column its design cannot separate from the others: the column is a linear combination of them, exactly or to within the tolerance the fit pivots at, so the fit has no unique solution for it and drops it.
      i `ipw()` rebuilds the numerator of the weights by multiplying the fitted coefficients against that design, so a column with no coefficient leaves every numerator undefined.
      i Refit `stabilize` without the redundant column, or combine it with the column it duplicates.

