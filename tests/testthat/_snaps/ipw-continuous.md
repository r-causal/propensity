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
      i Use `se_method = "bootstrap"`, passing the data the models were fit to in `.data`, which resamples the whole fit instead of stacking it.

# the MM refusal names the method and the psi it finishes on

    Code
      expr
    Condition <propensity_ipw_robust_psi_error>
      Error in `ipw()`:
      ! `ipw()` stacks only the Huber score of a <rlm/lm> propensity score model of a continuous exposure.
      x `wt_mod` was fit with `method = "MM"`, which starts from a high-breakdown fit and finishes on `MASS::psi.bisquare()`.
      i Refit `wt_mod` with `MASS::psi.huber()`, the default, whose threshold `ipw()` reads off the fit.
      i Use `se_method = "bootstrap"`, passing the data the models were fit to in `.data`, which resamples the whole fit instead of stacking it.

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
      ! `ipw()` supports only `stats::lm()`, gaussian `stats::glm()`, or `MASS::rlm()` propensity score models for a continuous exposure.
      x `wt_mod` has class <mymodel/lm>.
      i A <gam> is recognized and refused on its own terms; every other class reaches this refusal.
      i Refit `wt_mod` with `stats::lm()` or `stats::glm(family = gaussian())`.

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
      i Use `se_method = "bootstrap"`, passing the data the models were fit to in `.data`, which resamples the whole fit instead of stacking it.

# the kernel-density refusal names the bandwidth as the reason

    Code
      expr
    Condition <propensity_ipw_se_method_unavailable_error>
      Error in `ipw()`:
      ! `ipw()` cannot build a sandwich variance for weights built with a "kernel" density.
      x The bandwidth of a kernel estimate is chosen from the residuals of the propensity score model, so the weights are not a differentiable function of that model's parameters.
      i Use `se_method = "bootstrap"`, passing the data the models were fit to in `.data`, which resamples the whole fit instead of stacking it.

# ipw() continuous bootstrap refusals read well

    Code
      expr
    Condition <propensity_ipw_data_error>
      Error in `ipw()`:
      ! `ipw()` requires `.data` for "bootstrap" standard errors.
      x Each replicate refits both models on a resample of the rows, and there are no rows to resample without `.data`.
      i Supply the data the models were fit to in `.data`.

---

    Code
      expr
    Condition <propensity_unsupported_arg_error>
      Error in `ipw()`:
      ! `boot_reps` is not supported with "mestimation" standard errors.
      x `boot_reps` describes the resampling, and "mestimation" resamples nothing.
      i Use `se_method = "bootstrap"`, or drop `boot_reps`.

---

    Code
      expr
    Condition <propensity_ipw_bootstrap_error>
      Error in `ipw()`:
      ! `ipw()` needs at least 50 successful bootstrap replicates.
      x 40 of 40 replicates succeeded, and 0 failed.
      i Raise `boot_reps`, or refit the models on data every resample of them can be fit to.

# ipw() refuses the bootstrap for a binary exposure in so many words

    Code
      expr
    Condition <propensity_method_error>
      Error in `ipw()`:
      ! `ipw()` supports "bootstrap" standard errors only for a continuous exposure.
      x `wt_mod` is a propensity score model of a binary exposure.
      i Use `se_method = "mestimation"`, which builds a sandwich variance for every fit this exposure type accepts.

# ipw() says how many bootstrap replicates it dropped

    Code
      out <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat, se_method = "bootstrap",
      boot_reps = 100L, boot_seed = 1L)
    Condition <propensity_ipw_bootstrap_warning>
      Warning in `ipw()`:
      20 of 100 bootstrap replicates failed and were dropped.
      x A replicate fails when a refit errors or leaves a coefficient that is not finite, so the spread is the spread of the resamples the models could be fit to.
      i Check whether either model is fit to something a resample of the rows can break, such as a factor covariate with a level so rare that a resample can leave it out.

