# ipw() print output is stable per SE method

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Estimates:
              estimate  std.err      z ci.lower ci.upper conf.level   p.value    
      rd      0.243086 0.048326 5.0302  0.14837   0.3378       0.95 4.900e-07 ***
      log(rr) 0.597322 0.126778 4.7116  0.34884   0.8458       0.95 2.458e-06 ***
      log(or) 1.021974 0.211405 4.8342  0.60763   1.4363       0.95 1.337e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Estimates:
              estimate  std.err      z ci.lower ci.upper conf.level   p.value    
      rd      0.243086 0.048386 5.0239  0.14825  0.33792       0.95 5.064e-07 ***
      log(rr) 0.597322 0.126936 4.7057  0.34853  0.84611       0.95 2.530e-06 ***
      log(or) 1.021974 0.211670 4.8281  0.60711  1.43684       0.95 1.378e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the no-intercept rejection names the intercept and the SE method

    Code
      ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `ipw()` supports "linearization" standard errors only for an outcome model with an intercept.
      x `outcome_mod` was fit without an intercept.
      i This covers every no-intercept coding, including a saturated factor coding of the exposure, which the M-estimation path accepts.
      i Include an intercept in `outcome_mod`, or use `se_method = "mestimation"`.

# the linearization atu rejection names the SE method

    Code
      ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support "linearization" standard errors for the "atu" estimand.
      i Use `se_method = "mestimation"` for the "atu" estimand.

# the matrix-response propensity model error names the matrix response

    Code
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support a matrix response in the propensity score model.
      x `wt_mod` has a matrix response, such as `cbind(successes, failures)`; the exposure must be a single-column response.
      i Fit `wt_mod` with the exposure itself as the response, adding it to the data as its own column first if it has to be computed.

# the matrix-response outcome model error names the outcome model

    Code
      ipw(m$ps_mod, m$outcome_mod)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support a matrix response in the outcome model.
      x `outcome_mod` has a matrix response, such as `cbind(successes, failures)`; the marginal means are estimated from a single-column response.
      i Refit `outcome_mod` on one row per observation, weighted by the propensity score weights. Aggregated counts have to be expanded first, because the outcome model must carry exactly those weights.

# the ps_link mismatch error names both links

    Code
      ipw(mods$ps_mod, mods$outcome_mod, ps_link = "probit", se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `ps_link` must match the link `wt_mod` was fit with for "linearization" standard errors.
      x `ps_link` is "probit"; `wt_mod` was fit with a "logit" link.
      i Omit `ps_link` to use the model's own link, or refit `wt_mod` with the link you intend.

# the focal level rejection names both the recorded and coded levels

    Code
      ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! The weights that fit `outcome_mod` target a different exposure level than `ipw()` does.
      x The weights record "control" as the focal level, but `ipw()` treats "treat" as focal: it takes the second sorted level of a binary exposure as the exposed group.
      i Relevel `zf` so that "control" sorts second, then refit both `wt_mod` and `outcome_mod`.

# the linearization separation error names the count and the model

    Code
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `wt_mod` must not separate the exposure.
      x Putting the fitted linear predictors through the link's inverse gives a probability of exactly 0 or 1 for 389 observations, whose weights are then undefined.
      i This is usually separation: some covariate pattern predicts the exposure without error, so the fit has no finite maximum likelihood estimate. An extreme covariate pattern can saturate the scores even where the estimate is finite.
      i Check overlap in `wt_mod` rather than the weights. Dropping or combining the covariate that separates, or penalizing the fit, gives a model with finite coefficients.

# the call-form propensity response error reads in the user's terms

    Code
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support a transformed response in the propensity score model.
      x `wt_mod` reads the exposure through `factor(z)`, an expression rather than a single column.
      i Fit `wt_mod` with the exposure itself as the response, adding it to the data as its own column first if it has to be computed.

# the intercept-only outcome rejection reads in the user's terms

    Code
      ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `ipw()` supports "linearization" standard errors only for an outcome model of the exposure alone.
      x `outcome_mod` does not include the exposure "z".
      i Use `se_method = "mestimation"` for a covariate-adjusted outcome model.

