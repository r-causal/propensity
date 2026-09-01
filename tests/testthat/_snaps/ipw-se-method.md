# ipw() print output is stable per SE method

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Marginal estimates:
                     estimate  std.err       z ci.lower ci.upper conf.level   p.value    
      mean 0         0.297445 0.032228  9.2294  0.23428  0.36061       0.95 < 2.2e-16 ***
      mean 1         0.540531 0.036528 14.7976  0.46894  0.61213       0.95 < 2.2e-16 ***
      rd 1 vs 0      0.243086 0.048326  5.0302  0.14837  0.33780       0.95 4.900e-07 ***
      log(rr) 1 vs 0 0.597322 0.126778  4.7116  0.34884  0.84580       0.95 2.458e-06 ***
      log(or) 1 vs 0 1.021974 0.211405  4.8342  0.60763  1.43632       0.95 1.337e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Marginal estimates:
                     estimate  std.err       z ci.lower ci.upper conf.level   p.value    
      mean 0         0.297445 0.032268  9.2178  0.23420  0.36069       0.95 < 2.2e-16 ***
      mean 1         0.540531 0.036574 14.7790  0.46885  0.61222       0.95 < 2.2e-16 ***
      rd 1 vs 0      0.243086 0.048386  5.0239  0.14825  0.33792       0.95 5.064e-07 ***
      log(rr) 1 vs 0 0.597322 0.126936  4.7057  0.34853  0.84611       0.95 2.530e-06 ***
      log(or) 1 vs 0 1.021974 0.211670  4.8281  0.60711  1.43684       0.95 1.378e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust"))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Marginal estimates:
                     estimate  std.err       z ci.lower ci.upper conf.level   p.value    
      mean 0         0.297445 0.032318  9.2036  0.23410  0.36079       0.95 < 2.2e-16 ***
      mean 1         0.540531 0.037526 14.4041  0.46698  0.61408       0.95 < 2.2e-16 ***
      rd 1 vs 0      0.243086 0.049525  4.9084  0.14602  0.34015       0.95 9.183e-07 ***
      log(rr) 1 vs 0 0.597322 0.128939  4.6326  0.34461  0.85004       0.95 3.611e-06 ***
      log(or) 1 vs 0 1.021974 0.216214  4.7267  0.59820  1.44575       0.95 2.282e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Standard errors: robust, a diagnostic that treats the weights as known

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

# the linearization separation error reads in the user's terms

    Code
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
    Condition
      Error in `ipw()`:
      ! `wt_mod` must not separate the exposure.
      x Putting the fitted linear predictors through the link's inverse gives a probability of exactly 0 or 1 for <n> observations, whose weights are then undefined.
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

# a robust result has no conditional reading

    Code
      expr
    Condition <propensity_no_conditional_vcov_error>
      Error in `tidy()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result records none.
      x The "robust" standard error method treats the weights as known and reports what the outcome model computes for itself, so it stores no covariance for the conditional reading.
      i Fit with `se_method = "mestimation"`, which solves the two models as one system and stores that covariance.

# ipw() refuses .by with robust standard errors

    Code
      expr
    Condition <propensity_ipw_by_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` with "robust" standard errors.
      x The stratum effects and their contrasts are parameters of the stacked estimating equations, and the robust path solves no such system.
      i Use `se_method = "mestimation"` to report effects by `.by`.

# robust rejects the atu estimand as linearization does

    Code
      expr
    Condition <propensity_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support "robust" standard errors for the "atu" estimand.
      i Use `se_method = "mestimation"` for the "atu" estimand.

# robust rejects a covariate-adjusted outcome model as linearization does

    Code
      expr
    Condition <propensity_method_error>
      Error in `ipw()`:
      ! `ipw()` supports "robust" standard errors only for an outcome model of the exposure alone.
      x `outcome_mod` is adjusted for terms beyond "z".
      i Use `se_method = "mestimation"` for a covariate-adjusted outcome model.

