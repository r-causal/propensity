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
      x `wt_mod` has a matrix response, such as `cbind(successes, failures)`; a binary exposure must be a single-column response.
      i Fit `wt_mod` with a single binary response column.

