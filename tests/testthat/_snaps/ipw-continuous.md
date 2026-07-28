# ipw() continuous print output is stable

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
        Call: lm(formula = A ~ x1 + x2, data = dat) 
      
      Outcome Model:
        Call: lm(formula = msm_fmla, data = dat, weights = wts) 
      
      Estimates:
            estimate  std.err      z ci.lower ci.upper conf.level   p.value    
      slope 0.679126 0.037942 17.899  0.60476  0.75349       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the continuous propensity-link error names the unsupported link

    Code
      ipw(ps_mod, msm)
    Condition
      Error in `ipw()`:
      ! `ipw()` supports only an identity-link propensity score model for a continuous exposure.
      x `ps_mod` is a gaussian model with a "log" link.
      i Refit `ps_mod` as an `lm()` or a gaussian glm with an identity link.

