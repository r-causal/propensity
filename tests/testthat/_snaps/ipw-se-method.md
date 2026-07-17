# ipw() print output is stable per SE method

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Estimates:
              estimate  std.err        z ci.lower ci.upper conf.level   p.value    
      rd       0.24309 0.048386 5.023886   0.1483  0.33792       0.95 5.064e-07 ***
      log(rr)  0.59732 0.126936 4.705681   0.3485  0.84611       0.95 2.530e-06 ***
      log(or)  1.02197 0.211670 4.828145   0.6071  1.43684       0.95 1.378e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      print(ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"))
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Estimates:
              estimate  std.err        z ci.lower ci.upper conf.level   p.value    
      rd       0.24309 0.048386 5.023886   0.1483  0.33792       0.95 5.064e-07 ***
      log(rr)  0.59732 0.126936 4.705681   0.3485  0.84611       0.95 2.530e-06 ***
      log(or)  1.02197 0.211670 4.828145   0.6071  1.43684       0.95 1.378e-06 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

