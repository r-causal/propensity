# ipw works for binary outcome with a confounder, using logistic ps, logistic outcome

    Code
      res
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
      
      Estimates:
              estimate  std.err        z ci.lower ci.upper conf.level      p.value
      rd       0.19988 0.092425 2.162637   0.0187  0.38103       0.95    0.0305691
      log(rr)  0.56041 0.156172 3.588443   0.2543  0.86651       0.95    0.0003327
      log(or)  0.87831 0.173946 5.049330   0.5374  1.21924       0.95 0.0000004434
                 
      rd      *  
      log(rr) ***
      log(or) ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# ipw works for continuous outcome with a confounder, using logistic ps, linear outcome

    Code
      res
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: lm(formula = y ~ z, data = dat, weights = wts) 
      
      Estimates:
           estimate  std.err        z ci.lower ci.upper conf.level
      diff   2.2526  0.17524 12.85419   1.9091    2.596       0.95
                         p.value    
      diff < 0.00000000000000022 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

