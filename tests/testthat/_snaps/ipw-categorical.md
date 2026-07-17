# ipw() categorical print output is stable

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: nnet::multinom(formula = a ~ x1 + x2, data = dat, trace = FALSE, 
          reltol = 1e-14, maxit = 2000) 
      
      Outcome Model:
        Call: glm(formula = fmla, family = quasibinomial(), data = dat, weights = wts, 
          control = glm.control(epsilon = 1e-14, maxit = 200)) 
      
      Estimates:
                     estimate  std.err        z ci.lower ci.upper conf.level
      rd b vs a      0.081945 0.050387 1.626317  -0.0168  0.18070       0.95
      log(rr) b vs a 0.168870 0.104633 1.613927  -0.0362  0.37395       0.95
      log(or) b vs a 0.328762 0.203058 1.619051  -0.0692  0.72675       0.95
      rd c vs a      0.166939 0.045182 3.694815   0.0784  0.25549       0.95
      log(rr) c vs a 0.318293 0.091898 3.463534   0.1382  0.49841       0.95
      log(or) c vs a 0.676435 0.185786 3.640925   0.3123  1.04057       0.95
                       p.value    
      rd b vs a      0.1038822    
      log(rr) b vs a 0.1065433    
      log(or) b vs a 0.1054363    
      rd c vs a      0.0002200 ***
      log(rr) c vs a 0.0005331 ***
      log(or) c vs a 0.0002717 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

