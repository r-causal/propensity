# ipw() categorical print output is stable

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
        Call: nnet::multinom(formula = a ~ x1 + x2, data = dat, trace = FALSE, 
          reltol = 1e-14, maxit = 2000) 
      
      Outcome Model:
        Call: glm(formula = fmla, family = quasibinomial(), data = dat, weights = wts, 
          control = glm.control(epsilon = 1e-14, maxit = 200)) 
      
      Estimates:
                     estimate  std.err      z  ci.lower ci.upper conf.level   p.value
      rd b vs a      0.081945 0.050387 1.6263 -0.016811  0.18070       0.95 0.1038822
      log(rr) b vs a 0.168870 0.104633 1.6139 -0.036207  0.37395       0.95 0.1065433
      log(or) b vs a 0.328762 0.203058 1.6191 -0.069225  0.72675       0.95 0.1054363
      rd c vs a      0.166939 0.045182 3.6948  0.078384  0.25549       0.95 0.0002200
      log(rr) c vs a 0.318293 0.091898 3.4635  0.138176  0.49841       0.95 0.0005331
      log(or) c vs a 0.676435 0.185786 3.6409  0.312300  1.04057       0.95 0.0002717
                        
      rd b vs a         
      log(rr) b vs a    
      log(or) b vs a    
      rd c vs a      ***
      log(rr) c vs a ***
      log(or) c vs a ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the categorical degenerate-design error names the pinned levels

    Code
      ipw(mods$ps_mod, collapsed, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must be able to represent the outcome at every exposure level.
      x Setting `a` to "a" and "b" leaves the counterfactual design identically zero, which pins the marginal mean there to the outcome link's zero point instead of estimating it.
      i Include an intercept in `outcome_mod`, or code the exposure as a factor, whose no-intercept coding is saturated and represents every level.

