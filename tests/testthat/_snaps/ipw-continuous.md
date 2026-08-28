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
      x A continuous propensity score model has no link for `ps_link` to override.
      i Omit `ps_link`; it applies only to a binomial glm propensity score model.

# the transformed continuous propensity response error reads in the user's terms

    Code
      ipw(ps_mod, outcome_mod, .data = dat)
    Condition
      Error in `ipw()`:
      ! `ipw()` does not support a transformed response in the propensity score model.
      x `wt_mod` reads the exposure through `log(Apos)`, an expression rather than a single column.
      i Fit `wt_mod` with the exposure itself as the response, adding it to the data as its own column first if it has to be computed.

