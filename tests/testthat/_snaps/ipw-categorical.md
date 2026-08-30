# ipw() categorical print output is stable

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: nnet::multinom(formula = a ~ x1 + x2, data = dat, trace = FALSE, 
          reltol = 1e-14, maxit = 2000) 
      
      Outcome Model:
        Call: glm(formula = fmla, family = quasibinomial(), data = dat, weights = wts, 
          control = glm.control(epsilon = 1e-14, maxit = 200)) 
      
      Marginal estimates:
                     estimate std.err     z ci.lower ci.upper conf.level p.value    
      mean a           0.4454  0.0358 12.43   0.3752    0.516       0.95 < 2e-16 ***
      mean b           0.5274  0.0362 14.58   0.4565    0.598       0.95 < 2e-16 ***
      mean c           0.6124  0.0284 21.58   0.5568    0.668       0.95 < 2e-16 ***
      rd b vs a        0.0819  0.0504  1.63  -0.0168    0.181       0.95 0.10388    
      log(rr) b vs a   0.1689  0.1046  1.61  -0.0362    0.374       0.95 0.10654    
      log(or) b vs a   0.3288  0.2031  1.62  -0.0692    0.727       0.95 0.10544    
      rd c vs a        0.1669  0.0452  3.69   0.0784    0.255       0.95 0.00022 ***
      log(rr) c vs a   0.3183  0.0919  3.46   0.1382    0.498       0.95 0.00053 ***
      log(or) c vs a   0.6764  0.1858  3.64   0.3123    1.041       0.95 0.00027 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the skewed-design error names the column and both types

    Code
      ipw(mods$ps_mod, mods$outcome_mod, .data = dat_skew)
    Condition
      Error in `ipw()`:
      ! `.data` must supply "x2" as the numeric column the models were fit with.
      x `.data` has "x2" as a factor.
      x `wt_mod` recorded "x2" as a numeric vector, and the designs rebuilt from `.data` use that coding.
      i Supply "x2" as that numeric column, or refit the models on the factor.

# ipw() rejects a two-level multinom propensity score model

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ipw()`:
      ! `ipw()` needs a propensity score model of three or more exposure levels for a categorical exposure.
      x `wt_mod` was fit to 2 levels ("control" and "treated"), so it fits the single probability of a binary exposure.
      i For a binary exposure, fit the propensity score model with `stats::glm()` and `family = binomial()`; for a categorical one, fit a `nnet::multinom()` to its three or more levels.

# ipw() rejects a multinom fit to a matrix response

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ipw()`:
      ! A propensity score model must be fit to the levels of the exposure.
      x `wt_mod` was fit to a matrix response, which `nnet::multinom()` reads as counts rather than as levels.
      i Refit the model with the exposure factor on the left-hand side, as in `nnet::multinom(exposure ~ x)`.

# the categorical degenerate-design error names the pinned levels

    Code
      ipw(mods$ps_mod, collapsed, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must be able to represent the outcome at every exposure level.
      x Setting `a` to "a" and "b" leaves the counterfactual designs identically zero, which pins the marginal means there to the outcome link's zero point instead of estimating them.
      i Include an intercept in `outcome_mod`, or code the exposure as a factor, whose no-intercept coding is saturated and represents every level.

# the categorical indistinguishable-design error names the pinned levels

    Code
      ipw(mods$ps_mod, indistinguishable, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must be able to distinguish every pair of exposure levels.
      x Setting `a` to "a" and "b" produces identical counterfactual designs, so the model predicts the same outcome at both levels and the contrast between them is degenerate.
      i Code the exposure so `outcome_mod` separates every level, for example as a factor rather than an indicator for a single level.

# the categorical numeric-exposure error names the exposure and the remedy

    Code
      ipw(mods$ps_mod, mods$outcome_mod)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must hold a categorical exposure as a factor.
      x "z" enters `outcome_mod` as a numeric term, which gives it one coefficient instead of one per level.
      i Refit `outcome_mod` after converting "z" to a factor in the data, rather than wrapping it in the formula.

# the categorical level-order error names both orders

    Code
      ipw(mods$ps_mod, releveled)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` and `wt_mod` must code the exposure on the same levels in the same level order.
      x `outcome_mod` was fit on "c", "a", and "b"; `wt_mod` was fit on "a", "b", and "c".
      i Refit `outcome_mod` with "a" factored in the propensity score model's order. A character column is factored alphabetically, so convert it to a factor with that order first.

