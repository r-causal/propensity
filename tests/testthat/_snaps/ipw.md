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
              estimate  std.err        z ci.lower ci.upper conf.level   p.value    
      rd       0.19988 0.092425 2.162637   0.0187  0.38103       0.95 0.0305691 *  
      log(rr)  0.56041 0.156172 3.588443   0.2543  0.86651       0.95 0.0003327 ***
      log(or)  0.87831 0.173946 5.049330   0.5374  1.21924       0.95 4.434e-07 ***
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
           estimate  std.err        z ci.lower ci.upper conf.level   p.value    
      diff   2.2526  0.17524 12.85419   1.9091    2.596       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# ps_mod must be glm, outcome_mod must be glm or lm

    Code
      expr
    Condition <propensity_class_error>
      Error in `ipw()`:
      ! `ps_mod` must be of class "glm".
      x It has class "lm".

---

    Code
      expr
    Condition <propensity_class_error>
      Error:
      ! `"a"` must be of class "character" and length 2.
      x It has length 1.

---

    Code
      expr
    Condition <propensity_class_error>
      Error:
      ! `"a"` must be one of class "numeric" and "character" and length 2.
      x It has length 1.

---

    Code
      expr
    Condition <propensity_class_error>
      Error in `ipw()`:
      ! `outcome_mod` must be one of class "glm" and "lm".
      x It has class "list".

---

    Code
      expr
    Condition <propensity_columns_exist_error>
      Error in `ipw()`:
      ! The data frame `.data` is missing the "z" and "y" columns.

# ipw handles various errors correctly

    Code
      expr
    Condition <propensity_error>
      Error in `check_estimand()`:
      ! Can't determine the estimand from weights.
      i Please specify `estimand`.

# Estimand mismatch triggers an error if outcome weights differ from user-specified

    Code
      expr
    Condition <propensity_error>
      Error in `ipw()`:
      ! Estimand in weights different from `estimand`: "ate" vs. "att"

# ipw works for cloglog link in the propensity score model

    Code
      expr
    Condition <propensity_error>
      Error in `ipw()`:
      ! `exposure` and `outcome` must be the same length.
      x `exposure` is length 400
      x `outcome` is length 100

---

    Code
      expr
    Condition <propensity_columns_exist_error>
      Error in `ipw()`:
      ! "z" not found in `model.frame(outcome_mod)`.
      i The outcome model may have transformations in the formula.
      i Please specify `.data`

