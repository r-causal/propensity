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
              estimate  std.err      z ci.lower ci.upper conf.level p.value  
      rd      0.199882 0.092425 2.1626 0.018732  0.38103       0.95 0.03057 *
      log(rr) 0.560414 0.273519 2.0489 0.024326  1.09650       0.95 0.04047 *
      log(or) 0.878313 0.418661 2.0979 0.057753  1.69887       0.95 0.03591 *
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
           estimate std.err      z ci.lower ci.upper conf.level   p.value    
      diff  2.25255 0.17524 12.854   1.9091    2.596       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# ps_mod must be glm, outcome_mod must be glm or lm

    Code
      expr
    Condition <propensity_method_error>
      Error in `ipw()`:
      ! `ipw()` does not know how to handle `ps_mod` of class <not_a_model>.
      i `ps_mod` must be a fitted propensity score model of class <glm>, such as a logistic regression.

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
      Error in `ipw()`:
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

# the cannot-determine-estimand error is attributed to ipw()

    Code
      ipw(ps_mod, outcome_mod, .data = dat)
    Condition
      Error in `ipw()`:
      ! Can't determine the estimand from weights.
      i Please specify `estimand`.

