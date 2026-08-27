# the rank-deficient design errors read in the user's terms

    Code
      ipw(ps_mod, out, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must have a coefficient for every column of its design.
      x `outcome_mod` has no fitted coefficient for "x1_copy".
      i A model reports that for a column its design cannot separate from the others: the column is a linear combination of them, exactly or to within the tolerance the fit pivots at, so the fit has no unique solution for it and drops it.
      i `ipw()` estimates the marginal means by multiplying the fitted coefficients against that design, so a column with no coefficient leaves every predicted outcome undefined.
      i Refit `outcome_mod` without the redundant column, or combine it with the column it duplicates.

# the propensity rank error reads in the user's terms

    Code
      ipw(ps_mod, out, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! `wt_mod` must have a coefficient for every column of its design.
      x `wt_mod` has no fitted coefficient for "x1_copy".
      i A model reports that for a column its design cannot separate from the others: the column is a linear combination of them, exactly or to within the tolerance the fit pivots at, so the fit has no unique solution for it and drops it.
      i `ipw()` rebuilds the propensity scores by multiplying the fitted coefficients against that design, so a column with no coefficient leaves every score undefined.
      i Refit `wt_mod` without the redundant column, or combine it with the column it duplicates.

# the unsolved-equations warning reads in the user's terms

    Code
      ipw(mods$ps_mod, mods$out, se_method = "mestimation")
    Condition
      Warning in `ipw()`:
      The estimating equations behind these estimates have no unique root at the values the solver returned.
      x The standard errors reported for "mean for 0", "mean for 1", "rd for 1 vs 0", "log(rr) for 1 vs 0", and "log(or) for 1 vs 0" are not meaningful: they collapsed to essentially zero, which makes the test statistics and the intervals built from them meaningless too.
      i At least one direction in the parameter space leaves the equations unchanged, so the sandwich variance along it is not identified.
      i An exposure level in which every outcome is an event, or none is, is one cause: the outcome model has no finite fit within it. Check the outcome within each level of the exposure, and both designs for columns that duplicate one another.
      i The estimates are reported as the solver returned them.
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y0 ~ z + x1, family = quasibinomial(), data = dat, 
          weights = wts) 
      
      Marginal estimates:
                       estimate    std.err          z   ci.lower   ci.upper conf.level   p.value    
      mean 0         2.6064e-09 <degenerate> <degenerate> 2.6064e-09 2.6064e-09       0.95 < 2.2e-16 ***
      mean 1         6.5625e-01 <degenerate> <degenerate> 6.5625e-01 6.5625e-01       0.95 < 2.2e-16 ***
      rd 1 vs 0      6.5625e-01 <degenerate> <degenerate> 6.5625e-01 6.5625e-01       0.95 < 2.2e-16 ***
      log(rr) 1 vs 0 1.9344e+01 <degenerate> <degenerate> 1.9344e+01 1.9344e+01       0.95 < 2.2e-16 ***
      log(or) 1 vs 0 2.0412e+01 <degenerate> <degenerate> 2.0412e+01 2.0412e+01       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# the no-variance error reads in the user's terms

    Code
      ipw(mods$ps_mod, mods$out, se_method = "mestimation")
    Condition
      Warning:
      The "log(rr)" effect for "c vs a" is undefined at the marginal means the solver reached.
      i "log(rr)" needs a positive marginal mean on each side of the comparison, and at least one side is outside that range.
      i An exposure level whose fitted outcomes are all events, or all non-events, drives its marginal mean to the boundary. Check the outcome within each level of the exposure.
      i Estimates and standard errors from this fit are not reliable.
      Warning:
      The "log(or)" effect for "c vs a" is undefined at the marginal means the solver reached.
      i "log(or)" needs a marginal mean strictly between 0 and 1 on each side of the comparison, and at least one side is outside that range.
      i An exposure level whose fitted outcomes are all events, or all non-events, drives its marginal mean to the boundary. Check the outcome within each level of the exposure.
      i Estimates and standard errors from this fit are not reliable.
      Error in `ipw()`:
      ! Can't compute a variance for this fit.
      x The outcome does not vary within the c level, which holds no events.
      i The outcome model has no finite fit within an arm whose outcome never varies, so the stacked estimating equations have no finite derivative at the values the solver returned and there is no variance to build from them.
      i Drop that arm and estimate the effect among the ones that are left, or model an outcome that varies within every arm.
      i `ipw()` reports estimates with inference or not at all, so no estimates are returned.

# the undiagnosable no-variance error reads in the user's terms

    Code
      ipw(mods$ps_mod, mods$out, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! Can't compute a variance for this fit.
      i The stacked estimating equations have no finite derivative at the values the solver returned, so there is no variance to build from them.
      i Check the outcome within each level of the exposure, and both designs for columns that all but duplicate one another: a column the fit kept a coefficient for can still leave the equations flat enough here that their derivatives are not finite.
      i `ipw()` reports estimates with inference or not at all, so no estimates are returned.

# the collapsed standard error report reads in the user's terms

    Code
      ipw(mods$ps_mod, mods$out, se_method = "mestimation")
    Condition
      Warning in `ipw()`:
      The standard errors reported for "mean for 0", "mean for 1", and "diff for 1 vs 0" are not meaningful.
      x They are zero, or so small beside the estimates they accompany that the test statistics and the intervals built from them carry no information.
      i An exposure group the outcome does not vary within is one cause: the contrast is then a fixed value rather than a quantity with any spread. Check the outcome within each level of the exposure.
      i The estimates are reported as they were computed.
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: lm(formula = yconst ~ z, data = dat, weights = wts) 
      
      Marginal estimates:
                     estimate     std.err           z    ci.lower  ci.upper conf.level   p.value    
      mean 0      -<degenerate>  <degenerate> -1.9959e+01 -<degenerate> -2.01e-16       0.95 < 2.2e-16 ***
      mean 1       1.0000e+00  <degenerate>  <degenerate>  1.0000e+00  1.00e+00       0.95 < 2.2e-16 ***
      diff 1 vs 0  1.0000e+00  <degenerate>  <degenerate>  1.0000e+00  1.00e+00       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

