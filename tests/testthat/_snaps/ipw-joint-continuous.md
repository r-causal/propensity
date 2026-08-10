# a marginal structural model reading a covariate is still refused

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a marginal structural model whose exposure terms read the exposure alone.
      x `e:x1` is not a term of "e" alone.
      i A term reading anything else contributes a coefficient that depends on what it reads, so there is no one effect for a row to report.
      i A term reading the exposure alone is admitted however it is written, so a curve such as `e + I(e^2)` reports one row per coefficient.
      i Read the full coefficient vector from the returned fit object for a model this surface cannot report.

---

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` requires a marginal structural model whose exposure terms read the exposure alone.
      x `a:e` is not a term of "e" alone.
      i A term reading anything else contributes a coefficient that depends on what it reads, so there is no one effect for a row to report.
      i A term reading the exposure alone is admitted however it is written, so a curve such as `e + I(e^2)` reports one row per coefficient.
      i Read the full coefficient vector from the returned fit object for a model this surface cannot report.

# ipw() refuses a treatment column coded some other way

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `a` and `a:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "a" enters as 0 for "no" and 1 for "yes", "e" enters as itself, and their interaction is the product of the two.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i Refit `outcome_mod` with "a" as a 0/1 numeric, or as an unordered factor under treatment contrasts.

---

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `a` and `a:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "a" enters as 0 for "no" and 1 for "yes", "e" enters as itself, and their interaction is the product of the two.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i Refit `outcome_mod` with "a" as a 0/1 numeric, or as an unordered factor under treatment contrasts.

# ipw() refuses .by on a joint continuous fit

    Code
      expr
    Condition <propensity_ipw_joint_by_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a joint exposure.
      x A joint exposure already reports an interaction between two treatments, and reporting it within the levels of a modifier is a three-way question this surface does not answer.
      i Drop `.by` to report the joint surface, or weight the crossing of the two treatments as one plain categorical exposure to report each cell against the reference cell within each subgroup.

# ipw() refuses linearization on a joint continuous fit

    Code
      expr
    Condition <propensity_ipw_joint_models_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support "linearization" standard errors for a joint treatment model.
      x The cell means and every contrast built from them are parameters of the stacked estimating equations, and the linearization path solves no such system.
      i Use `se_method = "mestimation"` for a joint treatment model.

