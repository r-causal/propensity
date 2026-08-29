# ipw() refuses an outcome model that does not read both treatments

    Code
      expr
    Condition <propensity_ipw_exposure_error>
      Error in `ipw()`:
      ! `outcome_mod` must contain the exposure "e".
      x "e" does not appear in the formula of `outcome_mod`.
      i Refit `outcome_mod` with "e" on the right-hand side of the formula.

# ipw() refuses weights that are not a joint product

    Code
      expr
    Condition <propensity_ipw_joint_models_weights_error>
      Error in `ipw()`:
      ! `outcome_mod` must be fit with the product weights the two treatment models imply.
      x Its weights carry no record of being a product, so they weight one treatment rather than the crossing of the two.
      i A single treatment's weights are an ordinary <psw> of the right length and the right estimand, so nothing else would notice.
      i Build the weights with `wt_joint()` from the two components, and refit `outcome_mod` with them.

# ipw() refuses .by on the two-model route

    Code
      expr
    Condition <propensity_ipw_joint_by_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a joint exposure.
      x A joint exposure already reports an interaction between two treatments, and reporting it within the levels of a modifier is a three-way question this surface does not answer.
      i Drop `.by` to report the joint surface, or weight the crossing of the two treatments as one plain categorical exposure to report each cell against the reference cell within each subgroup.

# ipw() refuses linearization standard errors on the two-model route

    Code
      expr
    Condition <propensity_ipw_joint_models_method_error>
      Error in `ipw()`:
      ! `ipw()` does not support "linearization" standard errors for a joint treatment model.
      x The cell means and every contrast built from them are parameters of the stacked estimating equations, and the linearization path solves no such system.
      i Use `se_method = "mestimation"` for a joint treatment model.

