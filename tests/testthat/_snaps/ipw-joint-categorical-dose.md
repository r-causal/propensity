# the three-level coding refusal names each term once

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `z` and `z:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "z" enters as one 0/1 indicator per non-reference level, "mid" and "hi" in that order and each against "lo", "e" enters as itself, and each interaction column is the product of one indicator with it.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i A model with no intercept, written `- 1` or `+ 0`, expands a factor treatment to an indicator for every level rather than for the non-reference levels alone, so its columns are shifted by one against the indicators the rows describe. Keep the intercept.
      i Refit `outcome_mod` with "z" as an unordered factor under treatment contrasts. A treatment with more than two levels has no numeric coding whose bare term contributes those columns.

# a dose column rescaled after the treatment models are fit says so

    Code
      expr
    Condition <propensity_ipw_msm_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint intervention with a dose from a marginal structural model whose treatment columns are the treatments themselves.
      x `e` and `z:e` in `outcome_mod` contribute a column coded some other way.
      i The reported rows name the coefficients of a model in which "z" enters as one 0/1 indicator per non-reference level, "mid" and "hi" in that order and each against "lo", "e" enters as itself, and each interaction column is the product of one indicator with it.
      i A contrast coding other than treatment contrasts rescales or recenters those columns without changing what the formula says. An ordered factor carries polynomial contrasts, and `options(contrasts = )` sets a coding for every factor in the session.
      i A model with no intercept, written `- 1` or `+ 0`, expands a factor treatment to an indicator for every level rather than for the non-reference levels alone, so its columns are shifted by one against the indicators the rows describe. Keep the intercept.
      i Refit `outcome_mod` with "z" as an unordered factor under treatment contrasts. A treatment with more than two levels has no numeric coding whose bare term contributes those columns.

