# augment(data = ) rejects data holding a column augment adds

    Code
      expr
    Condition <propensity_augment_column_error>
      Error in `augment()`:
      ! `augment()` adds a column the frame carrying it already holds.
      x The frame already holds ".fitted".
      i Rename or drop it in `data`, or leave `data` as `NULL` to carry the fit's columns on the outcome model's own frame.

# augment(data = ) names every column of the data it clashes with

    Code
      expr
    Condition <propensity_augment_column_error>
      Error in `augment()`:
      ! `augment()` adds columns the frame carrying them already holds.
      x The frame already holds ".propensity", ".weights", and ".fitted".
      i Rename or drop them in `data`, or leave `data` as `NULL` to carry the fit's columns on the outcome model's own frame.

# augment(data = ) clashes on the level columns of a categorical fit

    Code
      expr
    Condition <propensity_augment_column_error>
      Error in `augment()`:
      ! `augment()` adds a column the frame carrying it already holds.
      x The frame already holds ".propensity_b".
      i Rename or drop it in `data`, or leave `data` as `NULL` to carry the fit's columns on the outcome model's own frame.

# augment() rejects a result whose models saw different counts

    Code
      expr
    Condition <propensity_augment_alignment_error>
      Error in `augment()`:
      ! `x` must hold two models fit to the same rows.
      x The propensity score model produced 599 propensity scores.
      x The outcome model was fit to 600 observations.
      i Two models of the same data most often part this way over missing values, each dropping the rows a variable of its own is missing on.
      i Refit both models on the rows that are complete on every variable either model uses.

# augment() rejects a result whose models saw different rows

    Code
      expr
    Condition <propensity_augment_alignment_error>
      Error in `augment()`:
      ! `x` must hold two models fit to the same rows.
      x Both models were fit to 599 observations, and the "z" values in the two model frames disagree.
      i Two models of the same data most often part this way over missing values, each dropping the rows a variable of its own is missing on.
      i Refit both models on the rows that are complete on every variable either model uses.

# augment() rejects a result whose outcome model has no weights

    Code
      expr
    Condition <propensity_ipw_weights_missing_error>
      Error in `augment()`:
      ! `x` must hold an outcome model fit with propensity score weights.
      x The outcome model `x` holds reports no weights.
      i An `ipw()` outcome model is weighted by construction. Refit it with the weights the propensity score model implies, such as those `wt_ate()` returns.

# tidy() has no conditional reading of a linearization fit

    Code
      expr
    Condition <propensity_no_conditional_vcov_error>
      Error in `tidy()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result records none.
      x The "linearization" standard error method corrects the marginal estimates only, so it stores no covariance for the conditional reading.
      i Fit with `se_method = "mestimation"`, which solves the two models as one system and stores that covariance.

