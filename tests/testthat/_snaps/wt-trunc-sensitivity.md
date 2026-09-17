# a lower grid of any other length is refused

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_trunc_sensitivity()`:
      ! `lower` must have length 1 or the length of `upper`.
      x `lower` has length 2 and `upper` has length 4.

# the absolute and count grids need upper, as wt_trunc() does

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `wt_trunc_sensitivity()`:
      ! For `method = 'wt'`, `upper` is required.
      i Supply the grid of upper bounds, or use `method = 'pctl'` to read a default grid of quantiles.

# one invalid bound refuses the whole grid with wt_trunc()'s class

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc_sensitivity()`:
      ! For `method = 'pctl'`, `upper` must be between 0.5 and 1.
      x `upper` is 1.5.
      i The refused bound is at position 2 of the grid.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc_sensitivity()`:
      ! `lower` must be smaller than `upper`.
      x `lower` is 3 and `upper` is 2.5.
      i The refused bound is at position 2 of the grid.

---

    Code
      expr
    Condition <propensity_missing_value_error>
      Error in `wt_trunc_sensitivity()`:
      ! `upper` must not be missing.
      i A missing bound decides nothing about which weights to move. Supply a value, or leave `lower` unset to bound only the upper tail.
      i The refused bound is at position 2 of the grid.

# an empty grid is refused

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_trunc_sensitivity()`:
      ! `upper` must have at least one bound.
      x `upper` has length 0, so there is no grid to tabulate.

# the adaptive method is not a grid method

    Code
      expr
    Condition <rlang_error>
      Error in `wt_trunc_sensitivity()`:
      ! `method` must be one of "pctl", "wt", or "count", not "adaptive".

# a grid over no weights is refused under every method

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc_sensitivity()`:
      ! `.weights` must have at least one weight present.
      x Of the 0 weights given, none are present, so there is no range or mean to report.

# propensity scores are refused, as wt_trunc() refuses them

    Code
      expr
    Condition <propensity_type_error>
      Error in `wt_trunc_sensitivity()`:
      ! `.weights` must be weights, not propensity scores.
      x It is <ps_trim>, a modified propensity score.
      i Build weights from it with a weight function such as `wt_ate()` and truncate those, or bound the scores themselves with `ps_trunc()`.

# a matrix of weights is refused

    Code
      expr
    Condition <propensity_type_error>
      Error in `wt_trunc_sensitivity()`:
      ! `.weights` must be a vector of weights, one per unit.
      x It has dimensions 2 x 2.
      i Pass a <psw> vector, such as one `wt_ate()` returns, or a plain numeric vector.

# a grid that is not an atomic vector is refused

    Code
      expr
    Condition <propensity_type_error>
      Error in `wt_trunc_sensitivity()`:
      ! `upper` must be a numeric vector of bounds.
      x It is <list>.

# weights that are already truncated are refused

    Code
      expr
    Condition <propensity_already_modified_error>
      Error in `wt_trunc_sensitivity()`:
      ! `.weights` must not already be truncated.
      x The weights have already been bounded with `wt_trunc()`, so the first row of the grid would describe weights that were already bounded rather than the weights before truncation.
      i Pass the weights as they were before `wt_trunc()`.

