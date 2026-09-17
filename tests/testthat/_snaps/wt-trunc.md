# bounds supplied to the adaptive method are ignored with a warning

    Code
      out <- wt_trunc(w, method = "adaptive", upper = 3)
    Condition <propensity_warning>
      Warning in `wt_trunc()`:
      For `method = 'adaptive'`, `lower` and `upper` are ignored.
      i The adaptive bound is set by the number of weights. To give a bound yourself, use `method = 'wt'`.

# method = 'wt' refuses a bound it cannot apply

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `wt_trunc()`:
      ! For `method = 'wt'`, `upper` is required.
      i Supply the bound, or use `method = 'adaptive'` to have it set from the number of weights.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'wt'`, `upper` must be a positive weight.
      x `upper` is 0.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! `lower` must be smaller than `upper`.
      x `lower` is 4 and `upper` is 3.

# method = 'pctl' refuses an upper probability below one half and names both readings

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'pctl'`, `upper` must be between 0.5 and 1.
      x `upper` is 0.1, which is in the lower tail.
      i To bound the upper tail write `upper = 0.9`; to bound the lower tail write `lower = 0.1`.

# method = 'pctl' refuses probabilities outside their half of the unit interval

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'pctl'`, `upper` must be between 0.5 and 1.
      x `upper` is 1.5.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'pctl'`, `lower` must be between 0 and 0.5.
      x `lower` is 0.6.

# method = 'count' refuses counts that meet or cross

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'count'`, `lower` and `upper` must leave at least two weights between them.
      x `lower` is 4, `upper` is 4, and 9 weights are present.

# method = 'count' refuses a count it cannot apply

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'count'`, `upper` must be a whole number of at least 1.
      x `upper` is 1.5.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'count'`, `upper` must leave a weight below the ones it bounds.
      x `upper` is 9 and 9 weights are present.

# a subset that drops the record prints a truncation line without counts

    Code
      wt_trunc(w, method = "wt", upper = 3)[1:3]
    Output
      <psw{estimand = ate; weights truncated}[3]>
      [1] 0.5 1.0 3.0
      truncation: weights truncated

# wt_trunc() refuses input that is not weights

    Code
      expr
    Condition <propensity_type_error>
      Error in `wt_trunc()`:
      ! `.weights` must be a numeric vector of weights.
      x It is <character>, which holds no weights to bound.
      i Pass a <psw> vector, such as one `wt_ate()` returns.

# wt_trunc() refuses propensity scores

    Code
      expr
    Condition <propensity_type_error>
      Error in `wt_trunc()`:
      ! `.weights` must be weights, not propensity scores.
      x It is <ps_trim>, a modified propensity score.
      i Build weights from it with a weight function such as `wt_ate()` and truncate those, or bound the scores themselves with `ps_trunc()`.

# the data-driven methods refuse fewer than two weights

    Code
      expr
    Condition <propensity_range_error>
      Error in `wt_trunc()`:
      ! For `method = 'adaptive'`, at least two weights must be present.
      x 1 weight is present.

# truncating truncated weights warns and returns them unchanged

    Code
      out <- wt_trunc(once, method = "wt", upper = 2)
    Condition <propensity_already_modified_warning>
      Warning in `wt_trunc()`:
      These weights have already been truncated. Returning them unchanged.
      i To truncate at a different bound, call `wt_trunc()` on the original weights.

# truncated weights print their bound and count

    Code
      wt_trunc(w, method = "wt", upper = 20)
    Output
      <psw{estimand = ate; weights truncated}[7]>
      [1]  0.50  1.00  4.00  2.00 20.00    NA  0.02
      truncation: wt (upper 20), 1 of 7 weights truncated
    Code
      wt_trunc(w, method = "wt", lower = 0.05, upper = 20)
    Output
      <psw{estimand = ate; weights truncated}[7]>
      [1]  0.50  1.00  4.00  2.00 20.00    NA  0.05
      truncation: wt (lower 0.05, upper 20), 2 of 7 weights truncated
    Code
      wt_trunc(c(0.5, 1, 4, 2, 25), method = "count", upper = 1)
    Output
      <psw{estimand = unknown}[5]>
      [1] 0.5 1.0 4.0 2.0 4.0
      truncation: count 1 (upper 4), 1 of 5 weights truncated

---

    Code
      wt_trunc(weights, method = "pctl", upper = 0.9)
    Output
      <psw{estimand = ate; weights truncated; stabilized}[12]>
       [1] 1.2601667 0.8158751 0.2876066 0.4445126 0.8514472 0.8062735 0.8305737
       [8] 1.2601667 0.6685880 0.7727881 0.6170081 0.8400276
      density:   normal
      numerator: marginal
      sigma:     pooled
      truncation: pctl 0.9 (upper 1.26), 2 of 12 weights truncated

