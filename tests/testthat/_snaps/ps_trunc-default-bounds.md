# a lone bound whose mirror crosses it is refused

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! `lower` must be smaller than `upper`
      x `lower` is 0.6 and `upper` is 0.4
      i `upper` was not supplied, so it is the mirror of `lower`, 1 - `lower`. Supply both to bound the scores asymmetrically.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! `lower` must be smaller than `upper`
      x `lower` is 0.6 and `upper` is 0.4
      i `lower` was not supplied, so it is the mirror of `upper`, 1 - `upper`. Supply both to bound the scores asymmetrically.

# the default floor meets the 1/k refusal at ten or more levels

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! The truncation threshold must fall below 1/k, for k columns of propensity scores.
      x `lower` is 0.1, and 1/k is 0.1 for the 10 columns the scores hold.
      i No row summing to one can hold every score above 1/k, so a threshold there leaves no rule to apply.
      i `lower` was not supplied, and 0.1 is its default. Supply a `lower` below 0.1.

