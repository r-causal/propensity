# the refusal of a density reads the same way for either type

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A binary exposure has a probability rather than a density, so leave `.density` unset for one.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as categorical.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A categorical exposure has a probability rather than a density, so leave `.density` unset for one.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_cens()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A binary exposure has a probability rather than a density, so leave `.density` unset for one.

