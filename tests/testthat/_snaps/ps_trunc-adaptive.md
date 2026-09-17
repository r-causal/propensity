# the adaptive bound on a vector prints its method and bounds

    Code
      print(ps_trunc(adaptive_ps, method = "adaptive"))
    Output
      <ps_trunc{[0.373,0.627], method=adaptive}[20]>
       [1] 0.3732089 0.3732089 0.3732089 0.5000000 0.6000000 0.6200000 0.4500000
       [8] 0.5500000 0.4000000 0.5000000 0.6267911 0.6267911 0.6267911 0.5000000
      [15] 0.4000000 0.3800000 0.4500000 0.5500000 0.6000000 0.3900000

# bounds supplied to the adaptive method are ignored with a warning

    Code
      out <- ps_trunc(adaptive_ps, method = "adaptive", lower = 0.2, upper = 0.8)
    Condition <propensity_warning>
      Warning in `ps_trunc()`:
      For `method = 'adaptive'`, `lower` and `upper` are ignored.
      i The adaptive bounds are set by the number of propensity scores. To give bounds yourself, use `method = 'ps'`.

# the adaptive bound needs at least 15 scores

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! For `method = 'adaptive'`, at least 15 propensity scores must be present.
      x 14 scores are present.
      i The bounds are 1/c and 1 - 1/c with c = sqrt(n) log(n) / 5, and below 15 scores c is under 2, so the lower bound would lie above the upper one.
      i Supply bounds yourself with `method = 'ps'`.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trunc()`:
      ! For `method = 'adaptive'`, at least 15 propensity scores must be present.
      x 6 scores are present.
      i The bounds are 1/c and 1 - 1/c with c = sqrt(n) log(n) / 5, and below 15 scores c is under 2, so the lower bound would lie above the upper one.
      i Supply bounds yourself with `method = 'ps'`.

# the adaptive bound refuses a categorical score matrix

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trunc()`:
      ! Method "adaptive" is not supported for categorical exposures.
      i Bounding a categorical score renormalizes each row, which moves the scores the bound never reached, so a bound on the scores is no longer a bound on the weights.
      i Use `ps_trim()` with `method = "optimal"`, or build the weights and bound them with `wt_trunc()`.

# the adaptive bound refuses a multinomial fit of three levels

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trunc()`:
      ! Method "adaptive" is not supported for categorical exposures.
      i Bounding a categorical score renormalizes each row, which moves the scores the bound never reached, so a bound on the scores is no longer a bound on the weights.
      i Use `ps_trim()` with `method = "optimal"`, or build the weights and bound them with `wt_trunc()`.

