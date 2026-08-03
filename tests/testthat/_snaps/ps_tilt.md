# ps_tilt() sends a modified propensity score to its plain scores

    Code
      ps_tilt(trimmed, "att")
    Condition
      Error in `ps_tilt()`:
      ! No method for objects of class <ps_trim/vctrs_vctr/double>.
      i `ps_tilt()` takes plain propensity scores. A score modified by `ps_trim()`, `ps_trunc()`, or `ps_calibrate()` carries a class of its own; pass the scores underneath it, with `as.numeric(x)` for a binary exposure or `as.matrix(x)` for a categorical one.

