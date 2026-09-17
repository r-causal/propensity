# assigning weights with a different score explains the refusal

    Code
      y[1:3] <- rev(w)
    Condition
      Error in `vec_cast.psw.psw()`:
      ! Can't convert <psw{estimand = ate; stabilized}> to <psw{estimand = ate; stabilized}>.
      different stabilization scores
      The value's per-observation `stabilization_score` does not match the target's, and a score describes the units of the weights that carry it.
      Assign `vec_data()` of the value to keep the target's score, or rebuild the weights once the values are assigned.

