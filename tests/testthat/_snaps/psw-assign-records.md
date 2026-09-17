# assigning weights truncated at another bound is refused

    Code
      target[1:3] <- value[1:3]
    Condition
      Error in `vec_cast.psw.psw()`:
      ! Can't convert <psw{estimand = unknown}> to <psw{estimand = unknown}>.
      different weight truncation bounds
      Assign `vec_data()` of the value to keep the target's record, or rebuild the weights from one modification.

