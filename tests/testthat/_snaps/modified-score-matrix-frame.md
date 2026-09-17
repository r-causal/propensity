# ps_refit() refuses a single column of a trimmed matrix's data frame

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_refit()`:
      ! `ps_refit()` cannot refit one column of a trimmed score matrix.
      x These scores carry the record of a categorical trim, which chose the retained rows from every column at once.
      i Refit the trimmed matrix with `ps_refit()` before converting it with `as.data.frame()`.

