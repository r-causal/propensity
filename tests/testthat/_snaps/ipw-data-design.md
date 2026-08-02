# the factored-exposure error names the term and the route that rebuilds it

    Code
      ipw(ps_mod, out)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must read the exposure from a column `ipw()` can set to counterfactual values.
      x `outcome_mod` reads "z" through the term `factor(z)`.
      x Without `.data` the designs come from `outcome_mod`'s own model frame, which holds `factor(z)` at the values it was fit on, so the counterfactual value `ipw()` sets is ignored and every level is given the fitted design.
      i Supply `.data`, which rebuilds `factor(z)` at each value under the levels "0" and "1" it was fit with, or refit `outcome_mod` on the plain "z" column.

# the mirrored type error names the column and both types

    Code
      ipw(ps_mod, out, .data = dat_fac)
    Condition
      Error in `ipw()`:
      ! `.data` must supply "num" as the numeric column the models were fit with.
      x `.data` has "num" as a factor.
      x `outcome_mod` recorded "num" as a numeric vector, and the designs rebuilt from `.data` use that coding.
      i Supply "num" as that numeric column, or refit the models on the factor.

# the character outcome error names the column and both types

    Code
      ipw(ps_mod, out, .data = dat_chr, se_method = "mestimation")
    Condition
      Error in `ipw()`:
      ! `.data` must supply the outcome "y" as the numeric column `outcome_mod` was fit on.
      x `.data` has "y" as a character vector.
      x `ipw()` reads the outcome values from `.data`: a factor response becomes an indicator for its non-first levels and any other response is used as its own values, so the two are not interchangeable.
      i Supply "y" as that numeric column.

# the new-level error names the column and the values the fit never saw

    Code
      ipw(ps_mod, out, .data = dat_new)
    Condition
      Error in `ipw()`:
      ! `.data` must hold only the levels the models were fit on.
      x `cov` takes the value "extra" in `.data`.
      x The models were fit with `cov` at "low", "mid", and "high", and the designs rebuilt from `.data` carry one column per non-reference level of that set, so a value outside it has no coefficient to multiply.
      i Supply the data the models were fit to, or refit them on data that holds every level you want represented.

# the re-leveled response error names the column and the coding

    Code
      ipw(ps_mod, out, .data = dat_flip)
    Condition
      Error in `ipw()`:
      ! `.data` must supply the outcome "yf" with the level `outcome_mod` was fit to treat as the failure first.
      x `.data` declares "yes" first; `outcome_mod` was fit with "no" first.
      x `ipw()` codes a factor outcome as an indicator for its non-first levels, following `glm()`, so the two orders describe opposite outcomes and every estimate reverses.
      i Re-level "yf" so "no" comes first, or supply the data the models were fit to.

