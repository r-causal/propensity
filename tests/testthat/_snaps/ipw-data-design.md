# the factored-exposure error names the term and the rebuild it blocks

    Code
      ipw(ps_mod, out, .data = dat)
    Condition
      Error in `ipw()`:
      ! `outcome_mod` must read the exposure from a column `ipw()` can set to counterfactual values.
      x `outcome_mod` reads "z" through the term `factor(z)`.
      x `ipw()` estimates the marginal means by setting "z" to one value at a time and rebuilding the outcome design, and `factor(z)` does not carry the levels "0" and "1" it was fit with once the column is constant.
      i Refit `outcome_mod` on the plain "z" column, converting it to a factor in the data first if you want a factor coding.

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
      i Supply "y" as that numeric column, or refit `outcome_mod` on the character vector.

