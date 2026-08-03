# Identify which units were trimmed

`is_unit_trimmed()` returns a logical vector indicating which
observations were removed by trimming. This is a per-unit query, as
opposed to
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md),
which tests whether the object has been trimmed at all.

The answer comes from the trimming record, which is written for a fixed
set of observations and can both be lost and outlive them. On a
`ps_trim`,
[`vctrs::vec_slice()`](https://vctrs.r-lib.org/reference/vec_slice.html)
and [`c()`](https://rdrr.io/r/base/c.html) drop it, and subassignment
that grows the vector carries it across a length change; see
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for the whole contract. On a
[psw](https://r-causal.github.io/propensity/reference/psw.md) vector
built from trimmed propensity scores, a subset drops it, while
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) re-attaches
it to the shortened weights column of an outcome model fit on these
weights.

`is_unit_trimmed()` therefore checks that the record covers the object
it is given, and raises an error of class
`propensity_missing_meta_error` when it does not, or when an object
marked as trimmed carries no record at all, rather than name trimmed
units at stale positions. Query the `ps_trim` object the record was
written for instead.

That check counts observations, which a reordering does not change, so
it does not catch one. A `ps_trim` reordered through vctrs rather than
through `[`, by `vctrs::vec_slice(x, 5:1)` or
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
keeps a record written for the old order, and a `psw` keeps one through
any same-length operation, a reordering included. `is_unit_trimmed()`
answers from those positions and names the wrong units. See
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
and [psw](https://r-causal.github.io/propensity/reference/psw.md) for
the whole contract.

## Usage

``` r
is_unit_trimmed(x)
```

## Arguments

- x:

  A `ps_trim` object created by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  or a [psw](https://r-causal.github.io/propensity/reference/psw.md)
  vector built from one.

## Value

A logical vector the same length as `x`, where `TRUE` marks a trimmed
unit.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for trimming propensity scores,
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)
to test whether an object has been trimmed,
[`ps_trim_meta()`](https://r-causal.github.io/propensity/reference/ps_trim_meta.md)
to retrieve full trimming metadata.

## Examples

``` r
ps <- c(0.05, 0.3, 0.6, 0.95)
trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)

is_unit_trimmed(trimmed)
#> [1]  TRUE FALSE FALSE  TRUE

# Use to subset data to retained observations
kept <- !is_unit_trimmed(trimmed)
ps[kept]
#> [1] 0.3 0.6
```
