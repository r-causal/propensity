# Identify which units were truncated

`is_unit_truncated()` returns a logical vector indicating which
observations had their propensity scores modified by truncation. Use
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
to test whether an object has been truncated at all.

The answer comes from the truncation record, which is written for a
fixed set of observations and can both be lost and outlive them. On a
`ps_trunc`,
[`vctrs::vec_slice()`](https://vctrs.r-lib.org/reference/vec_slice.html)
and [`c()`](https://rdrr.io/r/base/c.html) drop it, and subassignment
that grows the vector carries it across a length change; see
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
for the whole contract. On a
[psw](https://r-causal.github.io/propensity/reference/psw.md) vector
built from truncated propensity scores, a subset drops it, while
subassignment that grows the weights carries it across the length
change.

`is_unit_truncated()` therefore checks that the record covers the object
it is given, and raises an error of class
`propensity_missing_meta_error` when it does not, or when an object
marked as truncated carries no record at all, rather than name truncated
units at stale positions. Query the `ps_trunc` object the record was
written for instead.

That check counts observations, which a reordering does not change, so
it does not catch one. An operation that reorders through vctrs rather
than through `[`, such as `vctrs::vec_slice(x, 5:1)` or
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
keeps a record written for the old order, and a `psw` keeps one through
any same-length operation, a reordering included. `is_unit_truncated()`
answers from those positions and names the wrong units. See
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
and [psw](https://r-causal.github.io/propensity/reference/psw.md) for
the whole contract.

## Usage

``` r
is_unit_truncated(x)
```

## Arguments

- x:

  A `ps_trunc` object created by
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  or a [psw](https://r-causal.github.io/propensity/reference/psw.md)
  vector built from one.

## Value

A logical vector the same length as `x` (or number of rows for matrix
input). `TRUE` marks observations whose values were winsorized.

## See also

[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md),
[`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)

## Examples

``` r
ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
is_unit_truncated(ps_t)
#> [1]  TRUE FALSE FALSE FALSE  TRUE
```
