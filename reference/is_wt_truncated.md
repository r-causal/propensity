# Test whether the weights themselves have been truncated

`is_wt_truncated()` returns `TRUE` for a
[psw](https://r-causal.github.io/propensity/reference/psw.md) vector
whose own values were bounded, and `is_unit_wt_truncated()` returns
which of its weights the bound moved.

Truncating the weights is recorded apart from truncating the propensity
scores they were built from, which
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
and
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md)
report. The two are different operations, and a set of weights can carry
either, both, or neither.

## Usage

``` r
is_wt_truncated(x)

is_unit_wt_truncated(x)
```

## Arguments

- x:

  An object. `is_unit_wt_truncated()` accepts only a
  [psw](https://r-causal.github.io/propensity/reference/psw.md) vector.

## Value

- `is_wt_truncated()`: a single `TRUE` or `FALSE`. It is `FALSE` for
  anything that is not a
  [psw](https://r-causal.github.io/propensity/reference/psw.md) vector.

- `is_unit_wt_truncated()`: a logical vector the same length as `x`,
  `TRUE` for each weight the bound moved.

## Details

The truncation leaves a record on the weights, the `psw_trunc_meta`
attribute, holding the method, the bound as given and as applied, and
the positions of the weights it moved. `is_unit_wt_truncated()` answers
from those positions, so the record follows the rules
[psw](https://r-causal.github.io/propensity/reference/psw.md) describes
for the records a modified propensity score leaves: it is kept through
arithmetic and subassignment, re-indexed through the subscript of `[`,
and loses its positions to any other slice and to a combine, which keeps
the method and bounds. The `wt_truncated` flag describes the weights as
a whole and is kept through all of those, so `is_wt_truncated()` keeps
its answer where `is_unit_wt_truncated()` has none to give.

`is_unit_wt_truncated()` therefore checks that the record covers the
weights it is given, and raises an error of class
`propensity_missing_meta_error` when it does not, or when weights marked
as truncated carry no record at all, rather than name truncated weights
at stale positions. Weights that were not truncated have no record to
read, and every weight is reported as untouched.

## See also

[psw](https://r-causal.github.io/propensity/reference/psw.md) for the
weight vector class, and
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
for truncation of the propensity scores.

## Examples

``` r
w <- psw(c(1.2, 0.8, 1.5), estimand = "ate")
is_wt_truncated(w)
#> [1] FALSE
is_unit_wt_truncated(w)
#> [1] FALSE FALSE FALSE
```
