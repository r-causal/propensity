# Propensity Score Weight Vectors

`psw` objects are numeric vectors that carry metadata about propensity
score weights, including the target estimand and whether the underlying
propensity scores were trimmed, truncated, or calibrated.

Most users will encounter `psw` objects as return values from
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and related weight functions. These constructor and helper functions are
useful for inspecting weight objects or for package developers extending
propensity.

## Usage

``` r
new_psw(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE,
  ...
)

psw(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE
)

is_psw(x)

is_stabilized(wt)

is_causal_wt(x)

as_psw(x, estimand = NULL)

estimand(wt)

estimand(wt) <- value
```

## Arguments

- x:

  For `psw()` and `new_psw()`: a numeric vector of weights (default:
  [`double()`](https://rdrr.io/r/base/double.html)). For `is_psw()`,
  `is_causal_wt()`, and `as_psw()`: an object to test or coerce.

- estimand:

  A character string identifying the target estimand (e.g., `"ate"`,
  `"att"`, `"ato"`). Defaults to `NULL`.

- stabilized:

  Logical. Were the weights stabilized? Defaults to `FALSE`.

- trimmed:

  Logical. Were the weights derived from trimmed propensity scores?
  Defaults to `FALSE`.

- truncated:

  Logical. Were the weights derived from truncated propensity scores?
  Defaults to `FALSE`.

- calibrated:

  Logical. Were the weights derived from calibrated propensity scores?
  Defaults to `FALSE`.

- ...:

  Additional attributes stored on the object (developer use only).

- wt:

  A `psw` or `causal_wts` object.

- value:

  A character string: the new estimand to assign.

## Value

- `new_psw()`, `psw()`, `as_psw()`: A `psw` vector.

- `is_psw()`, `is_causal_wt()`, `is_stabilized()`: A single logical
  value.

- `estimand()`: A character string, or `NULL` if no estimand is set.

- `estimand<-`: The modified `psw` object (called for its side effect).

## Details

### Constructors

- `psw()` is the **user-facing** constructor. It coerces `x` to double
  and validates inputs before creating the object.

- `new_psw()` is the **low-level** constructor intended for developers.
  It assumes `x` is already a double vector and performs minimal
  validation.

- `as_psw()` coerces an existing numeric vector to a `psw` object.

### Queries

- `is_psw()` tests whether an object is a `psw` vector.

- `is_causal_wt()` tests whether an object inherits from the broader
  `causal_wts` class (which includes `psw` objects).

- `estimand()` and `estimand<-` get and set the estimand attribute.

- `is_stabilized()` returns `TRUE` if the weights are stabilized.

### Arithmetic and combining

Arithmetic operations on `psw` objects preserve the class and
attributes, so operations like normalization (`weights / sum(weights)`)
retain metadata. Combining `psw` objects with
[`c()`](https://rdrr.io/r/base/c.html) preserves the class only when all
metadata matches; mismatched metadata produces a warning and falls back
to a plain numeric vector.

Subsetting with `[` preserves class and attributes. Summary functions
([`sum()`](https://rdrr.io/r/base/sum.html),
[`mean()`](https://rdrr.io/r/base/mean.html), etc.) return plain numeric
values.

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for calculating propensity score weights (which return `psw` objects).

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
and
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
for modifying propensity scores before weight calculation.

## Examples

``` r
# Create psw objects directly
w <- psw(c(1.2, 0.8, 1.5), estimand = "ate")
w
#> <psw{estimand = ate}[3]>
#> [1] 1.2 0.8 1.5

# Query metadata
is_psw(w)
#> [1] TRUE
estimand(w)
#> [1] "ate"
is_stabilized(w)
#> [1] FALSE

# Coerce a plain numeric vector
as_psw(c(1.0, 2.0), estimand = "att")
#> <psw{estimand = att}[2]>
#> [1] 1 2

# Arithmetic preserves the psw class
w / sum(w)
#> <psw{estimand = ate}[3]>
#> [1] 0.3428571 0.2285714 0.4285714

# Combining: compatible metadata is preserved
x <- psw(c(1.2, 0.8), estimand = "ate")
y <- psw(c(1.1, 0.9), estimand = "ate")
c(x, y)
#> <psw{estimand = ate}[4]>
#> [1] 1.2 0.8 1.1 0.9

# Combining: incompatible metadata warns and returns numeric
x <- psw(c(1.2, 0.8), estimand = "ate")
y <- psw(c(1.1, 0.9), estimand = "att")
c(x, y)
#> Warning: Converting psw to numeric: incompatible estimands 'ate' and 'att'
#> ℹ Metadata cannot be preserved when combining incompatible objects
#> ℹ Use identical objects or explicitly cast to numeric to avoid this warning
#> [1] 1.2 0.8 1.1 0.9
```
