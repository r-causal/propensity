# Create and Manipulate `psw` Objects

Functions to create and manipulate `psw` objects, which are specialized
vectors for propensity score weights with optional `estimand`
attributes. Most users should use
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and friends, but these functions can help extend the functionality of
`psw` objects.

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

  A numeric vector (default:
  [`double()`](https://rdrr.io/r/base/double.html)).

- estimand:

  A character string representing the estimand (e.g., "ate", "att",
  "ato"). Default is `NULL`.

- stabilized:

  A logical `TRUE`

- trimmed:

  Logical, whether these weights came from a trimmed PS.

- truncated:

  Logical, whether these weights came from a truncated PS.

- calibrated:

  Logical, whether these weights came from a calibrated PS.

- ...:

  Additional attributes to track in the weights.

- wt:

  An object to check or convert.

- value:

  The value to add to the attribute.

## Value

- `new_psw()`: A `psw` object.

- `psw()`: A `psw` object.

- `is_psw()`: `TRUE` if the object is a `psw`, otherwise `FALSE`.

- `as_psw()`: A `psw` object.

- `estimand()`: The `estimand` attribute of a `psw` object.

- `is_stabilized()`: The `stabilized` attribute of a `psw` object.

## Details

The `psw` class is a vctrs-based S3 class that represents propensity
score weights. It extends numeric vectors with additional metadata
tracking the estimand type, stabilization status, and source
transformations.

**Arithmetic behavior**: Unlike `ps_trim` and `ps_trunc` objects,
arithmetic operations on `psw` objects preserve the class and
attributes. This allows weight manipulations like normalization
(`weights / sum(weights)`) while maintaining metadata.

**Combining behavior**: When combining `psw` objects with
[`c()`](https://rdrr.io/r/base/c.html), the class is preserved only if
all metadata matches. Mismatched metadata triggers a warning and returns
a numeric vector.

**Base R compatibility**: Most base R operations work seamlessly:

- Subsetting with `[` preserves class and attributes

- Summary functions ([`sum()`](https://rdrr.io/r/base/sum.html),
  [`mean()`](https://rdrr.io/r/base/mean.html), etc.) return numeric
  values

- Comparison operators return logical vectors

- Works in data frames and with tidyverse functions

## Examples

``` r
psw_weights <- new_psw(c(0.1, 0.2, 0.3), estimand = "ate")
is_psw(psw_weights)
#> [1] TRUE
estimand(psw_weights)
#> [1] "ate"

psw_helper <- psw(c(0.5, 0.7), estimand = "att")
as_psw(c(0.1, 0.2), estimand = "ato")
#> <psw{estimand = ato}[2]>
#> [1] 0.1 0.2

# Coercion behavior - compatible objects combine silently
x <- psw(c(0.5, 0.7), estimand = "ate")
y <- psw(c(0.3, 0.8), estimand = "ate")
c(x, y)  # Returns psw object
#> <psw{estimand = ate}[4]>
#> [1] 0.5 0.7 0.3 0.8

# Incompatible metadata triggers warning and returns numeric
x <- psw(c(0.5, 0.7), estimand = "ate")
y <- psw(c(0.3, 0.8), estimand = "att")
c(x, y)  # Warning: returns numeric
#> Warning: Converting psw to numeric: incompatible estimands 'ate' and 'att'
#> ℹ Metadata cannot be preserved when combining incompatible objects
#> ℹ Use identical objects or explicitly cast to numeric to avoid this warning
#> [1] 0.5 0.7 0.3 0.8

# Works with tidyr::pivot_longer for plotting
if (requireNamespace("tidyr", quietly = TRUE)) {
  df <- data.frame(
    id = 1:4,
    ate_wts = psw(c(0.5, 0.7, 0.3, 0.8), estimand = "ate"),
    att_wts = psw(c(0.4, 0.6, 0.2, 0.9), estimand = "att")
  )
  # This will warn but succeed, returning numeric in the pivoted column
  tidyr::pivot_longer(df, cols = c(ate_wts, att_wts))
}
#> Warning: Converting psw to numeric: incompatible estimands 'ate' and 'att'
#> ℹ Metadata cannot be preserved when combining incompatible objects
#> ℹ Use identical objects or explicitly cast to numeric to avoid this warning
#> # A tibble: 8 × 3
#>      id name    value
#>   <int> <chr>   <dbl>
#> 1     1 ate_wts   0.5
#> 2     1 att_wts   0.4
#> 3     2 ate_wts   0.7
#> 4     2 att_wts   0.6
#> 5     3 ate_wts   0.3
#> 6     3 att_wts   0.2
#> 7     4 ate_wts   0.8
#> 8     4 att_wts   0.9
```
