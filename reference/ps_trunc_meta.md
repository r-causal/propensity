# Extract truncation metadata from a `ps_trunc` object

Returns the metadata list attached to a
[`ps_trunc`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
object. The list includes fields such as `method`, `lower_bound`,
`upper_bound`, and `truncated_idx`.

## Usage

``` r
ps_trunc_meta(x)
```

## Arguments

- x:

  A `ps_trunc` object created by
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md).

## Value

A named list with truncation metadata, including:

- `method` – the truncation method used (`"ps"`, `"pctl"`, or `"cr"`)

- `lower_bound`, `upper_bound` – the applied bounds

- `truncated_idx` – integer positions of values that were winsorized

- `n_obs` – the number of observations those positions describe

`truncated_idx` and `n_obs` are absent when an operation that changed
the number of observations dropped the record; see
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md).

## See also

[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md),
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md)

## Examples

``` r
ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
ps_trunc_meta(ps_t)
#> $method
#> [1] "ps"
#> 
#> $lower_bound
#> [1] 0.05
#> 
#> $upper_bound
#> [1] 0.95
#> 
#> $truncated_idx
#> [1] 1 5
#> 
#> $n_obs
#> [1] 5
#> 
```
