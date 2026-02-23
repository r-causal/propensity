# Extract trimming metadata from a `ps_trim` object

`ps_trim_meta()` returns the metadata list attached to a `ps_trim`
object by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md).

## Usage

``` r
ps_trim_meta(x)
```

## Arguments

- x:

  A `ps_trim` object.

## Value

A named list with elements:

- `method`:

  Character string indicating the trimming method used.

- `keep_idx`:

  Integer vector of retained observation indices.

- `trimmed_idx`:

  Integer vector of trimmed observation indices.

- `lower`, `upper`:

  Numeric cutoffs, when applicable.

- `refit`:

  Logical, `TRUE` if the model was refit via
  [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md).

Additional method-specific elements (e.g. `cutoff`, `delta`, `lambda`)
may also be present.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for trimming propensity scores,
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)
and
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md)
for predicate queries.

## Examples

``` r
ps <- c(0.05, 0.3, 0.6, 0.95)
trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)

ps_trim_meta(trimmed)
#> $method
#> [1] "ps"
#> 
#> $lower
#> [1] 0.1
#> 
#> $upper
#> [1] 0.9
#> 
#> $keep_idx
#> [1] 2 3
#> 
#> $trimmed_idx
#> [1] 1 4
#> 
```
