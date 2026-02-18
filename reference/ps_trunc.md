# Truncate (Winsorize) Propensity Scores

**`ps_trunc()`** sets out‐of‐range propensity scores to fixed bounding
values (a form of *winsorizing*). This is an alternative to
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
which removes (sets `NA`) instead of bounding and is then refit with
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)

## Usage

``` r
ps_trunc(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)
```

## Arguments

- ps:

  The propensity score, either a numeric vector between 0 and 1 for
  binary exposures, or a matrix/data.frame where each column represents
  propensity scores for each level of a categorical exposure.

- method:

  One of `"ps"`, `"pctl"`, or `"cr"`.

  - `"ps"`: directly cut on `[lower, upper]` of `ps`. For categorical,
    uses symmetric truncation with `lower` as the threshold.

  - `"pctl"`: use quantiles of `ps` as bounding values. For categorical,
    calculates quantiles across all propensity score values.

  - `"cr"`: the common range of `ps` given `.exposure`, bounding
    `[min(ps[treated]), max(ps[untreated])]` (binary only)

- lower, upper:

  Numeric or quantile bounds. If `NULL`, defaults vary by method. For
  categorical exposures with method `"ps"`, `lower` represents the
  truncation threshold (delta).

- .exposure:

  For method "cr", a binary exposure vector. For categorical exposures,
  must be a factor or character vector.

- .focal_level:

  For binary exposures, the value representing the focal group
  (typically the treatment group). For categorical exposures with ATT or
  ATU estimands, specifies the focal category. Must be one of the levels
  of the exposure variable. Required for
  [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  and
  [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  with categorical exposures.

- .reference_level:

  For binary exposures, the value representing the reference group
  (typically the control group). If not provided, it is automatically
  detected.

- ...:

  Additional arguments passed to methods

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

## Value

A **`ps_trunc`** object (numeric vector or matrix). It has an attribute
`ps_trunc_meta` storing fields like `method`, `lower_bound`, and
`upper_bound`.

## Details

For binary exposures with each \\ps\[i\]\\:

- If \\ps\[i\] \< lower\\bound\\, we set \\ps\[i\] = lower\\bound\\.

- If \\ps\[i\] \> upper\\bound\\, we set \\ps\[i\] = upper\\bound\\.

For categorical exposures:

- Each value below the threshold is set to the threshold

- Rows are renormalized to sum to 1

This approach is often called *winsorizing*.

**Arithmetic behavior**: Like `ps_trim`, arithmetic operations on
`ps_trunc` objects return numeric vectors. The reasoning is the same -
transformed propensity scores (e.g., weights) are no longer propensity
scores.

**No NA values**: Unlike `ps_trim`, truncation doesn't create `NA`
values. Out-of-range values are set to the boundaries, so all values
remain finite and valid for calculations.

**Metadata tracking**: The `truncated_idx` tracks which positions had
their values modified (winsorized to boundaries):

- Subsetting with `[` updates indices to new positions

- [`sort()`](https://rdrr.io/r/base/sort.html) reorders data and updates
  indices accordingly

- Operations preserve finite values (no `NA` handling needed)

**Boundary detection**: Values exactly at the boundaries (after
truncation) may indicate truncation, but aren't necessarily truncated -
they could have been at the boundary originally.

**Combining behavior**: When combining `ps_trunc` objects with
[`c()`](https://rdrr.io/r/base/c.html), truncation parameters must
match. Mismatched parameters trigger a warning and return a numeric
vector.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
and
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
for removing extreme values vs. bounding

## Examples

``` r
set.seed(2)
n <- 30
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x))
fit <- glm(z ~ x, family = binomial)
ps <- predict(fit, type = "response")

# truncate just the 99th percentile
ps_trunc(ps, method = "pctl", lower = 0, upper = .99)
#> <ps_trunc{[0.341443426776033,0.805793268892769], method=pctl}[30]>
#>         1         2         3         4         5         6         7         8 
#> 0.5149714 0.6361298 0.7694837 0.4880712 0.6073989 0.6305169 0.6899234 0.5897388 
#>         9        10        11        12        13        14        15        16 
#> 0.8003122 0.6009455 0.6605909 0.7162599 0.5725720 0.4985231 0.7849940 0.3561684 
#>        17        18        19        20        21        22        23        24 
#> 0.7064972 0.6200818 0.7191624 0.6620999 0.8057933 0.4800637 0.7696302 0.7981060 
#>        25        26        27        28        29        30 
#> 0.6167236 0.3414434 0.6667225 0.5494305 0.6981704 0.6472363 

# Coercion behavior with ps_trunc objects
ps_trunc1 <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
ps_trunc2 <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

# Compatible objects combine silently
c(ps_trunc1[1:15], ps_trunc2[16:30])  # Returns ps_trunc object
#> <ps_trunc{[0.1,0.9], method=ps}[30]>
#>         1         2         3         4         5         6         7         8 
#> 0.5149714 0.6361298 0.7694837 0.4880712 0.6073989 0.6305169 0.6899234 0.5897388 
#>         9        10        11        12        13        14        15        16 
#> 0.8003122 0.6009455 0.6605909 0.7162599 0.5725720 0.4985231 0.7849940 0.3561684 
#>        17        18        19        20        21        22        23        24 
#> 0.7064972 0.6200818 0.7191624 0.6620999 0.8080320 0.4800637 0.7696302 0.7981060 
#>        25        26        27        28        29        30 
#> 0.6167236 0.3414434 0.6667225 0.5494305 0.6981704 0.6472363 

# Different truncation parameters trigger warning
ps_trunc3 <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
c(ps_trunc1[1:15], ps_trunc3[16:30])  # Warning: returns numeric
#> Warning: Converting ps_trunc to numeric: different truncation parameters
#> ℹ Metadata cannot be preserved when combining incompatible objects
#> ℹ Use identical objects or explicitly cast to numeric to avoid this warning
#>         1         2         3         4         5         6         7         8 
#> 0.5149714 0.6361298 0.7694837 0.4880712 0.6073989 0.6305169 0.6899234 0.5897388 
#>         9        10        11        12        13        14        15        16 
#> 0.8003122 0.6009455 0.6605909 0.7162599 0.5725720 0.4985231 0.7849940 0.3561684 
#>        17        18        19        20        21        22        23        24 
#> 0.7064972 0.6200818 0.7191624 0.6620999 0.8000000 0.4800637 0.7696302 0.7981060 
#>        25        26        27        28        29        30 
#> 0.6167236 0.3414434 0.6667225 0.5494305 0.6981704 0.6472363 

# Mixing with other propensity classes warns
ps_trim_obj <- ps_trim(ps[1:15], method = "ps", lower = 0.1)
c(ps_trunc1[1:15], ps_trim_obj)  # Warning: returns numeric
#> Warning: Converting ps_trunc and ps_trim to numeric
#> ℹ Class-specific attributes and metadata have been dropped
#> ℹ Use explicit casting to numeric to avoid this warning
#>         1         2         3         4         5         6         7         8 
#> 0.5149714 0.6361298 0.7694837 0.4880712 0.6073989 0.6305169 0.6899234 0.5897388 
#>         9        10        11        12        13        14        15         1 
#> 0.8003122 0.6009455 0.6605909 0.7162599 0.5725720 0.4985231 0.7849940 0.5149714 
#>         2         3         4         5         6         7         8         9 
#> 0.6361298 0.7694837 0.4880712 0.6073989 0.6305169 0.6899234 0.5897388 0.8003122 
#>        10        11        12        13        14        15 
#> 0.6009455 0.6605909 0.7162599 0.5725720 0.4985231 0.7849940 
```
