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
  stabilization_score = NULL,
  ...
)

psw(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE,
  stabilization_score = NULL
)

is_psw(x)

is_stabilized(wt)

stabilization_score(wt)

as_psw(x, estimand = NULL)
```

## Arguments

- x:

  For `psw()` and `new_psw()`: a numeric vector of weights (default:
  [`double()`](https://rdrr.io/r/base/double.html)). For `is_psw()` and
  `as_psw()`: an object to test or coerce.

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

- stabilization_score:

  Optional numeric stabilization score to record on the object, either a
  single value or one value per observation. Every value must be
  positive and finite. Defaults to `NULL`, meaning no fixed score was
  supplied.

- ...:

  Additional attributes stored on the object (developer use only).

- wt:

  A `psw` or `causal_wts` object.

## Value

- `new_psw()`, `psw()`, `as_psw()`: A `psw` vector.

- `is_psw()`, `is_stabilized()`: A single logical value.

- `stabilization_score()`: A numeric value or vector, or `NULL` if none
  was recorded or a per-observation score was dropped.

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

- `is_stabilized()` returns `TRUE` if the weights are stabilized.

- `stabilization_score()` returns the user-supplied stabilization score,
  or `NULL` when none was recorded or when a per-observation score was
  dropped because an operation changed the length of the weights.

### The stabilization score

A `stabilization_score` is the multiplier the weights were stabilized
on. It holds either a single value, which scales every weight, or one
value per observation, which scales each weight by its own. Both forms
are checked where the score is recorded: a score must be numeric,
positive, and finite, and hold a length the weights can use, and one
that does not is refused with an error of class
`propensity_stabilization_score_error`. A score is recorded as a plain
double vector, with its storage type normalized and any names dropped,
so a score written `1L` and one written `1` are the same score and
combine without conflict.

A zero-length `psw` is a prototype: it records what a result built from
it will carry rather than describing observations of its own, so a
per-observation score on one is recorded as given and checked for length
when the observations arrive. A cast whose data arrives at a length the
score does not match drops the score rather than refusing the cast,
since the score describes the prototype's observations and says nothing
about the data being cast. A cast to no observations keeps it, having
brought no observations for it to contradict.

`psw` objects also inherit the broader `causal_wts` class. The accessors
that class carries,
[`causalgenerics::is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html),
[`causalgenerics::estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html),
and `estimand<-`, are re-exported by propensity and documented at
[`causalgenerics::estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html).

### Arithmetic and combining

Arithmetic operations on `psw` objects preserve the class and
attributes, so operations like normalization (`weights / sum(weights)`)
retain metadata.

An operation between two `psw` objects merges what each of them records.
Two different estimands are pasted together, and an estimand only one
operand names stands for the result; the result is stabilized only when
both operands are, and it is marked as trimmed, truncated, or calibrated
when either operand is. The remaining attributes, the
`stabilization_score`, the records left by a modified propensity score,
and the attributes describing a categorical exposure, are carried by
agreement: one only a single operand records carries, and one both
record with the same value carries. One they record differently is
dropped, since neither value describes the result, and a warning of
class `propensity_metadata_conflict_warning` names it, once for each
attribute dropped that way and whatever order the inputs were given in.
The rule is applied per operation, so an attribute one operation drops
for a disagreement can be carried again by a later operation whose
operands agree.

A `stabilization_score` is carried only when the result is stabilized,
which takes both operands. A result that is not stabilized drops the
score without comment.

Combining `psw` objects with [`c()`](https://rdrr.io/r/base/c.html)
preserves the class only when all metadata matches; mismatched metadata
produces a warning and falls back to a plain numeric vector.
Concatenation appends one set of observations to another, so the
positions a modification record names would describe units from the
other input; those records are dropped from the result whether or not
the inputs agree on them. The categorical attributes name exposure
levels rather than positions and carry when the inputs agree.

Subsetting with `[` preserves class and attributes for vector
subscripts. Two kinds of attribute hold one value per observation and so
cannot be re-indexed for a subset: a `stabilization_score` with more
than one value, and the records left by a modified propensity score
(`ps_trim_meta`, `ps_trunc_meta`, and `ps_calib_meta`). Where an
operation goes through vctrs, these are carried when the result comes
back at the length they were recorded on and dropped when it does not.
Any same-length operation keeps them, a reordering or a subscript with
duplicates included, so the positions a record names can end up
describing different observations than they did.

Dropping the `stabilization_score` warns, because the score was supplied
by the user and the weights can be recomputed on the subset. Dropping a
modification record is silent, because these records also travel by
routes vctrs does not see, and a warning on the one route it controls
would be neither complete nor about anything the user wrote.
Subassignment that grows the vector carries a record across the length
change, since `[<-` casts the replacement and then leaves base R to
preserve the attributes; and
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) drops the
`NA`-weighted rows from a weights column in C and re-attaches the
original variable's attributes to the shortened result, so weights built
on trimmed propensity scores come back out of every outcome model fit on
them still carrying a record written for rows that are no longer there.

Honesty therefore lives at query time.
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md)
answers by position, so it checks that the record covers the vector it
is given and raises an error of class `propensity_missing_meta_error`
when it does not, or when weights marked as trimmed carry no record at
all, rather than name trimmed units at stale positions.
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
reads a single flag rather than a position, so it answers from any
record present and refuses only when the record is absent entirely.

The result of any of these operations stays a `psw` and keeps every
other attribute, including its stabilized, trimmed, truncated, and
calibrated status, and the attributes describing a categorical exposure,
which name the exposure levels rather than the units and so mean the
same thing at any length.

Matrix or array subscripts intentionally drop the `psw` class and return
a plain numeric vector via base R linear indexing; this is required so
[`glm.fit()`](https://rdrr.io/r/stats/glm.html)-style internal indexing
works on `psw`-weighted GLMs. Summary functions
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
