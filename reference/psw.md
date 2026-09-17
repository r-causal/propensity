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
  wt_truncated = FALSE,
  ...
)

psw(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE,
  stabilization_score = NULL,
  wt_truncated = FALSE
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

- wt_truncated:

  Logical. Were the weights themselves truncated? This is separate from
  `truncated`, which describes the propensity scores the weights were
  built from. Defaults to `FALSE`.

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

- [`is_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md)
  returns `TRUE` if the weights themselves were truncated, which is
  recorded apart from whether the propensity scores they were built from
  were.

- `stabilization_score()` returns the user-supplied stabilization score,
  or `NULL` when none was recorded or when a per-observation score was
  dropped because an operation changed the length of the weights.

- [`numerator_model()`](https://r-causal.github.io/propensity/reference/numerator_model.md)
  returns the fitted model supplied to `stabilize` that estimated the
  numerator of the weights, or `NULL` when no model did.

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

### The exposure records

Weights built by
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and its relatives also record the exposure they were built for:
`exposure_type` names its type, and `density_meta` describes the ratio
of densities that weights for a continuous exposure are. Both are read
back with
[`exposure_type()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
and
[`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md),
and both describe the exposure rather than the units, so neither carries
a length of its own.

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
both operands are, and it is marked as trimmed, truncated, calibrated,
or weight-truncated when either operand is. The remaining attributes,
the `stabilization_score`, the records left by a modified propensity
score, the attributes describing a categorical exposure, and the
exposure records, are carried by agreement: one only a single operand
records carries, and one both record with the same value carries. One
they record differently is dropped, since neither value describes the
result, and a warning of class `propensity_metadata_conflict_warning`
names it, once for each attribute dropped that way and whatever order
the inputs were given in. The rule is applied per operation, so an
attribute one operation drops for a disagreement can be carried again by
a later operation whose operands agree. Two density records agree when
they name the same family with the same parameters and the same
numerator and residual spread, so weights built the same way in two
calls agree even though each call builds its own copy of the function
that evaluates the density. A numerator estimated by a fitted model is
read the same way: two records agree when their models write the same
formula and were fit to the same coefficients, and a model is not
compared against the marginal density or against no numerator at all,
both of which are disagreements.

A `density_meta` record describes a continuous exposure, so a result
that drops `exposure_type` for a disagreement drops the density record
with it, without a second warning: the record has nothing left to be
about.

A `stabilization_score` is carried only when the result is stabilized,
which takes both operands. A result that is not stabilized drops the
score without comment.

Combining `psw` objects with [`c()`](https://rdrr.io/r/base/c.html)
preserves the class only when all metadata matches; mismatched metadata
produces a warning and falls back to a plain numeric vector. The
trimming, truncation, and weight truncation records are compared by what
they say about the modification: weights built from scores trimmed or
truncated differently, from a trim that was refit and one that was not,
or truncated at different bounds describe different estimands and are
combined only as numbers. Weights flagged as modified with no record
agree with any record. Concatenation appends one set of observations to
another, so the positions a record names would describe units from the
other input; the result keeps each record without its positions, which
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md),
the printed footer, and a later combine still read. The calibration
record (`ps_calib_meta`) names the curve the scores were calibrated with
rather than any position, the categorical attributes name exposure
levels, and the exposure records describe the exposure rather than any
unit, so all of them carry by the agreement rule above: two calibration
records agree when they name the same method and smoothing.

[`c()`](https://rdrr.io/r/base/c.html) of a single `psw` returns it
unchanged, with every record. A combine through vctrs drops the
positional records even when it is handed a single input:
`vctrs::vec_c(x)`, `dplyr::bind_rows(df)`, `vctrs::vec_rbind(df)`, and
an ungrouped
[`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
all rebuild the column, so a later
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md),
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md),
or
[`is_unit_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md)
on the result refuses it.

Subsetting with `[` preserves class and attributes for vector
subscripts. Two kinds of attribute hold one value per observation. The
records left by a trimmed or truncated propensity score (`ps_trim_meta`
and `ps_trunc_meta`) or by truncating the weights themselves
(`psw_trunc_meta`) name units by position, and `[` is handed the
subscript, so it re-indexes each record onto the result:
[`rev()`](https://rdrr.io/r/base/rev.html),
[`sort()`](https://rdrr.io/r/base/sort.html), `x[order(x)]`, and a
shorter subset all return records naming the units at their new
positions. Any other operation through vctrs is not handed a subscript,
so it keeps the positions only where every unit stays at its place, as
in elementwise arithmetic, and otherwise drops them, at any length:
[`vctrs::vec_slice()`](https://vctrs.r-lib.org/reference/vec_slice.html),
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`unique()`](https://rdrr.io/r/base/unique.html),
[`rep_len()`](https://rdrr.io/r/base/rep.html), and
[`vctrs::vec_c()`](https://vctrs.r-lib.org/reference/vec_c.html) of a
single input all return weights whose records name no units. A dropped
record keeps its method, its bounds, and whether the model was refit, so
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md),
the printed footer, and the bound comparison of a later combine still
read it. Subassignment with `[<-` and `is.na<-` moves no unit and keeps
the records whole; `[<-` refuses a value whose records describe a
different modification, and accepts one with no record. A `[<-` of a
value modified the same way that writes a missing weight leaves that
unit listed as retained in the record, as `is.na<-` on a `psw` does;
`is.na<-` on a `ps_trim` score instead removes the unit from its
retained set. vctrs' own assignment,
[`vctrs::vec_assign()`](https://vctrs.r-lib.org/reference/vec_slice.html),
reaches the same restore a slice does and cannot be told apart from one,
so it drops the positions, and so do the helpers built on it or on a
combine, such as
[`tidyr::replace_na()`](https://tidyr.tidyverse.org/reference/replace_na.html),
[`dplyr::coalesce()`](https://dplyr.tidyverse.org/reference/coalesce.html),
[`dplyr::if_else()`](https://dplyr.tidyverse.org/reference/if_else.html),
and
[`dplyr::case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html).

A `stabilization_score` with more than one value, on the weights or on a
component of a
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
product, is in the order of the units too. `[` subsets it with the
weights, elementwise arithmetic keeps it, and any other route through
vctrs that brings back observations drops it with a warning of class
`propensity_stabilization_score_warning`, since the score may no longer
line up with the weights.

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
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md),
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md),
and
[`is_unit_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md)
answer by position, so each checks that the record covers the vector it
is given and raises an error of class `propensity_missing_meta_error`
when it does not, or when weights marked as modified carry no record at
all, rather than name modified units at stale positions.
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
reads a single flag rather than a position, so it answers from any
record present and refuses only when the record is absent entirely.

The result of any of these operations stays a `psw` and keeps every
other attribute, including its stabilized, trimmed, truncated, and
calibrated status, the calibration record, which names a calibration
curve rather than any unit, the attributes describing a categorical
exposure, which name the exposure levels rather than the units, and the
exposure records, which describe the exposure rather than any unit, so
all of them mean the same thing at any length.

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
