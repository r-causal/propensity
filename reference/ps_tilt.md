# Propensity score tilting functions

Every estimand this package supports targets a population whose
covariate distribution is the study population reweighted by a tilting
function \\h\\ of the propensity score. `ps_tilt()` evaluates \\h\\ at
each observation's propensity score.

The tilt is the numerator of every weight: a weight is \\h\\ divided by
the propensity score of the exposure level the unit actually received.
It is also what standardizes a g-computation estimate to a target
population, which for `"atm"`, `"ato"`, and `"entropy"` is the only
route to the estimand, since those populations are not subsets of the
rows and cannot be reached by filtering.

## Usage

``` r
ps_tilt(
  .propensity,
  estimand,
  ...,
  .focal_level = NULL,
  ps = lifecycle::deprecated()
)
```

## Arguments

- .propensity:

  Propensity scores. A numeric vector of \\P(Z = \text{focal} \mid X)\\
  for a binary exposure, or a matrix or data frame with one column per
  exposure level, named for that level, for a categorical exposure. A
  data frame holding a `.pred_class` column, which a fitted tidymodels
  classification model returns when no prediction type is named, carries
  predicted levels rather than probabilities and is refused with an
  error of class `propensity_df_class_column_error`.

- estimand:

  One of `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`.

- ...:

  These dots are for future extensions and must be empty.

- .focal_level:

  The exposure level the `"att"` and `"atu"` tilts target, matched
  against the column names of `.propensity`. A column named for the
  level exactly wins; failing that, a `.pred_` prefix is stripped, so
  `.pred_a` matches the level `a`. Required for those two estimands with
  a categorical `.propensity`, and accepted nowhere else: a numeric
  `.propensity` is already the probability of the focal level, and the
  remaining tilts treat every level alike.

- ps:

  **\[deprecated\]** Use `.propensity` instead. A call that names `ps`
  must name the arguments after it as well, since a positional argument
  binds to `.propensity`.

## Value

A plain double vector, unnamed, with one element per observation: the
length of `.propensity` for the numeric method, and the number of rows
of `.propensity` for the matrix and data frame methods.

## Tilting functions

For a binary exposure with propensity score \\e = P(Z = \text{focal}
\mid X)\\:

|  |  |  |
|----|----|----|
| estimand | \\h(e)\\ | target population |
| `"ate"` | \\1\\ | everyone |
| `"att"` | \\e\\ | the focal group |
| `"atu"` | \\1 - e\\ | the reference group |
| `"atm"` | \\\min(e, 1 - e)\\ | the evenly matchable |
| `"ato"` | \\e(1 - e)\\ | the overlap population |
| `"entropy"` | \\-e \log(e) - (1 - e) \log(1 - e)\\ | the entropy-tilted population |

For a categorical exposure with propensity score vector \\(e_1, \ldots,
e_K)\\ and focal level \\f\\: `"ate"` is \\1\\, `"att"` is \\e_f\\,
`"atu"` is \\1 - e_f\\, `"atm"` is \\\min_k e_k\\, `"ato"` is \\(\sum_k
1 / e_k)^{-1}\\, and `"entropy"` is \\-\sum_k e_k \log(e_k)\\.

Stabilization and censoring weights are not tilts. Stabilization
multiplies a weight by a marginal quantity that does not depend on the
covariates, and
[`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
reuses the `"ate"` formula, so neither has an entry here.

## Standardizing model predictions

A g-computation estimate averages an outcome model's per-row
counterfactual predictions, and where those predictions vary from row to
row it is the weights of that average that standardize the estimate to a
target population. The tilt is those weights, which is what the second
example below supplies by hand. An average taken with equal weight
standardizes a covariate-adjusted outcome model to the whole sample, so
such a model targets the ATE however its own fitting weights were built:
fit for a non-`"ate"` estimand and averaged that way, it reports the
full-sample contrast rather than the one it was weighted for. Tools that
average per-row model predictions, such as the `avg_comparisons()`
function in the marginaleffects package, average with equal weight
unless they are told otherwise, and the remedy is to hand them the tilt
as the weights of the average: `ps_tilt(ps, "att")` for a binary
exposure and `ps_tilt(ps, "att", .focal_level = "b")` for a categorical
one.

The requirement is easy to miss because an outcome model saturated in
the exposure hides it. Such a model predicts one value per exposure
level, so every average of its predictions returns the same contrast,
and there the estimand is settled by the weights the model was fit with
rather than by the weights of the average: a model fit with
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
weights reports the ATT whether the average is taken with equal weight
or with the tilt. The two averages come apart at the first covariate the
outcome model adjusts for, where the predictions vary from row to row
and the weights the average is taken under are what decides the
population the estimate describes.

## Propensity score range

Every propensity score in `.propensity` must lie strictly inside \\(0,
1)\\, the same requirement
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and the rest of the weight family impose before they tilt. The bound
holds for a matrix or data frame `.propensity` entry by entry, and each
row must sum to one on top of it. A fitted model that separates the
exposure can return a probability of exactly zero or one; those scores
have no weight to divide and are rejected here rather than tilted.

A missing propensity score gives a missing tilt under every estimand,
`"ate"` included, so an observation whose propensity score is unknown
never counts toward a tilted mean. A numeric `.propensity` propagates
`NA` position by position, and a matrix `.propensity` gives `NA` for any
row holding one: a probability vector with a missing entry is not one
this can tilt on, whichever level the tilt reads.

## Fitted models

`.propensity` can be the fitted propensity score model instead of the
scores it reports. A binomial
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and a two-level
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) report
one score per unit and are tilted as a vector; a
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) of
three or more levels reports one column per level and is tilted as a
matrix, with `.focal_level` naming the column `"att"` and `"atu"` read.
Those are the shapes `predict(fit, type = "response")` and `fitted(fit)`
give, and tilting a fit tilts exactly what tilting those values would. A
tilt reads no exposure, so nothing but the scores is taken off the
model.

## Modified propensity scores

`ps_tilt()` takes plain propensity scores. A score modified by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
carries a class of its own and has no method here; pass the scores
underneath it, with `as.numeric(x)` for a binary exposure or
`as.matrix(x)` for a categorical one. Units
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
set to `NA` stay `NA` through that extraction and take an `NA` tilt.

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and the rest of the weight family, which divide the tilt by the
propensity score of the received exposure level.

## Examples

``` r
set.seed(1)
n <- 500
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.6 * x))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + x))
sim <- data.frame(x = x, z = z, y = y)
ps <- unname(
  predict(glm(z ~ x, data = sim, family = binomial()), type = "response")
)

# a weight is the tilt over the propensity of the received exposure level
received <- z * ps + (1 - z) * (1 - ps)
all.equal(
  as.numeric(wt_ato(ps, z, exposure_type = "binary")),
  ps_tilt(ps, "ato") / received
)
#> [1] TRUE

# tilted g-computation: standardize counterfactual predictions to the
# overlap population with an h-weighted mean, no weights in the outcome model
fit <- glm(y ~ z * x, data = sim, family = binomial())
m1 <- predict(fit, transform(sim, z = 1), type = "response")
m0 <- predict(fit, transform(sim, z = 0), type = "response")

h <- ps_tilt(ps, "ato")
weighted.mean(m1, h) - weighted.mean(m0, h)
#> [1] 0.1231355

# the same predictions standardized to everyone give the ATE instead
mean(m1) - mean(m0)
#> [1] 0.1210093

# the tilt of the scores a fitted model reports
ps_tilt(glm(z ~ x, data = sim, family = binomial()), "ato")
#>   [1] 0.23534442 0.24996155 0.22688346 0.20297217 0.24912843 0.22754947
#>   [7] 0.24701618 0.24118477 0.24530010 0.24487234 0.20761371 0.24846825
#>  [13] 0.23553414 0.14807733 0.22690699 0.24907816 0.24933571 0.23432405
#>  [19] 0.23863039 0.24490199 0.23524615 0.23987179 0.24987314 0.16185673
#>  [25] 0.24430551 0.24896655 0.24769579 0.19342860 0.24030505 0.24809853
#>  [31] 0.21573118 0.24843344 0.24849519 0.24899024 0.19888986 0.24213142
#>  [37] 0.24269139 0.24893364 0.22799655 0.24045032 0.24756098 0.24597673
#>  [43] 0.24234645 0.24570324 0.23299457 0.23225857 0.24876694 0.24028844
#>  [49] 0.24831076 0.23660357 0.24836362 0.23586681 0.24901562 0.21258879
#>  [55] 0.21185829 0.18022205 0.24339409 0.21698722 0.24542981 0.24800105
#>  [61] 0.15434193 0.24913250 0.24254146 0.24964987 0.23081698 0.24995065
#>  [67] 0.17321384 0.21012183 0.24999822 0.16842293 0.24721872 0.23216132
#>  [73] 0.24451846 0.22237432 0.20586842 0.24945033 0.24133495 0.24947049
#>  [79] 0.24987234 0.23666496 0.23738605 0.24799928 0.22451235 0.19029911
#>  [85] 0.24490098 0.24909569 0.22957253 0.24489936 0.24870536 0.24961765
#>  [91] 0.23826488 0.22313096 0.22531928 0.24225794 0.20344781 0.24566555
#>  [97] 0.20459136 0.23722863 0.20746765 0.24044851 0.23556584 0.24972924
#> [103] 0.22346335 0.24999559 0.23430212 0.19303347 0.24180145 0.23556683
#> [109] 0.24853793 0.19800785 0.23500393 0.24079937 0.21189758 0.23444804
#> [115] 0.24684468 0.24273073 0.24453928 0.24544619 0.24689823 0.24735637
#> [121] 0.23944548 0.21652810 0.24671556 0.24731998 0.24846598 0.24191443
#> [127] 0.24878014 0.24914753 0.23326976 0.24443984 0.24981591 0.23668688
#> [133] 0.24620728 0.19060706 0.24933132 0.18953087 0.24497098 0.23873143
#> [139] 0.23439563 0.24895866 0.16647657 0.22458135 0.18178069 0.24074350
#> [145] 0.21329468 0.23050693 0.17368270 0.24958353 0.20404829 0.18326059
#> [151] 0.24762588 0.24931606 0.24458374 0.22259819 0.19244221 0.21540593
#> [157] 0.23215149 0.23553319 0.19846509 0.18695391 0.24799805 0.24626561
#> [163] 0.22976620 0.23641655 0.23560655 0.16635897 0.24594338 0.19614054
#> [169] 0.24786615 0.24989949 0.16008534 0.24996078 0.24751946 0.24873987
#> [175] 0.24421043 0.24917441 0.23970096 0.17441520 0.23105044 0.22312906
#> [181] 0.20709937 0.23278756 0.24985581 0.19363493 0.24640806 0.24765051
#> [187] 0.21017382 0.22987351 0.24170755 0.22275155 0.24736006 0.24831297
#> [193] 0.23128651 0.23833033 0.20837075 0.21679260 0.21142640 0.21840492
#> [199] 0.24818034 0.24303861 0.24821506 0.19762018 0.20346159 0.24428383
#> [205] 0.14382389 0.14849372 0.24313814 0.24601401 0.24935848 0.24661164
#> [211] 0.24756329 0.24806019 0.24253227 0.19928412 0.23263302 0.20717701
#> [217] 0.24479678 0.20588746 0.24376446 0.24908034 0.17761584 0.24947801
#> [223] 0.23520374 0.24404345 0.21114685 0.19090962 0.24427853 0.18538394
#> [229] 0.24992995 0.24964179 0.21988514 0.10960368 0.23482856 0.24541304
#> [235] 0.24892937 0.24849096 0.24561706 0.20954333 0.22813710 0.24942201
#> [241] 0.24206302 0.23077602 0.24984182 0.22494911 0.22520300 0.16119123
#> [247] 0.23818970 0.24593045 0.24753590 0.23133185 0.24999806 0.24824495
#> [253] 0.24882327 0.24608982 0.24238477 0.22595871 0.13680936 0.24536541
#> [259] 0.24865086 0.24184640 0.23405334 0.24282510 0.24533484 0.23742339
#> [265] 0.19583144 0.24959894 0.24193245 0.20939993 0.24428087 0.22210245
#> [271] 0.24586465 0.24841124 0.22616195 0.13939912 0.24999684 0.22667347
#> [277] 0.14359052 0.24110698 0.20236222 0.23521589 0.24836331 0.24233558
#> [283] 0.21747637 0.23250604 0.23697514 0.21913662 0.23378734 0.23427333
#> [289] 0.24787395 0.23194717 0.24280186 0.24863154 0.24974763 0.19603774
#> [295] 0.19237513 0.24999719 0.24037725 0.23389688 0.24902280 0.24486274
#> [301] 0.23615955 0.21682732 0.18077445 0.24297207 0.19962289 0.20759008
#> [307] 0.24990162 0.24548277 0.21797117 0.24918861 0.23038481 0.24994614
#> [313] 0.24102026 0.23425707 0.24916341 0.22931711 0.24012777 0.24819562
#> [319] 0.20360847 0.24689604 0.21829449 0.20841849 0.23884195 0.16922318
#> [325] 0.24710879 0.24753307 0.24373993 0.24998323 0.22561475 0.24282092
#> [331] 0.24437356 0.18757717 0.24338828 0.21543633 0.24420376 0.24134579
#> [337] 0.23422058 0.24949443 0.24376653 0.23868544 0.24114806 0.21600639
#> [343] 0.24973712 0.24522326 0.14504696 0.19705098 0.23535425 0.24712435
#> [349] 0.23920824 0.18585898 0.20967351 0.24287255 0.24858889 0.24709850
#> [355] 0.20394985 0.24484961 0.21023696 0.24769812 0.16619612 0.24705641
#> [361] 0.12584872 0.21799022 0.23501106 0.24171412 0.24748532 0.24448381
#> [367] 0.24284437 0.24546730 0.23725344 0.19968137 0.24283866 0.24954705
#> [373] 0.22743525 0.24883215 0.21055419 0.24939901 0.24999339 0.24784445
#> [379] 0.24757056 0.19325383 0.24046804 0.22750020 0.22288813 0.24999033
#> [385] 0.22557167 0.24896252 0.15326782 0.24897798 0.16705614 0.22795409
#> [391] 0.21748915 0.24440402 0.22835737 0.24933093 0.24833924 0.22283613
#> [397] 0.20310558 0.24974426 0.23195499 0.23715564 0.22909358 0.18536545
#> [403] 0.23618948 0.24278208 0.24209761 0.24317869 0.24340905 0.24508823
#> [409] 0.21139115 0.23265128 0.24285326 0.24350816 0.22691399 0.22962101
#> [415] 0.23937830 0.24635875 0.23144142 0.24602052 0.19581970 0.19644406
#> [421] 0.21174981 0.23214445 0.24887289 0.17600742 0.24542974 0.20200732
#> [427] 0.18346218 0.22930688 0.23480279 0.23329020 0.15915456 0.24677777
#> [433] 0.18980847 0.24926124 0.24492250 0.24700577 0.23621876 0.24925500
#> [439] 0.23456156 0.24366253 0.24160517 0.19271909 0.24931856 0.23757966
#> [445] 0.24990697 0.10344777 0.19951946 0.24187898 0.24978366 0.14038672
#> [451] 0.23364664 0.23613865 0.23042200 0.18838497 0.19442030 0.24979891
#> [457] 0.24662522 0.15519232 0.21897423 0.24612380 0.24105241 0.16886923
#> [463] 0.22133433 0.24486616 0.24949792 0.24953715 0.23232149 0.24411050
#> [469] 0.20933151 0.22871074 0.22786409 0.18457652 0.24834575 0.24776745
#> [475] 0.21612358 0.20222840 0.24876915 0.24979893 0.22377977 0.24923575
#> [481] 0.24364327 0.21166598 0.23908146 0.24352267 0.15746999 0.15160059
#> [487] 0.24752673 0.21701081 0.16286750 0.24652723 0.21461351 0.16151988
#> [493] 0.22463320 0.24997069 0.07964278 0.21366102 0.24932296 0.21376613
#> [499] 0.24894915 0.22519687

if (rlang::is_installed("nnet")) {
  sim$trt <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  multinomial_fit <- nnet::multinom(trt ~ x, data = sim, trace = FALSE)
  head(ps_tilt(multinomial_fit, "att", .focal_level = "b"))
}
#> [1] 0.3497955 0.3525547 0.3490816 0.3573406 0.3530506 0.3491334
```
