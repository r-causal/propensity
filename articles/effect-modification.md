# Effect modification and interaction

Two questions get asked with the same words and answered with different
machinery. This vignette separates them, shows a dataset on which they
give different answers, and explains why that is the correct behavior
rather than a discrepancy to be reconciled.

``` r

library(propensity)

# The weight functions report which exposure they found and how they read it.
# That is useful once and repetitive here, so quiet it for the vignette.
options(propensity.quiet = TRUE)
```

## Two different questions

**Effect modification** asks whether the effect of a treatment `A`
differs across strata of some variable `V`. It is a question about
*conditioning*. You split the population by `V`, estimate the effect of
`A` in each part, and compare. `V` is never intervened on and needs no
causal role at all. It can be a biomarker, a region, a calendar year, or
anything else you can group people by.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
answers this with `.by`.

**Interaction** asks what happens under a joint intervention on two
treatments `A` and `E`. It is a question about *setting both*. Both
treatments have to be identifiable: each needs positivity and a
sufficient adjustment set, and the second needs a model conditioning on
the first. `propensity` answers this with `joint_exposure()` for a
declared crossing, or with
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
and
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
for a pair of treatment models.

The distinction is not a technicality about which function to call. The
two questions have different answers, and a variable can modify an
effect while being useless as a target of intervention.

There is a genuine coincidence result, but it is about one variable
analyzed both ways. Conditioning on `V` and intervening on `V` answer
the same question when `V` is independent of whatever carries the effect
heterogeneity, or when there is no heterogeneity to carry. That says
nothing about how effect modification by one variable relates to
interaction with a different one. Those two are free to take any pair of
values, as the sections below show.

## A dataset where they diverge

The design is a surrogate modifier. A covariate `l` carries genuine
heterogeneity in the effect of `a`. A second variable `v` is strongly
associated with `l` but does nothing at all: it does not appear in the
outcome, and setting it would change nothing. The second treatment `e`
depends on `l` and on `a`.

``` r

set.seed(2024)
n <- 2000

l <- rbinom(n, 1, 0.5)
v <- factor(
  ifelse(rbinom(n, 1, 0.15 + 0.7 * l) == 1, "high", "low"),
  levels = c("low", "high")
)
a <- rbinom(n, 1, plogis(-0.4 + 0.9 * l))
e <- rbinom(n, 1, plogis(-0.5 + 1.0 * l + 0.8 * a))
y <- 1 + 0.4 * a + 0.3 * e + 0.5 * l + 1.2 * a * l + 0.5 * a * e + rnorm(n)

dat <- data.frame(l, v, a, e, y)
```

Two facts about this data generating process matter later. The effect of
`a` is larger when `l` is 1, because of the `1.2 * a * l` term. And `v`
is a marker for `l` and nothing more:

``` r

round(prop.table(table(dat$v, dat$l), 1), 3)
#>       
#>            0     1
#>   low  0.842 0.158
#>   high 0.158 0.842
```

Each stratum of `v` is about 84% one value of `l`. That composition is
the whole of what `v` contributes.

By construction the additive interaction between `a` and `e` is 0.50,
the coefficient on `a * e`, since the outcome is linear and `l` averages
out. The effect modification by `v` works out to about 0.91: the effect
of `a` is 0.75 where `l` is 0 and 2.04 where `l` is 1, and the two
strata of `v` mix those in opposite proportions. Two different numbers,
both true, describing two different things.

## Effect modification with `.by`

Weights handle the confounding of `a` by `l`, and the outcome model is
the marginal structural model in `a` and `v`. `.by` then reports the
effect within each stratum along with the difference between them.

``` r

ps_a <- glm(a ~ l, data = dat, family = binomial())
w_a <- wt_ate(ps_a)

em_mod <- lm(y ~ a * v, data = dat, weights = w_a)
ipw(ps_a, em_mod, .by = v)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = a ~ l, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y ~ a * v, data = dat, weights = w_a) 
#> 
#> Marginal estimates:
#>                                 estimate  std.err      z ci.lower ci.upper conf.level   p.value    
#> mean 0 overall                  1.375819 0.034348 40.055  1.30850   1.4431       0.95 < 2.2e-16 ***
#> mean 1 overall                  2.810577 0.039774 70.663  2.73262   2.8885       0.95 < 2.2e-16 ***
#> diff 1 vs 0 overall             1.434758 0.051168 28.040  1.33447   1.5350       0.95 < 2.2e-16 ***
#> mean 0 v = low                  1.239646 0.041393 29.948  1.15852   1.3208       0.95 < 2.2e-16 ***
#> mean 1 v = low                  2.160110 0.056806 38.026  2.04877   2.2714       0.95 < 2.2e-16 ***
#> mean 0 v = high                 1.513086 0.054746 27.638  1.40579   1.6204       0.95 < 2.2e-16 ***
#> mean 1 v = high                 3.466269 0.049330 70.267  3.36958   3.5630       0.95 < 2.2e-16 ***
#> diff 1 vs 0 v = low             0.920464 0.069938 13.161  0.78339   1.0575       0.95 < 2.2e-16 ***
#> diff 1 vs 0 v = high            1.953183 0.073289 26.650  1.80954   2.0968       0.95 < 2.2e-16 ***
#> diff 1 vs 0 v = high vs v = low 1.032719 0.102924 10.034  0.83099   1.2344       0.95 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The effect of `a` is about 0.92 where `v` is low and about 1.95 where
`v` is high, and the difference between the strata is about 1.03.
Nothing here is wrong. The effect of `a` really is larger among people
with `v = high`, and if you wanted to know where the treatment does the
most good, this is the answer.

What it does not say is that `v` causes anything. Every bit of the
difference comes from the strata having different mixes of `l`. Give
someone `v = high` without changing their `l` and their response to `a`
is unchanged. This is why `.by` requires no assumptions about `V`: it
makes no claims about `V` either.

Stabilizing the weights on the modifier is worth doing here, and
[`?ipw`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
describes the recipe under **Effect modification**.

## Interaction with a joint intervention

Now the other question. Both treatments are intervened on, so both need
treatment models. The second conditions on the first, which is the
factorization a joint weight actually has.

``` r

ps_e <- glm(e ~ a * l, data = dat, family = binomial())

models <- joint_wt_models(a = ps_a, e = ps_e)
w_joint <- wt_joint(
  w_a,
  wt_ate(ps_e),
  exposure_type = c("binary", "binary")
)

joint_mod <- lm(y ~ a * e, data = dat, weights = w_joint)
joint_fit <- ipw(models, joint_mod)
joint_fit
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: <joint_wt_models> 
#> 
#> Outcome Model:
#>   Call: lm(formula = y ~ a * e, data = dat, weights = w_joint) 
#> 
#> Marginal estimates:
#>                               estimate  std.err       z ci.lower ci.upper conf.level   p.value    
#> mean a = 0, e = 0 overall     1.198370 0.052887 22.6591  1.09471  1.30203       0.95 < 2.2e-16 ***
#> mean a = 1, e = 0 overall     2.227287 0.054681 40.7327  2.12011  2.33446       0.95 < 2.2e-16 ***
#> mean a = 0, e = 1 overall     1.524161 0.045885 33.2172  1.43423  1.61409       0.95 < 2.2e-16 ***
#> mean a = 1, e = 1 overall     3.053515 0.044054 69.3132  2.96717  3.13986       0.95 < 2.2e-16 ***
#> diff a: 1 vs 0 e = 0          1.028917 0.075456 13.6360  0.88103  1.17681       0.95 < 2.2e-16 ***
#> diff a: 1 vs 0 e = 1          1.529354 0.062107 24.6244  1.40763  1.65108       0.95 < 2.2e-16 ***
#> diff e: 1 vs 0 a = 0          0.325791 0.069820  4.6661  0.18895  0.46264       0.95 3.069e-06 ***
#> diff e: 1 vs 0 a = 1          0.826228 0.065501 12.6140  0.69785  0.95461       0.95 < 2.2e-16 ***
#> diff a: 1 vs 0 e = 1 vs e = 0 0.500437 0.095685  5.2300  0.31290  0.68798       0.95 1.695e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The surface reports the four counterfactual cell means, each treatment’s
simple effects within the levels of the other, and one interaction row.
The interaction is about 0.50, which is the truth the data was built
with.

A declared crossing gives the same rows through one model over the cells
rather than two models of the treatments:

``` r

joint <- causalgenerics::joint_exposure(a = dat$a, e = dat$e)
ps_cells <- nnet::multinom(joint ~ l, trace = FALSE)
w_cells <- wt_ate(predict(ps_cells, type = "probs"), joint)
ipw(ps_cells, lm(y ~ a * e, data = dat, weights = w_cells))
```

Prefer the two-model route when the treatments call for different
adjustment sets, or when the dependence of the second on the first is
what you want to model. Prefer the declared crossing when one model over
the cells is the natural specification.

## Why the two numbers differ

The effect modification contrast is about 1.03. The interaction is about
0.50. These are not two estimates of one quantity that ought to agree.

The first compares the effect of `a` between two groups of people you
did not create, whose difference is their composition. The second
compares the effect of `a` between two worlds you did create, one in
which everyone gets `e` and one in which nobody does. The first is a
difference of two conditional effects, each computed in its own
covariate distribution. The second is a double difference within a
single population.

Neither number constrains the other, and this data generating process
shows that in both directions if you adjust one piece of it at a time.

Make `v` independent of `l` and change nothing else. The effect
modification contrast goes to zero, because the strata of `v` now have
the same mix of `l` and there is nothing left for stratifying on `v` to
find. The interaction stays at 0.50, because `l` had nothing to do with
it. Someone who read the effect modification result as the interaction
would conclude there was no interaction when there is one of 0.50, which
is why this is the case where conflating the two misleads most, not the
case where they agree. At n = 2000 a single run scatters around that
zero with a standard error of about 0.11, so read it as an expectation
rather than as a number one sample will land on.

Now put `v` back and drop the `a * e` term instead, keeping `a * l`. The
reverse happens: the interaction goes to zero while the effect
modification contrast stays near 0.95. The effect of `a` still differs
across strata of `v`, and setting `e` still does nothing to it.

So the pair can be about (0, 0.50), or about (0.95, 0), or, as in the
data above, (1.03, 0.50). Asking which one is “the” interaction is
asking the wrong question.

The practical consequence is about what each result licenses. A large
effect modification contrast tells you where to look for people who
benefit. Only the interaction tells you what happens if you act on the
second thing.

## Interaction is symmetric

The interaction is reported once, under the first treatment’s framing.
That is not a choice about which treatment matters more; the two
framings are the same number. The change in `a`’s effect when `e` goes
from 0 to 1 equals the change in `e`’s effect when `a` goes from 0 to 1:

``` r

est <- joint_fit$estimates
est[est$effect == "diff", c("contrast", "group", "estimate")]
#>    contrast          group  estimate
#> 5 a: 1 vs 0          e = 0 1.0289168
#> 6 a: 1 vs 0          e = 1 1.5293535
#> 7 e: 1 vs 0          a = 0 0.3257911
#> 8 e: 1 vs 0          a = 1 0.8262278
#> 9 a: 1 vs 0 e = 1 vs e = 0 0.5004367

pick <- function(contrast, group) {
  est$estimate[
    est$effect == "diff" & est$contrast == contrast & est$group == group
  ]
}

c(
  # a's effect at e = 1 minus a's effect at e = 0
  a_framing = pick("a: 1 vs 0", "e = 1") - pick("a: 1 vs 0", "e = 0"),
  # e's effect at a = 1 minus e's effect at a = 0
  e_framing = pick("e: 1 vs 0", "a = 1") - pick("e: 1 vs 0", "a = 0"),
  reported = pick("a: 1 vs 0", "e = 1 vs e = 0")
)
#> a_framing e_framing  reported 
#> 0.5004367 0.5004367 0.5004367
```

The three agree exactly rather than approximately, because all of them
are contrasts over the same block of cell means in one stacked system.
Reporting both framings would be one quantity under two labels, so only
the first treatment’s is shown. The equality is pinned in the package’s
tests.

## Odds ratios and effect modification

With a binary outcome, `.by` reports the risk difference and the log
risk ratio for the strata and their contrast, but the log odds ratio
only for the overall rows:

``` r

dat$yb <- rbinom(
  n,
  1,
  plogis(-1.2 + 0.5 * a + 0.3 * e + 0.6 * l + 0.9 * a * l)
)

emb_mod <- glm(
  yb ~ a * v,
  data = dat,
  family = quasibinomial(),
  weights = w_a
)
ipw(ps_a, emb_mod, .by = v)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = a ~ l, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = yb ~ a * v, family = quasibinomial(), data = dat, 
#>     weights = w_a) 
#> 
#> Marginal estimates:
#>                                    estimate  std.err       z  ci.lower ci.upper conf.level
#> mean 0 overall                     0.338597 0.016109 21.0189  0.307023  0.37017       0.95
#> mean 1 overall                     0.563799 0.015292 36.8694  0.533828  0.59377       0.95
#> rd 1 vs 0 overall                  0.225203 0.021977 10.2473  0.182129  0.26828       0.95
#> log(rr) 1 vs 0 overall             0.509889 0.054266  9.3960  0.403528  0.61625       0.95
#> log(or) 1 vs 0 overall             0.926150 0.094087  9.8436  0.741743  1.11056       0.95
#> mean 0 v = low                     0.287234 0.019854 14.4673  0.248321  0.32615       0.95
#> mean 1 v = low                     0.422081 0.023254 18.1507  0.376503  0.46766       0.95
#> mean 0 v = high                    0.390372 0.025330 15.4113  0.340726  0.44002       0.95
#> mean 1 v = high                    0.706656 0.019007 37.1779  0.669402  0.74391       0.95
#> rd 1 vs 0 v = low                  0.134847 0.030509  4.4199  0.075051  0.19464       0.95
#> log(rr) 1 vs 0 v = low             0.384900 0.088198  4.3640  0.212035  0.55777       0.95
#> rd 1 vs 0 v = high                 0.316284 0.031639  9.9968  0.254273  0.37829       0.95
#> log(rr) 1 vs 0 v = high            0.593443 0.070192  8.4545  0.455869  0.73102       0.95
#> rd 1 vs 0 v = high vs v = low      0.181437 0.044136  4.1108  0.094931  0.26794       0.95
#> log(rr) 1 vs 0 v = high vs v = low 0.208543 0.113040  1.8449 -0.013011  0.43010       0.95
#>                                      p.value    
#> mean 0 overall                     < 2.2e-16 ***
#> mean 1 overall                     < 2.2e-16 ***
#> rd 1 vs 0 overall                  < 2.2e-16 ***
#> log(rr) 1 vs 0 overall             < 2.2e-16 ***
#> log(or) 1 vs 0 overall             < 2.2e-16 ***
#> mean 0 v = low                     < 2.2e-16 ***
#> mean 1 v = low                     < 2.2e-16 ***
#> mean 0 v = high                    < 2.2e-16 ***
#> mean 1 v = high                    < 2.2e-16 ***
#> rd 1 vs 0 v = low                  9.873e-06 ***
#> log(rr) 1 vs 0 v = low             1.277e-05 ***
#> rd 1 vs 0 v = high                 < 2.2e-16 ***
#> log(rr) 1 vs 0 v = high            < 2.2e-16 ***
#> rd 1 vs 0 v = high vs v = low      3.943e-05 ***
#> log(rr) 1 vs 0 v = high vs v = low   0.06506 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

This is deliberate. The odds ratio is noncollapsible: a stratum-specific
odds ratio can differ from the population odds ratio, and from the odds
ratio in another stratum, purely because the outcome is distributed
differently in each stratum. Two strata with identical risk differences
and identical risk ratios can still have different odds ratios. A
difference of two odds ratios therefore moves with the outcome’s
baseline risk whether or not the effect is modified at all, and
reporting one under the heading of effect modification would invite
exactly the reading it does not support. The same reasoning removes the
odds ratio from the joint intervention surface entirely.

If you want an odds ratio, take the overall one. If you want effect
modification, take it on a collapsible scale.

## Guardrails on joint weights

A joint weight is not the product of any two weights. The weight for a
joint intervention on `A` and `E` is `1 / [f(A | L) f(E | A, L)]`, and
the second factor conditions on the first treatment. The product of two
*marginal* weights, `1 / [f(A | L) f(E | L)]`, is a different quantity.
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
refuses it:

``` r

marginal_e <- glm(e ~ l, data = dat, family = binomial())
joint_wt_models(a = ps_a, e = marginal_e)
#> Error in `joint_wt_models()`:
#> ! `joint_wt_models()` requires the second model to condition on the first treatment.
#> ✖ `e` does not read "a" on its right-hand side.
#> ℹ A joint weight factorizes as f(a | L) f(e | a, L). The product of two marginal models, f(a | L)
#>   f(e | L), is a different quantity, and it is not the joint weight wherever "e" depends on "a".
#> ℹ Nothing downstream can tell the two apart: the product is an ordinary vector of positive numbers
#>   either way.
#> ℹ Add "a" to the formula of `e`, and model that dependence flexibly rather than as a single
#>   additive term.
```

The refusal exists because nothing downstream can catch the mistake.
Both products are ordinary vectors of positive numbers, both fit an
outcome model without complaint, and both return estimates with standard
errors that look entirely reasonable.

How much damage it does depends on the data, and in a way that is worth
knowing. When the confounder enters the outcome additively, the marginal
product is very nearly unbiased for the interaction: the errors in the
four cell means cancel in the double difference. The failure appears
once the confounder also modifies the treatment effects, which is the
situation an interaction analysis is usually run in. In simulations for
this guardrail, with a confounder modifying both effects and an odds
ratio of about 4.5 for the first treatment’s effect on the second, a
true interaction of 0.75 came back as 0.39, an attenuation of roughly
half. Under stronger dependence and stronger confounding the estimate
changes sign altogether, returning a negative interaction where the
truth is positive. The magnitude scales with the outcome, so no single
number characterizes it; the direction of the failure is the part to
remember.

Two further requirements follow from the same factorization. The
dependence of the second treatment on the first should be modeled
flexibly rather than as a single additive term. An additive term
satisfies the check, since the check reads the formula, but a badly
modeled dependence biases the joint weight just as invisibly as omitting
the term. Above, `e ~ a * l` lets the dependence differ by `l` rather
than forcing one shift. And a continuous component must be stabilized:
the unstabilized density ratio has a heavy right tail on its own, and a
product inherits it.

The check lives on the container.
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
multiplies whatever it is given, so building weights by hand and fitting
an outcome model directly bypasses the guardrail. Route joint analyses
through
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md).

## Two notes on labels

Grouped rows are labeled with treatment levels, as in `"a: 1 vs 0"` and
`"e = 1"`. Those labels come from the data as it was coded, so they are
not portable across two codings of one intervention: the same analysis
with `a` as a factor of `"no"` and `"yes"` produces different label text
for identical estimates. This matters when pooling, since pooled results
are keyed by their labels. Pooling across multiply imputed datasets is
safe whenever the analyses share one data preparation, which they
normally do.

On the joint dose route, where the second treatment is continuous, the
coefficient surface accepts any coding, and the coding shows up in the
coefficient names the rows carry. Each row reports the fit’s own
coefficient under whatever coding it was given, so this is not the
previous paragraph’s situation of one estimate under two labels: an
ordered factor contributes a polynomial contrast column, and the row
named for it reports that coefficient, which is a different quantity
from the indicator’s. The vocabulary surface, whose rows claim to be the
effect at a dose of zero and the slope at the reference level, refuses a
non-indicator coding outright rather than making those claims of a
column that does not support them.

## Learning more

- [`?ipw`](https://r-causal.github.io/causalgenerics/reference/ipw.html),
  sections **Effect modification** and **Joint exposures**, for the
  reported row sets and their restrictions
- [`?wt_joint`](https://r-causal.github.io/propensity/reference/wt_joint.md)
  and
  [`?joint_wt_models`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
  for the factorization and its checks
- [`?wt_ate`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  section **Stabilization**, for stabilizing on a modifier
