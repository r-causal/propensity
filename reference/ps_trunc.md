# Truncate (Winsorize) Propensity Scores

`ps_trunc()` bounds extreme propensity scores to fixed limits, replacing
out-of-range values with the boundary value (a form of *winsorizing*).
The result is a vector or matrix of the same length and dimensions as
`ps`, with no observations removed. This contrasts with
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
which sets extreme values to `NA` (effectively removing those
observations from analysis).

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

  A numeric vector of propensity scores between 0 and 1 (binary
  exposures), or a matrix/data.frame where each column contains
  propensity scores for one level of a categorical exposure.

- method:

  One of `"ps"`, `"pctl"`, or `"cr"`:

  - `"ps"` (default): Truncate directly on propensity score values.
    Values outside `[lower, upper]` are set to the nearest bound. For
    categorical exposures, applies symmetric truncation using `lower` as
    the threshold (delta) and renormalizes rows to sum to 1.

  - `"pctl"`: Truncate at quantiles of the propensity score
    distribution. The `lower` and `upper` arguments specify quantile
    probabilities. For categorical exposures, quantiles are computed
    across all columns.

  - `"cr"`: Truncate to the common range of propensity scores across
    exposure groups (binary exposures only). Bounds are
    `[min(ps[focal]), max(ps[reference])]`. Requires `.exposure`.

- lower, upper:

  Bounds for truncation. Interpretation depends on `method`:

  - `method = "ps"`: Propensity score values (defaults: 0.1 and 0.9).
    For categorical exposures, `lower` is the truncation threshold delta
    (default: 0.01); `upper` is ignored.

  - `method = "pctl"`: Quantile probabilities (defaults: 0.05 and 0.95;
    categorical defaults: 0.01 and 0.99).

  - `method = "cr"`: Ignored; bounds are determined by the data.

- .exposure:

  An exposure vector. Required for method `"cr"` (binary exposure
  vector) and for categorical exposures (factor or character vector)
  with any method.

- .focal_level:

  The value of `.exposure` representing the focal (treated) group. For
  binary exposures, defaults to the higher value. Required for
  [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  and
  [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  with categorical exposures.

- .reference_level:

  The value of `.exposure` representing the reference (control) group.
  Automatically detected if not supplied.

- ...:

  Additional arguments passed to methods.

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

## Value

A `ps_trunc` object (a numeric vector for binary exposures, or a matrix
for categorical exposures). Use
[`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)
to inspect metadata including `method`, `lower_bound`, `upper_bound`,
and `truncated_idx` (positions of modified values).

## Details

Unlike
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
truncation preserves all observations. No `NA` values are introduced;
out-of-range scores are replaced with the boundary value.

For **binary exposures**, each propensity score \\e_i\\ is bounded:

- If \\e_i \< l\\, set \\e_i = l\\ (the lower bound).

- If \\e_i \> u\\, set \\e_i = u\\ (the upper bound).

For **categorical exposures**, values below the threshold are set to the
threshold and each row is renormalized to sum to 1.

**Arithmetic behavior**: Arithmetic operations on `ps_trunc` objects
return plain numeric vectors. Once propensity scores are transformed
(e.g., into weights), the result is no longer a propensity score.

**Combining behavior**: Combining `ps_trunc` objects with
[`c()`](https://rdrr.io/r/base/c.html) requires matching truncation
parameters. Mismatched parameters produce a warning and return a plain
numeric vector.

## References

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).
Dealing with limited overlap in estimation of average treatment effects.
*Biometrika*, 96(1), 187–199.

Walker, A. M., Patrick, A. R., Lauer, M. S., et al. (2013). A tool for
assessing the feasibility of comparative effectiveness research.
*Comparative Effectiveness Research*, 3, 11–20.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for removing (rather than bounding) extreme values,
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
for refitting the propensity model after trimming,
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md),
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md),
[`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)

## Examples

``` r
set.seed(2)
n <- 200
x <- rnorm(n)
z <- rbinom(n, 1, plogis(1.2 * x))
fit <- glm(z ~ x, family = binomial)
ps <- predict(fit, type = "response")

# Truncate to [0.1, 0.9]
ps_t <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
ps_t
#> <ps_trunc{[0.1,0.9], method=ps}[200]>
#>         1         2         3         4         5         6         7         8 
#> 0.2913095 0.5707511 0.8590561 0.2418820 0.4993060 0.5567591 0.7011122 0.4561604 
#>         9        10        11        12        13        14        15        16 
#> 0.9000000 0.4834312 0.6312438 0.7594616 0.4153613 0.2603881 0.8827176 0.1000000 
#>        17        18        19        20        21        22        23        24 
#> 0.7384237 0.5307543 0.7655683 0.6349278 0.9000000 0.2283122 0.8592915 0.9000000 
#>        25        26        27        28        29        30        31        32 
#> 0.5224033 0.1000000 0.6461643 0.3628349 0.7199125 0.5983593 0.7081103 0.6059826 
#>        33        34        35        36        37        38        39        40 
#> 0.7776782 0.4442186 0.3189618 0.3630602 0.1432262 0.2900408 0.3722937 0.4543264 
#>        41        42        43        44        45        46        47        48 
#> 0.4177638 0.1148898 0.3038315 0.8956760 0.6813225 0.9000000 0.4385124 0.4964323 
#>        49        50        51        52        53        54        55        56 
#> 0.4711469 0.2285337 0.3046166 0.9000000 0.3714863 0.8128663 0.2587396 0.1141443 
#>        57        58        59        60        61        62        63        64 
#> 0.4338453 0.7502468 0.7892860 0.8697086 0.1351329 0.9000000 0.3365403 0.5636424 
#>        65        66        67        68        69        70        71        72 
#> 0.6533256 0.3088378 0.1105762 0.3927369 0.5438044 0.2916295 0.2858818 0.6089556 
#>        73        74        75        76        77        78        79        80 
#> 0.4826524 0.6355772 0.5065031 0.2890305 0.8174117 0.7154239 0.7732114 0.1906359 
#>        81        82        83        84        85        86        87        88 
#> 0.7622717 0.1472971 0.3788316 0.1970403 0.1000000 0.8871259 0.3486998 0.4440783 
#>        89        90        91        92        93        94        95        96 
#> 0.4168763 0.6233901 0.8606965 0.8708768 0.2314475 0.1994226 0.1740398 0.2183041 
#>        97        98        99       100       101       102       103       104 
#> 0.9000000 0.5231365 0.3036226 0.3616811 0.7773582 0.5907629 0.4361655 0.3253709 
#>       105       106       107       108       109       110       111       112 
#> 0.2991480 0.9000000 0.7510710 0.9000000 0.4078242 0.4264336 0.2629645 0.4532487 
#>       113       114       115       116       117       118       119       120 
#> 0.6448288 0.8262184 0.6674247 0.6408723 0.8053646 0.7907095 0.5498332 0.3173982 
#>       121       122       123       124       125       126       127       128 
#> 0.8071017 0.5584826 0.8744315 0.4053976 0.2594361 0.6609892 0.3447196 0.6851535 
#>       129       130       131       132       133       134       135       136 
#> 0.1434914 0.1410494 0.6969682 0.6090883 0.7368408 0.1087328 0.8022199 0.8001309 
#>       137       138       139       140       141       142       143       144 
#> 0.7692948 0.7186431 0.9000000 0.1834141 0.3662171 0.6292392 0.3118604 0.5441734 
#>       145       146       147       148       149       150       151       152 
#> 0.7097460 0.3486308 0.6894218 0.6639810 0.3119192 0.2693231 0.7582975 0.4751339 
#>       153       154       155       156       157       158       159       160 
#> 0.7043397 0.3032090 0.8131267 0.2020948 0.7139969 0.6429234 0.5927018 0.6918371 
#>       161       162       163       164       165       166       167       168 
#> 0.6263848 0.3524856 0.4486289 0.6165339 0.2074392 0.2942181 0.9000000 0.1003085 
#>       169       170       171       172       173       174       175       176 
#> 0.2210197 0.7611782 0.7800143 0.7302196 0.5364399 0.6072561 0.2895753 0.3489979 
#>       177       178       179       180       181       182       183       184 
#> 0.4500408 0.2829251 0.7262051 0.1573127 0.2623292 0.2166743 0.6247878 0.2416706 
#>       185       186       187       188       189       190       191       192 
#> 0.6625838 0.7959533 0.5278943 0.6555095 0.3485232 0.6526880 0.2148033 0.5002503 
#>       193       194       195       196       197       198       199       200 
#> 0.2017086 0.4490035 0.7798096 0.6994296 0.4022315 0.3161761 0.3003833 0.3261363 

# Truncate at the 1st and 99th percentiles
ps_trunc(ps, method = "pctl", lower = 0.01, upper = 0.99)
#> <ps_trunc{[0.0900660672917143,0.912020808232711], method=pctl}[200]>
#>          1          2          3          4          5          6          7 
#> 0.29130945 0.57075108 0.85905610 0.24188201 0.49930600 0.55675915 0.70111217 
#>          8          9         10         11         12         13         14 
#> 0.45616038 0.90360175 0.48343119 0.63124383 0.75946160 0.41536130 0.26038808 
#>         15         16         17         18         19         20         21 
#> 0.88271762 0.09006607 0.73842370 0.53075430 0.76556830 0.63492785 0.91202081 
#>         22         23         24         25         26         27         28 
#> 0.22831221 0.85929151 0.90074563 0.52240328 0.09006607 0.64616425 0.36283492 
#>         29         30         31         32         33         34         35 
#> 0.71991248 0.59835929 0.70811028 0.60598256 0.77767823 0.44421858 0.31896179 
#>         36         37         38         39         40         41         42 
#> 0.36306016 0.14322616 0.29004079 0.37229372 0.45432638 0.41776382 0.11488981 
#>         43         44         45         46         47         48         49 
#> 0.30383148 0.89567605 0.68132246 0.90420941 0.43851237 0.49643227 0.47114693 
#>         50         51         52         53         54         55         56 
#> 0.22853370 0.30461661 0.91106441 0.37148631 0.81286633 0.25873963 0.11414425 
#>         57         58         59         60         61         62         63 
#> 0.43384531 0.75024682 0.78928598 0.86970860 0.13513287 0.90793306 0.33654034 
#>         64         65         66         67         68         69         70 
#> 0.56364237 0.65332558 0.30883775 0.11057616 0.39273691 0.54380444 0.29162947 
#>         71         72         73         74         75         76         77 
#> 0.28588182 0.60895559 0.48265241 0.63557725 0.50650310 0.28903045 0.81741170 
#>         78         79         80         81         82         83         84 
#> 0.71542388 0.77321135 0.19063586 0.76227170 0.14729713 0.37883165 0.19704028 
#>         85         86         87         88         89         90         91 
#> 0.09015377 0.88712592 0.34869977 0.44407832 0.41687627 0.62339015 0.86069650 
#>         92         93         94         95         96         97         98 
#> 0.87087680 0.23144748 0.19942258 0.17403983 0.21830407 0.90120123 0.52313647 
#>         99        100        101        102        103        104        105 
#> 0.30362260 0.36168114 0.77735817 0.59076288 0.43616553 0.32537093 0.29914802 
#>        106        107        108        109        110        111        112 
#> 0.90944557 0.75107099 0.90586644 0.40782418 0.42643359 0.26296452 0.45324867 
#>        113        114        115        116        117        118        119 
#> 0.64482877 0.82621843 0.66742471 0.64087231 0.80536462 0.79070955 0.54983317 
#>        120        121        122        123        124        125        126 
#> 0.31739822 0.80710169 0.55848264 0.87443145 0.40539762 0.25943605 0.66098916 
#>        127        128        129        130        131        132        133 
#> 0.34471960 0.68515350 0.14349137 0.14104940 0.69696824 0.60908833 0.73684079 
#>        134        135        136        137        138        139        140 
#> 0.10873282 0.80221987 0.80013088 0.76929476 0.71864309 0.91202081 0.18341413 
#>        141        142        143        144        145        146        147 
#> 0.36621714 0.62923918 0.31186040 0.54417339 0.70974601 0.34863085 0.68942181 
#>        148        149        150        151        152        153        154 
#> 0.66398099 0.31191915 0.26932313 0.75829753 0.47513393 0.70433975 0.30320896 
#>        155        156        157        158        159        160        161 
#> 0.81312672 0.20209485 0.71399686 0.64292343 0.59270175 0.69183710 0.62638476 
#>        162        163        164        165        166        167        168 
#> 0.35248559 0.44862887 0.61653394 0.20743922 0.29421813 0.91200893 0.10030845 
#>        169        170        171        172        173        174        175 
#> 0.22101971 0.76117824 0.78001428 0.73021961 0.53643993 0.60725607 0.28957526 
#>        176        177        178        179        180        181        182 
#> 0.34899794 0.45004081 0.28292513 0.72620510 0.15731275 0.26232917 0.21667429 
#>        183        184        185        186        187        188        189 
#> 0.62478778 0.24167062 0.66258378 0.79595325 0.52789428 0.65550948 0.34852324 
#>        190        191        192        193        194        195        196 
#> 0.65268801 0.21480329 0.50025030 0.20170859 0.44900347 0.77980959 0.69942959 
#>        197        198        199        200 
#> 0.40223149 0.31617609 0.30038326 0.32613630 

# Use truncated scores to calculate weights
wt_ate(ps_t, .exposure = z)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ate; truncated}[200]>
#>   [1] 1.411053 2.329651 1.164068 1.319056 2.002780 2.256110 1.426305 1.838777
#>   [9] 1.111111 2.068547 2.711819 1.316722 1.710458 1.352060 1.132865 1.111111
#>  [17] 1.354236 1.884111 1.306219 2.739185 1.111111 1.295861 1.163749 1.111111
#>  [25] 2.093817 1.111111 1.547594 1.569452 1.389058 1.671237 1.412209 2.537959
#>  [33] 4.497985 1.799268 1.468346 2.754364 1.167169 3.447791 2.686051 1.832597
#>  [41] 2.393697 1.129803 3.291298 1.116475 1.467734 1.111111 2.280437 2.014373
#>  [49] 1.890884 1.296233 1.438056 1.111111 1.591055 1.230215 3.864889 1.128852
#>  [57] 1.766302 4.003953 4.745769 1.149810 7.400124 1.111111 1.507251 1.774175
#>  [65] 1.530630 3.237946 1.124323 2.546234 1.838896 1.411691 1.400328 1.642156
#>  [73] 1.932936 1.573373 1.974322 3.459843 1.223374 1.397773 1.293307 1.235538
#>  [81] 1.311868 1.172741 2.639695 1.245392 1.111111 1.127236 2.867797 1.798815
#>  [89] 2.398793 2.655268 1.161850 1.148268 1.301147 5.014477 1.210712 1.279270
#>  [97] 1.111111 2.097036 1.436003 1.566615 1.286408 1.692727 1.773570 1.482296
#> [105] 1.426835 1.111111 1.331432 1.111111 2.452037 2.345031 1.356787 1.828985
#> [113] 1.550799 1.210334 1.498296 1.560373 5.137812 1.264687 1.818733 3.150616
#> [121] 5.184079 1.790566 1.143600 2.466714 1.350322 1.512884 1.526064 3.176151
#> [129] 1.167531 7.089715 1.434786 2.558123 3.799981 1.121998 1.246541 1.249796
#> [137] 1.299892 1.391511 1.111111 5.452143 1.577827 2.697157 1.453194 2.193817
#> [145] 1.408955 1.535228 1.450491 2.976022 1.453318 3.713012 4.137318 1.905248
#> [153] 1.419769 1.435150 5.351220 1.253282 1.400566 1.555395 1.687189 3.245037
#> [161] 2.676550 1.544367 2.229014 1.621971 4.820689 1.416868 1.111111 1.111492
#> [169] 4.524483 1.313753 4.545750 1.369451 2.157218 2.546188 1.407609 2.865346
#> [177] 2.222021 1.394555 1.377021 1.186680 1.355618 1.276608 2.665158 1.318688
#> [185] 1.509243 1.256355 1.894319 1.525531 2.869249 2.879256 1.273566 2.001002
#> [193] 1.252675 1.814893 1.282364 3.327007 1.672888 1.462365 1.429354 1.483980

# Inspect which observations were truncated
is_unit_truncated(ps_t)
#>   [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
#>  [13] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE
#>  [25] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [37] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
#>  [49] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [61] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [73] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [85]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [97]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE
#> [109] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [121] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [133] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
#> [145] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [157] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#> [169] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [181] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [193] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
