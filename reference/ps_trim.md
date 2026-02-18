# Trim Propensity Scores

`ps_trim()` applies trimming methods to a propensity-score vector or
matrix, returning a new vector/matrix of the *same length/dimensions*,
with trimmed entries replaced by `NA.` You can inspect further metadata
in `ps_trim_meta(x)`. After running `ps_trim()`, you should refit the
model with
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md).

## Usage

``` r
ps_trim(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
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

  One of `c("ps", "adaptive", "pctl", "pref", "cr", "optimal")`. For
  categorical exposures, only `"ps"` and `"optimal"` are supported.

- lower, upper:

  Numeric cutoffs or quantiles. If `NULL`, defaults vary by method. For
  categorical exposures with method `"ps"`, `lower` represents the
  symmetric trimming threshold (delta).

- .exposure:

  For methods like `"pref"` or `"cr"`, a vector for a binary exposure.
  For categorical exposures with method `"optimal"`, must be a factor or
  character vector.

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

A `ps_trim` object (numeric vector or matrix). The attribute
`ps_trim_meta` stores metadata.

## Details

The returned object is a **`ps_trim`** vector/matrix of the same
length/dimensions as `ps`, but with trimmed entries replaced by `NA`. An
attribute `ps_trim_meta` contains:

- `method`: Which trimming method was used

- `keep_idx`: Indices retained

- `trimmed_idx`: Indices replaced by `NA`

- Possibly other fields such as final cutoffs, etc.

For categorical exposures:

- **Symmetric trimming** (`method = "ps"`): Removes observations where
  any propensity score falls below the threshold delta (specified via
  `lower`).

- **Optimal trimming** (`method = "optimal"`): Uses the Yang et
  al. (2016) approach for multi-category treatments.

**Arithmetic behavior**: Arithmetic operations on `ps_trim` objects
return numeric vectors, not `ps_trim` objects. This is intentional -
once you transform propensity scores (e.g., `1/ps` for weights), the
result is no longer a propensity score.

**NA handling**: Trimmed values are set to `NA`. Operations that don't
handle `NA` values will propagate them (e.g.,
[`sum()`](https://rdrr.io/r/base/sum.html) returns `NA` unless
`na.rm = TRUE`).

**Metadata tracking**: The `trimmed_idx` and `keep_idx` are updated when
subsetting or reordering:

- Subsetting with `[` updates indices to new positions

- [`sort()`](https://rdrr.io/r/base/sort.html) reorders data and updates
  indices accordingly

- [`unique()`](https://rdrr.io/r/base/unique.html) may change lengths
  but preserves the class

- [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) removes trimmed
  values and updates indices

**Combining behavior**: When combining `ps_trim` objects with
[`c()`](https://rdrr.io/r/base/c.html), metadata must match (same
trimming parameters). Mismatched metadata triggers a warning and returns
a numeric vector.

## See also

[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
for bounding/winsorizing instead of discarding,
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md),
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)

## Examples

``` r
set.seed(2)
n <- 300
x <- rnorm(n)
z <- rbinom(n, 1, plogis(1.3 * x))
fit <- glm(z ~ x, family = binomial)
ps <- predict(fit, type = "response")

ps_trim(ps, method = "adaptive")
#> <ps_trim; trimmed 44 of [300]>
#>         1         2         3         4         5         6         7         8 
#> 0.1780112 0.4934483 0.8725819 0.1353577 0.4025855 0.4752489 0.6683913 0.3506150 
#>         9        10        11        12        13        14        15        16 
#>        NA 0.3831809 0.5738044 0.7467782 0.3038553 0.1508037 0.8997256        NA 
#>        17        18        19        20        21        22        23        24 
#> 0.7187207 0.4419184 0.7548592 0.5787647        NA 0.1244366 0.8728587        NA 
#>        25        26        27        28        29        30        31        32 
#> 0.4313639        NA 0.5939254 0.2474276 0.6938170 0.5298266 0.6778670 0.5399668 
#>        33        34        35        36        37        38        39        40 
#> 0.7707828 0.3366773 0.2037945 0.2476600        NA 0.1768609 0.2572600 0.3484614 
#>        41        42        43        44        45        46        47        48 
#> 0.3065403        NA 0.1895191        NA 0.6415565        NA 0.3300895 0.3990495 
#>        49        50        51        52        53        54        55        56 
#> 0.3683877 0.1246121 0.1902500        NA 0.2564150 0.8160957 0.1494022        NA 
#>        57        58        59        60        61        62        63        64 
#> 0.3247367 0.7345271 0.7859021 0.8849770        NA        NA 0.2208817 0.4841803 
#>        65        66        67        68        69        70        71        72 
#> 0.6036086 0.1941979        NA 0.2790100 0.4585602 0.1783019 0.1731103 0.5439312 
#>        73        74        75        76        77        78        79        80 
#> 0.3822372 0.5796397 0.4114856 0.1759469 0.8218240 0.6877562 0.7649259        NA 
#>        81        82        83        84        85        86        87        88 
#> 0.7505008        NA 0.2641422 0.1005949        NA        NA 0.2330120 0.3365148 
#>        89        90        91        92        93        94        95        96 
#> 0.3055473 0.5632496 0.8745082 0.8863194 0.1269292 0.1023453        NA 0.1166039 
#>        97        98        99       100       101       102       103       104 
#>        NA 0.4322875 0.1893249 0.2462384 0.7703638 0.5197606 0.3273939 0.2099624 
#>       105       106       107       108       109       110       111       112 
#> 0.1851823        NA 0.7356255        NA 0.2954896 0.3163021 0.1530042 0.3471981 
#>       113       114       115       116       117       118       119       120 
#> 0.5921212 0.8328274 0.6227067 0.5867797 0.8065734 0.7877456 0.4663064 0.2023006 
#>       121       122       123       124       125       126       127       128 
#> 0.8087856 0.4774812 0.8903829 0.2928150 0.1499937 0.6139848 0.2290136 0.6467536 
#>       129       130       131       132       133       134       135       136 
#>        NA        NA 0.6627758 0.5441083 0.7165979        NA 0.8025574 0.7998822 
#>       137       138       139       140       141       142       143       144 
#> 0.7597742 0.6921037        NA        NA 0.2509264 0.5711077 0.1970442 0.4590332 
#>       145       146       147       148       149       150       151       152 
#> 0.6800801 0.2329425 0.6525433 0.6180388 0.1970997 0.1584870 0.7452343 0.3731672 
#>       153       154       155       156       157       158       159       160 
#> 0.6727629 0.1889404 0.8164247 0.1043218 0.6858279 0.5895482 0.5223260 0.6558189 
#>       161       162       163       164       165       166       167       168 
#> 0.5672709 0.2368400 0.3418011 0.5540597 0.1083159 0.1806594        NA        NA 
#>       169       170       171       172       173       174       175       176 
#> 0.1187105 0.7490531 0.7738375 0.7077039 0.4491471 0.5416643 0.1764396 0.2333126 
#>       177       178       179       180       181       182       183       184 
#> 0.3434474 0.1704628 0.7023006        NA 0.1524604 0.1153463 0.5651259 0.1351849 
#>       185       186       187       188       189       190       191       192 
#> 0.6161453 0.7945146 0.4382952 0.6065642 0.2328341 0.6027459 0.1139088 0.4037497 
#>       193       194       195       196       197       198       199       200 
#> 0.1040353 0.3422376 0.7735701 0.6661116 0.2893391 0.2011360 0.1863223 0.2107038 
#>       201       202       203       204       205       206       207       208 
#> 0.5327159 0.1544197        NA 0.5052145 0.1642856 0.5622725 0.3887452 0.3163883 
#>       209       210       211       212       213       214       215       216 
#> 0.6341619 0.5095796 0.7582940 0.2665601        NA        NA 0.4774214 0.5852567 
#>       217       218       219       220       221       222       223       224 
#> 0.8030404 0.1067663 0.1338595 0.8855459 0.5671613 0.6707177 0.2000938        NA 
#>       225       226       227       228       229       230       231       232 
#> 0.4538222 0.5912550 0.3993802 0.7230688 0.3094663        NA 0.8702859        NA 
#>       233       234       235       236       237       238       239       240 
#> 0.2415827 0.4195113        NA 0.1483657 0.2541910 0.5957602 0.4626159 0.2496911 
#>       241       242       243       244       245       246       247       248 
#> 0.4020747        NA 0.7936546 0.6179515 0.6899355 0.2556513 0.4650462 0.5270157 
#>       249       250       251       252       253       254       255       256 
#> 0.2803714 0.2850014 0.6020472 0.3002118 0.3705819 0.3257820 0.7087374 0.5960756 
#>       257       258       259       260       261       262       263       264 
#> 0.3306472 0.3363861 0.7464321 0.3725508 0.8297224 0.3965473 0.3250342        NA 
#>       265       266       267       268       269       270       271       272 
#> 0.2345937        NA        NA 0.5886632 0.6106387 0.4172303        NA 0.4137080 
#>       273       274       275       276       277       278       279       280 
#> 0.4312953 0.1206438 0.5164038 0.8752400 0.4176263 0.7591012 0.3016993 0.4527210 
#>       281       282       283       284       285       286       287       288 
#>        NA 0.6501483 0.8168531 0.2393292 0.8311549 0.8851143 0.7936328 0.4317872 
#>       289       290       291       292       293       294       295       296 
#> 0.8222307 0.3983355 0.1356224 0.6302455 0.4602807 0.3527491 0.8562989 0.3153233 
#>       297       298       299       300 
#> 0.5741435        NA 0.1008201 0.2253829 

# Coercion behavior with ps_trim objects
ps_trim1 <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
ps_trim2 <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)

# Compatible objects combine silently
c(ps_trim1[1:50], ps_trim2[51:100])  # Returns ps_trim object
#> <ps_trim; trimmed 20 of [100]>
#>         1         2         3         4         5         6         7         8 
#> 0.1780112 0.4934483 0.8725819 0.1353577 0.4025855 0.4752489 0.6683913 0.3506150 
#>         9        10        11        12        13        14        15        16 
#>        NA 0.3831809 0.5738044 0.7467782 0.3038553 0.1508037 0.8997256        NA 
#>        17        18        19        20        21        22        23        24 
#> 0.7187207 0.4419184 0.7548592 0.5787647        NA 0.1244366 0.8728587        NA 
#>        25        26        27        28        29        30        31        32 
#> 0.4313639        NA 0.5939254 0.2474276 0.6938170 0.5298266 0.6778670 0.5399668 
#>        33        34        35        36        37        38        39        40 
#> 0.7707828 0.3366773 0.2037945 0.2476600        NA 0.1768609 0.2572600 0.3484614 
#>        41        42        43        44        45        46        47        48 
#> 0.3065403        NA 0.1895191        NA 0.6415565        NA 0.3300895 0.3990495 
#>        49        50        51        52        53        54        55        56 
#> 0.3683877 0.1246121 0.1902500        NA 0.2564150 0.8160957 0.1494022        NA 
#>        57        58        59        60        61        62        63        64 
#> 0.3247367 0.7345271 0.7859021 0.8849770        NA        NA 0.2208817 0.4841803 
#>        65        66        67        68        69        70        71        72 
#> 0.6036086 0.1941979        NA 0.2790100 0.4585602 0.1783019 0.1731103 0.5439312 
#>        73        74        75        76        77        78        79        80 
#> 0.3822372 0.5796397 0.4114856 0.1759469 0.8218240 0.6877562 0.7649259        NA 
#>        81        82        83        84        85        86        87        88 
#> 0.7505008        NA 0.2641422 0.1005949        NA        NA 0.2330120 0.3365148 
#>        89        90        91        92        93        94        95        96 
#> 0.3055473 0.5632496 0.8745082 0.8863194 0.1269292 0.1023453        NA 0.1166039 
#>        97        98        99       100 
#>        NA 0.4322875 0.1893249 0.2462384 

# Different trim parameters trigger warning
ps_trim3 <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
c(ps_trim1[1:50], ps_trim3[51:100])  # Warning: returns numeric
#> Warning: Converting ps_trim to numeric: different trimming parameters
#> ℹ Metadata cannot be preserved when combining incompatible objects
#> ℹ Use identical objects or explicitly cast to numeric to avoid this warning
#>         1         2         3         4         5         6         7         8 
#> 0.1780112 0.4934483 0.8725819 0.1353577 0.4025855 0.4752489 0.6683913 0.3506150 
#>         9        10        11        12        13        14        15        16 
#>        NA 0.3831809 0.5738044 0.7467782 0.3038553 0.1508037 0.8997256        NA 
#>        17        18        19        20        21        22        23        24 
#> 0.7187207 0.4419184 0.7548592 0.5787647        NA 0.1244366 0.8728587        NA 
#>        25        26        27        28        29        30        31        32 
#> 0.4313639        NA 0.5939254 0.2474276 0.6938170 0.5298266 0.6778670 0.5399668 
#>        33        34        35        36        37        38        39        40 
#> 0.7707828 0.3366773 0.2037945 0.2476600        NA 0.1768609 0.2572600 0.3484614 
#>        41        42        43        44        45        46        47        48 
#> 0.3065403        NA 0.1895191        NA 0.6415565        NA 0.3300895 0.3990495 
#>        49        50        51        52        53        54        55        56 
#> 0.3683877 0.1246121        NA        NA 0.2564150        NA        NA        NA 
#>        57        58        59        60        61        62        63        64 
#> 0.3247367 0.7345271 0.7859021        NA        NA        NA 0.2208817 0.4841803 
#>        65        66        67        68        69        70        71        72 
#> 0.6036086        NA        NA 0.2790100 0.4585602        NA        NA 0.5439312 
#>        73        74        75        76        77        78        79        80 
#> 0.3822372 0.5796397 0.4114856        NA        NA 0.6877562 0.7649259        NA 
#>        81        82        83        84        85        86        87        88 
#> 0.7505008        NA 0.2641422        NA        NA        NA 0.2330120 0.3365148 
#>        89        90        91        92        93        94        95        96 
#> 0.3055473 0.5632496        NA        NA        NA        NA        NA        NA 
#>        97        98        99       100 
#>        NA 0.4322875        NA 0.2462384 

# Cross-class combinations warn and return numeric
psw_obj <- psw(ps[1:50], estimand = "ate")
c(ps_trim1[1:50], psw_obj)  # Warning: returns numeric
#> Warning: Converting ps_trim and psw to numeric
#> ℹ Class-specific attributes and metadata have been dropped
#> ℹ Use explicit casting to numeric to avoid this warning
#>          1          2          3          4          5          6          7 
#> 0.17801124 0.49344831 0.87258189 0.13535765 0.40258554 0.47524889 0.66839132 
#>          8          9         10         11         12         13         14 
#> 0.35061504         NA 0.38318089 0.57380442 0.74677823 0.30385529 0.15080370 
#>         15         16         17         18         19         20         21 
#> 0.89972561         NA 0.71872070 0.44191837 0.75485924 0.57876473         NA 
#>         22         23         24         25         26         27         28 
#> 0.12443658 0.87285871         NA 0.43136394         NA 0.59392536 0.24742762 
#>         29         30         31         32         33         34         35 
#> 0.69381700 0.52982663 0.67786695 0.53996676 0.77078278 0.33667732 0.20379449 
#>         36         37         38         39         40         41         42 
#> 0.24766004         NA 0.17686094 0.25726004 0.34846140 0.30654026         NA 
#>         43         44         45         46         47         48         49 
#> 0.18951912         NA 0.64155648         NA 0.33008953 0.39904946 0.36838768 
#>         50                                                                   
#> 0.12461207 0.17801124 0.49344831 0.87258189 0.13535765 0.40258554 0.47524889 
#>                                                                              
#> 0.66839132 0.35061504 0.92239228 0.38318089 0.57380442 0.74677823 0.30385529 
#>                                                                              
#> 0.15080370 0.89972561 0.02943822 0.71872070 0.44191837 0.75485924 0.57876473 
#>                                                                              
#> 0.93233516 0.12443658 0.87285871 0.91937240 0.43136394 0.02433818 0.59392536 
#>                                                                              
#> 0.24742762 0.69381700 0.52982663 0.67786695 0.53996676 0.77078278 0.33667732 
#>                                                                              
#> 0.20379449 0.24766004 0.06402613 0.17686094 0.25726004 0.34846140 0.30654026 
#>                                                                              
#> 0.04714018 0.18951912 0.91394741 0.64155648 0.92303132 0.33008953 0.39904946 
#>                       
#> 0.36838768 0.12461207 
```
