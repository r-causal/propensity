# Tests for ps_tilt(), the exported estimand tilting function, and for the
# identity harness that pins every consumer of the tilt.
#
# The oracles below inline the weight arithmetic as each consumer wrote it while
# the tilt had no single source, so a shared implementation is measured against
# the formulas it replaced rather than against itself. The ipw() pins are the
# estimates and standard errors those same formulas produced from the fixtures
# in this file, recorded to 17 significant digits and compared at 1e-12.

tilt_estimands <- c("ate", "att", "atu", "atm", "ato", "entropy")

tilt_weight_fn <- function(estimand) {
  switch(
    estimand,
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )
}

# The binary weight formulas as R/weights.R wrote them, one per estimand.
oracle_binary_weights <- function(ps, z, estimand) {
  switch(
    estimand,
    ate = (z / ps) + ((1 - z) / (1 - ps)),
    att = ((ps * z) / ps) + ((ps * (1 - z)) / (1 - ps)),
    atu = (((1 - ps) * z) / ps) + (((1 - ps) * (1 - z)) / (1 - ps)),
    atm = pmin(ps, 1 - ps) / (z * ps + (1 - z) * (1 - ps)),
    ato = (1 - ps) * z + ps * (1 - z),
    entropy = (-ps * log(ps) - (1 - ps) * log(1 - ps)) /
      (z * ps + (1 - z) * (1 - ps))
  )
}

# The categorical tilting switch as calculate_categorical_weights() wrote it,
# divided by the propensity score of the level each unit was assigned.
oracle_categorical_weights <- function(
  ps,
  exposure,
  estimand,
  focal_level = NULL
) {
  levels_exp <- levels(exposure)
  k <- nlevels(exposure)
  n <- length(exposure)
  z <- matrix(0, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    z[exposure == levels_exp[j], j] <- 1
  }
  e_actual <- rowSums(z * ps)

  h_e <- switch(
    estimand,
    ate = rep(1, n),
    att = ps[, which(levels_exp == focal_level)],
    atu = 1 - ps[, which(levels_exp == focal_level)],
    ato = 1 / rowSums(1 / ps),
    atm = apply(ps, 1, min),
    entropy = {
      ps_safe <- ps
      ps_safe[ps == 0] <- .Machine$double.eps
      -rowSums(ps_safe * log(ps_safe))
    }
  )

  as.double(h_e / e_actual)
}

# ---- fixtures ---------------------------------------------------------------

tilt_ps_matrix <- function() {
  matrix(
    c(0.2, 0.5, 0.3, 0.25, 0.5, 0.25),
    nrow = 2,
    dimnames = list(NULL, c("a", "b", "c"))
  )
}

tilt_grid_binary <- function() {
  withr::local_seed(20240817)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  data.frame(x1, x2, z, y)
}

# The outcome model carries the exposure alone so the same fit serves both
# standard error methods; the linearization rejects a covariate-adjusted one.
tilt_grid_binary_models <- function(dat, estimand) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- tilt_weight_fn(estimand)(ps_mod)
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

tilt_grid_categorical <- function() {
  withr::local_seed(20240817)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  lab <- ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c"))
  a <- factor(lab, levels = c("a", "b", "c"))
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  data.frame(x1, x2, a, y)
}

tilt_grid_ps_multinom <- function(dat) {
  nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

tilt_grid_ps_matrix <- function(ps_mod, dat) {
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$a)
  ps
}

tilt_grid_categorical_weights <- function(ps, dat, estimand) {
  args <- list(ps, dat$a, exposure_type = "categorical")
  if (estimand %in% c("att", "atu")) {
    args$.focal_level <- "b"
  }
  do.call(tilt_weight_fn(estimand), args)
}

tilt_grid_categorical_models <- function(dat, estimand) {
  ps_mod <- tilt_grid_ps_multinom(dat)
  ps <- tilt_grid_ps_matrix(ps_mod, dat)
  wts <- tilt_grid_categorical_weights(ps, dat, estimand)
  outcome_mod <- glm(
    y ~ a + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

tilt_grid_continuous <- function() {
  withr::local_seed(20240817)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  a <- 0.5 + 0.8 * x1 - 0.3 * x2 + rnorm(n)
  y <- 1 + 0.6 * a + 0.4 * x1 + rnorm(n)
  data.frame(x1, x2, a, y)
}

tilt_grid_continuous_models <- function(dat) {
  ps_mod <- lm(a ~ x1 + x2, data = dat)
  wts <- wt_ate(
    fitted(ps_mod),
    dat$a,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  outcome_mod <- lm(y ~ a, data = dat, weights = wts)
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

# The pins below are recorded values, so the tolerance is set by what varies
# between the platforms that record them and the platforms that read them: the
# standard errors move in their eighth significant digit with the BLAS the
# sandwich is solved through. A change to the tilt itself moves these numbers by
# orders of magnitude more than that, so 1e-6 still catches what the pins exist
# to catch.
expect_ipw_pin <- function(result, pin) {
  expect_equal(result$estimates$estimate, pin$estimate, tolerance = 1e-6)
  expect_equal(result$estimates$std.err, pin$std.err, tolerance = 1e-6)
}

# ---- recorded ipw() estimates and standard errors ---------------------------

ipw_pins_binary_mestimation <- list(
  ate = list(
    estimate = c(0.29010061662792452, 0.65355441957308824, 1.2038273389659033),
    std.err = c(0.055559761142279257, 0.13544099333129908, 0.24371600012893724)
  ),
  att = list(
    estimate = c(0.28723748102862745, 0.62799353354706233, 1.186409382241653),
    std.err = c(0.056462238438017391, 0.13563518181106302, 0.246937056443264)
  ),
  atu = list(
    estimate = c(0.29251832822841661, 0.6765610946076448, 1.220189837432075),
    std.err = c(0.056379563698840056, 0.1394880601464632, 0.24798846083983966)
  ),
  atm = list(
    estimate = c(0.27934370422833743, 0.61922897371939933, 1.1541992959088954),
    std.err = c(0.056484381389992684, 0.13564959899654924, 0.24601551064917018)
  ),
  ato = list(
    estimate = c(0.28675427499435874, 0.64395454051039847, 1.1886349793765985),
    std.err = c(0.055708828887160958, 0.13539854701581117, 0.24387530467582558)
  ),
  entropy = list(
    estimate = c(0.28764741159939494, 0.64653820066298817, 1.1926924294331949),
    std.err = c(0.055650833369477321, 0.13536608297866395, 0.24375269962329513)
  )
)

ipw_pins_binary_linearization <- list(
  ate = list(
    estimate = c(0.29010061662792452, 0.65355441957308824, 1.2038273389659033),
    std.err = c(0.05565259546448844, 0.13566729677145253, 0.24412324304708569)
  ),
  att = list(
    estimate = c(0.28723748102862778, 0.62799353354706255, 1.1864093822416544),
    std.err = c(0.056556577605207589, 0.13586179535924031, 0.24734966568287917)
  ),
  atm = list(
    estimate = c(0.27934370422833715, 0.61922897371939944, 1.1541992959088943),
    std.err = c(0.056578757722559297, 0.13587623911747534, 0.24642657828266731)
  ),
  ato = list(
    estimate = c(0.28675427499435863, 0.64395454051039813, 1.1886349793765982),
    std.err = c(0.0558019126757743, 0.13562477345093912, 0.2442828079879289)
  )
)

ipw_pins_categorical <- list(
  ate = list(
    estimate = c(
      0.18518609022915994,
      0.38549082685120101,
      0.74998531518016087,
      0.24190058580946253,
      0.47895058532323886,
      0.98810957294555379
    ),
    std.err = c(
      0.07370085660243654,
      0.15827889763094391,
      0.30537142678528306,
      0.063028450560498342,
      0.14085937058807954,
      0.26704183639484924
    )
  ),
  att = list(
    estimate = c(
      0.23503487569004813,
      0.47864549227058389,
      0.95804910021072553,
      0.28797324844140643,
      0.56084730270408889,
      1.1893836755320653
    ),
    std.err = c(
      0.074189387468617224,
      0.17045750525166206,
      0.31407718897353298,
      0.077306882728530216,
      0.1687508796869126,
      0.33781405448735818
    )
  ),
  atu = list(
    estimate = c(
      0.16783773199537358,
      0.35401541305883766,
      0.67901527073108536,
      0.22566208931154194,
      0.451791065001444,
      0.91870965824704498
    ),
    std.err = c(
      0.079166411285920404,
      0.16717109161662033,
      0.32559748860276222,
      0.061377163530248746,
      0.13906871180386729,
      0.25812484291513893
    )
  ),
  atm = list(
    estimate = c(
      0.13419576228219657,
      0.27156935089552076,
      0.54005821707841173,
      0.2313883024517899,
      0.43048717393012403,
      0.9513708311615946
    ),
    std.err = c(
      0.074267820896721209,
      0.15436847983108382,
      0.30255157951499645,
      0.065013903067243278,
      0.13518175816142475,
      0.27426271138517139
    )
  ),
  ato = list(
    estimate = c(
      0.16302174139019077,
      0.33488416395741549,
      0.65816394753796714,
      0.24455601263873136,
      0.46795452639637442,
      1.003043233379884
    ),
    std.err = c(
      0.071839798207293637,
      0.15286603982546662,
      0.2953462742299377,
      0.063073657017100243,
      0.13654576698609627,
      0.26728072736989195
    )
  ),
  entropy = list(
    estimate = c(
      0.1766960656574352,
      0.36487588395278214,
      0.71460527020915499,
      0.24472375523203643,
      0.47614942654947556,
      1.0016528254102726
    ),
    std.err = c(
      0.072382294121891719,
      0.15469084299699334,
      0.29897349936836826,
      0.062576947248800266,
      0.13804654083112936,
      0.26518819630200902
    )
  )
)

ipw_pins_continuous <- list(
  ate = list(
    estimate = 0.37465254468061776,
    std.err = 0.067395457827484198
  )
)

# ---- ps_tilt(): the binary tilt ---------------------------------------------

test_that("ps_tilt() returns the binary tilt of every estimand", {
  ps <- c(0.2, 0.5, 0.8)

  expect_equal(ps_tilt(ps, "ate"), c(1, 1, 1))
  expect_equal(ps_tilt(ps, "att"), c(0.2, 0.5, 0.8))
  expect_equal(ps_tilt(ps, "atu"), c(0.8, 0.5, 0.2))
  expect_equal(ps_tilt(ps, "atm"), c(0.2, 0.5, 0.2))
  expect_equal(ps_tilt(ps, "ato"), c(0.16, 0.25, 0.16))
  expect_equal(
    ps_tilt(ps, "entropy"),
    c(0.50040242353818787, 0.69314718055994529, 0.50040242353818787)
  )
})

test_that("ps_tilt() returns a plain double the length of the propensity score", {
  ps <- c(0.2, 0.5, 0.8)

  for (estimand in tilt_estimands) {
    tilt <- ps_tilt(ps, estimand)
    expect_type(tilt, "double")
    expect_length(tilt, 3)
    expect_null(attributes(tilt))
  }
})

test_that("ps_tilt() passes missing propensity scores through", {
  for (estimand in tilt_estimands) {
    expect_equal(ps_tilt(c(0.2, NA), estimand)[[2]], NA_real_)
  }

  expect_equal(ps_tilt(c(0.2, NA), "att"), c(0.2, NA))
  expect_equal(ps_tilt(c(0.2, NA), "ato"), c(0.16, NA))
  # ate has a constant tilt and would otherwise report 1 for an observation
  # whose propensity score is unknown.
  expect_equal(ps_tilt(c(0.2, NA), "ate"), c(1, NA))
})

test_that("ps_tilt() returns a plain double from named input", {
  named_ps <- c(a = 0.2, b = 0.5, c = 0.8)

  rownamed <- tilt_ps_matrix()
  rownames(rownamed) <- c("r1", "r2")

  one_row <- rownamed[1, , drop = FALSE]

  for (estimand in tilt_estimands) {
    focal <- if (estimand %in% c("att", "atu")) "b" else NULL

    expect_null(names(ps_tilt(named_ps, estimand)))
    expect_null(names(ps_tilt(rownamed, estimand, .focal_level = focal)))
    expect_null(names(ps_tilt(one_row, estimand, .focal_level = focal)))
    expect_null(names(ps_tilt(
      as.data.frame(rownamed),
      estimand,
      .focal_level = focal
    )))
  }
})

# ---- ps_tilt(): the categorical tilt ----------------------------------------

test_that("ps_tilt() returns the categorical tilt of every estimand", {
  ps <- tilt_ps_matrix()

  expect_equal(ps_tilt(ps, "ate"), c(1, 1))
  expect_equal(ps_tilt(ps, "att", .focal_level = "b"), c(0.3, 0.25))
  expect_equal(ps_tilt(ps, "atu", .focal_level = "b"), c(0.7, 0.75))
  expect_equal(ps_tilt(ps, "atm"), c(0.2, 0.25))
  expect_equal(ps_tilt(ps, "ato"), c(0.096774193548387094, 0.1))
  expect_equal(
    ps_tilt(ps, "entropy"),
    c(1.0296530140645737, 1.0397207708399179)
  )
})

test_that("ps_tilt() gives a data frame the same tilt as the matrix", {
  ps <- tilt_ps_matrix()
  ps_df <- as.data.frame(ps)

  for (estimand in tilt_estimands) {
    focal <- if (estimand %in% c("att", "atu")) "b" else NULL
    expect_equal(
      ps_tilt(ps_df, estimand, .focal_level = focal),
      ps_tilt(ps, estimand, .focal_level = focal)
    )
  }
})

test_that("ps_tilt() gives a missing tilt for a row with a missing score", {
  ps <- tilt_ps_matrix()
  ps[1, "c"] <- NA

  for (estimand in tilt_estimands) {
    focal <- if (estimand %in% c("att", "atu")) "b" else NULL
    tilt <- ps_tilt(ps, estimand, .focal_level = focal)
    # The focal column of row one is known, but the row is no longer a
    # probability vector, so no estimand reports a tilt for it.
    expect_equal(tilt[[1]], NA_real_)
    expect_false(is.na(tilt[[2]]))
  }
})

test_that("ps_tilt() resolves .focal_level against the column names", {
  ps <- tilt_ps_matrix()
  expect_equal(ps_tilt(ps, "att", .focal_level = "c"), c(0.5, 0.25))

  # Column order is irrelevant: the focal level names the column to read.
  reordered <- ps[, c("c", "b", "a"), drop = FALSE]
  expect_equal(ps_tilt(reordered, "att", .focal_level = "c"), c(0.5, 0.25))

  # Prediction columns keep the level after their `.pred_` prefix.
  prefixed <- ps
  colnames(prefixed) <- paste0(".pred_", colnames(ps))
  expect_equal(ps_tilt(prefixed, "att", .focal_level = "c"), c(0.5, 0.25))
})

test_that("ps_tilt() prefers an exact column name over a prefix-stripped one", {
  ps <- matrix(
    c(0.5, 0.25, 0.3, 0.25, 0.2, 0.5),
    nrow = 2,
    dimnames = list(NULL, c(".pred_a", "b", "a"))
  )

  # Both `.pred_a` and `a` match the level `a` once the prefix comes off. The
  # exact name wins rather than whichever column sits first.
  expect_equal(ps_tilt(ps, "att", .focal_level = "a"), c(0.2, 0.5))
  expect_equal(ps_tilt(ps, "atu", .focal_level = "a"), c(0.8, 0.5))
})

test_that("ps_tilt() rejects a focal level two columns tie for", {
  duplicated_exact <- matrix(
    c(0.2, 0.5, 0.3, 0.25, 0.5, 0.25),
    nrow = 2,
    dimnames = list(NULL, c("a", "b", "a"))
  )
  expect_error(
    ps_tilt(duplicated_exact, "att", .focal_level = "a"),
    class = "propensity_tilt_focal_error"
  )

  duplicated_prefix <- matrix(
    c(0.2, 0.5, 0.3, 0.25, 0.5, 0.25),
    nrow = 2,
    dimnames = list(NULL, c(".pred_a", "b", ".pred_a"))
  )
  expect_error(
    ps_tilt(duplicated_prefix, "att", .focal_level = "a"),
    class = "propensity_tilt_focal_error"
  )
})

# ---- ps_tilt(): focal level errors ------------------------------------------

test_that("ps_tilt() requires a focal level for the categorical att and atu tilts", {
  ps <- tilt_ps_matrix()

  expect_error(ps_tilt(ps, "att"), class = "propensity_focal_required_error")
  expect_error(ps_tilt(ps, "atu"), class = "propensity_focal_required_error")
})

test_that("ps_tilt() rejects a focal level that names no column", {
  ps <- tilt_ps_matrix()

  expect_error(
    ps_tilt(ps, "att", .focal_level = "d"),
    class = "propensity_tilt_focal_error"
  )
  expect_error(
    ps_tilt(ps, "att", .focal_level = c("a", "b")),
    class = "propensity_tilt_focal_error"
  )
})

test_that("ps_tilt() rejects an unnamed matrix when a focal level is needed", {
  ps <- unname(tilt_ps_matrix())

  expect_error(
    ps_tilt(ps, "att", .focal_level = "b"),
    class = "propensity_matrix_names_error"
  )
})

test_that("ps_tilt() rejects a focal level the tilt does not use", {
  expect_error(
    ps_tilt(c(0.2, 0.5), "att", .focal_level = "a"),
    class = "propensity_tilt_focal_error"
  )
  expect_error(
    ps_tilt(tilt_ps_matrix(), "ato", .focal_level = "b"),
    class = "propensity_tilt_focal_error"
  )
})

# ---- ps_tilt(): argument validation -----------------------------------------

test_that("ps_tilt() rejects an estimand that is not a tilt", {
  expect_error(
    ps_tilt(c(0.2, 0.5), "cens"),
    class = "propensity_unknown_estimand_error"
  )
  expect_error(
    ps_tilt(tilt_ps_matrix(), "overlap"),
    class = "propensity_unknown_estimand_error"
  )
  expect_error(
    ps_tilt(c(0.2, 0.5), c("ate", "att")),
    class = "propensity_unknown_estimand_error"
  )
})

test_that("ps_tilt() rejects extra arguments", {
  expect_error(
    ps_tilt(c(0.2, 0.5), "ate", TRUE),
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    ps_tilt(tilt_ps_matrix(), "att", focal_level = "b"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ps_tilt() requires an estimand", {
  expect_error(
    ps_tilt(c(0.2, 0.5)),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    ps_tilt(tilt_ps_matrix()),
    class = "propensity_missing_arg_error"
  )
})

test_that("ps_tilt() has no method for a propensity score it cannot read", {
  expect_error(ps_tilt("0.5", "ate"), class = "propensity_method_error")
})

test_that("ps_tilt() sends a modified propensity score to its plain scores", {
  withr::local_options(propensity.quiet = TRUE)
  trimmed <- ps_trim(c(0.05, 0.4, 0.6, 0.95), method = "ps", lower = 0.1)

  expect_error(ps_tilt(trimmed, "att"), class = "propensity_method_error")
  expect_snapshot(error = TRUE, ps_tilt(trimmed, "att"))

  # The route the message names works, and the units trimming set to NA keep an
  # NA tilt through it.
  expect_equal(
    ps_tilt(as.numeric(trimmed), "att"),
    c(NA, 0.4, 0.6, NA)
  )
})

# ---- ps_tilt(): the range contract ------------------------------------------

test_that("ps_tilt() rejects a binary propensity score at the bounds", {
  # The same range check the binary weight functions apply before they tilt.
  expect_error(ps_tilt(c(0, 0.5), "att"), class = "propensity_range_error")
  expect_error(ps_tilt(c(0.5, 1), "att"), class = "propensity_range_error")
})

test_that("ps_tilt() rejects a categorical propensity score at the bounds", {
  # The same open interval the categorical weight path requires: a score of
  # exactly zero at a unit's own level divides its weight by zero, so the
  # matrix is rejected before any tilt is evaluated.
  ps <- matrix(
    c(0, 0.2, 1, 0.3, 0, 0.5),
    nrow = 2,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  expect_error(ps_tilt(ps, "atm"), class = "propensity_range_error")
  expect_error(ps_tilt(ps, "ato"), class = "propensity_range_error")
  expect_error(
    ps_tilt(ps, "att", .focal_level = "a"),
    class = "propensity_range_error"
  )
  expect_error(ps_tilt(ps, "entropy"), class = "propensity_range_error")
})

test_that("ps_tilt() rejects a matrix whose rows are not probability vectors", {
  expect_error(
    ps_tilt(cbind(a = c(0.2, 0.5), b = c(0.2, 0.5)), "ate"),
    class = "propensity_matrix_sum_error"
  )
  expect_error(
    ps_tilt(cbind(a = c(-0.2, 0.5), b = c(1.2, 0.5)), "ate"),
    class = "propensity_range_error"
  )
  expect_error(
    ps_tilt(matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "a")), "ate"),
    class = "propensity_matrix_dims_error"
  )
  expect_error(
    ps_tilt(matrix("a", nrow = 2, ncol = 2), "ate"),
    class = "propensity_matrix_type_error"
  )
})

# ---- ps_tilt(): the weight identity -----------------------------------------

test_that("the documented weight identity holds as the example writes it", {
  # The identity the `ps_tilt()` example prints. It compares a plain numeric
  # against a quotient, so it holds only when neither side carries names.
  withr::local_options(propensity.quiet = TRUE)
  set.seed(1)
  n <- 500
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.6 * x))
  sim <- data.frame(x = x, z = z)
  ps <- unname(
    predict(glm(z ~ x, data = sim, family = binomial()), type = "response")
  )

  received <- z * ps + (1 - z) * (1 - ps)
  expect_true(isTRUE(all.equal(
    as.numeric(wt_ato(ps, z, exposure_type = "binary")),
    ps_tilt(ps, "ato") / received
  )))
})

test_that("ps_tilt() divided by the received propensity score gives the binary weights", {
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_binary()
  ps <- as.double(predict(
    glm(z ~ x1 + x2, data = dat, family = binomial()),
    type = "response"
  ))
  received <- dat$z * ps + (1 - dat$z) * (1 - ps)

  for (estimand in tilt_estimands) {
    weights <- tilt_weight_fn(estimand)(ps, dat$z, exposure_type = "binary")
    expect_equal(
      as.double(weights),
      ps_tilt(ps, estimand) / received,
      tolerance = 1e-12
    )
  }
})

# ---- the identity grid: weights ---------------------------------------------

test_that("binary weight functions reproduce the pre-refactor formulas", {
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_binary()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  ps <- as.double(predict(ps_mod, type = "response"))

  for (estimand in tilt_estimands) {
    expected <- oracle_binary_weights(ps, dat$z, estimand)
    wt_fun <- tilt_weight_fn(estimand)

    expect_equal(
      as.double(wt_fun(ps, dat$z, exposure_type = "binary")),
      expected,
      tolerance = 1e-12
    )
    expect_equal(as.double(wt_fun(ps_mod)), expected, tolerance = 1e-12)
  }
})

test_that("categorical weight functions reproduce the pre-refactor formulas", {
  skip_if_not_installed("nnet")
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_categorical()
  ps <- tilt_grid_ps_matrix(tilt_grid_ps_multinom(dat), dat)

  for (estimand in tilt_estimands) {
    expect_equal(
      as.double(tilt_grid_categorical_weights(ps, dat, estimand)),
      oracle_categorical_weights(ps, dat$a, estimand, focal_level = "b"),
      tolerance = 1e-12
    )
  }
})

# ---- the identity grid: ipw() -----------------------------------------------

test_that("ipw() binary M-estimation results are unchanged for every estimand", {
  skip_if_not_installed("deli")
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_binary()

  for (estimand in tilt_estimands) {
    mods <- tilt_grid_binary_models(dat, estimand)
    expect_ipw_pin(
      ipw(mods$ps_mod, mods$outcome_mod),
      ipw_pins_binary_mestimation[[estimand]]
    )
  }
})

test_that("ipw() binary linearization results are unchanged for every estimand it supports", {
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_binary()

  for (estimand in names(ipw_pins_binary_linearization)) {
    mods <- tilt_grid_binary_models(dat, estimand)
    expect_ipw_pin(
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
      ipw_pins_binary_linearization[[estimand]]
    )
  }
})

test_that("ipw() categorical results are unchanged for every estimand", {
  skip_if_not_installed("deli")
  skip_if_not_installed("nnet")
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_categorical()

  for (estimand in tilt_estimands) {
    mods <- tilt_grid_categorical_models(dat, estimand)
    args <- list(mods$ps_mod, mods$outcome_mod)
    if (estimand %in% c("att", "atu")) {
      args$.focal_level <- "b"
    }
    expect_ipw_pin(do.call(ipw, args), ipw_pins_categorical[[estimand]])
  }
})

test_that("ipw() continuous results are unchanged", {
  skip_if_not_installed("deli")
  withr::local_options(propensity.quiet = TRUE)
  dat <- tilt_grid_continuous()
  mods <- tilt_grid_continuous_models(dat)

  expect_ipw_pin(ipw(mods$ps_mod, mods$outcome_mod), ipw_pins_continuous$ate)
})
