# Design mismatches between `.data` and the fitted models
#
# `ipw()` estimates the marginal means by setting the exposure column to one
# counterfactual value at a time and rebuilding the outcome design, and it
# rebuilds the propensity design from `.data` when one is supplied. Both
# rebuilds multiply the fitted coefficients positionally. So an outcome formula
# that reads the exposure through a call that cannot survive a constant column,
# and a `.data` column whose type differs from the one the models were fit on,
# both produce either a raw base R error or a design that is silently not the
# fitted one. These tests pin the classed errors that stand in for both, and the
# legitimate rebuilds the guards must leave alone.
#
# The companion `.data` guards, for a frame that cannot be reconstructed, a row
# count that disagrees, and a design width that disagrees, are pinned in
# test-ipw-mestimation.R.

design_data <- function(seed = 2024, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  cov <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  num <- sample(1:3, n, replace = TRUE)
  flag <- sample(c(TRUE, FALSE), n, replace = TRUE)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * (cov == "b")))
  y <- rbinom(n, 1, plogis(-0.4 + 1.0 * z + 0.5 * x1 + 0.3 * num))
  data.frame(x1, cov, num, flag, z, y)
}

design_weights <- function(ps_mod) {
  withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
}

design_outcome <- function(formula, data, wts) {
  glm(
    formula,
    data = data,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
}

design_msg <- function(err) {
  gsub("[[:space:]]+", " ", conditionMessage(err))
}

# ---- an outcome model that reads the exposure through a call ----------------
#
# `y ~ factor(z)` records the term under the label "factor(z)", so the fitted
# levels and the contrast coding belong to that label rather than to the column.
# The rebuild sets the column to one value and re-evaluates the label, which
# leaves `factor(z)` with a single level and no applicable contrasts. On the
# `.data` route that died inside `model.matrix()` as "contrasts can be applied
# only to factors with 2 or more levels"; without `.data` it died as a report
# that the exposure column is missing from the model frame, whose remedy is to
# supply `.data`, which walks into the first error. The guard reads the model
# alone, so both routes reach it and the diagnosis is the same on each.

test_that("ipw() rejects an outcome model that factors the exposure in the formula", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, .data = dat),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "factor(z)", fixed = TRUE)
  expect_false(grepl("contrasts can be applied", msg, fixed = TRUE))
})

test_that("the factored-exposure rejection is the same without .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  expect_match(design_msg(err), "factor(z)", fixed = TRUE)
})

test_that("the factored-exposure error names the term and the rebuild it blocks", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat))
})

test_that("ipw() rejects a factored exposure in a categorical outcome model", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$a <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_mod <- nnet::multinom(a ~ x1, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$a, exposure_type = "categorical")
  )
  out <- design_outcome(y ~ factor(a) + x1, dat, wts)

  for (supplied in list(NULL, dat)) {
    err <- expect_error(
      ipw(ps_mod, out, .data = supplied),
      class = "propensity_ipw_exposure_error"
    )
    expect_match(design_msg(err), "factor(a)", fixed = TRUE)
  }
})

# ---- exposure calls that do rebuild -----------------------------------------
#
# The rejection is not of transformed exposure terms in general. A numeric
# transformation is recomputed from the counterfactual column and rebuilds the
# fitted design, and so does a factor whose levels come from the call rather
# than from the data. Both are honest g-computation on the model as specified.
# A guard keyed on the term's spelling instead of on what it evaluates to would
# reject these, so they carry the discrimination the guard has to make.

test_that("ipw() accepts a numeric transformation of the exposure with .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  # z is 0/1, so I(z^2) is z and the two models are the same model
  out_sq <- design_outcome(y ~ I(z^2) + x1, dat, wts)
  out_plain <- design_outcome(y ~ z + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_sq, .data = dat))
  ref <- ipw(ps_mod, out_plain, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("a numeric transformation of the exposure still asks for .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_sq <- design_outcome(y ~ I(z^2) + x1, dat, wts)

  # The model frame holds `I(z^2)` and no `z` column to set, so the rebuild has
  # nowhere to write. That is the missing-column report, not the guard above.
  expect_error(
    ipw(ps_mod, out_sq),
    class = "propensity_columns_exist_error"
  )
})

test_that("ipw() accepts an exposure factored on levels the call fixes", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  # cut() takes its levels from the breaks, so a constant column still carries
  # both of them and the rebuilt design is the fitted one
  out_cut <- design_outcome(y ~ cut(z, c(-1, 0.5, 2)) + x1, dat, wts)
  out_plain <- design_outcome(y ~ z + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_cut, .data = dat))
  ref <- ipw(ps_mod, out_plain, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("ipw() accepts a factored covariate that is not the exposure", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  # `num` is untouched by the counterfactual rebuild, so factor(num) evaluates
  # to the fitted three-level factor every time
  out <- design_outcome(y ~ z + factor(num), dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("ipw() accepts a transformed covariate rebuilt from .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + log(x1 + 5), dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

# ---- a `.data` column of the wrong type -------------------------------------
#
# The fitted models record each variable's class, and every design rebuilt from
# `.data` is rebuilt under the coding that class implies. A numeric column takes
# one design column; a factor, character, or logical column expands into one per
# non-reference level. Supplying the other kind therefore either fails to
# multiply the coefficients at all or, when the two happen to be the same width,
# multiplies them by different numbers with nothing signaled.

test_that("ipw() rejects a factor .data column where the outcome model fit a numeric", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + num, dat, wts)

  dat_fac <- dat
  dat_fac$num <- factor(dat$num)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bnum\\b")
  expect_match(msg, "factor", fixed = TRUE)
  expect_match(msg, "numeric", fixed = TRUE)
  expect_false(grepl("non-conformable", msg, fixed = TRUE))
})

test_that("ipw() rejects a logical .data column where the outcome model fit a numeric", {
  skip_if_not_installed("deli")
  # The same-width case, and the reason a design-width check cannot stand in for
  # a type check: a logical column takes one design column, exactly as the
  # fitted numeric did, so the rebuild conforms and the marginal means are
  # computed from a design of ones and zeros the model was never fit on.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + num, dat, wts)

  dat_lgl <- dat
  dat_lgl$num <- dat$num > 1

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bnum\\b")
  expect_match(msg, "logical", fixed = TRUE)
  expect_match(msg, "numeric", fixed = TRUE)
})

test_that("ipw() rejects a numeric .data column where the outcome model fit a logical", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + flag, dat, wts)

  dat_num <- dat
  dat_num$flag <- as.numeric(dat$flag)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_num),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bflag\\b")
  expect_match(msg, "logical", fixed = TRUE)
})

test_that("ipw() rejects a factor .data column where the ps model fit a numeric", {
  skip_if_not_installed("deli")
  # The propensity design has its own width check, which catches a level count
  # that disagrees. The type check is more specific and reports the column.
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + num, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_fac <- dat
  dat_fac$num <- factor(dat$num)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bnum\\b")
  expect_match(msg, "wt_mod", fixed = TRUE)
})

test_that("ipw() rejects a logical .data column where the ps model fit a numeric", {
  skip_if_not_installed("deli")
  # Without a type check this reached the weight-consistency preflight, which
  # reports a disagreement about the estimand and the focal level.
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + num, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_lgl <- dat
  dat_lgl$num <- dat$num > 1

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bnum\\b")
  expect_match(msg, "logical", fixed = TRUE)
})

test_that("the mirrored type error names the column and both types", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + num, dat, wts)

  dat_fac <- dat
  dat_fac$num <- factor(dat$num)

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat_fac))
})

test_that("ipw() rejects a logical .data column where the models fit a factor", {
  skip_if_not_installed("deli")
  # Both are categorical codings, so the type check passes them and the rebuilt
  # counterfactual design is one column short of the coefficients it multiplies.
  # The outcome design has its own width check for exactly that.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_lgl <- dat
  dat_lgl$cov <- dat$cov == "a"

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "outcome_mod", fixed = TRUE)
})

test_that("ipw() keeps accepting the .data types the models were fit on", {
  skip_if_not_installed("deli")
  # The control for every rejection above: the same models, the same data, and
  # a factor covariate supplied as the factor it was fit as.
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov + num, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

# ---- the exposure column ----------------------------------------------------
#
# The exposure is the one column the rebuild overwrites, so what matters there
# is the value written rather than the column supplied. A numeric column hands
# the rebuild a number where the fit recorded a factor, and a character column
# hands it a single string, which carries no levels and collapses to one. Both
# died inside `model.matrix()`. Neither can arise without `.data`: the exposure
# then comes out of the propensity model's own frame.

test_that("ipw() rejects a numeric .data exposure where the models fit a factor", {
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf + x1, dat, wts)

  dat_num <- dat
  dat_num$zf <- dat$z

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_num),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bzf\\b")
  expect_match(msg, "factor", fixed = TRUE)
  expect_false(grepl("contrasts apply only to factors", msg, fixed = TRUE))
})

test_that("ipw() rejects a character .data exposure where the models fit a factor", {
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf + x1, dat, wts)

  dat_chr <- dat
  dat_chr$zf <- as.character(dat$z)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_chr),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bzf\\b")
  expect_false(grepl("contrasts can be applied", msg, fixed = TRUE))
})

test_that("a factor .data exposure supplied as the fitted factor still works", {
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("ipw() rejects a factor .data exposure where the models fit a numeric", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_fac <- dat
  dat_fac$z <- factor(dat$z)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\bz\\b")
})

test_that("ipw() rejects a factor .data exposure for a continuous exposure", {
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$a <- 0.5 + 0.8 * dat$x1 + rnorm(nrow(dat))
  ps_mod <- lm(a ~ x1, data = dat)
  wts <- wt_ate(
    fitted(ps_mod),
    dat$a,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  out <- lm(y ~ a + x1, data = dat, weights = wts)

  dat_fac <- dat
  dat_fac$a <- factor(round(dat$a))

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\ba\\b")
  expect_false(grepl("weights recomputed", msg, fixed = TRUE))
})

# ---- the outcome column -----------------------------------------------------
#
# A model can never be fit on a character response, so a character column under
# the fitted outcome name is always a `.data` mismatch. It reached the value
# conversion, which coerced it with an unclassed "NAs introduced by coercion"
# warning and then failed inside the M-estimator's own vocabulary; on the
# linearization path it left the standard errors non-finite instead. The guard
# sits at the shared extraction, before the conversion, so both paths reject it.

test_that("ipw() rejects a character outcome column under mestimation", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_chr <- dat
  dat_chr$y <- ifelse(dat$y == 1, "yes", "no")

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_chr, se_method = "mestimation"),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\by\\b")
  expect_match(msg, "character", fixed = TRUE)
  expect_match(msg, "numeric", fixed = TRUE)
})

test_that("ipw() rejects a character outcome column under linearization", {
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_chr <- dat
  dat_chr$y <- ifelse(dat$y == 1, "yes", "no")

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_chr, se_method = "linearization"),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\by\\b")
  expect_match(msg, "character", fixed = TRUE)
})

test_that("the character outcome error names the column and both types", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_chr <- dat
  dat_chr$y <- ifelse(dat$y == 1, "yes", "no")

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, out, .data = dat_chr, se_method = "mestimation")
  )
})

test_that("ipw() rejects a character outcome column for a categorical exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- design_data()
  dat$a <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_mod <- nnet::multinom(a ~ x1, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$a, exposure_type = "categorical")
  )
  out <- design_outcome(y ~ a + x1, dat, wts)

  dat_chr <- dat
  dat_chr$y <- ifelse(dat$y == 1, "yes", "no")

  expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_chr),
    class = "propensity_ipw_data_error"
  ))
})

test_that("ipw() rejects a factor outcome column where the models fit a numeric", {
  skip_if_not_installed("deli")
  # A factor response is converted to an indicator for its non-first levels and
  # any other response is used as its own values, so the two codings are not
  # interchangeable and a mismatch changes the outcome values silently.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_fac <- dat
  dat_fac$y <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\by\\b")
})

test_that("a logical outcome column is the numeric coding the fit used", {
  skip_if_not_installed("deli")
  # Logical and numeric responses are both converted by their own values, so
  # this direction is a working input rather than a mismatch.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_lgl <- dat
  dat_lgl$y <- as.logical(dat$y)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_lgl))
  ref <- ipw(ps_mod, out, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

# ---- the guards keep their place in the order -------------------------------
#
# The `.data` checks are ordered from the most specific diagnosis to the least,
# and the assertions run before any evaluation of the supplied frame. A type
# check that ran first would report a column type for a `.data` that is missing
# the column, or has the wrong number of rows, or belongs to a model whose frame
# cannot be reconstructed at all.

test_that("a missing .data column is still reported as the missing column", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_gone <- dat
  dat_gone$cov <- NULL

  expect_error(
    ipw(ps_mod, out, .data = dat_gone),
    class = "propensity_columns_exist_error"
  )
})

test_that("a wrong-sized .data is still reported as the row count", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_short <- dat[-1, ]
  dat_short$cov <- as.numeric(dat_short$cov)

  err <- expect_error(
    ipw(ps_mod, out, .data = dat_short),
    class = "propensity_ipw_data_error"
  )
  expect_match(design_msg(err), "one row per observation", fixed = TRUE)
})

test_that("an outcome model with no exposure term is still reported as such", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ x1, dat, wts)

  dat_bad <- dat
  dat_bad$z <- factor(dat$z)

  expect_error(
    ipw(ps_mod, out, .data = dat_bad),
    class = "propensity_ipw_exposure_error"
  )
})
