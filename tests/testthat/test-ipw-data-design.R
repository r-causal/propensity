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
# Whether the counterfactual designs can be built from it depends on which route
# builds them.
#
# Without `.data` they come from the outcome model's own model frame, which
# carries its terms attribute, so `model.matrix()` reads the `factor(z)` column
# already in it and evaluates nothing. The counterfactual column written beside
# it is dropped, both designs are the fitted design, and every contrast between
# them is zero. That route is rejected, and the remedy is to supply `.data`. It
# used to be rejected as a report that the exposure column was missing from the
# model frame, whose remedy was to supply `.data`, which then met an error of
# its own.
#
# With `.data` the design is rebuilt from the supplied frame, which re-evaluates
# the term at the counterfactual column and re-levels the result against the
# levels the fit recorded under the label. That route is accepted and is pinned,
# against the pre-converted factor model, at the end of this file.

test_that("ipw() rejects an outcome model that factors the exposure without .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "factor(z)", fixed = TRUE)
  expect_false(grepl("contrasts can be applied", msg, fixed = TRUE))
})

test_that("the factored-exposure error names the term and the route that rebuilds it", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, out))
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

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  expect_match(design_msg(err), "factor(a)", fixed = TRUE)
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
  # nowhere to write. It used to be reported as the missing column; the term is
  # now read off the formula and reported as the term it is.
  err <- expect_error(
    ipw(ps_mod, out_sq),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "I(z^2)", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
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

test_that("an exposure factored on levels the call fixes still asks for .data", {
  skip_if_not_installed("deli")
  # The levels rebuild, so the levels check passes it on to the check that asks
  # whether the variable is a column at all. The model frame holds the fitted
  # `cut()` values and no `z` column to set, so without `.data` the term is
  # reported as the term it is. It used to be reported as a missing column.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_cut <- design_outcome(y ~ cut(z, c(-1, 0.5, 2)) + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out_cut),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "cut(z, c(-1, 0.5, 2))", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
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
  # A logical column carries no levels for the rebuild to re-level against the
  # fit, so the type check rejects it whatever its width. Against a three-level
  # factor the rebuilt design is also one column short of the coefficients it
  # multiplies, which is what used to report it.
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

# ---- the levels the fit recorded --------------------------------------------
#
# Every design rebuilt from `.data` multiplies the fitted coefficients
# positionally, and for a categorical column the column a level lands in is
# fixed by the order the fit recorded rather than by the order `.data` happens
# to carry. A user who re-levels a factor for presentation after fitting, or who
# stores it as character, hands `ipw()` the same values under a different order,
# and the design rebuilt from them pairs each level with another level's
# coefficient. The fits record their levels, so the rebuilds are made to honor
# them; where honoring them is not possible, because the level order decides the
# meaning of the numbers rather than their arrangement, the mismatch is rejected.

design_level_data <- function(seed = 3025, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  # deliberately not alphabetical: a rebuild that re-levels on its own would
  # order these "high", "low", "mid" and pair every level with the wrong column
  cov <- factor(
    sample(c("low", "mid", "high"), n, replace = TRUE),
    levels = c("low", "mid", "high")
  )
  num <- sample(1:3, n, replace = TRUE)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * (cov == "mid")))
  y <- rbinom(
    n,
    1,
    plogis(-0.4 + 1.0 * z + 0.5 * x1 + 0.6 * (cov == "high") + 0.2 * num)
  )
  data.frame(x1, cov, num, z, y)
}

test_that("ipw() is unchanged when .data re-levels an outcome factor covariate", {
  skip_if_not_installed("deli")
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_relevel <- dat
  dat_relevel$cov <- factor(
    as.character(dat$cov),
    levels = c("high", "mid", "low")
  )

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_relevel))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

test_that("ipw() is unchanged when .data supplies an outcome factor covariate as character", {
  skip_if_not_installed("deli")
  # The fitted levels are not alphabetical, so a rebuild that re-levels the
  # character column on its own rebuilds a design the model was never fit on.
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_chr <- dat
  dat_chr$cov <- as.character(dat$cov)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_chr))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

test_that("ipw() is unchanged when .data supplies an outcome factor covariate as an ordered factor", {
  skip_if_not_installed("deli")
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_ord <- dat
  dat_ord$cov <- ordered(
    as.character(dat$cov),
    levels = c("low", "mid", "high")
  )

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_ord))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

test_that("ipw() is unchanged when .data declares a level the fit never saw", {
  skip_if_not_installed("deli")
  # The fits drop unused levels, so a level no observation carries contributes
  # no column to the fitted design and none to the rebuilt one either. The
  # design is the fitted design and the analysis is unchanged.
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_extra <- dat
  dat_extra$cov <- factor(
    as.character(dat$cov),
    levels = c("low", "mid", "high", "unused")
  )

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_extra))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

test_that("ipw() rejects a .data column holding a level the fit never saw", {
  skip_if_not_installed("deli")
  # A declared level with no observations is the design the fit used; an
  # observed one is not, and there is no coefficient for it. This died raw as
  # "factor cov has new levels extra", which names neither `.data` nor `ipw()`.
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_new <- dat
  levels(dat_new$cov) <- c("low", "mid", "high", "extra")
  dat_new$cov[1] <- "extra"

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_new),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bcov\\b")
  expect_match(msg, "extra", fixed = TRUE)
})

test_that("the new-level rejection covers a column read through a call", {
  skip_if_not_installed("deli")
  # `factor(num)` records its levels under the call rather than under `num`, so
  # the type check has no column to compare and the value reaches the rebuild.
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1 + factor(num), data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + factor(num), dat, wts)

  dat_new <- dat
  dat_new$num[1] <- 4L

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_new),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "factor(num)", fixed = TRUE)
  expect_match(msg, "\\b4\\b")
})

test_that("the new-level error names the column and the values the fit never saw", {
  skip_if_not_installed("deli")
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + cov, dat, wts)

  dat_new <- dat
  levels(dat_new$cov) <- c("low", "mid", "high", "extra")
  dat_new$cov[1] <- "extra"

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat_new))
})

# ---- the outcome response's levels ------------------------------------------
#
# A factor response is converted to an indicator for its non-first levels, so
# the level the fit treats as the failure has to be the first level of the
# column `.data` supplies as well. A response re-leveled the other way is the
# same values under the opposite coding: it silently reverses every estimate on
# the M-estimation path, and on the linearization path it leaves the point
# estimates alone, because they come from the outcome model's own predictions,
# and moves the standard errors, which read the response values.

design_level_response <- function(dat) {
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  dat
}

test_that("ipw() rejects a response re-leveled against the fit under mestimation", {
  skip_if_not_installed("deli")
  dat <- design_level_response(design_level_data())
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(yf ~ z + x1, dat, wts)

  dat_flip <- dat
  dat_flip$yf <- factor(as.character(dat$yf), levels = c("yes", "no"))

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_flip, se_method = "mestimation"),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\byf\\b")
  expect_match(msg, "no", fixed = TRUE)
})

test_that("ipw() rejects a response re-leveled against the fit under linearization", {
  dat <- design_level_response(design_level_data())
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(yf ~ z, dat, wts)

  dat_flip <- dat
  dat_flip$yf <- factor(as.character(dat$yf), levels = c("yes", "no"))

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_flip, se_method = "linearization"),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\byf\\b")
})

test_that("the re-leveled response error names the column and the coding", {
  skip_if_not_installed("deli")
  dat <- design_level_response(design_level_data())
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(yf ~ z + x1, dat, wts)

  dat_flip <- dat
  dat_flip$yf <- factor(as.character(dat$yf), levels = c("yes", "no"))

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat_flip))
})

test_that("a response supplied with the levels the fit recorded still works", {
  skip_if_not_installed("deli")
  dat <- design_level_response(design_level_data())
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(yf ~ z + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

# ---- the exposure column's levels -------------------------------------------
#
# The binary path takes the second sorted level of the exposure as the exposed
# group, which for a factor is the second level it declares. A `.data` exposure
# re-leveled against the fit therefore reverses the group the counterfactual
# designs and the influence functions treat as exposed. It errored, but at the
# weight-consistency preflight, whose message is about how the weights were
# built. The categorical path resolves the supplied column against the fitted
# levels instead and is unaffected.

test_that("ipw() rejects a binary factor exposure re-leveled against the fit", {
  skip_if_not_installed("deli")
  dat <- design_level_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf + x1, dat, wts)

  dat_flip <- dat
  dat_flip$zf <- factor(as.character(dat$zf), levels = c("1", "0"))

  err <- expect_error(
    ipw(ps_mod, out, .data = dat_flip),
    class = "propensity_ipw_data_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "\\bzf\\b")
  expect_false(grepl("weights recomputed", msg, fixed = TRUE))
})

test_that("the re-leveled binary exposure is rejected under linearization too", {
  dat <- design_level_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf, dat, wts)

  dat_flip <- dat
  dat_flip$zf <- factor(as.character(dat$zf), levels = c("1", "0"))

  err <- expect_error(
    ipw(ps_mod, out, .data = dat_flip, se_method = "linearization"),
    class = "propensity_ipw_data_error"
  )
  expect_match(design_msg(err), "\\bzf\\b")
})

test_that("a binary factor exposure declaring an unused level still works", {
  skip_if_not_installed("deli")
  # No observation carries the extra level, so the fitted coding is untouched.
  dat <- design_level_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ zf + x1, dat, wts)

  dat_extra <- dat
  dat_extra$zf <- factor(as.character(dat$zf), levels = c("0", "1", "2"))

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_extra))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

design_level_categorical <- function(seed = 4127, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  odds <- cbind(1, exp(0.5 * x1), exp(-0.4 * x1))
  probs <- odds / rowSums(odds)
  a <- factor(
    vapply(
      seq_len(n),
      function(i) sample(c("a", "b", "c"), 1, prob = probs[i, ]),
      character(1)
    ),
    levels = c("a", "b", "c")
  )
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.8 * (a == "b") + 1.2 * (a == "c") + 0.4 * x1)
  )
  data.frame(x1, a, y)
}

test_that("a re-leveled categorical exposure is resolved rather than rejected", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The categorical path resolves the exposure against the propensity model's
  # own level order before anything reads it, so the supplied order says
  # nothing and this stays an accepted input.
  dat <- design_level_categorical()
  ps_mod <- nnet::multinom(a ~ x1, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$a, exposure_type = "categorical")
  )
  out <- design_outcome(y ~ a + x1, dat, wts)

  dat_relevel <- dat
  dat_relevel$a <- factor(as.character(dat$a), levels = c("c", "b", "a"))

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_relevel))
  ref <- ipw(ps_mod, out, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

# ---- a factored exposure term rebuilds from .data ---------------------------
#
# `y ~ factor(z)` records the term's levels and contrast coding under the label
# "factor(z)", and the rebuild from `.data` re-evaluates that label at the
# counterfactual column and re-levels the result against the levels the fit
# recorded. The design that comes out is the design the fit used, so the
# marginal means are the ones the pre-converted factor model gives.
#
# Without `.data` there is nothing to re-evaluate: the model frame carries its
# own terms attribute, so the design is read from the `factor(z)` column already
# in it and the counterfactual column written beside it is dropped. Both
# counterfactual designs are then the fitted design, and every contrast between
# them is zero. That route keeps its rejection.

test_that("ipw() accepts a factored binary exposure term with .data", {
  skip_if_not_installed("deli")
  dat <- design_level_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  ps_fac <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_call <- design_outcome(y ~ factor(z) + x1, dat, wts)
  out_column <- design_outcome(y ~ zf + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_fac, out_column, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("a factored binary exposure term rebuilds under non-default contrasts", {
  skip_if_not_installed("deli")
  # The fit's own coding is carried into the rebuild rather than the default, so
  # a sum-coded exposure term reproduces the same marginal means.
  dat <- design_level_data()
  dat$zf <- factor(dat$z)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  ps_fac <- glm(zf ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_call <- glm(
    y ~ factor(z) + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    contrasts = list(`factor(z)` = "contr.sum"),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  out_column <- glm(
    y ~ zf + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    contrasts = list(zf = "contr.sum"),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_fac, out_column, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("ipw() accepts a factored categorical exposure term with .data", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- design_level_categorical()
  ps_mod <- nnet::multinom(a ~ x1, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$a, exposure_type = "categorical")
  )
  out_call <- design_outcome(y ~ factor(a) + x1, dat, wts)
  out_column <- design_outcome(y ~ a + x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_mod, out_column, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("the factored exposure term keeps its rejection without .data", {
  skip_if_not_installed("deli")
  # The route asymmetry: the same models are accepted with `.data` and rejected
  # without it, because only the `.data` route re-evaluates the term.
  dat <- design_level_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ factor(z) + x1, dat, wts)

  expect_no_error(ipw(ps_mod, out, .data = dat))

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "factor(z)", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
})

# ---- a term that mixes the exposure with a covariate ------------------------
#
# The rejection above is keyed on what a term evaluates to, which settles every
# term that carries levels. A term with no levels needs a different question,
# because the route decides it just as much: `y ~ z + x1 + I(z * x1)` fit on a
# numerically coded exposure rebuilds correctly from `.data` and gives the
# estimates `y ~ z * x1` gives, and without `.data` it reads `I(z * x1)` out of
# the model frame at the values it was fit on, so the counterfactual column is
# written and then ignored for that column alone. The interaction is dropped
# from the marginal means with nothing signaled, and the estimates move.
#
# A plain `*` or `:` interaction is a different construction: its model
# variables are the two plain columns, and the product is formed at design time
# from the column being set, so it survives the counterfactual write on both
# routes. What separates the two is the model frame's variables rather than the
# formula's term labels, so an exposure-reading call is caught wherever its
# terms sit while a plain interaction of the same order is accepted. Both are
# pinned here, since a guard keyed on the term label instead would either reject
# working models or let the call through inside an interaction.

test_that("ipw() rejects a call term mixing the exposure with a covariate without .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + I(z * x1), dat, wts)

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "I(z * x1)", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
})

test_that("the exposure-mixing term is accepted with .data and matches the interaction", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_call <- design_outcome(y ~ z + x1 + I(z * x1), dat, wts)
  out_star <- design_outcome(y ~ z * x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_mod, out_star, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("a star interaction with the exposure gives the same estimates on both routes", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z * x1, dat, wts)

  with_data <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  without_data <- expect_no_warning(ipw(ps_mod, out))
  expect_equal(with_data$estimates, without_data$estimates, tolerance = 1e-8)
})

test_that("a colon interaction with the exposure gives the same estimates on both routes", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + z:x1, dat, wts)

  with_data <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  without_data <- expect_no_warning(ipw(ps_mod, out))
  expect_equal(with_data$estimates, without_data$estimates, tolerance = 1e-8)
})

test_that("an interaction with a covariate that is not the exposure is left alone", {
  skip_if_not_installed("deli")
  # The guard asks whether a term reads the exposure, so a call term over other
  # covariates is untouched on either route.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + I(x1 * num), dat, wts)

  with_data <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  without_data <- expect_no_warning(ipw(ps_mod, out))
  expect_equal(with_data$estimates, without_data$estimates, tolerance = 1e-8)
})

test_that("the exposure-mixing rejection reads in the user's terms", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + I(z * x1), dat, wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, out))
})

test_that("ipw() rejects an exposure call inside an interaction without .data", {
  skip_if_not_installed("deli")
  # The term is `num:I(z * x1)`, an interaction, but the variable the frame
  # holds is `I(z * x1)`, which is as frozen at its fitted values as it is when
  # it stands alone. This ran silently and reported a risk difference the model
  # as written does not give.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + num + I(z * x1):num, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "I(z * x1)", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
})

test_that("the interacted exposure call rebuilds from .data", {
  skip_if_not_installed("deli")
  # The route that recomputes the variable at each counterfactual value gives
  # what the same interaction written over plain columns gives.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_call <- design_outcome(y ~ z + num + I(z * x1):num, dat, wts)
  out_plain <- design_outcome(y ~ z + num + z:x1:num, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_mod, out_plain, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("ipw() rejects an identity call on the exposure inside an interaction", {
  skip_if_not_installed("deli")
  # `I(z)` computes nothing at all, and the rejection is not of what the call
  # computes but of the frame column no counterfactual write reaches.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + I(z):x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out),
    class = "propensity_ipw_exposure_error"
  )
  expect_match(design_msg(err), "I(z)", fixed = TRUE)
})

test_that("the identity call on the exposure matches the plain interaction with .data", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_call <- design_outcome(y ~ z + x1 + I(z):x1, dat, wts)
  out_star <- design_outcome(y ~ z * x1, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out_call, .data = dat))
  ref <- ipw(ps_mod, out_star, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("a call over other covariates inside an interaction is left alone", {
  skip_if_not_installed("deli")
  # The question is whether the variable reads the exposure, not whether it is a
  # call, so a basis over another covariate is untouched at any order.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + poly(x1, 2):num, dat, wts)

  with_data <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  without_data <- expect_no_warning(ipw(ps_mod, out))
  expect_equal(with_data$estimates, without_data$estimates, tolerance = 1e-8)
})

# ---- a logical column where the fit recorded levels -------------------------
#
# A fit that recorded levels for a column rebuilds it by re-leveling the values
# against the ones it recorded. A logical column carries no levels to re-level,
# so it is coded by `model.matrix()` on its own terms, `FALSE` first, and the
# design lands on whichever mapping that happens to give. When the fitted factor
# has two levels the widths agree and nothing further notices: an aligned
# logical reproduces the fitted design by coincidence and a flipped one silently
# swaps the two levels' coefficients. The rejection is keyed on the type rather
# than on the values, because the values that would make it right are exactly
# the ones no check can recover once the levels are gone.

design_flag_data <- function() {
  dat <- design_data()
  dat$grp <- factor(
    ifelse(dat$num > 1, "high", "low"),
    levels = c("low", "high")
  )
  dat
}

test_that("ipw() rejects an aligned logical where the outcome model fit a two-level factor", {
  skip_if_not_installed("deli")
  dat <- design_flag_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + grp, dat, wts)

  dat_lgl <- dat
  dat_lgl$grp <- dat$grp == "high"

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bgrp\\b")
  expect_match(msg, "logical", fixed = TRUE)
  expect_match(msg, "factor", fixed = TRUE)
})

test_that("ipw() rejects a flipped logical where the outcome model fit a two-level factor", {
  skip_if_not_installed("deli")
  # The measured harm: this completed and reported estimates that differ from
  # the ones the factor column gives.
  dat <- design_flag_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + grp, dat, wts)

  dat_lgl <- dat
  dat_lgl$grp <- dat$grp == "low"

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\bgrp\\b")
})

test_that("ipw() rejects a logical where the propensity model fit a two-level factor", {
  skip_if_not_installed("deli")
  dat <- design_flag_data()
  ps_mod <- glm(z ~ x1 + grp, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_lgl <- dat
  dat_lgl$grp <- dat$grp == "high"

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_lgl),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "wt_mod", fixed = TRUE)
})

test_that("the logical-for-factor error names the column and both types", {
  skip_if_not_installed("deli")
  dat <- design_flag_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + grp, dat, wts)

  dat_lgl <- dat
  dat_lgl$grp <- dat$grp == "high"

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat_lgl))
})

test_that("a factor column supplied as character is still accepted", {
  skip_if_not_installed("deli")
  # The control for the rejection above: a character vector carries its values,
  # which the rebuild re-levels against the fit, so it rebuilds the fitted
  # design and stays accepted.
  dat <- design_flag_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + grp, dat, wts)

  dat_chr <- dat
  dat_chr$grp <- as.character(dat$grp)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_chr))
  ref <- ipw(ps_mod, out, .data = dat)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-12)
})

test_that("a logical column the models were fit on is still accepted", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + flag, dat, wts)

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

# ---- a column the fit recorded as a logical ---------------------------------
#
# The mirror of the rejection above. A fit that recorded a logical column
# recorded no levels for it either, so it took `FALSE` as the reference and the
# rebuilt design has to arrive on that coding. A factor, an ordered factor, and
# a character vector each bring a reference of their own: the widths agree, and
# which level the fitted coefficient is paired with is decided by whatever the
# column declares or sorts to. It is right by coincidence and wrong with nothing
# signaled, so all three are rejected on the type.

test_that("ipw() rejects an aligned factor where the models fit a logical", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + flag, dat, wts)

  dat_fac <- dat
  dat_fac$flag <- factor(dat$flag)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bflag\\b")
  expect_match(msg, "logical", fixed = TRUE)
  expect_match(msg, "factor", fixed = TRUE)
})

test_that("ipw() rejects a reversed factor where the models fit a logical", {
  skip_if_not_installed("deli")
  # The measured harm: this completed and reported estimates that differ from
  # the ones the logical column gives.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + flag, dat, wts)

  dat_rev <- dat
  dat_rev$flag <- factor(dat$flag, levels = c("TRUE", "FALSE"))

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_rev),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\bflag\\b")
})

test_that("ipw() rejects a character column where the models fit a logical", {
  skip_if_not_installed("deli")
  # `"FALSE"` sorts before `"TRUE"`, so this one reproduces the fitted design.
  # It is accepted only for as long as the values keep sorting that way, which
  # is nothing the column promises.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + flag, dat, wts)

  dat_chr <- dat
  dat_chr$flag <- as.character(dat$flag)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_chr),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bflag\\b")
  expect_match(msg, "character vector", fixed = TRUE)
  expect_match(msg, "logical", fixed = TRUE)
})

test_that("ipw() rejects a factor where the propensity model fit a logical", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1 + flag, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_fac <- dat
  dat_fac$flag <- factor(dat$flag)

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_fac),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "wt_mod", fixed = TRUE)
})

# ---- a response column with no levels ---------------------------------------

test_that("ipw() reports a factor outcome column that declares no levels", {
  skip_if_not_installed("deli")
  # A factor with no levels holds nothing but missing values, so the comparison
  # against the level the fit treats as the failure has no first level to read.
  # It reached that comparison and failed as a subscript out of bounds.
  dat <- design_data()
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(yf ~ z + x1, dat, wts)

  dat_empty <- dat
  dat_empty$yf <- factor(rep(NA_character_, nrow(dat)), levels = character(0))

  err <- expect_error(
    ipw(ps_mod, out, .data = dat_empty),
    class = "propensity_ipw_data_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "\\byf\\b")
  expect_false(grepl("subscript out of bounds", msg, fixed = TRUE))
})

# ---- missing values in `.data` ----------------------------------------------
#
# Every design is rebuilt with `model.frame()`, which drops a row holding a
# missing value in any column it reads, while the weights, the exposure, and the
# outcome values keep every row of `.data`. The two are then recycled against
# each other: base R warns twice that one object is not a multiple of the other,
# and the mismatch surfaces as weights that fail their consistency check, which
# reports how the weights were built when the mistake was the frame that was
# passed.

test_that("ipw() reports missing values in a .data covariate", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + num, dat, wts)

  dat_na <- dat
  dat_na$num[c(3, 17, 88)] <- NA

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_na),
    class = "propensity_ipw_data_error"
  ))
  msg <- design_msg(err)
  expect_match(msg, "\\bnum\\b")
  expect_match(msg, "3", fixed = TRUE)
  expect_false(grepl("weights recomputed", msg, fixed = TRUE))
})

test_that("ipw() reports missing values in the .data exposure", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_na <- dat
  dat_na$z[c(1, 2)] <- NA

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_na),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\bz\\b")
})

test_that("ipw() reports missing values in the .data outcome", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_na <- dat
  dat_na$y[5] <- NA

  err <- expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_na),
    class = "propensity_ipw_data_error"
  ))
  expect_match(design_msg(err), "\\by\\b")
})

test_that("ipw() reports missing values under linearization too", {
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z, dat, wts)

  dat_na <- dat
  dat_na$x1[c(4, 9)] <- NA

  expect_no_warning(expect_error(
    ipw(ps_mod, out, .data = dat_na, se_method = "linearization"),
    class = "propensity_ipw_data_error"
  ))
})

test_that("the missing-value report names the columns and their counts", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + num, dat, wts)

  dat_na <- dat
  dat_na$num[c(3, 17, 88)] <- NA
  dat_na$x1[c(3, 40)] <- NA

  expect_snapshot(error = TRUE, ipw(ps_mod, out, .data = dat_na))
})

test_that("a column .data holds but no model reads may still be missing", {
  skip_if_not_installed("deli")
  # The guard covers the columns the designs are rebuilt from, which is what
  # `model.frame()` drops rows on. A column neither model reads is never
  # consulted and is left alone.
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1, dat, wts)

  dat_na <- dat
  dat_na$flag[c(2, 6)] <- NA

  res <- expect_no_warning(ipw(ps_mod, out, .data = dat_na))
  ref <- ipw(ps_mod, out)
  expect_equal(res$estimates, ref$estimates, tolerance = 1e-8)
})

test_that("models fit on data with missing values still work without .data", {
  skip_if_not_installed("deli")
  # The model-frame route is unaffected: the fits dropped the incomplete rows
  # themselves, so every frame `ipw()` reads is already complete and aligned.
  dat <- design_data()
  dat$x1[c(3, 17, 88)] <- NA
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- glm(
    y ~ z + x1,
    data = dat[!is.na(dat$x1), ],
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  res <- expect_no_warning(ipw(ps_mod, out))
  expect_true(all(is.finite(res$estimates$std.err)))
})

# ---- the interaction hint belongs only where an interaction is involved ------
#
# The rebuild rejection closes with a note that a `*` or `:` interaction needs
# no `.data`. That is the right thing to say about `I(z * x1)`, which is an
# interaction the user wrote the long way. It is noise beside `I(z^2)`, which
# mixes in no covariate at all and leaves the reader looking for the
# interaction the message says they have.

test_that("the exposure-mixing rejection keeps the interaction hint", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out <- design_outcome(y ~ z + x1 + I(z * x1), dat, wts)

  err <- expect_error(ipw(ps_mod, out), class = "propensity_ipw_exposure_error")
  expect_match(design_msg(err), "An interaction written with", fixed = TRUE)
})

test_that("the single-variable transformation rejection drops the interaction hint", {
  skip_if_not_installed("deli")
  dat <- design_data()
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- design_weights(ps_mod)
  out_sq <- design_outcome(y ~ I(z^2) + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out_sq),
    class = "propensity_ipw_exposure_error"
  )
  msg <- design_msg(err)
  expect_match(msg, "I(z^2)", fixed = TRUE)
  expect_false(grepl("interaction", msg, fixed = TRUE))
})
