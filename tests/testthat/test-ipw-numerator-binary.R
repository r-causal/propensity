# `ipw()` with a binary exposure whose weights were stabilized on a fitted
# numerator model. The continuous route already estimates such a model in the
# stacked system, so the sandwich accounts for the numerator having been fitted
# rather than reading it as a constant. A binary exposure's numerator is a
# conditional probability rather than a conditional density, and the block that
# estimates it is the score of the binomial fit that reported it, but everything
# else about the arrangement is the same: the numerator parameters sit in the
# stabilization block, in place of the single marginal proportion the default
# stabilizer seeds there.
#
# The oracle for every point estimate is the weighted marginal structural model's
# own coefficients, and the oracle for the standard error is that it differs from
# the one the same weights get when the numerator is handed over as a fixed
# `stabilization_score`: the two routes weight the units identically, so anything
# that separates them is the sandwich reading the numerator differently.

# ---- data simulator ---------------------------------------------------------

# A binary exposure, a modifier that changes its effect, and a covariate the
# propensity score model reads and the numerator does not.
#
# The covariate is what the numerator has to be estimated against: it is the
# only thing the reported model reads that the numerator does not, so it is the
# only direction in which the numerator's own estimation can reach the reported
# standard errors at all. See the note above the fits below. It confounds the
# exposure strongly enough that the weights vary over it by more than rounding.
sim_binary_numerator <- function(seed = 6621, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  high <- as.numeric(v == "1")
  z <- rbinom(n, 1, plogis(1.5 * x1 - 0.6 * high))
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.7 * z + 0.5 * x1 - 0.4 * high + 0.9 * z * high)
  )
  data.frame(x1, v, z, y)
}

# ---- model fitting ----------------------------------------------------------

# The three fits every test below compares. All three weight the same propensity
# score model; what separates them is the numerator. `model` stabilizes on a
# fitted binomial model of the exposure on the modifier, `score` stabilizes on
# the very same numerator handed over as a vector of numbers, and `marginal`
# takes the default stabilizer. The first two build one and the same set of
# weights, so anything that differs between them is the sandwich rather than the
# weights.
#
# The reported model reads the covariate as well as the modifier, and that is
# what makes the numerator's estimation visible at all. A numerator of the
# exposure on the modifier is constant within each cell of the modifier and the
# exposure, so a model saturated in those two cells fits the very same
# coefficients weighted or unweighted by it: the numerator divides out of every
# cell mean, the estimator does not move when the numerator does, and no
# variance can describe an estimation the estimate does not depend on. The
# covariate is the direction the cells do not span, so it is where estimating
# the numerator has anything to say. The saturated case is pinned as its own
# test at the end of this file.
binary_numerator_fits <- function(dat, outcome_rhs = "z * v + x1") {
  ps_mod <- glm(z ~ x1 + v, data = dat, family = binomial())
  num_mod <- glm(z ~ v, data = dat, family = binomial())
  p <- as.numeric(fitted(num_mod))
  score <- dat$z * p + (1 - dat$z) * (1 - p)

  fit_one <- function(stabilize, stab_score) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(
        ps_mod,
        stabilize = stabilize,
        stabilization_score = stab_score
      )
    )
    outcome_mod <- glm(
      stats::reformulate(outcome_rhs, response = "y"),
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )

    list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
  }

  list(
    num_mod = num_mod,
    score = score,
    model = fit_one(num_mod, NULL),
    score_fit = fit_one(TRUE, score),
    marginal = fit_one(TRUE, NULL)
  )
}

# ---- the estimates a stacked binary numerator model reports -----------------

test_that("ipw() reports the same estimate from a binary numerator model and its score", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  # The two routes weight the units identically, which is what makes the
  # standard errors below comparable at all.
  expect_equal(
    as.numeric(fits$model$wts),
    as.numeric(fits$score_fit$wts),
    tolerance = 1e-12
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  expect_s3_class(res_model, "ipw")
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res_model$estimates$std.err)))
})

test_that("ipw() stacks the binary numerator model it was given", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)
  res_marginal <- ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod)

  # A supplied score is a constant of the stacked system and a supplied model is
  # not: the model's coefficients are parameters solved alongside everything
  # else. The default stabilizer estimates one parameter, the marginal
  # proportion; a numerator model replaces it with its own block rather than
  # sitting beside it.
  theta_model <- coef(res_model$fit)
  theta_score <- coef(res_score$fit)
  theta_marginal <- coef(res_marginal$fit)

  stab_names <- paste0("stab_", names(coef(fits$num_mod)))
  expect_true(all(stab_names %in% names(theta_model)))
  expect_false(any(grepl("^stab_", names(theta_score))))
  expect_identical(
    names(theta_marginal)[grepl("^stab_", names(theta_marginal))],
    "stab_pi"
  )
  expect_identical(
    length(theta_model),
    length(theta_score) + length(coef(fits$num_mod))
  )

  # The block is solved at the model it was given rather than carried at
  # whatever the seed was.
  expect_equal(
    unname(theta_model[stab_names]),
    unname(coef(fits$num_mod)),
    tolerance = 1e-6
  )
})

test_that("a stacked binary numerator model reports a different standard error", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  se_model <- res_model$estimates$std.err
  se_score <- res_score$estimates$std.err

  expect_true(all(is.finite(se_model) & se_model > 0))

  # The weights are the same weights, so the difference is the numerator having
  # been estimated rather than known: the score route treats it as a constant
  # and the model route carries its estimation into the sandwich.
  expect_false(isTRUE(all.equal(se_model, se_score, tolerance = 1e-6)))
})

# ---- effect modification end to end -----------------------------------------

test_that("a binary numerator model is the numerator `.by` asks for", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  # A fit reporting effects within the strata of a modifier reports a default
  # stabilizer that conditions on nothing, because a numerator conditioning on
  # the modifier would be consistent for the same stratum effects and tighter.
  # This numerator model is a model of the exposure on that modifier, which is
  # the numerator the report asks for, so the model route is silent.
  expect_warning(
    ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod, .by = v),
    class = "propensity_ipw_by_stabilizer_warning"
  )
  expect_no_warning(
    res <- ipw(fits$model$ps_mod, fits$model$outcome_mod, .by = v),
    class = "propensity_ipw_by_stabilizer_warning"
  )

  expect_s3_class(res, "ipw")
  expect_true("group" %in% names(res$estimates))
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("the stratum effects a binary numerator model reports are the score's", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  # The two routes build the same weights, so they fit the same outcome model
  # and report the same stratum effects and the same effect-modification rows.
  # Only the standard errors move, and they move because the numerator was
  # estimated rather than assumed.
  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod, .by = v)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod, .by = v)

  expect_identical(res_model$estimates$group, res_score$estimates$group)
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
  expect_false(isTRUE(all.equal(
    res_model$estimates$std.err,
    res_score$estimates$std.err,
    tolerance = 1e-6
  )))
})

# ---- what a saturated model reads a numerator as ----------------------------

test_that("a marginal structural model saturated in the numerator reads it as none", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat, outcome_rhs = "z * v")

  # A numerator of the exposure on the modifier takes one value in each cell of
  # the modifier and the exposure, so a model saturated in those cells fits the
  # same coefficients with it and without it: a constant within a cell divides
  # out of that cell's weighted mean. The estimator does not move when the
  # numerator does, so nothing about having estimated the numerator can reach
  # the standard error, and the model route reports exactly what the score
  # route reports.
  unstabilized <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fits$model$ps_mod)
  )
  saturated <- glm(
    y ~ z * v,
    data = dat,
    family = quasibinomial(),
    weights = unstabilized,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_equal(
    unname(coef(fits$model$outcome_mod)),
    unname(coef(saturated)),
    tolerance = 1e-8
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  expect_equal(
    res_model$estimates$std.err,
    res_score$estimates$std.err,
    tolerance = 1e-6
  )
})

# ---- a binary numerator design rebuilt from `.data` -------------------------
#
# A numerator model's design is one of the designs `ipw()` rebuilds when the
# caller supplies `.data`, and it is rebuilt on the same terms as the others:
# over the rows every model read, under the coding the fit recorded, and out of
# `.data` rather than out of a frame the fit may no longer keep. Read off the
# fit's own frame instead, it is a design over other rows than everything it is
# stacked with, and a column `.data` supplies as another type never reaches it
# at all.

# Whitespace in a cli-formatted message wraps where the console is narrow, so
# the message is flattened before anything is matched in it.
binary_numerator_ipw_message <- function(cnd) {
  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

# The binary fixtures above with a covariate only the outcome model reads, one
# of whose values is missing. The three fits then read one frame and keep
# different rows of it: the outcome model drops the incomplete row and the
# propensity score and numerator models, which never read the column, keep it.
# A `.data` holding the frame all three were given therefore has a row to drop
# before any design is built over it.
sim_binary_numerator_gap <- function(seed = 6621, n = 600) {
  dat <- sim_binary_numerator(seed = seed, n = n)
  dat$w <- rev(dat$x1)
  dat$w[[7]] <- NA

  dat
}

binary_numerator_gap_fits <- function(dat) {
  ps_mod <- glm(z ~ x1 + v, data = dat, family = binomial())
  num_mod <- glm(z ~ v, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod, stabilize = num_mod)
  )
  outcome_mod <- glm(
    y ~ z * v + x1 + w,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, num_mod = num_mod, outcome_mod = outcome_mod)
}

test_that("a binary numerator design is restricted to the rows .data keeps", {
  dat <- sim_binary_numerator_gap()
  fits <- binary_numerator_gap_fits(dat)
  kept <- !is.na(dat$w)

  # Supplying the frame the fits were given and supplying the rows `ipw()`
  # restricts it to are the same request, so they report the same thing. The
  # numerator's design is one of the designs that restriction is for.
  res_given <- ipw(fits$ps_mod, fits$outcome_mod, .data = dat)
  res_kept <- ipw(fits$ps_mod, fits$outcome_mod, .data = dat[kept, ])

  expect_s3_class(res_given, "ipw")
  expect_equal(
    res_given$estimates$estimate,
    res_kept$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    res_given$estimates$std.err,
    res_kept$estimates$std.err,
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res_given$estimates$std.err)))
})

test_that("a stacked binary numerator block solves over the rows .data keeps", {
  dat <- sim_binary_numerator_gap()
  fits <- binary_numerator_gap_fits(dat)
  kept <- !is.na(dat$w)

  # The block the numerator contributes is its own binomial score, which is
  # exactly identified, so the rows its design carries are the rows its solved
  # coefficients are the fit over. The numerator arrived fit to the whole frame
  # and the system reads it over the rows the outcome model kept, so what it
  # solves to is the refit on those rows rather than the coefficients it came
  # with. The two differ by more than the tolerance, which is what makes the pin
  # say anything.
  res <- ipw(fits$ps_mod, fits$outcome_mod, .data = dat)
  refit <- glm(z ~ v, data = dat[kept, ], family = binomial())

  stab_names <- paste0("stab_", names(coef(refit)))
  theta <- coef(res$fit)

  expect_true(all(stab_names %in% names(theta)))
  expect_equal(
    unname(theta[stab_names]),
    unname(coef(refit)),
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(
    unname(coef(refit)),
    unname(coef(fits$num_mod)),
    tolerance = 1e-6
  )))
})

# ---- a numerator covariate `.data` supplies as another type ------------------
#
# A numerator model is rebuilt from `.data` the way the other models are, so a
# column it alone reads is a column the rebuild can be given as the wrong type,
# or not given at all. The guards the other models meet report the column, the
# argument the model arrived in, and both types.

# A numerator of the exposure on one covariate, which nothing else reads.
binary_numerator_x2_fits <- function(dat) {
  ps_mod <- glm(z ~ x1 + v, data = dat, family = binomial())
  num_mod <- glm(z ~ x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod, stabilize = num_mod)
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, num_mod = num_mod, outcome_mod = outcome_mod)
}

# A three-level factor over the same rows, which nothing models. It is what a
# numeric numerator covariate is supplied as below, and what one is fit on when
# the mistake runs the other way.
binary_numerator_grouping <- function(dat) {
  factor(c("a", "b", "c")[1 + (rank(dat$x1) %% 3)], levels = c("a", "b", "c"))
}

test_that("a binary numerator covariate supplied as a factor where the fit read a number is refused", {
  dat <- sim_binary_numerator()
  dat$g <- binary_numerator_grouping(dat)
  dat$x2 <- round(dat$x1, 2)

  fits <- binary_numerator_x2_fits(dat)

  # The type sweep compares the class the fit recorded for the column with the
  # class `.data` supplies and refuses the pair before any design is rebuilt.
  # What it heads off is a factor of three levels taking two design columns
  # where the number it stands in for took one.
  supplied <- dat
  supplied$x2 <- dat$g

  err <- expect_error(
    ipw(fits$ps_mod, fits$outcome_mod, .data = supplied),
    class = "propensity_ipw_data_error"
  )

  message <- binary_numerator_ipw_message(err)
  expect_match(message, "x2", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a binary numerator covariate supplied as a number where the fit read a factor is refused", {
  dat <- sim_binary_numerator()
  dat$g <- binary_numerator_grouping(dat)
  dat$x2 <- dat$g

  fits <- binary_numerator_x2_fits(dat)

  # The other direction, which never reaches a width to compare: the rebuild
  # asks for the fit's contrast coding of a column that arrives with no levels
  # to code.
  supplied <- dat
  supplied$x2 <- as.numeric(dat$g)

  err <- expect_error(
    ipw(fits$ps_mod, fits$outcome_mod, .data = supplied),
    class = "propensity_ipw_data_error"
  )

  message <- binary_numerator_ipw_message(err)
  expect_match(message, "x2", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a binary numerator covariate absent from .data is refused", {
  dat <- sim_binary_numerator()
  dat$x2 <- round(dat$x1, 2)

  fits <- binary_numerator_x2_fits(dat)

  # The columns the rebuilds read are asked for before any of them runs, and a
  # numerator model's covariates are among them. Left out of the set, a column
  # only the numerator reads reaches `model.matrix()` as an object that is not
  # there.
  supplied <- dat
  supplied$x2 <- NULL

  err <- expect_error(
    ipw(fits$ps_mod, fits$outcome_mod, .data = supplied),
    class = "propensity_columns_exist_error"
  )

  expect_match(binary_numerator_ipw_message(err), "x2", fixed = TRUE)
})

# ---- recovering the data behind a binary numerator model --------------------
#
# A `glm` keeps its model frame, so its design is usually there to be read. One
# fit with `model = FALSE` keeps none, and the design is recovered by
# re-evaluating the fitting call; a fit made inside a wrapper whose frame is
# gone cannot be re-evaluated. What that owes the caller is the request the
# denominator's recovery already makes: name the model that cannot be rebuilt
# and ask for `.data`.

# The formula is written in the calling frame and the fitting call names a
# variable that lives only inside the wrapper, so nothing can rebuild the
# numerator's design once the wrapper has returned. Its fitted probabilities are
# still readable, and the weights were built from them.
binary_numerator_frame_gone <- function(dat) {
  fmla <- z ~ v
  fit_in_function <- function(fitting_data) {
    glm(fmla, data = fitting_data, family = binomial(), model = FALSE)
  }

  fit_in_function(dat)
}

binary_numerator_gone_fits <- function(dat, numerator) {
  ps_mod <- glm(z ~ x1 + v, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod, stabilize = numerator)
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

test_that("a binary numerator model whose data is gone asks for .data", {
  dat <- sim_binary_numerator()
  gone <- binary_numerator_frame_gone(dat)

  expect_error(model.matrix(gone))

  fits <- binary_numerator_gone_fits(dat, gone)

  err <- expect_error(
    ipw(fits$ps_mod, fits$outcome_mod),
    class = "propensity_ipw_data_error"
  )

  # The propensity score model here is rebuildable, so a message naming it would
  # send the caller to the wrong fit.
  message <- binary_numerator_ipw_message(err)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a binary numerator model whose data is gone is rebuilt from .data", {
  dat <- sim_binary_numerator()
  gone <- binary_numerator_frame_gone(dat)

  # The other half of the contract: the `.data` the refusal above asks for
  # rebuilds the numerator's design from the fit's own terms and contrasts, and
  # what comes back is what the same numerator reports when its frame was never
  # lost.
  kept <- glm(z ~ v, data = dat, family = binomial())
  expect_equal(unname(coef(gone)), unname(coef(kept)), tolerance = 1e-10)

  gone_fits <- binary_numerator_gone_fits(dat, gone)
  kept_fits <- binary_numerator_gone_fits(dat, kept)

  res_data <- ipw(gone_fits$ps_mod, gone_fits$outcome_mod, .data = dat)
  res_recovered <- ipw(kept_fits$ps_mod, kept_fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_recovered$estimates$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_data$estimates$std.err,
    res_recovered$estimates$std.err,
    tolerance = 1e-6
  )
})
