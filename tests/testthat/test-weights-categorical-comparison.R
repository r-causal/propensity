# Tests comparing propensity package categorical weights to WeightIt and PSweight

test_that("categorical weights match WeightIt for all estimands", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("nnet")
  skip_on_cran()

  # Note: WeightIt uses its own multinom_weightit() implementation for multi-category
  # treatments while we use nnet::multinom. Small differences are expected due to
  # different optimization algorithms.

  # Create test dataset with 3-category treatment
  set.seed(123)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)

  # Generate treatment with known propensity scores
  z1 <- 0.5 * x1 + 0.3 * x2 + 0.5 * x3 + rnorm(n, sd = 0.5)
  z2 <- -0.3 * x1 + 0.5 * x2 - 0.3 * x3 + rnorm(n, sd = 0.5)
  p1 <- exp(z1) / (1 + exp(z1) + exp(z2))
  p2 <- exp(z2) / (1 + exp(z1) + exp(z2))
  p3 <- 1 - p1 - p2

  # Sample treatment
  u <- runif(n)
  trt <- factor(ifelse(u < p1, "A", ifelse(u < p1 + p2, "B", "C")))

  # Create data frame
  test_data <- data.frame(
    trt = trt,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  # Fit propensity score model
  ps_model <- nnet::multinom(
    trt ~ x1 + x2 + x3,
    data = test_data,
    trace = FALSE
  )
  ps_matrix <- predict(ps_model, type = "probs")

  # Ensure matrix has columns for all levels in correct order
  ps_matrix <- ps_matrix[, levels(trt)]

  # Test ATE weights
  w_ate_propensity <- wt_ate(ps_matrix, trt, exposure_type = "categorical")
  w_ate_weightit <- WeightIt::weightit(
    trt ~ x1 + x2 + x3,
    data = test_data,
    method = "ps",
    estimand = "ATE"
  )$weights

  expect_equal(as.numeric(w_ate_propensity), w_ate_weightit, tolerance = 1e-5)

  # Test ATT weights for each focal category
  for (focal in levels(trt)) {
    w_att_propensity <- wt_att(
      ps_matrix,
      trt,
      focal = focal,
      exposure_type = "categorical"
    )
    w_att_weightit <- WeightIt::weightit(
      trt ~ x1 + x2 + x3,
      data = test_data,
      method = "ps",
      estimand = "ATT",
      focal = focal
    )$weights

    expect_equal(
      as.numeric(w_att_propensity),
      w_att_weightit,
      tolerance = 1e-5,
      label = paste("ATT weights with focal =", focal)
    )
  }

  # Test ATO weights
  w_ato_propensity <- wt_ato(ps_matrix, trt, exposure_type = "categorical")
  w_ato_weightit <- WeightIt::weightit(
    trt ~ x1 + x2 + x3,
    data = test_data,
    method = "ps",
    estimand = "ATO"
  )$weights

  expect_equal(as.numeric(w_ato_propensity), w_ato_weightit, tolerance = 1e-5)

  # Test ATM weights
  w_atm_propensity <- wt_atm(ps_matrix, trt, exposure_type = "categorical")
  w_atm_weightit <- WeightIt::weightit(
    trt ~ x1 + x2 + x3,
    data = test_data,
    method = "ps",
    estimand = "ATM"
  )$weights

  expect_equal(as.numeric(w_atm_propensity), w_atm_weightit, tolerance = 1e-5)
})

test_that("categorical weights produce same treatment effect estimates as PSweight", {
  skip_if_not_installed("PSweight")
  skip_if_not_installed("nnet")
  skip_on_cran()

  # Create test dataset with 3-category treatment
  set.seed(456)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)

  # Generate treatment
  z1 <- 0.5 * x1 + 0.3 * x2 + 0.5 * x3 + rnorm(n, sd = 0.5)
  z2 <- -0.3 * x1 + 0.5 * x2 - 0.3 * x3 + rnorm(n, sd = 0.5)
  p1 <- exp(z1) / (1 + exp(z1) + exp(z2))
  p2 <- exp(z2) / (1 + exp(z1) + exp(z2))
  p3 <- 1 - p1 - p2

  u <- runif(n)
  trt <- ifelse(u < p1, 1, ifelse(u < p1 + p2, 2, 3))

  # Generate outcome for PSweight
  y <- 2 + 1 * (trt == 1) + 3 * (trt == 2) + x1 + x2 + x3 + rnorm(n)

  # Create data frame - PSweight needs numeric treatment
  test_data <- data.frame(
    y = y,
    trt = trt,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  # Get propensity scores using our standard approach
  ps_model <- nnet::multinom(
    trt ~ x1 + x2 + x3,
    data = test_data,
    trace = FALSE
  )
  ps_matrix <- predict(ps_model, type = "probs")
  trt_factor <- factor(trt)

  # Test ATE estimation
  w_ate_propensity <- wt_ate(
    ps_matrix,
    trt_factor,
    exposure_type = "categorical"
  )

  # Calculate potential outcomes using our weights
  po_our <- numeric(3)
  for (g in 1:3) {
    po_our[g] <- sum(y[trt == g] * as.numeric(w_ate_propensity)[trt == g]) /
      sum(as.numeric(w_ate_propensity)[trt == g])
  }

  # PSweight ATE estimation
  psw_ate <- PSweight::PSweight(
    ps.formula = trt ~ x1 + x2 + x3,
    yname = "y",
    data = test_data,
    weight = "IPW"
  )

  # Compare potential outcome estimates
  expect_equal(po_our, as.numeric(psw_ate$muhat), tolerance = 1e-3)

  # Test entropy weights
  w_entropy_propensity <- wt_entropy(
    ps_matrix,
    trt_factor,
    exposure_type = "categorical"
  )

  # Calculate potential outcomes using entropy weights
  po_entropy_our <- numeric(3)
  for (g in 1:3) {
    po_entropy_our[g] <- sum(
      y[trt == g] * as.numeric(w_entropy_propensity)[trt == g]
    ) /
      sum(as.numeric(w_entropy_propensity)[trt == g])
  }

  psw_entropy <- PSweight::PSweight(
    ps.formula = trt ~ x1 + x2 + x3,
    yname = "y",
    data = test_data,
    weight = "entropy"
  )

  expect_equal(po_entropy_our, as.numeric(psw_entropy$muhat), tolerance = 1e-3)

  # Test overlap (ATO) weights
  w_ato_propensity <- wt_ato(
    ps_matrix,
    trt_factor,
    exposure_type = "categorical"
  )

  # Calculate potential outcomes using ATO weights
  po_ato_our <- numeric(3)
  for (g in 1:3) {
    po_ato_our[g] <- sum(y[trt == g] * as.numeric(w_ato_propensity)[trt == g]) /
      sum(as.numeric(w_ato_propensity)[trt == g])
  }

  psw_overlap <- PSweight::PSweight(
    ps.formula = trt ~ x1 + x2 + x3,
    yname = "y",
    data = test_data,
    weight = "overlap"
  )

  expect_equal(po_ato_our, as.numeric(psw_overlap$muhat), tolerance = 1e-3)
})

test_that("categorical weights handle 4+ categories correctly", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("nnet")
  skip_on_cran()

  # Create test dataset with 4-category treatment
  set.seed(789)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Generate 4-category treatment
  z1 <- 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
  z2 <- -0.3 * x1 + 0.5 * x2 + rnorm(n, sd = 0.5)
  z3 <- 0.2 * x1 - 0.4 * x2 + rnorm(n, sd = 0.5)

  denom <- 1 + exp(z1) + exp(z2) + exp(z3)
  p1 <- exp(z1) / denom
  p2 <- exp(z2) / denom
  p3 <- exp(z3) / denom
  p4 <- 1 - p1 - p2 - p3

  # Sample treatment
  u <- runif(n)
  trt <- factor(
    ifelse(
      u < p1,
      "A",
      ifelse(u < p1 + p2, "B", ifelse(u < p1 + p2 + p3, "C", "D"))
    )
  )

  test_data <- data.frame(trt = trt, x1 = x1, x2 = x2)

  # Fit propensity score model
  ps_model <- nnet::multinom(trt ~ x1 + x2, data = test_data, trace = FALSE)
  ps_matrix <- predict(ps_model, type = "probs")
  ps_matrix <- ps_matrix[, levels(trt)]

  # Test ATE weights
  w_ate_propensity <- wt_ate(ps_matrix, trt, exposure_type = "categorical")
  w_ate_weightit <- WeightIt::weightit(
    trt ~ x1 + x2,
    data = test_data,
    method = "ps",
    estimand = "ATE"
  )$weights

  expect_equal(as.numeric(w_ate_propensity), w_ate_weightit, tolerance = 1e-5)

  # Check that weights have correct attributes
  expect_equal(attr(w_ate_propensity, "n_categories"), 4)
  expect_equal(attr(w_ate_propensity, "category_names"), levels(trt))
})

test_that("stabilized categorical ATE weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("nnet")
  skip_on_cran()

  set.seed(101)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Generate 3-category treatment
  z1 <- 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
  z2 <- -0.3 * x1 + 0.5 * x2 + rnorm(n, sd = 0.5)
  p1 <- exp(z1) / (1 + exp(z1) + exp(z2))
  p2 <- exp(z2) / (1 + exp(z1) + exp(z2))
  p3 <- 1 - p1 - p2

  u <- runif(n)
  trt <- factor(ifelse(u < p1, "A", ifelse(u < p1 + p2, "B", "C")))

  test_data <- data.frame(trt = trt, x1 = x1, x2 = x2)

  # Fit propensity score model
  ps_model <- nnet::multinom(trt ~ x1 + x2, data = test_data, trace = FALSE)
  ps_matrix <- predict(ps_model, type = "probs")
  ps_matrix <- ps_matrix[, levels(trt)]

  # Test stabilized ATE weights
  w_ate_stab_propensity <- wt_ate(
    ps_matrix,
    trt,
    exposure_type = "categorical",
    stabilize = TRUE
  )

  # WeightIt stabilized weights
  # Note: WeightIt uses s.weights for stabilization
  w_obj <- WeightIt::weightit(
    trt ~ x1 + x2,
    data = test_data,
    method = "ps",
    estimand = "ATE",
    stabilize = TRUE
  )
  w_ate_stab_weightit <- w_obj$weights

  expect_equal(
    as.numeric(w_ate_stab_propensity),
    w_ate_stab_weightit,
    tolerance = 1e-5
  )
  expect_true(is_stabilized(w_ate_stab_propensity))
})

test_that("categorical weights work with parsnip models", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("nnet")
  skip_on_cran()

  set.seed(202)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)

  # Generate 3-category treatment
  z1 <- 0.5 * x1 + 0.3 * x2 + 0.5 * x3 + rnorm(n, sd = 0.5)
  z2 <- -0.3 * x1 + 0.5 * x2 - 0.3 * x3 + rnorm(n, sd = 0.5)
  p1 <- exp(z1) / (1 + exp(z1) + exp(z2))
  p2 <- exp(z2) / (1 + exp(z1) + exp(z2))
  p3 <- 1 - p1 - p2

  u <- runif(n)
  trt <- factor(ifelse(u < p1, "A", ifelse(u < p1 + p2, "B", "C")))

  test_data <- data.frame(
    trt = trt,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  # Fit propensity score model using parsnip
  ps_spec <- parsnip::multinom_reg() %>%
    parsnip::set_engine("nnet") %>%
    parsnip::set_mode("classification")

  ps_fit <- ps_spec %>%
    parsnip::fit(trt ~ x1 + x2 + x3, data = test_data)

  # Get predictions as probabilities (data frame)
  ps_probs <- predict(ps_fit, new_data = test_data, type = "prob")

  # Test that we can calculate weights directly with the data frame
  expect_no_error(
    w_ate <- wt_ate(ps_probs, trt, exposure_type = "categorical")
  )
  expect_no_error(
    w_att <- wt_att(ps_probs, trt, focal = "A", exposure_type = "categorical")
  )
  expect_no_error(
    w_ato <- wt_ato(ps_probs, trt, exposure_type = "categorical")
  )

  # Check basic properties
  expect_s3_class(w_ate, "psw")
  expect_equal(length(w_ate), n)
  expect_equal(attr(w_ate, "n_categories"), 3)

  # Compare to direct nnet fit
  ps_nnet <- nnet::multinom(trt ~ x1 + x2 + x3, data = test_data, trace = FALSE)
  ps_matrix_nnet <- predict(ps_nnet, type = "probs")
  ps_matrix_nnet <- ps_matrix_nnet[, levels(trt)]

  w_ate_nnet <- wt_ate(ps_matrix_nnet, trt, exposure_type = "categorical")

  # Should be very close (may have slight differences due to fitting algorithm)
  expect_equal(as.numeric(w_ate), as.numeric(w_ate_nnet), tolerance = 0.01)
})
