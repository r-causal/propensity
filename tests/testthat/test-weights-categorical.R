test_that("categorical exposure detection works correctly", {
  # Factor with 3 levels
  exposure_3 <- factor(c("A", "B", "C", "A", "B", "C"))
  set.seed(123)
  ps_matrix <- matrix(runif(18), ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix) # Normalize
  colnames(ps_matrix) <- levels(exposure_3)

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    wt_ate(ps_matrix, exposure_3),
    "Treating `.exposure` as categorical"
  )

  # Character vector converted to factor
  exposure_char <- c("low", "med", "high", "low", "med", "high")
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- unique(sort(exposure_char))

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    wt_ate(ps_matrix, exposure_char),
    "Treating `.exposure` as categorical"
  )
})

test_that("categorical exposure validation works", {
  # Only 2 levels should error
  exposure_2 <- factor(c("A", "B", "A", "B"))
  ps_matrix_2 <- matrix(c(0.7, 0.3, 0.4, 0.6, 0.8, 0.2, 0.5, 0.5), ncol = 2)

  expect_propensity_error(
    wt_ate(ps_matrix_2, exposure_2, exposure_type = "categorical")
  )

  # Invalid focal category should error
  exposure_3 <- factor(c("A", "B", "C", "A"))
  ps_matrix_3 <- matrix(runif(12), ncol = 3)

  expect_propensity_error(
    wt_att(ps_matrix_3, exposure_3, .focal_level = "D")
  )
})

test_that("propensity score matrix validation works", {
  exposure <- factor(c("A", "B", "C", "A", "B"))

  # Not a matrix or data.frame
  expect_propensity_error(
    wt_ate(c(0.3, 0.4, 0.3), exposure, exposure_type = "categorical")
  )

  # Wrong number of rows
  ps_wrong_rows <- matrix(runif(9), ncol = 3)
  expect_propensity_error(
    wt_ate(ps_wrong_rows, exposure, exposure_type = "categorical")
  )

  # Wrong number of columns
  ps_wrong_cols <- matrix(runif(10), ncol = 2)
  expect_propensity_error(
    wt_ate(ps_wrong_cols, exposure, exposure_type = "categorical")
  )

  # Rows don't sum to 1
  ps_bad_sum <- matrix(
    c(
      0.3,
      0.3,
      0.3, # Sums to 0.9
      0.4,
      0.4,
      0.2,
      0.2,
      0.3,
      0.5,
      0.3,
      0.3,
      0.4,
      0.25,
      0.25,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )

  expect_propensity_error(
    wt_ate(ps_bad_sum, exposure, exposure_type = "categorical")
  )

  # Invalid probabilities
  ps_invalid <- matrix(
    c(
      0.5,
      0.6,
      -0.1, # Negative value
      0.3,
      0.3,
      0.4,
      0.2,
      0.3,
      0.5,
      0.3,
      0.3,
      0.4,
      0.25,
      0.25,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )

  expect_propensity_error(
    wt_ate(ps_invalid, exposure, exposure_type = "categorical")
  )
})

test_that("categorical propensity scores at exactly 0 or 1 are rejected", {
  exposure <- factor(c("A", "B"), levels = c("A", "B", "C"))

  # A score of exactly 0 for a unit's observed level divides by zero in the
  # weight calculation, so it must be rejected rather than produce `Inf`.
  ps_zero <- data.frame(
    A = c(0, 0.3),
    B = c(0.5, 0.4),
    C = c(0.5, 0.3)
  )

  expect_error(
    wt_ate(ps_zero, exposure, exposure_type = "categorical"),
    class = "propensity_range_error"
  )

  # A score of exactly 1 is equally outside the open interval (0, 1).
  ps_one <- data.frame(
    A = c(1, 0.3),
    B = c(0, 0.4),
    C = c(0, 0.3)
  )

  expect_error(
    wt_ate(ps_one, exposure, exposure_type = "categorical"),
    class = "propensity_range_error"
  )

  expect_error(
    wt_ato(ps_zero, exposure, exposure_type = "categorical"),
    class = "propensity_range_error"
  )

  # Scores strictly inside (0, 1) are still accepted
  ps_interior <- data.frame(
    A = c(0.2, 0.3),
    B = c(0.5, 0.4),
    C = c(0.3, 0.3)
  )

  weights <- wt_ate(ps_interior, exposure, exposure_type = "categorical")

  expect_equal(as.numeric(weights), c(1 / 0.2, 1 / 0.4))
})

test_that("the categorical range refusal names `.propensity`", {
  exposure <- factor(c("A", "B"), levels = c("A", "B", "C"))
  ps_zero <- data.frame(
    A = c(0, 0.3),
    B = c(0.5, 0.4),
    C = c(0.5, 0.3)
  )

  cnd <- expect_error(
    wt_ate(ps_zero, exposure, exposure_type = "categorical"),
    class = "propensity_range_error"
  )

  # The binary route's range refusal names the argument the scores arrived in.
  # The matrix route reports a bare range, which says nothing about which
  # argument the caller has to correct.
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(cnd)),
    "The range of values in `.propensity`",
    fixed = TRUE
  )
})

test_that("ATE weights work for categorical exposures", {
  set.seed(123)
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights <- wt_ate(ps_matrix, exposure, exposure_type = "categorical")

  # Check structure
  expect_s3_class(weights, "psw")
  expect_equal(estimand(weights), "ate")
  expect_equal(length(weights), 6)

  # Check attributes
  expect_equal(attr(weights, "n_categories"), 3)
  expect_equal(attr(weights, "category_names"), c("A", "B", "C"))

  # Check weight calculations
  # For ATE, h(e) = 1, so w_i = 1 / e_{i,Z_i}
  expected_weights <- c(
    1 / 0.5, # A: 1/0.5 = 2
    1 / 0.5, # B: 1/0.5 = 2
    1 / 0.7, # C: 1/0.7 = 1.43
    1 / 0.6, # A: 1/0.6 = 1.67
    1 / 0.4, # B: 1/0.4 = 2.5
    1 / 0.6 # C: 1/0.6 = 1.67
  )
  expect_equal(as.numeric(weights), expected_weights, tolerance = 0.01)
})

test_that("ATE stabilization works for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights_stab <- wt_ate(
    ps_matrix,
    exposure,
    exposure_type = "categorical",
    stabilize = TRUE
  )

  expect_true(is_stabilized(weights_stab))

  # Check that stabilization is applied correctly
  # Marginal probabilities: A=2/6, B=2/6, C=2/6 = 1/3 each
  expected_weights <- c(
    (1 / 3) / 0.5, # A
    (1 / 3) / 0.5, # B
    (1 / 3) / 0.7, # C
    (1 / 3) / 0.6, # A
    (1 / 3) / 0.4, # B
    (1 / 3) / 0.6 # C
  )
  expect_equal(as.numeric(weights_stab), expected_weights, tolerance = 0.01)
})

test_that("ATT weights work for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  # ATT with .focal_level = "A"
  weights_att_a <- wt_att(
    ps_matrix,
    exposure,
    .focal_level = "A",
    exposure_type = "categorical"
  )

  expect_equal(estimand(weights_att_a), "att")
  expect_equal(attr(weights_att_a, "focal_category"), "A")

  # For ATT, h(e) = e_focal
  # So w_i = e_{i,A} / e_{i,Z_i}
  expected_weights_a <- c(
    0.5 / 0.5, # A: 0.5/0.5 = 1
    0.2 / 0.5, # B: 0.2/0.5 = 0.4
    0.1 / 0.7, # C: 0.1/0.7 = 0.143
    0.6 / 0.6, # A: 0.6/0.6 = 1
    0.3 / 0.4, # B: 0.3/0.4 = 0.75
    0.2 / 0.6 # C: 0.2/0.6 = 0.333
  )
  expect_equal(as.numeric(weights_att_a), expected_weights_a, tolerance = 0.01)

  # ATT requires focal for categorical
  expect_propensity_error(
    wt_att(ps_matrix, exposure, exposure_type = "categorical")
  )
})

test_that("ATU weights work for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  # ATU with .focal_level = "A" (weights for non-A)
  weights_atu_a <- wt_atu(
    ps_matrix,
    exposure,
    .focal_level = "A",
    exposure_type = "categorical"
  )

  expect_equal(estimand(weights_atu_a), "atu")
  expect_equal(attr(weights_atu_a, "focal_category"), "A")

  # For ATU, h(e) = 1 - e_focal
  # So w_i = (1 - e_{i,A}) / e_{i,Z_i}
  expected_weights_a <- c(
    (1 - 0.5) / 0.5, # A: 0.5/0.5 = 1
    (1 - 0.2) / 0.5, # B: 0.8/0.5 = 1.6
    (1 - 0.1) / 0.7, # C: 0.9/0.7 = 1.286
    (1 - 0.6) / 0.6, # A: 0.4/0.6 = 0.667
    (1 - 0.3) / 0.4, # B: 0.7/0.4 = 1.75
    (1 - 0.2) / 0.6 # C: 0.8/0.6 = 1.333
  )
  expect_equal(as.numeric(weights_atu_a), expected_weights_a, tolerance = 0.01)
})

test_that("ATM weights work for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights_atm <- wt_atm(ps_matrix, exposure, exposure_type = "categorical")

  expect_equal(estimand(weights_atm), "atm")

  # For ATM, h(e) = min(e_1, ..., e_K)
  # So w_i = min(e_i) / e_{i,Z_i}
  expected_weights <- c(
    0.2 / 0.5, # A: min(0.5,0.3,0.2)=0.2, 0.2/0.5 = 0.4
    0.2 / 0.5, # B: min(0.2,0.5,0.3)=0.2, 0.2/0.5 = 0.4
    0.1 / 0.7, # C: min(0.1,0.2,0.7)=0.1, 0.1/0.7 = 0.143
    0.1 / 0.6, # A: min(0.6,0.3,0.1)=0.1, 0.1/0.6 = 0.167
    0.3 / 0.4, # B: min(0.3,0.4,0.3)=0.3, 0.3/0.4 = 0.75
    0.2 / 0.6 # C: min(0.2,0.2,0.6)=0.2, 0.2/0.6 = 0.333
  )
  expect_equal(as.numeric(weights_atm), expected_weights, tolerance = 0.01)
})

test_that("ATO weights work for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights_ato <- wt_ato(ps_matrix, exposure, exposure_type = "categorical")

  expect_equal(estimand(weights_ato), "ato")

  # For ATO, h(e) = 1 / sum(1/e_k) - reciprocal of harmonic mean denominator
  # So w_i = h(e_i) / e_{i,Z_i}
  h_vals <- numeric(6)
  h_vals[1] <- 1 / (1 / 0.5 + 1 / 0.3 + 1 / 0.2) # 0.115
  h_vals[2] <- 1 / (1 / 0.2 + 1 / 0.5 + 1 / 0.3) # 0.095
  h_vals[3] <- 1 / (1 / 0.1 + 1 / 0.2 + 1 / 0.7) # 0.062
  h_vals[4] <- 1 / (1 / 0.6 + 1 / 0.3 + 1 / 0.1) # 0.067
  h_vals[5] <- 1 / (1 / 0.3 + 1 / 0.4 + 1 / 0.3) # 0.109
  h_vals[6] <- 1 / (1 / 0.2 + 1 / 0.2 + 1 / 0.6) # 0.086

  expected_weights <- c(
    h_vals[1] / 0.5,
    h_vals[2] / 0.5,
    h_vals[3] / 0.7,
    h_vals[4] / 0.6,
    h_vals[5] / 0.4,
    h_vals[6] / 0.6
  )
  expect_equal(as.numeric(weights_ato), expected_weights, tolerance = 0.01)
})

test_that("Entropy weights work for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights_entropy <- wt_entropy(
    ps_matrix,
    exposure,
    exposure_type = "categorical"
  )

  expect_equal(estimand(weights_entropy), "entropy")

  # For Entropy, h(e) = -sum(e_k * log(e_k))
  # Calculate entropy for each observation
  calc_entropy <- function(probs) {
    -sum(probs * log(probs))
  }

  h_vals <- apply(ps_matrix, 1, calc_entropy)

  expected_weights <- c(
    h_vals[1] / 0.5,
    h_vals[2] / 0.5,
    h_vals[3] / 0.7,
    h_vals[4] / 0.6,
    h_vals[5] / 0.4,
    h_vals[6] / 0.6
  )
  expect_equal(as.numeric(weights_entropy), expected_weights, tolerance = 0.01)
})

test_that("data.frame input works for categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))

  # Test with plain column names
  ps_df <- data.frame(
    A = c(0.5, 0.2, 0.1, 0.6, 0.3, 0.2),
    B = c(0.3, 0.5, 0.2, 0.3, 0.4, 0.2),
    C = c(0.2, 0.3, 0.7, 0.1, 0.3, 0.6)
  )

  weights_df <- wt_ate(ps_df, exposure, exposure_type = "categorical")

  # Test with parsnip-style column names
  ps_df_parsnip <- data.frame(
    .pred_A = c(0.5, 0.2, 0.1, 0.6, 0.3, 0.2),
    .pred_B = c(0.3, 0.5, 0.2, 0.3, 0.4, 0.2),
    .pred_C = c(0.2, 0.3, 0.7, 0.1, 0.3, 0.6)
  )

  weights_df_parsnip <- wt_ate(
    ps_df_parsnip,
    exposure,
    exposure_type = "categorical"
  )

  # Both should give same results
  expect_equal(as.numeric(weights_df), as.numeric(weights_df_parsnip))

  # Compare to matrix input
  ps_matrix <- as.matrix(ps_df)
  weights_matrix <- wt_ate(ps_matrix, exposure, exposure_type = "categorical")

  expect_equal(as.numeric(weights_df), as.numeric(weights_matrix))
})

test_that("stabilization works for ATE categorical exposures", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  set.seed(123)
  ps_matrix <- matrix(runif(18), ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix) # Normalize
  colnames(ps_matrix) <- levels(exposure)

  # Test that stabilization works for ATE
  expect_no_error(
    wt_ate(ps_matrix, exposure, exposure_type = "categorical", stabilize = TRUE)
  )
})

test_that("categorical weights handle different column orders correctly", {
  set.seed(456)
  n <- 100

  # Create treatment with 3 categories
  trt <- factor(sample(c("low", "medium", "high"), n, replace = TRUE))

  # Create propensity score matrix in CORRECT order (matching factor levels)
  ps_correct <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_correct <- ps_correct / rowSums(ps_correct)
  colnames(ps_correct) <- levels(trt) # "high", "low", "medium" (alphabetical)

  # Create same matrix with columns in WRONG order
  ps_wrong <- ps_correct[, c("medium", "high", "low")]

  # Calculate weights with both matrices
  w_ate_correct <- wt_ate(ps_correct, trt, exposure_type = "categorical")
  w_ate_wrong <- wt_ate(ps_wrong, trt, exposure_type = "categorical")

  # Weights should be identical after reordering
  expect_equal(as.numeric(w_ate_correct), as.numeric(w_ate_wrong))

  # Test with ATT
  w_att_correct <- wt_att(
    ps_correct,
    trt,
    exposure_type = "categorical",
    .focal_level = "medium"
  )
  w_att_wrong <- wt_att(
    ps_wrong,
    trt,
    exposure_type = "categorical",
    .focal_level = "medium"
  )

  expect_equal(as.numeric(w_att_correct), as.numeric(w_att_wrong))

  # Verify ATT weights are correct (focal group should have weight 1)
  expect_equal(unique(as.numeric(w_att_correct[trt == "medium"])), 1)
  expect_equal(unique(as.numeric(w_att_wrong[trt == "medium"])), 1)
})

test_that("categorical weights work with parsnip-style column names", {
  set.seed(789)
  n <- 50

  trt <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # Create matrix with parsnip-style names
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c(".pred_A", ".pred_B", ".pred_C")

  # Should work without error
  w_ate <- expect_no_error(
    wt_ate(ps_matrix, trt, exposure_type = "categorical")
  )

  # Test focal matching works correctly
  w_att <- expect_no_error(
    wt_att(
      ps_matrix,
      trt,
      exposure_type = "categorical",
      .focal_level = "B"
    )
  )

  # Focal group should have weight 1
  expect_equal(unique(as.numeric(w_att[trt == "B"])), 1)
})

test_that("categorical weights error on mismatched column names", {
  n <- 50
  trt <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # Matrix with wrong column names
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("X", "Y", "Z")

  expect_propensity_error(
    wt_ate(ps_matrix, trt, exposure_type = "categorical")
  )
})

test_that("categorical weights warn on unnamed columns", {
  n <- 50
  trt <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # Matrix with no column names
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)

  # Test warning for all weight functions
  expect_propensity_warning(
    wt_ate(ps_matrix, trt, exposure_type = "categorical"),
    print = TRUE
  )

  expect_propensity_warning(
    wt_att(ps_matrix, trt, exposure_type = "categorical", .focal_level = "A"),
    print = TRUE
  )

  expect_propensity_warning(
    wt_atu(ps_matrix, trt, exposure_type = "categorical", .focal_level = "A"),
    print = TRUE
  )

  expect_propensity_warning(
    wt_atm(ps_matrix, trt, exposure_type = "categorical"),
    print = TRUE
  )

  expect_propensity_warning(
    wt_ato(ps_matrix, trt, exposure_type = "categorical"),
    print = TRUE
  )

  expect_propensity_warning(
    wt_entropy(ps_matrix, trt, exposure_type = "categorical"),
    print = TRUE
  )
})

test_that("categorical ATE weights are missing where the exposure is missing", {
  exposure <- factor(c("A", "B", NA, "C"), levels = c("A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.2,
      0.5,
      0.3,
      0.3,
      0.4,
      0.3,
      0.5,
      0.2,
      0.3,
      0.3,
      0.4,
      0.3
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights <- wt_ate(ps_matrix, exposure, exposure_type = "categorical")

  # For ATE, h(e) = 1, so w_i = 1 / e_{i,Z_i}. The third unit has no observed
  # level, so it has no weight.
  expected_weights <- c(
    1 / 0.2, # A
    1 / 0.4, # B
    NA, # missing exposure
    1 / 0.3 # C
  )

  expect_equal(as.numeric(weights), expected_weights)
  expect_equal(which(is.na(weights)), 3L)
  expect_true(all(is.finite(as.numeric(weights)[-3])))
})

test_that("stabilized categorical ATE marginals use the complete cases", {
  exposure <- factor(c("A", "B", NA, "C"), levels = c("A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.2,
      0.5,
      0.3,
      0.3,
      0.4,
      0.3,
      0.5,
      0.2,
      0.3,
      0.3,
      0.4,
      0.3
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights <- wt_ate(
    ps_matrix,
    exposure,
    exposure_type = "categorical",
    stabilize = TRUE
  )

  expect_true(is_stabilized(weights))

  # The marginal probabilities are shares of the observed exposures, so they
  # sum to 1 over the three units with a level: one A, one B, and one C give
  # 1/3 each. Dividing by the number of rows instead would deflate every
  # weight by the share of units with a missing exposure.
  p_marginal <- 1 / 3

  expected_weights <- c(
    p_marginal / 0.2, # A
    p_marginal / 0.4, # B
    NA, # missing exposure
    p_marginal / 0.3 # C
  )

  expect_equal(as.numeric(weights), expected_weights)
  expect_equal(which(is.na(weights)), 3L)
  expect_true(all(is.finite(as.numeric(weights)[-3])))
})

test_that("stabilized categorical weights refuse an exposure with nothing observed", {
  exposure <- factor(rep(NA_character_, 4), levels = c("A", "B", "C"))
  ps_matrix <- matrix(
    1 / 3,
    nrow = 4,
    ncol = 3,
    dimnames = list(NULL, levels(exposure))
  )

  # The marginal probabilities the default stabilizer builds are shares of the
  # units with an observed level, so with none observed every share is 0 / 0.
  # Weights that are missing everywhere are not an answer to the call that was
  # made, so the degenerate stabilizer is reported rather than returned.
  expect_error(
    wt_ate(
      ps_matrix,
      exposure,
      exposure_type = "categorical",
      stabilize = TRUE
    ),
    class = "propensity_stabilize_categorical_error"
  )

  # Without stabilization there is no marginal to form, and a unit with no
  # observed level has no weight, which is already the answer.
  weights <- wt_ate(ps_matrix, exposure, exposure_type = "categorical")
  expect_true(all(is.na(as.numeric(weights))))
})

test_that("tilted categorical weights are missing where the exposure is missing", {
  exposure <- factor(c("A", "B", NA, "C"), levels = c("A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.2,
      0.5,
      0.3,
      0.3,
      0.4,
      0.3,
      0.5,
      0.2,
      0.3,
      0.3,
      0.4,
      0.3
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  weights_ato <- wt_ato(ps_matrix, exposure, exposure_type = "categorical")

  # For ATO, h(e) = 1 / sum(1 / e_k), so w_i = h(e_i) / e_{i,Z_i}
  h_vals <- numeric(4)
  h_vals[1] <- 1 / (1 / 0.2 + 1 / 0.5 + 1 / 0.3)
  h_vals[2] <- 1 / (1 / 0.3 + 1 / 0.4 + 1 / 0.3)
  h_vals[3] <- 1 / (1 / 0.5 + 1 / 0.2 + 1 / 0.3)
  h_vals[4] <- 1 / (1 / 0.3 + 1 / 0.4 + 1 / 0.3)

  expected_ato <- c(
    h_vals[1] / 0.2, # A
    h_vals[2] / 0.4, # B
    NA, # missing exposure
    h_vals[4] / 0.3 # C
  )

  expect_equal(as.numeric(weights_ato), expected_ato)
  expect_equal(which(is.na(weights_ato)), 3L)

  weights_att <- wt_att(
    ps_matrix,
    exposure,
    .focal_level = "A",
    exposure_type = "categorical"
  )

  # For ATT, h(e) = e_focal, so w_i = e_{i,A} / e_{i,Z_i}
  expected_att <- c(
    0.2 / 0.2, # A
    0.3 / 0.4, # B
    NA, # missing exposure
    0.3 / 0.3 # C
  )

  expect_equal(as.numeric(weights_att), expected_att)
  expect_equal(which(is.na(weights_att)), 3L)
})
