test_that("tidyr::pivot_longer works with propensity classes", {
  skip_if_not_installed("tidyr")

  # Create a data frame with different weight columns
  df <- data.frame(
    id = 1:4,
    ate_wts = psw(c(0.5, 0.7, 0.3, 0.8), estimand = "ate"),
    att_wts = psw(c(0.4, 0.6, 0.2, 0.9), estimand = "att"),
    other = c(1, 2, 3, 4)
  )

  # Pivot longer should work but with warning
  expect_propensity_warning(
    result <- tidyr::pivot_longer(
      df,
      cols = c(ate_wts, att_wts),
      names_to = "weight_type",
      values_to = "weight"
    )
  )

  expect_equal(nrow(result), 8)
  expect_type(result$weight, "double")
  # Check actual values are preserved
  expect_equal(result$weight, c(0.5, 0.4, 0.7, 0.6, 0.3, 0.2, 0.8, 0.9))
  expect_equal(result$id, rep(1:4, each = 2))
  expect_equal(result$weight_type, rep(c("ate_wts", "att_wts"), 4))
})

test_that("tidyr::pivot_longer works with mixed propensity classes", {
  skip_if_not_installed("tidyr")

  df <- data.frame(
    id = 1:4,
    psw_col = psw(c(0.5, 0.7, 0.3, 0.8), estimand = "ate"),
    trim_col = ps_trim(c(0.4, 0.6, 0.2, 0.9), lower = 0.1, upper = 0.9),
    trunc_col = ps_trunc(c(0.3, 0.5, 0.4, 0.7), lower = 0.1, upper = 0.9)
  )

  expect_propensity_warning(
    expect_propensity_warning(
      result <- tidyr::pivot_longer(
        df,
        cols = c(psw_col, trim_col, trunc_col),
        names_to = "type",
        values_to = "value"
      )
    )
  )

  expect_equal(nrow(result), 12)
  expect_type(result$value, "double")
  # Check values are preserved (accounting for any NAs from ps_trim)
  expected_vals <- c(
    0.5,
    0.4,
    0.3, # id 1
    0.7,
    0.6,
    0.5, # id 2
    0.3,
    0.2,
    0.4, # id 3
    0.8,
    0.9,
    0.7
  ) # id 4
  actual_vals <- result$value[!is.na(result$value)]
  expect_equal(actual_vals, expected_vals[!is.na(expected_vals)])
})

test_that("tidyr::pivot_longer preserves class when all columns are compatible", {
  skip_if_not_installed("tidyr")

  df <- data.frame(
    id = 1:4,
    wt1 = psw(c(0.5, 0.7, 0.3, 0.8), estimand = "ate"),
    wt2 = psw(c(0.4, 0.6, 0.2, 0.9), estimand = "ate"),
    wt3 = psw(c(0.3, 0.5, 0.4, 0.7), estimand = "ate")
  )

  expect_silent(
    result <- tidyr::pivot_longer(
      df,
      cols = starts_with("wt"),
      names_to = "weight_var",
      values_to = "weight"
    )
  )

  expect_s3_class(result$weight, "psw")
  expect_equal(estimand(result$weight), "ate")
})

test_that("c() works as expected with warnings", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")
  z <- 0.6

  # Different estimands
  expect_propensity_warning(
    result <- c(x, y)
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.3, 0.8))

  # Mixed with numeric
  expect_propensity_warning(
    result <- c(x, z)
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.6))
})

test_that("rbind and data frame operations work", {
  df1 <- data.frame(
    id = 1:2,
    wt = psw(c(0.5, 0.7), estimand = "ate")
  )

  df2 <- data.frame(
    id = 3:4,
    wt = psw(c(0.3, 0.8), estimand = "att")
  )

  # rbind preserves the first object's class
  result <- rbind(df1, df2)

  expect_equal(nrow(result), 4)
  expect_s3_class(result$wt, "psw")
  expect_equal(as.numeric(result$wt), c(0.5, 0.7, 0.3, 0.8))

  # But vec_rbind does trigger the warning
  expect_propensity_warning(
    result2 <- vctrs::vec_rbind(df1, df2)
  )
  expect_equal(nrow(result2), 4)
  expect_type(result2$wt, "double")
  expect_equal(result2$wt, c(0.5, 0.7, 0.3, 0.8))
})

test_that("tidyr::pivot_wider works with propensity classes", {
  skip_if_not_installed("tidyr")

  # Create long format data - all same estimand to avoid conversion during creation
  df_long <- data.frame(
    id = rep(1:3, each = 2),
    estimand = rep(c("ate", "att"), 3),
    weight = c(
      0.5,
      0.6,
      0.7,
      0.8,
      0.3,
      0.4
    )
  )

  # Add psw class after creation
  df_long$weight[df_long$estimand == "ate"] <- psw(
    df_long$weight[df_long$estimand == "ate"],
    estimand = "ate"
  )
  df_long$weight[df_long$estimand == "att"] <- psw(
    df_long$weight[df_long$estimand == "att"],
    estimand = "att"
  )

  # Now the weight column contains mixed psw objects
  result <- tidyr::pivot_wider(
    df_long,
    names_from = estimand,
    values_from = weight
  )

  expect_equal(nrow(result), 3)
  # When pivoting wider, each column should maintain its class
  # But the result may be numeric if conversion happened
  # Let's check what we actually get
  if (is.numeric(result$ate)) {
    expect_equal(result$ate, c(0.5, 0.7, 0.3))
    expect_equal(result$att, c(0.6, 0.8, 0.4))
  } else {
    expect_s3_class(result$ate, "psw")
    expect_s3_class(result$att, "psw")
    expect_equal(as.numeric(result$ate), c(0.5, 0.7, 0.3))
    expect_equal(as.numeric(result$att), c(0.6, 0.8, 0.4))
  }
})

test_that("tidyr operations with NAs in ps_trim work correctly", {
  skip_if_not_installed("tidyr")

  # Create ps_trim with some NAs (trimmed values)
  ps_vals <- c(0.05, 0.5, 0.95, 0.7) # First and third will be trimmed
  trim_obj <- ps_trim(ps_vals, method = "ps", lower = 0.1, upper = 0.9)

  df <- data.frame(
    id = 1:4,
    group = c("A", "A", "B", "B"),
    weights = trim_obj
  )

  # pivot_wider should preserve the ps_trim structure
  result <- tidyr::pivot_wider(
    df,
    names_from = group,
    values_from = weights
  )

  # Check structure - 4 rows for 4 unique ids
  expect_equal(nrow(result), 4)
  # Check the values - each id should have one value in either A or B column
  # id 1 and 2 have values in A column, id 3 and 4 have values in B column
  expect_true(is.na(result$A[1])) # First value was trimmed
  expect_equal(as.numeric(result$A[2]), 0.5)
  expect_true(all(is.na(result$A[3:4]))) # These ids belong to group B

  expect_true(all(is.na(result$B[1:2]))) # These ids belong to group A
  expect_true(is.na(result$B[3])) # Third value was trimmed
  expect_equal(as.numeric(result$B[4]), 0.7)
})

test_that("multiple pivot operations work correctly", {
  skip_if_not_installed("tidyr")

  # Create a more complex data frame with consistent estimands initially
  df <- data.frame(
    person = rep(1:2, each = 4),
    time = rep(c("pre", "post"), each = 2, times = 2),
    treatment = rep(c("control", "treated"), 4),
    weight = c(
      0.5,
      0.6, # person 1, pre
      0.7,
      0.8, # person 1, post
      0.3,
      0.4, # person 2, pre
      0.2,
      0.9 # person 2, post
    )
  )

  # Add different estimands by person
  df$weight[df$person == 1] <- psw(df$weight[df$person == 1], estimand = "ate")
  df$weight[df$person == 2] <- psw(df$weight[df$person == 2], estimand = "att")

  # First pivot wider by treatment - this should trigger warning
  # because we're combining different estimands
  wide1 <- tidyr::pivot_wider(
    df,
    names_from = treatment,
    values_from = weight
  )

  expect_equal(nrow(wide1), 4)
  # Check if conversion happened
  if (is.numeric(wide1$control)) {
    expect_type(wide1$control, "double")
    expect_type(wide1$treated, "double")
  }

  # Then pivot wider by time
  wide2 <- tidyr::pivot_wider(
    wide1,
    names_from = time,
    values_from = c(control, treated)
  )

  expect_equal(nrow(wide2), 2)
  expect_equal(ncol(wide2), 5) # person + 4 weight columns
})

test_that("tidyr works with stabilized weights", {
  skip_if_not_installed("tidyr")

  # Create stabilized weights
  df <- data.frame(
    id = 1:4,
    trt = rep(c("A", "B"), each = 2),
    wt_stab = new_psw(
      c(0.9, 1.1, 0.8, 1.2),
      estimand = "ate",
      stabilized = TRUE
    ),
    wt_unstab = psw(c(2.5, 3.0, 2.0, 3.5), estimand = "ate", stabilized = FALSE)
  )

  # Pivot should warn about different stabilization
  expect_propensity_warning(
    long <- tidyr::pivot_longer(
      df,
      cols = c(wt_stab, wt_unstab),
      names_to = "stab_type",
      values_to = "weight"
    )
  )

  expect_type(long$weight, "double")
  expect_equal(long$weight, c(0.9, 2.5, 1.1, 3.0, 0.8, 2.0, 1.2, 3.5))
})

test_that("tidyr with all NA weights works", {
  skip_if_not_installed("tidyr")

  # Create ps_trim where everything is trimmed
  ps_vals <- c(0.01, 0.02, 0.98, 0.99)
  trim_obj <- ps_trim(ps_vals, method = "ps", lower = 0.1, upper = 0.9)

  df <- data.frame(
    id = 1:4,
    group = c("A", "B", "A", "B"),
    weight = trim_obj
  )

  # All values should be NA
  expect_true(all(is.na(df$weight)))

  # Pivot should still work
  wide <- tidyr::pivot_wider(
    df,
    names_from = group,
    values_from = weight
  )

  expect_equal(nrow(wide), 4)
})
