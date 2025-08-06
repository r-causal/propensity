expect_propensity_error <- function(expr) {
  testthat::expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    expr
  )
}

expect_propensity_warning <- function(expr, class = NULL) {
  testthat::expect_snapshot(
    cnd_class = TRUE,
    expr
  )
}
