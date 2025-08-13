# Suppress partial match warnings and info messages during tests
op <- options(
  warnPartialMatchDollar = FALSE,
  warnPartialMatchArgs = FALSE,
  warnPartialMatchAttr = FALSE,
  propensity.quiet = TRUE
)

# Restore options on exit
withr::defer(options(op), teardown_env())
