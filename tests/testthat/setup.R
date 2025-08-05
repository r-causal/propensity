# Suppress partial match warnings during tests
op <- options(
  warnPartialMatchDollar = FALSE,
  warnPartialMatchArgs = FALSE,
  warnPartialMatchAttr = FALSE
)

# Restore options on exit
withr::defer(options(op), teardown_env())
