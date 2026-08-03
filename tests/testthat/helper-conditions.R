# Runs `expr` with the dropped-record warnings muffled and reports how many were
# raised, so a test can assert that an operation dropped a record without the
# warnings reaching the test run.
#
# The count is of warnings, not of drops. dplyr collects the warnings raised
# inside a verb and re-emits one summary of them per verb, so several drops in a
# grouped call arrive as one warning and count as one. Assert that the count is
# positive, or that it is zero, rather than that it equals a number of drops.
# Reports the name at the head of the call a condition was attributed to, which
# is what a reader sees after "Error in" or "Warning in". A condition raised
# from a dispatched method reports the generic, so the name is the one the
# caller wrote whichever method answered the call. `NA` stands for a condition
# that never arrived or that named no call at all, neither of which any
# attribution test accepts.
condition_call_name <- function(expr, classes = "error") {
  cnd <- rlang::catch_cnd(expr, classes = classes)
  if (is.null(cnd) || is.null(conditionCall(cnd))) {
    return(NA_character_)
  }

  paste(deparse(conditionCall(cnd)[[1]]), collapse = " ")
}

count_record_drops <- function(expr) {
  drops <- 0L
  count <- function(cnd) {
    drops <<- drops + 1L
    invokeRestart("muffleWarning")
  }

  # dplyr collects the warnings raised inside a verb and re-emits one plain
  # summary of them, so a drop that happens there arrives without its class and
  # is recognized by what it says instead.
  count_rewrapped <- function(cnd) {
    if (grepl("Dropping the record of which units", conditionMessage(cnd))) {
      count(cnd)
    }
  }

  value <- withCallingHandlers(
    expr,
    propensity_trim_record_warning = count,
    propensity_trunc_record_warning = count,
    warning = count_rewrapped
  )

  list(value = value, drops = drops)
}
