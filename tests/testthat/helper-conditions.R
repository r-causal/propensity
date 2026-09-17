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

# Runs `expr` with the finite-variance report muffled. Stabilized normal weights
# have no finite second moment once the marginal variance of the exposure
# reaches twice the conditional one, and a fixture whose fitted means track the
# exposure closely, or whose `.sigma` is a small number chosen to exercise the
# record, sits past that boundary. The report is true of those fixtures and
# beside the point of the tests written over them, so it is muffled where it is
# not the thing being asserted.
muffle_variance_warning <- function(expr) {
  withCallingHandlers(
    expr,
    propensity_density_variance_message = function(cnd) {
      invokeRestart("muffleMessage")
    }
  )
}

# Runs `expr` with the stabilizer coverage report muffled. A fixture built to
# exercise the machinery that rebuilds a numerator inside the stacked system
# stabilizes on whatever gives that machinery something to rebuild, which is
# often a covariate the marginal structural model beside it has no reason to
# carry. The report is true of those fixtures and beside the point of the tests
# written over them, so it is muffled where it is not the thing being asserted.
# Only that class is muffled: anything else the call raises still reaches the
# test run.
muffle_coverage_warning <- function(expr) {
  withCallingHandlers(
    expr,
    propensity_ipw_stabilizer_coverage_warning = function(cnd) {
      invokeRestart("muffleWarning")
    }
  )
}
