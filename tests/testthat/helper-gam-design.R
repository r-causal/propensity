# How many times a call evaluates an additive fit's design, which both the
# single-dose route and the two-model route are held to.
#
# `predict(type = "lpmatrix")` rebuilds every smooth basis over the whole data,
# and it is what an additive fit's design costs: at n = 2000 it is the larger
# part of the whole `ipw()` call. `mgcv::model.matrix.gam()` is the same call
# under another name, so counting `predict.gam()` at `type = "lpmatrix"` counts
# the design however the package asks for it.
count_gam_designs <- function(expr) {
  builds <- 0L
  bump <- function() builds <<- builds + 1L

  suppressMessages(trace(
    "predict.gam",
    tracer = bquote(if (identical(type, "lpmatrix")) .(bump)()),
    print = FALSE,
    where = asNamespace("mgcv")
  ))
  withr::defer(
    suppressMessages(untrace("predict.gam", where = asNamespace("mgcv")))
  )

  value <- expr
  list(builds = builds, value = value)
}
