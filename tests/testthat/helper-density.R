# The maximum likelihood scale of a Student's t with `df` degrees of freedom
# fit to `residuals`, for the tests that hold the package's own estimate against
# an oracle.
#
# The score for the scale, multiplied through by the scale so that it is a mean
# of bounded terms, is `(df + 1) * r^2 / (df * s^2 + r^2) - 1` at each residual.
# It falls from `df` at a scale of nothing to -1 as the scale grows without
# bound, so the root is bracketed either side of the root mean square and found
# rather than searched for: a search over the likelihood settles to about the
# fourth root of the machine epsilon, which is fewer digits than these tests
# read the two against each other to.
t_scale_mle <- function(residuals, df) {
  residuals <- residuals[!is.na(residuals)]
  rms <- sqrt(mean(residuals^2))

  score <- function(s) {
    mean((df + 1) * residuals^2 / (df * s^2 + residuals^2)) - 1
  }

  stats::uniroot(
    score,
    interval = c(rms / 1e4, rms * 50),
    tol = .Machine$double.eps^0.75
  )$root
}
