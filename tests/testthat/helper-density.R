# The maximum likelihood scale of a Student's t with `df` degrees of freedom
# fit to `residuals`, for the tests that hold the package's own estimate against
# an oracle.
#
# The oracle maximizes the likelihood itself rather than solving the equation
# the package solves, and it searches the log of the scale over a range wide
# enough to hold a fit that all but interpolates its exposure, so that neither
# the bracket the package brackets the root with nor the equation it reads the
# root off of is assumed here. A search settles to about the fourth root of the
# machine epsilon, which is fewer digits than these tests read the two against
# each other to, so the maximum is then polished by Newton steps on the score
# `(df + 1) * r^2 / (df * s^2 + r^2) - 1` averaged over the residuals, which is
# the derivative of the log likelihood multiplied through by the scale.
t_scale_mle <- function(residuals, df) {
  residuals <- residuals[!is.na(residuals)]
  rms <- sqrt(mean(residuals^2))

  loglik <- function(log_scale) {
    sum(stats::dt(residuals / exp(log_scale), df = df, log = TRUE) - log_scale)
  }

  log_scale <- stats::optimize(
    loglik,
    interval = log(rms) + c(-40, 5),
    maximum = TRUE,
    tol = .Machine$double.eps^0.5
  )$maximum

  score <- function(scale) {
    mean((df + 1) * residuals^2 / (df * scale^2 + residuals^2)) - 1
  }

  # The score with respect to the log of the scale, which is what the steps are
  # taken on so that no step can land on a scale of nothing or less.
  slope <- function(scale) {
    mean(
      -2 *
        df *
        (df + 1) *
        scale^2 *
        residuals^2 /
        (df * scale^2 + residuals^2)^2
    )
  }

  for (i in seq_len(3)) {
    scale <- exp(log_scale)
    log_scale <- log_scale - score(scale) / slope(scale)
  }

  exp(log_scale)
}
