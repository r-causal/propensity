# A joint product carrying the record a version of the package that kept no
# per-component stabilization score would have written, which is also the record
# a caller who assembled the attribute by hand would hold. The weights
# themselves are untouched, so a fit built on them is the fit the caller made
# and only the record the estimator reads has changed.
strip_joint_scores <- function(wts) {
  meta <- joint_wt_meta(wts)
  meta$stabilization_score <- NULL
  attr(wts, "joint_wt_meta") <- meta

  wts
}
