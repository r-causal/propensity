be_quiet <- function() {
  getOption("propensity.quiet", default = FALSE)
}

alert_info <- function(.message) {
  if (!be_quiet()) {
    cli::cli_alert_info(text = .message)
  }
}
