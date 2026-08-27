# The `ipw` result class is owned by causalgenerics: the `new_ipw()`
# constructor, the `print()` method, and the `as.data.frame()` method all live
# upstream. propensity's `ipw()` methods must build their result with the
# upstream constructor and let the result's own methods resolve upstream, and
# the behavior a user sees must be identical whichever package supplies them.
# Shared ownership of a class is the point of the arrangement, not a violation
# of it: propensity registers the broom tidiers for the same class, and those
# must resolve here.

# ---- registration provenance ------------------------------------------------

# The namespace whose S3 methods table records registrations for each generic.
# `print` and `as.data.frame` are both in `.knownS3Generics`, so a method for
# either registers into base's table whichever package declares it; `ipw`,
# `as_marginal`, `as_conditional`, and `estimand` are causalgenerics' own
# generics and `tidy` is the generics package's, so a method for any of them
# registers into the table of the package that defines it, again whichever
# package declares the method. The six accessors all reach stats' table by one
# route or the other: `coef`, `vcov`, `confint`, and `df.residual` are named in
# `.knownS3Generics` as stats', and `nobs` and `weights` are not named there but
# resolve to closures of stats' namespace, which is the same table. The two mode
# generics and `estimand()` reach no stats table at all, because stats defines
# nothing of any of those names; looking for them there would report every
# result class as having no method rather than reporting who owns one.
# Reading these tables reports which package owns a method independently
# of what is on the search path. It does not, on its own, prove a method is
# gone: `UseMethod()` searches the environment the generic was called from as
# well as the registration table, so a method left in
# `asNamespace("propensity")` but no longer registered still serves every call
# made from inside propensity. The absence checks below therefore pair a table
# read with a namespace read. A wrong entry in this map cannot pass silently,
# because the same map serves the presence assertions further down.
ipw_generic_namespaces <- c(
  print = "base",
  as.data.frame = "base",
  ipw = "causalgenerics",
  as_marginal = "causalgenerics",
  as_conditional = "causalgenerics",
  estimand = "causalgenerics",
  tidy = "generics",
  glance = "generics",
  augment = "generics",
  coef = "stats",
  vcov = "stats",
  confint = "stats",
  nobs = "stats",
  df.residual = "stats",
  weights = "stats"
)

# Package that owns the method registered for one generic and class, or NA when
# no such method is registered.
ipw_method_owner <- function(generic, cls) {
  methods_table <- get(
    ".__S3MethodsTable__.",
    envir = asNamespace(ipw_generic_namespaces[[generic]]),
    inherits = FALSE
  )
  method <- get0(
    paste0(generic, ".", cls),
    envir = methods_table,
    inherits = FALSE
  )
  if (is.null(method)) {
    return(NA_character_)
  }
  environmentName(environment(method))
}

# Owners for a set of methods given as classes named by their generic, keyed by
# full method name so a failure names the method that moved.
ipw_method_owners <- function(methods) {
  owners <- vapply(
    seq_along(methods),
    function(i) ipw_method_owner(names(methods)[[i]], methods[[i]]),
    character(1)
  )
  stats::setNames(owners, paste0(names(methods), ".", methods))
}

# ---- constructor provenance -------------------------------------------------

# Every call to `name` inside an unevaluated expression, reported as the package
# that qualifies the call or `NA` when the call is unqualified. A missing
# argument leaves an empty symbol behind, which cannot be walked.
ipw_call_qualifiers <- function(expr, name) {
  target <- as.name(name)
  qualifiers <- character()
  walk <- function(node) {
    if (!is.call(node)) {
      return(invisible(NULL))
    }
    head <- node[[1]]
    if (identical(head, target)) {
      qualifiers <<- c(qualifiers, NA_character_)
    } else if (
      is.call(head) &&
        identical(head[[1]], quote(`::`)) &&
        identical(head[[3]], target)
    ) {
      qualifiers <<- c(qualifiers, as.character(head[[2]]))
    }
    for (i in seq_along(node)) {
      element <- node[[i]]
      if (!rlang::is_missing(element)) {
        walk(element)
      }
    }
    invisible(NULL)
  }
  walk(expr)
  qualifiers
}

# The function object a call site evaluates. Both constructors return the same
# list, so the returned object cannot say which one built it; the function
# object can. R looks an unqualified name up when the call runs, starting from
# the environment the calling function was defined in. For a propensity function
# that is propensity's namespace, whose parent is the imports frame, so reading
# the name out of `environment(fn)` reads exactly what the call site evaluates:
# propensity's own definition while one exists in the namespace, and
# causalgenerics' once only the imports frame supplies it. A `pkg::name` call
# resolves in `pkg` instead.
ipw_call_target <- function(qualifier, name, fn) {
  env <- if (is.na(qualifier)) environment(fn) else asNamespace(qualifier)
  get0(name, envir = env)
}

# The `new_ipw()` calls one propensity function makes, as function objects.
ipw_constructors_called_by <- function(fn_name) {
  fn <- get(fn_name, envir = asNamespace("propensity"))
  lapply(
    ipw_call_qualifiers(body(fn), "new_ipw"),
    ipw_call_target,
    name = "new_ipw",
    fn = fn
  )
}

# ---- what moves and what stays ----------------------------------------------

# The objects causalgenerics now owns. `new_ipw()` is exported upstream and
# `print.ipw()` and `as.data.frame.ipw()` are registered upstream.
# `format_model_call()` is internal to causalgenerics and unexported, so
# propensity cannot import it; its only caller in propensity is `print.ipw()`,
# so once that method resolves upstream nothing reaches for the helper and it
# goes with no replacement.
ipw_moved_objects <- c(
  "new_ipw",
  "print.ipw",
  "format_model_call",
  "as.data.frame.ipw"
)

# Methods of the shared result class, as classes named by their generic.
ipw_result_methods <- c(print = "ipw", as.data.frame = "ipw")

# The accessors of the shared result class, which causalgenerics owns for the
# same reason it owns `print()`: they read the fields the class contracts to
# hold and must answer identically whichever package fitted the result.
# propensity supplies the covariance they report rather than methods of its own,
# so a method of any of these appearing here would be propensity's answer
# competing with the shared one. `vcov.ipw_model` belongs to the component model
# subclass the M-estimation path wraps its models in, and is upstream for the
# same reason.
#
# `as_marginal()` and `as_conditional()` move a result between the two readings
# the accessors report, and are upstream for the same reason again: two packages
# each registering `as_conditional.ipw()` would collide in the shared method
# table, and a caller would then get whichever was installed last rather than
# the contract. propensity records the reading a result is built in through the
# `effects` field rather than through a method of its own.
#
# The pooled result class is upstream for the same reasons. `pool_ipw()` builds
# it, pools both readings, and answers for the reading it holds, so the flips
# and the estimand it reports are the pooling's own answers rather than ones a
# downstream package supplies. propensity owns the tidiers of that class and
# nothing else about it, exactly as it does for a single result.
ipw_accessor_methods <- c(
  coef = "ipw",
  vcov = "ipw",
  confint = "ipw",
  nobs = "ipw",
  df.residual = "ipw",
  weights = "ipw",
  vcov = "ipw_model",
  as_marginal = "ipw",
  as_conditional = "ipw",
  as_marginal = "ipw_pooled",
  as_conditional = "ipw_pooled",
  estimand = "ipw_pooled"
)

# Methods of the shared result class that propensity owns rather than inherits.
# causalgenerics defines no tidier, so the broom methods for the class are
# propensity's, and a user reaching them through `broom::tidy()` must land in
# propensity. Each further tidier is one entry here and one in
# `ipw_generic_namespaces`.
ipw_tidier_methods <- c(tidy = "ipw", glance = "ipw", augment = "ipw")

# `ipw()` methods that must survive. causalgenerics has no `ipw.default` of its
# own, so losing propensity's would replace a message that names the supported
# propensity score model classes with `UseMethod()`'s generic complaint.
ipw_retained_methods <- c(
  ipw = "default",
  ipw = "glm",
  ipw = "lm",
  ipw = "multinom"
)

# Internal helpers that must survive. The moved objects are interleaved with
# them in `R/ipw.R`: `new_ipw()` sits between `ipw.glm()` and
# `ipw_continuous_estimate()`, and `print.ipw()`, `format_model_call()`, and
# `as.data.frame.ipw()` form a run between `ipw.default()` and
# `calculate_estimates()`. `R/ipw.R` is a 1400-line file holding the whole
# linearization machinery, so a deletion range that overshoots by one function
# is a realistic mistake with nothing else to catch it.
ipw_retained_internals <- c(
  "ipw_continuous_estimate",
  "check_ipw_weights",
  "check_ipw_linearization_outcome",
  "calculate_estimates"
)

# ---- fixtures ---------------------------------------------------------------

sim_ipw_binary <- function(seed = 920, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x1))
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
  yc <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
  data.frame(x1, z, y, yc)
}

# A binary-exposure result on the linearization path. It needs no M-estimation,
# so the printing and coercion surfaces are pinned without a solver in the way.
fit_ipw_binary <- function(dat, outcome = c("binary", "continuous")) {
  outcome <- match.arg(outcome)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  outcome_mod <- if (identical(outcome, "binary")) {
    glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
  } else {
    lm(yc ~ z, data = dat, weights = wts)
  }
  ipw(ps_mod, outcome_mod, se_method = "linearization")
}

sim_ipw_categorical <- function(seed = 921, n = 300) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  eta_b <- -0.2 + 0.5 * x1
  eta_c <- 0.1 - 0.4 * x1
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  a <- factor(
    ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c")),
    levels = c("a", "b", "c")
  )
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  data.frame(x1, a, y)
}

fit_ipw_categorical <- function(dat) {
  ps_mod <- nnet::multinom(a ~ x1, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$a)
  wts <- wt_ate(ps, dat$a, exposure_type = "categorical")
  outcome_mod <- glm(y ~ a, data = dat, family = quasibinomial(), weights = wts)
  ipw(ps_mod, outcome_mod)
}

# ---- provenance -------------------------------------------------------------

test_that("causalgenerics supplies the ipw result class methods", {
  owners <- ipw_method_owners(ipw_result_methods)

  expect_identical(
    owners,
    c(print.ipw = "causalgenerics", as.data.frame.ipw = "causalgenerics")
  )
})

test_that("causalgenerics supplies the ipw result class accessors", {
  owners <- ipw_method_owners(ipw_accessor_methods)

  expect_identical(
    owners,
    stats::setNames(rep("causalgenerics", length(owners)), names(owners))
  )

  # The registration tables above cannot see a method defined in propensity's
  # namespace but registered nowhere, and `UseMethod()` finds one of those for
  # every call made from inside propensity. Pairing the table read with a
  # namespace read is what the moved-objects test does, for the same reason.
  survivors <- vapply(
    names(owners),
    exists,
    logical(1),
    envir = asNamespace("propensity"),
    inherits = FALSE
  )

  expect_identical(names(owners)[survivors], character())
})

test_that("no definition of the moved ipw objects survives in propensity", {
  # An unregistered method still defined in the namespace is invisible to the
  # registration tables and to `R CMD check`, which only errors the other way
  # round, on a declared-but-missing method. Every call `ipw()` makes would keep
  # reaching propensity's copy while users reach causalgenerics'.
  survivors <- vapply(
    ipw_moved_objects,
    exists,
    logical(1),
    envir = asNamespace("propensity"),
    inherits = FALSE
  )

  expect_identical(ipw_moved_objects[survivors], character())
})

test_that("nothing in propensity calls the internal format_model_call helper", {
  # causalgenerics keeps `format_model_call()` unexported, so propensity has no
  # way to import it. Deleting propensity's copy outright is only safe while
  # `print.ipw()` is its sole caller, which is what this asserts.
  namespace <- asNamespace("propensity")
  callers <- Filter(
    function(nm) {
      object <- get0(nm, envir = namespace, inherits = FALSE)
      is.function(object) &&
        length(ipw_call_qualifiers(body(object), "format_model_call")) > 0
    },
    ls(namespace, all.names = TRUE)
  )

  expect_identical(callers, character())
})

test_that("propensity builds its ipw result with causalgenerics' constructor", {
  builders <- c("ipw.glm", "ipw_continuous_estimate", "ipw.multinom")
  constructors <- lapply(builders, ipw_constructors_called_by)
  names(constructors) <- builders

  # `ipw.glm()` builds a result on each of its two standard error paths; the
  # other two builders construct once each.
  expect_identical(
    lengths(constructors),
    c(ipw.glm = 2L, ipw_continuous_estimate = 1L, ipw.multinom = 1L)
  )

  upstream <- vapply(
    constructors,
    function(targets) {
      all(vapply(targets, identical, logical(1), causalgenerics::new_ipw))
    },
    logical(1)
  )

  expect_identical(
    upstream,
    c(ipw.glm = TRUE, ipw_continuous_estimate = TRUE, ipw.multinom = TRUE)
  )
})

test_that("propensity still owns the ipw methods and helpers that stay", {
  owners <- ipw_method_owners(ipw_retained_methods)

  expect_identical(
    owners,
    stats::setNames(rep("propensity", length(owners)), names(owners))
  )

  internals <- vapply(
    ipw_retained_internals,
    exists,
    logical(1),
    envir = asNamespace("propensity"),
    inherits = FALSE
  )

  expect_identical(
    internals,
    stats::setNames(
      rep(TRUE, length(ipw_retained_internals)),
      ipw_retained_internals
    )
  )
})

test_that("propensity supplies the broom tidiers for the ipw result class", {
  owners <- ipw_method_owners(ipw_tidier_methods)

  expect_identical(
    owners,
    stats::setNames(rep("propensity", length(owners)), names(owners))
  )
})

test_that("ipw() rejects an unsupported model with propensity's error class", {
  # `ipw.default()` is the upward neighbour of the methods that move, and its
  # only current coverage is a snapshot. Snapshots report as skips under
  # `--as-cran`, so the error class needs an ordinary assertion too.
  expect_error(
    ipw(structure(list(), class = "not_a_model"), NULL),
    regexp = "not_a_model",
    class = "propensity_method_error"
  )
})

# ---- behavior that must survive the change of owner -------------------------

test_that("an ipw result prints its labelled sections and model calls", {
  # testthat 3e pins the output width but not the number of significant digits,
  # and `printCoefmat()` wraps the estimates table past the console width, which
  # splits the rows this test reads by position. The row labels name the effect
  # measure and the contrast together, so the table needs more than the 80
  # columns testthat pins.
  withr::local_options(digits = 7, width = 120)

  res <- fit_ipw_binary(sim_ipw_binary())
  expect_identical(class(res), "ipw")

  out <- capture.output(print(res))
  expect_length(out, 19L)

  expect_identical(out[[1]], "Inverse Probability Weight Estimator")
  expect_identical(out[[2]], "Estimand: ATE ")

  # The reading the result presents is a fact about the result, so it is named
  # beside the estimand as well as over the table it decides.
  expect_identical(out[[3]], "Effects: marginal (population-averaged) ")
  expect_identical(out[[4]], "")
  expect_identical(out[[5]], "Weight Estimator:")
  expect_identical(
    out[[6]],
    paste0(
      "  Call: ",
      paste(deparse(stats::getCall(res$wt_mod)), collapse = "\n"),
      " "
    )
  )
  expect_identical(out[[7]], "")
  expect_identical(out[[8]], "Outcome Model:")
  expect_identical(
    out[[9]],
    paste0(
      "  Call: ",
      paste(deparse(stats::getCall(res$outcome_mod)), collapse = "\n"),
      " "
    )
  )
  expect_identical(out[[10]], "")
  expect_identical(out[[11]], "Marginal estimates:")

  # The printed call is the model's own, not a class label.
  expect_true(grepl("glm(formula = z ~ x1", out[[6]], fixed = TRUE))
  expect_true(grepl("glm(formula = y ~ z", out[[9]], fixed = TRUE))

  # `printCoefmat()` formats the seven numeric columns and takes the effect
  # names as row labels.
  expect_identical(
    strsplit(trimws(out[[12]]), " +")[[1]],
    c(
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(sub(" .*$", "", out[seq(13, 17)]), res$estimates$effect)

  expect_identical(out[[18]], "---")
  expect_true(startsWith(out[[19]], "Signif. codes:"))
})

test_that("an ipw result prints a weighting object with no call as a class label", {
  res <- fit_ipw_binary(sim_ipw_binary())

  # A weighting object that carries no call falls back to a class label so
  # print() works for any weighting object.
  res$wt_mod <- structure(list(), class = c("weighting_object", "list"))
  expect_identical(
    capture.output(print(res))[[6]],
    "  Call: <weighting_object/list> "
  )

  # An object that cannot be subset at all reaches the same fallback, by way of
  # the guard around `stats::getCall()` rather than a `NULL` return.
  res$wt_mod <- 1
  expect_identical(capture.output(print(res))[[6]], "  Call: <numeric> ")
})

test_that("printing an ipw result returns it invisibly", {
  res <- fit_ipw_binary(sim_ipw_binary())
  withr::local_output_sink(withr::local_tempfile())

  expect_invisible(print(res))
  expect_identical(print(res), res)
})

test_that("as.data.frame() on an ipw result returns the estimates in tidy columns", {
  res <- fit_ipw_binary(sim_ipw_binary())
  df <- as.data.frame(res)
  estimates <- res$estimates

  expect_s3_class(df, "data.frame")
  expect_named(
    df,
    c("term", "contrast", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(df$term, c("mean", "mean", "rd", "log(rr)", "log(or)"))
  expect_identical(df$contrast, c("0", "1", rep("1 vs 0", 3)))

  # The values are the result's own, renamed rather than recomputed.
  expect_identical(df$estimate, estimates$estimate)
  expect_identical(df$std.error, estimates$std.err)
  expect_identical(df$statistic, estimates$z)
  expect_identical(df$p.value, estimates$p.value)

  # The bounds are not columns of the default frame: they arrive, under the
  # names the broom column conventions use, only when `conf.int` asks for them.
  # The level they were built at is an argument rather than a column, and is
  # reported by neither frame.
  bounded <- as.data.frame(res, conf.int = TRUE)
  expect_named(
    bounded,
    c(
      "term",
      "contrast",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
  expect_identical(bounded$conf.low, estimates$ci.lower)
  expect_identical(bounded$conf.high, estimates$ci.upper)

  # The covariance of the reported effects travels with the frame, since it is
  # recoverable from nothing else the frame carries.
  expect_identical(
    attr(df, "ipw_vcov", exact = TRUE),
    attr(estimates, "ipw_vcov", exact = TRUE)
  )

  # `row.names` reaches `base::as.data.frame()`.
  named <- as.data.frame(res, row.names = c("a", "b", "c", "d", "e"))
  expect_identical(rownames(named), c("a", "b", "c", "d", "e"))

  # `optional` is part of the signature `base::as.data.frame()` requires. The
  # data frame method ignores it, so it must be accepted and leave the frame
  # untouched rather than error or reshape it.
  expect_identical(as.data.frame(res, optional = TRUE), df)

  # An argument the method does not use is absorbed by the dots, not rejected
  # as unused here.
  expect_identical(as.data.frame(res, stringsAsFactors = FALSE), df)
})

test_that("exponentiating an ipw result transforms only the estimate and bounds", {
  res <- fit_ipw_binary(sim_ipw_binary())
  log_scale <- as.data.frame(res, conf.int = TRUE)
  exp_scale <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  ratio <- log_scale$term %in% c("log(rr)", "log(or)")

  # The mean rows are counterfactual risks rather than ratios, so they are
  # neither exponentiated nor relabelled.
  expect_identical(exp_scale$term, c("mean", "mean", "rd", "rr", "or"))
  expect_equal(exp_scale$estimate[ratio], exp(log_scale$estimate[ratio]))
  expect_equal(exp_scale$conf.low[ratio], exp(log_scale$conf.low[ratio]))
  expect_equal(exp_scale$conf.high[ratio], exp(log_scale$conf.high[ratio]))

  # The standard error, test statistic, and p-value stay on the log scale.
  expect_identical(exp_scale$std.error, log_scale$std.error)
  expect_identical(exp_scale$statistic, log_scale$statistic)
  expect_identical(exp_scale$p.value, log_scale$p.value)

  # The covariance describes the effects on the log scale, so a frame that has
  # moved any of them off it carries none.
  expect_false(is.null(attr(log_scale, "ipw_vcov", exact = TRUE)))
  expect_null(attr(exp_scale, "ipw_vcov", exact = TRUE))

  # The risk difference row is untouched, label included. Subsetting rows
  # carries the covariance of the whole table along with the row it kept, so it
  # is set aside here rather than compared as though it described that row.
  on_log_scale <- log_scale
  attr(on_log_scale, "ipw_vcov") <- NULL
  expect_equal(exp_scale[!ratio, ], on_log_scale[!ratio, ])
})

test_that("exponentiating a difference-only ipw result changes nothing", {
  res <- fit_ipw_binary(sim_ipw_binary(), outcome = "continuous")
  plain <- as.data.frame(res, conf.int = TRUE)
  exponentiated <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(plain$term, c("mean", "mean", "diff"))

  # No row of this table is on a log scale, so every column arrives as it was.
  # The covariance is the one thing the exponentiated frame drops, and it drops
  # it whether or not any row moved.
  expect_null(attr(exponentiated, "ipw_vcov", exact = TRUE))
  attr(plain, "ipw_vcov") <- NULL
  expect_equal(exponentiated, plain)
})

test_that("a categorical ipw result keeps its contrast column through both surfaces", {
  skip_if_not_installed("nnet")
  res <- fit_ipw_categorical(sim_ipw_categorical())

  df <- as.data.frame(res)
  expect_identical(names(df)[seq(1, 2)], c("term", "contrast"))
  expect_identical(
    df$contrast,
    c(c("a", "b", "c"), rep(c("b vs a", "c vs a"), each = 3))
  )

  exp_scale <- as.data.frame(res, exponentiate = TRUE)
  ratio <- df$term %in% c("log(rr)", "log(or)")
  expect_identical(
    exp_scale$term,
    c(rep("mean", 3), rep(c("rd", "rr", "or"), times = 2))
  )
  expect_identical(exp_scale$contrast, df$contrast)
  expect_equal(exp_scale$estimate[ratio], exp(df$estimate[ratio]))
  expect_identical(exp_scale$std.error, df$std.error)

  # print() keys the rows by effect and contrast together and drops the
  # character contrast column from the matrix it formats.
  out <- capture.output(print(res))
  expect_true(any(grepl("rd b vs a", out, fixed = TRUE)))
  expect_true(any(grepl("log(or) c vs a", out, fixed = TRUE)))
  expect_false(any(grepl("contrast", out, fixed = TRUE)))
})
