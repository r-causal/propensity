# The `ipw` result class is owned by causalgenerics: the `new_ipw()`
# constructor, the `print()` method, and the `as.data.frame()` method all live
# upstream. propensity's `ipw()` methods must build their result with the
# upstream constructor and let the result's own methods resolve upstream, and
# the behavior a user sees must be identical whichever package supplies them.

# ---- registration provenance ------------------------------------------------

# The namespace whose S3 methods table records registrations for each generic.
# `print` and `as.data.frame` are both in `.knownS3Generics`, so a method for
# either registers into base's table whichever package declares it; `ipw` is
# causalgenerics' own generic. Reading these tables reports which package owns a
# method independently of what is on the search path. It does not, on its own,
# prove a method is gone: `UseMethod()` searches the environment the generic was
# called from as well as the registration table, so a method left in
# `asNamespace("propensity")` but no longer registered still serves every call
# made from inside propensity. The absence checks below therefore pair a table
# read with a namespace read. A wrong entry in this map cannot pass silently,
# because the same map serves the presence assertions further down.
ipw_generic_namespaces <- c(
  print = "base",
  as.data.frame = "base",
  ipw = "causalgenerics"
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
  # and `printCoefmat()` wraps the estimates table past 80 columns under a
  # larger `digits`, which splits the rows this test reads by position.
  withr::local_options(digits = 7)

  res <- fit_ipw_binary(sim_ipw_binary())
  expect_identical(class(res), "ipw")

  out <- capture.output(print(res))
  expect_length(out, 16L)

  expect_identical(out[[1]], "Inverse Probability Weight Estimator")
  expect_identical(out[[2]], "Estimand: ATE ")
  expect_identical(out[[3]], "")
  expect_identical(out[[4]], "Propensity Score Model:")
  expect_identical(
    out[[5]],
    paste0(
      "  Call: ",
      paste(deparse(stats::getCall(res$ps_mod)), collapse = "\n"),
      " "
    )
  )
  expect_identical(out[[6]], "")
  expect_identical(out[[7]], "Outcome Model:")
  expect_identical(
    out[[8]],
    paste0(
      "  Call: ",
      paste(deparse(stats::getCall(res$outcome_mod)), collapse = "\n"),
      " "
    )
  )
  expect_identical(out[[9]], "")
  expect_identical(out[[10]], "Estimates:")

  # The printed call is the model's own, not a class label.
  expect_true(grepl("glm(formula = z ~ x1", out[[5]], fixed = TRUE))
  expect_true(grepl("glm(formula = y ~ z", out[[8]], fixed = TRUE))

  # `printCoefmat()` formats the seven numeric columns and takes the effect
  # names as row labels.
  expect_identical(
    strsplit(trimws(out[[11]]), " +")[[1]],
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
  expect_identical(sub(" .*$", "", out[seq(12, 14)]), res$estimates$effect)

  expect_identical(out[[15]], "---")
  expect_true(startsWith(out[[16]], "Signif. codes:"))
})

test_that("an ipw result prints a weighting object with no call as a class label", {
  res <- fit_ipw_binary(sim_ipw_binary())

  # A weighting object that carries no call falls back to a class label so
  # print() works for any weighting object.
  res$ps_mod <- structure(list(), class = c("weighting_object", "list"))
  expect_identical(
    capture.output(print(res))[[5]],
    "  Call: <weighting_object/list> "
  )

  # An object that cannot be subset at all reaches the same fallback, by way of
  # the guard around `stats::getCall()` rather than a `NULL` return.
  res$ps_mod <- 1
  expect_identical(capture.output(print(res))[[5]], "  Call: <numeric> ")
})

test_that("printing an ipw result returns it invisibly", {
  res <- fit_ipw_binary(sim_ipw_binary())
  withr::local_output_sink(withr::local_tempfile())

  expect_invisible(print(res))
  expect_identical(print(res), res)
})

test_that("as.data.frame() on an ipw result returns the estimates component", {
  res <- fit_ipw_binary(sim_ipw_binary())
  df <- as.data.frame(res)

  expect_s3_class(df, "data.frame")
  expect_named(
    df,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_equal(df, res$estimates)
  expect_identical(df$effect, c("rd", "log(rr)", "log(or)"))

  # `row.names` reaches `base::as.data.frame()`.
  named <- as.data.frame(res, row.names = c("a", "b", "c"))
  expect_identical(rownames(named), c("a", "b", "c"))

  # `optional` is part of the signature `base::as.data.frame()` requires. The
  # data frame method ignores it, so it must be accepted and leave the
  # estimates untouched rather than error or reshape them.
  expect_equal(as.data.frame(res, optional = TRUE), res$estimates)

  # `...` passes through to `base::as.data.frame()`. An argument the data frame
  # method does not use is absorbed there, not rejected as unused here.
  expect_equal(as.data.frame(res, stringsAsFactors = FALSE), res$estimates)
})

test_that("exponentiating an ipw result transforms only the estimate and bounds", {
  res <- fit_ipw_binary(sim_ipw_binary())
  log_scale <- as.data.frame(res)
  exp_scale <- as.data.frame(res, exponentiate = TRUE)

  ratio <- log_scale$effect %in% c("log(rr)", "log(or)")

  expect_identical(exp_scale$effect, c("rd", "rr", "or"))
  expect_equal(exp_scale$estimate[ratio], exp(log_scale$estimate[ratio]))
  expect_equal(exp_scale$ci.lower[ratio], exp(log_scale$ci.lower[ratio]))
  expect_equal(exp_scale$ci.upper[ratio], exp(log_scale$ci.upper[ratio]))

  # The standard error, z statistic, and p-value stay on the log scale.
  expect_identical(exp_scale$std.err, log_scale$std.err)
  expect_identical(exp_scale$z, log_scale$z)
  expect_identical(exp_scale$p.value, log_scale$p.value)
  expect_identical(exp_scale$conf.level, log_scale$conf.level)

  # The risk difference row is untouched, label included.
  expect_equal(exp_scale[!ratio, ], log_scale[!ratio, ])
})

test_that("exponentiating a difference-only ipw result changes nothing", {
  res <- fit_ipw_binary(sim_ipw_binary(), outcome = "continuous")

  expect_identical(as.data.frame(res)$effect, "diff")
  expect_equal(as.data.frame(res, exponentiate = TRUE), as.data.frame(res))
})

test_that("a categorical ipw result keeps its comparison column through both surfaces", {
  skip_if_not_installed("nnet")
  res <- fit_ipw_categorical(sim_ipw_categorical())

  df <- as.data.frame(res)
  expect_identical(names(df)[seq(1, 2)], c("effect", "comparison"))
  expect_identical(df$comparison, rep(c("b vs a", "c vs a"), each = 3))

  exp_scale <- as.data.frame(res, exponentiate = TRUE)
  ratio <- df$effect %in% c("log(rr)", "log(or)")
  expect_identical(exp_scale$effect, rep(c("rd", "rr", "or"), times = 2))
  expect_identical(exp_scale$comparison, df$comparison)
  expect_equal(exp_scale$estimate[ratio], exp(df$estimate[ratio]))
  expect_identical(exp_scale$std.err, df$std.err)

  # print() keys the rows by effect and comparison together and drops the
  # character comparison column from the matrix it formats.
  out <- capture.output(print(res))
  expect_true(any(grepl("rd b vs a", out, fixed = TRUE)))
  expect_true(any(grepl("log(or) c vs a", out, fixed = TRUE)))
  expect_false(any(grepl("comparison", out, fixed = TRUE)))
})
