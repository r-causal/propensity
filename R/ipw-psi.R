# Internal building blocks for the M-estimation path of ipw(). Three targets:
# ipw_weight_fn() (the weight registry), ipw_theta_layout() (the theta
# partition and its root-valued seed), and build_ipw_psi() (the stacked
# estimating-function builder). Everything here is unexported and follows the
# M-estimation design contract: the psi matrix stacks blocks in the fixed order
# [ps score | stabilization | outcome score | marginal means | contrasts], and
# the weights that enter the outcome-score block are recomputed from the
# ps-parameter block of theta on every evaluation so the sandwich variance
# accounts for propensity score estimation.

# Inverse link functions. Covers the propensity score links (logit, probit,
# cloglog) and the outcome links (logit, identity, log) used when reconstructing
# fitted values from a coefficient block.
ipw_inv_link <- function(link) {
  switch(
    link,
    logit = stats::plogis,
    probit = stats::pnorm,
    cloglog = function(x) 1 - exp(-exp(x)),
    identity = function(x) x,
    log = exp,
    abort(
      "Unsupported link {.val {link}}.",
      error_class = "propensity_ipw_link_error"
    )
  )
}

# ---- weight registry --------------------------------------------------------

ipw_weight_fn <- function(exposure_type, estimand) {
  supported <- switch(
    exposure_type,
    binary = c("ate", "att", "atu", "atm", "ato", "entropy"),
    categorical = c("ate", "att", "atu", "atm", "ato", "entropy"),
    continuous = "ate",
    abort(
      "Unsupported exposure type {.val {exposure_type}}.",
      error_class = "propensity_ipw_estimand_error"
    )
  )

  if (!estimand %in% supported) {
    abort(
      c(
        "{.fun ipw} does not support the {.val {estimand}} estimand for \\
        {.val {exposure_type}} exposures.",
        i = "Supported estimands for {.val {exposure_type}} exposures: \\
        {.val {supported}}."
      ),
      error_class = "propensity_ipw_estimand_error"
    )
  }

  switch(
    exposure_type,
    binary = ipw_binary_weight_fn(estimand),
    categorical = ipw_categorical_weight_fn(estimand),
    continuous = ipw_continuous_weight_fn(estimand)
  )
}

# Binary weight registry. `ps` is the vector of propensity scores e_i,
# `exposure` is the 0/1 exposure z, and `extras` carries the stabilization
# probability and the fixed stabilization score. Formulas mirror the binary
# helpers in R/weights.R exactly.
ipw_binary_weight_fn <- function(estimand) {
  switch(
    estimand,
    ate = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      if (!is.null(extras$score)) {
        wt <- (z / e) + ((1 - z) / (1 - e))
        wt * extras$score
      } else if (!is.null(extras$stab_prob)) {
        p1 <- extras$stab_prob
        p0 <- 1 - p1
        (z * p1 / e) + ((1 - z) * p0 / (1 - e))
      } else {
        (z / e) + ((1 - z) / (1 - e))
      }
    },
    att = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      ((e * z) / e) + ((e * (1 - z)) / (1 - e))
    },
    atu = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      (((1 - e) * z) / e) + (((1 - e) * (1 - z)) / (1 - e))
    },
    atm = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      pmin(e, 1 - e) / (z * e + (1 - z) * (1 - e))
    },
    ato = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      (1 - e) * z + e * (1 - z)
    },
    entropy = function(ps, exposure, extras) {
      z <- exposure
      e <- ps
      h_e <- -e * log(e) - (1 - e) * log(1 - e)
      h_e / (z * e + (1 - z) * (1 - e))
    }
  )
}

# Categorical weight registry. `ps` is the n-by-K propensity score matrix,
# `exposure` the n-by-K reference-first indicator, and `extras` carries the
# focal column index, the length-K stabilization probabilities (column order),
# and the fixed stabilization score. Formulas mirror
# calculate_categorical_weights() in R/weights.R.
ipw_categorical_weight_fn <- function(estimand) {
  tilt <- switch(
    estimand,
    ate = function(ps, extras) rep(1, nrow(ps)),
    att = function(ps, extras) ps[, extras$focal_idx],
    atu = function(ps, extras) 1 - ps[, extras$focal_idx],
    ato = function(ps, extras) 1 / rowSums(1 / ps),
    atm = function(ps, extras) {
      do.call(pmin, lapply(seq_len(ncol(ps)), function(j) ps[, j]))
    },
    entropy = function(ps, extras) {
      ps_safe <- ps
      ps_safe[ps == 0] <- .Machine$double.eps
      -rowSums(ps_safe * log(ps_safe))
    }
  )

  function(ps, exposure, extras) {
    e_actual <- rowSums(exposure * ps)
    h_e <- tilt(ps, extras)
    weights <- h_e / e_actual

    if (estimand == "ate") {
      if (!is.null(extras$score)) {
        weights <- weights * extras$score
      } else if (!is.null(extras$stab_probs)) {
        stab_row <- matrix(
          extras$stab_probs,
          nrow = nrow(ps),
          ncol = ncol(ps),
          byrow = TRUE
        )
        weights <- weights * rowSums(exposure * stab_row)
      }
    }

    weights
  }
}

# Continuous weight registry (ate only). `ps` is the fitted mean X alpha,
# `exposure` the continuous A, and `extras` carries the conditional variance
# sigma2_d, the marginal moments mu_a and sigma2_a, the fixed stabilization
# score, and the stabilized flag. Replicates ate_continuous() exactly, including
# the z-score dnorm without the 1/sigma normalization and the silent
# unstabilized branch (the registry never emits the alert).
ipw_continuous_weight_fn <- function(estimand) {
  function(ps, exposure, extras) {
    z_den <- (exposure - ps) / sqrt(extras$sigma2_d)
    f_den <- stats::dnorm(z_den)
    wt <- 1 / f_den

    if (isTRUE(extras$stabilized) && is.null(extras$score)) {
      z_num <- (exposure - extras$mu_a) / sqrt(extras$sigma2_a)
      wt * stats::dnorm(z_num)
    } else if (isTRUE(extras$stabilized) && !is.null(extras$score)) {
      wt * extras$score
    } else {
      wt
    }
  }
}

# ---- theta layout -----------------------------------------------------------

ipw_theta_layout <- function(spec) {
  blocks <- switch(
    spec$exposure_type,
    binary = ipw_init_binary(spec),
    categorical = ipw_init_categorical(spec),
    continuous = ipw_init_continuous(spec)
  )

  block_order <- c("ps", "stab", "out", "mu", "contrast")
  idx <- vector("list", length(block_order))
  names(idx) <- block_order
  pos <- 0L
  for (b in block_order) {
    n_b <- length(blocks[[b]])
    idx[[b]] <- if (n_b > 0L) pos + seq_len(n_b) else integer(0)
    pos <- pos + n_b
  }

  init <- do.call(c, unname(blocks[block_order]))

  list(idx = idx, init = init)
}

# Value of a single contrast form from a pair of marginal means.
ipw_contrast_value <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo),
    "log(or)" = stats::qlogis(mu_hi) - stats::qlogis(mu_lo)
  )
}

# Plug-in contrast init values, named by form (with an optional level suffix for
# categorical exposures).
ipw_init_contrasts <- function(contrasts, mu_hi, mu_lo, suffix = NULL) {
  if (is.null(contrasts) || !length(contrasts)) {
    return(numeric(0))
  }
  vals <- vapply(
    contrasts,
    function(f) ipw_contrast_value(f, mu_hi, mu_lo),
    numeric(1)
  )
  nm <- if (is.null(suffix)) contrasts else paste0(contrasts, "_", suffix)
  stats::setNames(vals, nm)
}

ipw_init_binary <- function(spec) {
  ps_block <- spec$ps$coefs

  if (spec$stab$stabilized && is.null(spec$stab$score)) {
    stab_block <- c(stab_pi = mean(spec$exposure))
  } else {
    stab_block <- numeric(0)
  }

  beta <- spec$outcome$coefs
  inv_out <- ipw_inv_link(spec$outcome$link)
  mu1 <- mean(inv_out(as.vector(spec$outcome$X_counterfactual$X1 %*% beta)))
  mu0 <- mean(inv_out(as.vector(spec$outcome$X_counterfactual$X0 %*% beta)))
  mu_block <- c(mu1 = mu1, mu0 = mu0)

  con_block <- ipw_init_contrasts(spec$contrasts, mu1, mu0)

  list(
    ps = ps_block,
    stab = stab_block,
    out = beta,
    mu = mu_block,
    contrast = con_block
  )
}

# Generate readable names for the multinomial ps coefficient block. The block is
# as.vector(t(coef(multinom))): level-major, term-minor.
ipw_name_categorical_ps <- function(spec) {
  terms <- colnames(spec$ps$X)
  levs <- names(spec$outcome$X_counterfactual)
  nonref <- levs[-1]
  nm <- as.vector(vapply(
    nonref,
    function(l) paste0(l, ":", terms),
    character(length(terms))
  ))
  stats::setNames(spec$ps$coefs, nm)
}

ipw_init_categorical <- function(spec) {
  ps_block <- ipw_name_categorical_ps(spec)

  levs <- names(spec$outcome$X_counterfactual)

  if (spec$stab$stabilized && is.null(spec$stab$score)) {
    props <- colMeans(spec$exposure)
    stab_block <- stats::setNames(props[-1], paste0("stab_", levs[-1]))
  } else {
    stab_block <- numeric(0)
  }

  beta <- spec$outcome$coefs
  inv_out <- ipw_inv_link(spec$outcome$link)
  k <- spec$ps$k
  mu_vals <- vapply(
    seq_len(k),
    function(j) {
      mean(inv_out(as.vector(spec$outcome$X_counterfactual[[j]] %*% beta)))
    },
    numeric(1)
  )
  mu_block <- stats::setNames(mu_vals, paste0("mu_", levs))

  con_list <- lapply(seq_len(k)[-1], function(j) {
    ipw_init_contrasts(
      spec$contrasts,
      mu_vals[[j]],
      mu_vals[[1]],
      suffix = levs[[j]]
    )
  })
  con_block <- if (length(con_list)) do.call(c, con_list) else numeric(0)

  list(
    ps = ps_block,
    stab = stab_block,
    out = beta,
    mu = mu_block,
    contrast = con_block
  )
}

ipw_init_continuous <- function(spec) {
  alpha <- spec$ps$coefs
  fitted_ps <- as.vector(spec$ps$X %*% alpha)
  sigma2_d <- mean((spec$exposure - fitted_ps)^2)
  ps_block <- c(alpha, sigma2_d = sigma2_d)

  if (spec$stab$stabilized && is.null(spec$stab$score)) {
    mu_a <- mean(spec$exposure)
    sigma2_a <- mean((spec$exposure - mu_a)^2)
    stab_block <- c(mu_a = mu_a, sigma2_a = sigma2_a)
  } else {
    stab_block <- numeric(0)
  }

  list(
    ps = ps_block,
    stab = stab_block,
    out = spec$outcome$coefs,
    mu = numeric(0),
    contrast = numeric(0)
  )
}

# ---- psi builder ------------------------------------------------------------

build_ipw_psi <- function(spec, layout) {
  weight_fn <- ipw_weight_fn(spec$exposure_type, spec$estimand)
  switch(
    spec$exposure_type,
    binary = ipw_psi_binary(spec, layout, weight_fn),
    categorical = ipw_psi_categorical(spec, layout, weight_fn),
    continuous = ipw_psi_continuous(spec, layout, weight_fn)
  )
}

# Deterministic contrast rows: each form's residual repeated across the n
# observation columns so the stacked matrix stays rectangular.
ipw_contrast_rows <- function(contrasts, mu_hi, mu_lo, th_con, n) {
  if (is.null(contrasts) || !length(contrasts)) {
    return(NULL)
  }
  rows <- lapply(seq_along(contrasts), function(i) {
    val <- ipw_contrast_value(contrasts[[i]], mu_hi, mu_lo)
    matrix(val - th_con[[i]], nrow = 1, ncol = n)
  })
  do.call(rbind, rows)
}

ipw_stack <- function(blocks) {
  do.call(rbind, blocks[!vapply(blocks, is.null, logical(1))])
}

# Weighted outcome-score block. Branch on the outcome family so a gaussian
# outcome uses the weighted least-squares score and a binomial outcome uses the
# weighted GLM score with the outcome link.
ipw_outcome_rows <- function(theta, x, y, family, link, weights) {
  if (family == "binomial") {
    deli::ee_glm(
      theta,
      X = x,
      y = y,
      distribution = "binomial",
      link = link,
      weights = weights
    )
  } else {
    deli::ee_regression(
      theta,
      X = x,
      y = y,
      model = "linear",
      weights = weights
    )
  }
}

ipw_psi_binary <- function(spec, layout, weight_fn) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  ps_link <- spec$ps$link
  inv_ps <- ipw_inv_link(ps_link)
  z <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  inv_out <- ipw_inv_link(out_link)
  x1 <- spec$outcome$X_counterfactual$X1
  x0 <- spec$outcome$X_counterfactual$X0
  score <- spec$stab$score
  contrasts <- spec$contrasts
  n <- spec$n

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    ps_rows <- deli::ee_glm(
      th_ps,
      X = x_ps,
      y = z,
      distribution = "binomial",
      link = ps_link
    )
    e <- inv_ps(as.vector(x_ps %*% th_ps))

    stab_prob <- if (length(th_stab)) th_stab[[1]] else NULL
    w <- weight_fn(e, z, list(stab_prob = stab_prob, score = score))

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    stab_rows <- if (length(th_stab)) {
      matrix(z - th_stab[[1]], nrow = 1)
    } else {
      NULL
    }

    mu1 <- th_mu[[1]]
    mu0 <- th_mu[[2]]
    mu_rows <- rbind(
      inv_out(as.vector(x1 %*% th_out)) - mu1,
      inv_out(as.vector(x0 %*% th_out)) - mu0
    )

    con_rows <- ipw_contrast_rows(contrasts, mu1, mu0, th_con, n)

    ipw_stack(list(ps_rows, stab_rows, out_rows, mu_rows, con_rows))
  }
}

# Reconstruct the n-by-K propensity score matrix from the multinomial ps block,
# matching deli::ee_mlogit's internal softmax. Column order is reference-first.
ipw_categorical_ps <- function(x, theta, k) {
  n <- nrow(x)
  b <- ncol(x)
  denom <- rep(1, n)
  exp_pred <- matrix(0, nrow = n, ncol = k - 1)
  for (j in seq_len(k - 1)) {
    idx <- ((j - 1) * b + 1):(j * b)
    exp_pred[, j] <- exp(as.vector(x %*% theta[idx]))
    denom <- denom + exp_pred[, j]
  }
  ps <- matrix(0, nrow = n, ncol = k)
  ps[, 1] <- 1 / denom
  for (j in seq_len(k - 1)) {
    ps[, j + 1] <- exp_pred[, j] / denom
  }
  ps
}

ipw_psi_categorical <- function(spec, layout, weight_fn) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  k <- spec$ps$k
  z_ind <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  inv_out <- ipw_inv_link(out_link)
  x_cf <- spec$outcome$X_counterfactual
  levs <- names(x_cf)
  score <- spec$stab$score
  focal_idx <- if (!is.null(spec$focal_level)) {
    match(spec$focal_level, levs)
  } else {
    NULL
  }
  contrasts <- spec$contrasts
  n <- spec$n

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]
    th_mu <- theta[idx$mu]
    th_con <- theta[idx$contrast]

    ps_rows <- deli::ee_mlogit(th_ps, X = x_ps, y = z_ind)
    ps_mat <- ipw_categorical_ps(x_ps, th_ps, k)

    stab_probs <- if (length(th_stab)) {
      c(1 - sum(th_stab), th_stab)
    } else {
      NULL
    }
    w <- weight_fn(
      ps_mat,
      z_ind,
      list(focal_idx = focal_idx, stab_probs = stab_probs, score = score)
    )

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    stab_rows <- if (length(th_stab)) {
      do.call(
        rbind,
        lapply(seq_len(k - 1), function(j) {
          matrix(z_ind[, j + 1] - th_stab[[j]], nrow = 1)
        })
      )
    } else {
      NULL
    }

    mu_pred <- vapply(
      seq_len(k),
      function(j) inv_out(as.vector(x_cf[[j]] %*% th_out)),
      numeric(n)
    )
    mu_rows <- t(mu_pred) - th_mu

    con_rows <- if (!is.null(contrasts) && length(contrasts)) {
      mu_ref <- th_mu[[1]]
      con_index <- 0L
      row_list <- list()
      for (j in seq_len(k)[-1]) {
        mu_j <- th_mu[[j]]
        for (form in contrasts) {
          con_index <- con_index + 1L
          val <- ipw_contrast_value(form, mu_j, mu_ref)
          row_list[[con_index]] <- matrix(
            val - th_con[[con_index]],
            nrow = 1,
            ncol = n
          )
        }
      }
      do.call(rbind, row_list)
    } else {
      NULL
    }

    ipw_stack(list(ps_rows, stab_rows, out_rows, mu_rows, con_rows))
  }
}

ipw_psi_continuous <- function(spec, layout, weight_fn) {
  idx <- layout$idx
  x_ps <- spec$ps$X
  a <- spec$exposure
  x_out <- spec$outcome$X
  y <- spec$outcome$y
  family <- spec$outcome$family
  out_link <- spec$outcome$link
  score <- spec$stab$score
  stabilized <- spec$stab$stabilized
  n_alpha <- ncol(x_ps)

  function(theta) {
    th_ps <- theta[idx$ps]
    th_stab <- theta[idx$stab]
    th_out <- theta[idx$out]

    alpha <- th_ps[seq_len(n_alpha)]
    sigma2_d <- th_ps[[n_alpha + 1]]

    ps_score <- deli::ee_regression(alpha, X = x_ps, y = a, model = "linear")
    fitted_ps <- as.vector(x_ps %*% alpha)
    var_row <- matrix((a - fitted_ps)^2 - sigma2_d, nrow = 1)
    ps_rows <- rbind(ps_score, var_row)

    if (length(th_stab)) {
      mu_a <- th_stab[[1]]
      sigma2_a <- th_stab[[2]]
      stab_rows <- rbind(
        matrix(a - mu_a, nrow = 1),
        matrix((a - mu_a)^2 - sigma2_a, nrow = 1)
      )
    } else {
      mu_a <- NULL
      sigma2_a <- NULL
      stab_rows <- NULL
    }

    w <- weight_fn(
      fitted_ps,
      a,
      list(
        sigma2_d = sigma2_d,
        mu_a = mu_a,
        sigma2_a = sigma2_a,
        score = score,
        stabilized = stabilized
      )
    )

    out_rows <- ipw_outcome_rows(th_out, x_out, y, family, out_link, w)

    ipw_stack(list(ps_rows, stab_rows, out_rows))
  }
}
