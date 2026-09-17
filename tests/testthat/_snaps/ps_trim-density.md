# ps_trim() refuses an upper bound on a density floor

    Code
      expr
    Condition <propensity_unsupported_arg_error>
      Error in `ps_trim()`:
      ! `upper` is not used by `method = "density"`.
      x The method sets aside the units whose conditional density falls below a quantile floor, and a large conditional density is a plausible dose, so there is no upper tail to trim.
      i Drop `upper`, and set the floor with `lower`.

# ps_trim() refuses a lower bound on a residual bound, and needs an upper one

    Code
      expr
    Condition <propensity_unsupported_arg_error>
      Error in `ps_trim()`:
      ! `lower` is not used by `method = "resid"`.
      x The method bounds the absolute standardized residual, which has no lower end to trim.
      i Drop `lower`, and set the bound with `upper`.

---

    Code
      expr
    Condition <propensity_missing_arg_error>
      Error in `ps_trim()`:
      ! `method = "resid"` needs `upper`.
      x The bound on the absolute standardized residual has no default.
      i Supply the number of spreads from the predicted dose beyond which a unit is trimmed, such as `upper = 3`.

# ps_trim() refuses a density floor outside (0, 0.5) and a residual bound that is not positive

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trim()`:
      ! For `method = "density"`, `lower` must be a single probability between 0 and 0.5.
      x `lower` is 0.5.
      i `lower` is the quantile of the conditional density below which units are trimmed, such as the default `lower = 0.01`.

---

    Code
      expr
    Condition <propensity_range_error>
      Error in `ps_trim()`:
      ! For `method = "resid"`, `upper` must be a single positive, finite number.
      x `upper` is -2.
      i `upper` is the number of spreads from the predicted dose beyond which a unit is trimmed, such as `upper = 3`.

# ps_trim() refuses a residual bound under a density that need not be symmetric

    Code
      expr
    Condition <propensity_density_error>
      Error in `ps_trim()`:
      ! `method = "resid"` needs a symmetric unimodal density.
      x `.density` is a "kernel" density, which need not be symmetric or unimodal, so a bound on the absolute residual is not a floor on it.
      i Use `method = "density"`, which trims on the density itself under any family, or a `dens_normal()`, `dens_t()`, or `dens_laplace()` density.

# ps_trim() refuses a spread it cannot record on a dose model

    Code
      expr
    Condition <propensity_sigma_error>
      Error in `ps_trim()`:
      ! `.sigma` must be a single spread when trimming a dose model.
      x It holds 60 values.
      i A trim is recorded at one spread, which a refit of the dose model can keep; a spread for each unit cannot be refit. Supply one value, or leave `.sigma` unset to estimate it.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `ps_trim()`:
      ! `sigma_method = "mle"` cannot be used with `.sigma`.
      x `sigma_method = "mle"` estimates the scale of the "t" conditional density from the residuals of the model of the exposure, and `.sigma` is a spread of your own that replaces it.
      i A "t" density is spread by its own estimator unless you ask for another, so `sigma_method = "mle"` is what it carries whether or not you wrote it.
      i Drop `.sigma` to estimate the scale under the "t" density, or build the density with `sigma_method = "rms"` to spread the one you supplied.

# ps_trim() refuses a dose that does not line up with the fitted means

    Code
      expr
    Condition <propensity_length_error>
      Error in `ps_trim()`:
      ! `.exposure` must have one value for each fitted mean.
      x It has 60 values, and the model reports 58 fitted means.
      i Supply the dose for the rows the model was fit to, or refit it with `na.action = na.exclude`, which reports a mean for every row of its data.

---

    Code
      expr
    Condition <propensity_type_error>
      Error in `ps_trim()`:
      ! `.exposure` must be a numeric vector when trimming a dose model.
      x It is a character vector.
      i Supply the dose the model was fit to, one value per unit.

# ps_trim() refuses a dose model of several responses

    Code
      expr
    Condition <propensity_ps_shape_error>
      Error in `ps_trim()`:
      ! Trimming a dose model needs a model of one conditional mean for each unit.
      x `.propensity` is <mlm>, a fit of 2 responses, whose fitted values are 60 by 2.
      i Fit the dose on its own and trim that model.

# ps_trim() refuses a dose model fit with a spread that moves with its mean

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `ps_trim()`:
      ! Trimming a dose model needs a model of its conditional mean with a single spread.
      x `.propensity` was fit with `poisson()`, whose spread changes with its fitted values.
      i Fit the dose model with `gaussian()`, `lm()`, `mgcv::gam()`, or `MASS::rlm()`.

# ps_trim() refuses a score method on a dose model and names the dose methods

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "ps" cannot trim a dose model.
      x `.propensity` is <lm>, a model of a continuous exposure's conditional mean, which has no propensity score in (0, 1) for "ps" to bound.
      i Trim it on the scale of its conditional density with `method = "density"` or `method = "resid"`.

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "adaptive" cannot trim a dose model.
      x `.propensity` is <glm>, a model of a continuous exposure's conditional mean, which has no propensity score in (0, 1) for "adaptive" to bound.
      i Trim it on the scale of its conditional density with `method = "density"` or `method = "resid"`.

# ps_trim() refuses levels on a dose model

    Code
      expr
    Condition <propensity_focal_level_error>
      Error in `ps_trim()`:
      ! `.focal_level` cannot be used when trimming a dose model.
      x A continuous exposure has no levels, so none is focal and none is reference.
      i Drop `.focal_level`.

# ps_trim() refuses a density or a spread on the score routes

    Code
      expr
    Condition <propensity_density_error>
      Error in `ps_trim()`:
      ! `.density` applies only to trimming a dose model.
      x The trimming method bounds propensity scores, which have no conditional density to read.
      i Leave `.density` unset, or supply the model of a continuous exposure with `method = "density"` or `method = "resid"`.

---

    Code
      expr
    Condition <propensity_sigma_error>
      Error in `ps_trim()`:
      ! `.sigma` applies only to trimming a dose model.
      x The trimming method bounds propensity scores, which have no spread to read.
      i Leave `.sigma` unset, or supply the model of a continuous exposure with `method = "density"` or `method = "resid"`.

# ps_trim() refuses the density methods on a model of a probability

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "density" cannot trim a model of a probability.
      x `.propensity` is <glm>, whose fitted values are probabilities of the exposure, and a model of a probability has no conditional density to trim.
      i Trim its propensity scores with a score method such as "ps", or supply the model of a continuous exposure's conditional mean to trim with "density".

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "resid" cannot trim a model of a probability.
      x `.propensity` is <glm>, whose fitted values are probabilities of the exposure, and a model of a probability has no conditional density to trim.
      i Trim its propensity scores with a score method such as "ps", or supply the model of a continuous exposure's conditional mean to trim with "resid".

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "density" cannot trim a model of a probability.
      x `.propensity` is <multinom>, whose fitted values are probabilities of the exposure, and a model of a probability has no conditional density to trim.
      i Trim its propensity scores with a score method such as "ps", or supply the model of a continuous exposure's conditional mean to trim with "density".

---

    Code
      expr
    Condition <propensity_method_error>
      Error in `ps_trim()`:
      ! Method "resid" cannot trim a model of a probability.
      x `.propensity` is <multinom>, whose fitted values are probabilities of the exposure, and a model of a probability has no conditional density to trim.
      i Trim its propensity scores with a score method such as "ps", or supply the model of a continuous exposure's conditional mean to trim with "resid".

