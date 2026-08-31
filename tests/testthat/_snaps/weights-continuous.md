# the refusal of a density reads the same way for either type

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A binary exposure has a probability rather than a density, so leave `.density` unset for one.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as categorical.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A categorical exposure has a probability rather than a density, so leave `.density` unset for one.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_cens()`:
      ! `.density` applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i `.density` chooses the family of the conditional density a continuous exposure's weights are a ratio of. A binary exposure has a probability rather than a density, so leave `.density` unset for one.

# the refusal of a matrix of conditional means reads plainly

    Code
      expr
    Condition <propensity_ps_shape_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure need one conditional mean for each unit.
      x `.propensity` is 60 by 2 and of class <matrix>.
      i Pass a numeric vector of conditional means, such as the single column of `.propensity` that holds the mean of this exposure.

---

    Code
      expr
    Condition <propensity_ps_shape_error>
      Error in `wt_cens()`:
      ! Weights for a continuous exposure need one conditional mean for each unit.
      x `.propensity` is 60 by 2 and of class <matrix>.
      i Pass a numeric vector of conditional means, such as the single column of `.propensity` that holds the mean of this exposure.

---

    Code
      expr
    Condition <propensity_ps_shape_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure need a model of one conditional mean for each unit.
      x `.propensity` is <mlm>, a fit of 2 responses, whose fitted values are 60 by 2.
      i Fit the exposure on its own, or pass the conditional means of this exposure to `.propensity` as a numeric vector.

---

    Code
      expr
    Condition <propensity_ps_shape_error>
      Error in `wt_cens()`:
      ! Weights for a continuous exposure need a model of one conditional mean for each unit.
      x `.propensity` is <mlm>, a fit of 2 responses, whose fitted values are 60 by 2.
      i Fit the exposure on its own, or pass the conditional means of this exposure to `.propensity` as a numeric vector.

# an infinite exposure or fitted value is refused where it arrives

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure cannot be computed from an infinite value.
      x 1 value of `.exposure` is infinite.
      i An infinite value leaves the spread of the conditional density infinite, so every weight is missing rather than that unit's alone. Drop those observations before weighting.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure cannot be computed from an infinite value.
      x 1 value of `.propensity` is infinite.
      i An infinite value leaves the spread of the conditional density infinite, so every weight is missing rather than that unit's alone. Drop those observations before weighting.

# the zero-density refusal reads the way it is written

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! The conditional density is zero where the weight would divide by it.
      x It is zero at 1 observation: 33.
      i It failed at the standardized residual -2.19.
      i A weight built from a conditional density of zero is infinite. Use a family with heavier tails, such as `dens_t()` or `dens_kernel()`.

# the finite-variance report reads the way it is written

    Code
      out <- wt_ate(past$mu, past$exposure, exposure_type = "continuous", stabilize = TRUE)
    Message <propensity_density_variance_message>
      Stabilized normal weights have no finite variance for this exposure.
      x The marginal variance of the exposure is 2.6, which is at least twice the conditional variance 1.09.
      i The second moment of a stabilized normal density ratio exists only while the marginal variance stays below twice the conditional one. The weights are returned, but estimates built from them are erratic however large the sample is.
      i The boundary tightens as the model explains more of the exposure, since a better fit lowers the conditional variance while the marginal variance is fixed by the data. A family with heavier tails, such as `dens_t()`, has a different boundary, and an unstabilized weight has none.

# the refusals of an integrated numerator read the way they should

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i An integrated numerator averages the conditional density of a continuous exposure over the units. A binary exposure has a probability rather than a density, so leave `numerator` unset for one.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" applies only to continuous exposures.
      x `.exposure` is being treated as categorical.
      i An integrated numerator averages the conditional density of a continuous exposure over the units. A categorical exposure has a probability rather than a density, so leave `numerator` unset for one.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" needs stabilized weights.
      x `stabilize` is `FALSE`, so the weights carry no numerator at all.
      i Set `stabilize = TRUE` to stabilize the weights on the marginalized conditional density, or leave `numerator` unset.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" cannot be used with `stabilization_score`.
      x A score you supply is itself the numerator of the weights.
      i Drop `stabilization_score` to stabilize on the marginalized conditional density, or leave `numerator` unset to keep the numerator you wrote.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" cannot be used with `.sigma`.
      x The marginalization reads the conditional density the propensity model estimated at every unit's fitted mean, and `.sigma` replaces the spread of that density with one of your own.
      i Leave `.sigma` unset to spread the conditional density by the pooled residual root mean square, or use `numerator` = "marginal", which takes a spread you supply.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" cannot marginalize over an exposure that does not vary.
      x Every exposure with a fitted conditional mean is 2, so the grid the conditional density is averaged over has no width.
      i The integrated numerator reads the conditional density at points spanning `.exposure`. Use `numerator` = "marginal" for an exposure that takes one value.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! The integrated numerator must return finite, non-negative values.
      x It returned 3 values that are negative or not finite.
      i It failed at standardized residuals from -1.5 to -1.27.
      i The marginalized density is interpolated back to the exposure with a cubic spline, which can dip below zero where the density on the grid comes close to it. Use `numerator` = "marginal", or a density with a heavier tail, to stabilize on a density that is positive everywhere.

# a spread supplied and a spread estimated under the density are refused together

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! `sigma_method = "mle"` cannot be used with `.sigma`.
      x `sigma_method = "mle"` estimates the scale of the conditional density from the residuals of the propensity score model, and `.sigma` is a spread of your own that replaces it.
      i Drop `.sigma` to estimate the scale under the t density, or build the density with `sigma_method = "rms"` to spread the one you supplied.

# a likelihood with no maximum at a positive scale is refused

    Code
      expr
    Condition <propensity_density_error>
      Error in `wt_ate()`:
      ! The scale of a "t" density cannot be estimated by maximum likelihood from these residuals.
      x 9 of the 10 residuals of the propensity score model are exactly zero, and with 4 degrees of freedom that leaves the likelihood no maximum at a positive scale.
      i Use `sigma_method = "rms"`, or supply a spread with `.sigma`.

# the numerator model refusals read as the refusals they are

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! A model supplied to `stabilize` applies only to continuous exposures.
      x `.exposure` is being treated as binary.
      i A fitted model stabilizes the weights on the conditional density it estimates. A binary exposure has a probability rather than a density, so stabilize one with `stabilize = TRUE`, which reads the marginal probability of the exposure, or with a `stabilization_score` of your own.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! A model supplied to `stabilize` cannot be used with `stabilization_score`.
      x A score you supply is itself the numerator of the weights, and the model estimates a second one.
      i Drop `stabilization_score` to stabilize on the density the model estimates, or set `stabilize = TRUE` to keep the numerator you wrote.

---

    Code
      expr
    Condition <propensity_numerator_error>
      Error in `wt_ate()`:
      ! `numerator` = "integrated" cannot be used with a model supplied to `stabilize`.
      x A model you supply estimates the numerator of the weights itself.
      i Set `stabilize = TRUE` to stabilize on the marginalized conditional density, or leave `numerator` unset to keep the numerator the model estimates.

---

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a continuous exposure need a model of its conditional mean with a single spread.
      x `stabilize` was fit with `binomial()`, whose spread changes with its fitted values.
      i The numerator is a density read at the model's fitted mean and the spread of its residuals, so refit it with `stats::lm()` or `stats::glm(family = gaussian())`.

---

    Code
      expr
    Condition <propensity_length_error>
      Error in `wt_ate()`:
      ! The model supplied to `stabilize` must have one fitted value for each observation.
      x It has 20 fitted values and `.exposure` has 60 observations.
      i Fit the numerator model on the data the weights are being built for.

# stabilize takes one of the answers or a model and nothing else

    Code
      expr
    Condition <propensity_stabilize_error>
      Error in `wt_ate()`:
      ! `stabilize` must be `TRUE`, `FALSE`, `NULL`, or a fitted model of the exposure.
      x It is a string.
      i A fitted model stabilizes a continuous exposure's weights on the conditional density it estimates. To stabilize on a numerator you computed yourself, set `stabilize = TRUE` and pass it as `stabilization_score`.

---

    Code
      expr
    Condition <propensity_stabilize_error>
      Error in `wt_ate()`:
      ! `stabilize` must be `TRUE`, `FALSE`, `NULL`, or a fitted model of the exposure.
      x It is `NA`.
      i A fitted model stabilizes a continuous exposure's weights on the conditional density it estimates. To stabilize on a numerator you computed yourself, set `stabilize = TRUE` and pass it as `stabilization_score`.

