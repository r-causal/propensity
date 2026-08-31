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
    Condition <propensity_density_variance_warning>
      Warning in `wt_ate()`:
      Stabilized normal weights have no finite variance for this exposure.
      x The marginal variance of the exposure is 2.6, which is at least twice the conditional variance 1.09.
      i The second moment of a stabilized normal density ratio exists only while the marginal variance stays below twice the conditional one. The weights are returned, but estimates built from them are erratic however large the sample is.
      i A model that explains more of the exposure lowers the conditional variance, and a family with heavier tails, such as `dens_t()`, has a different boundary.

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

