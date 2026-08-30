# a multinomial fit reads its own exposure and says which it read

    Code
      from_model <- wt_ate(fit)
    Message
      i Using exposure variable "trt" from the propensity score model
      i Treating `.exposure` as categorical

# a two-level fit refuses a categorical exposure of three levels

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a categorical exposure need a probability for every level of the exposure being weighted.
      x `.propensity` was fit to 2 levels ("control" and "treated"), and the exposure has 3 ("a", "b", and "c").
      i Fit the propensity score model to the exposure being weighted.

# a four-level fit refuses a categorical exposure of three levels

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a categorical exposure need a probability for every level of the exposure being weighted.
      x `.propensity` was fit to 4 levels ("a", "b", "c", and "d"), and the exposure has 3 ("a", "b", and "c").
      i Fit the propensity score model to the exposure being weighted.

# an exposure carrying an unused level is sent to droplevels()

    Code
      expr
    Condition <propensity_model_family_error>
      Error in `wt_ate()`:
      ! Weights for a categorical exposure need a probability for every level of the exposure being weighted.
      x `.propensity` was fit to 3 levels ("a", "b", and "c"), and the exposure has 4 ("a", "b", "c", and "d").
      i `nnet::multinom()` drops a level no unit took, so an exposure declaring "d" with no unit there is fit one column short; call `droplevels()` on `.exposure` if that is the disagreement.
      i Fit the propensity score model to the exposure being weighted.

