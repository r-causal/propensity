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

